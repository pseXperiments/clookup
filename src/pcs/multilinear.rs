use crate::{
    poly::multilinear::MultilinearPolynomial,
    utils::{
        arithmetic::Field, end_timer, izip, parallel::parallelize, start_timer, Itertools,
        ProtocolError,
    },
};

pub mod kzg;

fn validate_input<'a, F: Field>(
    function: &str,
    param_num_vars: usize,
    polys: impl IntoIterator<Item = &'a MultilinearPolynomial<F>>,
    points: impl IntoIterator<Item = &'a Vec<F>>,
) -> Result<(), ProtocolError> {
    let polys = polys.into_iter().collect_vec();
    let points = points.into_iter().collect_vec();
    for poly in polys.iter() {
        if param_num_vars < poly.num_vars() {
            return Err(err_too_many_variates(
                function,
                param_num_vars,
                poly.num_vars(),
            ));
        }
    }
    let input_num_vars = polys
        .iter()
        .map(|poly| poly.num_vars())
        .chain(points.iter().map(|point| point.len()))
        .next()
        .expect("To have at least 1 poly or point");
    for point in points.into_iter() {
        if point.len() != input_num_vars {
            return Err(ProtocolError::InvalidPcsParam(format!(
                "Invalid point (expect point to have {input_num_vars} variates but got {})",
                point.len()
            )));
        }
    }
    Ok(())
}

fn err_too_many_variates(function: &str, upto: usize, got: usize) -> ProtocolError {
    ProtocolError::InvalidPcsParam(if function == "trim" {
        format!(
            "Too many variates to {function} (param supports variates up to {upto} but got {got})"
        )
    } else {
        format!(
            "Too many variates of poly to {function} (param supports variates up to {upto} but got {got})"
        )
    })
}

fn quotients<F: Field, T>(
    poly: &MultilinearPolynomial<F>,
    point: &[F],
    f: impl Fn(usize, Vec<F>) -> T,
) -> (Vec<T>, F) {
    assert_eq!(poly.num_vars(), point.len());

    let mut remainder = poly.evals().to_vec();
    let mut quotients = point
        .iter()
        .zip(0..poly.num_vars())
        .rev()
        .map(|(x_i, num_vars)| {
            let timer = start_timer(|| "quotients");
            let (remaimder_lo, remainder_hi) = remainder.split_at_mut(1 << num_vars);
            let mut quotient = vec![F::ZERO; remaimder_lo.len()];

            parallelize(&mut quotient, |(quotient, start)| {
                izip!(quotient, &remaimder_lo[start..], &remainder_hi[start..])
                    .for_each(|(q, r_lo, r_hi)| *q = *r_hi - r_lo);
            });
            parallelize(remaimder_lo, |(remaimder_lo, start)| {
                izip!(remaimder_lo, &remainder_hi[start..])
                    .for_each(|(r_lo, r_hi)| *r_lo += (*r_hi - r_lo as &_) * x_i);
            });

            remainder.truncate(1 << num_vars);
            end_timer(timer);

            f(num_vars, quotient)
        })
        .collect_vec();
    quotients.reverse();

    (quotients, remainder[0])
}

mod additive {
    use ff::Field;
    use itertools::Itertools;
    use std::{borrow::Cow, ops::Deref, ptr::addr_of};

    use crate::{
        pcs::{Additive, Evaluation, Point, PolynomialCommitmentScheme},
        poly::multilinear::MultilinearPolynomial,
        sumcheck::{
            classic::{ClassicSumcheck, ClassicSumcheckProverParam, ClassicSumcheckVerifierParam},
            eq_xy_eval, SumCheck as _, VirtualPolynomial,
        },
        utils::{
            arithmetic::{inner_product, PrimeField},
            end_timer, start_timer,
            transcript::{TranscriptRead, TranscriptWrite},
            ProtocolError,
        },
    };

    use super::validate_input;

    type SumCheck = ClassicSumcheck;

    fn generate_additive_comm_fn<'a, F: Field>(
        info: &'a Vec<(usize, &F, &MultilinearPolynomial<F>)>,
        eq_size: usize,
    ) -> impl Fn(&Vec<F>) -> F + 'a {
        move |evals: &Vec<F>| {
            let eq_res = evals
                .iter()
                .skip(evals.len() - eq_size)
                .take(eq_size)
                .collect_vec();
            info
                .iter()
                .zip(evals.iter())
                .map(|((idx, scalar, _), eval)| {
                    eval.clone() * (*scalar).clone() * eq_res[*idx]
                })
                .sum()
        }
    }

    pub fn batch_open<F, Pcs>(
        pp: &Pcs::ProverParam,
        num_vars: usize,
        polys: Vec<&Pcs::Polynomial>,
        comms: Vec<&Pcs::Commitment>,
        points: &[Point<F, Pcs::Polynomial>],
        evals: &[Evaluation<F>],
        transcript: &mut impl TranscriptWrite<Pcs::CommitmentChunk, F>,
    ) -> Result<(), ProtocolError>
    where
        F: PrimeField,
        Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
        Pcs::Commitment: Additive<F>,
    {
        validate_input("batch open", num_vars, polys.clone(), points)?;

        let ell = evals.len().next_power_of_two().ilog2() as usize;
        let t = transcript.squeeze_challenges(ell);

        let timer = start_timer(|| "merged_polys");
        let eq_xt = MultilinearPolynomial::eq_xy(&t);

        let merged_polys = evals.iter().zip(eq_xt.evals().iter()).fold(
            vec![(F::ONE, Cow::<MultilinearPolynomial<_>>::default()); points.len()],
            |mut merged_polys, (eval, eq_xt_i)| {
                if merged_polys[eval.point()].1.is_empty() {
                    merged_polys[eval.point()] = (*eq_xt_i, Cow::Borrowed(polys[eval.poly()]));
                } else {
                    let coeff = merged_polys[eval.point()].0;
                    if coeff != F::ONE {
                        merged_polys[eval.point()].0 = F::ONE;
                        *merged_polys[eval.point()].1.to_mut() *= &coeff;
                    }
                    *merged_polys[eval.point()].1.to_mut() += (eq_xt_i, polys[eval.poly()]);
                }
                merged_polys
            },
        );
        end_timer(timer);

        let used_polys = merged_polys
            .iter()
            .enumerate()
            .map(|(idx, (scalar, poly))| (idx, scalar, poly.deref()))
            .collect_vec();

        let eq_xys = points
            .iter()
            .map(|y| MultilinearPolynomial::eq_xy(y))
            .collect_vec();
        let combine_function = generate_additive_comm_fn(&used_polys, eq_xys.len());
        let mut virtual_polys = used_polys.iter().map(|(_, _, poly)| *poly).collect_vec();
        for p in eq_xys.iter() {
            virtual_polys.push(p);
        }
        let virtual_poly = VirtualPolynomial::new(num_vars, virtual_polys.as_slice());

        let tilde_gs_sum =
            inner_product(evals.iter().map(Evaluation::value), &eq_xt[..evals.len()]);
        let spp = ClassicSumcheckProverParam::new(num_vars, 2);
        let (challenges, _) = SumCheck::prove(
            &spp,
            &combine_function,
            tilde_gs_sum,
            virtual_poly,
            transcript,
        )?;

        let timer = start_timer(|| "g_prime");
        let eq_xy_evals = points
            .iter()
            .map(|point| eq_xy_eval(&challenges, point))
            .collect_vec();
        let g_prime = merged_polys
            .iter()
            .zip(eq_xy_evals.iter())
            .map(|((scalar, poly), eq_xy_eval)| {
                (scalar.clone() * eq_xy_eval, poly.clone().into_owned())
            })
            .sum::<MultilinearPolynomial<_>>();
        end_timer(timer);

        let (g_prime_comm, g_prime_eval) = if cfg!(feature = "sanity-check") {
            let scalars = evals
                .iter()
                .zip(eq_xt.evals())
                .map(|(eval, eq_xt_i)| eq_xy_evals[eval.point()] * eq_xt_i)
                .collect_vec();
            let bases = evals.iter().map(|eval| comms[eval.poly()]);
            let comm = Pcs::Commitment::msm(&scalars, bases);
            (comm, g_prime.evaluate(&challenges))
        } else {
            (Pcs::Commitment::default(), F::ZERO)
        };
        Pcs::open(
            pp,
            &g_prime,
            &g_prime_comm,
            &challenges,
            &g_prime_eval,
            transcript,
        )
    }

    pub fn batch_verify<F, Pcs>(
        vp: &Pcs::VerifierParam,
        num_vars: usize,
        comms: Vec<&Pcs::Commitment>,
        points: &[Point<F, Pcs::Polynomial>],
        evals: &[Evaluation<F>],
        transcript: &mut impl TranscriptRead<Pcs::CommitmentChunk, F>,
    ) -> Result<(), ProtocolError>
    where
        F: PrimeField,
        Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
        Pcs::Commitment: Additive<F>,
    {
        validate_input("batch verify", num_vars, [], points)?;

        let ell = evals.len().next_power_of_two().ilog2() as usize;
        let t = transcript.squeeze_challenges(ell);

        let eq_xt = MultilinearPolynomial::eq_xy(&t);
        let tilde_gs_sum =
            inner_product(evals.iter().map(Evaluation::value), &eq_xt[..evals.len()]);
        let svp = ClassicSumcheckVerifierParam::new(num_vars, 2);
        let num_polys = points.len() + evals.len();
        let (g_prime_eval, _, challenges) =
            SumCheck::verify(&svp, 2, tilde_gs_sum, num_polys, transcript)?;
        let eq_xy_evals = points
            .iter()
            .map(|point| eq_xy_eval(&challenges, point))
            .collect_vec();
        let g_prime_comm = {
            let scalars = evals
                .iter()
                .zip(eq_xt.evals())
                .map(|(eval, eq_xt_i)| eq_xy_evals[eval.point()] * eq_xt_i)
                .collect_vec();
            let bases = evals.iter().map(|eval| comms[eval.poly()]);
            Pcs::Commitment::msm(&scalars, bases)
        };
        Pcs::verify(vp, &g_prime_comm, &challenges, &g_prime_eval, transcript)
    }
}
