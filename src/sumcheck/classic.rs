use super::{SumCheck, VirtualPolynomial};
use crate::utils::transcript::{FieldTranscriptRead, FieldTranscriptWrite};
use crate::utils::{
    arithmetic::{barycentric_interpolate, barycentric_weights},
    ProtocolError,
};
use ff::PrimeField;
use std::fmt::Debug;

#[derive(Clone, Debug)]
pub struct ClassicSumcheck;

#[derive(Clone, Debug)]
pub struct ClassicSumcheckProverParam {
    num_vars: usize,
    max_degree: usize,
}

impl ClassicSumcheckProverParam {
    pub fn new(num_vars: usize, max_degree: usize) -> Self {
        ClassicSumcheckProverParam {
            num_vars,
            max_degree,
        }
    }
}

#[derive(Clone, Debug)]
pub struct ClassicSumcheckVerifierParam {
    num_vars: usize,
    max_degree: usize,
}

impl ClassicSumcheckVerifierParam {
    pub fn new(num_vars: usize, max_degree: usize) -> Self {
        ClassicSumcheckVerifierParam {
            num_vars,
            max_degree,
        }
    }
}

impl<F: PrimeField> SumCheck<F> for ClassicSumcheck {
    type ProverParam = ClassicSumcheckProverParam;
    type VerifierParam = ClassicSumcheckVerifierParam;

    fn prove(
        pp: &Self::ProverParam,
        combine_function: &impl Fn(&Vec<F>) -> F,
        _: F,
        mut virtual_poly: VirtualPolynomial<F>,
        transcript: &mut impl FieldTranscriptWrite<F>,
    ) -> Result<(Vec<F>, Vec<F>), ProtocolError> {
        // Declare r_polys and initialise it with 0s
        let r_degree = pp.max_degree;
        let mut r_polys: Vec<Vec<F>> = (0..pp.num_vars)
            .map(|_| vec![F::ZERO; r_degree + 1])
            .collect();

        let mut challenges = vec![];
        let mut evaluations = vec![];
        for round_index in 0..pp.num_vars {
            for k in 0..(r_degree + 1) {
                for i in 0..virtual_poly.polys()[0].size() {
                    let evaluations_at_k = virtual_poly
                        .polys()
                        .iter()
                        .map(|poly| {
                            let o = poly.table()[i].odd;
                            let e = poly.table()[i].even;
                            e + F::from(k as u64) * (o - e)
                        })
                        .collect::<Vec<F>>();

                    // apply combine function
                    r_polys[round_index][k] += combine_function(&evaluations_at_k);
                }
            }
            // append the round polynomial (i.e. prover message) to the transcript
            transcript.write_field_elements(&r_polys[round_index])?;

            // generate challenge α_i = H( transcript );
            let alpha = transcript.squeeze_challenge();
            challenges.push(alpha);

            if round_index == pp.num_vars - 1 {
                // last round
                evaluations = virtual_poly.evaluations(alpha);
                transcript.write_field_elements(&evaluations)?;
            } else {
                // update prover state polynomials
                virtual_poly.fold_into_half(alpha);
            }
        }

        for i in 0..virtual_poly.polys().len() {
            assert_eq!(virtual_poly.polys()[i].size(), 1);
        }
        challenges.reverse();

        Ok((challenges, evaluations))
    }

    fn verify(
        vp: &Self::VerifierParam,
        _: usize,
        sum: F,
        num_polys: usize,
        transcript: &mut impl FieldTranscriptRead<F>,
    ) -> Result<(F, Vec<F>, Vec<F>), ProtocolError> {
        let (msgs, mut challenges) = {
            let mut msgs = Vec::with_capacity(vp.num_vars);
            let mut challenges = Vec::with_capacity(vp.num_vars);
            for _ in 0..vp.num_vars {
                msgs.push(transcript.read_field_elements(vp.max_degree + 1)?);
                challenges.push(transcript.squeeze_challenge());
            }
            (msgs, challenges)
        };

        let evaluations = transcript.read_field_elements(num_polys)?;
        let mut expected_sum = sum.clone();
        let points_vec: Vec<F> = (0..vp.max_degree + 1)
            .map(|i| F::from_u128(i as u128))
            .collect();
        let weights = barycentric_weights(&points_vec);

        for round_index in 0..vp.num_vars {
            let round_poly_evaluations: &Vec<F> = &msgs[round_index];
            if round_poly_evaluations.len() != (vp.max_degree + 1) {
                return Err(ProtocolError::InvalidSumcheck(format!(
                    "incorrect number of evaluations of the {}-th round polynomial",
                    (round_index + 1)
                )));
            }

            let round_poly_evaluation_at_0 = round_poly_evaluations[0];
            let round_poly_evaluation_at_1 = round_poly_evaluations[1];
            let computed_sum = round_poly_evaluation_at_0 + round_poly_evaluation_at_1;

            // Check r_{i}(α_i) == r_{i+1}(0) + r_{i+1}(1)
            if computed_sum != expected_sum {
                return Err(ProtocolError::InvalidSumcheck(format!(
                    "computed sum != expected sum"
                )));
            }

            // Compute r_{i}(α_i) using barycentric interpolation
            expected_sum = barycentric_interpolate(
                &weights,
                &points_vec,
                round_poly_evaluations,
                &challenges[round_index],
            );
        }
        challenges.reverse();
        Ok((expected_sum, evaluations, challenges))
    }
}

#[cfg(test)]
mod test {
    use std::{borrow::Borrow, io::Cursor, iter};

    use crate::{
        poly::multilinear::MultilinearPolynomial,
        sumcheck::{EvalTable, SumCheck, VirtualPolynomial},
        utils::{
            transcript::{InMemoryTranscript, Keccak256Transcript},
            ProtocolError,
        },
    };
    use ff::Field;
    use halo2curves::bn256::Fr;
    use itertools::Itertools;

    use super::{ClassicSumcheck, ClassicSumcheckProverParam, ClassicSumcheckVerifierParam};

    #[test]
    fn test_fold_into_half() {
        let num_vars = 16;
        let evals = (0..1 << num_vars)
            .map(|_| crate::utils::random_fe())
            .collect::<Vec<Fr>>();
        let poly = MultilinearPolynomial::new(evals, vec![], num_vars);

        let mut eval_table = EvalTable::new(num_vars, &poly);
        let evals = eval_table.table().clone();
        let size_before = eval_table.size();

        let alpha: Fr = crate::utils::random_fe();
        eval_table.fold_into_half(alpha);
        let size_after = eval_table.size();
        assert_eq!(2 * size_after, size_before);

        for i in 0..eval_table.size() {
            let expected_even = (Fr::ONE - alpha) * evals[i].even + alpha * evals[i].odd;
            let expected_odd =
                (Fr::ONE - alpha) * evals[size_after + i].even + alpha * evals[size_after + i].odd;

            assert_eq!(eval_table.table()[i].even, expected_even);
            assert_eq!(eval_table.table()[i].odd, expected_odd);
        }
    }

    #[test]
    fn test_sumcheck() -> Result<(), ProtocolError> {
        // Take a simple polynomial
        let num_vars = 3;
        let evals = (0..1 << num_vars)
            .map(|_| crate::utils::random_fe::<Fr>())
            .collect_vec();
        let polys = iter::repeat(MultilinearPolynomial::new(evals, vec![], num_vars))
            .take(3)
            .collect_vec();
        let combine_function = |evals: &Vec<Fr>| evals.iter().product();
        let max_degree = 3;
        let claimed_sum = {
            (0..polys[0].evals().len())
                .map(|idx| {
                    combine_function(
                        polys
                            .iter()
                            .map(|poly| poly.evals()[idx])
                            .collect_vec()
                            .borrow(),
                    )
                })
                .sum()
        };

        // Prover
        let pp = ClassicSumcheckProverParam {
            num_vars,
            max_degree,
        };
        let mut transcript = Keccak256Transcript::<Cursor<Vec<u8>>>::default();
        let virtual_poly = VirtualPolynomial::new(num_vars, polys.iter().collect_vec().borrow());
        ClassicSumcheck::prove(
            &pp,
            &combine_function,
            claimed_sum,
            virtual_poly,
            &mut transcript,
        )?;
        let proof = transcript.into_proof();
        let vp = &ClassicSumcheckVerifierParam {
            num_vars,
            max_degree,
        };
        let mut transcript =
            Keccak256Transcript::<Cursor<Vec<u8>>>::from_proof((), proof.as_slice());
        ClassicSumcheck::verify(vp, max_degree, claimed_sum, polys.len(), &mut transcript)?;
        Ok(())
    }
}
