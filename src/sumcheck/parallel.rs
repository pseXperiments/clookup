use ff::PrimeField;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::utils::{
    arithmetic::{barycentric_interpolate, barycentric_weights},
    transcript::{FieldTranscriptRead, FieldTranscriptWrite},
    ProtocolError,
};

use super::{SumCheck, VirtualPolynomial};

#[derive(Clone, Debug)]
pub struct ParallelSumcheck;

#[derive(Clone, Debug)]
pub struct ParallelSumcheckProverParam {
    num_vars: usize,
    max_degree: usize,
}

#[derive(Clone, Debug)]
pub struct ParallelSumcheckVerifierParam {
    num_vars: usize,
    max_degree: usize,
}

impl<F: PrimeField> SumCheck<F> for ParallelSumcheck {
    type ProverParam = ParallelSumcheckProverParam;
    type VerifierParam = ParallelSumcheckVerifierParam;

    fn generate_pp(num_vars: usize, max_degree: usize) -> Result<Self::ProverParam, ProtocolError> {
        Ok(ParallelSumcheckProverParam {
            num_vars,
            max_degree,
        })
    }

    fn generate_vp(
        num_vars: usize,
        max_degree: usize,
    ) -> Result<Self::VerifierParam, ProtocolError> {
        Ok(ParallelSumcheckVerifierParam {
            num_vars,
            max_degree,
        })
    }

    fn prove(
        pp: &Self::ProverParam,
        combine_function: &impl Fn(&Vec<F>) -> F,
        sum: F,
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
                        .par_iter()
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
        degree: usize,
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
