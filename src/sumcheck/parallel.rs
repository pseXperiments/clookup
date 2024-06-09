use ff::PrimeField;

use crate::utils::{
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
        degree: usize,
        sum: F,
        num_polys: usize,
        transcript: &mut impl FieldTranscriptRead<F>,
    ) -> Result<(F, Vec<F>, Vec<F>), ProtocolError> {
        Ok((F::ZERO, vec![], vec![]))
    }
}
