use super::{SumCheck, VirtualPolynomial};
use crate::utils::transcript::{FieldTranscriptRead, FieldTranscriptWrite};
use crate::utils::ProtocolError;
use ff::{Field, PrimeField};

#[derive(Clone, Debug)]
struct ClassicSumcheck;

#[derive(Clone, Debug)]
struct ClassicSumcheckProverParam<F: Field> {
    num_vars: usize,
    max_degree: usize,
    combine_function: fn(&Vec<F>) -> F,
}

#[derive(Clone, Debug)]
struct ClassicSumcheckVerifierParam<F: Field> {
    combine_function: fn(Vec<F>) -> F,
}

impl<F: PrimeField> SumCheck<F> for ClassicSumcheck {
    type ProverParam = ClassicSumcheckProverParam<F>;
    type VerifierParam = ClassicSumcheckVerifierParam<F>;

    fn prove(
        pp: &Self::ProverParam,
        sum: F,
        mut virtual_poly: VirtualPolynomial<F>,
        transcript: &mut impl FieldTranscriptWrite<F>,
    ) -> Result<(Vec<F>, Vec<F>), ProtocolError> {
        // Declare r_polys and initialise it with 0s
        let r_degree = pp.max_degree;
        let mut r_polys: Vec<Vec<F>> = (0..pp.num_vars)
            .map(|_| vec![F::ZERO; r_degree + 1])
            .collect();

        transcript.write_field_element(&sum)?;

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
                    r_polys[round_index][k] += (pp.combine_function)(&evaluations_at_k);
                }
            }
            // append the round polynomial (i.e. prover message) to the transcript
            transcript.write_field_elements(&r_polys[round_index])?;

            // generate challenge α_i = H( transcript );
            let alpha = transcript.squeeze_challenge();

            // update prover state polynomials
            virtual_poly.fold_into_half(alpha);
        }

        for i in 0..virtual_poly.polys().len() {
            assert_eq!(virtual_poly.polys()[i].size(), 1);
        }

        Ok((vec![], vec![]))
    }

    fn verify(
        vp: &Self::VerifierParam,
        degree: usize,
        sum: F,
        transcript: &mut impl FieldTranscriptRead<F>,
    ) -> Result<(Vec<F>, Vec<F>), ProtocolError> {
        Ok((vec![], vec![]))
    }
}

#[cfg(test)]
mod test {
    use std::{borrow::Borrow, io::Cursor, iter};

    use crate::{
        poly::multilinear::MultilinearPolynomial,
        sumcheck::{EvalTable, SumCheck, VirtualPolynomial},
        utils::{transcript::Keccak256Transcript, ProtocolError},
    };
    use ff::Field;
    use halo2curves::bn256::Fr;
    use itertools::Itertools;

    use super::{ClassicSumcheck, ClassicSumcheckProverParam};

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
            max_degree: 3,
            combine_function,
        };
        let mut transcript = Keccak256Transcript::<Cursor<Vec<u8>>>::default();
        let virtual_poly = VirtualPolynomial::new(num_vars, &polys);
        ClassicSumcheck::prove(&pp, claimed_sum, virtual_poly, &mut transcript)?;

        Ok(())
    }
}
