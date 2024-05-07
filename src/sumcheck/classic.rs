use super::{SumCheck, VirtualPolynomial};
use crate::poly::multilinear::MultilinearPolynomial;
use crate::utils::transcript::{FieldTranscriptRead, FieldTranscriptWrite};
use crate::utils::ProtocolError;
use ff::{Field, PrimeField};

#[derive(Clone, Debug)]
struct ClassicSumcheck;

#[derive(Clone, Debug)]
struct ClassicSumcheckProverParam<F: Field> {
    num_vars: usize,
    round: usize,
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
        num_vars: usize,
        sum: F,
        mut virtual_polys: Vec<VirtualPolynomial<F>>,
        transcript: &mut impl FieldTranscriptWrite<F>,
    ) -> Result<(F, Vec<F>), ProtocolError> {
        // Declare r_polys and initialise it with 0s
        let r_degree = pp.max_degree;
        let mut r_polys: Vec<Vec<F>> = (0..pp.num_vars)
            .map(|_| vec![F::ZERO; r_degree + 1])
            .collect();

        transcript.write_field_element(&sum)?;

        for round_index in 0..pp.num_vars {
            let virtual_polynomial_len = virtual_polys[0].evals().len();
            for k in 0..(r_degree + 1) {
                for i in 0..virtual_polynomial_len {
                    let evaluations_at_k = virtual_polys
                        .iter()
                        .map(|virtual_poly| {
                            let o = virtual_poly.evals()[i].odd;
                            let e = virtual_poly.evals()[i].even;
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
            for j in 0..virtual_polys.len() {
                virtual_polys[j].fold_into_half(alpha);
            }
        }
        Ok((F::ONE, vec![]))
    }

    fn verify(
        vp: &Self::VerifierParam,
        num_vars: usize,
        degree: usize,
        sum: F,
        transcript: &mut impl FieldTranscriptRead<F>,
    ) -> Result<(F, Vec<F>), ProtocolError> {
        Ok((F::ONE, vec![]))
    }
}

#[cfg(test)]
mod test {
    use crate::sumcheck::classic::*;
    use ff::Field;
    use halo2curves::bn256::Fr;

    #[test]
    fn test_fold_into_half() {
        let num_vars = 3;
        let evals_list = (0..1 << num_vars)
            .map(|_| crate::utils::random_fe())
            .collect::<Vec<Fr>>();

        let mut virtual_poly: VirtualPolynomial<Fr> = VirtualPolynomial::new(num_vars, evals_list);
        let evals = virtual_poly.evals().clone();
        let size_before = virtual_poly.size();

        let alpha: Fr = crate::utils::random_fe();
        virtual_poly.fold_into_half(alpha);
        let size_after = virtual_poly.size();
        assert_eq!(2 * size_after, size_before);

        for i in 0..virtual_poly.size() {
            let expected_even = (Fr::ONE - alpha) * evals[i].even + alpha * evals[i].odd;
            let expected_odd =
                (Fr::ONE - alpha) * evals[size_after + i].even + alpha * evals[size_after + i].odd;

            assert_eq!(virtual_poly.evals()[i].even, expected_even);
            assert_eq!(virtual_poly.evals()[i].odd, expected_odd);
        }
    }
}
