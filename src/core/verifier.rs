use std::{hash::Hash, marker::PhantomData, iter};

use ff::PrimeField;
use itertools::Itertools;

use crate::{
    pcs::{PolynomialCommitmentScheme, Evaluation},
    poly::multilinear::MultilinearPolynomial,
    utils::{transcript::TranscriptRead, ProtocolError}, sumcheck::{classic::{ClassicSumcheck, ClassicSumcheckVerifierParam}, SumCheck},
};

#[derive(Clone, Debug)]
pub struct Verifier<
    F: PrimeField + Hash,
    Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
>(PhantomData<F>, PhantomData<Pcs>);

impl<
        F: PrimeField + Hash,
        Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
    > Verifier<F, Pcs>
{
    pub fn verify(
        vp: &Pcs::VerifierParam,
        transcript: &mut impl TranscriptRead<Pcs::CommitmentChunk, F>,
        num_polys: usize,
        table_dimension: usize,
        witness_num_vars: usize,
        max_degree: usize,
    ) -> Result<(), ProtocolError> {
        let table_comm = Pcs::read_commitment(vp, transcript)?;
        let witness_comm = Pcs::read_commitment(vp, transcript)?;
        let sigma_comm = Pcs::read_commitments(vp, table_dimension, transcript)?;
        
        let gamma = transcript.squeeze_challenge();
        let ys = transcript.squeeze_challenges(witness_num_vars);

        let svp = ClassicSumcheckVerifierParam::new(witness_num_vars, max_degree);
        let (evals, x) = ClassicSumcheck::verify(&svp, max_degree, F::ZERO, num_polys, transcript)?;
        let witness_poly_x = evals.first().unwrap();
        let table_poly_x = evals.get(1).unwrap();
        let sigma_poly_x = evals
            .iter()
            .skip(2)
            .take(table_dimension)
            .collect_vec();
        Pcs::verify(vp, &table_comm, &x, table_poly_x, transcript)?;
        
        let comms = iter::once(&witness_comm).chain(sigma_comm.iter());
        let points_vec = iter::repeat(x).take(1 + table_dimension).collect_vec();
        let points = points_vec.as_slice();
        let evals_vec = iter::once(witness_poly_x)
            .chain(sigma_poly_x)
            .enumerate()
            .map(|(poly, value)| Evaluation::new(poly, 0, *value))
            .collect_vec();
        let evals = evals_vec.as_slice();
        Pcs::batch_verify(vp, comms, points, evals, transcript)?;
        Ok(())
    }
}
