use std::{hash::Hash, iter, marker::PhantomData};

use ff::PrimeField;
use itertools::Itertools;

use crate::{
    pcs::{Evaluation, PolynomialCommitmentScheme},
    poly::multilinear::MultilinearPolynomial,
    sumcheck::SumCheck,
    utils::ProtocolError,
};
use transcript_utils::transcript::TranscriptRead;

#[derive(Clone, Debug)]
pub struct Verifier<
    F: PrimeField + Hash,
    Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
    Scs: SumCheck<F>,
>(PhantomData<F>, PhantomData<Pcs>, PhantomData<Scs>);

impl<
        F: PrimeField + Hash,
        Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
        Scs: SumCheck<F>,
    > Verifier<F, Pcs, Scs>
{
    pub fn verify(
        vp: &Pcs::VerifierParam,
        transcript: &mut impl TranscriptRead<Pcs::CommitmentChunk, F>,
        num_polys: usize,
        table_dimension: usize,
        witness_num_vars: usize,
        max_degree: usize,
    ) -> Result<(), ProtocolError> {
        let witness_comm = Pcs::read_commitment(vp, transcript)?;
        let sigma_comm = Pcs::read_commitments(vp, table_dimension, transcript)?;

        let gamma = transcript.squeeze_challenge();
        let ys = transcript.squeeze_challenges(witness_num_vars);

        let svp = Scs::generate_vp(witness_num_vars, max_degree)?;
        let (_, evals, x) = Scs::verify(&svp, max_degree, F::ZERO, num_polys, transcript)?;
        let witness_poly_x = evals.first().unwrap();

        let sigma_polys_x = evals.iter().skip(1).take(table_dimension).collect_vec();

        let comms = iter::once(&witness_comm).chain(sigma_comm.iter());
        let points_vec = iter::repeat(x).take(1 + table_dimension).collect_vec();
        let points = points_vec.as_slice();
        let evals_vec = iter::once(witness_poly_x)
            .chain(sigma_polys_x)
            .enumerate()
            .map(|(poly, value)| Evaluation::new(poly, 0, *value))
            .collect_vec();
        let evals = evals_vec.as_slice();

        Pcs::batch_verify(vp, comms, points, evals, transcript)?;
        Ok(())
    }
}
