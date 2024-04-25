use super::precomputation::Table;
use crate::{
    pcs::PolynomialCommitmentScheme, poly::multilinear::MultilinearPolynomial, utils::{transcript::TranscriptWrite, transpose, ProtocolError}
};
use ff::PrimeField;
use rand::RngCore;
use std::{cmp::max, hash::Hash, marker::PhantomData};

#[derive(Clone, Debug)]
struct DomainTransformation<F> {
    sigma: Vec<MultilinearPolynomial<F>>,
}

impl<F> DomainTransformation<F> {
    fn empty() -> Self {
        DomainTransformation { sigma: vec![] }
    }

    fn new(sigma: Vec<MultilinearPolynomial<F>>) -> Self {
        Self { sigma }
    }
}

#[derive(Clone, Debug)]
struct Prover<F: PrimeField + Hash, Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,> (PhantomData<F>, PhantomData<Pcs>);

impl<F: PrimeField + Hash, Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>> Prover<F, Pcs> {
    fn setup(
        table: &Table<F>,
        witness: &Vec<F>,
        rng: impl RngCore,
    ) -> Result<Pcs::Param, ProtocolError> {
        let poly_size = max(table.len(), witness.len());
        let batch_size = 1 + 1 + table.num_vars();
        Pcs::setup(poly_size, batch_size, rng)
    }

    fn set_domain_transformation(
        table: &Table<F>,
        witness: &Vec<F>,
    ) -> Result<Vec<MultilinearPolynomial<F>>, ProtocolError> {
        let indices = table.find_indices(&witness)?;
        let sigma: Vec<MultilinearPolynomial<F>> = transpose(indices)
            .iter()
            .map(|idx| MultilinearPolynomial::eval_to_coeff(idx, idx.len()))
            .collect();
        Ok(sigma)
    }

    pub fn prove(
        transcript: &mut impl TranscriptWrite<Pcs::CommitmentChunk, F>,
        table: &Table<F>,
        witness: &Vec<F>,
    ) -> Result<(), ProtocolError> {
        todo!()
    }
}
