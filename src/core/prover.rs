use super::precomputation::Table;
use crate::{
    pcs::PolynomialCommitmentScheme, poly::multilinear::MultilinearPolynomial, utils::{end_timer, start_timer, transcript::TranscriptWrite, transpose, ProtocolError}
};
use ff::PrimeField;
use rand::RngCore;
use std::{cmp::max, hash::Hash, marker::PhantomData};

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

    fn sigma_polys(
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
        pp: &Pcs::ProverParam,
        transcript: &mut impl TranscriptWrite<Pcs::CommitmentChunk, F>,
        table: &Table<F>,
        witness: &Vec<F>,
    ) -> Result<(), ProtocolError> {
        let witness_poly = MultilinearPolynomial::new(witness.clone(), vec![], witness.len().ilog2() as usize);
        let num_vars = witness_poly.num_vars();
        // get sigma_polys
        let timer = start_timer(|| "sigma_polys");
        let sigma_polys = Self::sigma_polys(table, witness)?;
        end_timer(timer);
        // commit to sigma_polys, witness polys, table polys
        let table_poly_comm = Pcs::commit_and_write(pp, &table.polynomial(), transcript)?;
        let witness_poly_comm = Pcs::commit_and_write(pp, &witness_poly, transcript)?;
        let sigma_polys_comms = Pcs::batch_commit_and_write(pp, &sigma_polys, transcript)?;

        // squeeze challenges
        let gamma = transcript.squeeze_challenge();
        let r = transcript.squeeze_challenges(num_vars);
        todo!()
    }
}
