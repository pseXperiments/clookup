use super::precomputation::Table;
use crate::{
    pcs::PolynomialCommitmentScheme,
    poly::multilinear::MultilinearPolynomial,
    sumcheck::{classic::ClassicSumcheckProverParam, VirtualPolynomial},
    utils::{
        arithmetic::powers, end_timer, start_timer, transcript::TranscriptWrite, transpose,
        ProtocolError,
    },
};
use ff::PrimeField;
use itertools::Itertools;
use rand::RngCore;
use std::{cmp::max, hash::Hash, iter, marker::PhantomData};

#[derive(Clone, Debug)]
struct Prover<
    F: PrimeField + Hash,
    Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
>(PhantomData<F>, PhantomData<Pcs>);

impl<
        F: PrimeField + Hash,
        Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
    > Prover<F, Pcs>
{
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

    fn h_function<'a>(table_poly: &MultilinearPolynomial<F>, gamma: F) -> impl Fn(&'a Vec<F>) -> F + '_ {
        move |evals: &'a Vec<F>| {
            let table_dim = table_poly.num_vars();
            let sigmas = &evals[1..1 + table_dim];
            (evals[0] - table_poly.evaluate(sigmas)
                + sigmas
                    .iter()
                    .zip(powers(gamma).skip(1).take(table_dim))
                    .map(|(sigma, gamma_power)| gamma_power * sigma * (sigma.clone() - F::ONE))
                    .sum::<F>())
                * evals.last().unwrap()
        }
    }

    pub fn prove(
        pp: &Pcs::ProverParam,
        transcript: &mut impl TranscriptWrite<Pcs::CommitmentChunk, F>,
        table: &Table<F>,
        witness: &Vec<F>,
    ) -> Result<(), ProtocolError> {
        let witness_poly =
            MultilinearPolynomial::new(witness.clone(), vec![], witness.len().ilog2() as usize);
        let table_poly = &table.polynomial();
        let num_vars = witness_poly.num_vars();
        let max_degree = 1 + max(2, table_poly.num_vars());
        // get sigma_polys
        let timer = start_timer(|| "sigma_polys");
        let sigma_polys = Self::sigma_polys(table, witness)?;
        end_timer(timer);
        // commit to sigma_polys, witness polys, table polys
        let table_poly_comm = Pcs::commit_and_write(pp, table_poly, transcript)?;
        let witness_poly_comm = Pcs::commit_and_write(pp, &witness_poly, transcript)?;
        let sigma_polys_comms = Pcs::batch_commit_and_write(pp, &sigma_polys, transcript)?;

        // squeeze challenges
        let gamma = transcript.squeeze_challenge();
        let ys = transcript.squeeze_challenges(num_vars);
        let eq = MultilinearPolynomial::eq_xy(&ys);
        let h_function = Self::h_function(table_poly, gamma);
        let virtual_poly = VirtualPolynomial::new(
            num_vars,
            iter::once(&witness_poly)
                .chain([table_poly])
                .chain(sigma_polys.iter())
                .chain(&[eq])
                .collect_vec()
                .as_ref(),
        );
        let pp = ClassicSumcheckProverParam::new(num_vars, max_degree, h_function);
        todo!()
    }
}
