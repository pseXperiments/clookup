use super::precomputation::Table;
use crate::{
    poly::multilinear::MultilinearPolynomial,
    utils::{transpose, ProtocolError},
};
use ff::PrimeField;
use std::hash::Hash;

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
struct Prover<F> {
    set: Vec<F>,
    poly: MultilinearPolynomial<F>,
    domain_trans: DomainTransformation<F>,
}

impl<F: PrimeField + Hash> Prover<F> {
    fn setup(set: Vec<F>) -> Self {
        let num_vars = set.len().ilog2() as usize;
        let poly = MultilinearPolynomial::eval_to_coeff(&set, num_vars);
        let domain_trans = DomainTransformation::empty();
        Prover {
            set,
            poly,
            domain_trans,
        }
    }

    fn set_domain_transformation(&mut self, table: Table<F>) -> Result<(), ProtocolError> {
        let indices = table.find_indices(&self.set)?;
        let sigma: Vec<MultilinearPolynomial<F>> = transpose(indices)
            .iter()
            .map(|idx| MultilinearPolynomial::eval_to_coeff(idx, idx.len()))
            .collect();
        self.domain_trans = DomainTransformation::new(sigma);
        Ok(())
    }
}
