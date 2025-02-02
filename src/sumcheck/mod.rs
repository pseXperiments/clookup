use std::fmt::Debug;

use ff::{Field, PrimeField};
use itertools::Itertools;

use crate::{
    poly::multilinear::MultilinearPolynomial,
    utils::{
        arithmetic::product,
        ProtocolError,
    },
};
use transcript_utils::transcript::{FieldTranscriptRead, FieldTranscriptWrite};

pub mod classic;
pub mod cuda;
pub mod parallel;

pub fn eq_xy_eval<F: PrimeField>(x: &[F], y: &[F]) -> F {
    assert!(!x.is_empty());
    assert_eq!(x.len(), y.len());

    product(
        x.iter()
            .zip(y)
            .map(|(x_i, y_i)| (*x_i * y_i).double() + F::ONE - x_i - y_i),
    )
}

#[derive(Clone, Debug)]
pub(super) struct EvalPair<F: Field> {
    even: F,
    odd: F,
}

#[derive(Clone, Debug)]
pub(super) struct EvalTable<F: Field> {
    num_vars: usize,
    table: Vec<EvalPair<F>>,
}

impl<F: Field> EvalTable<F> {
    pub fn new(num_vars: usize, poly: &MultilinearPolynomial<F>) -> Self {
        if poly.evals().len() == 0 {
            let table = vec![
                EvalPair {
                    even: F::ZERO,
                    odd: F::ZERO
                };
                1 << num_vars - 1
            ];
            return Self { num_vars, table };
        }
        assert_eq!(poly.evals().len(), 1 << num_vars);
        let table = poly
            .iter()
            .take(1 << (num_vars - 1))
            .copied()
            .zip(poly.iter().skip(1 << (num_vars - 1)).copied())
            .map(|(even, odd)| EvalPair { even, odd })
            .collect_vec();
        Self { num_vars, table }
    }

    pub fn size(&self) -> usize {
        self.table.len()
    }

    pub(super) fn table(&self) -> &Vec<EvalPair<F>> {
        &self.table
    }

    pub fn fold_into_half(&mut self, challenge: F) {
        assert_ne!(self.table.len(), 1);
        let len = self.table.len();
        let (pairs, pairs_distant) = self.table.split_at_mut(len / 2);
        let pairs_distant = &*pairs_distant;
        for i in 0..len / 2 {
            let eval_pair_distant = &pairs_distant[i];
            let eval_pair = &mut pairs[i];
            eval_pair.even = eval_pair.even + challenge * (eval_pair.odd - eval_pair.even);
            eval_pair.odd = eval_pair_distant.even
                + challenge * (eval_pair_distant.odd - eval_pair_distant.even);
        }
        self.num_vars -= 1;
        self.table.truncate(self.table.len() / 2);
    }

    pub fn to_evaluations(&self) -> Vec<F> {
        let mut evals = vec![];
        let mut odd_evals = vec![];
        for pair in &self.table {
            evals.push(pair.even);
            odd_evals.push(pair.odd);
        }
        evals.append(&mut odd_evals);
        evals
    }
}

pub struct VirtualPolynomial<F: Field> {
    polys: Vec<EvalTable<F>>,
}

impl<F: Field> VirtualPolynomial<F> {
    pub fn new(num_vars: usize, polys: &[&MultilinearPolynomial<F>]) -> Self {
        let polys = polys
            .iter()
            .map(|poly| EvalTable::new(num_vars, poly))
            .collect_vec();
        Self { polys }
    }

    pub(super) fn polys(&self) -> &Vec<EvalTable<F>> {
        &self.polys
    }

    pub fn fold_into_half(&mut self, challenge: F) {
        for poly in &mut self.polys {
            poly.fold_into_half(challenge);
        }
    }

    /// called at the last round of sumcheck
    pub fn evaluations(&self, challenge: F) -> Vec<F> {
        self.polys.iter().for_each(|poly| {
            assert_eq!(poly.size(), 1);
        });
        self.polys
            .iter()
            .map(|poly| {
                poly.table()[0].even + challenge * (poly.table()[0].odd - poly.table()[0].even)
            })
            .collect_vec()
    }
}

pub trait SumCheck<F: Field>: Clone + Debug {
    type ProverParam: Clone + Debug;
    type VerifierParam: Clone + Debug;

    fn generate_pp(num_vars: usize, max_degree: usize) -> Result<Self::ProverParam, ProtocolError>;
    fn generate_vp(
        num_vars: usize,
        max_degree: usize,
    ) -> Result<Self::VerifierParam, ProtocolError>;

    /// Returns the challenges and the evaluations of the polynomials at the challenges.
    fn prove(
        pp: &Self::ProverParam,
        combine_function: &impl Fn(&Vec<F>) -> F,
        sum: F,
        virtual_poly: VirtualPolynomial<F>,
        transcript: &mut impl FieldTranscriptWrite<F>,
    ) -> Result<(Vec<F>, Vec<F>), ProtocolError>;

    fn verify(
        vp: &Self::VerifierParam,
        degree: usize,
        sum: F,
        num_polys: usize,
        transcript: &mut impl FieldTranscriptRead<F>,
    ) -> Result<(F, Vec<F>, Vec<F>), ProtocolError>;
}
