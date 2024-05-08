use std::fmt::Debug;

use ff::Field;
use itertools::Itertools;

use crate::{
    poly::multilinear::MultilinearPolynomial,
    utils::{
        transcript::{FieldTranscriptRead, FieldTranscriptWrite},
        ProtocolError,
    },
};

pub mod classic;
pub mod parallel;

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
        for i in 0..self.table.len() / 2 {
            let eval_pair_distant = self.table[i + self.table.len() / 2].clone();
            let eval_pair = &mut self.table[i];
            eval_pair.even = eval_pair.even + challenge * (eval_pair.odd - eval_pair.even);
            eval_pair.odd = eval_pair_distant.even
                + challenge * (eval_pair_distant.odd - eval_pair_distant.even);
        }
        self.num_vars -= 1;
        self.table.truncate(self.table.len() / 2);
    }
}

pub struct VirtualPolynomial<F: Field> {
    polys: Vec<EvalTable<F>>,
}

impl<F: Field> VirtualPolynomial<F> {
    pub fn new(num_vars: usize, polys: &Vec<MultilinearPolynomial<F>>) -> Self {
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
}

pub trait SumCheck<F: Field>: Clone + Debug {
    type ProverParam: Clone + Debug;
    type VerifierParam: Clone + Debug;

    fn prove(
        pp: &Self::ProverParam,
        sum: F,
        virtual_poly: VirtualPolynomial<F>,
        transcript: &mut impl FieldTranscriptWrite<F>,
    ) -> Result<(F, Vec<F>), ProtocolError>;

    fn verify(
        vp: &Self::VerifierParam,
        degree: usize,
        sum: F,
        transcript: &mut impl FieldTranscriptRead<F>,
    ) -> Result<(F, Vec<F>), ProtocolError>;
}
