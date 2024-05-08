use std::fmt::Debug;

use ff::Field;
use itertools::Itertools;
use rand::RngCore;

use crate::utils::{
    random_fe,
    transcript::{FieldTranscriptRead, FieldTranscriptWrite},
    ProtocolError,
};

pub mod classic;
pub mod parallel;

#[derive(Clone, Debug)]
pub(super) struct EvalPair<F: Field> {
    even: F,
    odd: F,
}

#[derive(Clone, Debug)]
pub struct VirtualPolynomial<F: Field> {
    num_vars: usize,
    evals: Vec<EvalPair<F>>,
}

impl<F: Field> VirtualPolynomial<F> {
    pub fn new(num_vars: usize, evals: Vec<F>) -> Self {
        assert_eq!(evals.len(), 1 << num_vars);
        let evals = evals[..1 << (num_vars - 1)]
            .iter()
            .copied()
            .zip(evals[(1 << (num_vars - 1))..].iter().copied())
            .map(|(even, odd)| EvalPair { even, odd })
            .collect_vec();
        Self { num_vars, evals }
    }

    pub fn size(&self) -> usize {
        self.evals.len()
    }

    pub(super) fn evals(&self) -> &Vec<EvalPair<F>> {
        &self.evals
    }

    pub fn fold_into_half(&mut self, challenge: F) {
        for i in 0..self.evals.len() / 2 {
            let eval_pair_distant = self.evals[i + self.evals.len() / 2].clone();
            let eval_pair = &mut self.evals[i];
            eval_pair.even = eval_pair.even + challenge * (eval_pair.odd - eval_pair.even);
            eval_pair.odd = eval_pair_distant.even
                + challenge * (eval_pair_distant.odd - eval_pair_distant.even);
        }
        self.num_vars -= 1;
        self.evals.truncate(self.evals.len() / 2);
    }

    pub fn get_random_virtual(num_vars: usize) -> Self {
        Self {
            num_vars,
            evals: (0..1 << num_vars)
                .map(|_| EvalPair {
                    even: random_fe(),
                    odd: random_fe(),
                })
                .collect(),
        }
    }
}

pub struct VirtualPolynomial<F: Field> {
    num_vars: usize,
    polys: Vec<EvalTable<F>>,
}

impl<F: Field> VirtualPolynomial<F> {
    pub fn new(num_vars: usize, polys: &Vec<MultilinearPolynomial<F>>) -> Self {
        assert_eq!(polys[0].evals().len(), 1 << num_vars);
        let polys = polys
            .iter()
            .map(|poly| EvalTable::new(num_vars, poly))
            .collect_vec();
        Self { num_vars, polys }
    }

    pub fn polys(&self) -> &Vec<EvalTable<F>> {
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
