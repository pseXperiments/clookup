use crate::{
    poly::Polynomial,
    utils::{
        arithmetic::Field,
        transcript::{TranscriptRead, TranscriptWrite},
        DeserializeOwned, ProtocolError, Serialize,
    },
};
use rand::RngCore;
use std::fmt::Debug;

pub mod multilinear;

pub type Point<F, P> = <P as Polynomial<F>>::Point;

pub type Commitment<F, Pcs> = <Pcs as PolynomialCommitmentScheme<F>>::Commitment;

pub type CommitmentChunk<F, Pcs> = <Pcs as PolynomialCommitmentScheme<F>>::CommitmentChunk;

pub trait PolynomialCommitmentScheme<F: Field>: Clone + Debug {
    type Param: Clone + Debug + Serialize + DeserializeOwned;
    type ProverParam: Clone + Debug + Serialize + DeserializeOwned;
    type VerifierParam: Clone + Debug + Serialize + DeserializeOwned;
    type Polynomial: Polynomial<F> + Serialize + DeserializeOwned;
    type Commitment: Clone
        + Debug
        + Default
        + AsRef<[Self::CommitmentChunk]>
        + Serialize
        + DeserializeOwned;
    type CommitmentChunk: Clone + Debug + Default;

    fn setup(
        poly_size: usize,
        batch_size: usize,
        rng: impl RngCore,
    ) -> Result<Self::Param, ProtocolError>;

    fn trim(
        param: &Self::Param,
        poly_size: usize,
        batch_size: usize,
    ) -> Result<(Self::ProverParam, Self::VerifierParam), ProtocolError>;

    fn commit(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
    ) -> Result<Self::Commitment, ProtocolError>;

    fn commit_and_write(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
        transcript: &mut impl TranscriptWrite<Self::CommitmentChunk, F>,
    ) -> Result<Self::Commitment, ProtocolError> {
        let comm = Self::commit(pp, poly)?;
        transcript.write_commitments(comm.as_ref())?;
        Ok(comm)
    }

    fn batch_commit<'a>(
        pp: &Self::ProverParam,
        polys: impl IntoIterator<Item = &'a Self::Polynomial>,
    ) -> Result<Vec<Self::Commitment>, ProtocolError>
    where
        Self::Polynomial: 'a;

    fn batch_commit_and_write<'a>(
        pp: &Self::ProverParam,
        polys: impl IntoIterator<Item = &'a Self::Polynomial>,
        transcript: &mut impl TranscriptWrite<Self::CommitmentChunk, F>,
    ) -> Result<Vec<Self::Commitment>, ProtocolError>
    where
        Self::Polynomial: 'a,
    {
        let comms = Self::batch_commit(pp, polys)?;
        for comm in comms.iter() {
            transcript.write_commitments(comm.as_ref())?;
        }
        Ok(comms)
    }

    fn open(
        pp: &Self::ProverParam,
        poly: &Self::Polynomial,
        comm: &Self::Commitment,
        point: &Point<F, Self::Polynomial>,
        eval: &F,
        transcript: &mut impl TranscriptWrite<Self::CommitmentChunk, F>,
    ) -> Result<(), ProtocolError>;

    fn batch_open<'a>(
        pp: &Self::ProverParam,
        polys: impl IntoIterator<Item = &'a Self::Polynomial>,
        comms: impl IntoIterator<Item = &'a Self::Commitment>,
        points: &[Point<F, Self::Polynomial>],
        evals: &[Evaluation<F>],
        transcript: &mut impl TranscriptWrite<Self::CommitmentChunk, F>,
    ) -> Result<(), ProtocolError>
    where
        Self::Polynomial: 'a,
        Self::Commitment: 'a;

    fn read_commitment(
        vp: &Self::VerifierParam,
        transcript: &mut impl TranscriptRead<Self::CommitmentChunk, F>,
    ) -> Result<Self::Commitment, ProtocolError> {
        let comms = Self::read_commitments(vp, 1, transcript)?;
        assert_eq!(comms.len(), 1);
        Ok(comms.into_iter().next().unwrap())
    }

    fn read_commitments(
        vp: &Self::VerifierParam,
        num_polys: usize,
        transcript: &mut impl TranscriptRead<Self::CommitmentChunk, F>,
    ) -> Result<Vec<Self::Commitment>, ProtocolError>;

    fn verify(
        vp: &Self::VerifierParam,
        comm: &Self::Commitment,
        point: &Point<F, Self::Polynomial>,
        eval: &F,
        transcript: &mut impl TranscriptRead<Self::CommitmentChunk, F>,
    ) -> Result<(), ProtocolError>;

    fn batch_verify<'a>(
        vp: &Self::VerifierParam,
        comms: impl IntoIterator<Item = &'a Self::Commitment>,
        points: &[Point<F, Self::Polynomial>],
        evals: &[Evaluation<F>],
        transcript: &mut impl TranscriptRead<Self::CommitmentChunk, F>,
    ) -> Result<(), ProtocolError>
    where
        Self::Commitment: 'a;
}

#[derive(Clone, Debug)]
pub struct Evaluation<F> {
    poly: usize,
    point: usize,
    value: F,
}

impl<F> Evaluation<F> {
    pub fn new(poly: usize, point: usize, value: F) -> Self {
        Self { poly, point, value }
    }

    pub fn poly(&self) -> usize {
        self.poly
    }

    pub fn point(&self) -> usize {
        self.point
    }

    pub fn value(&self) -> &F {
        &self.value
    }
}

pub trait Additive<F: Field>: Clone + Debug + Default + PartialEq + Eq {
    fn msm<'a, 'b>(
        scalars: impl IntoIterator<Item = &'a F>,
        bases: impl IntoIterator<Item = &'b Self>,
    ) -> Self
    where
        Self: 'b;
}
