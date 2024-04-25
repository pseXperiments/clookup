use std::fmt::Debug;

use ff::Field;

use crate::utils::{transcript::{FieldTranscriptRead, FieldTranscriptWrite}, ProtocolError};

pub mod classic;
pub mod parallel;

pub trait SumCheck<F: Field>: Clone + Debug {
    type ProverParam: Clone + Debug;
    type VerifierParam: Clone + Debug;

    fn prove(
        pp: &Self::ProverParam,
        num_vars: usize,
        sum: F,
        transcript: &mut impl FieldTranscriptWrite<F>,
    ) -> Result<(F, Vec<F>), ProtocolError>;

    fn verify(
        vp: &Self::VerifierParam,
        num_vars: usize,
        degree: usize,
        sum: F,
        transcript: &mut impl FieldTranscriptRead<F>,
    ) -> Result<(F, Vec<F>), ProtocolError>;
}
