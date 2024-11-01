use super::{SumCheck, VirtualPolynomial};
use crate::utils::ProtocolError;
use ff::PrimeField;

#[derive(Clone, Debug)]
pub struct CudaSumcheck;

#[derive(Clone, Debug)]
pub struct CudaSumcheckProverParam {
    num_vars: usize,
    max_degree: usize,
}
#[derive(Clone, Debug)]
pub struct CudaSumcheckVerifierParam {
    num_vars: usize,
    max_degree: usize,
}

impl<F: PrimeField> SumCheck<F> for CudaSumcheck {
    type ProverParam = CudaSumcheckProverParam;
    type VerifierParam = CudaSumcheckVerifierParam;

    fn generate_pp(num_vars: usize, max_degree: usize) -> Result<Self::ProverParam, ProtocolError> {
        todo!()
    }

    fn generate_vp(
        num_vars: usize,
        max_degree: usize,
    ) -> Result<Self::VerifierParam, ProtocolError> {
        todo!()
    }

    fn prove(
        pp: &Self::ProverParam,
        combine_function: &impl Fn(&Vec<F>) -> F,
        sum: F,
        virtual_poly: VirtualPolynomial<F>,
        transcript: &mut impl crate::utils::transcript::FieldTranscriptWrite<F>,
    ) -> Result<(Vec<F>, Vec<F>), ProtocolError> {
        todo!()
    }

    fn verify(
        vp: &Self::VerifierParam,
        degree: usize,
        sum: F,
        num_polys: usize,
        transcript: &mut impl crate::utils::transcript::FieldTranscriptRead<F>,
    ) -> Result<(F, Vec<F>, Vec<F>), ProtocolError> {
        todo!()
    }
}
