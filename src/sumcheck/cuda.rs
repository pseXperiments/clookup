use std::cell::RefCell;

use super::{SumCheck, VirtualPolynomial};
use crate::utils::{
    arithmetic::{barycentric_interpolate, barycentric_weights},
    ProtocolError,
};
use cuda_sumcheck::{
    fieldbinding::{FromFieldBinding, ToFieldBinding},
    GPUApiWrapper,
};
use transcript_utils::transcript::{FieldTranscriptRead, FieldTranscriptWrite};
use cudarc::nvrtc::Ptx;
use ff::PrimeField;
use itertools::Itertools;
use sha3::Keccak256;

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
const MULTILINEAR_PTX: &str = include_str!(concat!(env!("OUT_DIR"), "/multilinear.ptx"));
const SUMCHECK_PTX: &str = include_str!(concat!(env!("OUT_DIR"), "/sumcheck.ptx"));

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

impl<F: PrimeField + FromFieldBinding<F> + ToFieldBinding<F>> SumCheck<F> for CudaSumcheck {
    type ProverParam = CudaSumcheckProverParam;
    type VerifierParam = CudaSumcheckVerifierParam;

    fn generate_pp(num_vars: usize, max_degree: usize) -> Result<Self::ProverParam, ProtocolError> {
        Ok(CudaSumcheckProverParam {
            num_vars,
            max_degree,
        })
    }

    fn generate_vp(
        num_vars: usize,
        max_degree: usize,
    ) -> Result<Self::VerifierParam, ProtocolError> {
        Ok(CudaSumcheckVerifierParam {
            num_vars,
            max_degree,
        })
    }

    fn prove(
        pp: &Self::ProverParam,
        combine_function: &impl Fn(&Vec<F>) -> F,
        sum: F,
        virtual_poly: VirtualPolynomial<F>,
        transcript: &mut impl FieldTranscriptWrite<F>,
    ) -> Result<(Vec<F>, Vec<F>), ProtocolError> {
        let mut challenges = vec![];
        let mut gpu_api_wrapper =
            GPUApiWrapper::<F>::setup().map_err(|e| ProtocolError::CudaLibraryError(e.to_string()))?;
        gpu_api_wrapper
            .gpu
            .load_ptx(
                Ptx::from_src(MULTILINEAR_PTX),
                "multilinear",
                &["convert_to_montgomery_form"],
            )
            .map_err(|e| ProtocolError::CudaLibraryError(e.to_string()))?;

        gpu_api_wrapper
            .gpu
            .load_ptx(
                Ptx::from_src(SUMCHECK_PTX),
                "sumcheck",
                &[
                    "fold_into_half",
                    "fold_into_half_in_place",
                    "combine",
                    "sum",
                ],
            )
            .map_err(|e| ProtocolError::CudaLibraryError(e.to_string()))?;
        let polys: Vec<Vec<F>> = virtual_poly
            .polys()
            .iter()
            .map(|table| table.to_evaluations())
            .collect();

        let mut gpu_polys = gpu_api_wrapper
            .copy_to_device(&polys.concat())
            .map_err(|e| ProtocolError::CudaLibraryError(e.to_string()))?;
        let device_ks = (0..pp.max_degree + 1)
            .map(|k| {
                gpu_api_wrapper
                    .gpu
                    .htod_copy(vec![F::to_montgomery_form(F::from(k as u64))])
            })
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| ProtocolError::CudaLibraryError(e.to_string()))?;

        let mut buf = gpu_api_wrapper
            .malloc_on_device(polys.len() << (pp.num_vars - 1))
            .map_err(|e| ProtocolError::CudaLibraryError(e.to_string()))?;
        let buf_view = RefCell::new(buf.slice_mut(..));

        let mut challenges_cuda = gpu_api_wrapper
            .malloc_on_device(1)
            .map_err(|e| ProtocolError::CudaLibraryError(e.to_string()))?;

        let mut round_evals = gpu_api_wrapper
            .malloc_on_device(pp.max_degree + 1)
            .map_err(|e| ProtocolError::CudaLibraryError(e.to_string()))?;
        let round_evals_view = RefCell::new(round_evals.slice_mut(..));

        // This is shit code
        let gamma = gpu_api_wrapper
            .gpu
            .htod_copy(vec![F::to_montgomery_form(combine_function(&vec![]))])
            .map_err(|e| ProtocolError::CudaLibraryError(e.to_string()))?;

        gpu_api_wrapper
            .prove_sumcheck(
                pp.num_vars,
                polys.len(),
                pp.max_degree,
                sum,
                &mut gpu_polys.slice_mut(..),
                &device_ks
                    .iter()
                    .map(|device_k| device_k.slice(..))
                    .collect_vec(),
                buf_view,
                &mut challenges_cuda,
                round_evals_view,
                transcript,
                &gamma.slice(..),
                &mut challenges,
            )
            .map_err(|e| ProtocolError::CudaLibraryError(String::from("library error")))?;

        gpu_api_wrapper
            .gpu
            .synchronize()
            .map_err(|e| ProtocolError::CudaLibraryError(e.to_string()))?;

        let evaluations = (0..polys.len())
            .map(|i| {
                gpu_api_wrapper
                    .dtoh_sync_copy(
                        &gpu_polys.slice(i << pp.num_vars..(i * 2 + 1) << (pp.num_vars - 1)),
                        true,
                    )
                    .map(|res| res.first().unwrap().clone())
                    .map_err(|e| ProtocolError::CudaLibraryError(e.to_string()))
            })
            .collect::<Result<Vec<F>, _>>()
            .map_err(|e| ProtocolError::CudaLibraryError(String::from("")))?;

        challenges.reverse();

        println!("{:?}", challenges);

        Ok((challenges, evaluations))
    }

    fn verify(
        vp: &Self::VerifierParam,
        _: usize,
        sum: F,
        num_polys: usize,
        transcript: &mut impl FieldTranscriptRead<F>,
    ) -> Result<(F, Vec<F>, Vec<F>), ProtocolError> {
        let (msgs, mut challenges) = {
            let mut msgs = Vec::with_capacity(vp.num_vars);
            let mut challenges = Vec::with_capacity(vp.num_vars);
            for _ in 0..vp.num_vars {
                msgs.push(transcript.read_field_elements(vp.max_degree + 1).map_err(|_| ProtocolError::Transcript)?);
                challenges.push(transcript.squeeze_challenge());
            }
            (msgs, challenges)
        };

        

        let evaluations = transcript.read_field_elements(num_polys).map_err(|_| ProtocolError::Transcript)?;
        let mut expected_sum = sum.clone();
        let points_vec: Vec<F> = (0..vp.max_degree + 1)
            .map(|i| F::from_u128(i as u128))
            .collect();
        let weights = barycentric_weights(&points_vec);

        for round_index in 0..vp.num_vars {
            let round_poly_evaluations: &Vec<F> = &msgs[round_index];
            if round_poly_evaluations.len() != (vp.max_degree + 1) {
                return Err(ProtocolError::InvalidSumcheck(format!(
                    "incorrect number of evaluations of the {}-th round polynomial",
                    (round_index + 1)
                )));
            }

            let round_poly_evaluation_at_0 = round_poly_evaluations[0];
            let round_poly_evaluation_at_1 = round_poly_evaluations[1];
            let computed_sum = round_poly_evaluation_at_0 + round_poly_evaluation_at_1;

            // Check r_{i}(α_i) == r_{i+1}(0) + r_{i+1}(1)
            if computed_sum != expected_sum {
                println!("round : {:?}", round_index);
                return Err(ProtocolError::InvalidSumcheck(format!(
                    "computed sum != expected sum"
                )));
            }

            // Compute r_{i}(α_i) using barycentric interpolation
            expected_sum = barycentric_interpolate(
                &weights,
                &points_vec,
                round_poly_evaluations,
                &challenges[round_index],
            );
        }
        challenges.reverse();
        Ok((expected_sum, evaluations, challenges))
    }
}
