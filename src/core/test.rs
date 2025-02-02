mod test {
    use crate::core::cuda_prover::CudaProver;
    use crate::core::{precomputation::Table, verifier::Verifier};
    use crate::pcs::multilinear::kzg::MultilinearKzg;
    use crate::pcs::PolynomialCommitmentScheme;
    use crate::poly::multilinear::MultilinearPolynomial;
    use crate::sumcheck::classic::ClassicSumcheck;
    use crate::sumcheck::cuda::CudaSumcheck;
    use crate::utils::{end_timer, start_timer};
    use crate::utils::{
        random_fe,
        ProtocolError,
    };
    use halo2curves::bn256::{Bn256, Fr};
    use itertools::Itertools;
    use std::cmp::max;
    use std::io::Cursor;
    use transcript_utils::transcript::{InMemoryTranscript, Keccak256Transcript};

    type ClookupProver = CudaProver<Fr, MultilinearKzg<Bn256>>;
    type ClookupVerifier = Verifier<Fr, MultilinearKzg<Bn256>, CudaSumcheck>;

    #[test]
    pub fn test_clookup() -> Result<(), ProtocolError> {
        let table_dim = 8;
        let witness_dim = 4;
        // Range table 0..1 << table_dim - 1
        let table_vec: Vec<Fr> = (0..1 << table_dim).map(|i| Fr::from(i)).collect_vec();
        let witness_vec = table_vec
            .iter()
            .take(1 << witness_dim)
            .cloned()
            .collect_vec();
        let table: Table<Fr> = table_vec.try_into()?;
        let max_degree = 3;
        let (pp, vp) = {
            let rng = rand::thread_rng();
            let param = ClookupProver::setup(&table, &witness_vec, rng)?;
            MultilinearKzg::trim(&param, 1 << witness_dim, 1).unwrap()
        };
        let timer = start_timer(|| "clookup prover");
        let proof = {
            let mut transcript = Keccak256Transcript::<Cursor<Vec<u8>>>::default();
            ClookupProver::prove(&pp, &mut transcript, &table, &witness_vec)?;
            transcript.into_proof()
        };
        end_timer(timer);
        let mut transcript =
            Keccak256Transcript::<Cursor<Vec<u8>>>::from_proof((), proof.as_slice());
        ClookupVerifier::verify(
            &vp,
            &mut transcript,
            table_dim + 2,
            table_dim,
            witness_dim,
            max_degree,
        )?;
        Ok(())
    }
}
