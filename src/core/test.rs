mod test {
    use crate::core::{precomputation::Table, prover::Prover, verifier::Verifier};
    use crate::pcs::multilinear::kzg::MultilinearKzg;
    use crate::pcs::PolynomialCommitmentScheme;
    use crate::poly::multilinear::MultilinearPolynomial;
    use crate::utils::{end_timer, start_timer};
    use crate::utils::{
        random_fe,
        transcript::{InMemoryTranscript, Keccak256Transcript},
        ProtocolError,
    };
    use halo2curves::bn256::{Bn256, Fr};
    use itertools::Itertools;
    use std::cmp::max;
    use std::io::Cursor;

    type ClookupProver = Prover<Fr, MultilinearKzg<Bn256>>;
    type ClookupVerifier = Verifier<Fr, MultilinearKzg<Bn256>>;

    #[test]
    pub fn test_clookup() -> Result<(), ProtocolError> {
        let table_dim = 8;
        let witness_dim = 4;
        let table_vec: Vec<Fr> = (0..1 << table_dim).map(|_| random_fe()).collect_vec();
        let witness_vec = table_vec
            .iter()
            .take(1 << witness_dim)
            .cloned()
            .collect_vec();
        let table: Table<Fr> = table_vec.try_into()?;
        let max_degree = 1 + max(2, table_dim);
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
