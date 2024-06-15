use std::{cmp::max, io::Cursor};

use clookup::{
    core::{precomputation::Table, prover::Prover, verifier::Verifier},
    pcs::{multilinear::kzg::MultilinearKzg, PolynomialCommitmentScheme},
    sumcheck::{classic::ClassicSumcheck, parallel::ParallelSumcheck},
    utils::{
        end_timer, random_fe, start_timer,
        transcript::{InMemoryTranscript, Keccak256Transcript},
        ProtocolError,
    },
};
use criterion::{criterion_group, criterion_main, Criterion};
use halo2curves::bn256::{Bn256, Fr};
use itertools::Itertools;

type ClookupProverClassic = Prover<Fr, MultilinearKzg<Bn256>, ClassicSumcheck>;
type ClookupVerifierClassic = Verifier<Fr, MultilinearKzg<Bn256>, ClassicSumcheck>;

type ClookupProverPar = Prover<Fr, MultilinearKzg<Bn256>, ParallelSumcheck>;
type ClookupVerifierPar = Verifier<Fr, MultilinearKzg<Bn256>, ParallelSumcheck>;

fn set_env(table_dim: usize, witness_dim: usize) -> (Table<Fr>, Vec<Fr>) {
    let table_vec: Vec<Fr> = (0..1 << table_dim).map(|_| random_fe()).collect_vec();
    let witness_vec = table_vec
        .iter()
        .take(1 << witness_dim)
        .cloned()
        .collect_vec();
    let table: Table<Fr> = table_vec.try_into().unwrap();
    (table, witness_vec)
}

pub fn test_clookup_classic() -> Result<(), ProtocolError> {
    let table_dim = 16;
    let witness_dim = 8;

    let (table, witness_vec) = set_env(table_dim, witness_dim);

    let max_degree = 1 + max(2, table_dim);
    let (pp, vp) = {
        let rng = rand::thread_rng();
        let param = ClookupProverClassic::setup(&table, &witness_vec, rng)?;
        MultilinearKzg::trim(&param, 1 << witness_dim, 1).unwrap()
    };
    let timer = start_timer(|| "clookup prover");
    let proof = {
        let mut transcript = Keccak256Transcript::<Cursor<Vec<u8>>>::default();
        ClookupProverClassic::prove(&pp, &mut transcript, &table, &witness_vec)?;
        transcript.into_proof()
    };
    end_timer(timer);
    let mut transcript = Keccak256Transcript::<Cursor<Vec<u8>>>::from_proof((), proof.as_slice());
    ClookupVerifierClassic::verify(
        &vp,
        &mut transcript,
        table_dim + 2,
        table_dim,
        witness_dim,
        max_degree,
    )?;
    Ok(())
}

pub fn test_clookup_par() -> Result<(), ProtocolError> {
    let table_dim = 16;
    let witness_dim = 8;

    let (table, witness_vec) = set_env(table_dim, witness_dim);

    let max_degree = 1 + max(2, table_dim);
    let (pp, vp) = {
        let rng = rand::thread_rng();
        let param = ClookupProverPar::setup(&table, &witness_vec, rng)?;
        MultilinearKzg::trim(&param, 1 << witness_dim, 1).unwrap()
    };
    let timer = start_timer(|| "clookup prover");
    let proof = {
        let mut transcript = Keccak256Transcript::<Cursor<Vec<u8>>>::default();
        ClookupProverPar::prove(&pp, &mut transcript, &table, &witness_vec)?;
        transcript.into_proof()
    };
    end_timer(timer);
    let mut transcript = Keccak256Transcript::<Cursor<Vec<u8>>>::from_proof((), proof.as_slice());
    ClookupVerifierPar::verify(
        &vp,
        &mut transcript,
        table_dim + 2,
        table_dim,
        witness_dim,
        max_degree,
    )?;
    Ok(())
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("e2e-clookup", |b| b.iter(|| test_clookup_classic()));
}
criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
