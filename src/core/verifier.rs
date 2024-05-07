use std::{hash::Hash, marker::PhantomData};

use ff::PrimeField;

use crate::{pcs::PolynomialCommitmentScheme, poly::multilinear::MultilinearPolynomial};

struct Verifier<
    F: PrimeField + Hash,
    Pcs: PolynomialCommitmentScheme<F, Polynomial = MultilinearPolynomial<F>>,
>(PhantomData<F>, PhantomData<Pcs>);
