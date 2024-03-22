use std::fmt::Debug;
use ff::Field;

// These are different form of multivariate polynomial
pub mod coefficient;
pub mod factored;
pub mod composition;

pub trait MultivariatePolynomial<F: Field>: Clone + Debug  {
    fn interpolate(points: Vec<F>) -> Self;
}
