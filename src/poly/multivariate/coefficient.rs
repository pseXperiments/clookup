use ff::Field;

use crate::poly::Polynomial;
use crate::utils::*;

use super::MultivariatePolynomial;

const TWO: usize = 2;

/// num_vars: number of variables
/// coefficients: [[coeff, deg(x_1), ..., deg(x_n)], ...] respectively
#[derive(Clone, Debug)]
pub struct CoefficientForm<F: Field> {
    num_vars: u32,
    coefficients: Vec<F>,
}

impl<F: Field> Polynomial<F> for CoefficientForm<F> {
    type Point = Vec<F>;

    fn evaluate(&self, point: &Self::Point) -> F {
        todo!()
    }
}

impl<F: Field> MultivariatePolynomial<F> for CoefficientForm<F> {
    fn interpolate(points: Vec<F>) -> Self {
        let num_vars = log(points.len());
        // TODO: interpolation can be done with 2^k values
        assert_eq!(TWO.pow(num_vars), points.len());
        
        Self { num_vars, coefficients: () }
    }
}
