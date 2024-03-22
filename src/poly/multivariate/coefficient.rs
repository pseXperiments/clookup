use ff::Field;

use crate::poly::Polynomial;

use super::MultivariatePolynomial;

#[derive(Clone, Debug)]
pub struct CoefficientForm<F: Field> {
    variables: u64,
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
        todo!()
    }
}
