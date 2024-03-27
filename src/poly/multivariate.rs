use std::fmt::Debug;
use ff::Field;

/// [coeff, deg(x_1), ..., deg(x_n)]
#[derive(Debug, Clone)]
struct Term<F: Field>(Vec<F>);

/// num_vars: number of variables
/// coefficients: [[coeff, deg(x_1), ..., deg(x_n)], ...] respectively
#[derive(Clone, Debug)]
struct CoefficientForm<F: Field> {
    num_vars: u32,
    coefficients: Vec<Term<F>>,
}

/// format: f(g(h(..)))
/// outer: f
/// inner: g(..)
#[derive(Debug, Clone)]
struct CompositionForm {
}

#[derive(Debug, Clone)]
struct FactoredForm<F: Field> {
    terms: Vec<Vec<Term<F>>>,
}

#[derive(Debug, Clone)]
pub enum MultivariatePolynomial<F: Field> {
    Coeff(CoefficientForm<F>),
    Comp(CompositionForm),
    Fac(FactoredForm<F>),
}

impl<F: Field> MultivariatePolynomial<F>  {
    fn interpolate(points: Vec<F>, k: u32) -> Self {
        todo!()
    }
}
