use ff::Field;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use std::fmt::Debug;

#[derive(Debug)]
pub struct MultilinearPolynomial<F> {
    coeffs: Vec<F>,
    num_vars: usize,
}

// #[derive(Debug, Clone)]
// pub enum MultivariatePolynomial<F: Field> {
//     Coeff(CoefficientForm<F>),
// }

impl<F: Field> MultilinearPolynomial<F> {
    fn new(coeffs: Vec<F>, num_vars: usize) -> Self {
        MultilinearPolynomial { coeffs, num_vars }
    }

    fn eval_to_coeff(poly: &[F], num_vars: usize) -> MultilinearPolynomial<F> {
        let mut result = poly.to_vec();
        for i in (0..num_vars).rev() {
            let chunk_size = 2usize.pow((i + 1) as u32);
            for chunk_offset in (0..poly.len()).step_by(chunk_size) {
                for j in chunk_offset..(chunk_offset + chunk_size / 2) {
                    result[j + chunk_size / 2] -= poly[j];
                }
            }
        }
        Self::new(result, num_vars)
    }
}
