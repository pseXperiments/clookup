use ff::Field;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use std::fmt::Debug;

use crate::core::precomputation::Table;

#[derive(Debug, Clone)]
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

    pub fn eval_to_coeff(eval: &Vec<F>, num_vars: usize) -> MultilinearPolynomial<F> {
        let mut result = eval.clone();
        for i in (0..num_vars).rev() {
            let chunk_size = 2usize.pow((i + 1) as u32);
            for chunk_offset in (0..eval.len()).step_by(chunk_size) {
                for j in chunk_offset..(chunk_offset + chunk_size / 2) {
                    result[j + chunk_size / 2] = result[j + chunk_size / 2] - result[j];
                }
            }
        }
        Self::new(result, num_vars)
    }
}

impl<F: Field> From<Table<F>> for MultilinearPolynomial<F> {
    fn from(value: Table<F>) -> Self {
        Self::eval_to_coeff(&value.table, value.num_vars)
    }
}

#[cfg(test)]
mod test {
    use super::MultilinearPolynomial;
    use halo2curves::bn256::Fr;
    #[test]
    fn test_conversion() {
        let poly = vec![Fr::from(1), Fr::from(3), Fr::from(5), Fr::from(7)];
        let res = vec![Fr::from(1), Fr::from(2), Fr::from(4), Fr::from(0)];
        assert_eq!(MultilinearPolynomial::eval_to_coeff(&poly, 2).coeffs, res);
    }
}
