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

    fn eval_to_coeff(poly: &[F], k: u32) -> MultilinearPolynomial<F> {
        let mut result = poly.to_vec();
        for i in (0..k).rev() {
            let chunk_size = 2usize.pow(i + 1);
            for chunk_offset in (0..poly.len()).step_by(chunk_size) {
                for j in chunk_offset..(chunk_offset + chunk_size / 2) {
                    result[j + chunk_size / 2] = result[j + chunk_size / 2] - result[j];
                }
            }
        }
        Self::new(result, 2usize.pow(k))
    }
}

#[cfg(test)]
mod test {
    use halo2curves::bn256::Fr;
    use super::MultilinearPolynomial;
    #[test]
    fn test_conversion() {
        let poly = [Fr::from(1), Fr::from(3), Fr::from(5), Fr::from(7)];
        let res = vec![Fr::from(1), Fr::from(2), Fr::from(4), Fr::from(0)];
        assert_eq!(MultilinearPolynomial::eval_to_coeff(&poly, 2).coeffs, res);
    }
}
