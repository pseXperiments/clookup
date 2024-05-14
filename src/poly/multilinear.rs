use crate::utils::arithmetic::usize_from_bits_le;
use crate::utils::parallel::parallelize;
use crate::{
    core::precomputation::Table,
    poly::Polynomial,
    utils::{
        arithmetic::div_ceil,
        impl_index,
        parallel::{num_threads, parallelize_iter},
    },
};
use ff::Field;
use num::Integer;
use serde::{Deserialize, Serialize};
use std::mem;
use std::{
    borrow::{Borrow, Cow},
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign},
};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultilinearPolynomial<F> {
    evals: Vec<F>,
    coeffs: Vec<F>,
    num_vars: usize,
}

impl<F: Field> Default for MultilinearPolynomial<F> {
    fn default() -> Self {
        MultilinearPolynomial::zero()
    }
}

impl<F: Field> MultilinearPolynomial<F> {
    pub const fn zero() -> Self {
        Self {
            evals: vec![],
            coeffs: vec![],
            num_vars: 0,
        }
    }

    pub fn new(evals: Vec<F>, coeffs: Vec<F>, num_vars: usize) -> Self {
        MultilinearPolynomial {
            evals,
            coeffs,
            num_vars,
        }
    }

    fn is_empty(&self) -> bool {
        self.evals.is_empty()
    }

    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn evals(&self) -> &[F] {
        &self.evals
    }

    pub fn into_evals(self) -> Vec<F> {
        self.evals
    }

    pub fn iter(&self) -> impl Iterator<Item = &F> {
        self.evals.iter()
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
        Self::new(eval.clone(), result, num_vars)
    }

    pub fn evaluate(&self, point: &[F]) -> F {
        assert_eq!(point.len(), self.num_vars);

        let mut evals = Cow::Borrowed(self.evals());
        let mut bits = Vec::new();
        let mut buf = Vec::with_capacity(self.evals.len() >> 1);
        for x_i in point.iter() {
            if x_i == &F::ZERO || x_i == &F::ONE {
                bits.push(x_i == &F::ONE);
                continue;
            }

            let distance = bits.len() + 1;
            let skip = usize_from_bits_le(&bits);
            merge_in_place(&mut evals, x_i, distance, skip, &mut buf);
            bits.clear();
        }

        evals[usize_from_bits_le(&bits)]
    }

    pub fn eq_xy(y: &[F]) -> Self {
        if y.is_empty() {
            return Self::zero();
        }

        let expand_serial = |next_evals: &mut [F], evals: &[F], y_i: &F| {
            for (next_evals, eval) in next_evals.chunks_mut(2).zip(evals.iter()) {
                next_evals[1] = *eval * y_i;
                next_evals[0] = *eval - &next_evals[1];
            }
        };

        let mut evals = vec![F::ONE];
        for y_i in y.iter().rev() {
            let mut next_evals = vec![F::ZERO; 2 * evals.len()];
            if evals.len() < 32 {
                expand_serial(&mut next_evals, &evals, y_i);
            } else {
                let mut chunk_size = div_ceil(evals.len(), num_threads());
                if chunk_size.is_odd() {
                    chunk_size += 1;
                }
                parallelize_iter(
                    next_evals
                        .chunks_mut(chunk_size)
                        .zip(evals.chunks(chunk_size >> 1)),
                    |(next_evals, evals)| expand_serial(next_evals, evals, y_i),
                );
            }
            evals = next_evals;
        }

        Self {
            evals,
            coeffs: vec![],
            num_vars: y.len(),
        }
    }
}

impl<F: Field> From<Table<F>> for MultilinearPolynomial<F> {
    fn from(value: Table<F>) -> Self {
        Self::eval_to_coeff(value.table(), value.num_vars())
    }
}

impl_index!(MultilinearPolynomial, evals);

impl<F: Field, P: Borrow<MultilinearPolynomial<F>>> Add<P> for &MultilinearPolynomial<F> {
    type Output = MultilinearPolynomial<F>;

    fn add(self, rhs: P) -> MultilinearPolynomial<F> {
        let mut output = self.clone();
        output += rhs;
        output
    }
}

impl<F: Field, P: Borrow<MultilinearPolynomial<F>>> AddAssign<P> for MultilinearPolynomial<F> {
    fn add_assign(&mut self, rhs: P) {
        let rhs = rhs.borrow();
        match (self.is_empty(), rhs.is_empty()) {
            (_, true) => {}
            (true, false) => *self = rhs.clone(),
            (false, false) => {
                assert_eq!(self.num_vars, rhs.num_vars);
                parallelize(&mut self.evals, |(lhs, start)| {
                    for (lhs, rhs) in lhs.iter_mut().zip(rhs[start..].iter()) {
                        *lhs += rhs;
                    }
                });
                parallelize(&mut self.coeffs, |(lhs, start)| {
                    for (lhs, rhs) in lhs.iter_mut().zip(rhs[start..].iter()) {
                        *lhs += rhs;
                    }
                });
            }
        }
    }
}

impl<F: Field, BF: Borrow<F>, P: Borrow<MultilinearPolynomial<F>>> AddAssign<(BF, P)>
    for MultilinearPolynomial<F>
{
    fn add_assign(&mut self, (scalar, rhs): (BF, P)) {
        let (scalar, rhs) = (scalar.borrow(), rhs.borrow());
        match (self.is_empty(), rhs.is_empty() | (scalar == &F::ZERO)) {
            (_, true) => {}
            (true, false) => {
                *self = rhs.clone();
                *self *= scalar;
            }
            (false, false) => {
                assert_eq!(self.num_vars, rhs.num_vars);

                if scalar == &F::ONE {
                    *self += rhs;
                } else if scalar == &-F::ONE {
                    *self -= rhs;
                } else {
                    parallelize(&mut self.evals, |(lhs, start)| {
                        for (lhs, rhs) in lhs.iter_mut().zip(rhs[start..].iter()) {
                            *lhs += &(*scalar * rhs);
                        }
                    });
                    parallelize(&mut self.coeffs, |(lhs, start)| {
                        for (lhs, rhs) in lhs.iter_mut().zip(rhs[start..].iter()) {
                            *lhs += &(*scalar * rhs);
                        }
                    });
                }
            }
        }
    }
}

impl<F: Field, P: Borrow<MultilinearPolynomial<F>>> Sub<P> for &MultilinearPolynomial<F> {
    type Output = MultilinearPolynomial<F>;

    fn sub(self, rhs: P) -> MultilinearPolynomial<F> {
        let mut output = self.clone();
        output -= rhs;
        output
    }
}

impl<F: Field, P: Borrow<MultilinearPolynomial<F>>> SubAssign<P> for MultilinearPolynomial<F> {
    fn sub_assign(&mut self, rhs: P) {
        let rhs = rhs.borrow();
        match (self.is_empty(), rhs.is_empty()) {
            (_, true) => {}
            (true, false) => {
                *self = rhs.clone();
                *self *= &-F::ONE;
            }
            (false, false) => {
                assert_eq!(self.num_vars, rhs.num_vars);

                parallelize(&mut self.evals, |(lhs, start)| {
                    for (lhs, rhs) in lhs.iter_mut().zip(rhs[start..].iter()) {
                        *lhs -= rhs;
                    }
                });

                parallelize(&mut self.coeffs, |(lhs, start)| {
                    for (lhs, rhs) in lhs.iter_mut().zip(rhs[start..].iter()) {
                        *lhs -= rhs;
                    }
                });
            }
        }
    }
}

impl<F: Field, BF: Borrow<F>, P: Borrow<MultilinearPolynomial<F>>> SubAssign<(BF, P)>
    for MultilinearPolynomial<F>
{
    fn sub_assign(&mut self, (scalar, rhs): (BF, P)) {
        *self += (-*scalar.borrow(), rhs);
    }
}

impl<F: Field, BF: Borrow<F>> Mul<BF> for &MultilinearPolynomial<F> {
    type Output = MultilinearPolynomial<F>;

    fn mul(self, rhs: BF) -> MultilinearPolynomial<F> {
        let mut output = self.clone();
        output *= rhs;
        output
    }
}

impl<F: Field, BF: Borrow<F>> MulAssign<BF> for MultilinearPolynomial<F> {
    fn mul_assign(&mut self, rhs: BF) {
        let rhs = rhs.borrow();
        if rhs == &F::ZERO {
            self.evals = vec![F::ZERO; self.evals.len()];
            self.coeffs = vec![F::ZERO; self.coeffs.len()];
        } else if rhs == &-F::ONE {
            parallelize(&mut self.evals, |(evals, _)| {
                for eval in evals.iter_mut() {
                    *eval = -*eval;
                }
            });
            parallelize(&mut self.coeffs, |(coeffs, _)| {
                for coeff in coeffs.iter_mut() {
                    *coeff = -*coeff;
                }
            });
        } else if rhs != &F::ONE {
            parallelize(&mut self.evals, |(lhs, _)| {
                for lhs in lhs.iter_mut() {
                    *lhs *= rhs;
                }
            });
            parallelize(&mut self.coeffs, |(lhs, _)| {
                for lhs in lhs.iter_mut() {
                    *lhs *= rhs;
                }
            });
        }
    }
}

impl<F: Field> Polynomial<F> for MultilinearPolynomial<F> {
    type Point = Vec<F>;

    fn coeffs(&self) -> &[F] {
        self.evals.as_slice()
    }

    fn evaluate(&self, point: &Self::Point) -> F {
        MultilinearPolynomial::evaluate(self, point)
    }
}

pub(crate) fn merge_into<F: Field>(
    target: &mut Vec<F>,
    evals: &[F],
    x_i: &F,
    distance: usize,
    skip: usize,
) {
    assert!(target.capacity() >= evals.len() >> distance);
    target.resize(evals.len() >> distance, F::ZERO);

    let step = 1 << distance;
    parallelize(target, |(target, start)| {
        let start = (start << distance) + skip;
        for (target, (eval_0, eval_1)) in
            target.iter_mut().zip(zip_self!(evals.iter(), step, start))
        {
            *target = (*eval_1 - eval_0) * x_i + eval_0;
        }
    });
}

fn merge_in_place<F: Field>(
    evals: &mut Cow<[F]>,
    x_i: &F,
    distance: usize,
    skip: usize,
    buf: &mut Vec<F>,
) {
    merge_into(buf, evals, x_i, distance, skip);
    if let Cow::Owned(_) = evals {
        mem::swap(evals.to_mut(), buf);
    } else {
        *evals = mem::replace(buf, Vec::with_capacity(buf.len() >> 1)).into();
    }
}

macro_rules! zip_self {
    (@ $iter:expr, $step:expr, $skip:expr) => {
        $iter.skip($skip).step_by($step).zip($iter.skip($skip + ($step >> 1)).step_by($step))
    };
    ($iter:expr) => {
        zip_self!(@ $iter, 2, 0)
    };
    ($iter:expr, $step:expr) => {
        zip_self!(@ $iter, $step, 0)
    };
    ($iter:expr, $step:expr, $skip:expr) => {
        zip_self!(@ $iter, $step, $skip)
    };
}

pub(crate) use zip_self;

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
