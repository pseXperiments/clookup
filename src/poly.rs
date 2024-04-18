use crate::ff::Field;
use std::{fmt::Debug, ops::AddAssign};

pub mod multilinear;

pub trait Polynomial<F: Field>:
    Clone + Debug + Default + for<'a> AddAssign<(&'a F, &'a Self)>
{
    type Point: Clone + Debug;

    fn coeffs(&self) -> &[F];

    fn evaluate(&self, point: &Self::Point) -> F;
}
