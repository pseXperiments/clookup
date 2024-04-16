pub mod multilinear;

use crate::ff::Field;
use std::fmt::Debug;

pub trait Polynomial<F: Field>: Clone + Debug {
    type Point: Clone + Debug;

    fn evaluate(&self, point: &Self::Point) -> F;
}
