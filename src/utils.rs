pub use itertools::{chain, izip, Either, Itertools};
pub use num_bigint::BigUint;
pub use serde::{de::DeserializeOwned, Deserialize, Deserializer, Serialize, Serializer};
pub use timer::{end_timer, start_timer, start_unit_timer};

pub mod arithmetic;
pub mod parallel;
pub mod timer;

#[derive(Debug)]
pub enum ProtocolError {
    InvalidSumcheck(String),
    InvalidPcsParam(String),
    InvalidPcsOpen(String),
    SizeError,
    NotInclusion,
}

pub fn transpose<T>(v: Vec<Vec<T>>) -> Vec<Vec<T>> {
    let len = v[0].len();
    let mut iters: Vec<_> = v.into_iter().map(|n| n.into_iter()).collect();
    (0..len)
        .map(|_| {
            iters
                .iter_mut()
                .map(|n| n.next().unwrap())
                .collect::<Vec<T>>()
        })
        .collect()
}

macro_rules! impl_index {
    (@ $name:ty, $field:tt, [$($range:ty => $output:ty),*$(,)?]) => {
        $(
            impl<F> std::ops::Index<$range> for $name {
                type Output = $output;

                fn index(&self, index: $range) -> &$output {
                    self.$field.index(index)
                }
            }

            impl<F> std::ops::IndexMut<$range> for $name {
                fn index_mut(&mut self, index: $range) -> &mut $output {
                    self.$field.index_mut(index)
                }
            }
        )*
    };
    (@ $name:ty, $field:tt) => {
        impl_index!(
            @ $name, $field,
            [
                usize => F,
                std::ops::Range<usize> => [F],
                std::ops::RangeFrom<usize> => [F],
                std::ops::RangeFull => [F],
                std::ops::RangeInclusive<usize> => [F],
                std::ops::RangeTo<usize> => [F],
                std::ops::RangeToInclusive<usize> => [F],
            ]
        );
    };
    ($name:ident, $field:tt) => {
        impl_index!(@ $name<F>, $field);
    };
}

macro_rules! izip_eq {
    (@closure $p:pat => $tup:expr) => {
        |$p| $tup
    };
    (@closure $p:pat => ($($tup:tt)*) , $_iter:expr $(, $tail:expr)*) => {
        $crate::utils::izip_eq!(@closure ($p, b) => ($($tup)*, b) $(, $tail)*)
    };
    ($first:expr $(,)*) => {
        itertools::__std_iter::IntoIterator::into_iter($first)
    };
    ($first:expr, $second:expr $(,)*) => {
        $crate::utils::izip_eq!($first).zip_eq($second)
    };
    ($first:expr $(, $rest:expr)* $(,)*) => {
        $crate::util::izip_eq!($first)
            $(.zip_eq($rest))*
            .map($crate::utils::izip_eq!(@closure a => (a) $(, $rest)*))
    };
}

pub(crate) use {impl_index, izip_eq};
