use crate::{poly::multilinear::MultilinearPolynomial, utils::ProtocolError};
use ff::{Field, PrimeField};
use num::Integer;
use std::{collections::HashMap, hash::Hash};

#[derive(Clone, Debug)]
pub struct Table<F> {
    pub table: Vec<F>,
    /// Table size should be 2^k
    pub num_vars: usize,
    index_map: HashMap<F, usize>,
}

impl<F: Field + Hash> TryFrom<Vec<F>> for Table<F> {
    type Error = ProtocolError;
    fn try_from(table: Vec<F>) -> Result<Self, Self::Error> {
        let num_vars = table.len().ilog2() as usize;
        if 1 << num_vars != table.len() {
            Err(ProtocolError::SizeError)
        } else {
            let mut index_map = HashMap::new();
            for (index, &element) in table.iter().enumerate() {
                index_map.insert(element, index);
            }
            Ok(Self {
                table,
                num_vars,
                index_map,
            })
        }
    }
}

impl<F: PrimeField + Hash> Table<F> {
    pub fn find_indices(&self, elements: &Vec<F>) -> Result<Vec<Vec<F>>, ProtocolError> {
        let length = self.num_vars;
        let mut indices = Vec::new();
        for elem in elements {
            if let Some(idx) = self.index_map.get(elem) {
                let binary_string = format!("{:0width$b}", idx, width = length);
                let mut binary_vec = Vec::new();
                for c in binary_string.chars() {
                    binary_vec.push(F::from(c.to_digit(10).unwrap() as u64));
                }
                indices.push(binary_vec);
            } else {
                return Err(ProtocolError::NotInclusion);
            }
        }
        Ok(indices)
    }
}

#[derive(Clone, Debug)]
pub struct TablePolynomial<F> {
    poly: MultilinearPolynomial<F>,
}
