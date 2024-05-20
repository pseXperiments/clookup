use crate::{poly::multilinear::MultilinearPolynomial, utils::ProtocolError};
use ff::{Field, PrimeField};
use std::{collections::HashMap, hash::Hash};

#[derive(Clone, Debug)]
pub struct Table<F> {
    table: Vec<F>,
    /// Table size should be 2^k
    num_vars: usize,
    index_map: HashMap<F, usize>,
}

impl<F> Table<F> {
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn len(&self) -> usize {
        assert_eq!(self.table.len(), 1 << self.num_vars);
        self.table.len()
    }

    pub fn table(&self) -> &Vec<F> {
        &self.table
    }
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
                binary_vec.reverse();
                indices.push(binary_vec);
            } else {
                return Err(ProtocolError::NotInclusion);
            }
        }
        Ok(indices)
    }

    pub fn polynomial(&self) -> MultilinearPolynomial<F> {
        MultilinearPolynomial::eval_to_coeff(&self.table, self.num_vars)
    }
}

#[derive(Clone, Debug)]
pub struct TablePolynomial<F> {
    poly: MultilinearPolynomial<F>,
}

#[cfg(test)]
mod test {
    use super::Table;
    use crate::utils::{transpose, ProtocolError};
    use ff::Field;
    use halo2curves::bn256::Fr;

    #[test]
    fn test_find_indices() -> Result<(), ProtocolError> {
        let table_vec = vec![
            Fr::from(11),
            Fr::from(5),
            Fr::from(3),
            Fr::from(17),
            Fr::from(2),
            Fr::from(13),
            Fr::from(7),
            Fr::from(19),
        ];
        let table: Table<Fr> = table_vec.try_into()?;
        let witness = vec![Fr::from(2), Fr::from(3), Fr::from(5), Fr::from(7)];
        let indices = table.find_indices(&witness)?;

        let res = vec![
            vec![Fr::ZERO, Fr::ZERO, Fr::ONE],
            vec![Fr::ZERO, Fr::ONE, Fr::ZERO],
            vec![Fr::ONE, Fr::ZERO, Fr::ZERO],
            vec![Fr::ZERO, Fr::ONE, Fr::ONE],
        ];
        assert_eq!(indices, res);
        Ok(())
    }
}
