use ff::Field;
use crate::utils::{ProtocolError, TWO, log};

#[derive(Clone, Debug)]
pub struct Table<F: Field> {
    table: Vec<F>,
    /// Table size should be 2^k
    /// pow_vars = k
    exp_vars: u32
}

impl<F: Field> TryFrom<Vec<F>> for Table<F> {
    type Error = ProtocolError;
    fn try_from(table: Vec<F>) -> Result<Self, Self::Error> {
        let exp_vars = log(table.len());
        if TWO.pow(exp_vars) != table.len() {
            Err(ProtocolError::SizeError)
        } else {
            Ok(Self { table, exp_vars })
        }
    }
}
