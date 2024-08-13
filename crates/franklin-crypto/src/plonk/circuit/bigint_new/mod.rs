use bellman::plonk::better_better_cs::cs::LookupTableApplication;

use crate::bellman::pairing::Engine;
use crate::bellman::pairing::ff::{Field, PrimeField, PrimeFieldRepr, BitIterator};
use crate::bellman::SynthesisError;
use crate::bellman::plonk::better_better_cs::cs::{
    Variable, ConstraintSystem, ArithmeticTerm, MainGateTerm, Width4MainGateWithDNext, MainGate, GateInternal,
    Gate, LinearCombinationOfTerms, PolynomialMultiplicativeTerm, PolynomialInConstraint, TimeDilation,
    Coefficient, PlonkConstraintSystemParams, PlonkCsWidth4WithNextStepParams, TrivialAssembly
};
use bellman::plonk::better_better_cs::data_structures::*;
use crate::plonk::circuit::Assignment;
use crate::plonk::circuit::utils::u64_to_fe;
use super::allocated_num::*;
use super::boolean::*;
use std::iter::FromIterator;
use std::sync::Arc;

pub mod bigint;
pub mod range_check_custom_gate2;
pub mod range_checks;
pub mod range_check_table2;

pub use self::bigint::*;
pub use self::range_check_custom_gate2::*;
pub use self::range_checks::*;
pub use self::range_check_table2::*;

pub mod amplified_linear_combination;
pub mod field;
pub use self::amplified_linear_combination::*;
pub use self::field::*;


pub const BITWISE_LOGICAL_OPS_TABLE_NAME: &'static str = "Table for bitwise logical ops";
pub const DEFAULT_RANGE_TABLE_GRANULARITY: usize = 8;


// splits an element into slices of fixed bit widths in LE order
#[track_caller]
pub fn split_into_slices<F: PrimeField>(el: &F, slice_width: usize, num_slices: usize) -> Vec<F> {
    let mut repr = el.into_repr();
    assert!(repr.num_bits() as usize <= slice_width * num_slices);
    let mut slices = Vec::with_capacity(num_slices);
    if slice_width < 64 {    
        let mask = (1u64 << slice_width) - 1u64;
        for _ in 0..num_slices {
            let slice = repr.as_ref()[0] & mask;

            let mut r = F::Repr::default();
            r.as_mut()[0] = slice;

            let slice = F::from_repr(r).unwrap();
            slices.push(slice);

            repr.shr(slice_width as u32);
        }
    }
    else {
        let it = repr.as_ref().iter().map(|x| u64_to_fe::<F>(*x)).take(num_slices);
        slices.extend(it);
    };

    slices
}

#[track_caller]
pub fn split_some_into_slices<F: PrimeField>(
    el: Option<F>,
    slice_width: usize,
    num_slices: usize
) -> Vec<Option<F>> {
    if let Some(v) = el.as_ref() {
        split_into_slices(v, slice_width, num_slices).into_iter().map(|el| Some(el)).collect()
    } else {
        vec![None; num_slices]
    }
}


#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum RangeConstraintStrategy {
    NaiveSingleBit,
    CustomTwoBitGate,
    WithBitwiseOpTable(usize) // parameter here is the chunk width    
}

impl RangeConstraintStrategy {
    pub fn get_range_check_granularity(&self) -> usize {
        match self {
            RangeConstraintStrategy::NaiveSingleBit => 1,
            RangeConstraintStrategy::CustomTwoBitGate => 2,
            RangeConstraintStrategy::WithBitwiseOpTable(n) => *n,
        }
    }
}

pub fn get_optimal_strategy<E: Engine, CS: ConstraintSystem<E>>(cs: &CS) -> RangeConstraintStrategy {
    if let Ok(table) = cs.get_table(BITWISE_LOGICAL_OPS_TABLE_NAME) {
        let width = crate::log2_floor(table.size())/2;
        return RangeConstraintStrategy::WithBitwiseOpTable(width as usize);
    }  
    if CS::Params::STATE_WIDTH == 4 && CS::Params::HAS_CUSTOM_GATES {
        return RangeConstraintStrategy::CustomTwoBitGate
    }
    RangeConstraintStrategy::NaiveSingleBit
}

pub fn inscribe_default_bitop_range_table<E, CS>(cs: &mut CS) -> Result<Arc<LookupTableApplication<E>>, SynthesisError> 
where E: Engine, CS: ConstraintSystem<E>
{
    use crate::plonk::circuit::hashes_with_tables::get_or_create_table;

    let columns3 = vec![
        PolyIdentifier::VariablesPolynomial(0), 
        PolyIdentifier::VariablesPolynomial(1), 
        PolyIdentifier::VariablesPolynomial(2)
    ];

    get_or_create_table(
        cs, BITWISE_LOGICAL_OPS_TABLE_NAME, || {
            LookupTableApplication::new(
                BITWISE_LOGICAL_OPS_TABLE_NAME, CombinedBitwiseLogicRangeTable::new(
                    BITWISE_LOGICAL_OPS_TABLE_NAME, DEFAULT_RANGE_TABLE_GRANULARITY,
                ),
                columns3, None, true
            )
        }
    )
}


pub(crate) fn compute_shifts<F: PrimeField>() -> Vec<F> {
    let mut result = Vec::with_capacity(F::CAPACITY as usize);
    let mut el = F::one();
    result.push(el);
    for _ in 1..F::CAPACITY {
        el.double();
        result.push(el);
    }

    result
}


pub(crate) fn round_up(x: usize, granularity: usize) -> usize {
    let rem = x % granularity;
    let to_add = if rem == 0 { 0 } else { granularity - rem };
    x + to_add
}


