use super::*;

const TABLE_NAME_16: &str = "Range check 16 bits table";

#[derive(Derivative)]
#[derivative(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RangeCheck16BitsTable<const N: usize>;

pub fn create_range_check_16_bits_table<const N: usize, F: SmallField>() -> LookupTable<F, N> {
    let mut all_keys = Vec::with_capacity(1 << 16);
    for a in 0..=u16::MAX {
        let key = smallvec::smallvec![F::from_u64_unchecked(a as u64),];
        all_keys.push(key);
    }
    LookupTable::new_from_keys_and_generation_function(
        &all_keys,
        TABLE_NAME_16.to_string(),
        1,
        |_keys| smallvec::smallvec![F::ZERO; N-1],
    )
}

const TABLE_NAME_15: &str = "Range check 15 bits table";

#[derive(Derivative)]
#[derivative(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RangeCheck15BitsTable<const N: usize>;

pub fn create_range_check_15_bits_table<const N: usize, F: SmallField>() -> LookupTable<F, N> {
    let mut all_keys = Vec::with_capacity(1 << 15);
    for a in 0..(1 << 15) {
        let key = smallvec::smallvec![F::from_u64_unchecked(a as u64),];
        all_keys.push(key);
    }
    LookupTable::new_from_keys_and_generation_function(
        &all_keys,
        TABLE_NAME_15.to_string(),
        1,
        |_keys| smallvec::smallvec![F::ZERO; N-1],
    )
}
