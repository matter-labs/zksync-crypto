#![allow(dead_code, unused_imports, unused_macros)]
#![allow(macro_expanded_macro_exports_accessed_by_absolute_paths)]
#![warn(unused_assignments)]
#![feature(array_chunks)]

pub extern crate bellman;
extern crate blake2;
extern crate blake2_rfc_bellman_edition as blake2_rfc;
pub extern crate boojum;
extern crate byteorder;
extern crate core;
extern crate derivative;
extern crate digest;
extern crate indexmap;
extern crate itertools;
extern crate num_bigint;
extern crate num_derive;
extern crate num_integer;
extern crate num_traits;
extern crate rand;
extern crate serde;
extern crate sha2;
extern crate sha3;
extern crate splitmut;
extern crate tiny_keccak;

use bellman::pairing;
use bellman::pairing::ff;

#[macro_use]
extern crate lazy_static;

#[macro_use]
extern crate arr_macro;

#[cfg(test)]
extern crate hex;

pub mod alt_babyjubjub;
pub mod as_waksman;
pub mod constants;
pub mod generic_twisted_edwards;
pub mod group_hash;
pub mod interpolation;
pub mod jubjub;
pub mod pedersen_hash;
pub mod plonk;
pub mod primitives;
pub mod redjubjub;
pub mod rescue;
pub mod util;

pub fn log2_floor(num: usize) -> u32 {
    assert!(num > 0);

    let mut pow = 0;

    while (1 << (pow + 1)) <= num {
        pow += 1;
    }

    pow
}
