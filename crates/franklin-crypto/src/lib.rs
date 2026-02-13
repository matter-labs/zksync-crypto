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
extern crate rand as rand_crate;
extern crate serde;
extern crate sha2;
extern crate sha3;
extern crate splitmut;
extern crate tiny_keccak;

use bellman::pairing;
use bellman::pairing::ff;

pub mod rand {
    pub use crate::bellman::pairing::ff::rand::Rng;
    pub use crate::bellman::pairing::ff::Rand;
    pub use crate::rand_crate::{distributions, random, rngs, seq, thread_rng, CryptoRng, RngCore, SeedableRng};

    #[derive(Clone, Debug)]
    pub struct XorShiftRng(pub rand_xorshift::XorShiftRng);

    impl XorShiftRng {
        pub fn from_seed(seed: [u32; 4]) -> Self {
            let mut seed_bytes = [0u8; 16];
            for (chunk, word) in seed_bytes.chunks_exact_mut(4).zip(seed.iter()) {
                chunk.copy_from_slice(&word.to_le_bytes());
            }

            <Self as SeedableRng>::from_seed(seed_bytes)
        }
    }

    impl RngCore for XorShiftRng {
        fn next_u32(&mut self) -> u32 {
            self.0.next_u32()
        }

        fn next_u64(&mut self) -> u64 {
            self.0.next_u64()
        }

        fn fill_bytes(&mut self, dest: &mut [u8]) {
            self.0.fill_bytes(dest)
        }

        fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), crate::rand_crate::Error> {
            self.0.try_fill_bytes(dest)
        }
    }

    impl SeedableRng for XorShiftRng {
        type Seed = <rand_xorshift::XorShiftRng as SeedableRng>::Seed;

        fn from_seed(seed: Self::Seed) -> Self {
            Self(rand_xorshift::XorShiftRng::from_seed(seed))
        }
    }

    #[derive(Clone, Debug)]
    pub struct ChaChaRng(pub rand_chacha::ChaChaRng);

    impl ChaChaRng {
        pub fn from_seed(seed: &[u32]) -> Self {
            let mut seed_bytes = [0u8; 32];
            for (chunk, word) in seed_bytes.chunks_exact_mut(4).zip(seed.iter().take(8)) {
                chunk.copy_from_slice(&word.to_le_bytes());
            }

            <Self as SeedableRng>::from_seed(seed_bytes)
        }
    }

    impl RngCore for ChaChaRng {
        fn next_u32(&mut self) -> u32 {
            self.0.next_u32()
        }

        fn next_u64(&mut self) -> u64 {
            self.0.next_u64()
        }

        fn fill_bytes(&mut self, dest: &mut [u8]) {
            self.0.fill_bytes(dest)
        }

        fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), crate::rand_crate::Error> {
            self.0.try_fill_bytes(dest)
        }
    }

    impl CryptoRng for ChaChaRng {}

    impl SeedableRng for ChaChaRng {
        type Seed = <rand_chacha::ChaChaRng as SeedableRng>::Seed;

        fn from_seed(seed: Self::Seed) -> Self {
            Self(rand_chacha::ChaChaRng::from_seed(seed))
        }
    }

    pub mod chacha {
        pub use super::ChaChaRng;
    }
}

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
