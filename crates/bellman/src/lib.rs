#![allow(dead_code, unused_imports, unused_mut, unused_variables, unused_macros, unused_assignments, unreachable_patterns)]
#![cfg_attr(feature = "allocator", feature(allocator_api))]
#![allow(clippy::needless_borrow, clippy::needless_borrows_for_generic_args)]

#[macro_use]
extern crate cfg_if;
extern crate bit_vec;
extern crate byteorder;
pub extern crate pairing;
pub extern crate rand as rand_crate;

pub use pairing::*;
pub use smallvec;

use crate::pairing::ff;
pub use ff::*;

pub mod rand {
    pub use crate::pairing::ff::rand::Rng;
    pub use crate::pairing::ff::Rand;
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
mod log;

pub mod domain;
pub mod groth16;

#[cfg(feature = "gm17")]
pub mod gm17;

#[cfg(feature = "sonic")]
pub mod sonic;

#[cfg(feature = "plonk")]
pub mod plonk;

#[macro_use]
#[cfg(feature = "plonk")]
extern crate lazy_static;

#[cfg(feature = "marlin")]
pub mod marlin;

#[cfg(any(feature = "marlin", feature = "plonk"))]
pub mod kate_commitment;

pub mod constants;
mod group;
mod multiexp;
mod prefetch;
mod source;

#[cfg(test)]
mod tests;

cfg_if! {
    if #[cfg(feature = "multicore")] {
        #[cfg(feature = "wasm")]
        compile_error!("Multicore feature is not yet compatible with wasm target arch");

        mod multicore;
        pub mod worker {
            pub use super::multicore::*;
        }
    } else {
        mod singlecore;
        pub mod worker {
            pub use super::singlecore::*;
        }
    }
}

mod cs;
pub use self::cs::*;

use std::env;
use std::str::FromStr;

cfg_if! {
    if #[cfg(any(not(feature = "nolog"), feature = "sonic"))] {
        fn verbose_flag() -> bool {
            option_env!("BELLMAN_VERBOSE").unwrap_or("0") == "1"
        }
    }
}
