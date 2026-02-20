#![feature(array_chunks)]
#![feature(allocator_api)]
#![feature(type_changing_struct_update)]

pub mod traits;
pub mod verifier;
pub mod verifier_structs;

pub mod implementations;

pub extern crate rescue_poseidon;
pub use franklin_crypto::boojum;
pub use rescue_poseidon::franklin_crypto;

pub mod rand {
    pub use crate::franklin_crypto::bellman::pairing::ff::rand::Rng;
    pub use crate::franklin_crypto::bellman::pairing::ff::Rand;
    pub use rand::{distributions, random, rngs, seq, thread_rng, RngCore, SeedableRng};
}
