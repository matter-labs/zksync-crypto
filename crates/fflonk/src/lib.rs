#![feature(generic_const_exprs)]
#![feature(allocator_api)]
pub use franklin_crypto;
pub use franklin_crypto::bellman;
use franklin_crypto::bellman::plonk::better_better_cs::{cs::PlonkCsWidth3Params, gates::naive_main_gate::NaiveMainGate};

use bellman::{
    bn256::{Bn256, Fr},
    kate_commitment::{commit_using_monomials, Crs, CrsForMonomialForm},
    pairing::ff::{Field, PrimeField},
    pairing::{CurveAffine, CurveProjective},
    plonk::{
        better_better_cs::{
            cs::{
                ensure_in_map_or_create, get_from_map_unchecked, AssembledPolynomialStorage, AssembledPolynomialStorageForMonomialForms, Assembly, Circuit, Gate, GateInternal, MainGate,
                PlonkConstraintSystemParams, PolyIdentifier, PolynomialInConstraint, PolynomialProxy, Setup, SynthesisMode, SynthesisModeTesting,
            },
            gates::selector_optimized_with_d_next::SelectorOptimizedWidth4MainGateWithDNext,
            utils::BinopAddAssignScaled,
        },
        better_cs::generator::make_non_residues,
        commitments::transcript::Transcript,
        domains::Domain,
        fft::cooley_tukey_ntt::{BitReversedOmegas, CTPrecomputations, OmegasInvBitreversed},
        polynomials::{Coefficients, Polynomial, Values},
    },
    worker::Worker,
    Engine, ScalarEngine, SynthesisError,
};
pub use franklin_crypto::plonk::circuit::custom_rescue_gate::Rescue5CustomGate;

pub mod rand {
    pub use crate::franklin_crypto::bellman::pairing::ff::rand::Rng;
    pub use crate::franklin_crypto::bellman::pairing::ff::Rand;
    pub use rand::{distributions, random, rngs, seq, thread_rng, RngCore, SeedableRng};

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

        fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), ::rand::Error> {
            self.0.try_fill_bytes(dest)
        }
    }

    impl SeedableRng for XorShiftRng {
        type Seed = <rand_xorshift::XorShiftRng as SeedableRng>::Seed;

        fn from_seed(seed: Self::Seed) -> Self {
            Self(rand_xorshift::XorShiftRng::from_seed(seed))
        }
    }
}

mod definitions;
pub use definitions::*;
pub mod prover;
use prover::*;
pub mod utils;
pub use utils::*;
pub mod verifier;
pub use verifier::*;

#[cfg(test)]
mod test;

pub const L1_VERIFIER_DOMAIN_SIZE_LOG: usize = 23;
pub const MAX_COMBINED_DEGREE_FACTOR: usize = 9;
pub(crate) const SANITY_CHECK: bool = true;
#[cfg(feature = "allocator")]
pub type FflonkAssembly<E, S, A = std::alloc::Global> = Assembly<E, PlonkCsWidth3Params, NaiveMainGate, S, A>;
#[cfg(not(feature = "allocator"))]
pub type FflonkAssembly<E, S> = Assembly<E, PlonkCsWidth3Params, NaiveMainGate, S>;
