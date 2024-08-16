#![feature(generic_const_exprs)]
#![feature(allocator_api)]
pub use circuit_definitions::snark_wrapper::franklin_crypto as franklin_crypto;
pub use franklin_crypto::bellman;

use bellman::{
    bn256::{Bn256, Fr},
    kate_commitment::{commit_using_monomials, Crs, CrsForMonomialForm},
    pairing::ff::{Field, PrimeField},
    pairing::{CurveAffine, CurveProjective},
    plonk::{
        better_better_cs::{
            cs::{
                ensure_in_map_or_create, get_from_map_unchecked, AssembledPolynomialStorage,
                AssembledPolynomialStorageForMonomialForms, Assembly, Circuit, Gate, GateInternal,
                MainGate, PlonkConstraintSystemParams, PolyIdentifier, PolynomialInConstraint,
                PolynomialProxy, Setup, SynthesisMode, SynthesisModeTesting,
            },
            gates::naive_main_gate::NaiveMainGate,
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

mod definitions;
pub use definitions::*;
mod prover;
use prover::*;
pub mod verifier;
pub use utils::*;
pub use verifier::*;
pub mod convenience;
pub use convenience::*;

#[cfg(test)]
mod test;
mod utils;
use utils::*;
pub use utils::{
    compute_generators, compute_power_of_two_root_of_generator, num_system_polys_from_vk,
};

pub(crate) const SANITY_CHECK: bool = true;
// pub use shivini::circuit_definitions as circuit_definitions;
