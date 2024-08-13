use super::sbox::sbox;
use super::sponge::circuit_generic_hash_num;
use super::matrix::{matrix_vector_product, mul_by_sparse_matrix};
use crate::{DomainStrategy, poseidon::params::PoseidonParams};
use crate::poseidon2::Poseidon2Params;
use crate::traits::{HashFamily, HashParams};
use franklin_crypto::bellman::plonk::better_better_cs::cs::ConstraintSystem;
use franklin_crypto::bellman::{Field, SynthesisError};
use franklin_crypto::{
    bellman::Engine,
    plonk::circuit::{allocated_num::Num, linear_combination::LinearCombination},
};

/// Receives inputs whose length `known` prior(fixed-length).
/// Also uses custom domain strategy which basically sets value of capacity element to
/// length of input and applies a padding rule which makes input size equals to multiple of
/// rate parameter.
/// Uses pre-defined state-width=3 and rate=2.
pub fn circuit_poseidon2_hash<E: Engine, CS: ConstraintSystem<E>, const L: usize>(
    cs: &mut CS,
    input: &[Num<E>; L],
    domain_strategy: Option<DomainStrategy>,
) -> Result<[Num<E>; 2], SynthesisError> {
    const WIDTH: usize = 3;
    const RATE: usize = 2;
    let params = Poseidon2Params::<E, RATE, WIDTH>::default();
    circuit_generic_hash_num(cs, input, &params, domain_strategy)
}

pub fn circuit_poseidon2_round_function<
    E: Engine,
    CS: ConstraintSystem<E>,
    const RATE: usize,
    const WIDTH: usize,
>(
    cs: &mut CS,
    params: &Poseidon2Params<E, RATE, WIDTH>,
    state: &mut [LinearCombination<E>; WIDTH],
) -> Result<(), SynthesisError> {
    assert!(params.number_of_full_rounds() % 2 == 0);

    let half_of_full_rounds = params.number_of_full_rounds() / 2;

    // Linear layer at beginning
    matrix_vector_product(&params.mds_external_matrix, state)?;

    // first full rounds
    for round in 0..half_of_full_rounds {
        let round_constants = &params.round_constants[round];

        // add round constatnts
        for (s, c) in state.iter_mut().zip(round_constants.iter()) {
            s.add_assign_constant(*c);
        }
        // non linear sbox
        sbox(
            cs,
            params.alpha(),
            state,
            Some(0..WIDTH),
            params.custom_gate(),
        )?;

        // mul state by mds
        matrix_vector_product(&params.mds_external_matrix, state)?;
    }

    let mut diag_internal_matrix_decreased = params.diag_internal_matrix.clone();
    for coeff in diag_internal_matrix_decreased.iter_mut() {
        coeff.sub_assign(&E::Fr::one());
    }

    for round in half_of_full_rounds..(params.partial_rounds + half_of_full_rounds) {
        // add round constatnt
        let round_constant = params.round_constants[round][0];
        state[0].add_assign_constant(round_constant);

        // non linear sbox
        sbox(cs, params.alpha(), state, Some(0..1), params.custom_gate())?;

        // mul state by internal matrix
        let mut sum = state[0].clone();
        for s in state.iter().skip(1) {
            sum.add_assign(s);
        }

        for (s, coeff) in state.iter_mut().zip(diag_internal_matrix_decreased.iter()) {
            s.scale(coeff);
            s.add_assign(&sum);
        }
    }

    // second full round
    for round in (params.number_of_partial_rounds() + half_of_full_rounds)
        ..(params.number_of_partial_rounds() + params.number_of_full_rounds())
    {
        let round_constants = &params.round_constants[round];

        // add round constatnts
        for (s, c) in state.iter_mut().zip(round_constants.iter()) {
            s.add_assign_constant(*c);
        }
        // non linear sbox
        sbox(
            cs,
            params.alpha(),
            state,
            Some(0..WIDTH),
            params.custom_gate(),
        )?;

        // mul state by mds
        matrix_vector_product(&params.mds_external_matrix, state)?;
    }

    Ok(())
}
