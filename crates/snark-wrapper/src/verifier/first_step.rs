use super::*;

/// First of all we should absorb commitments to transcript
/// and get challenges in the following order:
/// - absorb  setup commitment
/// - absorb  public inputs
/// - absorb  witness commitment
/// - get     beta and gamma challenges
/// - get     lookup_beta and lookup_gamma challenges
/// - absorb  stage_2 commitment
/// - get     alpha challenge
/// - absorb  quotient commitment
/// - get     z challenge
/// - absorb  evaluations at z
pub(crate) fn verify_first_step<E: Engine, CS: ConstraintSystem<E> + 'static, H: CircuitGLTreeHasher<E>, TR: CircuitGLTranscript<E, CircuitCompatibleCap = H::CircuitOutput>>(
    cs: &mut CS,
    proof: &AllocatedProof<E, H>,
    vk: &AllocatedVerificationKey<E, H>,
    challenges: &mut ChallengesHolder<E, CS>,
    transcript: &mut TR,
    // parameters
    verifier: &WrapperVerifier<E, CS>,
    fixed_parameters: &VerificationKeyCircuitGeometry,
    constants: &ConstantsHolder,
) -> Result<Vec<(GL, Vec<(usize, GoldilocksAsFieldWrapper<E, CS>)>)>, SynthesisError> {
    // allocate everything
    let setup_tree_cap = &vk.setup_merkle_tree_cap;
    assert_eq!(fixed_parameters.cap_size, setup_tree_cap.len());
    transcript.witness_merkle_tree_cap(cs, &setup_tree_cap)?;

    if proof.public_inputs.len() != fixed_parameters.public_inputs_locations.len() {
        panic!("Invalid number of public inputs");
    }

    let num_public_inputs = proof.public_inputs.len();
    let mut public_inputs_with_values = Vec::with_capacity(num_public_inputs);
    let mut public_input_allocated = Vec::with_capacity(num_public_inputs);

    // commit public inputs
    for ((column, row), value) in fixed_parameters.public_inputs_locations.iter().copied().zip(proof.public_inputs.iter().copied()) {
        transcript.witness_field_elements(cs, &[value])?;
        public_input_allocated.push(value);
        let value = value.into();
        public_inputs_with_values.push((column, row, value));
    }

    // commit witness
    assert_eq!(fixed_parameters.cap_size, proof.witness_oracle_cap.len());
    transcript.witness_merkle_tree_cap(cs, &proof.witness_oracle_cap)?;

    // draw challenges for stage 2
    challenges.get_beta_gamma_challenges(cs, transcript, verifier)?;

    assert_eq!(fixed_parameters.cap_size, proof.stage_2_oracle_cap.len());
    transcript.witness_merkle_tree_cap(cs, &proof.stage_2_oracle_cap)?;

    challenges.get_alpha_powers(cs, transcript, constants)?;

    // commit quotient
    assert_eq!(fixed_parameters.cap_size, proof.quotient_oracle_cap.len());
    transcript.witness_merkle_tree_cap(cs, &proof.quotient_oracle_cap)?;

    challenges.get_z_challenge(cs, transcript, fixed_parameters)?;

    // commit claimed values at z, and form our poly storage
    for set in proof.values_at_z.iter().chain(proof.values_at_z_omega.iter()).chain(proof.values_at_0.iter()) {
        transcript.witness_field_elements(cs, set)?;
    }

    // and public inputs should also go into quotient
    let mut public_input_opening_tuples: Vec<(GL, Vec<(usize, GoldilocksAsFieldWrapper<E, CS>)>)> = vec![];
    {
        let omega = domain_generator_for_size::<GL>(fixed_parameters.domain_size as u64);

        for (column, row, value) in public_inputs_with_values.into_iter() {
            let open_at = BoojumField::pow_u64(&omega, row as u64);
            let pos = public_input_opening_tuples.iter().position(|el| el.0 == open_at);
            if let Some(pos) = pos {
                public_input_opening_tuples[pos].1.push((column, value));
            } else {
                public_input_opening_tuples.push((open_at, vec![(column, value)]));
            }
        }
    }

    assert_eq!(proof.values_at_z.len(), constants.num_poly_values_at_z);
    assert_eq!(proof.values_at_z_omega.len(), constants.num_poly_values_at_z_omega);
    assert_eq!(proof.values_at_0.len(), constants.num_poly_values_at_zero);

    Ok(public_input_opening_tuples)
}
