use super::*;

use crate::traits::pow::RecursivePoWRunner;
use crate::traits::transcript::BoolsBuffer;
use crate::verifier_structs::allocated_queries::AllocatedSingleRoundQueries;

pub struct GatesInfo<E: Engine> {
    pub variable_polynomials: Vec<PolyIdentifier>,
    pub witness_polynomials: Vec<PolyIdentifier>,
    pub variable_assignments: Vec<Variable>,
    pub witness_assignments: Vec<E::Fr>,
}

pub struct MyCS<E: Engine> {
    pub num_variables: usize,
    pub num_aux_gates: usize,
    pub aux_assingments: Vec<E::Fr>,
    pub gates: Vec<GatesInfo<E>>,
}

impl<E: Engine> MyCS<E> {
    pub fn new(num_variables: usize) -> Self {
        Self {
            num_variables,
            num_aux_gates: 0,
            aux_assingments: Vec::new(),
            gates: Vec::new(),
        }
    }
}

impl<E: Engine> ConstraintSystem<E> for MyCS<E> {
    type Params = PlonkCsWidth4WithNextStepAndCustomGatesParams;

    type MainGate = SelectorOptimizedWidth4MainGateWithDNext;

    fn alloc<F>(&mut self, value: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>,
    {
        self.num_variables += 1;
        self.aux_assingments.push(value().unwrap());
        Ok(Variable::new_unchecked(Index::Aux(self.num_variables)))
    }

    fn alloc_input<F>(&mut self, value: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>,
    {
        todo!()
    }

    fn get_main_gate(&self) -> &Self::MainGate {
        todo!()
    }

    fn begin_gates_batch_for_step(&mut self) -> Result<(), SynthesisError> {
        self.num_aux_gates += 1;
        Ok(())
    }

    fn new_gate_in_batch<G: Gate<E>>(&mut self, equation: &G, coefficients_assignments: &[E::Fr], variables_assignments: &[Variable], witness_assignments: &[E::Fr]) -> Result<(), SynthesisError> {
        let gate_info = GatesInfo {
            variable_polynomials: equation.variable_polynomials().to_vec(),
            witness_polynomials: equation.witness_polynomials().to_vec(),
            variable_assignments: variables_assignments.to_vec(),
            witness_assignments: witness_assignments.to_vec(),
        };
        self.gates.push(gate_info);

        Ok(())
    }

    fn end_gates_batch_for_step(&mut self) -> Result<(), SynthesisError> {
        Ok(())
    }

    fn allocate_variables_without_gate(&mut self, variables_assignments: &[Variable], witness_assignments: &[E::Fr]) -> Result<(), SynthesisError> {
        todo!()
    }

    fn get_dummy_variable() -> Variable {
        Variable::new_unchecked(Index::Aux(0))
    }

    fn get_explicit_zero(&mut self) -> Result<Variable, SynthesisError> {
        todo!()
    }

    fn get_explicit_one(&mut self) -> Result<Variable, SynthesisError> {
        todo!()
    }

    fn add_table(&mut self, table: LookupTableApplication<E>) -> Result<std::sync::Arc<LookupTableApplication<E>>, SynthesisError> {
        todo!()
    }

    fn get_table(&self, functional_name: &str) -> Result<std::sync::Arc<LookupTableApplication<E>>, SynthesisError> {
        todo!()
    }

    fn add_multitable(&mut self, table: MultiTableApplication<E>) -> Result<(), SynthesisError> {
        todo!()
    }

    fn get_multitable(&self, functional_name: &str) -> Result<std::sync::Arc<MultiTableApplication<E>>, SynthesisError> {
        todo!()
    }

    fn apply_single_lookup_gate(&mut self, variables: &[Variable], gate: std::sync::Arc<LookupTableApplication<E>>) -> Result<(), SynthesisError> {
        todo!()
    }

    fn apply_multi_lookup_gate(&mut self, variables: &[Variable], gate: std::sync::Arc<MultiTableApplication<E>>) -> Result<(), SynthesisError> {
        todo!()
    }

    fn get_current_step_number(&self) -> usize {
        todo!()
    }

    fn get_current_aux_gate_number(&self) -> usize {
        todo!()
    }
}

pub(crate) fn verify_fri_part<
    E: Engine,
    CS: ConstraintSystem<E> + 'static,
    H: CircuitGLTreeHasher<E>,
    TR: CircuitGLTranscript<E, CircuitCompatibleCap = H::CircuitOutput>,
    POW: RecursivePoWRunner<E>,
>(
    cs: &mut CS,
    proof: &AllocatedProof<E, H>,
    vk: &AllocatedVerificationKey<E, H>,
    challenges: &mut ChallengesHolder<E, CS>,
    transcript: &mut TR,
    public_input_opening_tuples: Vec<(GL, Vec<(usize, GoldilocksAsFieldWrapper<E, CS>)>)>,
    // parameters
    verifier: &WrapperVerifier<E, CS>,
    fixed_parameters: &VerificationKeyCircuitGeometry,
    constants: &ConstantsHolder,
) -> Result<Vec<Boolean>, SynthesisError> {
    let now = std::time::Instant::now();
    let mut validity_flags = vec![];

    // get challenges
    challenges.get_challenges_for_fri_quotiening(cs, transcript, constants.total_num_challenges_for_fri_quotiening)?;
    challenges.get_fri_intermediate_challenges(cs, transcript, proof, fixed_parameters, constants)?;

    assert_eq!(constants.final_expected_degree as usize, proof.final_fri_monomials[0].len());
    assert_eq!(constants.final_expected_degree as usize, proof.final_fri_monomials[1].len());
    assert!(proof.final_fri_monomials[0].len() > 0);

    // witness monomial coeffs
    transcript.witness_field_elements(cs, &proof.final_fri_monomials[0])?;
    transcript.witness_field_elements(cs, &proof.final_fri_monomials[1])?;

    println!("  Witness field elements took {:?}", now.elapsed());
    let now = std::time::Instant::now();

    if constants.new_pow_bits != 0 {
        const SEED_BITS: usize = 256;
        // pull enough challenges from the transcript
        let num_challenges = SEED_BITS.next_multiple_of(GL::CHAR_BITS) / GL::CHAR_BITS;
        let challenges: Vec<_> = transcript.get_multiple_challenges(cs, num_challenges as usize)?;
        let (is_valid, pow_challenge_limbs) = POW::verify_from_field_elements(cs, challenges, proof.pow_challenge_le, constants.new_pow_bits)?;
        match is_valid.get_value() {
            Some(is_valid) => {
                if is_valid == false {
                    println!("PoW challenge is invalid")
                }
            }
            None => (),
        }
        transcript.witness_field_elements(cs, &pow_challenge_limbs)?;
    }
    println!("  PoW challenge limbs took {:?}", now.elapsed());
    let now = std::time::Instant::now();

    let max_needed_bits = (fixed_parameters.domain_size * fixed_parameters.fri_lde_factor as u64).trailing_zeros() as usize;

    let mut bools_buffer = BoolsBuffer {
        available: vec![],
        max_needed: max_needed_bits,
    };

    let num_bits_for_in_coset_index = max_needed_bits - fixed_parameters.fri_lde_factor.trailing_zeros() as usize;
    let base_tree_index_shift = fixed_parameters.domain_size.trailing_zeros();
    assert_eq!(num_bits_for_in_coset_index, base_tree_index_shift as usize);

    assert_eq!(constants.num_fri_repetitions, proof.queries_per_fri_repetition.len());

    let multiplicative_generator = GoldilocksAsFieldWrapper::constant(GL::multiplicative_generator(), cs);

    // precompute once, will be handy later
    let mut precomputed_powers = vec![];
    let mut precomputed_powers_inversed = vec![];
    for i in 0..=(fixed_parameters.domain_size * fixed_parameters.fri_lde_factor as u64).trailing_zeros() {
        let omega = domain_generator_for_size::<GL>(1u64 << i);
        precomputed_powers.push(omega);
        precomputed_powers_inversed.push(BoojumPrimeField::inverse(&omega).unwrap());
    }

    println!("  Precomputation took {:?}", now.elapsed());
    let now = std::time::Instant::now();

    // we also want to precompute "steps" for different interpolation degrees
    // e.g. if we interpolate 8 elements,
    // then those will be ordered as bitreverses of [0..=7], namely
    // [0, 4, 2, 6, 1, 5, 3, 7]

    // so we want to have exactly half of it, because separation by 4
    // is exactly -1, so we need [1, sqrt4(1), sqrt8(1), sqrt4(1)*sqrt8(1)]

    let mut interpolation_steps = vec![GL::ONE; 4]; // max size

    for idx in [1, 3].into_iter() {
        BoojumField::mul_assign(&mut interpolation_steps[idx], &precomputed_powers_inversed[2]);
    }
    for idx in [2, 3].into_iter() {
        BoojumField::mul_assign(&mut interpolation_steps[idx], &precomputed_powers_inversed[3]);
    }

    assert_eq!(interpolation_steps[0], GL::ONE);
    assert_eq!(BoojumField::pow_u64(&interpolation_steps[1], 4), GL::ONE);
    assert_eq!(BoojumField::pow_u64(&interpolation_steps[2], 8), GL::ONE);

    let precomputed_powers: Vec<_> = precomputed_powers.into_iter().map(|el| GoldilocksAsFieldWrapper::constant(el, cs)).collect();
    let precomputed_powers_inversed: Vec<_> = precomputed_powers_inversed.into_iter().map(|el| GoldilocksAsFieldWrapper::constant(el, cs)).collect();
    let interpolation_steps: Vec<_> = interpolation_steps.into_iter().map(|el| GoldilocksAsFieldWrapper::constant(el, cs)).collect();

    let base_oracle_depth = fixed_parameters.base_oracles_depth();

    println!("  Before queries took {:?}. Queries count: {}", now.elapsed(), proof.queries_per_fri_repetition.len());
    let now = std::time::Instant::now();

    for (idx, queries) in proof.queries_per_fri_repetition.iter().enumerate() {
        let now = std::time::Instant::now();
        let query_index_lsb_first_bits = bools_buffer.get_bits(cs, transcript, max_needed_bits)?;

        // we consider it to be some convenient for us encoding of coset + inner index.

        // Small note on indexing: when we commit to elements we use bitreversal enumeration everywhere.
        // So index `i` in the tree corresponds to the element of `omega^bitreverse(i)`.
        // This gives us natural separation of LDE cosets, such that subtrees form independent cosets,
        // and if cosets are in the form of `{1, gamma, ...} x {1, omega, ...} where gamma^lde_factor == omega,
        // then subtrees are enumerated by bitreverse powers of gamma

        // let inner_idx = &query_index_lsb_first_bits[0..num_bits_for_in_coset_index];
        // let coset_idx = &query_index_lsb_first_bits[num_bits_for_in_coset_index..];
        let base_tree_idx = query_index_lsb_first_bits.clone();

        if idx < 10 {
            println!("    - {} just before proof verification {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();
        // first verify basic inclusion proofs
        validity_flags.extend(verify_inclusion_proofs(cs, queries, proof, vk, &base_tree_idx, constants, base_oracle_depth)?);
        if idx < 10 {
            println!("    - {} inclusion proofs verification {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();

        // now perform the quotiening operation
        let zero_ext = GoldilocksExtAsFieldWrapper::<E, CS>::zero(cs);
        let mut simulated_ext_element = zero_ext;

        assert_eq!(query_index_lsb_first_bits.len(), precomputed_powers.len() - 1);

        let domain_element = pow_from_precomputations(cs, &precomputed_powers[1..], &query_index_lsb_first_bits);

        if idx < 10 {
            println!("    - {} pow from precomputation took {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();

        // we will find it handy to have power of the generator with some bits masked to be zero
        let mut power_chunks = vec![];
        let mut skip_highest_powers = 0;
        // TODO: we may save here (in circuits case especially) if we compute recursively
        for interpolation_degree_log2 in constants.fri_folding_schedule.iter() {
            let domain_element = pow_from_precomputations(
                cs,
                &precomputed_powers_inversed[(1 + interpolation_degree_log2)..],
                &query_index_lsb_first_bits[(skip_highest_powers + interpolation_degree_log2)..],
            );

            skip_highest_powers += *interpolation_degree_log2;
            power_chunks.push(domain_element);
        }
        if idx < 10 {
            println!("    - {} interpolation took {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();

        // don't forget that we are shifted
        let mut domain_element_for_quotiening = domain_element;
        domain_element_for_quotiening.mul_assign(&multiplicative_generator, cs);

        let mut domain_element_for_interpolation = domain_element_for_quotiening;

        if idx < 10 {
            println!("    - {} before quotening took {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();

        verify_quotening_operations(
            cs,
            &mut simulated_ext_element,
            queries,
            proof,
            &public_input_opening_tuples,
            domain_element_for_quotiening,
            challenges,
            verifier,
            constants,
        )?;
        if idx < 10 {
            println!("    - {} verifying quotening took {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();

        let base_coset_inverse = BoojumPrimeField::inverse(&GL::multiplicative_generator()).unwrap();

        let mut current_folded_value: GoldilocksExtAsFieldWrapper<E, CS> = simulated_ext_element;
        let mut subidx = base_tree_idx;
        let mut coset_inverse = base_coset_inverse;

        let mut expected_fri_query_len = base_oracle_depth;

        if idx < 10 {
            println!("    - {} before preparing took {:?} for attempts: {}", idx, now.elapsed(), constants.fri_folding_schedule.len());
        }
        let now = std::time::Instant::now();

        for (idx, (interpolation_degree_log2, fri_query)) in constants.fri_folding_schedule.iter().zip(queries.fri_queries.iter()).enumerate() {
            expected_fri_query_len -= *interpolation_degree_log2;
            let interpolation_degree = 1 << *interpolation_degree_log2;
            let subidx_in_leaf = &subidx[..*interpolation_degree_log2];
            let tree_idx = &subidx[*interpolation_degree_log2..];

            assert_eq!(fri_query.leaf_elements.len(), interpolation_degree * 2);

            let [c0, c1] = current_folded_value.into_coeffs_in_base();

            let c0_from_leaf = binary_select(cs, &fri_query.leaf_elements[..interpolation_degree], subidx_in_leaf)?;
            let c1_from_leaf = binary_select(cs, &fri_query.leaf_elements[interpolation_degree..], subidx_in_leaf)?;

            let c0_is_valid = GoldilocksField::equals(cs, &c0, &c0_from_leaf)?;
            let c1_is_valid = GoldilocksField::equals(cs, &c1, &c1_from_leaf)?;

            validity_flags.push(c0_is_valid);
            validity_flags.push(c1_is_valid);

            // verify query itself
            let cap = if idx == 0 { &proof.fri_base_oracle_cap } else { &proof.fri_intermediate_oracles_caps[idx - 1] };
            assert_eq!(fri_query.proof.len(), expected_fri_query_len);
            validity_flags.push(check_if_included::<E, CS, H>(cs, &fri_query.leaf_elements, &fri_query.proof, &cap, tree_idx)?);

            // interpolate
            let mut elements_to_interpolate = Vec::with_capacity(interpolation_degree);
            for (c0, c1) in fri_query.leaf_elements[..interpolation_degree].iter().zip(fri_query.leaf_elements[interpolation_degree..].iter()) {
                let as_ext = GoldilocksExtAsFieldWrapper::<E, CS>::from_coeffs_in_base([*c0, *c1]);
                elements_to_interpolate.push(as_ext);
            }

            let mut next = Vec::with_capacity(interpolation_degree / 2);
            let challenges = &challenges.fri_intermediate_challenges[idx];
            assert_eq!(challenges.len(), *interpolation_degree_log2);

            let mut base_pow = power_chunks[idx];

            for challenge in challenges.iter() {
                for (i, [a, b]) in elements_to_interpolate.array_chunks::<2>().enumerate() {
                    let mut result = *a;
                    result.add_assign(b, cs);

                    let mut diff = *a;
                    diff.sub_assign(&b, cs);
                    diff.mul_assign(&challenge, cs);
                    // divide by corresponding power
                    let mut pow = base_pow;
                    pow.mul_assign(&interpolation_steps[i], cs);
                    let coset_inverse = GoldilocksAsFieldWrapper::constant(coset_inverse, cs);
                    pow.mul_assign(&coset_inverse, cs);

                    GoldilocksExtAsFieldWrapper::<E, CS>::mul_by_base_and_accumulate_into(&mut result, &pow, &diff, cs)?;

                    // diff.mul_assign_by_base(&pow, cs);
                    // result.add_assign(&diff, cs);

                    next.push(result);
                }

                std::mem::swap(&mut next, &mut elements_to_interpolate);
                next.clear();
                base_pow.square(cs);
                BoojumField::square(&mut coset_inverse);
            }

            for _ in 0..*interpolation_degree_log2 {
                domain_element_for_interpolation.square(cs);
            }

            // recompute the index
            subidx = tree_idx.to_vec();
            current_folded_value = elements_to_interpolate[0];
        }
        if idx < 10 {
            println!("    - {} fri folding took {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();

        // and we should evaluate monomial form and compare

        let mut result_from_monomial = zero_ext;
        // horner rule
        for (c0, c1) in proof.final_fri_monomials[0].iter().zip(proof.final_fri_monomials[1].iter()).rev() {
            let coeff = GoldilocksExtAsFieldWrapper::<E, CS>::from_coeffs_in_base([*c0, *c1]);

            // result_from_monomial = result_from_monomial * z + coeff

            let mut tmp = coeff;
            GoldilocksExtAsFieldWrapper::<E, CS>::mul_by_base_and_accumulate_into(&mut tmp, &domain_element_for_interpolation, &result_from_monomial, cs)?;

            result_from_monomial = tmp;

            // result_from_monomial.mul_assign_by_base(&domain_element_for_interpolation, cs);
            // result_from_monomial.add_assign(&coeff, cs);
        }

        if idx < 10 {
            println!("    - {} horner took {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();

        let result_from_monomial = result_from_monomial.into_coeffs_in_base();
        let current_folded_value = current_folded_value.into_coeffs_in_base();

        let c0_is_valid = GoldilocksField::equals(cs, &result_from_monomial[0], &current_folded_value[0])?;
        let c1_is_valid = GoldilocksField::equals(cs, &result_from_monomial[1], &current_folded_value[1])?;

        validity_flags.push(c0_is_valid);
        validity_flags.push(c1_is_valid);
        if idx < 10 {
            println!("    - {} last part took {:?}", idx, now.elapsed());
        }
    }
    println!("  Queries took {:?}", now.elapsed());

    Ok(validity_flags)
}

fn verify_inclusion_proofs<E: Engine, CS: ConstraintSystem<E> + 'static, H: CircuitGLTreeHasher<E>>(
    cs: &mut CS,
    queries: &AllocatedSingleRoundQueries<E, H>,
    proof: &AllocatedProof<E, H>,
    vk: &AllocatedVerificationKey<E, H>,
    base_tree_idx: &Vec<Boolean>,
    constants: &ConstantsHolder,
    base_oracle_depth: usize,
) -> Result<Vec<Boolean>, SynthesisError> {
    let mut validity_flags = Vec::new();

    assert_eq!(constants.witness_leaf_size, queries.witness_query.leaf_elements.len());
    assert_eq!(base_oracle_depth, queries.witness_query.proof.len());

    let mut mycs = MyCS::new(0);
    check_if_included::<E, _, H>(&mut mycs, &queries.witness_query.leaf_elements, &queries.witness_query.proof, &proof.witness_oracle_cap, &base_tree_idx).unwrap();

    println!("MyCS Stats: num variables: {} num gates: {} ", mycs.num_variables, mycs.num_aux_gates);

    validity_flags.push(check_if_included::<E, CS, H>(
        cs,
        &queries.witness_query.leaf_elements,
        &queries.witness_query.proof,
        &proof.witness_oracle_cap,
        &base_tree_idx,
    )?);

    assert_eq!(constants.stage_2_leaf_size, queries.stage_2_query.leaf_elements.len());
    assert_eq!(base_oracle_depth, queries.stage_2_query.proof.len());
    validity_flags.push(check_if_included::<E, CS, H>(
        cs,
        &queries.stage_2_query.leaf_elements,
        &queries.stage_2_query.proof,
        &proof.stage_2_oracle_cap,
        &base_tree_idx,
    )?);

    assert_eq!(constants.quotient_leaf_size, queries.quotient_query.leaf_elements.len());
    assert_eq!(base_oracle_depth, queries.quotient_query.proof.len());
    validity_flags.push(check_if_included::<E, CS, H>(
        cs,
        &queries.quotient_query.leaf_elements,
        &queries.quotient_query.proof,
        &proof.quotient_oracle_cap,
        &base_tree_idx,
    )?);

    assert_eq!(constants.setup_leaf_size, queries.setup_query.leaf_elements.len());
    assert_eq!(base_oracle_depth, queries.setup_query.proof.len());
    validity_flags.push(check_if_included::<E, CS, H>(
        cs,
        &queries.setup_query.leaf_elements,
        &queries.setup_query.proof,
        &vk.setup_merkle_tree_cap,
        &base_tree_idx,
    )?);

    Ok(validity_flags)
}

fn verify_quotening_operations<E: Engine, CS: ConstraintSystem<E> + 'static, H: CircuitGLTreeHasher<E>>(
    cs: &mut CS,
    simulated_ext_element: &mut GoldilocksExtAsFieldWrapper<E, CS>,
    queries: &AllocatedSingleRoundQueries<E, H>,
    proof: &AllocatedProof<E, H>,
    public_input_opening_tuples: &Vec<(GL, Vec<(usize, GoldilocksAsFieldWrapper<E, CS>)>)>,
    domain_element_for_quotiening: GoldilocksAsFieldWrapper<E, CS>,
    challenges: &ChallengesHolder<E, CS>,
    verifier: &WrapperVerifier<E, CS>,
    constants: &ConstantsHolder,
) -> Result<(), SynthesisError> {
    let zero_num = GoldilocksField::zero();
    let zero_base = GoldilocksAsFieldWrapper::<E, CS>::zero(cs);
    let zero_ext = GoldilocksExtAsFieldWrapper::<E, CS>::zero(cs);

    let mut challenge_offset = 0;

    let z_polys_offset = 0;
    let intermediate_polys_offset = 2;
    let lookup_witness_encoding_polys_offset = intermediate_polys_offset + constants.num_intermediate_partial_product_relations * 2;
    let lookup_multiplicities_encoding_polys_offset = lookup_witness_encoding_polys_offset + constants.num_lookup_subarguments * 2;
    let copy_permutation_polys_offset = 0;
    let constants_offset = 0 + constants.num_copy_permutation_polys;
    let lookup_tables_values_offset = 0 + constants.num_copy_permutation_polys + constants.num_constant_polys;
    let variables_offset = 0;
    let witness_columns_offset = constants.num_variable_polys;
    let lookup_multiplicities_offset = witness_columns_offset + constants.num_witness_polys;

    let evaluations = EvaluationsHolder::from_proof(proof);

    {
        let z = challenges.z;
        let z_omega = challenges.z_omega;

        let cast_from_base = move |el: &[GoldilocksField<E>]| {
            el.iter()
                .map(|el| GoldilocksExtAsFieldWrapper::<E, CS>::from_coeffs_in_base([*el, GoldilocksField::zero()]))
                .collect::<Vec<_>>()
        };

        let cast_from_extension = move |el: &[GoldilocksField<E>]| {
            assert_eq!(el.len() % 2, 0);

            el.array_chunks::<2>()
                .map(|[c0, c1]| GoldilocksExtAsFieldWrapper::<E, CS>::from_coeffs_in_base([*c0, *c1]))
                .collect::<Vec<_>>()
        };

        let mut sources = vec![];
        // witness
        sources.extend(cast_from_base(
            &queries.witness_query.leaf_elements[variables_offset..(variables_offset + constants.num_variable_polys)],
        ));
        sources.extend(cast_from_base(
            &queries.witness_query.leaf_elements[witness_columns_offset..(witness_columns_offset + constants.num_witness_polys)],
        ));
        // normal setup
        sources.extend(cast_from_base(&queries.setup_query.leaf_elements[constants_offset..(constants_offset + constants.num_constant_polys)]));
        sources.extend(cast_from_base(
            &queries.setup_query.leaf_elements[copy_permutation_polys_offset..(copy_permutation_polys_offset + constants.num_copy_permutation_polys)],
        ));
        // copy-permutation
        sources.extend(cast_from_extension(&queries.stage_2_query.leaf_elements[z_polys_offset..intermediate_polys_offset]));
        sources.extend(cast_from_extension(
            &queries.stage_2_query.leaf_elements[intermediate_polys_offset..lookup_witness_encoding_polys_offset],
        ));
        // lookup if exists
        sources.extend(cast_from_base(
            &queries.witness_query.leaf_elements[lookup_multiplicities_offset..(lookup_multiplicities_offset + constants.num_multiplicities_polys)],
        ));
        sources.extend(cast_from_extension(
            &queries.stage_2_query.leaf_elements[lookup_witness_encoding_polys_offset..lookup_multiplicities_encoding_polys_offset],
        ));
        sources.extend(cast_from_extension(&queries.stage_2_query.leaf_elements[lookup_multiplicities_encoding_polys_offset..]));
        // lookup setup
        if verifier.lookup_parameters.lookup_is_allowed() {
            let num_lookup_setups = verifier.lookup_parameters.lookup_width() + 1;
            sources.extend(cast_from_base(
                &queries.setup_query.leaf_elements[lookup_tables_values_offset..(lookup_tables_values_offset + num_lookup_setups)],
            ));
        }
        // quotient
        sources.extend(cast_from_extension(&queries.quotient_query.leaf_elements));

        assert_eq!(sources.len(), evaluations.all_values_at_z.len());
        // log!("Making quotiening at Z");
        quotening_operation(
            cs,
            simulated_ext_element,
            &sources,
            &evaluations.all_values_at_z,
            domain_element_for_quotiening,
            z,
            &challenges.challenges_for_fri_quotiening[challenge_offset..(challenge_offset + sources.len())],
        );
        challenge_offset += sources.len();

        // now z*omega
        let mut sources = vec![];
        sources.extend(cast_from_extension(&queries.stage_2_query.leaf_elements[z_polys_offset..intermediate_polys_offset]));

        assert_eq!(sources.len(), evaluations.all_values_at_z_omega.len());
        // log!("Making quotiening at Z*omega");
        quotening_operation(
            cs,
            simulated_ext_element,
            &sources,
            &evaluations.all_values_at_z_omega,
            domain_element_for_quotiening,
            z_omega,
            &challenges.challenges_for_fri_quotiening[challenge_offset..(challenge_offset + sources.len())],
        );

        challenge_offset += sources.len();
        // now at 0 if lookup is needed
        if verifier.lookup_parameters.lookup_is_allowed() {
            let mut sources = vec![];
            // witness encoding
            sources.extend(cast_from_extension(
                &queries.stage_2_query.leaf_elements[lookup_witness_encoding_polys_offset..lookup_multiplicities_encoding_polys_offset],
            ));
            // multiplicities encoding
            sources.extend(cast_from_extension(&queries.stage_2_query.leaf_elements[lookup_multiplicities_encoding_polys_offset..]));

            assert_eq!(sources.len(), evaluations.all_values_at_0.len());
            // log!("Making quotiening at 0 for lookups sumchecks");
            quotening_operation(
                cs,
                simulated_ext_element,
                &sources,
                &evaluations.all_values_at_0,
                domain_element_for_quotiening,
                zero_ext,
                &challenges.challenges_for_fri_quotiening[challenge_offset..(challenge_offset + sources.len())],
            );

            challenge_offset += sources.len();
        }
    }

    // and public inputs
    for (open_at, set) in public_input_opening_tuples.iter() {
        let mut sources = Vec::with_capacity(set.len());
        let mut values = Vec::with_capacity(set.len());
        for (column, expected_value) in set.into_iter() {
            let c0 = queries.witness_query.leaf_elements[*column];
            let el = GoldilocksExtAsFieldWrapper::<E, CS>::from_coeffs_in_base([c0, zero_num]);
            sources.push(el);

            let value = GoldilocksExtAsFieldWrapper::<E, CS>::from_wrapper_coeffs_in_base([*expected_value, zero_base]);
            values.push(value);
        }
        let num_challenges_required = sources.len();
        assert_eq!(values.len(), num_challenges_required);

        // log!("Making quotiening at {} for public inputs", open_at);

        let open_at = GoldilocksAsFieldWrapper::constant(*open_at, cs);
        let open_at = GoldilocksExtAsFieldWrapper::<E, CS>::from_wrapper_coeffs_in_base([open_at, zero_base]);

        quotening_operation(
            cs,
            simulated_ext_element,
            &sources,
            &values,
            domain_element_for_quotiening,
            open_at,
            &challenges.challenges_for_fri_quotiening[challenge_offset..(challenge_offset + sources.len())],
        );

        challenge_offset += num_challenges_required;
    }

    assert_eq!(challenge_offset, challenges.challenges_for_fri_quotiening.len());

    Ok(())
}

fn check_if_included<E: Engine, CS: ConstraintSystem<E>, H: CircuitGLTreeHasher<E>>(
    cs: &mut CS,
    leaf_elements: &Vec<GoldilocksField<E>>,
    proof: &Vec<H::CircuitOutput>,
    tree_cap: &Vec<H::CircuitOutput>,
    path: &[Boolean],
) -> Result<Boolean, SynthesisError> {
    let leaf_hash = <H as CircuitGLTreeHasher<E>>::hash_into_leaf(cs, leaf_elements)?;
    verify_proof_over_cap::<E, H, CS>(cs, proof, tree_cap, &leaf_hash, path)
}

pub fn verify_proof_over_cap<E: Engine, H: CircuitGLTreeHasher<E>, CS: ConstraintSystem<E>>(
    cs: &mut CS,
    proof: &[H::CircuitOutput],
    cap: &[H::CircuitOutput],
    leaf_hash: &H::CircuitOutput,
    path: &[Boolean],
) -> Result<Boolean, SynthesisError> {
    assert!(path.len() >= proof.len());

    let mut current = leaf_hash.clone();
    let path_bits = &path[..proof.len()];
    let cap_bits = &path[proof.len()..];

    for (proof_el, path_bit) in proof.iter().zip(path_bits.iter()) {
        let (left, right) = H::swap_nodes(cs, *path_bit, &current, &proof_el, 0)?;
        current = <H as CircuitGLTreeHasher<E>>::hash_into_node(cs, &left, &right, 0)?;
    }

    let selected_cap_el = H::select_cap_node(cs, cap_bits, cap)?;

    H::compare_output(cs, &current, &selected_cap_el)
    //Ok(Boolean::Constant(true))
}

fn pow_from_precomputations<E: Engine, CS: ConstraintSystem<E> + 'static>(cs: &mut CS, bases: &[GoldilocksAsFieldWrapper<E, CS>], bits: &[Boolean]) -> GoldilocksAsFieldWrapper<E, CS> {
    let mut result = GoldilocksAsFieldWrapper::<E, CS>::one(cs);

    for (base, bit) in bases.iter().zip(bits.iter()) {
        let mut tmp = result;
        tmp.mul_assign(base, cs);
        result = GoldilocksAsFieldWrapper::conditionally_select(cs, *bit, &tmp, &result);
    }

    result
}

fn quotening_operation<E: Engine, CS: ConstraintSystem<E> + 'static>(
    cs: &mut CS,
    dst: &mut GoldilocksExtAsFieldWrapper<E, CS>,
    polynomial_values: &Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    values_at: &Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    domain_element: GoldilocksAsFieldWrapper<E, CS>,
    at: GoldilocksExtAsFieldWrapper<E, CS>,
    challenges: &[GoldilocksExtAsFieldWrapper<E, CS>],
) {
    // we precompute challenges outside to avoid any manual extension ops here
    assert_eq!(polynomial_values.len(), values_at.len());
    assert_eq!(polynomial_values.len(), challenges.len());

    let zero_base = GoldilocksAsFieldWrapper::zero(cs);

    let mut denom = GoldilocksExtAsFieldWrapper::<E, CS>::from_wrapper_coeffs_in_base([domain_element, zero_base]);
    denom.sub_assign(&at, cs);
    denom = denom.inverse(cs);

    let mut acc = GoldilocksExtAsFieldWrapper::<E, CS>::zero(cs);

    for ((poly_value, value_at), challenge) in polynomial_values.iter().zip(values_at.iter()).zip(challenges.iter()) {
        // (f(x) - f(z))/(x - z)
        let mut tmp = *poly_value;
        tmp.sub_assign(&value_at, cs);

        GoldilocksExtAsFieldWrapper::<E, CS>::mul_and_accumulate_into(&mut acc, &tmp, &challenge, cs);

        // let mut as_ext = *challenge;
        // as_ext.mul_assign(&tmp, cs);
        // acc.add_assign(&as_ext, cs);
    }

    GoldilocksExtAsFieldWrapper::<E, CS>::mul_and_accumulate_into(dst, &acc, &denom, cs);

    // acc.mul_assign(&denom, cs);
    // dst.add_assign(&acc, cs);
}
