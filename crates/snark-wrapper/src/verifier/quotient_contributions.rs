use super::*;

use std::collections::HashMap;

use crate::boojum::cs::gates::lookup_marker::LookupFormalGate;
use crate::boojum::cs::gates::lookup_marker::LookupGateMarkerFormalEvaluator;
use crate::boojum::cs::implementations::copy_permutation::non_residues_for_copy_permutation;
use crate::boojum::cs::implementations::verifier::*;
use crate::boojum::cs::traits::gate::GatePlacementStrategy;
use std::alloc::Global;

/// Run verifier at z.
/// We should check:
/// - lookup contribution
/// - specialized gates contribution
/// - general purpose gates contribution
/// - copy permutation contribution
pub(crate) fn check_quotient_contributions_in_z<E: Engine, CS: ConstraintSystem<E> + 'static, H: CircuitGLTreeHasher<E>>(
    cs: &mut CS,
    proof: &AllocatedProof<E, H>,
    challenges: &ChallengesHolder<E, CS>,
    // parameters
    verifier: &WrapperVerifier<E, CS>,
    fixed_parameters: &VerificationKeyCircuitGeometry,
    constants: &ConstantsHolder,
) -> Result<Vec<Boolean>, SynthesisError> {
    let mut validity_flags = vec![];

    let zero_ext = GoldilocksExtAsFieldWrapper::<E, CS>::zero(cs);
    let one_ext = GoldilocksExtAsFieldWrapper::<E, CS>::one(cs);

    let non_residues_for_copy_permutation = non_residues_for_copy_permutation::<GL, Global>(fixed_parameters.domain_size as usize, constants.num_variable_polys);

    let non_residues_for_copy_permutation: Vec<_> = non_residues_for_copy_permutation.into_iter().map(|el| GoldilocksAsFieldWrapper::constant(el, cs)).collect();

    let evaluations = EvaluationsHolder::from_proof(proof);

    let mut source_it = evaluations.all_values_at_z.iter();
    // witness
    let variables_polys_values: Vec<_> = (&mut source_it).take(constants.num_variable_polys).copied().collect();
    let witness_polys_values: Vec<_> = (&mut source_it).take(constants.num_witness_polys).copied().collect();
    // normal setup
    let constant_poly_values: Vec<_> = (&mut source_it).take(constants.num_constant_polys).copied().collect();
    let sigmas_values: Vec<_> = (&mut source_it).take(constants.num_copy_permutation_polys).copied().collect();
    let copy_permutation_z_at_z = *source_it.next().unwrap();
    let grand_product_intermediate_polys: Vec<_> = (&mut source_it).take(constants.num_intermediate_partial_product_relations).copied().collect();
    // lookup if exists
    let multiplicities_polys_values: Vec<_> = (&mut source_it).take(constants.num_multiplicities_polys).copied().collect();
    let lookup_witness_encoding_polys_values: Vec<_> = (&mut source_it).take(constants.num_lookup_subarguments).copied().collect();
    let multiplicities_encoding_polys_values: Vec<_> = (&mut source_it).take(constants.num_multiplicities_polys).copied().collect();
    // lookup setup
    let lookup_tables_columns: Vec<_> = (&mut source_it).take(constants.num_lookup_table_setup_polys).copied().collect();
    // quotient
    let quotient_chunks: Vec<_> = source_it.copied().collect();

    assert_eq!(quotient_chunks.len(), constants.quotient_degree);

    let mut source_it = evaluations.all_values_at_z_omega.iter();
    let copy_permutation_z_at_z_omega = *source_it.next().unwrap();

    let mut t_accumulator = GoldilocksExtAsFieldWrapper::<E, CS>::zero(cs);
    // precompute selectors at z

    let mut selectors_buffer = HashMap::new();
    for (gate_idx, evaluator) in verifier.evaluators_over_general_purpose_columns.iter().enumerate() {
        if let Some(path) = fixed_parameters.selectors_placement.output_placement(gate_idx) {
            if selectors_buffer.contains_key(&path) {
                panic!("same selector for different gates");
            }

            compute_selector_subpath_at_z(path, &mut selectors_buffer, &constant_poly_values, cs);
        } else {
            assert!(evaluator.num_quotient_terms == 0);
        }
    }

    validity_flags.extend(check_lookup_contribution(
        cs,
        &evaluations,
        challenges,
        &mut t_accumulator,
        &variables_polys_values,
        &lookup_witness_encoding_polys_values,
        &lookup_tables_columns,
        &constant_poly_values,
        &mut selectors_buffer,
        &multiplicities_encoding_polys_values,
        &multiplicities_polys_values,
        verifier,
        fixed_parameters,
        constants,
    )?);

    let constants_for_gates_over_general_purpose_columns = fixed_parameters.extra_constant_polys_for_selectors + verifier.parameters.num_constant_columns;

    let src = VerifierPolyStorage::new(variables_polys_values.clone(), witness_polys_values, constant_poly_values);

    check_specialized_gates_contribution(cs, challenges, &mut t_accumulator, &src, verifier, constants, constants_for_gates_over_general_purpose_columns)?;

    // log!("Evaluating general purpose gates");

    check_general_purpose_gates_contribution(
        cs,
        challenges,
        &mut t_accumulator,
        &src,
        &mut selectors_buffer,
        verifier,
        fixed_parameters,
        constants,
        constants_for_gates_over_general_purpose_columns,
    )?;

    // then copy_permutation algorithm

    let z_in_domain_size = challenges.z.pow_u64(fixed_parameters.domain_size as u64, cs);

    let mut vanishing_at_z = z_in_domain_size;
    vanishing_at_z.sub_assign(&one_ext, cs);

    check_copy_permutation_contribution(
        cs,
        challenges,
        &mut t_accumulator,
        &variables_polys_values,
        copy_permutation_z_at_z,
        copy_permutation_z_at_z_omega,
        &grand_product_intermediate_polys,
        &sigmas_values,
        &non_residues_for_copy_permutation,
        vanishing_at_z,
        constants,
    )?;

    let mut t_from_chunks = zero_ext;
    let mut pow = one_ext;
    for el in quotient_chunks.into_iter() {
        GoldilocksExtAsFieldWrapper::<E, CS>::mul_and_accumulate_into(&mut t_from_chunks, &el, &pow, cs);

        // let mut tmp = el;
        // tmp.mul_assign(&pow, cs);
        // t_from_chunks.add_assign(&tmp, cs);

        pow.mul_assign(&z_in_domain_size, cs);
    }

    t_from_chunks.mul_assign(&vanishing_at_z, cs);

    let t_accumulator = t_accumulator.into_coeffs_in_base();
    let t_from_chunks = t_from_chunks.into_coeffs_in_base();

    let c0_is_valid = GoldilocksField::equals(cs, &t_accumulator[0], &t_from_chunks[0])?;
    let c1_is_valid = GoldilocksField::equals(cs, &t_accumulator[1], &t_from_chunks[1])?;

    validity_flags.push(c0_is_valid);
    validity_flags.push(c1_is_valid);

    Ok(validity_flags)
}

pub(crate) fn check_lookup_contribution<E: Engine, CS: ConstraintSystem<E> + 'static>(
    cs: &mut CS,
    evaluations: &EvaluationsHolder<E, CS>,
    challenges: &ChallengesHolder<E, CS>,
    t_accumulator: &mut GoldilocksExtAsFieldWrapper<E, CS>,
    // polynomial values
    variables_polys_values: &Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    lookup_witness_encoding_polys_values: &Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    lookup_tables_columns: &Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    constant_poly_values: &Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    selectors_buffer: &mut HashMap<Vec<bool>, GoldilocksExtAsFieldWrapper<E, CS>>,
    multiplicities_encoding_polys_values: &Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    multiplicities_polys_values: &Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    // parameters
    verifier: &WrapperVerifier<E, CS>,
    fixed_parameters: &VerificationKeyCircuitGeometry,
    constants: &ConstantsHolder,
) -> Result<Vec<Boolean>, SynthesisError> {
    let mut validity_flags = vec![];

    let one_ext = GoldilocksExtAsFieldWrapper::<E, CS>::one(cs);
    let lookup_challenges = &challenges.pregenerated_challenges_for_lookup;

    // first we do the lookup
    if verifier.lookup_parameters != LookupParameters::NoLookup {
        // immediatelly do sumchecks
        let lookup_witness_encoding_polys_polys_at_0 = &evaluations.all_values_at_0[..constants.num_lookup_subarguments];
        let multiplicities_encoding_polys_at_0 = &evaluations.all_values_at_0[constants.num_lookup_subarguments..];

        let mut witness_subsum = GoldilocksExtAsFieldWrapper::<E, CS>::zero(cs);
        for a in lookup_witness_encoding_polys_polys_at_0.iter() {
            witness_subsum.add_assign(a, cs);
        }

        let mut multiplicities_subsum = GoldilocksExtAsFieldWrapper::<E, CS>::zero(cs);
        for b in multiplicities_encoding_polys_at_0.iter() {
            multiplicities_subsum.add_assign(b, cs);
        }

        let witness_subsum = witness_subsum.into_coeffs_in_base();
        let multiplicities_subsum = multiplicities_subsum.into_coeffs_in_base();

        let c0_is_valid = GoldilocksField::equals(cs, &witness_subsum[0], &multiplicities_subsum[0])?;
        let c1_is_valid = GoldilocksField::equals(cs, &witness_subsum[1], &multiplicities_subsum[1])?;

        validity_flags.push(c0_is_valid);
        validity_flags.push(c1_is_valid);

        // lookup argument related parts
        match verifier.lookup_parameters {
            LookupParameters::TableIdAsVariable { width: _, share_table_id: _ } | LookupParameters::TableIdAsConstant { width: _, share_table_id: _ } => {
                // exists by our setup
                let lookup_evaluator_id = 0;
                let selector_subpath = fixed_parameters.selectors_placement.output_placement(lookup_evaluator_id).expect("lookup gate must be placed");
                let selector = selectors_buffer.remove(&selector_subpath).expect("path must be unique and precomputed");

                let column_elements_per_subargument = verifier.lookup_parameters.columns_per_subargument() as usize;
                assert!(fixed_parameters.table_ids_column_idxes.len() == 0 || fixed_parameters.table_ids_column_idxes.len() == 1);

                // this is our lookup width, either counted by number of witness columns only, or if one includes setup
                let num_lookup_columns = column_elements_per_subargument + ((fixed_parameters.table_ids_column_idxes.len() == 1) as usize);
                assert_eq!(lookup_tables_columns.len(), num_lookup_columns);

                let capacity = column_elements_per_subargument + ((fixed_parameters.table_ids_column_idxes.len() == 1) as usize);
                let mut powers_of_gamma = Vec::with_capacity(capacity);
                let mut tmp = GoldilocksExtAsFieldWrapper::<E, CS>::one(cs);
                powers_of_gamma.push(tmp);
                for _idx in 1..capacity {
                    if _idx == 1 {
                        tmp = challenges.lookup_gamma;
                    } else {
                        tmp.mul_assign(&challenges.lookup_gamma, cs);
                    }

                    powers_of_gamma.push(tmp);
                }

                // precompute aggregation of lookup table polys
                assert_eq!(powers_of_gamma.len(), capacity);
                let mut lookup_table_columns_aggregated = challenges.lookup_beta;
                for (gamma, column) in powers_of_gamma.iter().zip(lookup_tables_columns.iter()) {
                    GoldilocksExtAsFieldWrapper::<E, CS>::mul_and_accumulate_into(&mut lookup_table_columns_aggregated, gamma, column, cs);
                }

                let mut challenges_it = lookup_challenges.iter();

                // first A polys
                let variables_columns_for_lookup = &variables_polys_values[..(column_elements_per_subargument * constants.num_lookup_subarguments)];
                assert_eq!(
                    lookup_witness_encoding_polys_values.len(),
                    variables_columns_for_lookup.chunks_exact(column_elements_per_subargument as usize).len()
                );

                for (a_poly, witness_columns) in lookup_witness_encoding_polys_values
                    .iter()
                    .zip(variables_columns_for_lookup.chunks_exact(column_elements_per_subargument as usize))
                {
                    let alpha = *challenges_it.next().expect("challenge for lookup A poly contribution");
                    let mut contribution = challenges.lookup_beta;

                    let table_id = if let Some(table_id_poly) = fixed_parameters.table_ids_column_idxes.get(0).copied() {
                        vec![constant_poly_values[table_id_poly]]
                    } else {
                        vec![]
                    };

                    for (gamma, column) in powers_of_gamma.iter().zip(witness_columns.iter().chain(table_id.iter())) {
                        GoldilocksExtAsFieldWrapper::<E, CS>::mul_and_accumulate_into(&mut contribution, gamma, column, cs);
                    }

                    // mul by A(x)
                    contribution.mul_assign(a_poly, cs);
                    // sub selector
                    contribution.sub_assign(&selector, cs);

                    // mul by power of challenge and accumulate
                    GoldilocksExtAsFieldWrapper::<E, CS>::mul_and_accumulate_into(t_accumulator, &alpha, &contribution, cs);

                    // contribution.mul_assign(&alpha, cs);
                    // t_accumulator.add_assign(&contribution, cs);
                }

                // then B polys
                assert_eq!(multiplicities_encoding_polys_values.len(), multiplicities_polys_values.len());
                for (b_poly, multiplicities_poly) in multiplicities_encoding_polys_values.iter().zip(multiplicities_polys_values.iter()) {
                    let alpha = *challenges_it.next().expect("challenge for lookup B poly contribution");
                    let mut contribution = lookup_table_columns_aggregated;
                    // mul by B(x)
                    contribution.mul_assign(b_poly, cs);
                    // sub multiplicity
                    contribution.sub_assign(multiplicities_poly, cs);

                    // mul by power of challenge and accumulate
                    GoldilocksExtAsFieldWrapper::<E, CS>::mul_and_accumulate_into(t_accumulator, &alpha, &contribution, cs);

                    // contribution.mul_assign(&alpha, cs);
                    // t_accumulator.add_assign(&contribution, cs);
                }
            }
            LookupParameters::UseSpecializedColumnsWithTableIdAsConstant {
                width: _,
                num_repetitions: _,
                share_table_id: _,
            }
            | LookupParameters::UseSpecializedColumnsWithTableIdAsVariable {
                width: _,
                num_repetitions: _,
                share_table_id: _,
            } => {
                let column_elements_per_subargument = verifier.lookup_parameters.specialized_columns_per_subargument() as usize;
                assert!(fixed_parameters.table_ids_column_idxes.len() == 0 || fixed_parameters.table_ids_column_idxes.len() == 1);

                // this is our lookup width, either counted by number of witness columns only, or if one includes setup
                let num_lookup_columns = column_elements_per_subargument + ((fixed_parameters.table_ids_column_idxes.len() == 1) as usize);
                assert_eq!(lookup_tables_columns.len(), num_lookup_columns);

                let capacity = column_elements_per_subargument + ((fixed_parameters.table_ids_column_idxes.len() == 1) as usize);
                let mut powers_of_gamma = Vec::with_capacity(capacity);
                let mut tmp = GoldilocksExtAsFieldWrapper::<E, CS>::one(cs);
                powers_of_gamma.push(tmp);
                for _idx in 1..capacity {
                    if _idx == 1 {
                        tmp = challenges.lookup_gamma;
                    } else {
                        tmp.mul_assign(&challenges.lookup_gamma, cs);
                    }

                    powers_of_gamma.push(tmp);
                }

                // precompute aggregation of lookup table polys
                assert_eq!(powers_of_gamma.len(), capacity);
                let mut lookup_table_columns_aggregated = challenges.lookup_beta;
                for (gamma, column) in powers_of_gamma.iter().zip(lookup_tables_columns.iter()) {
                    GoldilocksExtAsFieldWrapper::<E, CS>::mul_and_accumulate_into(&mut lookup_table_columns_aggregated, gamma, column, cs);
                }

                let mut challenges_it = lookup_challenges.iter();

                // first A polys
                let variables_columns_for_lookup = &variables_polys_values[verifier.parameters.num_columns_under_copy_permutation
                    ..(verifier.parameters.num_columns_under_copy_permutation + column_elements_per_subargument * constants.num_lookup_subarguments)];
                assert_eq!(
                    lookup_witness_encoding_polys_values.len(),
                    variables_columns_for_lookup.chunks_exact(column_elements_per_subargument as usize).len()
                );

                for (a_poly, witness_columns) in lookup_witness_encoding_polys_values
                    .iter()
                    .zip(variables_columns_for_lookup.chunks_exact(column_elements_per_subargument as usize))
                {
                    let alpha = *challenges_it.next().expect("challenge for lookup A poly contribution");
                    let mut contribution = challenges.lookup_beta;

                    let table_id = if let Some(table_id_poly) = fixed_parameters.table_ids_column_idxes.get(0).copied() {
                        vec![constant_poly_values[table_id_poly]]
                    } else {
                        vec![]
                    };

                    for (gamma, column) in powers_of_gamma.iter().zip(witness_columns.iter().chain(table_id.iter())) {
                        GoldilocksExtAsFieldWrapper::<E, CS>::mul_and_accumulate_into(&mut contribution, gamma, column, cs);
                    }

                    // mul by A(x)
                    contribution.mul_assign(a_poly, cs);
                    // sub numerator
                    contribution.sub_assign(&one_ext, cs);

                    // mul by power of challenge and accumulate
                    GoldilocksExtAsFieldWrapper::<E, CS>::mul_and_accumulate_into(t_accumulator, &alpha, &contribution, cs);

                    // contribution.mul_assign(&alpha, cs);
                    // t_accumulator.add_assign(&contribution, cs);
                }

                // then B polys
                assert_eq!(multiplicities_encoding_polys_values.len(), multiplicities_polys_values.len());
                for (b_poly, multiplicities_poly) in multiplicities_encoding_polys_values.iter().zip(multiplicities_polys_values.iter()) {
                    let alpha = *challenges_it.next().expect("challenge for lookup B poly contribution");
                    let mut contribution = lookup_table_columns_aggregated;
                    // mul by B(x)
                    contribution.mul_assign(b_poly, cs);
                    // sub multiplicity
                    contribution.sub_assign(multiplicities_poly, cs);

                    // mul by power of challenge and accumulate
                    GoldilocksExtAsFieldWrapper::<E, CS>::mul_and_accumulate_into(t_accumulator, &alpha, &contribution, cs);

                    // contribution.mul_assign(&alpha, cs);
                    // t_accumulator.add_assign(&contribution, cs);
                }
            }
            _ => {
                unreachable!()
            }
        }
    }

    Ok(validity_flags)
}

pub(crate) fn check_specialized_gates_contribution<E: Engine, CS: ConstraintSystem<E> + 'static>(
    cs: &mut CS,
    challenges: &ChallengesHolder<E, CS>,
    t_accumulator: &mut GoldilocksExtAsFieldWrapper<E, CS>,
    src: &VerifierPolyStorage<GL, GoldilocksExtAsFieldWrapper<E, CS>>,
    // parameters
    verifier: &WrapperVerifier<E, CS>,
    constants: &ConstantsHolder,
    constants_for_gates_over_general_purpose_columns: usize,
) -> Result<(), SynthesisError> {
    let zero_ext = GoldilocksExtAsFieldWrapper::<E, CS>::zero(cs);
    let one_ext = GoldilocksExtAsFieldWrapper::<E, CS>::one(cs);

    let specialized_evaluators_challenges = &challenges.pregenerated_challenges_for_gates_over_specialized_columns;

    let mut specialized_placement_data = vec![];
    let mut evaluation_functions = vec![];

    for (idx, (gate_type_id, evaluator)) in verifier
        .gate_type_ids_for_specialized_columns
        .iter()
        .zip(verifier.evaluators_over_specialized_columns.iter())
        .enumerate()
    {
        if gate_type_id == &std::any::TypeId::of::<LookupFormalGate>() {
            continue;
        }

        assert!(
            evaluator.total_quotient_terms_over_all_repetitions != 0,
            "evaluator {} has no contribution to quotient",
            &evaluator.debug_name,
        );
        // log!(
        //     "Will be evaluating {} over specialized columns",
        //     &evaluator.debug_name
        // );

        let num_terms = evaluator.num_quotient_terms;
        let placement_strategy = verifier.placement_strategies.get(gate_type_id).copied().expect("gate must be allowed");
        let GatePlacementStrategy::UseSpecializedColumns { num_repetitions, share_constants } = placement_strategy else {
            unreachable!();
        };

        let total_terms = num_terms * num_repetitions;

        let (initial_offset, per_repetition_offset, total_constants_available) = verifier.offsets_for_specialized_evaluators[idx];

        let placement_data = (num_repetitions, share_constants, initial_offset, per_repetition_offset, total_constants_available, total_terms);

        specialized_placement_data.push(placement_data);
        let t = &**evaluator.columnwise_satisfiability_function.as_ref().expect("must be properly configured");
        evaluation_functions.push(t);
    }

    let mut challenges_offset = 0;

    for (placement_data, evaluation_fn) in specialized_placement_data.iter().zip(evaluation_functions.iter()) {
        let (num_repetitions, share_constants, initial_offset, per_repetition_offset, _total_constants_available, total_terms) = *placement_data;

        // we self-check again
        if share_constants {
            assert_eq!(per_repetition_offset.constants_offset, 0);
        }
        let mut final_offset = initial_offset;
        for _ in 0..num_repetitions {
            final_offset.add_offset(&per_repetition_offset);
        }

        let mut dst = VerifierRelationDestination {
            accumulator: zero_ext,
            selector_value: one_ext,
            challenges: specialized_evaluators_challenges.clone(),
            current_challenge_offset: challenges_offset,
            _marker: std::marker::PhantomData,
        };

        let mut src = src.subset(
            initial_offset.variables_offset..final_offset.variables_offset,
            initial_offset.witnesses_offset..final_offset.witnesses_offset,
            (constants_for_gates_over_general_purpose_columns + initial_offset.constants_offset)..(constants_for_gates_over_general_purpose_columns + final_offset.constants_offset),
        );

        evaluation_fn.evaluate_over_columns(&mut src, &mut dst, cs);

        t_accumulator.add_assign(&dst.accumulator, cs);

        challenges_offset += total_terms;
    }

    assert_eq!(challenges_offset, constants.total_num_gate_terms_for_specialized_columns);

    Ok(())
}

pub(crate) fn check_general_purpose_gates_contribution<E: Engine, CS: ConstraintSystem<E> + 'static>(
    cs: &mut CS,
    challenges: &ChallengesHolder<E, CS>,
    t_accumulator: &mut GoldilocksExtAsFieldWrapper<E, CS>,
    src: &VerifierPolyStorage<GL, GoldilocksExtAsFieldWrapper<E, CS>>,
    selectors_buffer: &mut HashMap<Vec<bool>, GoldilocksExtAsFieldWrapper<E, CS>>,
    // parameters
    verifier: &WrapperVerifier<E, CS>,
    fixed_parameters: &VerificationKeyCircuitGeometry,
    constants: &ConstantsHolder,
    constants_for_gates_over_general_purpose_columns: usize,
) -> Result<(), SynthesisError> {
    let src = src.subset(
        0..verifier.parameters.num_columns_under_copy_permutation,
        0..verifier.parameters.num_witness_columns,
        0..constants_for_gates_over_general_purpose_columns,
    );

    let mut challenges_offset = 0;
    let zero_ext = GoldilocksExtAsFieldWrapper::<E, CS>::zero(cs);

    let general_purpose_challenges = &challenges.pregenerated_challenges_for_gates_over_general_purpose_columns;

    for (gate_idx, evaluator) in verifier.evaluators_over_general_purpose_columns.iter().enumerate() {
        if &evaluator.evaluator_type_id == &std::any::TypeId::of::<LookupGateMarkerFormalEvaluator>() {
            continue;
        }

        if evaluator.total_quotient_terms_over_all_repetitions == 0 {
            // we MAY formally have NOP gate in the set here, but we should not evaluate it.
            // NOP gate will affect selectors placement, but not the rest
            continue;
        }

        if let Some(path) = fixed_parameters.selectors_placement.output_placement(gate_idx) {
            let selector = selectors_buffer.remove(&path).expect("path must be unique and precomputed");
            let constant_placement_offset = path.len();

            let mut dst = VerifierRelationDestination {
                accumulator: zero_ext,
                selector_value: selector,
                challenges: general_purpose_challenges.clone(),
                current_challenge_offset: challenges_offset,
                _marker: std::marker::PhantomData,
            };

            let mut source = src.clone();

            let evaluation_fn = &**evaluator.rowwise_satisfiability_function.as_ref().expect("gate must be allowed");

            evaluation_fn.evaluate_over_general_purpose_columns(&mut source, &mut dst, constant_placement_offset, cs);

            t_accumulator.add_assign(&dst.accumulator, cs);
            challenges_offset += evaluator.total_quotient_terms_over_all_repetitions;
        } else {
            assert!(evaluator.num_quotient_terms == 0);
        }
    }

    assert_eq!(challenges_offset, constants.total_num_gate_terms_for_general_purpose_columns);

    Ok(())
}

pub(crate) fn check_copy_permutation_contribution<E: Engine, CS: ConstraintSystem<E> + 'static>(
    cs: &mut CS,
    challenges: &ChallengesHolder<E, CS>,
    t_accumulator: &mut GoldilocksExtAsFieldWrapper<E, CS>,
    // polynomial values
    variables_polys_values: &Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    copy_permutation_z_at_z: GoldilocksExtAsFieldWrapper<E, CS>,
    copy_permutation_z_at_z_omega: GoldilocksExtAsFieldWrapper<E, CS>,
    grand_product_intermediate_polys: &Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    sigmas_values: &Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    non_residues_for_copy_permutation: &Vec<GoldilocksAsFieldWrapper<E, CS>>,
    vanishing_at_z: GoldilocksExtAsFieldWrapper<E, CS>,
    // parameters
    constants: &ConstantsHolder,
) -> Result<(), SynthesisError> {
    let one_ext = GoldilocksExtAsFieldWrapper::<E, CS>::one(cs);

    let mut challenges_it = challenges.remaining_challenges.iter();

    {
        // (x^n - 1) / (x - 1),
        let mut z_minus_one = challenges.z;
        z_minus_one.sub_assign(&one_ext, cs);

        let mut unnormalized_l1_inverse_at_z = vanishing_at_z;
        let z_minus_one_inversed = z_minus_one.inverse(cs);
        unnormalized_l1_inverse_at_z.mul_assign(&z_minus_one_inversed, cs);

        let alpha = *challenges_it.next().expect("challenge for z(1) == 1");
        // (z(x) - 1) * l(1)
        let mut contribution = copy_permutation_z_at_z;
        contribution.sub_assign(&one_ext, cs);
        contribution.mul_assign(&unnormalized_l1_inverse_at_z, cs);

        // mul by power of challenge and accumulate
        GoldilocksExtAsFieldWrapper::<E, CS>::mul_and_accumulate_into(t_accumulator, &alpha, &contribution, cs);

        // contribution.mul_assign(&alpha, cs);
        // t_accumulator.add_assign(&contribution, cs);
    }

    // partial products

    let lhs = grand_product_intermediate_polys.iter().chain(std::iter::once(&copy_permutation_z_at_z_omega));

    let rhs = std::iter::once(&copy_permutation_z_at_z).chain(grand_product_intermediate_polys.iter());

    for (((((lhs, rhs), alpha), non_residues), variables), sigmas) in lhs
        .zip(rhs)
        .zip(&mut challenges_it)
        .zip(non_residues_for_copy_permutation.chunks(constants.quotient_degree))
        .zip(variables_polys_values.chunks(constants.quotient_degree))
        .zip(sigmas_values.chunks(constants.quotient_degree))
    {
        let mut lhs = *lhs;
        for (variable, sigma) in variables.iter().zip(sigmas.iter()) {
            // denominator is w + beta * sigma(x) + gamma
            let mut subres = *sigma;
            subres.mul_assign(&challenges.beta, cs);
            subres.add_assign(&variable, cs);
            subres.add_assign(&challenges.gamma, cs);
            lhs.mul_assign(&subres, cs);
        }

        let mut rhs = *rhs;
        let x_poly_value = challenges.z;
        for (non_res, variable) in non_residues.iter().zip(variables.iter()) {
            // numerator is w + beta * non_res * x + gamma
            let mut subres = x_poly_value;
            subres.mul_assign_by_base(cs, non_res)?;
            subres.mul_assign(&challenges.beta, cs);
            subres.add_assign(&variable, cs);
            subres.add_assign(&challenges.gamma, cs);
            rhs.mul_assign(&subres, cs);
        }

        let mut contribution = lhs;
        contribution.sub_assign(&rhs, cs);

        // mul by power of challenge and accumulate
        GoldilocksExtAsFieldWrapper::<E, CS>::mul_and_accumulate_into(t_accumulator, &alpha, &contribution, cs);

        // contribution.mul_assign(&alpha, cs);
        // t_accumulator.add_assign(&contribution, cs);
    }

    assert_eq!(challenges_it.len(), 0, "must exhaust all the challenges");

    Ok(())
}
