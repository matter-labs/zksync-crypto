use crate::cs::implementations::copy_permutation::num_intermediate_partial_product_relations;
use crate::cs::implementations::prover::compute_fri_schedule;
use crate::cs::implementations::prover::ProofConfig;
use crate::cs::implementations::verifier::Verifier;
use crate::cs::implementations::verifier::{SizeCalculator, VerificationKeyCircuitGeometry};
use crate::cs::traits::gate::GatePlacementStrategy;
use crate::field::goldilocks::{GoldilocksExt2 as GLExt2, GoldilocksField as GL};

/// This module is responsible for calculating security levels for different stages of the protocol.
/// Here the query phase is skipped, because number of queries is already calculated to achieve security level.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SecurityLevels {
    pub copy_permutation_security_bits: usize,
    pub lookup_security_bits: usize,
    pub quotient_alpha_security_bits: usize,
    pub deep_z_security_bits: usize,
    pub deep_poly_alpha_security_bits: usize,
    pub folding_round_security_bits: usize,
}

// We use Goldilocks field 2nd extension
const CHALLENGE_FIELD_SIZE_LOG2: usize = 128;

impl SecurityLevels {
    pub fn from_parameters(
        domain_size_log2: usize,
        columns_under_copy_permutation: usize,
        number_of_lookup_constraints: usize,
        number_of_quotient_terms: usize,
        lde_factor_log2: usize,
        max_constraint_degree_log2: usize,
        number_of_evaluations: usize,
        folding_factor_log2: usize,
    ) -> Self {
        let field_size_log2 = CHALLENGE_FIELD_SIZE_LOG2;
        let lde_size_log2 = domain_size_log2 + lde_factor_log2;

        Self {
            copy_permutation_security_bits: pow_bits_for_permutation(
                domain_size_log2,
                columns_under_copy_permutation,
                field_size_log2,
            ),
            lookup_security_bits: pow_bits_for_lookup(
                domain_size_log2,
                number_of_lookup_constraints,
                field_size_log2,
            ),
            quotient_alpha_security_bits: pow_bits_for_quotient(
                field_size_log2,
                number_of_quotient_terms,
                lde_factor_log2,
            ),
            deep_z_security_bits: pow_bits_for_deep_z(
                field_size_log2,
                lde_size_log2,
                max_constraint_degree_log2,
            ),
            deep_poly_alpha_security_bits: pow_bits_for_deep_poly_alpha(
                field_size_log2,
                lde_size_log2,
                number_of_evaluations,
            ),
            folding_round_security_bits: pow_bits_for_folding_round(
                field_size_log2,
                domain_size_log2,
                folding_factor_log2,
            ),
        }
    }

    pub fn from_verifier(
        verifier: &Verifier<GL, GLExt2>,
        fixed_parameters: &VerificationKeyCircuitGeometry,
        proof_config: &ProofConfig,
    ) -> Self {
        assert_eq!(verifier.parameters, fixed_parameters.parameters);
        assert_eq!(
            verifier.lookup_parameters,
            fixed_parameters.lookup_parameters
        );
        assert!(proof_config.fri_folding_schedule.is_none());
        assert_eq!(fixed_parameters.cap_size, proof_config.merkle_tree_cap_size);
        assert_eq!(fixed_parameters.fri_lde_factor, proof_config.fri_lde_factor);
        assert_eq!(
            verifier.evaluators_over_specialized_columns.len(),
            verifier.gate_type_ids_for_specialized_columns.len()
        );

        let quotient_degree = SizeCalculator::<GL, 2, GLExt2>::quotient_degree(fixed_parameters);
        let num_lookup_subarguments = SizeCalculator::<GL, 2, GLExt2>::num_sublookup_arguments(
            &verifier.parameters,
            &verifier.lookup_parameters,
        );
        let num_variable_polys = SizeCalculator::<GL, 2, GLExt2>::num_variable_polys(
            &verifier.parameters,
            verifier.total_num_variables_for_specialized_columns,
        );
        let num_witness_polys = SizeCalculator::<GL, 2, GLExt2>::num_witness_polys(
            &verifier.parameters,
            verifier.total_num_witnesses_for_specialized_columns,
        );
        let num_constant_polys = SizeCalculator::<GL, 2, GLExt2>::num_constant_polys(
            &verifier.parameters,
            fixed_parameters,
            verifier.total_num_constants_for_specialized_columns,
        );
        let num_multiplicities_polys = SizeCalculator::<GL, 2, GLExt2>::num_multipicities_polys(
            &verifier.lookup_parameters,
            fixed_parameters.total_tables_len as usize,
            fixed_parameters.domain_size,
        );
        let num_copy_permutation_polys = num_variable_polys;

        let total_num_lookup_argument_terms = num_lookup_subarguments + num_multiplicities_polys;

        let num_intermediate_partial_product_relations =
            num_intermediate_partial_product_relations(num_copy_permutation_polys, quotient_degree);

        let total_num_gate_terms_for_specialized_columns = verifier
            .evaluators_over_specialized_columns
            .iter()
            .zip(verifier.gate_type_ids_for_specialized_columns.iter())
            .map(|(evaluator, gate_type_id)| {
                let placement_strategy = verifier
                    .placement_strategies
                    .get(gate_type_id)
                    .copied()
                    .expect("gate must be allowed");
                let num_repetitions = match placement_strategy {
                    GatePlacementStrategy::UseSpecializedColumns {
                        num_repetitions, ..
                    } => num_repetitions,
                    _ => unreachable!(),
                };
                assert_eq!(evaluator.num_repetitions_on_row, num_repetitions);
                let terms_per_repetition = evaluator.num_quotient_terms;

                terms_per_repetition * num_repetitions
            })
            .sum::<usize>();

        let total_num_gate_terms_for_general_purpose_columns = verifier
            .evaluators_over_general_purpose_columns
            .iter()
            .map(|evaluator| evaluator.total_quotient_terms_over_all_repetitions)
            .sum::<usize>();

        let number_of_quotient_terms = total_num_lookup_argument_terms // and lookup is first
            + total_num_gate_terms_for_specialized_columns // then gates over specialized columns
            + total_num_gate_terms_for_general_purpose_columns // all getes terms over general purpose columns 
            + 1 // z(1) == 1 copy permutation
            + 1 // z(x * omega) = ...
            + num_intermediate_partial_product_relations; // chunking copy permutation part;

        let lde_factor_log2 = fixed_parameters.fri_lde_factor.trailing_zeros() as usize;

        // Number of openings
        let expected_lookup_polys_total = if fixed_parameters.lookup_parameters.lookup_is_allowed()
        {
            num_lookup_subarguments + // lookup witness encoding polys
            num_multiplicities_polys * 2 + // multiplicity and multiplicity encoding
            fixed_parameters.lookup_parameters.lookup_width() + // encode tables itself
            1 // encode table IDs
        } else {
            0
        };
        let num_poly_values_at_z = num_variable_polys + num_witness_polys +
            num_constant_polys + num_copy_permutation_polys +
            1 + // z_poly
            num_intermediate_partial_product_relations + // partial products in copy-permutation
            expected_lookup_polys_total + // everything from lookup
            quotient_degree; // chunks of quotient poly
        let num_public_inputs = fixed_parameters.public_inputs_locations.len();

        let total_openings_num =
            num_poly_values_at_z + 1 + total_num_lookup_argument_terms + num_public_inputs;

        let (_new_pow_bits, _num_queries, interpolation_log2s_schedule, _final_expected_degree) =
            compute_fri_schedule(
                proof_config.security_level as u32,
                proof_config.merkle_tree_cap_size,
                proof_config.pow_bits,
                fixed_parameters.fri_lde_factor.trailing_zeros(),
                fixed_parameters.domain_size.trailing_zeros(),
            );

        let max_folding_factor_log2 = interpolation_log2s_schedule
            .iter()
            .max()
            .cloned()
            .unwrap_or(0);

        Self::from_parameters(
            fixed_parameters.domain_size.trailing_zeros() as usize,
            num_copy_permutation_polys,
            num_lookup_subarguments,
            number_of_quotient_terms,
            lde_factor_log2,
            quotient_degree.next_power_of_two().trailing_zeros() as usize,
            total_openings_num,
            max_folding_factor_log2,
        )
    }
}

// TODO: check
pub fn pow_bits_for_permutation(
    domain_size_log2: usize,
    number_of_columns: usize,
    field_size_log2: usize,
) -> usize {
    let number_of_columns_log2 = number_of_columns.next_power_of_two().trailing_zeros() as usize;
    field_size_log2 - domain_size_log2 - number_of_columns_log2 - 2
}

pub fn pow_bits_for_lookup(
    domain_size_log2: usize,
    number_of_constraints: usize,
    field_size_log2: usize,
) -> usize {
    let number_of_columns_log2 =
        number_of_constraints.next_power_of_two().trailing_zeros() as usize;
    field_size_log2 - number_of_columns_log2 - domain_size_log2 - 2
}

// https://eprint.iacr.org/2022/1216.pdf
// We can bound L^+ as 4
pub fn pow_bits_for_quotient(
    challenge_field_size_log2: usize,
    powers_of_alpha: usize,
    lde_factor_log2: usize,
) -> usize {
    let powers_of_alpha_log2 = powers_of_alpha.next_power_of_two().trailing_zeros() as usize;
    challenge_field_size_log2 - powers_of_alpha_log2 - 4 - lde_factor_log2.div_ceil(2)
}

// https://eprint.iacr.org/2022/1216.pdf
// We can bound L^+ as 4
pub fn pow_bits_for_deep_z(
    challenge_field_size_log2: usize,
    lde_domain_size_log2: usize,
    max_constraint_degree_log2: usize,
) -> usize {
    challenge_field_size_log2 - lde_domain_size_log2 - max_constraint_degree_log2 - 4
}

// https://hackmd.io/@pgaf/HkKs_1ytT
pub fn pow_bits_for_deep_poly_alpha(
    challenge_field_size_log2: usize,
    domain_size_log2: usize,
    powers_of_alpha: usize,
) -> usize {
    let powers_of_alpha_log2 = powers_of_alpha.next_power_of_two().trailing_zeros() as usize;
    challenge_field_size_log2 - powers_of_alpha_log2 - domain_size_log2
}

// https://hackmd.io/@pgaf/HkKs_1ytT
pub fn pow_bits_for_folding_round(
    challenge_field_size_log2: usize,
    domain_size_log2: usize,
    folding_factor_log2: usize,
) -> usize {
    challenge_field_size_log2 - folding_factor_log2 - domain_size_log2
}
