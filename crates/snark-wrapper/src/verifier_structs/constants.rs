use super::*;

use crate::boojum::cs::implementations::verifier::{SizeCalculator, VerificationKeyCircuitGeometry};
use crate::boojum::field::goldilocks::{GoldilocksExt2 as GLExt2, GoldilocksField as GL};

#[derive(Clone, Default, Debug)]
pub(crate) struct ConstantsHolder {
    // quotient parameters
    pub(crate) quotient_degree: usize,
    pub(crate) num_lookup_subarguments: usize,
    pub(crate) num_variable_polys: usize,
    pub(crate) num_witness_polys: usize,
    pub(crate) num_constant_polys: usize,
    pub(crate) num_multiplicities_polys: usize,
    pub(crate) num_copy_permutation_polys: usize,
    pub(crate) num_lookup_table_setup_polys: usize,
    pub(crate) num_intermediate_partial_product_relations: usize,
    pub(crate) total_num_gate_terms_for_specialized_columns: usize,
    pub(crate) total_num_gate_terms_for_general_purpose_columns: usize,
    pub(crate) total_num_lookup_argument_terms: usize,
    pub(crate) total_num_terms: usize,

    // commitments parameters
    pub(crate) witness_leaf_size: usize,
    pub(crate) stage_2_leaf_size: usize,
    pub(crate) quotient_leaf_size: usize,
    pub(crate) setup_leaf_size: usize,

    // opening parameters
    pub(crate) num_poly_values_at_z: usize,
    pub(crate) num_poly_values_at_z_omega: usize,
    pub(crate) num_poly_values_at_zero: usize,
    pub(crate) num_public_inputs: usize,

    // fri parameters
    pub(crate) new_pow_bits: usize,
    pub(crate) num_fri_repetitions: usize,
    pub(crate) fri_folding_schedule: Vec<usize>,
    pub(crate) final_expected_degree: usize,
    pub(crate) total_num_challenges_for_fri_quotiening: usize,
}

impl ConstantsHolder {
    pub fn generate<E: Engine, CS: ConstraintSystem<E>>(proof_config: &ProofConfig, verifier: &WrapperVerifier<E, CS>, fixed_parameters: &VerificationKeyCircuitGeometry) -> Self {
        assert_eq!(verifier.parameters, fixed_parameters.parameters);
        assert_eq!(verifier.lookup_parameters, fixed_parameters.lookup_parameters);
        assert!(proof_config.fri_folding_schedule.is_none());
        assert_eq!(fixed_parameters.cap_size, proof_config.merkle_tree_cap_size);
        assert_eq!(fixed_parameters.fri_lde_factor, proof_config.fri_lde_factor,);

        let mut result = Self::default();

        result.quotient_degree = SizeCalculator::<GL, 2, GLExt2>::quotient_degree(fixed_parameters);
        result.num_lookup_subarguments = SizeCalculator::<GL, 2, GLExt2>::num_sublookup_arguments(&verifier.parameters, &verifier.lookup_parameters);
        result.num_variable_polys = SizeCalculator::<GL, 2, GLExt2>::num_variable_polys(&verifier.parameters, verifier.total_num_variables_for_specialized_columns);
        result.num_witness_polys = SizeCalculator::<GL, 2, GLExt2>::num_witness_polys(&verifier.parameters, verifier.total_num_witnesses_for_specialized_columns);
        result.num_constant_polys = SizeCalculator::<GL, 2, GLExt2>::num_constant_polys(&verifier.parameters, fixed_parameters, verifier.total_num_constants_for_specialized_columns);
        result.num_multiplicities_polys =
            SizeCalculator::<GL, 2, GLExt2>::num_multipicities_polys(&verifier.lookup_parameters, fixed_parameters.total_tables_len as usize, fixed_parameters.domain_size);
        result.num_copy_permutation_polys = result.num_variable_polys;

        result.num_lookup_table_setup_polys = SizeCalculator::<GL, 2, GLExt2>::num_lookup_table_setup_polys(&verifier.lookup_parameters);

        result.witness_leaf_size = SizeCalculator::<GL, 2, GLExt2>::witness_leaf_size(
            &verifier.parameters,
            &verifier.lookup_parameters,
            fixed_parameters,
            verifier.total_num_variables_for_specialized_columns,
            verifier.total_num_witnesses_for_specialized_columns,
        );
        result.stage_2_leaf_size = SizeCalculator::<GL, 2, GLExt2>::stage_2_leaf_size(
            &verifier.parameters,
            &verifier.lookup_parameters,
            fixed_parameters,
            verifier.total_num_variables_for_specialized_columns,
        );
        result.quotient_leaf_size = SizeCalculator::<GL, 2, GLExt2>::quotient_leaf_size(fixed_parameters);
        result.setup_leaf_size = SizeCalculator::<GL, 2, GLExt2>::setup_leaf_size(
            &verifier.parameters,
            &verifier.lookup_parameters,
            fixed_parameters,
            verifier.total_num_variables_for_specialized_columns,
            verifier.total_num_constants_for_specialized_columns,
        );

        result.total_num_lookup_argument_terms = result.num_lookup_subarguments + result.num_multiplicities_polys;

        use crate::boojum::cs::implementations::copy_permutation::num_intermediate_partial_product_relations;
        result.num_intermediate_partial_product_relations = num_intermediate_partial_product_relations(result.num_copy_permutation_polys, result.quotient_degree);

        result.compute_num_gate_terms(verifier);
        result.compute_total_num_terms();

        result.compute_num_poly_values_at_z(fixed_parameters);
        result.compute_num_poly_values_at_z_omega();
        result.compute_num_poly_values_at_zero();
        result.num_public_inputs = fixed_parameters.public_inputs_locations.len();

        result.compute_fri_parameters(fixed_parameters, proof_config);

        result.compute_total_num_challenges_for_fri_quotiening(verifier);

        result
    }

    fn compute_num_gate_terms<E: Engine, CS: ConstraintSystem<E>>(&mut self, verifier: &WrapperVerifier<E, CS>) {
        assert_eq!(verifier.evaluators_over_specialized_columns.len(), verifier.gate_type_ids_for_specialized_columns.len());

        self.total_num_gate_terms_for_specialized_columns = verifier
            .evaluators_over_specialized_columns
            .iter()
            .zip(verifier.gate_type_ids_for_specialized_columns.iter())
            .map(|(evaluator, gate_type_id)| {
                let placement_strategy = verifier.placement_strategies.get(gate_type_id).copied().expect("gate must be allowed");
                let num_repetitions = match placement_strategy {
                    GatePlacementStrategy::UseSpecializedColumns { num_repetitions, .. } => num_repetitions,
                    _ => unreachable!(),
                };
                assert_eq!(evaluator.num_repetitions_on_row, num_repetitions);
                let terms_per_repetition = evaluator.num_quotient_terms;

                terms_per_repetition * num_repetitions
            })
            .sum();

        self.total_num_gate_terms_for_general_purpose_columns = verifier
            .evaluators_over_general_purpose_columns
            .iter()
            .map(|evaluator| evaluator.total_quotient_terms_over_all_repetitions)
            .sum();
    }

    fn compute_total_num_terms(&mut self) {
        self.total_num_terms = self.total_num_lookup_argument_terms // and lookup is first
        + self.total_num_gate_terms_for_specialized_columns // then gates over specialized columns
        + self.total_num_gate_terms_for_general_purpose_columns // all getes terms over general purpose columns 
        + 1 // z(1) == 1 copy permutation
        + 1 // z(x * omega) = ...
        + self.num_intermediate_partial_product_relations; // chunking copy permutation part;
    }

    fn compute_num_poly_values_at_z(&mut self, fixed_parameters: &VerificationKeyCircuitGeometry) {
        let expected_lookup_polys_total = if fixed_parameters.lookup_parameters.lookup_is_allowed() {
            self.num_lookup_subarguments + // lookup witness encoding polys
            self.num_multiplicities_polys * 2 + // multiplicity and multiplicity encoding
            fixed_parameters.lookup_parameters.lookup_width() + // encode tables itself
            1 // encode table IDs
        } else {
            0
        };

        self.num_poly_values_at_z = self.num_variable_polys + self.num_witness_polys +
            self.num_constant_polys + self.num_copy_permutation_polys +
            1 + // z_poly
            self.num_intermediate_partial_product_relations + // partial products in copy-permutation
            expected_lookup_polys_total + // everything from lookup
            self.quotient_degree; // chunks of quotient poly
    }

    fn compute_num_poly_values_at_z_omega(&mut self) {
        self.num_poly_values_at_z_omega = 1;
    }

    fn compute_num_poly_values_at_zero(&mut self) {
        self.num_poly_values_at_zero = self.num_lookup_subarguments + self.num_multiplicities_polys;
    }

    fn compute_total_num_challenges_for_fri_quotiening<E: Engine, CS: ConstraintSystem<E>>(&mut self, verifier: &WrapperVerifier<E, CS>) {
        let expected_lookup_polys_total = if verifier.lookup_parameters.lookup_is_allowed() {
            self.num_lookup_subarguments + // lookup witness encoding polys
            self.num_multiplicities_polys * 2 + // multiplicity and multiplicity encoding
            verifier.lookup_parameters.lookup_width() + // encode tables itself
            1 // encode table IDs
        } else {
            0
        };

        let num_poly_values_at_z = self.num_variable_polys + self.num_witness_polys +
        self.num_constant_polys + self.num_copy_permutation_polys +
            1 + // z_poly
            self.num_intermediate_partial_product_relations + // partial products in copy-permutation
            expected_lookup_polys_total + // everything from lookup
            self.quotient_degree; // chunks of quotient poly

        let mut total_num_challenges = 0;
        total_num_challenges += num_poly_values_at_z;
        total_num_challenges += 1;
        total_num_challenges += self.total_num_lookup_argument_terms;
        total_num_challenges += self.num_public_inputs;

        self.total_num_challenges_for_fri_quotiening = total_num_challenges;
    }

    fn compute_fri_parameters(&mut self, fixed_parameters: &VerificationKeyCircuitGeometry, proof_config: &ProofConfig) {
        let (
            new_pow_bits,                 // updated POW bits if needed
            num_queries,                  // num queries
            interpolation_log2s_schedule, // folding schedule
            final_expected_degree,
        ) = crate::boojum::cs::implementations::prover::compute_fri_schedule(
            proof_config.security_level as u32,
            proof_config.merkle_tree_cap_size,
            proof_config.pow_bits,
            fixed_parameters.fri_lde_factor.trailing_zeros(),
            fixed_parameters.domain_size.trailing_zeros(),
        );

        let mut expected_degree = fixed_parameters.domain_size;

        for interpolation_degree_log2 in interpolation_log2s_schedule.iter() {
            expected_degree >>= interpolation_degree_log2;
        }

        assert_eq!(final_expected_degree, expected_degree as usize);

        self.new_pow_bits = new_pow_bits as usize;
        self.num_fri_repetitions = num_queries;
        self.fri_folding_schedule = interpolation_log2s_schedule;
        self.final_expected_degree = final_expected_degree;
    }
}
