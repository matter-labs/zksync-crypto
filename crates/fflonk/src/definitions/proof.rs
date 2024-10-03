use bellman::PrimeFieldRepr;

use super::*;

#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct FflonkProof<E: Engine, C: Circuit<E>> {
    pub n: usize,
    pub inputs: Vec<E::Fr>,
    pub commitments: Vec<E::G1Affine>,
    pub evaluations: Vec<E::Fr>,
    pub lagrange_basis_inverses: Vec<E::Fr>,
    _c: std::marker::PhantomData<C>,
}

pub fn fe_slice_into_be_byte_array<F: PrimeField>(values: &[F]) -> Vec<[u8; 32]> {
    let mut out = vec![];

    for value in values {
        let mut buf = [0; 32];
        value.into_repr().write_be(&mut buf[..]).unwrap();
        out.push(buf);
    }

    out
}

impl<E: Engine, C: Circuit<E>> FflonkProof<E, C> {
    pub fn empty() -> Self {
        Self {
            n: 0,
            inputs: vec![],
            commitments: vec![],
            evaluations: vec![],
            lagrange_basis_inverses: vec![],
            _c: std::marker::PhantomData,
        }
    }
    pub fn as_evm_format(&self) -> (Vec<[u8; 32]>, Vec<[u8; 32]>) {
        let Self {
            inputs,
            commitments,
            evaluations,
            lagrange_basis_inverses,
            ..
        } = self;

        let serialized_inputs = fe_slice_into_be_byte_array(inputs);

        let mut flattened_commitments = vec![];
        for p in commitments.into_iter() {
            let (x, y) = p.as_xy();
            flattened_commitments.push(x.clone());
            flattened_commitments.push(y.clone());
        }
        let mut serialized_proof = fe_slice_into_be_byte_array(&flattened_commitments);

        serialized_proof.extend(fe_slice_into_be_byte_array(evaluations));
        serialized_proof.extend(fe_slice_into_be_byte_array(lagrange_basis_inverses));

        (serialized_inputs, serialized_proof)
    }
}

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct EvaluationOffsets {
    pub trace: TraceAndGateEvaluationOffsets,
    pub copy_permutation: CopyPermutationEvaluationOffsets,
    pub lookup: Option<LookupEvaluationOffsets>,
    pub setup: SetupEvaluationOffsets,
    pub has_custom_gate: bool,
    pub c0: usize,
    pub c0_shifted: usize,
    pub c1: usize,
    pub c1_shifted: usize,
    pub c2: usize,
    pub c2_shifted: usize,
}
impl EvaluationOffsets {
    pub fn from_setup<E: Engine, C: Circuit<E>>(setup: &FflonkSetup<E, C>, crs: &Crs<E, CrsForMonomialForm>) -> Self {
        let vk = FflonkVerificationKey::from_setup(setup, crs).unwrap();

        Self::from_vk(&vk)
    }

    pub fn from_vk<E: Engine, C: Circuit<E>>(vk: &FflonkVerificationKey<E, C>) -> Self {
        let (num_setup_polys, num_first_round_polys, num_second_round_polys, max_num_polys) = num_system_polys_from_vk(vk);
        let gates = sorted_gates_from_circuit_definitions::<_, C>();
        let has_custom_gate = gates.len() > 1;
        let main_gate = gates[0].clone();
        let has_lookup = vk.total_lookup_entries_length > 0;

        let some_lookup = match has_lookup {
            true => Some(true),
            false => None,
        };
        let some_custom_gate = match has_custom_gate {
            true => Some(true),
            false => None,
        };
        let num_state_polys = vk.num_state_polys;
        assert_eq!(vk.num_witness_polys, 0);
        let num_gate_setups = gates[0].setup_polynomials().len();
        let num_lookup_table_cols = if has_lookup { 4 } else { 0 };
        let num_gate_selectors = some_custom_gate.map(|_| 2).unwrap_or(0);
        let num_permutations = num_state_polys;
        let num_all_polys = num_setup_polys + num_first_round_polys + num_second_round_polys;
        let main_gate_name = main_gate.name();
        let (a_term, b_term, c_term, d_term, q_ab_term, q_ac_term, q_const_term, q_dnext_term) = if main_gate_name == STD_MAIN_GATE_NAME {
            (0, 1, 2, 0, 3, 0, 4, 0)
        } else if main_gate_name == STD_MAIN_GATE_NAME_WITH_DNEXT {
            (
                0,
                1,
                2,
                3,
                SelectorOptimizedWidth4MainGateWithDNext::AB_MULTIPLICATION_TERM_COEFF_INDEX,
                SelectorOptimizedWidth4MainGateWithDNext::AC_MULTIPLICATION_TERM_COEFF_INDEX,
                SelectorOptimizedWidth4MainGateWithDNext::CONSTANT_TERM_COEFF_INDEX,
                SelectorOptimizedWidth4MainGateWithDNext::D_NEXT_TERM_COEFF_INDEX,
            )
        } else if main_gate_name == SELECTOR_OPTIMIZED_MAIN_GATE_NAME {
            (
                0,
                1,
                2,
                3,
                SelectorOptimizedWidth4MainGateWithDNext::AB_MULTIPLICATION_TERM_COEFF_INDEX,
                SelectorOptimizedWidth4MainGateWithDNext::AC_MULTIPLICATION_TERM_COEFF_INDEX,
                SelectorOptimizedWidth4MainGateWithDNext::CONSTANT_TERM_COEFF_INDEX,
                SelectorOptimizedWidth4MainGateWithDNext::D_NEXT_TERM_COEFF_INDEX,
            )
        } else {
            unreachable!("only 3 main gate types are allowed");
        };
        let num_first_round_evals = num_first_round_polys - 1;
        let num_copy_permutation_evals = 1;

        let num_evals = num_setup_polys + num_first_round_evals + num_copy_permutation_evals;
        let mut setup = SetupEvaluationOffsets::default();
        let shifted_first_round_pos = {
            setup.num_gate_setups = num_gate_setups;
            setup.num_permutations = num_permutations;
            setup.num_lookup_table_cols = num_lookup_table_cols;
            setup.num_gate_selectors = num_gate_selectors;
            setup.q_ab_term = q_ab_term;
            setup.q_ac_term = q_ac_term;
            setup.q_const_term = q_const_term;
            setup.q_dnext_term = q_dnext_term;

            setup.gate_setups_at_z = 0;
            setup.gate_selectors_at_z = num_gate_setups;
            setup.permutations_at_z = num_gate_setups + num_gate_selectors;
            setup.lookup_selector_at_z = num_gate_setups + num_gate_selectors + num_permutations;
            setup.lookup_tables_at_z = num_gate_setups + num_gate_selectors + num_permutations + some_lookup.map(|_| 1).unwrap_or(0);
            setup.lookup_table_type_at_z = num_gate_setups + num_gate_selectors + num_permutations + some_lookup.map(|_| 1).unwrap_or(0) + some_lookup.map(|_| num_lookup_table_cols).unwrap_or(0);

            let shifted_first_round_pos = if requires_setup_polys_opening_at_shifted_point(vk) {
                setup.gate_setups_at_z_omega = num_evals;
                setup.gate_selectors_at_z_omega = num_evals + num_gate_setups;
                setup.permutations_at_z_omega = num_evals + num_gate_setups + num_gate_selectors;
                setup.lookup_selector_at_z_omega = num_evals + num_gate_setups + num_gate_selectors + num_permutations;
                setup.lookup_tables_at_z_omega = num_evals + num_gate_setups + num_gate_selectors + num_permutations + some_lookup.map(|_| 1).unwrap_or(0);
                setup.lookup_table_type_at_z_omega =
                    num_evals + num_gate_setups + num_gate_selectors + num_permutations + some_lookup.map(|_| 1).unwrap_or(0) + some_lookup.map(|_| num_lookup_table_cols).unwrap_or(0);

                assert_eq!(num_evals + num_setup_polys, setup.lookup_table_type_at_z_omega + some_lookup.map(|_| 1).unwrap_or(0));

                setup.lookup_table_type_at_z_omega + some_lookup.map(|_| 1).unwrap_or(0)
            } else {
                assert_eq!(num_setup_polys, setup.lookup_table_type_at_z + some_lookup.map(|_| 1).unwrap_or(0));

                num_evals
            };

            shifted_first_round_pos
        };

        let mut trace = TraceAndGateEvaluationOffsets::default();
        trace.a_term = a_term;
        trace.b_term = b_term;
        trace.c_term = c_term;
        trace.d_term = d_term;

        trace.trace_evaluations_at_z = num_setup_polys;
        // shifted openings for first round
        let shifted_second_round_pos = if requires_trace_polys_opening_at_shifted_point(main_gate) {
            trace.trace_evaluations_at_z_omega = shifted_first_round_pos;
            trace.main_gate_quotient_at_z_omega = shifted_first_round_pos + num_state_polys;
            if has_custom_gate {
                trace.custom_gate_quotient_at_z_omega = shifted_first_round_pos + num_state_polys + 1;
                shifted_first_round_pos + num_state_polys + 1 + 1
            } else {
                shifted_first_round_pos + num_state_polys + 1
            }
        } else {
            shifted_first_round_pos
        };

        let mut copy_permutation = CopyPermutationEvaluationOffsets::default();
        copy_permutation.grand_product_at_z = num_setup_polys + num_first_round_evals;

        let mut lookup = LookupEvaluationOffsets::default();
        if has_lookup {
            lookup.s_poly_at_z = num_setup_polys + num_first_round_evals + num_copy_permutation_evals;
            lookup.grand_product_at_z = num_setup_polys + num_first_round_evals + num_copy_permutation_evals + 1;
        } else {
            assert_eq!(num_setup_polys + num_first_round_evals + num_copy_permutation_evals, copy_permutation.grand_product_at_z + 1);
        };

        // shifted openings for second round
        copy_permutation.grand_product_at_z_omega = shifted_second_round_pos;
        copy_permutation.first_quotient_at_z_omega = shifted_second_round_pos + 1;
        copy_permutation.second_quotient_at_z_omega = shifted_second_round_pos + 1 + 1;
        let shifted_lookup_pos = copy_permutation.second_quotient_at_z_omega + 1;
        let lookup = if has_lookup {
            lookup.s_poly_at_z_omega = shifted_lookup_pos;
            lookup.grand_product_at_z_omega = shifted_lookup_pos + 1;
            lookup.first_quotient_at_z_omega = shifted_lookup_pos + 1 + 1;
            lookup.second_quotient_at_z_omega = shifted_lookup_pos + 1 + 1 + 1;
            lookup.third_quotient_at_z_omega = shifted_lookup_pos + 1 + 1 + 1 + 1;

            Some(lookup)
        } else {
            None
        };

        Self {
            trace,
            copy_permutation,
            lookup,
            setup,
            c0: 0,
            c1: num_setup_polys,
            c2: num_setup_polys + num_first_round_evals,
            c0_shifted: num_setup_polys + num_first_round_evals + num_copy_permutation_evals,
            c1_shifted: shifted_first_round_pos,
            c2_shifted: shifted_second_round_pos,
            has_custom_gate,
        }
    }
}

pub fn flatten_all_evaluations<F: PrimeField>(
    setup_evaluations: &SetupEvaluations<F>,
    first_round_evaluations: &FirstRoundEvaluations<F>,
    second_round_evaluations: &SecondRoundEvaluations<F>,
) -> Vec<F> {
    let (c0_evals, c0_evals_shifted) = setup_evaluations.flatten();
    let (c1_evals, c1_evals_shifted) = first_round_evaluations.flatten();
    let (c2_evals, c2_evals_shifted) = second_round_evaluations.flatten();
    if setup_evaluations.requires_opening_at_shifted_point() {
        assert_eq!(c0_evals.len(), c0_evals_shifted.len());
    }
    if first_round_evaluations.requires_opening_at_shifted_point() {
        assert_eq!(c1_evals.len() + 1, c1_evals_shifted.len());
    }
    assert_eq!(c2_evals.len() + 2, c2_evals_shifted.len());

    // this is verification friendly representation of the evaluations
    // c0 || c1 || c2 ||  Option<c1 shifted> || c2 shifted
    let mut flattened_evaluations = vec![];
    flattened_evaluations.extend(c0_evals);
    flattened_evaluations.extend(c1_evals);
    flattened_evaluations.extend(c2_evals);
    flattened_evaluations.extend(c0_evals_shifted);
    flattened_evaluations.extend(c1_evals_shifted);
    flattened_evaluations.extend(c2_evals_shifted);

    flattened_evaluations
}
