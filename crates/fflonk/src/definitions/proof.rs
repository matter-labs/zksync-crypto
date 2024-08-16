use std::io::{Read, Write};

use super::*;

#[derive(Clone, Default, serde::Serialize, serde::Deserialize, PartialEq, Eq)]
#[serde(bound = "")]
pub struct FflonkProof<E: Engine, C: Circuit<E>> {
    pub n: usize,
    pub inputs: Vec<<E as ScalarEngine>::Fr>,
    pub c1: E::G1Affine,
    pub c2: E::G1Affine,
    pub w: E::G1Affine,
    pub w_prime: E::G1Affine,

    pub setup_evaluations: SetupEvaluations<E::Fr>,
    pub first_round_evaluations: FirstRoundEvaluations<E::Fr>,
    pub second_round_evaluations: SecondRoundEvaluations<E::Fr>,
    pub lagrange_basis_inverses: Vec<E::Fr>,

    _c: std::marker::PhantomData<C>,
}

#[derive(Clone,serde::Serialize, serde::Deserialize)]
pub struct FlattenedFflonkProof<E: Engine, C: Circuit<E>> {
    pub inputs: Vec<E::Fr>,
    pub commitments: Vec<E::G1Affine>,
    pub evaluations: Vec<E::Fr>,
    pub lagrange_basis_inverses: Vec<E::Fr>,
    _c: std::marker::PhantomData<C>,
}

impl<E: Engine, C: Circuit<E>> FlattenedFflonkProof<E, C> {
    pub fn into_original_proof(&self, vk: &FflonkVerificationKey<E, C>) -> FflonkProof<E, C> {
        let has_lookup = vk.total_lookup_entries_length > 0;
        let gates = sorted_gates_from_circuit_definitions::<_, C>();
        let has_custom_gate = gates.len() > 1;
        let num_state_polys = vk.num_state_polys;

        let evaluations = &self.evaluations;
        let offsets = EvaluationOffsets::from_vk(vk);

        let mut setup_evaluations = SetupEvaluations::default();
        {
            setup_evaluations.gate_setups_at_z = evaluations[offsets.setup.gate_setups_at_z
                ..offsets.setup.gate_setups_at_z + offsets.setup.num_gate_setups]
                .to_vec();

            setup_evaluations.gate_selectors_at_z = evaluations[offsets.setup.gate_selectors_at_z
                ..offsets.setup.gate_selectors_at_z + offsets.setup.num_gate_selectors]
                .to_vec();

            setup_evaluations.permutations_at_z = evaluations[offsets.setup.permutations_at_z
                ..offsets.setup.permutations_at_z + offsets.setup.num_permutations]
                .to_vec();
            if has_lookup {
                setup_evaluations.lookup_selector_at_z =
                    Some(evaluations[offsets.setup.lookup_selector_at_z]);

                setup_evaluations.lookup_tables_at_z = Some(
                    evaluations[offsets.setup.lookup_tables_at_z
                        ..offsets.setup.lookup_tables_at_z + offsets.setup.num_lookup_table_cols]
                        .to_vec(),
                );
                setup_evaluations.lookup_tables_at_z_omega = Some(
                    evaluations[offsets.setup.lookup_tables_at_z_omega
                        ..offsets.setup.lookup_tables_at_z_omega
                            + offsets.setup.num_lookup_table_cols]
                        .to_vec(),
                );
                setup_evaluations.lookup_table_type_at_z =
                    Some(evaluations[offsets.setup.lookup_table_type_at_z]);
            }
        }

        let mut trace_and_gate_evaluations = TraceAndGateEvaluations::default();
        let main_gate = gates[0].clone();
        {
            trace_and_gate_evaluations.trace_evaluations_at_z =
                evaluations[offsets.trace.trace_evaluations_at_z
                    ..offsets.trace.trace_evaluations_at_z + num_state_polys]
                    .to_vec();
            trace_and_gate_evaluations.main_gate_quotient_at_z =
                evaluations[offsets.trace.main_gate_quotient_at_z];

            if has_custom_gate {
                trace_and_gate_evaluations.custom_gate_quotient_at_z =
                    Some(evaluations[offsets.trace.custom_gate_quotient_at_z]);
            }
            if requires_trace_polys_opening_at_shifted_point(main_gate) {
                trace_and_gate_evaluations.trace_evaluations_at_z_omega = Some(
                    evaluations[offsets.trace.trace_evaluations_at_z_omega
                        ..offsets.trace.trace_evaluations_at_z_omega + num_state_polys]
                        .to_vec(),
                );
                trace_and_gate_evaluations.main_gate_quotient_at_z_omega =
                    Some(evaluations[offsets.trace.main_gate_quotient_at_z_omega]);
                if has_custom_gate {
                    trace_and_gate_evaluations.custom_gate_quotient_at_z_omega =
                        Some(evaluations[offsets.trace.custom_gate_quotient_at_z_omega]);
                }
            }
        }

        let first_round_evaluations = FirstRoundEvaluations {
            trace_and_gate_evaluations,
        };

        let mut copy_permutation_evaluations = CopyPermutationEvaluations::default();
        {
            copy_permutation_evaluations.grand_product_at_z =
                evaluations[offsets.copy_permutation.grand_product_at_z];
            copy_permutation_evaluations.grand_product_at_z_omega =
                evaluations[offsets.copy_permutation.grand_product_at_z_omega];
            copy_permutation_evaluations.first_quotient_at_z =
                evaluations[offsets.copy_permutation.first_quotient_at_z];
            copy_permutation_evaluations.first_quotient_at_z_omega =
                evaluations[offsets.copy_permutation.first_quotient_at_z_omega];
            copy_permutation_evaluations.second_quotient_at_z =
                evaluations[offsets.copy_permutation.second_quotient_at_z];
            copy_permutation_evaluations.second_quotient_at_z_omega =
                evaluations[offsets.copy_permutation.second_quotient_at_z_omega];
        }

        let lookup_evaluations = if has_lookup {
            let lookup_offsets = offsets.lookup.expect("lookup offsets");
            let mut lookup_evaluations = LookupEvaluations::default();
            lookup_evaluations.s_poly_at_z = evaluations[lookup_offsets.s_poly_at_z];
            lookup_evaluations.s_poly_at_z_omega = evaluations[lookup_offsets.s_poly_at_z_omega];
            lookup_evaluations.grand_product_at_z = evaluations[lookup_offsets.grand_product_at_z];
            lookup_evaluations.grand_product_at_z_omega =
                evaluations[lookup_offsets.grand_product_at_z_omega];

            lookup_evaluations.first_quotient_at_z =
                evaluations[lookup_offsets.first_quotient_at_z];
            lookup_evaluations.first_quotient_at_z_omega =
                evaluations[lookup_offsets.first_quotient_at_z_omega];
            lookup_evaluations.second_quotient_at_z =
                evaluations[lookup_offsets.second_quotient_at_z];
            lookup_evaluations.second_quotient_at_z_omega =
                evaluations[lookup_offsets.second_quotient_at_z_omega];
            lookup_evaluations.third_quotient_at_z =
                evaluations[lookup_offsets.third_quotient_at_z];
            lookup_evaluations.third_quotient_at_z_omega =
                evaluations[lookup_offsets.third_quotient_at_z_omega];
            Some(lookup_evaluations)
        } else {
            None
        };

        let second_round_evaluations = SecondRoundEvaluations {
            copy_permutation_evaluations,
            lookup_evaluations,
        };

        let public_inputs = self.inputs.clone();
        let c1 = self.commitments[0];
        let c2 = self.commitments[1];
        let w = self.commitments[2];
        let w_prime = self.commitments[3];
        FflonkProof {
            n: vk.n,
            inputs: public_inputs,
            c1,
            c2,
            w,
            w_prime,
            setup_evaluations,
            first_round_evaluations,
            second_round_evaluations,
            lagrange_basis_inverses: self.lagrange_basis_inverses.clone(),
            _c: std::marker::PhantomData,
        }
    }

    pub fn serialize_into_evm_format<W: std::io::Write>(
        &self,
        writer: W,
    ) -> Result<(), serde_json::error::Error> {
        use franklin_crypto::bellman::to_hex;

        #[derive(serde::Serialize)]
        struct G1Point {
            x: String,
            y: String,
        }
        #[derive(serde::Serialize)]
        struct JsonProof {
            pub inputs: Vec<String>,
            pub commitments: Vec<G1Point>,
            pub evaluations: Vec<String>,
            pub lagrange_basis_inverses: Vec<String>,
        }

        let mut json_proof = JsonProof {
            inputs: vec![],
            commitments: vec![],
            evaluations: vec![],
            lagrange_basis_inverses: vec![],
        };

        for el in self.inputs.iter() {
            let hex = to_hex(el);
            json_proof.inputs.push(hex);
        }

        for el in self.commitments.iter() {
            let (x, y) = el.into_xy_unchecked();
            let x = to_hex(&x);
            let y = to_hex(&y);
            json_proof.commitments.push(G1Point { x, y });
        }

        for el in self.evaluations.iter() {
            let hex = to_hex(el);
            json_proof.evaluations.push(hex);
        }

        for el in self.lagrange_basis_inverses.iter() {
            let hex = to_hex(el);
            json_proof.lagrange_basis_inverses.push(hex);
        }

        serde_json::to_writer(writer, &json_proof)
    }
}

impl<E: Engine, C: Circuit<E>> FflonkProof<E, C> {
    pub fn empty() -> Self {
        Self {
            n: 0,
            inputs: vec![],
            c1: E::G1Affine::zero(),
            c2: E::G1Affine::zero(),
            w: E::G1Affine::zero(),
            w_prime: E::G1Affine::zero(),
            setup_evaluations: Default::default(),
            first_round_evaluations: Default::default(),
            second_round_evaluations: Default::default(),
            lagrange_basis_inverses: Default::default(),
            _c: std::marker::PhantomData,
        }
    }

    pub fn flatten(self) -> FlattenedFflonkProof<E, C> {
        let flattened_evaluations = flatten_all_evaluations(
            &self.setup_evaluations,
            &self.first_round_evaluations,
            &self.second_round_evaluations,
        );
        let mut commitments = vec![self.c1, self.c2];
        commitments.push(self.w);
        commitments.push(self.w_prime);

        FlattenedFflonkProof {
            inputs: self.inputs,
            commitments,
            evaluations: flattened_evaluations,
            lagrange_basis_inverses: self.lagrange_basis_inverses,
            _c: std::marker::PhantomData,
        }
    }

    pub fn read<R: Read>(mut reader: R) -> std::io::Result<Self> {
        // serde serialize doesn't work for this type so this is a dirty way of serialization
        let n = serde_json::from_reader(&mut reader)?;
        let inputs = serde_json::from_reader(&mut reader)?;
        let c1 = serde_json::from_reader(&mut reader)?;
        let c2 = serde_json::from_reader(&mut reader)?;
        let w = serde_json::from_reader(&mut reader)?;
        let w_prime = serde_json::from_reader(&mut reader)?;
        let setup_evaluations = serde_json::from_reader(&mut reader)?;
        let first_round_evaluations = serde_json::from_reader(&mut reader)?;
        let second_round_evaluations = serde_json::from_reader(&mut reader)?;
        let lagrange_basis_inverses = serde_json::from_reader(&mut reader)?;

        Ok(Self {
            n,
            inputs,
            c1,
            c2,
            w,
            w_prime,
            setup_evaluations,
            first_round_evaluations,
            second_round_evaluations,
            lagrange_basis_inverses,
            _c: std::marker::PhantomData,
        })
    }

    pub fn write<W: Write>(&self, mut writer: W) -> std::io::Result<()> {
        serde_json::to_writer(&mut writer, &self.n)?;
        serde_json::to_writer(&mut writer, &self.inputs)?;
        serde_json::to_writer(&mut writer, &self.c1)?;
        serde_json::to_writer(&mut writer, &self.c2)?;
        serde_json::to_writer(&mut writer, &self.w)?;
        serde_json::to_writer(&mut writer, &self.w_prime)?;
        serde_json::to_writer(&mut writer, &self.setup_evaluations)?;
        serde_json::to_writer(&mut writer, &self.first_round_evaluations)?;
        serde_json::to_writer(&mut writer, &self.second_round_evaluations)?;

        Ok(())
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
    pub fn from_vk<E: Engine, C: Circuit<E>>(vk: &FflonkVerificationKey<E, C>) -> Self {
        let (num_setup_polys, num_first_round_polys, num_second_round_polys, max_num_polys) =
            num_system_polys_from_vk(vk);
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
        let (a_term, b_term, c_term, d_term, q_ab_term, q_ac_term, q_const_term, q_dnext_term) =
            if main_gate_name == STD_MAIN_GATE_NAME {
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
            setup.lookup_tables_at_z = num_gate_setups
                + num_gate_selectors
                + num_permutations
                + some_lookup.map(|_| 1).unwrap_or(0);
            setup.lookup_table_type_at_z = num_gate_setups
                + num_gate_selectors
                + num_permutations
                + some_lookup.map(|_| 1).unwrap_or(0)
                + some_lookup.map(|_| num_lookup_table_cols).unwrap_or(0);

            let shifted_first_round_pos = if requires_setup_polys_opening_at_shifted_point(vk) {
                setup.gate_setups_at_z_omega = num_all_polys;
                setup.gate_selectors_at_z_omega = num_all_polys + num_gate_setups;
                setup.permutations_at_z_omega =
                    num_all_polys + num_gate_setups + num_gate_selectors;
                setup.lookup_selector_at_z_omega =
                    num_all_polys + num_gate_setups + num_gate_selectors + num_permutations;
                setup.lookup_tables_at_z_omega = num_all_polys
                    + num_gate_setups
                    + num_gate_selectors
                    + num_permutations
                    + some_lookup.map(|_| 1).unwrap_or(0);
                setup.lookup_table_type_at_z_omega = num_all_polys
                    + num_gate_setups
                    + num_gate_selectors
                    + num_permutations
                    + some_lookup.map(|_| 1).unwrap_or(0)
                    + some_lookup.map(|_| num_lookup_table_cols).unwrap_or(0);

                assert_eq!(
                    num_all_polys + num_setup_polys,
                    setup.lookup_table_type_at_z_omega + some_lookup.map(|_| 1).unwrap_or(0)
                );

                setup.lookup_table_type_at_z_omega + some_lookup.map(|_| 1).unwrap_or(0)
            } else {
                assert_eq!(
                    num_setup_polys,
                    setup.lookup_table_type_at_z + some_lookup.map(|_| 1).unwrap_or(0)
                );

                num_all_polys
            };

            shifted_first_round_pos
        };

        let mut trace = TraceAndGateEvaluationOffsets::default();
        trace.a_term = a_term;
        trace.b_term = b_term;
        trace.c_term = c_term;
        trace.d_term = d_term;

        trace.trace_evaluations_at_z = num_setup_polys;
        trace.main_gate_quotient_at_z = num_setup_polys + num_state_polys;
        if has_custom_gate {
            trace.custom_gate_quotient_at_z = num_setup_polys + num_state_polys + 1;
        }

        // shifted openings for first round
        let shifted_second_round_pos = if requires_trace_polys_opening_at_shifted_point(main_gate) {
            trace.trace_evaluations_at_z_omega = shifted_first_round_pos;
            trace.main_gate_quotient_at_z_omega = shifted_first_round_pos + num_state_polys;
            if has_custom_gate {
                trace.custom_gate_quotient_at_z_omega =
                    shifted_first_round_pos + num_state_polys + 1;
                shifted_first_round_pos + num_state_polys + 1 + 1
            } else {
                shifted_first_round_pos + num_state_polys + 1
            }
        } else {
            shifted_first_round_pos
        };

        let mut copy_permutation = CopyPermutationEvaluationOffsets::default();
        copy_permutation.grand_product_at_z = num_setup_polys + num_first_round_polys;
        copy_permutation.first_quotient_at_z = num_setup_polys + num_first_round_polys + 1;
        copy_permutation.second_quotient_at_z = num_setup_polys + num_first_round_polys + 1 + 1;

        let num_copy_permutation_polys = 1 + 1 + 1;

        let mut lookup = LookupEvaluationOffsets::default();
        if has_lookup {
            lookup.s_poly_at_z =
                num_setup_polys + num_first_round_polys + num_copy_permutation_polys;
            lookup.grand_product_at_z =
                num_setup_polys + num_first_round_polys + num_copy_permutation_polys + 1;
            lookup.first_quotient_at_z =
                num_setup_polys + num_first_round_polys + num_copy_permutation_polys + 1 + 1;
            lookup.second_quotient_at_z =
                num_setup_polys + num_first_round_polys + num_copy_permutation_polys + 1 + 1 + 1;
            lookup.third_quotient_at_z = num_setup_polys
                + num_first_round_polys
                + num_copy_permutation_polys
                + 1
                + 1
                + 1
                + 1;
            assert_eq!(
                num_setup_polys + num_first_round_polys + num_second_round_polys,
                lookup.third_quotient_at_z + 1
            );
        } else {
            assert_eq!(
                num_setup_polys + num_first_round_polys + num_second_round_polys,
                copy_permutation.second_quotient_at_z + 1
            );
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
            c2: num_setup_polys + num_first_round_polys,
            c0_shifted: num_all_polys,
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
        assert_eq!(c1_evals.len(), c1_evals_shifted.len());
    }
    assert_eq!(c2_evals.len(), c2_evals_shifted.len());

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
