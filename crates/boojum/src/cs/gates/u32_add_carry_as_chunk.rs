use super::*;
use crate::cs::cs_builder::*;
use crate::gadgets::traits::castable::WitnessCastable;

// for a, b, c being u8x4 we output a + b mod 2^32 as u8x4 and a + b / 2^32
// `carry_out` is boolean constrainted
// but `c` parts are NOT. We will use separate range checks

#[derive(Derivative)]
#[derivative(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct U32AddCarryAsChunkConstraintEvaluator;

impl U32AddCarryAsChunkConstraintEvaluator {
    const fn principal_width() -> usize {
        4 * 2 + // input
        4 + 1 // output
    }
}

impl<F: PrimeField> GateConstraintEvaluator<F> for U32AddCarryAsChunkConstraintEvaluator {
    type UniqueParameterizationParams = ();

    #[inline(always)]
    fn new_from_parameters(_params: Self::UniqueParameterizationParams) -> Self {
        Self
    }

    #[inline(always)]
    fn unique_params(&self) -> Self::UniqueParameterizationParams {}

    #[inline]
    fn type_name() -> std::borrow::Cow<'static, str> {
        Cow::Borrowed(std::any::type_name::<Self>())
    }

    #[inline]
    fn instance_width(&self) -> GatePrincipalInstanceWidth {
        GatePrincipalInstanceWidth {
            num_variables: Self::principal_width(),
            num_witnesses: 0,
            num_constants: 0,
        }
    }

    #[inline]
    fn gate_purpose() -> GatePurpose {
        GatePurpose::Evaluatable {
            max_constraint_degree: 2,
            num_quotient_terms: 2,
        }
    }

    #[inline]
    fn placement_type(&self) -> GatePlacementType {
        GatePlacementType::MultipleOnRow {
            per_chunk_offset: PerChunkOffset {
                variables_offset: Self::principal_width(),
                witnesses_offset: 0,
                constants_offset: 0,
            },
        }
    }

    #[inline]
    fn num_repetitions_in_geometry(&self, geometry: &CSGeometry) -> usize {
        debug_assert!(geometry.num_columns_under_copy_permutation >= Self::principal_width());

        geometry.num_columns_under_copy_permutation / Self::principal_width()
    }

    #[inline]
    fn num_required_constants_in_geometry(&self, _geometry: &CSGeometry) -> usize {
        0
    }

    type GlobalConstants<P: field::traits::field_like::PrimeFieldLike<Base = F>> = [P; 4];

    #[inline(always)]
    fn create_global_constants<P: field::traits::field_like::PrimeFieldLike<Base = F>>(
        &self,
        ctx: &mut P::Context,
    ) -> Self::GlobalConstants<P> {
        [
            P::constant(F::from_u64_with_reduction(1u64 << 8), ctx),
            P::constant(F::from_u64_with_reduction(1u64 << 16), ctx),
            P::constant(F::from_u64_with_reduction(1u64 << 24), ctx),
            P::constant(F::from_u64_with_reduction(1u64 << 32), ctx),
        ]
    }

    type RowSharedConstants<P: field::traits::field_like::PrimeFieldLike<Base = F>> = [P; 0];

    #[inline(always)]
    fn load_row_shared_constants<
        P: field::traits::field_like::PrimeFieldLike<Base = F>,
        S: TraceSource<F, P>,
    >(
        &self,
        _trace_source: &S,
        _ctx: &mut P::Context,
    ) -> Self::RowSharedConstants<P> {
        []
    }

    #[inline(always)]
    fn evaluate_once<
        P: field::traits::field_like::PrimeFieldLike<Base = F>,
        S: TraceSource<F, P>,
        D: EvaluationDestination<F, P>,
    >(
        &self,
        trace_source: &S,
        destination: &mut D,
        _shared_constants: &Self::RowSharedConstants<P>,
        global_constants: &Self::GlobalConstants<P>,
        ctx: &mut P::Context,
    ) {
        let [shift8, shift16, shift24, shift32] = global_constants;

        let one = P::one(ctx);

        let a0 = trace_source.get_variable_value(0);
        let a1 = trace_source.get_variable_value(1);
        let a2 = trace_source.get_variable_value(2);
        let a3 = trace_source.get_variable_value(3);

        let b0 = trace_source.get_variable_value(4 + 0);
        let b1 = trace_source.get_variable_value(4 + 1);
        let b2 = trace_source.get_variable_value(4 + 2);
        let b3 = trace_source.get_variable_value(4 + 3);

        let out0 = trace_source.get_variable_value(8 + 0);
        let out1 = trace_source.get_variable_value(8 + 1);
        let out2 = trace_source.get_variable_value(8 + 2);
        let out3 = trace_source.get_variable_value(8 + 3);

        let carry = trace_source.get_variable_value(12);

        let mut contribution = P::zero(ctx);

        let coeffs = [&one, shift8, shift16, shift24];

        for (a, coeff) in [a0, a1, a2, a3].into_iter().zip(coeffs.iter()) {
            P::mul_and_accumulate_into(&mut contribution, &a, *coeff, ctx);
        }

        for (b, coeff) in [b0, b1, b2, b3].into_iter().zip(coeffs.iter()) {
            P::mul_and_accumulate_into(&mut contribution, &b, *coeff, ctx);
        }

        contribution.sub_assign(&out0, ctx);

        let mut tmp = out1;
        tmp.mul_assign(shift8, ctx);
        contribution.sub_assign(&tmp, ctx);

        let mut tmp = out2;
        tmp.mul_assign(shift16, ctx);
        contribution.sub_assign(&tmp, ctx);

        let mut tmp = out3;
        tmp.mul_assign(shift24, ctx);
        contribution.sub_assign(&tmp, ctx);

        let mut tmp = carry;
        tmp.mul_assign(shift32, ctx);
        contribution.sub_assign(&tmp, ctx);

        destination.push_evaluation_result(contribution, ctx);

        // carry is boolean

        let mut contribution = carry;
        contribution.mul_assign(&carry, ctx);
        contribution.sub_assign(&carry, ctx);

        destination.push_evaluation_result(contribution, ctx);
    }
}

#[derive(Derivative)]
#[derivative(Clone, Debug, PartialEq, Eq, Hash)]
pub struct U32AddCarryAsChunkGate {
    pub a: [Variable; 4],
    pub b: [Variable; 4],
    pub out: [Variable; 4],
    pub carry_out: Variable,
}

// individual for each width

impl<F: SmallField> Gate<F> for U32AddCarryAsChunkGate {
    #[inline(always)]
    fn check_compatible_with_cs<CS: ConstraintSystem<F>>(&self, cs: &CS) -> bool {
        let geometry = cs.get_params();
        geometry.max_allowed_constraint_degree >= 2
            && geometry.num_columns_under_copy_permutation
                >= U32AddCarryAsChunkConstraintEvaluator::principal_width()
    }

    // even though we make width a type parameter here, they share the evaluator
    type Evaluator = U32AddCarryAsChunkConstraintEvaluator;

    #[inline]
    fn evaluator(&self) -> Self::Evaluator {
        U32AddCarryAsChunkConstraintEvaluator
    }
}

impl U32AddCarryAsChunkGate {
    pub fn configure_builder<
        F: SmallField,
        GC: GateConfigurationHolder<F>,
        TB: StaticToolboxHolder,
        TImpl: CsBuilderImpl<F, TImpl>,
    >(
        builder: CsBuilder<TImpl, F, GC, TB>,
        placement_strategy: GatePlacementStrategy,
        // ) -> CsBuilder<TImpl, F, GC::DescendantHolder<Self, NextGateCounterWithoutParams>, TB> {
    ) -> CsBuilder<TImpl, F, (GateTypeEntry<F, Self, NextGateCounterWithoutParams>, GC), TB> {
        builder.allow_gate(placement_strategy, (), None)
    }

    pub fn add_to_cs<F: SmallField, CS: ConstraintSystem<F>>(self, cs: &mut CS) {
        debug_assert!(F::CAPACITY_BITS >= 33);
        debug_assert!(cs.gate_is_allowed::<Self>());

        if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == false {
            return;
        }

        match cs.get_gate_placement_strategy::<Self>() {
            GatePlacementStrategy::UseGeneralPurposeColumns => {
                let offered_row_idx = cs.next_available_row();
                let capacity_per_row = self.capacity_per_row(&*cs);
                let tooling: &mut NextGateCounterWithoutParams = cs
                    .get_gates_config_mut()
                    .get_aux_data_mut::<Self, _>()
                    .expect("gate must be allowed");
                let (row, num_instances_already_placed) =
                    find_next_gate_without_params(tooling, capacity_per_row, offered_row_idx);
                drop(tooling);

                // now we can use methods of CS to inform it of low level operations
                let mut offset = num_instances_already_placed
                    * U32AddCarryAsChunkConstraintEvaluator::principal_width();
                if offered_row_idx == row {
                    cs.place_gate(&self, row);
                }
                cs.place_multiple_variables_into_row(&self.a, row, offset);
                offset += 4;
                cs.place_multiple_variables_into_row(&self.b, row, offset);
                offset += 4;
                cs.place_multiple_variables_into_row(&self.out, row, offset);
                offset += 4;
                cs.place_variable(self.carry_out, row, offset);
            }
            GatePlacementStrategy::UseSpecializedColumns {
                num_repetitions,
                share_constants: _,
            } => {
                // gate knows how to place itself
                let capacity_per_row = num_repetitions;
                let tooling: &mut NextGateCounterWithoutParams = cs
                    .get_gates_config_mut()
                    .get_aux_data_mut::<Self, _>()
                    .expect("gate must be allowed");
                let (row, num_instances_already_placed) =
                    find_next_specialized_gate_without_params(tooling, capacity_per_row);
                cs.place_gate_specialized(&self, num_instances_already_placed, row);
                let mut starting_column = 0;
                cs.place_multiple_variables_into_row_specialized::<Self, 4>(
                    &self.a,
                    num_instances_already_placed,
                    row,
                    starting_column,
                );
                starting_column += 4;
                cs.place_multiple_variables_into_row_specialized::<Self, 4>(
                    &self.b,
                    num_instances_already_placed,
                    row,
                    starting_column,
                );
                starting_column += 4;
                cs.place_multiple_variables_into_row_specialized::<Self, 4>(
                    &self.out,
                    num_instances_already_placed,
                    row,
                    starting_column,
                );
                starting_column += 4;
                cs.place_variable_specialized::<Self>(
                    self.carry_out,
                    num_instances_already_placed,
                    row,
                    starting_column,
                );
            }
        }
    }

    pub fn perform_addition<F: SmallField, CS: ConstraintSystem<F>>(
        cs: &mut CS,
        a: [Variable; 4],
        b: [Variable; 4],
    ) -> ([Variable; 4], Variable) {
        debug_assert!(F::CAPACITY_BITS >= 33);
        debug_assert!(cs.gate_is_allowed::<Self>());

        let output_variables = cs.alloc_multiple_variables_without_values::<5>();

        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 8]| {
                let inputs = inputs.map(|el| <u8 as WitnessCastable<F, F>>::cast_from_source(el));
                let [a0, a1, a2, a3, b0, b1, b2, b3] = inputs;

                let a = u32::from_le_bytes([a0, a1, a2, a3]);
                let b = u32::from_le_bytes([b0, b1, b2, b3]);

                let tmp = (a as u64) + (b as u64);
                let out = tmp as u32;
                let carry = (tmp >> 32) as u32;
                let [out0, out1, out2, out3] = out.to_le_bytes();

                [
                    F::from_u64_with_reduction(out0 as u64),
                    F::from_u64_with_reduction(out1 as u64),
                    F::from_u64_with_reduction(out2 as u64),
                    F::from_u64_with_reduction(out3 as u64),
                    F::from_u64_with_reduction(carry as u64),
                ]
            };

            let mut dependencies = [Place::placeholder(); 8];
            dependencies[0..4].copy_from_slice(&Place::from_variables(a));
            dependencies[4..8].copy_from_slice(&Place::from_variables(b));

            cs.set_values_with_dependencies(
                &dependencies,
                &Place::from_variables(output_variables),
                value_fn,
            );
        }

        let [out0, out1, out2, out3, carry] = output_variables;

        if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP {
            let gate = Self {
                a,
                b,
                out: [out0, out1, out2, out3],
                carry_out: carry,
            };

            gate.add_to_cs(cs);
        }

        ([out0, out1, out2, out3], carry)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::cs::gates::testing_tools::test_evaluator;
    use crate::field::goldilocks::GoldilocksField;
    type F = GoldilocksField;

    #[test]
    fn test_properties() {
        // particular geometry is not important
        let _geometry = CSGeometry {
            num_columns_under_copy_permutation: 80,
            num_witness_columns: 80,
            num_constant_columns: 10,
            max_allowed_constraint_degree: 8,
        };

        let evaluator = <U32AddCarryAsChunkConstraintEvaluator as GateConstraintEvaluator<F>>::new_from_parameters(
            ()
        );

        test_evaluator::<F, _>(evaluator);
    }
}
