
#[cfg(test)]
mod tests {
    use crate::{config::DevCSConfig, cs::{cs_builder::{new_builder, CsBuilder, CsBuilderImpl}, cs_builder_reference::CsReferenceImplementationBuilder, gates::{BooleanConstraintGate, ConstantsAllocatorGate, DotProductGate, FmaGateInBaseFieldWithoutConstant, NopGate, ReductionGate, SelectionGate, U8x4FMAGate, UIntXAddGate, ZeroCheckGate}, traits::{cs::ConstraintSystem, gate::GatePlacementStrategy}, CSGeometry, GateConfigurationHolder, LookupParameters, StaticToolboxHolder}, field::SmallField, gadgets::{tables::{create_and8_table, create_byte_split_table, create_xor8_table, And8Table, Xor8Table}, tower_extension::tests::utils::cs::create_test_cs, traits::{allocatable::CSAllocatable, witnessable::WitnessHookable}, u256::UInt256}, worker::Worker};
    use ethereum_types::U256;
    use crate::field::goldilocks::GoldilocksField;

    type F = GoldilocksField;
    type P = GoldilocksField;
    #[test]
    fn test_modmul() {
        let mut owned_cs = create_test_cs(1 << 21);
        let cs = &mut owned_cs;

        let a_value = U256::from(123456789u64);
        let b_value = U256::from(987654321u64);
        let modulo_value = U256::from(11111u64);


        let a = UInt256::allocate(cs, a_value);
        let b = UInt256::allocate(cs, b_value);
        let modulo = UInt256::allocate(cs, modulo_value);


        let result = a.modmul(cs, &b, &modulo);


        let expected_value = (a_value.full_mul(b_value) % modulo_value).try_into().unwrap();


        let result_value = result.witness_hook(cs)().unwrap();
        assert_eq!(result_value, expected_value, "modmul result is incorrect");
    }
    #[test]
    fn test_modmul_2() {
        let geometry = CSGeometry {
            num_columns_under_copy_permutation: 200,
            num_witness_columns: 0,
            num_constant_columns: 8,
            max_allowed_constraint_degree: 4,
        };
        let max_variables = 1 << 27;
    
        fn configure<
            F: SmallField,
            T: CsBuilderImpl<F, T>,
            GC: GateConfigurationHolder<F>,
            TB: StaticToolboxHolder,
        >(
            builder: CsBuilder<T, F, GC, TB>,
        ) -> CsBuilder<T, F, impl GateConfigurationHolder<F>, impl StaticToolboxHolder> {
            let builder = builder.allow_lookup(
                LookupParameters::UseSpecializedColumnsWithTableIdAsConstant {
                    width: 3,
                    num_repetitions: 20,
                    share_table_id: true,
                },
            );
            let builder = U8x4FMAGate::configure_builder(
                builder,
                GatePlacementStrategy::UseGeneralPurposeColumns,
            );
            let builder = ConstantsAllocatorGate::configure_builder(
                builder,
                GatePlacementStrategy::UseGeneralPurposeColumns,
            );
            let builder = FmaGateInBaseFieldWithoutConstant::configure_builder(
                builder,
                GatePlacementStrategy::UseGeneralPurposeColumns,
            );
            let builder = ReductionGate::<F, 4>::configure_builder(
                builder,
                GatePlacementStrategy::UseGeneralPurposeColumns,
            );
            let builder = UIntXAddGate::<32>::configure_builder(
                builder,
                GatePlacementStrategy::UseGeneralPurposeColumns,
            );
            let builder = SelectionGate::configure_builder(
                builder,
                GatePlacementStrategy::UseGeneralPurposeColumns,
            );
            let builder = ZeroCheckGate::configure_builder(
                builder,
                GatePlacementStrategy::UseGeneralPurposeColumns,
                false,
            );

            let builder =
                NopGate::configure_builder(builder, GatePlacementStrategy::UseGeneralPurposeColumns);
    
            builder
        }
    
        let builder_impl =
            CsReferenceImplementationBuilder::<F, P, DevCSConfig>::new(geometry, 1 << 21);
        let builder = new_builder::<_, F>(builder_impl);
    
        let builder = configure(builder);
        let mut owned_cs = builder.build(max_variables);
    
        let table = create_xor8_table();
        owned_cs.add_lookup_table::<Xor8Table, 3>(table);


        let cs = &mut owned_cs;

        let a_value = U256::max_value() - U256::from(1u64);
        let b_value = U256::max_value() - U256::from(2u64);
        let modulo_value = U256::from(1u64);

        let a = UInt256::allocate(cs, a_value);
        let b = UInt256::allocate(cs, b_value);
        let modulo = UInt256::allocate(cs, modulo_value);

        let result = a.modmul(cs, &b, &modulo);

        let expected_value = (a_value.full_mul(b_value) % modulo_value).try_into().unwrap();
        let result_value = result.witness_hook(cs)().unwrap();
        assert_eq!(result_value, expected_value, "modmul result is incorrect");

        cs.pad_and_shrink();
        let worker = Worker::new();
        let mut owned_cs = owned_cs.into_assembly::<std::alloc::Global>();
        owned_cs.print_gate_stats();
        assert!(owned_cs.check_if_satisfied(&worker));
    }
}