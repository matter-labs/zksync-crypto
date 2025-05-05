use crate::{
    config::DevCSConfig,
    cs::{
        cs_builder::{new_builder, CsBuilder, CsBuilderImpl},
        cs_builder_reference::CsReferenceImplementationBuilder,
        gates::{
            BooleanConstraintGate, ConstantsAllocatorGate, DotProductGate,
            FmaGateInBaseFieldWithoutConstant, NopGate, ReductionGate, SelectionGate, U8x4FMAGate,
            UIntXAddGate, ZeroCheckGate,
        },
        implementations::reference_cs::CSReferenceImplementation,
        traits::{cs::ConstraintSystem, gate::GatePlacementStrategy},
        CSGeometry, GateConfigurationHolder, LookupParameters, StaticToolboxHolder,
    },
    field::{goldilocks::GoldilocksField, SmallField},
    gadgets::non_native_field::implementations::NonNativeFieldOverU16Params,
};

use crate::gadgets::tower_extension::tests::json::types::{
    BN256BaseNNFieldParams, BN256ScalarNNFieldParams,
};

type F = GoldilocksField;
type P = GoldilocksField;

/// Creates a test constraint system for testing purposes that includes the
/// majority (even possibly unneeded) of the gates and tables.
pub fn create_test_cs(
    max_trace_len: usize,
) -> CSReferenceImplementation<
    F,
    P,
    DevCSConfig,
    impl GateConfigurationHolder<F>,
    impl StaticToolboxHolder,
> {
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
                width: 1,
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
        // let owned_cs = ReductionGate::<F, 4>::configure_for_cs(owned_cs, GatePlacementStrategy::UseSpecializedColumns { num_repetitions: 8, share_constants: true });
        let builder = BooleanConstraintGate::configure_builder(
            builder,
            GatePlacementStrategy::UseGeneralPurposeColumns,
        );
        let builder = UIntXAddGate::<32>::configure_builder(
            builder,
            GatePlacementStrategy::UseGeneralPurposeColumns,
        );
        let builder = UIntXAddGate::<16>::configure_builder(
            builder,
            GatePlacementStrategy::UseGeneralPurposeColumns,
        );
        let builder = UIntXAddGate::<8>::configure_builder(
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
        let builder = DotProductGate::<4>::configure_builder(
            builder,
            GatePlacementStrategy::UseGeneralPurposeColumns,
        );
        let builder =
            NopGate::configure_builder(builder, GatePlacementStrategy::UseGeneralPurposeColumns);

        builder
    }

    let builder_impl =
        CsReferenceImplementationBuilder::<F, P, DevCSConfig>::new(geometry, max_trace_len);
    let builder = new_builder::<_, F>(builder_impl);

    let builder = configure(builder);
    let mut owned_cs = builder.build(max_variables);

    use crate::gadgets::tables::create_range_check_16_bits_table;
    use crate::gadgets::tables::RangeCheck16BitsTable;

    let table = create_range_check_16_bits_table();
    owned_cs.add_lookup_table::<RangeCheck16BitsTable<1>, 1>(table);

    owned_cs
}

/// Returns BN254 base field parameters
pub fn bn254_base_field_params() -> BN256BaseNNFieldParams {
    NonNativeFieldOverU16Params::create()
}

/// Returns BN254 scalar field parameters
pub fn bn254_scalar_field_params() -> BN256ScalarNNFieldParams {
    NonNativeFieldOverU16Params::create()
}
