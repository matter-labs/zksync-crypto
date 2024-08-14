use std::any::TypeId;
use std::collections::HashMap;

use crate::franklin_crypto::bellman::pairing::Engine;
use crate::franklin_crypto::bellman::plonk::better_better_cs::cs::ConstraintSystem;

use crate::boojum::cs::cs_builder::new_builder;
use crate::boojum::cs::cs_builder::CsBuilderImpl;
use crate::boojum::cs::traits::circuit::{CircuitBuilder, CircuitBuilderProxy};
use crate::boojum::cs::traits::evaluator::GateBatchEvaluationComparisonFunction;
use crate::boojum::cs::traits::evaluator::GatePlacementType;
use crate::boojum::cs::traits::evaluator::PerChunkOffset;
use crate::boojum::cs::{CSGeometry, LookupParameters};
use crate::boojum::field::goldilocks::GoldilocksField as GL;

use crate::traits::circuit::ErasedBuilderForWrapperVerifier;
use crate::verifier_structs::gate_evaluator::TypeErasedGateEvaluationWrapperVerificationFunction;
use crate::verifier_structs::WrapperVerifier;

impl<E: Engine, CS: ConstraintSystem<E> + 'static, T: CircuitBuilder<GL>> ErasedBuilderForWrapperVerifier<E, CS> for CircuitBuilderProxy<GL, T> {
    fn geometry(&self) -> CSGeometry {
        T::geometry()
    }

    fn lookup_parameters(&self) -> LookupParameters {
        T::lookup_parameters()
    }

    fn create_wrapper_verifier(&self, cs: &mut CS) -> WrapperVerifier<E, CS> {
        let geometry = T::geometry();
        let builder_impl = CsWrapperVerifierBuilder::<'_, E, CS>::new_from_parameters(cs, geometry);
        let builder = new_builder::<_, GL>(builder_impl);

        let builder = T::configure_builder(builder);
        let verifier = builder.build(());

        verifier
    }
}

pub struct CsWrapperVerifierBuilder<'a, E: Engine, CS: ConstraintSystem<E> + 'static> {
    pub(crate) cs: &'a mut CS,

    pub parameters: CSGeometry,
    pub lookup_parameters: LookupParameters,

    pub(crate) gate_type_ids_for_specialized_columns: Vec<TypeId>,
    pub(crate) evaluators_over_specialized_columns: Vec<TypeErasedGateEvaluationWrapperVerificationFunction<E, CS>>,
    pub(crate) offsets_for_specialized_evaluators: Vec<(PerChunkOffset, PerChunkOffset, usize)>,

    pub(crate) evaluators_over_general_purpose_columns: Vec<TypeErasedGateEvaluationWrapperVerificationFunction<E, CS>>,
    pub(crate) general_purpose_evaluators_comparison_functions: HashMap<TypeId, Vec<(GateBatchEvaluationComparisonFunction, usize)>>,

    pub(crate) total_num_variables_for_specialized_columns: usize,
    pub(crate) total_num_witnesses_for_specialized_columns: usize,
    pub(crate) total_num_constants_for_specialized_columns: usize,
}

impl<'a, E: Engine, CS: ConstraintSystem<E> + 'static> CsWrapperVerifierBuilder<'a, E, CS> {
    pub fn new_from_parameters(cs: &'a mut CS, parameters: CSGeometry) -> Self {
        Self {
            cs,

            parameters: parameters,
            lookup_parameters: LookupParameters::NoLookup,

            gate_type_ids_for_specialized_columns: Vec::with_capacity(16),
            evaluators_over_specialized_columns: Vec::with_capacity(16),
            offsets_for_specialized_evaluators: Vec::with_capacity(16),

            evaluators_over_general_purpose_columns: Vec::with_capacity(16),
            general_purpose_evaluators_comparison_functions: HashMap::with_capacity(16),

            total_num_variables_for_specialized_columns: 0,
            total_num_witnesses_for_specialized_columns: 0,
            total_num_constants_for_specialized_columns: 0,
        }
    }
}

use crate::boojum::cs::cs_builder::CsBuilder;
use crate::boojum::cs::gates::lookup_marker::LookupFormalGate;
use crate::boojum::cs::gates::LookupTooling;
use crate::boojum::cs::traits::gate::GatePlacementStrategy;
use crate::boojum::cs::traits::{evaluator::GateConstraintEvaluator, gate::Gate};
use crate::boojum::cs::GateConfigurationHolder;
use crate::boojum::cs::GateTypeEntry;
use crate::boojum::cs::StaticToolboxHolder;
use crate::boojum::cs::Tool;

impl<'a, E: Engine, CS: ConstraintSystem<E> + 'static> CsBuilderImpl<GL, CsWrapperVerifierBuilder<'a, E, CS>> for CsWrapperVerifierBuilder<'a, E, CS> {
    type Final<GC: GateConfigurationHolder<GL>, TB: StaticToolboxHolder> = WrapperVerifier<E, CS>;

    type BuildParams<'b> = ();

    fn parameters<GC: GateConfigurationHolder<GL>, TB: StaticToolboxHolder>(builder: &CsBuilder<Self, GL, GC, TB>) -> CSGeometry {
        builder.implementation.parameters
    }

    fn allow_gate<GC: GateConfigurationHolder<GL>, TB: StaticToolboxHolder, G: Gate<GL>, TAux: 'static + Send + Sync + Clone>(
        mut builder: CsBuilder<Self, GL, GC, TB>,
        placement_strategy: GatePlacementStrategy,
        params: <<G as Gate<GL>>::Evaluator as GateConstraintEvaluator<GL>>::UniqueParameterizationParams,
        aux_data: TAux,
    ) -> CsBuilder<Self, GL, (GateTypeEntry<GL, G, TAux>, GC), TB> {
        // log!("Adding gate {:?}", std::any::type_name::<G>());

        let this = &mut builder.implementation;

        let new_configuration = builder.gates_config.add_gate::<G, _>(placement_strategy, params.clone(), aux_data);
        let evaluator_type_id = TypeId::of::<G::Evaluator>();
        let gate_type_id = TypeId::of::<G>();
        let evaluator = G::Evaluator::new_from_parameters(params.clone());

        // // depending on the configuration we should place it into corresponding set,
        // // and create some extra staff

        match placement_strategy {
            GatePlacementStrategy::UseGeneralPurposeColumns => {
                // we should batch gates that have the same evaluator
                if let Some(comparison_fns) = this.general_purpose_evaluators_comparison_functions.get_mut(&evaluator_type_id) {
                    let (dynamic_evaluator, comparator) = TypeErasedGateEvaluationWrapperVerificationFunction::from_evaluator(this.cs, evaluator, &this.parameters, placement_strategy);

                    let mut existing_idx = None;
                    for (other_comparator, idx) in comparison_fns.iter() {
                        if other_comparator.equals_to(&comparator) {
                            existing_idx = Some(*idx);
                            break;
                        }
                    }

                    if let Some(_existing_idx) = existing_idx {
                        // nothing, same evaluator
                    } else {
                        if comparison_fns.len() > 0 {
                            panic!("not yet supported");
                        }
                        let idx = this.evaluators_over_general_purpose_columns.len();
                        this.evaluators_over_general_purpose_columns.push(dynamic_evaluator);
                        // evaluator_type_id_into_evaluator_index_over_general_purpose_columns.insert(evaluator_type_id, idx);
                        comparison_fns.push((comparator, idx));
                    }

                    // gate_type_ids_for_general_purpose_columns.push(gate_type_id);
                } else {
                    // new one
                    let idx = this.evaluators_over_general_purpose_columns.len();
                    let (dynamic_evaluator, comparator) = TypeErasedGateEvaluationWrapperVerificationFunction::from_evaluator(this.cs, evaluator, &this.parameters, placement_strategy);
                    this.evaluators_over_general_purpose_columns.push(dynamic_evaluator);
                    // gate_type_ids_for_general_purpose_columns.push(gate_type_id);
                    // evaluator_type_id_into_evaluator_index_over_general_purpose_columns.insert(evaluator_type_id, idx);
                    this.general_purpose_evaluators_comparison_functions.insert(evaluator_type_id, vec![(comparator, idx)]);
                }
            }
            GatePlacementStrategy::UseSpecializedColumns { num_repetitions, share_constants } => {
                // we always add an evaluator

                let (dynamic_evaluator, _comparator) = TypeErasedGateEvaluationWrapperVerificationFunction::from_evaluator(this.cs, evaluator.clone(), &this.parameters, placement_strategy);

                // we need to extend copy-permutation data and witness placement data,
                // as well as keep track on offsets into them

                let _idx = this.evaluators_over_specialized_columns.len();
                this.gate_type_ids_for_specialized_columns.push(gate_type_id);
                this.evaluators_over_specialized_columns.push(dynamic_evaluator);
                // gate_type_id_into_evaluator_index_over_specialized_columns.insert(gate_type_id, idx);

                let principal_width = evaluator.instance_width();
                let mut total_width = principal_width;

                for _ in 1..num_repetitions {
                    total_width.num_variables += principal_width.num_variables;
                    total_width.num_witnesses += principal_width.num_witnesses;
                    if share_constants == false {
                        total_width.num_constants += principal_width.num_constants;
                    }
                }

                let total_constants_available = principal_width.num_constants;

                if share_constants {
                    match evaluator.placement_type() {
                        GatePlacementType::MultipleOnRow { per_chunk_offset: _ } => {}
                        GatePlacementType::UniqueOnRow => {
                            panic!("Can not share constants if placement type is unique");
                        }
                    }
                }

                let initial_offset = PerChunkOffset {
                    variables_offset: this.parameters.num_columns_under_copy_permutation + this.total_num_variables_for_specialized_columns,
                    witnesses_offset: this.parameters.num_witness_columns + this.total_num_witnesses_for_specialized_columns,
                    constants_offset: this.total_num_constants_for_specialized_columns, // we use separate vector for them
                };

                let offset_per_repetition = if share_constants == false {
                    PerChunkOffset {
                        variables_offset: principal_width.num_variables,
                        witnesses_offset: principal_width.num_witnesses,
                        constants_offset: principal_width.num_constants,
                    }
                } else {
                    let offset_per_repetition = match evaluator.placement_type() {
                        GatePlacementType::MultipleOnRow { per_chunk_offset } => per_chunk_offset,
                        GatePlacementType::UniqueOnRow => {
                            panic!("Can not share constants if placement type is unique");
                        }
                    };

                    assert_eq!(offset_per_repetition.variables_offset, principal_width.num_variables);
                    assert_eq!(offset_per_repetition.witnesses_offset, principal_width.num_witnesses);

                    // offset_per_repetition.variables_offset = principal_width.num_variables;
                    // offset_per_repetition.witnesses_offset = principal_width.num_witnesses;
                    // and we only leave constants untouched

                    offset_per_repetition
                };

                this.offsets_for_specialized_evaluators.push((initial_offset, offset_per_repetition, total_constants_available));

                this.total_num_variables_for_specialized_columns += total_width.num_variables;
                this.total_num_witnesses_for_specialized_columns += total_width.num_witnesses;
                this.total_num_constants_for_specialized_columns += total_width.num_constants;
            }
        }

        CsBuilder {
            gates_config: new_configuration,
            ..builder
        }
    }

    fn add_tool<GC: GateConfigurationHolder<GL>, TB: StaticToolboxHolder, M: 'static + Send + Sync + Clone, T: 'static + Send + Sync>(
        builder: CsBuilder<Self, GL, GC, TB>,
        tool: T,
        // ) -> CsBuilder<Self, F, GC, TB::DescendantHolder<M, T>> {
    ) -> CsBuilder<Self, GL, GC, (Tool<M, T>, TB)> {
        // TODO: toolbox is not used in the verifier, so perhaps it should be
        // moved out to the builder impl so it would not get in the way and just
        // hold the type in the phantom.
        let new_toolbox = builder.toolbox.add_tool(tool);

        CsBuilder { toolbox: new_toolbox, ..builder }
    }

    // type GcWithLookup<GC: GateConfigurationHolder<F>> = GC;
    type GcWithLookup<GC: GateConfigurationHolder<GL>> = (GateTypeEntry<GL, LookupFormalGate, LookupTooling>, GC);
    // GC::DescendantHolder<LookupFormalGate, LookupTooling>;

    fn allow_lookup<GC: GateConfigurationHolder<GL>, TB: StaticToolboxHolder>(
        builder: CsBuilder<Self, GL, GC, TB>,
        lookup_parameters: LookupParameters,
    ) -> CsBuilder<Self, GL, Self::GcWithLookup<GC>, TB> {
        let mut builder = builder;
        builder.implementation.lookup_parameters = lookup_parameters;

        let (placement_strategy, num_variables, num_constants, share_table_id) = match lookup_parameters {
            LookupParameters::NoLookup => {
                // this is formal

                (
                    GatePlacementStrategy::UseSpecializedColumns {
                        num_repetitions: 0,
                        share_constants: false,
                    },
                    0,
                    0,
                    false,
                )
            }
            LookupParameters::TableIdAsVariable { width, share_table_id } => {
                assert!(share_table_id == false, "other option is not yet implemented");
                // we need to resize multiplicities
                assert!(builder.implementation.parameters.num_columns_under_copy_permutation >= (width + 1) as usize);

                (GatePlacementStrategy::UseGeneralPurposeColumns, (width + 1) as usize, 0, share_table_id)
            }
            LookupParameters::TableIdAsConstant { width, share_table_id } => {
                assert!(share_table_id == true, "other option is not yet implemented");
                assert!(builder.implementation.parameters.num_columns_under_copy_permutation >= width as usize);

                (GatePlacementStrategy::UseGeneralPurposeColumns, width as usize, 1, share_table_id)
            }
            LookupParameters::UseSpecializedColumnsWithTableIdAsVariable {
                width,
                num_repetitions,
                share_table_id,
            } => {
                assert!(share_table_id == false, "other option is not yet implemented");

                (
                    GatePlacementStrategy::UseSpecializedColumns {
                        num_repetitions,
                        share_constants: false,
                    },
                    (width + 1) as usize,
                    0,
                    share_table_id,
                )
            }
            LookupParameters::UseSpecializedColumnsWithTableIdAsConstant {
                width,
                num_repetitions,
                share_table_id,
            } => {
                assert!(share_table_id == true, "other option is not yet implemented");

                (
                    GatePlacementStrategy::UseSpecializedColumns {
                        num_repetitions,
                        share_constants: share_table_id,
                    },
                    width as usize,
                    1,
                    share_table_id,
                )
            }
        };

        Self::allow_gate(builder, placement_strategy, (num_variables, num_constants, share_table_id), (Vec::with_capacity(32), 0))
    }

    fn build<'b, GC: GateConfigurationHolder<GL>, TB: StaticToolboxHolder, ARG: Into<Self::BuildParams<'b>>>(builder: CsBuilder<Self, GL, GC, TB>, _params: ARG) -> Self::Final<GC, TB> {
        let this: CsWrapperVerifierBuilder<E, CS> = builder.implementation;

        // capture small pieces of information from the gate configuration
        assert_eq!(this.evaluators_over_specialized_columns.len(), this.gate_type_ids_for_specialized_columns.len());

        let capacity = this.evaluators_over_specialized_columns.len();
        let mut placement_strategies = HashMap::with_capacity(capacity);

        for gate_type_id in this.gate_type_ids_for_specialized_columns.iter() {
            let placement_strategy = builder.gates_config.placement_strategy_for_type_id(*gate_type_id).expect("gate must be allowed");
            placement_strategies.insert(*gate_type_id, placement_strategy);
        }

        WrapperVerifier {
            parameters: this.parameters,
            lookup_parameters: this.lookup_parameters,
            gate_type_ids_for_specialized_columns: this.gate_type_ids_for_specialized_columns,
            evaluators_over_specialized_columns: this.evaluators_over_specialized_columns,
            offsets_for_specialized_evaluators: this.offsets_for_specialized_evaluators,
            evaluators_over_general_purpose_columns: this.evaluators_over_general_purpose_columns,
            total_num_variables_for_specialized_columns: this.total_num_variables_for_specialized_columns,
            total_num_witnesses_for_specialized_columns: this.total_num_witnesses_for_specialized_columns,
            total_num_constants_for_specialized_columns: this.total_num_constants_for_specialized_columns,
            placement_strategies,
        }
    }
}
