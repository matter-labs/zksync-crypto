use super::*;
use derivative::*;

use crate::boojum::cs::implementations::verifier::VerifierPolyStorage;
use crate::boojum::cs::implementations::verifier::VerifierRelationDestination;
use crate::boojum::cs::traits::evaluator::GatePlacementType;
use crate::boojum::cs::traits::evaluator::GatePurpose;
use crate::boojum::cs::traits::evaluator::GenericColumnwiseEvaluator;
use crate::boojum::cs::traits::evaluator::GenericDynamicEvaluatorOverGeneralPurposeColumns;
use crate::boojum::cs::traits::evaluator::GenericDynamicEvaluatorOverSpecializedColumns;
use crate::boojum::cs::traits::evaluator::GenericRowwiseEvaluator;

use crate::franklin_crypto::plonk::circuit::goldilocks::prime_field_like::GoldilocksExtAsFieldWrapper;

#[derive(Derivative)]
#[derivative(Debug)]
pub(crate) struct TypeErasedGateEvaluationWrapperVerificationFunction<E: Engine, CS: ConstraintSystem<E> + 'static> {
    pub(crate) debug_name: String,
    pub(crate) evaluator_type_id: TypeId,
    pub(crate) gate_purpose: GatePurpose,
    pub(crate) max_constraint_degree: usize,
    pub(crate) num_quotient_terms: usize,
    pub(crate) num_required_constants: usize,
    pub(crate) total_quotient_terms_over_all_repetitions: usize,
    pub(crate) num_repetitions_on_row: usize,
    pub(crate) placement_type: GatePlacementType,
    #[derivative(Debug = "ignore")]
    pub(crate) columnwise_satisfiability_function: Option<
        Box<
            dyn GenericDynamicEvaluatorOverSpecializedColumns<
                    GL,
                    GoldilocksExtAsFieldWrapper<E, CS>,
                    VerifierPolyStorage<GL, GoldilocksExtAsFieldWrapper<E, CS>>,
                    VerifierRelationDestination<GL, GoldilocksExtAsFieldWrapper<E, CS>>,
                >
                + 'static
                + Send
                + Sync,
        >,
    >,
    #[derivative(Debug = "ignore")]
    pub(crate) rowwise_satisfiability_function: Option<
        Box<
            dyn GenericDynamicEvaluatorOverGeneralPurposeColumns<
                    GL,
                    GoldilocksExtAsFieldWrapper<E, CS>,
                    VerifierPolyStorage<GL, GoldilocksExtAsFieldWrapper<E, CS>>,
                    VerifierRelationDestination<GL, GoldilocksExtAsFieldWrapper<E, CS>>,
                >
                + 'static
                + Send
                + Sync,
        >,
    >,
}

use crate::boojum::cs::traits::evaluator::GateBatchEvaluationComparisonFunction;
use crate::boojum::cs::traits::evaluator::GateConstraintEvaluator;
impl<E: Engine, CS: ConstraintSystem<E> + 'static> TypeErasedGateEvaluationWrapperVerificationFunction<E, CS> {
    pub fn from_evaluator<EV: GateConstraintEvaluator<GL>>(
        cs: &mut CS,
        evaluator: EV,
        geometry: &CSGeometry,
        placement_strategy: GatePlacementStrategy,
    ) -> (Self, GateBatchEvaluationComparisonFunction) {
        let debug_name = evaluator.instance_name();
        let evaluator_type_id = std::any::TypeId::of::<EV>();
        let gate_purpose = EV::gate_purpose();
        let max_constraint_degree = EV::max_constraint_degree();
        let num_quotient_terms = EV::num_quotient_terms();
        let num_required_constants = evaluator.num_required_constants_in_geometry(geometry);
        let placement_type = evaluator.placement_type();
        let mut final_per_chunk_offset = PerChunkOffset::zero();
        let (num_repetitions_on_row, total_quotient_terms_over_all_repetitions) = match placement_strategy {
            GatePlacementStrategy::UseGeneralPurposeColumns => {
                let num_repetitions_on_row = evaluator.num_repetitions_in_geometry(geometry);
                if let GatePlacementType::MultipleOnRow { per_chunk_offset } = &placement_type {
                    debug_assert!(num_repetitions_on_row > 0, "invalid for evaluator {}", std::any::type_name::<EV>());
                    final_per_chunk_offset = *per_chunk_offset;
                } else {
                    debug_assert_eq!(num_repetitions_on_row, 1);
                }

                let total_quotient_terms_in_geometry = evaluator.total_quotient_terms_in_geometry(geometry);

                (num_repetitions_on_row, total_quotient_terms_in_geometry)
            }
            GatePlacementStrategy::UseSpecializedColumns { num_repetitions, share_constants } => {
                let principal_width = evaluator.instance_width();
                final_per_chunk_offset = PerChunkOffset {
                    variables_offset: principal_width.num_variables,
                    witnesses_offset: principal_width.num_witnesses,
                    constants_offset: principal_width.num_constants,
                };
                if share_constants {
                    final_per_chunk_offset.constants_offset = 0;
                }

                (num_repetitions, num_repetitions * num_quotient_terms)
            }
        };

        let (specialized_satisfiability_evaluator, general_purpose_satisfiability_evaluator) = match placement_strategy {
            GatePlacementStrategy::UseSpecializedColumns { .. } => {
                let specialized_evaluator = GenericColumnwiseEvaluator {
                    evaluator: evaluator.clone(),
                    global_constants: evaluator.create_global_constants(cs),
                    num_repetitions: num_repetitions_on_row,
                    per_chunk_offset: final_per_chunk_offset,
                };

                // dbg!(&specialized_evaluator);

                (
                    Some(Box::new(specialized_evaluator)
                        as Box<
                            dyn GenericDynamicEvaluatorOverSpecializedColumns<
                                    GL,
                                    GoldilocksExtAsFieldWrapper<E, CS>,
                                    VerifierPolyStorage<GL, GoldilocksExtAsFieldWrapper<E, CS>>,
                                    VerifierRelationDestination<GL, GoldilocksExtAsFieldWrapper<E, CS>>,
                                >
                                + 'static
                                + Send
                                + Sync,
                        >),
                    None,
                )
            }
            GatePlacementStrategy::UseGeneralPurposeColumns => {
                let general_purpose_evaluator = GenericRowwiseEvaluator {
                    evaluator: evaluator.clone(),
                    global_constants: evaluator.create_global_constants(cs),
                    num_repetitions: num_repetitions_on_row,
                    per_chunk_offset: final_per_chunk_offset,
                };

                // dbg!(&general_purpose_evaluator);

                (
                    None,
                    Some(Box::new(general_purpose_evaluator)
                        as Box<
                            dyn GenericDynamicEvaluatorOverGeneralPurposeColumns<
                                    GL,
                                    GoldilocksExtAsFieldWrapper<E, CS>,
                                    VerifierPolyStorage<GL, GoldilocksExtAsFieldWrapper<E, CS>>,
                                    VerifierRelationDestination<GL, GoldilocksExtAsFieldWrapper<E, CS>>,
                                >
                                + 'static
                                + Send
                                + Sync,
                        >),
                )
            }
        };

        let this_params = evaluator.unique_params();

        let comparison_fn = move |other_evaluator: &dyn std::any::Any| -> bool {
            assert_eq!(other_evaluator.type_id(), evaluator_type_id);
            let other_evaluator: &EV = other_evaluator.downcast_ref().expect("must downcast");

            this_params == other_evaluator.unique_params()
        };

        let comparator = GateBatchEvaluationComparisonFunction {
            type_id: evaluator_type_id,
            evaluator_dyn: Box::new(evaluator),
            equality_fn: Box::new(comparison_fn),
        };

        let new = Self {
            debug_name,
            evaluator_type_id,
            gate_purpose,
            max_constraint_degree,
            num_quotient_terms,
            num_required_constants,
            total_quotient_terms_over_all_repetitions,
            num_repetitions_on_row,
            placement_type,
            columnwise_satisfiability_function: specialized_satisfiability_evaluator,
            rowwise_satisfiability_function: general_purpose_satisfiability_evaluator,
        };

        (new, comparator)
    }
}
