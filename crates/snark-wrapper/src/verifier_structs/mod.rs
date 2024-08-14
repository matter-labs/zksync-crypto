use crate::franklin_crypto::plonk::circuit::goldilocks::GoldilocksField;

use crate::boojum::cs::implementations::proof::Proof;
use crate::boojum::cs::implementations::proof::{OracleQuery, SingleRoundQueries};
use crate::boojum::cs::implementations::prover::ProofConfig;
use crate::boojum::cs::implementations::verifier::VerificationKeyCircuitGeometry;
use crate::boojum::cs::oracle::TreeHasher;
use crate::boojum::cs::traits::evaluator::PerChunkOffset;
use crate::boojum::cs::traits::gate::GatePlacementStrategy;
use crate::boojum::cs::CSGeometry;
use crate::boojum::cs::LookupParameters;
use crate::boojum::field::goldilocks::{GoldilocksExt2 as GLExt2, GoldilocksField as GL};

use crate::franklin_crypto::bellman::pairing::Engine;
use crate::franklin_crypto::bellman::plonk::better_better_cs::cs::ConstraintSystem;
use crate::franklin_crypto::bellman::SynthesisError;
use crate::franklin_crypto::plonk::circuit::allocated_num::Num;
use crate::franklin_crypto::plonk::circuit::boolean::Boolean;

use crate::traits::tree_hasher::CircuitGLTreeHasher;
use crate::verifier_structs::constants::ConstantsHolder;
use crate::verifier_structs::gate_evaluator::TypeErasedGateEvaluationWrapperVerificationFunction;

use std::any::TypeId;
use std::collections::HashMap;

pub mod allocated_proof;
pub mod allocated_queries;
pub mod allocated_vk;
pub mod challenges;
pub mod constants;
pub mod gate_evaluator;

pub struct WrapperVerifier<E: Engine, CS: ConstraintSystem<E> + 'static> {
    // when we init we get the following from VK
    pub parameters: CSGeometry,
    pub lookup_parameters: LookupParameters,

    pub(crate) gate_type_ids_for_specialized_columns: Vec<TypeId>,
    pub(crate) evaluators_over_specialized_columns: Vec<TypeErasedGateEvaluationWrapperVerificationFunction<E, CS>>,
    pub(crate) offsets_for_specialized_evaluators: Vec<(PerChunkOffset, PerChunkOffset, usize)>,

    pub(crate) evaluators_over_general_purpose_columns: Vec<TypeErasedGateEvaluationWrapperVerificationFunction<E, CS>>,

    pub(crate) total_num_variables_for_specialized_columns: usize,
    pub(crate) total_num_witnesses_for_specialized_columns: usize,
    pub(crate) total_num_constants_for_specialized_columns: usize,

    pub(crate) placement_strategies: HashMap<TypeId, GatePlacementStrategy>,
}

pub fn allocate_num_elements<T, R, E: Engine, CS: ConstraintSystem<E>>(
    cs: &mut CS,
    num_elements: usize,
    mut source: Option<impl Iterator<Item = T>>,
    allocating_function: impl Fn(&mut CS, Option<T>) -> Result<R, SynthesisError>,
) -> Result<Vec<R>, SynthesisError> {
    let mut result = Vec::with_capacity(num_elements);

    for _ in 0..num_elements {
        let el = source.as_mut().map(|el| el.next().expect("Should be enough elements in the source"));
        result.push(allocating_function(cs, el)?);
    }
    debug_assert!(source.as_mut().map(|el| el.next().is_none()).unwrap_or(true));

    Ok(result)
}
