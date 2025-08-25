use crate::boojum::cs::implementations::proof::Proof;
use crate::boojum::cs::implementations::prover::ProofConfig;
use crate::boojum::cs::implementations::utils::domain_generator_for_size;
use crate::boojum::cs::implementations::verifier::VerificationKey;
use crate::boojum::cs::implementations::verifier::VerificationKeyCircuitGeometry;
use crate::boojum::cs::oracle::TreeHasher;
use crate::boojum::cs::LookupParameters;
use crate::boojum::field::goldilocks::{GoldilocksExt2 as GLExt2, GoldilocksField as GL};
use crate::boojum::field::traits::field_like::PrimeFieldLike;
use crate::boojum::field::Field as BoojumField;
use crate::boojum::field::PrimeField as BoojumPrimeField;

use crate::franklin_crypto::bellman::pairing::Engine;
use crate::franklin_crypto::bellman::plonk::better_better_cs::cs::ConstraintSystem;
use crate::franklin_crypto::bellman::plonk::better_better_cs::cs::*;
use crate::franklin_crypto::bellman::plonk::better_better_cs::gates::selector_optimized_with_d_next::SelectorOptimizedWidth4MainGateWithDNext;
use crate::franklin_crypto::bellman::{Field, PrimeField, PrimeFieldRepr, SynthesisError};
use crate::franklin_crypto::plonk::circuit::allocated_num::{AllocatedNum, Num};
use crate::franklin_crypto::plonk::circuit::bigint_new::BITWISE_LOGICAL_OPS_TABLE_NAME;
use crate::franklin_crypto::plonk::circuit::boolean::Boolean;
use crate::franklin_crypto::plonk::circuit::custom_rescue_gate::Rescue5CustomGate;
use crate::franklin_crypto::plonk::circuit::goldilocks::prime_field_like::{GoldilocksAsFieldWrapper, GoldilocksExtAsFieldWrapper};
use crate::franklin_crypto::plonk::circuit::goldilocks::GoldilocksField;
use crate::franklin_crypto::plonk::circuit::linear_combination::LinearCombination;
use crate::franklin_crypto::plonk::circuit::Assignment;

use crate::implementations::poseidon2::pow::ConcretePoseidon2SpongeGadget;
use crate::traits::circuit::*;
use crate::traits::transcript::CircuitGLTranscript;
use crate::traits::tree_hasher::CircuitGLTreeHasher;
use crate::traits::*;
use crate::verifier_structs::allocated_vk::AllocatedVerificationKey;
use crate::verifier_structs::challenges::{ChallengesHolder, EvaluationsHolder};
use crate::verifier_structs::constants::ConstantsHolder;
use crate::verifier_structs::{allocated_proof::*, *};

mod first_step;
mod fri;
mod quotient_contributions;
pub(crate) mod utils;

use first_step::*;
use fri::*;
use pow::RecursivePoWRunner;
use quotient_contributions::*;
use utils::*;

#[derive(Clone, Debug, serde::Serialize)]
pub struct WrapperCircuit<
    E: Engine,
    HS: TreeHasher<GL, Output = E::Fr>,
    H: CircuitGLTreeHasher<E, CircuitOutput = Num<E>, NonCircuitSimulator = HS>,
    TR: CircuitGLTranscript<E, CircuitCompatibleCap = H::CircuitOutput>,
    PWF: ProofWrapperFunction<E>,
> {
    pub witness: Option<Proof<GL, HS, GLExt2>>,
    pub vk: VerificationKey<GL, H::NonCircuitSimulator>,
    pub fixed_parameters: VerificationKeyCircuitGeometry,
    pub transcript_params: TR::TranscriptParameters,
    pub wrapper_function: PWF,
}

impl<
        E: Engine,
        HS: TreeHasher<GL, Output = E::Fr>,
        H: CircuitGLTreeHasher<E, CircuitOutput = Num<E>, NonCircuitSimulator = HS>,
        TR: CircuitGLTranscript<E, CircuitCompatibleCap = H::CircuitOutput>,
        PWF: ProofWrapperFunction<E>,
    > Circuit<E> for WrapperCircuit<E, HS, H, TR, PWF>
{
    type MainGate = SelectorOptimizedWidth4MainGateWithDNext;

    fn declare_used_gates() -> Result<Vec<Box<dyn GateInternal<E>>>, SynthesisError> {
        Ok(vec![Self::MainGate::default().into_internal(), Rescue5CustomGate::default().into_internal()])
    }

    fn synthesize<CS: ConstraintSystem<E> + 'static>(&self, cs: &mut CS) -> Result<(), SynthesisError> {
        let now = std::time::Instant::now();
        // Add table for range check
        let columns3 = vec![PolyIdentifier::VariablesPolynomial(0), PolyIdentifier::VariablesPolynomial(1), PolyIdentifier::VariablesPolynomial(2)];

        let name = BITWISE_LOGICAL_OPS_TABLE_NAME;
        let bitwise_logic_table = LookupTableApplication::new(name, TwoKeysOneValueBinopTable::<E, XorBinop>::new(8, name), columns3.clone(), None, true);
        println!("Table created in {:?}", now.elapsed());
        let now = std::time::Instant::now();
        cs.add_table(bitwise_logic_table).unwrap();
        println!("Table added in {:?}", now.elapsed());
        let now = std::time::Instant::now();

        // Prepare for proof verification
        let verifier_builder = self.wrapper_function.builder_for_wrapper();
        let verifier = verifier_builder.create_wrapper_verifier(cs);

        let proof_config = self.wrapper_function.proof_config_for_compression_step();
        let fixed_parameters = self.fixed_parameters.clone();

        let vk = AllocatedVerificationKey::<E, H>::allocate_constant(&self.vk, &fixed_parameters);
        println!("Constant allocation took {:?}", now.elapsed());
        let now = std::time::Instant::now();

        let proof: AllocatedProof<E, H> = AllocatedProof::allocate_from_witness(cs, &self.witness, &verifier, &fixed_parameters, &proof_config)?;

        println!("Proof from witness took {:?}", now.elapsed());
        let now = std::time::Instant::now();

        // Verify proof
        let correct = crate::verifier::verify::<E, CS, H, TR, ConcretePoseidon2SpongeGadget<E>>(cs, self.transcript_params.clone(), &proof_config, &proof, &verifier, &fixed_parameters, &vk)?;
        Boolean::enforce_equal(cs, &correct, &Boolean::constant(true))?;
        println!("proof verify took {:?}", now.elapsed());
        let now = std::time::Instant::now();

        // Aggregate PI
        let _pi = aggregate_public_inputs(cs, &proof.public_inputs)?;
        println!("PI aggregation took {:?}", now.elapsed());

        Ok(())
    }
}

#[derive(Clone)]
pub struct WrapperCircuitWidth3NoLookupNoCustomGate<
    E: Engine,
    HS: TreeHasher<GL, Output = E::Fr>,
    H: CircuitGLTreeHasher<E, CircuitOutput = Num<E>, NonCircuitSimulator = HS>,
    TR: CircuitGLTranscript<E, CircuitCompatibleCap = H::CircuitOutput>,
    PWF: ProofWrapperFunction<E>,
> {
    pub witness: Option<Proof<GL, HS, GLExt2>>,
    pub vk: VerificationKey<GL, H::NonCircuitSimulator>,
    pub fixed_parameters: VerificationKeyCircuitGeometry,
    pub transcript_params: TR::TranscriptParameters,
    pub wrapper_function: PWF,
}

impl<
        E: Engine,
        HS: TreeHasher<GL, Output = E::Fr>,
        H: CircuitGLTreeHasher<E, CircuitOutput = Num<E>, NonCircuitSimulator = HS>,
        TR: CircuitGLTranscript<E, CircuitCompatibleCap = H::CircuitOutput>,
        PWF: ProofWrapperFunction<E>,
    > Circuit<E> for WrapperCircuitWidth3NoLookupNoCustomGate<E, HS, H, TR, PWF>
{
    type MainGate = rescue_poseidon::franklin_crypto::bellman::plonk::better_better_cs::gates::naive_main_gate::NaiveMainGate;

    fn declare_used_gates() -> Result<Vec<Box<dyn GateInternal<E>>>, SynthesisError> {
        Ok(vec![Self::MainGate::default().into_internal()])
    }

    fn synthesize<CS: ConstraintSystem<E> + 'static>(&self, cs: &mut CS) -> Result<(), SynthesisError> {
        // Prepare for proof verification
        let verifier_builder = self.wrapper_function.builder_for_wrapper();
        let verifier = verifier_builder.create_wrapper_verifier(cs);

        let proof_config = self.wrapper_function.proof_config_for_compression_step();
        let fixed_parameters = self.fixed_parameters.clone();

        let vk = AllocatedVerificationKey::<E, H>::allocate_constant(&self.vk, &fixed_parameters);
        let proof: AllocatedProof<E, H> = AllocatedProof::allocate_from_witness(cs, &self.witness, &verifier, &fixed_parameters, &proof_config)?;
        // Verify proof
        let correct = crate::verifier::verify::<E, CS, H, TR, ConcretePoseidon2SpongeGadget<E>>(cs, self.transcript_params.clone(), &proof_config, &proof, &verifier, &fixed_parameters, &vk)?;
        Boolean::enforce_equal(cs, &correct, &Boolean::constant(true))?;

        // Aggregate PI
        let _pi = aggregate_public_inputs(cs, &proof.public_inputs)?;

        Ok(())
    }
}

pub fn verify<E: Engine, CS: ConstraintSystem<E> + 'static, H: CircuitGLTreeHasher<E>, TR: CircuitGLTranscript<E, CircuitCompatibleCap = H::CircuitOutput>, POW: RecursivePoWRunner<E>>(
    cs: &mut CS,
    transcript_params: TR::TranscriptParameters,
    proof_config: &ProofConfig,
    proof: &AllocatedProof<E, H>,
    verifier: &WrapperVerifier<E, CS>,
    fixed_parameters: &VerificationKeyCircuitGeometry,
    vk: &AllocatedVerificationKey<E, H>,
) -> Result<Boolean, SynthesisError> {
    let now = std::time::Instant::now();
    let mut validity_flags = Vec::with_capacity(256);

    let mut transcript = TR::new(cs, transcript_params)?;
    let mut challenges = ChallengesHolder::new(cs);

    println!("Before constants took {:?}", now.elapsed());
    let now = std::time::Instant::now();
    // prepare constants
    let constants = ConstantsHolder::generate(proof_config, verifier, fixed_parameters);
    assert_eq!(fixed_parameters.cap_size, vk.setup_merkle_tree_cap.len());

    println!("Constants took {:?}", now.elapsed());
    let now = std::time::Instant::now();

    let public_input_opening_tuples = verify_first_step(cs, proof, vk, &mut challenges, &mut transcript, verifier, fixed_parameters, &constants)?;
    println!("First step took {:?}", now.elapsed());
    let now = std::time::Instant::now();

    validity_flags.extend(check_quotient_contributions_in_z(cs, proof, &challenges, verifier, fixed_parameters, &constants)?);
    println!("Quotient contributions took {:?}", now.elapsed());
    let now = std::time::Instant::now();
    validity_flags.extend(verify_fri_part::<E, CS, H, TR, POW>(
        cs,
        proof,
        vk,
        &mut challenges,
        &mut transcript,
        public_input_opening_tuples,
        verifier,
        fixed_parameters,
        &constants,
    )?);
    println!("Verify fri took {:?}", now.elapsed());

    let correct = smart_and(cs, &validity_flags)?;

    Ok(correct)
}

/// aggregate public inputs to one scalar field element
fn aggregate_public_inputs<E: Engine, CS: ConstraintSystem<E>>(cs: &mut CS, public_inputs: &[GoldilocksField<E>]) -> Result<Num<E>, SynthesisError> {
    let chunk_bit_size = (GL::CAPACITY_BITS / 8) * 8;
    assert!(
        public_inputs.len() * chunk_bit_size <= E::Fr::CAPACITY as usize,
        "scalar field capacity is not enough to fit all public inputs"
    );

    // Firstly we check that public inputs have correct size
    for pi in public_inputs.iter() {
        if let Ok(_) = cs.get_table(BITWISE_LOGICAL_OPS_TABLE_NAME) {
            range_check_with_lookup(cs, &pi.into_num(), chunk_bit_size)?;
        } else {
            range_check_with_naive(cs, &pi.into_num(), chunk_bit_size)?;
        }
    }

    // compute aggregated pi value
    let mut tmp = E::Fr::one();
    let mut shift_repr = <E::Fr as PrimeField>::Repr::from(1);
    shift_repr.shl(chunk_bit_size as u32);
    let shift = E::Fr::from_repr(shift_repr).unwrap();

    let mut lc = LinearCombination::<E>::zero();
    for pi in public_inputs.iter().rev() {
        lc.add_assign_number_with_coeff(&pi.into_num(), tmp);
        tmp.mul_assign(&shift);
    }

    // allocate as pi
    let pi = Num::Variable(AllocatedNum::alloc_input(cs, || Ok(*lc.get_value().get()?))?);

    // check sum
    let mut minus_one = E::Fr::one();
    minus_one.negate();
    lc.add_assign_number_with_coeff(&pi, minus_one);
    lc.enforce_zero(cs)?;

    Ok(pi)
}

pub fn range_check_with_naive<E: Engine, CS: ConstraintSystem<E>>(cs: &mut CS, num: &Num<E>, num_bits: usize) -> Result<(), SynthesisError> {
    use rescue_poseidon::franklin_crypto::plonk::circuit::goldilocks::range_check_for_num_bits;
    range_check_for_num_bits(cs, num, num_bits)?;

    Ok(())
}

pub fn range_check_with_lookup<E: Engine, CS: ConstraintSystem<E>>(cs: &mut CS, num: &Num<E>, num_bits: usize) -> Result<(), SynthesisError> {
    let table = cs.get_table(BITWISE_LOGICAL_OPS_TABLE_NAME).unwrap();
    use rescue_poseidon::franklin_crypto::plonk::circuit::bigint_new::enforce_range_check_using_bitop_table;
    enforce_range_check_using_bitop_table(cs, &num.get_variable(), num_bits, table, false)?;
    Ok(())
}
