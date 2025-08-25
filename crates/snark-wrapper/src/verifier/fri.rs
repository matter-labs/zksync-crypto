use std::collections::HashMap;
use std::os::unix::thread;

use super::*;

use crate::traits::pow::RecursivePoWRunner;
use crate::traits::transcript::BoolsBuffer;
use crate::verifier_structs::allocated_queries::AllocatedSingleRoundQueries;
use crossbeam::channel::{unbounded, Sender};
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};

#[derive(Clone, Debug)]
pub struct GatesInfo<E: Engine> {
    pub coefficient_assignments: Vec<E::Fr>,
    pub variable_assignments: Vec<Variable>,
    pub witness_assignments: Vec<E::Fr>,
    pub gate_type: GateType,
}

#[derive(Clone, Debug)]
pub enum GateType {
    Selector4,
    Rescue5Custom,
}

impl GateType {
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "main gate of width 4 with D_next and selector optimization" => Some(GateType::Selector4),
            "Alpha=5 custom gate for Rescue/Poseidon" => Some(GateType::Rescue5Custom),
            _ => None,
        }
    }
}

#[derive(Clone)]
pub struct MyCS<E: Engine> {
    pub num_variables: usize,
    pub started_variable: usize,
    pub num_aux_gates: usize,
    pub aux_assingments: Vec<E::Fr>,
    pub gates: Vec<GatesInfo<E>>,
    pub old_nums: HashMap<Variable, E::Fr>,
}

impl<E: Engine> Default for MyCS<E> {
    fn default() -> Self {
        Self {
            num_variables: 0,
            started_variable: 0,
            num_aux_gates: 0,
            aux_assingments: Vec::new(),
            gates: Vec::new(),
            old_nums: HashMap::new(),
        }
    }
}

impl<E: Engine> MyCS<E> {
    pub fn new(num_variables: usize) -> Self {
        Self {
            num_variables,
            started_variable: num_variables,
            num_aux_gates: 0,
            aux_assingments: Vec::new(),
            gates: Vec::new(),
            old_nums: HashMap::new(),
        }
    }

    pub fn add_allocated_num(&mut self, allocated_num: &AllocatedNum<E>) {
        self.old_nums.insert(allocated_num.get_variable(), allocated_num.get_value().unwrap());
    }

    pub fn dump_to_existing_cs<CS: ConstraintSystem<E>>(&self, other_cs: &mut CS, start_position: (usize, usize)) {
        // move pointer to the proper location.
        let current_pos = other_cs.reposition(start_position.0, start_position.1).unwrap();

        for (i, entry) in self.aux_assingments.iter().enumerate() {
            let tmp = other_cs.alloc(|| Ok(*entry)).unwrap();
            assert_eq!(tmp, Variable::new_unchecked(Index::Aux(i + self.started_variable + 1)));
        }
        for entry in &self.gates {
            match entry.gate_type {
                GateType::Selector4 => {
                    let tmp = SelectorOptimizedWidth4MainGateWithDNext;
                    other_cs
                        .new_single_gate_for_trace_step(&tmp, &entry.coefficient_assignments, &entry.variable_assignments, &entry.witness_assignments)
                        .unwrap();
                }
                GateType::Rescue5Custom => {
                    let tmp = Rescue5CustomGate;
                    other_cs
                        .new_single_gate_for_trace_step(&tmp, &entry.coefficient_assignments, &entry.variable_assignments, &entry.witness_assignments)
                        .unwrap();
                }
            }
        }
        // move pointers back.
        other_cs.reposition(current_pos.0, current_pos.1).unwrap();
    }
}

impl<E: Engine> ConstraintSystem<E> for MyCS<E> {
    type Params = PlonkCsWidth4WithNextStepAndCustomGatesParams;

    type MainGate = SelectorOptimizedWidth4MainGateWithDNext;

    fn get_value(&self, variable: Variable) -> Result<E::Fr, SynthesisError> {
        let previous = self.old_nums.get(&variable);
        if let Some(value) = previous {
            return Ok(*value);
        }
        match variable.get_unchecked() {
            Index::Aux(index) => {
                if index <= self.started_variable {
                    panic!("aux variable index {} is less than started_variable {}", index, self.started_variable);
                }

                return Ok(self.aux_assingments[index - self.started_variable - 1]);
            }
            Index::Input(index) => panic!("input not supported"),
        }
    }

    /*fn allocate_main_gate(&mut self, term: MainGateTerm<E>) -> Result<(), SynthesisError> {
        todo!();
    }*/

    fn get_current_aux_assignment_number(&self) -> usize {
        todo!();
    }

    fn reserve(&mut self, aux_variables: usize, aux_gates: usize) -> Result<(usize, usize), SynthesisError> {
        todo!()
    }
    fn reposition(&mut self, aux_variable_pos: usize, aux_gates_pos: usize) -> Result<(usize, usize), SynthesisError> {
        todo!()
    }

    fn alloc<F>(&mut self, value: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>,
    {
        self.num_variables += 1;

        let value = value().unwrap();
        self.aux_assingments.push(value);

        Ok(Variable::new_unchecked(Index::Aux(self.num_variables)))
    }

    fn alloc_input<F>(&mut self, value: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>,
    {
        todo!()
    }

    fn get_main_gate(&self) -> &Self::MainGate {
        todo!()
    }

    fn begin_gates_batch_for_step(&mut self) -> Result<(), SynthesisError> {
        self.num_aux_gates += 1;
        Ok(())
    }

    fn new_gate_in_batch<G: Gate<E>>(&mut self, equation: &G, coefficients_assignments: &[E::Fr], variables_assignments: &[Variable], witness_assignments: &[E::Fr]) -> Result<(), SynthesisError> {
        let gate_info = GatesInfo {
            coefficient_assignments: coefficients_assignments.to_vec(),
            variable_assignments: variables_assignments.to_vec(),
            witness_assignments: witness_assignments.to_vec(),
            gate_type: GateType::from_str(equation.name()).expect(&format!("Failed to detect gate {} ", equation.name())),
        };
        self.gates.push(gate_info);

        Ok(())
    }

    fn end_gates_batch_for_step(&mut self) -> Result<(), SynthesisError> {
        Ok(())
    }

    fn allocate_variables_without_gate(&mut self, variables_assignments: &[Variable], witness_assignments: &[E::Fr]) -> Result<(), SynthesisError> {
        todo!()
    }

    fn get_dummy_variable() -> Variable {
        Variable::new_unchecked(Index::Aux(0))
    }

    fn get_explicit_zero(&mut self) -> Result<Variable, SynthesisError> {
        todo!()
    }

    fn get_explicit_one(&mut self) -> Result<Variable, SynthesisError> {
        todo!()
    }

    fn add_table(&mut self, table: LookupTableApplication<E>) -> Result<std::sync::Arc<LookupTableApplication<E>>, SynthesisError> {
        todo!()
    }

    fn get_table(&self, functional_name: &str) -> Result<std::sync::Arc<LookupTableApplication<E>>, SynthesisError> {
        todo!()
    }

    fn add_multitable(&mut self, table: MultiTableApplication<E>) -> Result<(), SynthesisError> {
        todo!()
    }

    fn get_multitable(&self, functional_name: &str) -> Result<std::sync::Arc<MultiTableApplication<E>>, SynthesisError> {
        todo!()
    }

    fn apply_single_lookup_gate(&mut self, variables: &[Variable], gate: std::sync::Arc<LookupTableApplication<E>>) -> Result<(), SynthesisError> {
        todo!()
    }

    fn apply_multi_lookup_gate(&mut self, variables: &[Variable], gate: std::sync::Arc<MultiTableApplication<E>>) -> Result<(), SynthesisError> {
        todo!()
    }

    fn get_current_step_number(&self) -> usize {
        todo!()
    }

    fn get_current_aux_gate_number(&self) -> usize {
        todo!()
    }
}

pub(crate) fn verify_fri_part<
    E: Engine,
    CS: ConstraintSystem<E> + 'static,
    H: CircuitGLTreeHasher<E>,
    TR: CircuitGLTranscript<E, CircuitCompatibleCap = H::CircuitOutput>,
    POW: RecursivePoWRunner<E>,
>(
    cs: &mut CS,
    proof: &AllocatedProof<E, H>,
    vk: &AllocatedVerificationKey<E, H>,
    challenges: &mut ChallengesHolder<E, CS>,
    transcript: &mut TR,
    public_input_opening_tuples: Vec<(GL, Vec<(usize, GoldilocksAsFieldWrapper<E, CS>)>)>,
    // parameters
    verifier: &WrapperVerifier<E, CS>,
    fixed_parameters: &VerificationKeyCircuitGeometry,
    constants: &ConstantsHolder,
) -> Result<Vec<Boolean>, SynthesisError> {
    let thread_pool = ThreadPoolBuilder::new().num_threads(32).build().unwrap();
    let now = std::time::Instant::now();
    let mut validity_flags = vec![];

    // get challenges
    challenges.get_challenges_for_fri_quotiening(cs, transcript, constants.total_num_challenges_for_fri_quotiening)?;
    challenges.get_fri_intermediate_challenges(cs, transcript, proof, fixed_parameters, constants)?;

    assert_eq!(constants.final_expected_degree as usize, proof.final_fri_monomials[0].len());
    assert_eq!(constants.final_expected_degree as usize, proof.final_fri_monomials[1].len());
    assert!(proof.final_fri_monomials[0].len() > 0);

    // witness monomial coeffs
    transcript.witness_field_elements(cs, &proof.final_fri_monomials[0])?;
    transcript.witness_field_elements(cs, &proof.final_fri_monomials[1])?;

    println!("  Witness field elements took {:?}", now.elapsed());
    let now = std::time::Instant::now();

    if constants.new_pow_bits != 0 {
        const SEED_BITS: usize = 256;
        // pull enough challenges from the transcript
        let num_challenges = SEED_BITS.next_multiple_of(GL::CHAR_BITS) / GL::CHAR_BITS;
        let challenges: Vec<_> = transcript.get_multiple_challenges(cs, num_challenges as usize)?;
        let (is_valid, pow_challenge_limbs) = POW::verify_from_field_elements(cs, challenges, proof.pow_challenge_le, constants.new_pow_bits)?;
        match is_valid.get_value() {
            Some(is_valid) => {
                if is_valid == false {
                    println!("PoW challenge is invalid")
                }
            }
            None => (),
        }
        transcript.witness_field_elements(cs, &pow_challenge_limbs)?;
    }
    println!("  PoW challenge limbs took {:?}", now.elapsed());
    let now = std::time::Instant::now();

    let max_needed_bits = (fixed_parameters.domain_size * fixed_parameters.fri_lde_factor as u64).trailing_zeros() as usize;

    let mut bools_buffer = BoolsBuffer {
        available: vec![],
        max_needed: max_needed_bits,
    };

    let num_bits_for_in_coset_index = max_needed_bits - fixed_parameters.fri_lde_factor.trailing_zeros() as usize;
    let base_tree_index_shift = fixed_parameters.domain_size.trailing_zeros();
    assert_eq!(num_bits_for_in_coset_index, base_tree_index_shift as usize);

    assert_eq!(constants.num_fri_repetitions, proof.queries_per_fri_repetition.len());

    let multiplicative_generator = GoldilocksAsFieldWrapper::constant(GL::multiplicative_generator(), cs);

    // precompute once, will be handy later
    let mut precomputed_powers = vec![];
    let mut precomputed_powers_inversed = vec![];
    for i in 0..=(fixed_parameters.domain_size * fixed_parameters.fri_lde_factor as u64).trailing_zeros() {
        let omega = domain_generator_for_size::<GL>(1u64 << i);
        precomputed_powers.push(omega);
        precomputed_powers_inversed.push(BoojumPrimeField::inverse(&omega).unwrap());
    }

    println!("  Precomputation took {:?}", now.elapsed());
    let now = std::time::Instant::now();

    // we also want to precompute "steps" for different interpolation degrees
    // e.g. if we interpolate 8 elements,
    // then those will be ordered as bitreverses of [0..=7], namely
    // [0, 4, 2, 6, 1, 5, 3, 7]

    // so we want to have exactly half of it, because separation by 4
    // is exactly -1, so we need [1, sqrt4(1), sqrt8(1), sqrt4(1)*sqrt8(1)]

    let mut interpolation_steps = vec![GL::ONE; 4]; // max size

    for idx in [1, 3].into_iter() {
        BoojumField::mul_assign(&mut interpolation_steps[idx], &precomputed_powers_inversed[2]);
    }
    for idx in [2, 3].into_iter() {
        BoojumField::mul_assign(&mut interpolation_steps[idx], &precomputed_powers_inversed[3]);
    }

    assert_eq!(interpolation_steps[0], GL::ONE);
    assert_eq!(BoojumField::pow_u64(&interpolation_steps[1], 4), GL::ONE);
    assert_eq!(BoojumField::pow_u64(&interpolation_steps[2], 8), GL::ONE);

    let precomputed_powers: Vec<_> = precomputed_powers.into_iter().map(|el| GoldilocksAsFieldWrapper::constant(el, cs)).collect();
    let precomputed_powers_inversed: Vec<_> = precomputed_powers_inversed.into_iter().map(|el| GoldilocksAsFieldWrapper::constant(el, cs)).collect();
    let interpolation_steps: Vec<_> = interpolation_steps.into_iter().map(|el| GoldilocksAsFieldWrapper::constant(el, cs)).collect();

    let base_oracle_depth = fixed_parameters.base_oracles_depth();

    println!("  Before queries took {:?}. Queries count: {}", now.elapsed(), proof.queries_per_fri_repetition.len());
    let now = std::time::Instant::now();

    let (tx, rx) = unbounded::<InclusionResponse<E>>();

    for (idx, queries) in proof.queries_per_fri_repetition.iter().enumerate() {
        let now = std::time::Instant::now();
        let query_index_lsb_first_bits = bools_buffer.get_bits(cs, transcript, max_needed_bits)?;

        // we consider it to be some convenient for us encoding of coset + inner index.

        // Small note on indexing: when we commit to elements we use bitreversal enumeration everywhere.
        // So index `i` in the tree corresponds to the element of `omega^bitreverse(i)`.
        // This gives us natural separation of LDE cosets, such that subtrees form independent cosets,
        // and if cosets are in the form of `{1, gamma, ...} x {1, omega, ...} where gamma^lde_factor == omega,
        // then subtrees are enumerated by bitreverse powers of gamma

        // let inner_idx = &query_index_lsb_first_bits[0..num_bits_for_in_coset_index];
        // let coset_idx = &query_index_lsb_first_bits[num_bits_for_in_coset_index..];
        let base_tree_idx = query_index_lsb_first_bits.clone();

        if idx < 10 {
            println!("    - {} just before proof verification {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();
        // first verify basic inclusion proofs
        verify_inclusion_proofs(&thread_pool, cs, queries, proof, vk, &base_tree_idx, constants, base_oracle_depth, validity_flags.len(), tx.clone())?;

        let tmp: Vec<Boolean> = vec![Default::default(); 4];
        validity_flags.extend(tmp);

        if idx < 10 {
            println!("    - {} inclusion proofs verification {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();

        // now perform the quotiening operation
        let zero_ext = GoldilocksExtAsFieldWrapper::<E, CS>::zero(cs);
        let mut simulated_ext_element = zero_ext;

        assert_eq!(query_index_lsb_first_bits.len(), precomputed_powers.len() - 1);

        let domain_element = pow_from_precomputations(cs, &precomputed_powers[1..], &query_index_lsb_first_bits);

        if idx < 10 {
            println!("    - {} pow from precomputation took {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();

        // we will find it handy to have power of the generator with some bits masked to be zero
        let mut power_chunks = vec![];
        let mut skip_highest_powers = 0;
        // TODO: we may save here (in circuits case especially) if we compute recursively
        for interpolation_degree_log2 in constants.fri_folding_schedule.iter() {
            let domain_element = pow_from_precomputations(
                cs,
                &precomputed_powers_inversed[(1 + interpolation_degree_log2)..],
                &query_index_lsb_first_bits[(skip_highest_powers + interpolation_degree_log2)..],
            );

            skip_highest_powers += *interpolation_degree_log2;
            power_chunks.push(domain_element);
        }
        if idx < 10 {
            println!("    - {} interpolation took {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();

        // don't forget that we are shifted
        let mut domain_element_for_quotiening = domain_element;
        domain_element_for_quotiening.mul_assign(&multiplicative_generator, cs);

        let mut domain_element_for_interpolation = domain_element_for_quotiening;

        if idx < 10 {
            println!("    - {} before quotening took {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();

        verify_quotening_operations(
            cs,
            &mut simulated_ext_element,
            queries,
            proof,
            &public_input_opening_tuples,
            domain_element_for_quotiening,
            challenges,
            verifier,
            constants,
        )?;
        if idx < 10 {
            println!("    - {} verifying quotening took {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();

        let base_coset_inverse = BoojumPrimeField::inverse(&GL::multiplicative_generator()).unwrap();

        let mut current_folded_value: GoldilocksExtAsFieldWrapper<E, CS> = simulated_ext_element;
        let mut subidx = base_tree_idx;
        let mut coset_inverse = base_coset_inverse;

        let mut expected_fri_query_len = base_oracle_depth;

        if idx < 10 {
            println!("    - {} before preparing took {:?} for attempts: {}", idx, now.elapsed(), constants.fri_folding_schedule.len());
        }
        let now = std::time::Instant::now();

        for (idx, (interpolation_degree_log2, fri_query)) in constants.fri_folding_schedule.iter().zip(queries.fri_queries.iter()).enumerate() {
            expected_fri_query_len -= *interpolation_degree_log2;
            let interpolation_degree = 1 << *interpolation_degree_log2;
            let subidx_in_leaf = &subidx[..*interpolation_degree_log2];
            let tree_idx = &subidx[*interpolation_degree_log2..];

            assert_eq!(fri_query.leaf_elements.len(), interpolation_degree * 2);

            let [c0, c1] = current_folded_value.into_coeffs_in_base();

            let c0_from_leaf = binary_select(cs, &fri_query.leaf_elements[..interpolation_degree], subidx_in_leaf)?;
            let c1_from_leaf = binary_select(cs, &fri_query.leaf_elements[interpolation_degree..], subidx_in_leaf)?;

            let c0_is_valid = GoldilocksField::equals(cs, &c0, &c0_from_leaf)?;
            let c1_is_valid = GoldilocksField::equals(cs, &c1, &c1_from_leaf)?;

            validity_flags.push(c0_is_valid);
            validity_flags.push(c1_is_valid);

            // verify query itself
            let cap = if idx == 0 { &proof.fri_base_oracle_cap } else { &proof.fri_intermediate_oracles_caps[idx - 1] };
            assert_eq!(fri_query.proof.len(), expected_fri_query_len);
            validity_flags.push(check_if_included::<E, CS, H>(cs, &fri_query.leaf_elements, &fri_query.proof, &cap, tree_idx)?);

            // interpolate
            let mut elements_to_interpolate = Vec::with_capacity(interpolation_degree);
            for (c0, c1) in fri_query.leaf_elements[..interpolation_degree].iter().zip(fri_query.leaf_elements[interpolation_degree..].iter()) {
                let as_ext = GoldilocksExtAsFieldWrapper::<E, CS>::from_coeffs_in_base([*c0, *c1]);
                elements_to_interpolate.push(as_ext);
            }

            let mut next = Vec::with_capacity(interpolation_degree / 2);
            let challenges = &challenges.fri_intermediate_challenges[idx];
            assert_eq!(challenges.len(), *interpolation_degree_log2);

            let mut base_pow = power_chunks[idx];

            for challenge in challenges.iter() {
                for (i, [a, b]) in elements_to_interpolate.array_chunks::<2>().enumerate() {
                    let mut result = *a;
                    result.add_assign(b, cs);

                    let mut diff = *a;
                    diff.sub_assign(&b, cs);
                    diff.mul_assign(&challenge, cs);
                    // divide by corresponding power
                    let mut pow = base_pow;
                    pow.mul_assign(&interpolation_steps[i], cs);
                    let coset_inverse = GoldilocksAsFieldWrapper::constant(coset_inverse, cs);
                    pow.mul_assign(&coset_inverse, cs);

                    GoldilocksExtAsFieldWrapper::<E, CS>::mul_by_base_and_accumulate_into(&mut result, &pow, &diff, cs)?;

                    // diff.mul_assign_by_base(&pow, cs);
                    // result.add_assign(&diff, cs);

                    next.push(result);
                }

                std::mem::swap(&mut next, &mut elements_to_interpolate);
                next.clear();
                base_pow.square(cs);
                BoojumField::square(&mut coset_inverse);
            }

            for _ in 0..*interpolation_degree_log2 {
                domain_element_for_interpolation.square(cs);
            }

            // recompute the index
            subidx = tree_idx.to_vec();
            current_folded_value = elements_to_interpolate[0];
        }
        if idx < 10 {
            println!("    - {} fri folding took {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();

        // and we should evaluate monomial form and compare

        let mut result_from_monomial = zero_ext;
        // horner rule
        for (c0, c1) in proof.final_fri_monomials[0].iter().zip(proof.final_fri_monomials[1].iter()).rev() {
            let coeff = GoldilocksExtAsFieldWrapper::<E, CS>::from_coeffs_in_base([*c0, *c1]);

            // result_from_monomial = result_from_monomial * z + coeff

            let mut tmp = coeff;
            GoldilocksExtAsFieldWrapper::<E, CS>::mul_by_base_and_accumulate_into(&mut tmp, &domain_element_for_interpolation, &result_from_monomial, cs)?;

            result_from_monomial = tmp;

            // result_from_monomial.mul_assign_by_base(&domain_element_for_interpolation, cs);
            // result_from_monomial.add_assign(&coeff, cs);
        }

        if idx < 10 {
            println!("    - {} horner took {:?}", idx, now.elapsed());
        }
        let now = std::time::Instant::now();

        let result_from_monomial = result_from_monomial.into_coeffs_in_base();
        let current_folded_value = current_folded_value.into_coeffs_in_base();

        let c0_is_valid = GoldilocksField::equals(cs, &result_from_monomial[0], &current_folded_value[0])?;
        let c1_is_valid = GoldilocksField::equals(cs, &result_from_monomial[1], &current_folded_value[1])?;

        validity_flags.push(c0_is_valid);
        validity_flags.push(c1_is_valid);
        if idx < 10 {
            println!("    - {} last part took {:?}", idx, now.elapsed());
        }
    }

    {
        let recv_start = std::time::Instant::now();
        // random order.
        for i in 0..400 {
            let data = rx.recv().unwrap();

            data.mycs.dump_to_existing_cs(cs, data.position);
            //println!("Setting {:?} at {}", data.result, data.request_id);
            validity_flags[data.request_id] = data.result;
        }
        println!("  Receives async took {:?}", recv_start.elapsed());
    }
    println!("  Queries took {:?}", now.elapsed());

    Ok(validity_flags)
}

pub const QUERY_VAR_COUNT: usize = 32003;
pub const STAGE2_VAR_COUNT: usize = 14615;
pub const QUOTIENT_VAR_COUNT: usize = 14618;
pub const SETUP_VAR_COUNT: usize = 21935;

pub const TOTAL_VAR_COUNT: usize = QUERY_VAR_COUNT + STAGE2_VAR_COUNT + QUOTIENT_VAR_COUNT + SETUP_VAR_COUNT;

pub const QUERY_GATES_COUNT: usize = 26403;
pub const STAGE2_GATES_COUNT: usize = 12055;
pub const QUOTIENT_GATES_COUNT: usize = 12058;
pub const SETUP_GATES_COUNT: usize = 18095;

pub const TOTAL_GATES_COUNT: usize = QUERY_GATES_COUNT + STAGE2_GATES_COUNT + QUOTIENT_GATES_COUNT + SETUP_GATES_COUNT;

#[derive(Clone)]
pub struct InclusionResponse<E: Engine> {
    pub result: Boolean,
    pub mycs: MyCS<E>,
    pub position: (usize, usize),
    pub request_id: usize,
}

fn check_if_included_async<E: Engine, CS: ConstraintSystem<E>, H: CircuitGLTreeHasher<E>>(
    request_id: usize,
    start_offset: (usize, usize),
    expected_var_count: usize,
    expected_gates_count: usize,
    leaf_elements: &Vec<GoldilocksField<E>>,
    proof: &Vec<H::CircuitOutput>,
    tree_cap: &Vec<H::CircuitOutput>,
    path: &[Boolean],
) -> InclusionResponse<E> {
    let mut mycs = MyCS::new(start_offset.0);
    use crate::traits::tree_hasher::ToAllocatedNum;
    for elem in leaf_elements {
        let var = elem.into_num().get_variable();
        mycs.add_allocated_num(&var);
    }
    for elem in proof {
        let var = elem.into_allocated_num();
        if let Some(var) = var {
            mycs.add_allocated_num(&var);
        }
    }
    for elem in tree_cap {
        let var = elem.into_allocated_num();
        if let Some(var) = var {
            mycs.add_allocated_num(&var);
        }
    }

    let leaf_hash = <H as CircuitGLTreeHasher<E>>::hash_into_leaf(&mut mycs, leaf_elements).unwrap();
    let result = verify_proof_over_cap::<E, H, _>(&mut mycs, proof, tree_cap, &leaf_hash, path).unwrap();

    let created_vars = mycs.num_variables - mycs.started_variable;
    assert_eq!(created_vars, expected_var_count);
    assert_eq!(mycs.gates.len(), expected_gates_count);

    //mycs.dump_to_existing_cs(cs);
    InclusionResponse {
        result,
        mycs,
        position: start_offset,
        request_id,
    }
}

fn verify_inclusion_proofs<E: Engine, CS: ConstraintSystem<E> + 'static, H: CircuitGLTreeHasher<E>>(
    thread_pool: &ThreadPool,
    cs: &mut CS,
    queries: &AllocatedSingleRoundQueries<E, H>,
    proof: &AllocatedProof<E, H>,
    vk: &AllocatedVerificationKey<E, H>,
    base_tree_idx: &Vec<Boolean>,
    constants: &ConstantsHolder,
    base_oracle_depth: usize,
    validity_flags_index: usize,
    tx: Sender<InclusionResponse<E>>,
) -> Result<(), SynthesisError> {
    //let (tx, rx) = unbounded::<InclusionResponse<E>>();
    let start_offset = cs.reserve(TOTAL_VAR_COUNT, TOTAL_GATES_COUNT)?;

    assert_eq!(constants.witness_leaf_size, queries.witness_query.leaf_elements.len());
    assert_eq!(base_oracle_depth, queries.witness_query.proof.len());

    let mut start_offset = start_offset;
    let tx2 = tx.clone();

    let witness_leaves = queries.witness_query.leaf_elements.clone();
    let witness_proof = queries.witness_query.proof.clone();

    let witness_oracle_cap = proof.witness_oracle_cap.clone();
    let bb = base_tree_idx.clone();

    thread_pool.spawn(move || {
        tx2.send(check_if_included_async::<E, CS, H>(
            validity_flags_index,
            start_offset,
            QUERY_VAR_COUNT,
            QUERY_GATES_COUNT,
            &witness_leaves,
            &witness_proof,
            &witness_oracle_cap,
            &bb,
        ))
        .unwrap();
    });
    start_offset.0 += QUERY_VAR_COUNT;
    start_offset.1 += QUERY_GATES_COUNT;

    assert_eq!(constants.stage_2_leaf_size, queries.stage_2_query.leaf_elements.len());
    assert_eq!(base_oracle_depth, queries.stage_2_query.proof.len());

    let tx2 = tx.clone();

    let stage2_oracle_cap = proof.stage_2_oracle_cap.clone();
    let bb = base_tree_idx.clone();
    let stage2_leaves = queries.stage_2_query.leaf_elements.clone();
    let stage2_query_proof = queries.stage_2_query.proof.clone();

    thread_pool.spawn(move || {
        tx2.send(check_if_included_async::<E, CS, H>(
            validity_flags_index + 1,
            start_offset,
            STAGE2_VAR_COUNT,
            STAGE2_GATES_COUNT,
            &stage2_leaves,
            &stage2_query_proof,
            &stage2_oracle_cap,
            &bb,
        ))
        .unwrap();
    });
    start_offset.0 += STAGE2_VAR_COUNT;
    start_offset.1 += STAGE2_GATES_COUNT;

    assert_eq!(constants.quotient_leaf_size, queries.quotient_query.leaf_elements.len());
    assert_eq!(base_oracle_depth, queries.quotient_query.proof.len());
    let tx2 = tx.clone();

    let quotient_oracle_cap = proof.quotient_oracle_cap.clone();
    let bb = base_tree_idx.clone();
    let quotient_leaves = queries.quotient_query.leaf_elements.clone();
    let quotient_query_proof = queries.quotient_query.proof.clone();

    thread_pool.spawn(move || {
        tx2.send(check_if_included_async::<E, CS, H>(
            validity_flags_index + 2,
            start_offset,
            QUOTIENT_VAR_COUNT,
            QUOTIENT_GATES_COUNT,
            &quotient_leaves,
            &quotient_query_proof,
            &quotient_oracle_cap,
            &bb,
        ))
        .unwrap();
    });
    start_offset.0 += QUOTIENT_VAR_COUNT;
    start_offset.1 += QUOTIENT_GATES_COUNT;

    assert_eq!(constants.setup_leaf_size, queries.setup_query.leaf_elements.len());
    assert_eq!(base_oracle_depth, queries.setup_query.proof.len());
    let tx2 = tx.clone();

    let setup_merkle_tree_cap = vk.setup_merkle_tree_cap.clone();
    let bb = base_tree_idx.clone();
    let setup_leaves = queries.setup_query.leaf_elements.clone();
    let setup_query_proof = queries.setup_query.proof.clone();

    thread_pool.spawn(move || {
        tx2.send(check_if_included_async::<E, CS, H>(
            validity_flags_index + 3,
            start_offset,
            SETUP_VAR_COUNT,
            SETUP_GATES_COUNT,
            &setup_leaves,
            &setup_query_proof,
            &setup_merkle_tree_cap,
            &bb,
        ))
        .unwrap();
    });

    /*
    let mut responses = vec![None; 4];

    let mut validity_flags = vec![Default::default(); 4];

    println!("Validity flag index: {}", validity_flags_index);

    if true {
        //validity_flags_index < 400 {
        for _ in 0..4 {
            let data = rx.recv().unwrap();
            let pos = data.request_id - validity_flags_index;
            responses[pos] = Some(data);
        }

        for i in [0, 1, 3, 2] {
            let data = responses[i].as_ref().unwrap();

            data.mycs.dump_to_existing_cs(cs, data.position);
            println!("Setting {:?} at {}", data.result, data.request_id - validity_flags_index);
            validity_flags[data.request_id - validity_flags_index] = data.result;
        }
    } else {
        // random order.
        for i in 0..4 {
            let data = rx.recv().unwrap();

            data.mycs.dump_to_existing_cs(cs, data.position);
            println!("Setting {:?} at {}", data.result, data.request_id - validity_flags_index);
            validity_flags[data.request_id - validity_flags_index] = data.result;
        }
    }
    */

    /*let validity_flags = validity_flags
    .into_iter()
    .map(|x| x.unwrap())
    .map(|x| {
        x.mycs.dump_to_existing_cs(cs, x.position);
        x.result
    })
    .collect();*/

    Ok(())
}

fn verify_quotening_operations<E: Engine, CS: ConstraintSystem<E> + 'static, H: CircuitGLTreeHasher<E>>(
    cs: &mut CS,
    simulated_ext_element: &mut GoldilocksExtAsFieldWrapper<E, CS>,
    queries: &AllocatedSingleRoundQueries<E, H>,
    proof: &AllocatedProof<E, H>,
    public_input_opening_tuples: &Vec<(GL, Vec<(usize, GoldilocksAsFieldWrapper<E, CS>)>)>,
    domain_element_for_quotiening: GoldilocksAsFieldWrapper<E, CS>,
    challenges: &ChallengesHolder<E, CS>,
    verifier: &WrapperVerifier<E, CS>,
    constants: &ConstantsHolder,
) -> Result<(), SynthesisError> {
    let zero_num = GoldilocksField::zero();
    let zero_base = GoldilocksAsFieldWrapper::<E, CS>::zero(cs);
    let zero_ext = GoldilocksExtAsFieldWrapper::<E, CS>::zero(cs);

    let mut challenge_offset = 0;

    let z_polys_offset = 0;
    let intermediate_polys_offset = 2;
    let lookup_witness_encoding_polys_offset = intermediate_polys_offset + constants.num_intermediate_partial_product_relations * 2;
    let lookup_multiplicities_encoding_polys_offset = lookup_witness_encoding_polys_offset + constants.num_lookup_subarguments * 2;
    let copy_permutation_polys_offset = 0;
    let constants_offset = 0 + constants.num_copy_permutation_polys;
    let lookup_tables_values_offset = 0 + constants.num_copy_permutation_polys + constants.num_constant_polys;
    let variables_offset = 0;
    let witness_columns_offset = constants.num_variable_polys;
    let lookup_multiplicities_offset = witness_columns_offset + constants.num_witness_polys;

    let evaluations = EvaluationsHolder::from_proof(proof);

    {
        let z = challenges.z;
        let z_omega = challenges.z_omega;

        let cast_from_base = move |el: &[GoldilocksField<E>]| {
            el.iter()
                .map(|el| GoldilocksExtAsFieldWrapper::<E, CS>::from_coeffs_in_base([*el, GoldilocksField::zero()]))
                .collect::<Vec<_>>()
        };

        let cast_from_extension = move |el: &[GoldilocksField<E>]| {
            assert_eq!(el.len() % 2, 0);

            el.array_chunks::<2>()
                .map(|[c0, c1]| GoldilocksExtAsFieldWrapper::<E, CS>::from_coeffs_in_base([*c0, *c1]))
                .collect::<Vec<_>>()
        };

        let mut sources = vec![];
        // witness
        sources.extend(cast_from_base(
            &queries.witness_query.leaf_elements[variables_offset..(variables_offset + constants.num_variable_polys)],
        ));
        sources.extend(cast_from_base(
            &queries.witness_query.leaf_elements[witness_columns_offset..(witness_columns_offset + constants.num_witness_polys)],
        ));
        // normal setup
        sources.extend(cast_from_base(&queries.setup_query.leaf_elements[constants_offset..(constants_offset + constants.num_constant_polys)]));
        sources.extend(cast_from_base(
            &queries.setup_query.leaf_elements[copy_permutation_polys_offset..(copy_permutation_polys_offset + constants.num_copy_permutation_polys)],
        ));
        // copy-permutation
        sources.extend(cast_from_extension(&queries.stage_2_query.leaf_elements[z_polys_offset..intermediate_polys_offset]));
        sources.extend(cast_from_extension(
            &queries.stage_2_query.leaf_elements[intermediate_polys_offset..lookup_witness_encoding_polys_offset],
        ));
        // lookup if exists
        sources.extend(cast_from_base(
            &queries.witness_query.leaf_elements[lookup_multiplicities_offset..(lookup_multiplicities_offset + constants.num_multiplicities_polys)],
        ));
        sources.extend(cast_from_extension(
            &queries.stage_2_query.leaf_elements[lookup_witness_encoding_polys_offset..lookup_multiplicities_encoding_polys_offset],
        ));
        sources.extend(cast_from_extension(&queries.stage_2_query.leaf_elements[lookup_multiplicities_encoding_polys_offset..]));
        // lookup setup
        if verifier.lookup_parameters.lookup_is_allowed() {
            let num_lookup_setups = verifier.lookup_parameters.lookup_width() + 1;
            sources.extend(cast_from_base(
                &queries.setup_query.leaf_elements[lookup_tables_values_offset..(lookup_tables_values_offset + num_lookup_setups)],
            ));
        }
        // quotient
        sources.extend(cast_from_extension(&queries.quotient_query.leaf_elements));

        assert_eq!(sources.len(), evaluations.all_values_at_z.len());
        // log!("Making quotiening at Z");
        quotening_operation(
            cs,
            simulated_ext_element,
            &sources,
            &evaluations.all_values_at_z,
            domain_element_for_quotiening,
            z,
            &challenges.challenges_for_fri_quotiening[challenge_offset..(challenge_offset + sources.len())],
        );
        challenge_offset += sources.len();

        // now z*omega
        let mut sources = vec![];
        sources.extend(cast_from_extension(&queries.stage_2_query.leaf_elements[z_polys_offset..intermediate_polys_offset]));

        assert_eq!(sources.len(), evaluations.all_values_at_z_omega.len());
        // log!("Making quotiening at Z*omega");
        quotening_operation(
            cs,
            simulated_ext_element,
            &sources,
            &evaluations.all_values_at_z_omega,
            domain_element_for_quotiening,
            z_omega,
            &challenges.challenges_for_fri_quotiening[challenge_offset..(challenge_offset + sources.len())],
        );

        challenge_offset += sources.len();
        // now at 0 if lookup is needed
        if verifier.lookup_parameters.lookup_is_allowed() {
            let mut sources = vec![];
            // witness encoding
            sources.extend(cast_from_extension(
                &queries.stage_2_query.leaf_elements[lookup_witness_encoding_polys_offset..lookup_multiplicities_encoding_polys_offset],
            ));
            // multiplicities encoding
            sources.extend(cast_from_extension(&queries.stage_2_query.leaf_elements[lookup_multiplicities_encoding_polys_offset..]));

            assert_eq!(sources.len(), evaluations.all_values_at_0.len());
            // log!("Making quotiening at 0 for lookups sumchecks");
            quotening_operation(
                cs,
                simulated_ext_element,
                &sources,
                &evaluations.all_values_at_0,
                domain_element_for_quotiening,
                zero_ext,
                &challenges.challenges_for_fri_quotiening[challenge_offset..(challenge_offset + sources.len())],
            );

            challenge_offset += sources.len();
        }
    }

    // and public inputs
    for (open_at, set) in public_input_opening_tuples.iter() {
        let mut sources = Vec::with_capacity(set.len());
        let mut values = Vec::with_capacity(set.len());
        for (column, expected_value) in set.into_iter() {
            let c0 = queries.witness_query.leaf_elements[*column];
            let el = GoldilocksExtAsFieldWrapper::<E, CS>::from_coeffs_in_base([c0, zero_num]);
            sources.push(el);

            let value = GoldilocksExtAsFieldWrapper::<E, CS>::from_wrapper_coeffs_in_base([*expected_value, zero_base]);
            values.push(value);
        }
        let num_challenges_required = sources.len();
        assert_eq!(values.len(), num_challenges_required);

        // log!("Making quotiening at {} for public inputs", open_at);

        let open_at = GoldilocksAsFieldWrapper::constant(*open_at, cs);
        let open_at = GoldilocksExtAsFieldWrapper::<E, CS>::from_wrapper_coeffs_in_base([open_at, zero_base]);

        quotening_operation(
            cs,
            simulated_ext_element,
            &sources,
            &values,
            domain_element_for_quotiening,
            open_at,
            &challenges.challenges_for_fri_quotiening[challenge_offset..(challenge_offset + sources.len())],
        );

        challenge_offset += num_challenges_required;
    }

    assert_eq!(challenge_offset, challenges.challenges_for_fri_quotiening.len());

    Ok(())
}

fn check_if_included<E: Engine, CS: ConstraintSystem<E>, H: CircuitGLTreeHasher<E>>(
    cs: &mut CS,
    leaf_elements: &Vec<GoldilocksField<E>>,
    proof: &Vec<H::CircuitOutput>,
    tree_cap: &Vec<H::CircuitOutput>,
    path: &[Boolean],
) -> Result<Boolean, SynthesisError> {
    let leaf_hash = <H as CircuitGLTreeHasher<E>>::hash_into_leaf(cs, leaf_elements)?;
    verify_proof_over_cap::<E, H, CS>(cs, proof, tree_cap, &leaf_hash, path)
}

pub fn verify_proof_over_cap<E: Engine, H: CircuitGLTreeHasher<E>, CS: ConstraintSystem<E>>(
    cs: &mut CS,
    proof: &[H::CircuitOutput],
    cap: &[H::CircuitOutput],
    leaf_hash: &H::CircuitOutput,
    path: &[Boolean],
) -> Result<Boolean, SynthesisError> {
    assert!(path.len() >= proof.len());

    let mut current = leaf_hash.clone();
    let path_bits = &path[..proof.len()];
    let cap_bits = &path[proof.len()..];

    for (proof_el, path_bit) in proof.iter().zip(path_bits.iter()) {
        let (left, right) = H::swap_nodes(cs, *path_bit, &current, &proof_el, 0)?;
        current = <H as CircuitGLTreeHasher<E>>::hash_into_node(cs, &left, &right, 0)?;
    }

    let selected_cap_el = H::select_cap_node(cs, cap_bits, cap)?;

    H::compare_output(cs, &current, &selected_cap_el)
    //Ok(Boolean::Constant(true))
}

fn pow_from_precomputations<E: Engine, CS: ConstraintSystem<E> + 'static>(cs: &mut CS, bases: &[GoldilocksAsFieldWrapper<E, CS>], bits: &[Boolean]) -> GoldilocksAsFieldWrapper<E, CS> {
    let mut result = GoldilocksAsFieldWrapper::<E, CS>::one(cs);

    for (base, bit) in bases.iter().zip(bits.iter()) {
        let mut tmp = result;
        tmp.mul_assign(base, cs);
        result = GoldilocksAsFieldWrapper::conditionally_select(cs, *bit, &tmp, &result);
    }

    result
}

fn quotening_operation<E: Engine, CS: ConstraintSystem<E> + 'static>(
    cs: &mut CS,
    dst: &mut GoldilocksExtAsFieldWrapper<E, CS>,
    polynomial_values: &Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    values_at: &Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    domain_element: GoldilocksAsFieldWrapper<E, CS>,
    at: GoldilocksExtAsFieldWrapper<E, CS>,
    challenges: &[GoldilocksExtAsFieldWrapper<E, CS>],
) {
    // we precompute challenges outside to avoid any manual extension ops here
    assert_eq!(polynomial_values.len(), values_at.len());
    assert_eq!(polynomial_values.len(), challenges.len());

    let zero_base = GoldilocksAsFieldWrapper::zero(cs);

    let mut denom = GoldilocksExtAsFieldWrapper::<E, CS>::from_wrapper_coeffs_in_base([domain_element, zero_base]);
    denom.sub_assign(&at, cs);
    denom = denom.inverse(cs);

    let mut acc = GoldilocksExtAsFieldWrapper::<E, CS>::zero(cs);

    for ((poly_value, value_at), challenge) in polynomial_values.iter().zip(values_at.iter()).zip(challenges.iter()) {
        // (f(x) - f(z))/(x - z)
        let mut tmp = *poly_value;
        tmp.sub_assign(&value_at, cs);

        GoldilocksExtAsFieldWrapper::<E, CS>::mul_and_accumulate_into(&mut acc, &tmp, &challenge, cs);

        // let mut as_ext = *challenge;
        // as_ext.mul_assign(&tmp, cs);
        // acc.add_assign(&as_ext, cs);
    }

    GoldilocksExtAsFieldWrapper::<E, CS>::mul_and_accumulate_into(dst, &acc, &denom, cs);

    // acc.mul_assign(&denom, cs);
    // dst.add_assign(&acc, cs);
}
