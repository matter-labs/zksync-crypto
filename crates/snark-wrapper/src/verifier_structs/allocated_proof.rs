use super::*;

use super::allocated_queries::AllocatedSingleRoundQueries;
pub struct AllocatedProof<E: Engine, H: CircuitGLTreeHasher<E>> {
    pub public_inputs: Vec<GoldilocksField<E>>,

    pub witness_oracle_cap: Vec<H::CircuitOutput>,
    pub stage_2_oracle_cap: Vec<H::CircuitOutput>,
    pub quotient_oracle_cap: Vec<H::CircuitOutput>,
    pub final_fri_monomials: [Vec<GoldilocksField<E>>; 2],

    pub values_at_z: Vec<[GoldilocksField<E>; 2]>,
    pub values_at_z_omega: Vec<[GoldilocksField<E>; 2]>,
    pub values_at_0: Vec<[GoldilocksField<E>; 2]>,

    pub fri_base_oracle_cap: Vec<H::CircuitOutput>,
    pub fri_intermediate_oracles_caps: Vec<Vec<H::CircuitOutput>>,

    pub queries_per_fri_repetition: Vec<AllocatedSingleRoundQueries<E, H>>,

    pub pow_challenge_le: [Boolean; 64],
}

impl<E: Engine, HS: TreeHasher<GL, Output = E::Fr>, H: CircuitGLTreeHasher<E, CircuitOutput = Num<E>, NonCircuitSimulator = HS>> AllocatedProof<E, H> {
    pub fn allocate_from_witness<CS: ConstraintSystem<E>>(
        cs: &mut CS,
        witness: &Option<Proof<GL, HS, GLExt2>>,
        verifier: &WrapperVerifier<E, CS>,
        fixed_parameters: &VerificationKeyCircuitGeometry,
        proof_config: &ProofConfig,
    ) -> Result<Self, SynthesisError> {
        if let Some(config) = witness.as_ref().map(|el| &el.proof_config) {
            assert_eq!(config, proof_config);
        }

        let constants = ConstantsHolder::generate(proof_config, verifier, fixed_parameters);

        let num_elements = fixed_parameters.num_public_inputs();
        let public_inputs = witness.as_ref().map(|el| el.public_inputs.iter().copied());
        let public_inputs = allocate_num_elements(cs, num_elements, public_inputs, GoldilocksField::alloc_from_field)?;

        let num_elements = fixed_parameters.cap_size;
        let witness_oracle_cap = witness.as_ref().map(|el| el.witness_oracle_cap.iter().cloned());
        let witness_oracle_cap = allocate_num_elements(cs, num_elements, witness_oracle_cap, Num::alloc)?;

        let num_elements = fixed_parameters.cap_size;
        let stage_2_oracle_cap = witness.as_ref().map(|el| el.stage_2_oracle_cap.iter().cloned());
        let stage_2_oracle_cap = allocate_num_elements(cs, num_elements, stage_2_oracle_cap, Num::alloc)?;

        let num_elements = fixed_parameters.cap_size;
        let quotient_oracle_cap = witness.as_ref().map(|el| el.quotient_oracle_cap.iter().cloned());
        let quotient_oracle_cap = allocate_num_elements(cs, num_elements, quotient_oracle_cap, Num::alloc)?;

        let num_elements = constants.final_expected_degree;
        let final_fri_monomials_c0 = witness.as_ref().map(|el| el.final_fri_monomials[0].iter().cloned());
        let final_fri_monomials_c0 = allocate_num_elements(cs, num_elements, final_fri_monomials_c0, GoldilocksField::alloc_from_field)?;

        let num_elements = constants.final_expected_degree;
        let final_fri_monomials_c1 = witness.as_ref().map(|el| el.final_fri_monomials[1].iter().cloned());
        let final_fri_monomials_c1 = allocate_num_elements(cs, num_elements, final_fri_monomials_c1, GoldilocksField::alloc_from_field)?;

        let num_elements = constants.num_poly_values_at_z;
        let values_at_z = witness.as_ref().map(|el| el.values_at_z.iter().map(|el| el.into_coeffs_in_base()));
        let values_at_z = allocate_num_elements(cs, num_elements, values_at_z, allocate_gl_array)?;

        let num_elements = constants.num_poly_values_at_z_omega;
        let values_at_z_omega = witness.as_ref().map(|el| el.values_at_z_omega.iter().map(|el| el.into_coeffs_in_base()));
        let values_at_z_omega = allocate_num_elements(cs, num_elements, values_at_z_omega, allocate_gl_array)?;

        let num_elements = constants.num_poly_values_at_zero;
        let values_at_0 = witness.as_ref().map(|el| el.values_at_0.iter().map(|el| el.into_coeffs_in_base()));
        let values_at_0 = allocate_num_elements(cs, num_elements, values_at_0, allocate_gl_array)?;

        let num_elements = fixed_parameters.cap_size;
        let fri_base_oracle_cap = witness.as_ref().map(|el| el.fri_base_oracle_cap.iter().cloned());
        let fri_base_oracle_cap = allocate_num_elements(cs, num_elements, fri_base_oracle_cap, Num::alloc)?;

        let fri_folding_schedule = constants.fri_folding_schedule;
        assert!(fri_folding_schedule.len() > 0);
        let mut fri_intermediate_oracles_caps = Vec::with_capacity(fri_folding_schedule.len() - 1);
        for idx in 0..(fri_folding_schedule.len() - 1) {
            let num_elements = fixed_parameters.cap_size;
            let fri_intermediate_cap = witness.as_ref().map(|el| el.fri_intermediate_oracles_caps[idx].iter().cloned());
            let fri_intermediate_cap = allocate_num_elements(cs, num_elements, fri_intermediate_cap, Num::alloc)?;
            fri_intermediate_oracles_caps.push(fri_intermediate_cap);
        }

        let num_items = constants.num_fri_repetitions;
        let mut queries_per_fri_repetition = Vec::with_capacity(num_items);
        for idx in 0..num_items {
            let wit = witness.as_ref().map(|el| el.queries_per_fri_repetition[idx].clone());
            let queries = AllocatedSingleRoundQueries::allocate_from_witness(cs, wit, verifier, fixed_parameters, proof_config)?;
            queries_per_fri_repetition.push(queries);
        }

        let mut pow_challenge_boolean = [Boolean::Constant(true); 64];
        let pow_challenge = witness.as_ref().map(|el| vec![el.pow_challenge]).unwrap_or(vec![]);

        let mut lsb_iter = crate::boojum::utils::LSBIterator::new(&pow_challenge);

        for i in 0..64 {
            pow_challenge_boolean[i] = Boolean::alloc(cs, lsb_iter.next())?;
        }

        let final_fri_monomials = [final_fri_monomials_c0, final_fri_monomials_c1];

        Ok(Self {
            public_inputs,

            witness_oracle_cap,
            stage_2_oracle_cap,
            quotient_oracle_cap,
            final_fri_monomials,

            values_at_z,
            values_at_z_omega,
            values_at_0,

            fri_base_oracle_cap,
            fri_intermediate_oracles_caps,

            queries_per_fri_repetition,

            pow_challenge_le: pow_challenge_boolean,
        })
    }
}

pub fn allocate_gl_array<E: Engine, CS: ConstraintSystem<E>, const N: usize>(cs: &mut CS, source: Option<[GL; N]>) -> Result<[GoldilocksField<E>; N], SynthesisError> {
    let mut result = [GoldilocksField::zero(); N];

    let mut source_it = source.map(|s| s.into_iter());

    for i in 0..N {
        let el = source_it.as_mut().map(|el| el.next().expect("Should be enough elements in the source"));
        result[i] = GoldilocksField::alloc_from_field(cs, el)?;
    }

    Ok(result)
}
