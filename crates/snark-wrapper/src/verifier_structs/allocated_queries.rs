use super::*;

pub struct AllocatedSingleRoundQueries<E: Engine, H: CircuitGLTreeHasher<E>> {
    // we need query for witness, setup, stage 2 and quotient
    pub witness_query: AllocatedOracleQuery<E, H>,
    pub stage_2_query: AllocatedOracleQuery<E, H>,
    pub quotient_query: AllocatedOracleQuery<E, H>,
    pub setup_query: AllocatedOracleQuery<E, H>,

    pub fri_queries: Vec<AllocatedOracleQuery<E, H>>,
}

impl<E: Engine, HS: TreeHasher<GL, Output = E::Fr>, H: CircuitGLTreeHasher<E, CircuitOutput = Num<E>, NonCircuitSimulator = HS>> AllocatedSingleRoundQueries<E, H> {
    pub fn allocate_from_witness<CS: ConstraintSystem<E>>(
        cs: &mut CS,
        witness: Option<SingleRoundQueries<GL, H::NonCircuitSimulator>>,
        verifier: &WrapperVerifier<E, CS>,
        fixed_parameters: &VerificationKeyCircuitGeometry,
        proof_config: &ProofConfig,
    ) -> Result<Self, SynthesisError> {
        let base_oracle_depth = fixed_parameters.base_oracles_depth();
        let constants = ConstantsHolder::generate(proof_config, verifier, fixed_parameters);

        let witness_leaf_size = constants.witness_leaf_size;
        let witness_query = AllocatedOracleQuery::allocate_from_witness(cs, witness.as_ref().map(|el| el.witness_query.clone()), witness_leaf_size, base_oracle_depth)?;

        let stage_2_leaf_size = constants.stage_2_leaf_size;
        let stage_2_query = AllocatedOracleQuery::allocate_from_witness(cs, witness.as_ref().map(|el| el.stage_2_query.clone()), stage_2_leaf_size, base_oracle_depth)?;

        let quotient_leaf_size = constants.quotient_leaf_size;
        let quotient_query = AllocatedOracleQuery::allocate_from_witness(cs, witness.as_ref().map(|el| el.quotient_query.clone()), quotient_leaf_size, base_oracle_depth)?;

        let setup_leaf_size = constants.setup_leaf_size;
        let setup_query = AllocatedOracleQuery::allocate_from_witness(cs, witness.as_ref().map(|el| el.setup_query.clone()), setup_leaf_size, base_oracle_depth)?;

        // fri is a little bit more involved
        let mut expected_fri_query_len = base_oracle_depth;
        let interpolation_schedule = constants.fri_folding_schedule;
        let mut fri_queries = Vec::with_capacity(interpolation_schedule.len());
        for (idx, interpolation_log_2) in interpolation_schedule.into_iter().enumerate() {
            expected_fri_query_len -= interpolation_log_2;
            let leaf_size = (1 << interpolation_log_2) * 2; // in extension
            let wit = witness.as_ref().map(|el| el.fri_queries[idx].clone());
            let query = AllocatedOracleQuery::allocate_from_witness(cs, wit, leaf_size, expected_fri_query_len)?;
            fri_queries.push(query);
        }

        Ok(Self {
            witness_query,
            stage_2_query,
            quotient_query,
            setup_query,
            fri_queries,
        })
    }
}

pub struct AllocatedOracleQuery<E: Engine, H: CircuitGLTreeHasher<E>> {
    pub leaf_elements: Vec<GoldilocksField<E>>,
    pub proof: Vec<H::CircuitOutput>,
}

impl<E: Engine, HS: TreeHasher<GL, Output = E::Fr>, H: CircuitGLTreeHasher<E, CircuitOutput = Num<E>, NonCircuitSimulator = HS>> AllocatedOracleQuery<E, H> {
    pub fn allocate_from_witness<CS: ConstraintSystem<E>>(cs: &mut CS, witness: Option<OracleQuery<GL, H::NonCircuitSimulator>>, leaf_size: usize, proof_depth: usize) -> Result<Self, SynthesisError> {
        let num_elements = leaf_size;
        let leaf_elements = witness.as_ref().map(|el| el.leaf_elements.iter().copied());
        let leaf_elements = allocate_num_elements(cs, num_elements, leaf_elements, GoldilocksField::alloc_from_field)?;

        let num_elements = proof_depth;
        let proof = witness.as_ref().map(|el| el.proof.iter().cloned());
        let proof = allocate_num_elements(cs, num_elements, proof, Num::alloc)?;

        Ok(Self { leaf_elements, proof })
    }
}
