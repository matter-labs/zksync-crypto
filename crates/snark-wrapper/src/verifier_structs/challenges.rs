use crate::traits::transcript::CircuitGLTranscript;

use super::allocated_proof::AllocatedProof;
use super::*;
use crate::boojum::field::traits::field_like::PrimeFieldLike;

use crate::franklin_crypto::bellman::plonk::better_better_cs::cs::ConstraintSystem;
use crate::franklin_crypto::plonk::circuit::goldilocks::prime_field_like::GoldilocksExtAsFieldWrapper;
use crate::franklin_crypto::plonk::circuit::goldilocks::prime_field_like::*;

use crate::verifier::utils::materialize_powers_serial;

pub(crate) struct ChallengesHolder<E: Engine, CS: ConstraintSystem<E>> {
    pub(crate) beta: GoldilocksExtAsFieldWrapper<E, CS>,
    pub(crate) gamma: GoldilocksExtAsFieldWrapper<E, CS>,
    pub(crate) lookup_beta: GoldilocksExtAsFieldWrapper<E, CS>,
    pub(crate) lookup_gamma: GoldilocksExtAsFieldWrapper<E, CS>,

    pub(crate) alpha: GoldilocksExtAsFieldWrapper<E, CS>,
    pub(crate) pregenerated_challenges_for_lookup: Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    pub(crate) pregenerated_challenges_for_gates_over_specialized_columns: Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    pub(crate) pregenerated_challenges_for_gates_over_general_purpose_columns: Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    pub(crate) remaining_challenges: Vec<GoldilocksExtAsFieldWrapper<E, CS>>,

    pub(crate) z: GoldilocksExtAsFieldWrapper<E, CS>,
    pub(crate) z_omega: GoldilocksExtAsFieldWrapper<E, CS>,

    pub(crate) challenges_for_fri_quotiening: Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    pub(crate) fri_intermediate_challenges: Vec<Vec<GoldilocksExtAsFieldWrapper<E, CS>>>,
    // pub(crate) challenges: Vec<GoldilocksFieldExt<E>>,
}

impl<E: Engine, CS: ConstraintSystem<E> + 'static> ChallengesHolder<E, CS> {
    pub fn new(cs: &mut CS) -> Self {
        Self {
            beta: GoldilocksExtAsFieldWrapper::zero(cs),
            gamma: GoldilocksExtAsFieldWrapper::zero(cs),
            lookup_beta: GoldilocksExtAsFieldWrapper::zero(cs),
            lookup_gamma: GoldilocksExtAsFieldWrapper::zero(cs),

            alpha: GoldilocksExtAsFieldWrapper::zero(cs),
            pregenerated_challenges_for_lookup: vec![],
            pregenerated_challenges_for_gates_over_specialized_columns: vec![],
            pregenerated_challenges_for_gates_over_general_purpose_columns: vec![],
            remaining_challenges: vec![],

            z: GoldilocksExtAsFieldWrapper::zero(cs),
            z_omega: GoldilocksExtAsFieldWrapper::zero(cs),

            challenges_for_fri_quotiening: vec![],
            fri_intermediate_challenges: vec![],
        }
    }

    pub fn get_beta_gamma_challenges<T: CircuitGLTranscript<E>>(&mut self, cs: &mut CS, transcript: &mut T, verifier: &WrapperVerifier<E, CS>) -> Result<(), SynthesisError> {
        let beta = transcript.get_multiple_challenges_fixed::<_, 2>(cs)?;
        self.beta = GoldilocksExtAsFieldWrapper::from_coeffs_in_base(beta);

        let gamma = transcript.get_multiple_challenges_fixed::<_, 2>(cs)?;
        self.gamma = GoldilocksExtAsFieldWrapper::from_coeffs_in_base(gamma);

        (self.lookup_beta, self.lookup_gamma) = if verifier.lookup_parameters != LookupParameters::NoLookup {
            // lookup argument related parts
            let lookup_beta = transcript.get_multiple_challenges_fixed::<_, 2>(cs)?;
            let lookup_beta = GoldilocksExtAsFieldWrapper::from_coeffs_in_base(lookup_beta);
            let lookup_gamma = transcript.get_multiple_challenges_fixed::<_, 2>(cs)?;
            let lookup_gamma = GoldilocksExtAsFieldWrapper::from_coeffs_in_base(lookup_gamma);

            (lookup_beta, lookup_gamma)
        } else {
            let zero_ext = GoldilocksExtAsFieldWrapper::zero(cs);
            (zero_ext, zero_ext)
        };

        Ok(())
    }

    pub fn get_alpha_powers<T: CircuitGLTranscript<E>>(&mut self, cs: &mut CS, transcript: &mut T, constants: &ConstantsHolder) -> Result<(), SynthesisError> {
        let alpha = transcript.get_multiple_challenges_fixed::<_, 2>(cs)?;
        self.alpha = GoldilocksExtAsFieldWrapper::from_coeffs_in_base(alpha);

        let powers: Vec<_> = materialize_powers_serial(cs, self.alpha, constants.total_num_terms);
        let rest = &powers[..];
        let (take, rest) = rest.split_at(constants.total_num_lookup_argument_terms);
        self.pregenerated_challenges_for_lookup = take.to_vec();
        let (take, rest) = rest.split_at(constants.total_num_gate_terms_for_specialized_columns);
        self.pregenerated_challenges_for_gates_over_specialized_columns = take.to_vec();
        let (take, rest) = rest.split_at(constants.total_num_gate_terms_for_general_purpose_columns);
        self.pregenerated_challenges_for_gates_over_general_purpose_columns = take.to_vec();
        self.remaining_challenges = rest.to_vec();

        Ok(())
    }

    pub fn get_z_challenge<T: CircuitGLTranscript<E>>(&mut self, cs: &mut CS, transcript: &mut T, fixed_parameters: &VerificationKeyCircuitGeometry) -> Result<(), SynthesisError> {
        let z = transcript.get_multiple_challenges_fixed::<_, 2>(cs)?;
        self.z = GoldilocksExtAsFieldWrapper::from_coeffs_in_base(z);

        use crate::boojum::cs::implementations::utils::domain_generator_for_size;
        let omega = domain_generator_for_size::<GL>(fixed_parameters.domain_size as u64);
        let omega_cs_constant = GoldilocksAsFieldWrapper::constant(omega, cs);
        self.z_omega = self.z;
        self.z_omega.mul_assign_by_base(cs, &omega_cs_constant)?;

        Ok(())
    }

    pub fn get_challenges_for_fri_quotiening<T: CircuitGLTranscript<E>>(&mut self, cs: &mut CS, transcript: &mut T, total_num_challenges: usize) -> Result<(), SynthesisError> {
        // get challenges
        let c0 = transcript.get_challenge(cs)?;
        let c1 = transcript.get_challenge(cs)?;

        let challenge = GoldilocksExtAsFieldWrapper::from_coeffs_in_base([c0, c1]);

        self.challenges_for_fri_quotiening = crate::verifier::utils::materialize_powers_serial(cs, challenge, total_num_challenges);

        Ok(())
    }

    pub fn get_fri_intermediate_challenges<H: CircuitGLTreeHasher<E>, TR: CircuitGLTranscript<E, CircuitCompatibleCap = H::CircuitOutput>>(
        &mut self,
        cs: &mut CS,
        transcript: &mut TR,
        proof: &AllocatedProof<E, H>,
        fixed_parameters: &VerificationKeyCircuitGeometry,
        constants: &ConstantsHolder,
    ) -> Result<(), SynthesisError> {
        {
            // now witness base FRI oracle
            assert_eq!(fixed_parameters.cap_size, proof.fri_base_oracle_cap.len());
            transcript.witness_merkle_tree_cap(cs, &proof.fri_base_oracle_cap)?;

            let reduction_degree_log_2 = constants.fri_folding_schedule[0];

            let c0 = transcript.get_challenge(cs)?;
            let c1 = transcript.get_challenge(cs)?;

            let mut challenge_powers = Vec::with_capacity(reduction_degree_log_2);
            let as_extension = GoldilocksExtAsFieldWrapper::from_coeffs_in_base([c0, c1]);
            challenge_powers.push(as_extension);

            let mut current = as_extension;

            for _ in 1..reduction_degree_log_2 {
                current.square(cs);
                challenge_powers.push(current);
            }

            self.fri_intermediate_challenges.push(challenge_powers);
        }

        assert_eq!(constants.fri_folding_schedule[1..].len(), proof.fri_intermediate_oracles_caps.len());

        for (interpolation_degree_log2, cap) in constants.fri_folding_schedule[1..].iter().zip(proof.fri_intermediate_oracles_caps.iter()) {
            // commit new oracle
            assert_eq!(fixed_parameters.cap_size, cap.len());
            transcript.witness_merkle_tree_cap(cs, &cap)?;

            // get challenge
            let reduction_degree_log_2 = *interpolation_degree_log2;
            let c0 = transcript.get_challenge(cs)?;
            let c1 = transcript.get_challenge(cs)?;

            let mut challenge_powers = Vec::with_capacity(reduction_degree_log_2);
            let as_extension = GoldilocksExtAsFieldWrapper::from_coeffs_in_base([c0, c1]);
            challenge_powers.push(as_extension);

            let mut current = as_extension;

            for _ in 1..reduction_degree_log_2 {
                current.square(cs);
                challenge_powers.push(current);
            }

            self.fri_intermediate_challenges.push(challenge_powers);
        }

        Ok(())
    }
}

pub(crate) struct EvaluationsHolder<E: Engine, CS: ConstraintSystem<E>> {
    pub(crate) all_values_at_z: Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    pub(crate) all_values_at_z_omega: Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
    pub(crate) all_values_at_0: Vec<GoldilocksExtAsFieldWrapper<E, CS>>,
}

impl<E: Engine, CS: ConstraintSystem<E>> EvaluationsHolder<E, CS> {
    pub(crate) fn from_proof<H: CircuitGLTreeHasher<E>>(proof: &AllocatedProof<E, H>) -> Self {
        Self {
            all_values_at_z: proof.values_at_z.iter().map(|el| GoldilocksExtAsFieldWrapper::<E, CS>::from_coeffs_in_base(*el)).collect(),
            all_values_at_z_omega: proof.values_at_z_omega.iter().map(|el| GoldilocksExtAsFieldWrapper::<E, CS>::from_coeffs_in_base(*el)).collect(),
            all_values_at_0: proof.values_at_0.iter().map(|el| GoldilocksExtAsFieldWrapper::<E, CS>::from_coeffs_in_base(*el)).collect(),
        }
    }
}
