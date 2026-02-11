use super::*;

use crate::traits::transcript::CircuitGLTranscript;

#[derive(Derivative)]
#[derivative(Clone, Debug)]
pub struct CircuitPoseidon2Transcript<E: Engine, const RATE: usize, const WIDTH: usize, const CHUNK_BY: usize, const ABSORB_BY_REPLACEMENT: bool> {
    buffer: Vec<LinearCombination<E>>,
    last_filled: usize,
    available_challenges: Vec<GoldilocksField<E>>,
    #[derivative(Debug = "ignore")]
    sponge: CircuitPoseidon2Sponge<E, RATE, WIDTH, CHUNK_BY, ABSORB_BY_REPLACEMENT>,
}

impl<E: Engine, const RATE: usize, const WIDTH: usize, const CHUNK_BY: usize, const ABSORB_BY_REPLACEMENT: bool> CircuitPoseidon2Transcript<E, RATE, WIDTH, CHUNK_BY, ABSORB_BY_REPLACEMENT> {
    pub fn new() -> Self {
        Self {
            buffer: Vec::new(),
            last_filled: 0,
            available_challenges: Vec::new(),
            sponge: CircuitPoseidon2Sponge::<E, RATE, WIDTH, CHUNK_BY, ABSORB_BY_REPLACEMENT>::new(),
        }
    }
}

impl<E: Engine, const RATE: usize, const WIDTH: usize, const CHUNK_BY: usize, const ABSORB_BY_REPLACEMENT: bool> CircuitGLTranscript<E>
    for CircuitPoseidon2Transcript<E, RATE, WIDTH, CHUNK_BY, ABSORB_BY_REPLACEMENT>
{
    type CircuitCompatibleCap = Num<E>;
    type TranscriptParameters = ();

    const IS_ALGEBRAIC: bool = true;

    fn new<CS: ConstraintSystem<E>>(_cs: &mut CS, _params: Self::TranscriptParameters) -> Result<Self, SynthesisError> {
        Ok(Self::new())
    }

    fn witness_field_elements<CS: ConstraintSystem<E>>(&mut self, _cs: &mut CS, field_els: &[GoldilocksField<E>]) -> Result<(), SynthesisError> {
        debug_assert!(self.last_filled < CHUNK_BY);

        let add_to_last = field_els.len().min((CHUNK_BY - self.last_filled) % CHUNK_BY);

        if add_to_last != 0 {
            for (i, el) in field_els[..add_to_last].iter().enumerate() {
                let mut coeff = <E::Fr as PrimeField>::Repr::from(1);
                coeff.shl(((i + self.last_filled) * GL::CHAR_BITS) as u32);

                self.buffer.last_mut().unwrap().add_assign_number_with_coeff(&el.into_num(), E::Fr::from_repr(coeff).unwrap());
            }
        }

        for chunk in field_els[add_to_last..].chunks(CHUNK_BY) {
            let mut new = LinearCombination::zero();
            let mut coeff = <E::Fr as PrimeField>::Repr::from(1);
            for el in chunk.iter() {
                new.add_assign_number_with_coeff(&el.into_num(), E::Fr::from_repr(coeff).unwrap());
                coeff.shl(GL::CHAR_BITS as u32);
            }
            self.buffer.push(new);
        }

        self.last_filled = (self.last_filled + field_els.len()) % CHUNK_BY;

        Ok(())
    }

    fn witness_merkle_tree_cap<CS: ConstraintSystem<E>>(&mut self, _cs: &mut CS, cap: &Vec<Self::CircuitCompatibleCap>) -> Result<(), SynthesisError> {
        self.last_filled = 0;
        self.buffer.extend(cap.iter().map(|&el| el.into()));

        Ok(())
    }

    fn get_challenge<CS: ConstraintSystem<E>>(&mut self, cs: &mut CS) -> Result<GoldilocksField<E>, SynthesisError> {
        assert_eq!(self.sponge.filled, 0);

        if self.buffer.is_empty() {
            if self.available_challenges.len() > 0 {
                let first_el = self.available_challenges.first().unwrap().clone();
                self.available_challenges.drain(..1);
                return Ok(first_el);
            } else {
                self.sponge.run_round_function(cs)?;

                {
                    let commitment = self.sponge.try_get_commitment(cs)?.expect("must have no pending elements in the buffer");
                    for &el in commitment.iter() {
                        self.available_challenges.extend(get_challenges_from_num(cs, el)?);
                    }
                }

                return self.get_challenge(cs);
            }
        }

        let to_absorb = std::mem::replace(&mut self.buffer, vec![]);
        self.sponge.absorb(cs, &to_absorb)?;
        self.last_filled = 0;

        self.available_challenges = vec![];
        let commitment = self.sponge.finalize(cs)?;
        for &el in commitment.iter() {
            self.available_challenges.extend(get_challenges_from_num(cs, el)?);
        }

        // to avoid duplication
        self.get_challenge(cs)
    }
}

fn get_challenges_from_num<E: Engine, CS: ConstraintSystem<E>>(cs: &mut CS, num: Num<E>) -> Result<Vec<GoldilocksField<E>>, SynthesisError> {
    Ok(GoldilocksField::from_num_to_multiple_with_reduction::<_, 3>(cs, num)?.to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boojum::cs::implementations::transcript::Transcript;
    use crate::boojum::field::{SmallField, U64Representable};
    use rand::{Rand, Rng};

    use crate::franklin_crypto::bellman::pairing::bn256::{Bn256, Fr};
    use crate::franklin_crypto::bellman::plonk::better_better_cs::cs::*;
    use crate::franklin_crypto::plonk::circuit::bigint_new::BITWISE_LOGICAL_OPS_TABLE_NAME;

    use rescue_poseidon::poseidon2::transcript::Poseidon2Transcript;

    use crate::implementations::poseidon2::tree_hasher::AbsorptionModeReplacement;

    #[test]
    fn test_poseidon2_transcript() {
        let mut assembly = TrivialAssembly::<Bn256, PlonkCsWidth4WithNextStepParams, Width4MainGateWithDNext>::new();
        let _before = assembly.n();

        let mut rng = rand::thread_rng();
        let buffer_u64 = [0; 100].map(|_| rng.gen_range(0, GL::CHAR));

        let buffer_circuit = buffer_u64.map(|x| GoldilocksField::alloc_from_u64(&mut assembly, Some(x)).unwrap());

        let buffer_gl = buffer_u64.map(|x| GL::from_u64_unchecked(x));

        // add table for range check
        let columns3 = vec![PolyIdentifier::VariablesPolynomial(0), PolyIdentifier::VariablesPolynomial(1), PolyIdentifier::VariablesPolynomial(2)];

        let name = BITWISE_LOGICAL_OPS_TABLE_NAME;
        let bitwise_logic_table = LookupTableApplication::new(name, TwoKeysOneValueBinopTable::<Bn256, XorBinop>::new(8, name), columns3.clone(), None, true);
        assembly.add_table(bitwise_logic_table).unwrap();

        let mut transcript = Poseidon2Transcript::<Bn256, GL, AbsorptionModeReplacement<Fr>, 2, 3>::new();
        let mut circuit_transcript = CircuitPoseidon2Transcript::<Bn256, 2, 3, 3, true>::new();

        transcript.witness_field_elements(&buffer_gl);
        circuit_transcript.witness_field_elements(&mut assembly, &buffer_circuit).unwrap();

        for _ in 0..5 {
            let chal = transcript.get_challenge();
            let chal_circuit = circuit_transcript.get_challenge(&mut assembly).unwrap();

            assert_eq!(chal, chal_circuit.into_num().get_value().unwrap().into_repr().as_ref()[0]);
        }

        transcript.witness_field_elements(&buffer_gl);
        circuit_transcript.witness_field_elements(&mut assembly, &buffer_circuit).unwrap();

        for _ in 0..10 {
            let chal = transcript.get_challenge();
            let chal_circuit = circuit_transcript.get_challenge(&mut assembly).unwrap();

            assert_eq!(chal, chal_circuit.into_num().get_value().unwrap().into_repr().as_ref()[0]);
        }

        let rand_fr: Vec<_> = (0..10).map(|_| Fr::rand(&mut rng)).collect();
        let num: Vec<_> = rand_fr.iter().map(|x| Num::alloc(&mut assembly, Some(*x)).unwrap()).collect();

        transcript.witness_merkle_tree_cap(&rand_fr);
        circuit_transcript.witness_merkle_tree_cap(&mut assembly, &num).unwrap();

        for _ in 0..5 {
            let chal = transcript.get_challenge();
            let chal_circuit = circuit_transcript.get_challenge(&mut assembly).unwrap();

            assert_eq!(chal, chal_circuit.into_num().get_value().unwrap().into_repr().as_ref()[0]);
        }
    }
}
