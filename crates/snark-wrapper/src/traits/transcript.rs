use super::*;

pub trait CircuitGLTranscript<E: Engine>: Clone + Send + Sync + std::fmt::Debug {
    type CircuitCompatibleCap: Clone;
    type TranscriptParameters: Clone + Send + Sync;

    const IS_ALGEBRAIC: bool = true;

    fn new<CS: ConstraintSystem<E>>(cs: &mut CS, params: Self::TranscriptParameters) -> Result<Self, SynthesisError>;

    fn witness_field_elements<CS: ConstraintSystem<E>>(&mut self, cs: &mut CS, field_els: &[GoldilocksField<E>]) -> Result<(), SynthesisError>;

    fn witness_merkle_tree_cap<CS: ConstraintSystem<E>>(&mut self, cs: &mut CS, cap: &Vec<Self::CircuitCompatibleCap>) -> Result<(), SynthesisError>;

    fn get_challenge<CS: ConstraintSystem<E>>(&mut self, cs: &mut CS) -> Result<GoldilocksField<E>, SynthesisError>;

    fn get_multiple_challenges_fixed<CS: ConstraintSystem<E>, const N: usize>(&mut self, cs: &mut CS) -> Result<[GoldilocksField<E>; N], SynthesisError> {
        let mut result = [GoldilocksField::zero(); N];
        for res in result.iter_mut() {
            *res = self.get_challenge(cs)?;
        }

        Ok(result)
    }

    fn get_multiple_challenges<CS: ConstraintSystem<E>>(&mut self, cs: &mut CS, num_challenges: usize) -> Result<Vec<GoldilocksField<E>>, SynthesisError> {
        let mut result = Vec::with_capacity(num_challenges);
        for _ in 0..num_challenges {
            let chal = self.get_challenge(cs)?;
            result.push(chal);
        }

        Ok(result)
    }
}

pub(crate) struct BoolsBuffer {
    pub(crate) available: Vec<Boolean>,
    pub(crate) max_needed: usize,
}

impl BoolsBuffer {
    pub fn get_bits<E: Engine, CS: ConstraintSystem<E>, T: CircuitGLTranscript<E>>(&mut self, cs: &mut CS, transcript: &mut T, num_bits: usize) -> Result<Vec<Boolean>, SynthesisError> {
        if self.available.len() >= num_bits {
            let give: Vec<_> = self.available.drain(..num_bits).collect();

            Ok(give)
        } else {
            let bits_avaiable = GoldilocksField::<E>::ORDER_BITS - self.max_needed;

            // get 1 field element form transcript
            let field_el = transcript.get_challenge(cs)?;
            let el_bits = field_el.spread_into_bits::<CS, 64>(cs)?;
            let mut lsb_iterator = el_bits.iter();

            for _ in 0..bits_avaiable {
                let bit = lsb_iterator.next().unwrap();
                self.available.push(*bit);
            }

            self.get_bits(cs, transcript, num_bits)
        }
    }
}
