use crate::franklin_crypto::bellman::pairing::Engine;
use crate::franklin_crypto::bellman::plonk::better_better_cs::cs::ConstraintSystem;
use crate::franklin_crypto::bellman::PrimeFieldRepr;
use crate::franklin_crypto::bellman::{PrimeField, SynthesisError};
use crate::franklin_crypto::plonk::circuit::allocated_num::Num;
use crate::franklin_crypto::plonk::circuit::boolean::Boolean;
use crate::franklin_crypto::plonk::circuit::goldilocks::GoldilocksField;
use crate::franklin_crypto::plonk::circuit::linear_combination::LinearCombination;

use rescue_poseidon::circuit::poseidon2::circuit_poseidon2_round_function;
use rescue_poseidon::poseidon2::Poseidon2Params;

use crate::boojum::field::goldilocks::GoldilocksField as GL;
use crate::boojum::field::PrimeField as BoojumPrimeField;

use derivative::*;

pub mod pow;
pub mod transcript;
pub mod tree_hasher;

#[derive(Derivative)]
#[derivative(Clone, Debug)]
pub struct CircuitPoseidon2Sponge<E: Engine, const RATE: usize, const WIDTH: usize, const CHUNK_BY: usize, const ABSORB_BY_REPLACEMENT: bool> {
    pub(crate) state: [LinearCombination<E>; WIDTH],
    pub(crate) buffer: [LinearCombination<E>; RATE],
    pub(crate) gl_buffer: [GoldilocksField<E>; CHUNK_BY],
    pub(crate) filled: usize,
    #[derivative(Debug = "ignore")]
    pub(crate) params: Poseidon2Params<E, RATE, WIDTH>,
}

impl<E: Engine, const RATE: usize, const WIDTH: usize, const CHUNK_BY: usize, const ABSORB_BY_REPLACEMENT: bool> CircuitPoseidon2Sponge<E, RATE, WIDTH, CHUNK_BY, ABSORB_BY_REPLACEMENT> {
    pub fn new() -> Self {
        Self::new_from_params(Poseidon2Params::default())
    }

    pub fn new_from_params(params: Poseidon2Params<E, RATE, WIDTH>) -> Self {
        assert!(CHUNK_BY == (E::Fr::CAPACITY as usize) / (GL::CHAR_BITS as usize));
        assert!(ABSORB_BY_REPLACEMENT, "Only replacement mode is implemented");

        Self {
            state: [(); WIDTH].map(|_| LinearCombination::zero()),
            buffer: [(); RATE].map(|_| LinearCombination::zero()),
            gl_buffer: [GoldilocksField::zero(); CHUNK_BY],
            filled: 0,
            params,
        }
    }

    pub fn run_round_function<CS: ConstraintSystem<E>>(&mut self, cs: &mut CS) -> Result<(), SynthesisError> {
        circuit_poseidon2_round_function(cs, &self.params, &mut self.state)
    }

    pub fn try_get_commitment<CS: ConstraintSystem<E>>(&self, cs: &mut CS) -> Result<Option<[Num<E>; RATE]>, SynthesisError> {
        if self.filled != 0 {
            return Ok(None);
        }

        let mut result = [Num::zero(); RATE];
        for (dst, src) in result.iter_mut().zip(self.state.iter()) {
            *dst = src.clone().into_num(cs)?;
        }

        Ok(Some(result))
    }

    pub fn absorb_buffer_to_state<CS: ConstraintSystem<E>>(&mut self, cs: &mut CS) -> Result<(), SynthesisError> {
        for (dst, src) in self.state.iter_mut().zip(self.buffer.iter_mut()) {
            *dst = std::mem::replace(src, LinearCombination::zero());
        }

        self.run_round_function(cs)?;
        self.filled = 0;

        Ok(())
    }

    pub fn absorb_single_gl<CS: ConstraintSystem<E>>(&mut self, cs: &mut CS, value: &GoldilocksField<E>) -> Result<(), SynthesisError> {
        debug_assert!(self.filled < RATE * CHUNK_BY);
        let pos = self.filled / CHUNK_BY;
        let exp = self.filled % CHUNK_BY;

        let mut coeff = <E::Fr as PrimeField>::Repr::from(1);
        coeff.shl((exp * GL::CHAR_BITS) as u32);

        self.buffer[pos].add_assign_number_with_coeff(&value.into_num(), E::Fr::from_repr(coeff).unwrap());
        self.filled += 1;

        if self.filled == RATE * CHUNK_BY {
            self.absorb_buffer_to_state(cs)?;
        }

        Ok(())
    }

    pub fn absorb_single<CS: ConstraintSystem<E>>(&mut self, cs: &mut CS, value: Num<E>) -> Result<(), SynthesisError> {
        debug_assert!(self.filled < RATE * CHUNK_BY);
        let pos = self.filled / CHUNK_BY;
        let exp = self.filled % CHUNK_BY;

        match exp {
            0 => {
                self.filled += CHUNK_BY;
                self.buffer[pos] = value.into();
            }
            _ => {
                self.filled = (pos + 1) * CHUNK_BY;

                if self.filled == RATE * CHUNK_BY {
                    self.absorb_buffer_to_state(cs)?;

                    self.buffer[0] = value.into();
                    self.filled = CHUNK_BY;
                } else {
                    self.filled += CHUNK_BY;
                    self.buffer[pos + 1] = value.into();
                }
            }
        }

        if self.filled == RATE * CHUNK_BY {
            self.absorb_buffer_to_state(cs)?;
        }

        Ok(())
    }

    pub fn absorb<T: Into<LinearCombination<E>> + Clone, CS: ConstraintSystem<E>>(&mut self, cs: &mut CS, values: &[T]) -> Result<(), SynthesisError> {
        debug_assert!(self.filled < RATE * CHUNK_BY);
        let mut pos = self.filled / CHUNK_BY;
        let exp = self.filled % CHUNK_BY;
        let len = values.len();

        if exp != 0 {
            pos += 1;
        }

        if len + pos < RATE {
            for (dst, src) in self.buffer[pos..pos + len].iter_mut().zip(values.iter()) {
                *dst = src.clone().into();
            }

            self.filled += len * CHUNK_BY;

            return Ok(());
        }

        let chunks_start = RATE - pos;
        let num_chunks = (len - chunks_start) / RATE;
        let chunk_finish = chunks_start + num_chunks * RATE;

        for (i, value) in values[..chunks_start].iter().enumerate() {
            self.buffer[pos + i] = value.clone().into();
        }
        self.absorb_buffer_to_state(cs)?;

        for chunk in values[chunks_start..chunk_finish].chunks_exact(RATE) {
            for (j, value) in chunk.iter().enumerate() {
                self.state[j] = value.clone().into();
            }
            self.run_round_function(cs)?;
        }

        let new_pos = len - chunk_finish;
        for (dst, src) in self.buffer[..new_pos].iter_mut().zip(values[chunk_finish..].iter()) {
            *dst = src.clone().into();
        }
        self.filled = new_pos * CHUNK_BY;

        Ok(())
    }

    pub fn finalize<CS: ConstraintSystem<E>>(&mut self, cs: &mut CS) -> Result<[Num<E>; RATE], SynthesisError> {
        // padding
        self.absorb_single_gl(cs, &GoldilocksField::one())?;

        if self.filled > 0 {
            self.absorb_buffer_to_state(cs)?;
        }

        let mut result = [Num::zero(); RATE];

        for (dst, src) in result.iter_mut().zip(self.state.iter()) {
            *dst = src.clone().into_num(cs)?;
        }

        Ok(result)
    }

    pub fn finalize_reset<CS: ConstraintSystem<E>>(&mut self, cs: &mut CS) -> Result<[Num<E>; RATE], SynthesisError> {
        // padding
        self.absorb_single_gl(cs, &GoldilocksField::one())?;

        // reset
        let mut state = std::mem::replace(&mut self.state, [(); WIDTH].map(|_| LinearCombination::zero()));

        let filled = self.filled;
        self.filled = 0;

        // run round function if necessary
        if filled > 0 {
            for (dst, src) in state.iter_mut().zip(self.buffer.iter_mut()) {
                *dst = std::mem::replace(src, LinearCombination::zero());
            }

            circuit_poseidon2_round_function(cs, &self.params, &mut state)?;
        }

        let mut result = [Num::zero(); RATE];

        for (dst, src) in result.iter_mut().zip(state.into_iter()) {
            *dst = src.into_num(cs)?;
        }

        Ok(result)
    }
}
