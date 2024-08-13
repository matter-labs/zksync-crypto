use super::*;

use derivative::*;
use franklin_crypto::boojum::field::SmallField;
use franklin_crypto::boojum::cs::oracle::TreeHasher;
use franklin_crypto::bellman::{Engine, Field, PrimeField, PrimeFieldRepr};
use franklin_crypto::boojum::algebraic_props::round_function::AbsorptionModeTrait;

use typemap_rev::{TypeMap, TypeMapKey};
use std::sync::{Arc, RwLock};

impl<E: Engine, const RATE: usize, const WIDTH: usize> TypeMapKey for Poseidon2Params::<E, RATE, WIDTH> {
    type Value = Arc<Poseidon2Params::<E, RATE, WIDTH>>;
}

#[derive(Derivative)]
#[derivative(Clone, Debug)]
pub struct Poseidon2Sponge<
    E: Engine,
    F: SmallField,
    M: AbsorptionModeTrait<E::Fr>,
    const RATE: usize,
    const WIDTH: usize
>{
    pub(crate) state: [E::Fr; WIDTH],
    pub(crate) buffer: [E::Fr; RATE],
    pub(crate) filled: usize,
    #[derivative(Debug = "ignore")]
    pub(crate) params: Arc<Poseidon2Params<E, RATE, WIDTH>>,
    _marker: std::marker::PhantomData<(F, M)>,
}

impl<
    E: Engine,
    F: SmallField,
    M: AbsorptionModeTrait<E::Fr>,
    const RATE: usize,
    const WIDTH: usize,
> Poseidon2Sponge<E, F, M, RATE, WIDTH> {
    pub fn new() -> Self {
        assert!(Self::capasity_per_element() > 0);

        lazy_static::lazy_static!{
            static ref POSEIDON_PARAMS: RwLock<TypeMap> = RwLock::new(TypeMap::new());
        };

        let static_params = POSEIDON_PARAMS.read().unwrap();
        let params = static_params.get::<Poseidon2Params<E, RATE, WIDTH>>().map(|p| p.clone());
        drop(static_params);

        let params = if let Some(params) = params {
            params
        } else {
            let params = Arc::new(Poseidon2Params::<E, RATE, WIDTH>::default());
            let mut static_params = POSEIDON_PARAMS.write().unwrap();
            static_params.insert::<Poseidon2Params<E, RATE, WIDTH>>(params.clone());
            params
        };

        Self {
            params,
            state: [E::Fr::zero(); WIDTH],
            buffer: [E::Fr::zero(); RATE],
            filled: 0,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn capasity_per_element() -> usize {
        (E::Fr::CAPACITY as usize) / (F::CHAR_BITS as usize)
    }

    pub fn run_round_function(&mut self) {
        poseidon2_round_function(&mut self.state, &self.params);
    }

    pub fn try_get_committment(&mut self) -> Option<[E::Fr; RATE]> {
        if self.filled != 0 {
            return None;
        }

        Some(self.state[..RATE].try_into().unwrap())
    }

    pub fn absorb_buffer_to_state(&mut self) {
        for (dst, src) in self.state.iter_mut()
            .zip(self.buffer.iter_mut())
        {
            M::absorb(dst, src);
            *src = E::Fr::zero();
        }

        self.run_round_function();
        self.filled = 0;
    }

    pub fn absorb_single_small_field(&mut self, value: &F) {
        let capasity_per_element = Self::capasity_per_element();
        debug_assert!(self.filled < RATE * capasity_per_element);
        let pos = self.filled / capasity_per_element;
        let exp = self.filled % capasity_per_element;

        let mut value_repr = <E::Fr as PrimeField>::Repr::from(value.as_u64_reduced());
        value_repr.shl((exp * F::CHAR_BITS) as u32);

        self.buffer[pos].add_assign(&E::Fr::from_repr(value_repr).unwrap());
        self.filled += 1;

        if self.filled == RATE * capasity_per_element {
            self.absorb_buffer_to_state();
        }
    }

    pub fn absorb_single(&mut self, value: &E::Fr) {
        let capasity_per_element = Self::capasity_per_element();
        debug_assert!(self.filled < RATE * capasity_per_element);
        let pos = self.filled / capasity_per_element;
        let exp = self.filled % capasity_per_element;

        match exp {
            0 => {
                self.filled += capasity_per_element;
                self.buffer[pos] = *value;
            },
            _ => {
                self.filled = (pos + 1) * capasity_per_element;

                if self.filled == RATE * capasity_per_element {
                    self.absorb_buffer_to_state();

                    self.buffer[0] = *value;
                    self.filled = capasity_per_element;
                } else {
                    self.filled += capasity_per_element;
                    self.buffer[pos + 1] = *value;
                }
            }
        }

        if self.filled == RATE * capasity_per_element {
            self.absorb_buffer_to_state();
        }
    }

    pub fn absorb(&mut self, values: &[E::Fr]) {
        let capasity_per_element = Self::capasity_per_element();
        debug_assert!(self.filled < RATE * capasity_per_element);
        let mut pos = self.filled / capasity_per_element;
        let exp = self.filled % capasity_per_element;
        let len = values.len();

        if exp != 0 {
            pos += 1;
        }

        if len + pos < RATE {
            self.buffer[pos..pos+len].copy_from_slice(values);

            self.filled += len * capasity_per_element;

            return;
        }

        let chunks_start = RATE - pos;
        let num_chunks = (len - chunks_start) / RATE;
        let chunk_finish = chunks_start + num_chunks * RATE;

        for (i, value) in values[..chunks_start].iter().enumerate() {
            self.buffer[pos + i] = *value;
        }
        self.absorb_buffer_to_state();

        for chunk in values[chunks_start..chunk_finish].chunks_exact(RATE) {
            for (j, value) in chunk.iter().enumerate() {
                M::absorb(&mut self.state[j], value);
            }
            self.run_round_function();
        }

        let new_pos = len - chunk_finish;
        self.buffer[..new_pos].copy_from_slice(&values[chunk_finish..]);
        self.filled = new_pos * capasity_per_element;
    }

    pub fn finalize(&mut self) -> [E::Fr; RATE] {
        // padding
        self.absorb_single_small_field(&F::ONE);

        if self.filled > 0 {
            self.absorb_buffer_to_state();
        }

        self.state[..RATE].try_into().unwrap()
    }

    pub fn finalize_reset(&mut self) -> [E::Fr; RATE] {
        // padding
        self.absorb_single_small_field(&F::ONE);

        // reset
        let mut state = std::mem::replace(&mut self.state, [E::Fr::zero(); WIDTH]);
        let filled = self.filled;
        self.filled = 0;

        // run round function if necessary
        if filled > 0 {
            for (dst, src) in state.iter_mut().zip(self.buffer.iter_mut()) {
                M::absorb(dst, src);
                *src = E::Fr::zero();
            }

            poseidon2_round_function(&mut state, &self.params);
        }

        self.state[..RATE].try_into().unwrap()
    }
}

impl<
    E: Engine,
    F: SmallField,
    M: AbsorptionModeTrait<E::Fr>,
    const RATE: usize,
    const WIDTH: usize,
> TreeHasher<F> for Poseidon2Sponge<E, F, M, RATE, WIDTH> {
    type Output = E::Fr;

    #[inline]
    fn new() -> Self {
        Self::new()
    }

    #[inline]
    fn placeholder_output() -> Self::Output {
        E::Fr::zero()
    }

    #[inline]
    fn accumulate_into_leaf(&mut self, value: &F) {
        self.absorb_single_small_field(value);
    }

    #[inline]
    fn finalize_into_leaf_hash_and_reset(&mut self) -> Self::Output {
        self.finalize_reset()[0]
    }

    #[inline]
    fn hash_into_leaf<'a, S: IntoIterator<Item = &'a F>>(source: S) -> Self::Output
    where
        F: 'a 
    {
        let mut hasher = Self::new();

        for el in source.into_iter() {
            hasher.absorb_single_small_field(el);
        }
        hasher.finalize()[0]
    }

    #[inline]
    fn hash_into_leaf_owned<S: IntoIterator<Item = F>>(source: S) -> Self::Output {
        let mut hasher = Self::new();

        for el in source.into_iter() {
            hasher.absorb_single_small_field(&el);
        }
        hasher.finalize()[0]
    }

    #[inline]
    fn hash_into_node(left: &Self::Output, right: &Self::Output, _depth: usize) -> Self::Output {
        lazy_static::lazy_static!{
            static ref POSEIDON_PARAMS: RwLock<TypeMap> = RwLock::new(TypeMap::new());
        };

        let static_params = POSEIDON_PARAMS.read().unwrap();
        let params = static_params.get::<Poseidon2Params<E, RATE, WIDTH>>().map(|p| p.clone());
        drop(static_params);

        let params = if let Some(params) = params {
            params
        } else {
            let params = Arc::new(Poseidon2Params::<E, RATE, WIDTH>::default());
            let mut static_params = POSEIDON_PARAMS.write().unwrap();
            static_params.insert::<Poseidon2Params<E, RATE, WIDTH>>(params.clone());
            params
        };

        let mut state = [E::Fr::zero(); WIDTH];
        M::absorb(&mut state[0], left);
        M::absorb(&mut state[1], right);

        poseidon2_round_function(&mut state, params.as_ref());

        state[0]
    }
}
