use super::*;

use derivative::*;

use franklin_crypto::boojum::algebraic_props::round_function::AbsorptionModeTrait;
use franklin_crypto::boojum::cs::implementations::transcript::Transcript;
use franklin_crypto::boojum::cs::traits::GoodAllocator;
use franklin_crypto::boojum::field::SmallField;
use std::collections::VecDeque;

use franklin_crypto::bellman::{Engine, Field, PrimeField, PrimeFieldRepr};

#[derive(Derivative)]
#[derivative(Clone, Debug)]
pub struct Poseidon2Transcript<E: Engine, F: SmallField, M: AbsorptionModeTrait<E::Fr>, const RATE: usize, const WIDTH: usize> {
    buffer: Vec<E::Fr>,
    last_filled: usize,
    available_challenges: VecDeque<F>,
    #[derivative(Debug = "ignore")]
    sponge: Poseidon2Sponge<E, F, M, RATE, WIDTH>,
}

impl<E: Engine, F: SmallField, M: AbsorptionModeTrait<E::Fr>, const RATE: usize, const WIDTH: usize> Poseidon2Transcript<E, F, M, RATE, WIDTH> {
    pub fn new() -> Self {
        Self {
            buffer: Vec::new(),
            last_filled: 0,
            available_challenges: VecDeque::new(),
            sponge: Poseidon2Sponge::<E, F, M, RATE, WIDTH>::new(),
        }
    }
}

impl<E: Engine, F: SmallField, M: AbsorptionModeTrait<E::Fr>, const RATE: usize, const WIDTH: usize> Transcript<F> for Poseidon2Transcript<E, F, M, RATE, WIDTH> {
    type CompatibleCap = E::Fr;
    type TransciptParameters = ();

    const IS_ALGEBRAIC: bool = true;

    fn new(_params: Self::TransciptParameters) -> Self {
        Self {
            buffer: Vec::new(),
            last_filled: 0,
            available_challenges: VecDeque::new(),
            sponge: Poseidon2Sponge::<E, F, M, RATE, WIDTH>::new(),
        }
    }

    fn witness_field_elements(&mut self, field_els: &[F]) {
        let capasity_per_element = Poseidon2Sponge::<E, F, M, RATE, WIDTH>::capasity_per_element();
        debug_assert!(self.last_filled < capasity_per_element);

        let add_to_last = field_els.len().min((capasity_per_element - self.last_filled) % capasity_per_element);

        if add_to_last != 0 {
            let mut repr_to_add = <E::Fr as PrimeField>::Repr::default();
            for (i, el) in field_els[..add_to_last].iter().enumerate() {
                let mut value_repr = <E::Fr as PrimeField>::Repr::from(el.as_u64_reduced());
                value_repr.shl((i * F::CHAR_BITS) as u32);
                repr_to_add.add_nocarry(&value_repr);
            }
            repr_to_add.shl((self.last_filled * F::CHAR_BITS) as u32);
            self.buffer.last_mut().unwrap().add_assign(&E::Fr::from_repr(repr_to_add).unwrap());
        }

        for chunk in field_els[add_to_last..].chunks(capasity_per_element) {
            let mut repr = <E::Fr as PrimeField>::Repr::default();
            for (i, el) in chunk.iter().enumerate() {
                let mut value_repr = <E::Fr as PrimeField>::Repr::from(el.as_u64_reduced());
                value_repr.shl((i * F::CHAR_BITS) as u32);
                repr.add_nocarry(&value_repr);
            }
            self.buffer.push(E::Fr::from_repr(repr).unwrap());
        }

        self.last_filled = (self.last_filled + field_els.len()) % capasity_per_element;

        self.available_challenges = VecDeque::new();
    }

    fn witness_merkle_tree_cap(&mut self, cap: &[Self::CompatibleCap]) {
        self.last_filled = 0;
        self.buffer.extend_from_slice(cap);

        self.available_challenges = VecDeque::new();
    }

    fn get_challenge(&mut self) -> F {
        assert_eq!(self.sponge.filled, 0);

        if self.buffer.is_empty() {
            if self.available_challenges.len() > 0 {
                return self.available_challenges.pop_front().unwrap();
            } else {
                self.sponge.run_round_function();

                {
                    let commitment = self.sponge.try_get_committment().expect("must have no pending elements in the buffer");
                    for &el in commitment.iter() {
                        self.available_challenges.extend(get_challenges_from_fr::<E, F>(el));
                    }
                }

                return self.get_challenge();
            }
        }

        let to_absorb = std::mem::replace(&mut self.buffer, vec![]);
        self.sponge.absorb(&to_absorb);
        self.last_filled = 0;

        self.available_challenges = VecDeque::new();
        let commitment = self.sponge.finalize();
        for &el in commitment.iter() {
            self.available_challenges.extend(get_challenges_from_fr::<E, F>(el));
        }

        // to avoid duplication
        self.get_challenge()
    }
}

fn get_challenges_from_fr<E: Engine, F: SmallField>(scalar_element: E::Fr) -> Vec<F> {
    assert!(F::CHAR_BITS <= 64, "Goldilocks has less than 64 bits per element");
    let num_challenges = (E::Fr::CAPACITY as usize) / (F::CHAR_BITS as usize);

    scalar_element.into_repr().as_ref()[..num_challenges].iter().map(|x| F::from_u64_with_reduction(*x)).collect()
}
