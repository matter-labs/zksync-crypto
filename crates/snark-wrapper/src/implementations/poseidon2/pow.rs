use crate::traits::pow::RecursivePoWRunner;

use super::*;
use rescue_poseidon::franklin_crypto::bellman::Field;

pub type ConcretePoseidon2SpongeGadget<E> = CircuitPoseidon2Sponge<E, 2, 3, 3, true>;

impl<E: Engine> RecursivePoWRunner<E> for ConcretePoseidon2SpongeGadget<E> {
    fn verify_from_field_elements<CS: ConstraintSystem<E>>(
        cs: &mut CS,
        seed: Vec<GoldilocksField<E>>,
        pow_challenge_le_bits: [Boolean; 64],
        pow_bits: usize,
    ) -> Result<(Boolean, [GoldilocksField<E>; 2]), SynthesisError> {
        let mut sponge = ConcretePoseidon2SpongeGadget::new();

        for el in seed.iter() {
            sponge.absorb_single_gl(cs, el)?;
        }

        // commit nonce
        let mut lc = LinearCombination::zero();
        let mut coeff = E::Fr::one();
        for bit in pow_challenge_le_bits[..32].iter() {
            lc.add_assign_boolean_with_coeff(bit, coeff.clone());
            coeff.double();
        }
        let low = lc.into_num(cs)?;
        let low = GoldilocksField::from_num(cs, low)?;
        sponge.absorb_single_gl(cs, &low)?;

        let mut lc = LinearCombination::zero();
        let mut coeff = E::Fr::one();
        for bit in pow_challenge_le_bits[32..].iter() {
            lc.add_assign_boolean_with_coeff(bit, coeff.clone());
            coeff.double();
        }
        let high = lc.into_num(cs)?;
        let high = GoldilocksField::from_num(cs, high)?;
        sponge.absorb_single_gl(cs, &high)?;

        // get final pow challenge
        let result = sponge.finalize(cs)?[0];

        // verify that  pow challenge has enough zeroes
        let allocated_bools = result.into_bits_le(cs, None)?;

        let mut lc = LinearCombination::zero();
        let coeff = E::Fr::one();
        for b in allocated_bools.iter().take(pow_bits) {
            lc.add_assign_boolean_with_coeff(b, coeff);
        }
        let num_zeroes = lc.into_num(cs)?;
        let result = num_zeroes.is_zero(cs)?;
        Boolean::enforce_equal(cs, &result, &Boolean::constant(true))?;

        Ok((result, [low, high]))
    }
}
