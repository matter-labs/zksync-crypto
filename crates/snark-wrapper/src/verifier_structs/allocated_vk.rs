use super::*;

use crate::boojum::cs::implementations::verifier::VerificationKey;

#[derive(Clone, Debug)]
pub struct AllocatedVerificationKey<E: Engine, H: CircuitGLTreeHasher<E>> {
    pub setup_merkle_tree_cap: Vec<H::CircuitOutput>,
}

impl<E: Engine, HS: TreeHasher<GL, Output = E::Fr>, H: CircuitGLTreeHasher<E, CircuitOutput = Num<E>, NonCircuitSimulator = HS>> AllocatedVerificationKey<E, H> {
    pub fn allocate_from_witness<CS: ConstraintSystem<E>>(
        cs: &mut CS,
        witness: Option<VerificationKey<GL, H::NonCircuitSimulator>>,
        vk_fixed_parameters: &VerificationKeyCircuitGeometry,
    ) -> Result<Self, SynthesisError> {
        if let Some(VerificationKey {
            setup_merkle_tree_cap,
            fixed_parameters,
        }) = witness
        {
            assert_eq!(vk_fixed_parameters, &fixed_parameters);

            // allocate fixed length
            assert!(fixed_parameters.cap_size > 0);
            let cap = allocate_num_elements(cs, fixed_parameters.cap_size, Some(setup_merkle_tree_cap.into_iter()), Num::alloc)?;

            Ok(Self { setup_merkle_tree_cap: cap })
        } else {
            let cap = allocate_num_elements(cs, vk_fixed_parameters.cap_size, None::<core::option::IntoIter<E::Fr>>, Num::alloc)?;

            Ok(Self { setup_merkle_tree_cap: cap })
        }
    }

    pub fn allocate_constant(witness: &VerificationKey<GL, H::NonCircuitSimulator>, vk_fixed_parameters: &VerificationKeyCircuitGeometry) -> Self {
        let VerificationKey {
            setup_merkle_tree_cap,
            fixed_parameters,
        } = witness;

        assert_eq!(vk_fixed_parameters, fixed_parameters);

        // allocate fixed length
        assert!(fixed_parameters.cap_size > 0);
        let cap = setup_merkle_tree_cap.iter().map(|x| Num::Constant(*x)).collect();

        Self { setup_merkle_tree_cap: cap }
    }
}
