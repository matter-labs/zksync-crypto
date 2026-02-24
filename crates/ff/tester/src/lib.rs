#![feature(llvm_asm)]
#![feature(asm)]

extern crate ff;
extern crate rand as rand_crate;

pub mod assembly_4;
mod check_assembly_4;
mod check_cios;
pub mod mul_variant0;
mod test_large_cios_field;
mod test_large_field;
mod test_short_field;

pub mod adx_4;

pub mod rand {
    pub use crate::rand_crate::{distributions, random, rngs, seq, thread_rng, RngCore, SeedableRng};
    pub use ff::rand::Rng;
    pub use ff::Rand;

    #[derive(Clone, Debug)]
    pub struct XorShiftRng(pub rand_xorshift::XorShiftRng);

    impl XorShiftRng {
        pub fn from_seed(seed: [u32; 4]) -> Self {
            let mut seed_bytes = [0u8; 16];
            for (chunk, word) in seed_bytes.chunks_exact_mut(4).zip(seed.iter()) {
                chunk.copy_from_slice(&word.to_le_bytes());
            }

            <Self as SeedableRng>::from_seed(seed_bytes)
        }
    }

    impl RngCore for XorShiftRng {
        fn next_u32(&mut self) -> u32 {
            self.0.next_u32()
        }

        fn next_u64(&mut self) -> u64 {
            self.0.next_u64()
        }

        fn fill_bytes(&mut self, dest: &mut [u8]) {
            self.0.fill_bytes(dest)
        }

        fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), crate::rand_crate::Error> {
            self.0.try_fill_bytes(dest)
        }
    }

    impl SeedableRng for XorShiftRng {
        type Seed = <rand_xorshift::XorShiftRng as SeedableRng>::Seed;

        fn from_seed(seed: Self::Seed) -> Self {
            Self(rand_xorshift::XorShiftRng::from_seed(seed))
        }
    }
}
