use super::*;

use franklin_crypto::boojum::{worker::Worker, pairing::bls12_381::Fr};
use franklin_crypto::boojum::algebraic_props::round_function::AbsorptionModeTrait;
use franklin_crypto::boojum::field::SmallField;
use franklin_crypto::boojum::cs::implementations::pow::PoWRunner;

use franklin_crypto::bellman::{Engine, Field, PrimeField, PrimeFieldRepr};


const BN256_POSEIDON2_NO_RESULT: u64 = u64::MAX;
const BN256_POSEIDON2_ROUNDS_PER_INVOCAITON: usize = 1 << 16u32;

impl<
    E: Engine,
    F: SmallField,
    M: AbsorptionModeTrait<E::Fr>,
    const RATE: usize,
    const WIDTH: usize,
> PoWRunner for Poseidon2Sponge<E, F, M, RATE, WIDTH> {
    fn run_from_bytes(_seed: Vec<u8>, _pow_bits: u32, _worker: &Worker) -> u64 {
        unimplemented!()
    }

    fn verify_from_bytes(_seed: Vec<u8>, _pow_bits: u32, _challenge: u64) -> bool {
        unimplemented!()
    }

    fn run_from_field_elements<FF: SmallField>(seed: Vec<FF>, pow_bits: u32, worker: &Worker) -> u64 {
        assert!(pow_bits <= 32);

        let mut base_transcript = Self::new();

        // We expect that F == FF == Goldilocks
        if F::CHAR >= FF::CHAR {
            for el in seed.iter() {
                base_transcript.absorb_single_small_field(
                    &F::from_u64(el.as_u64_reduced()).expect("Should be in range")
                );
            }
        } else {
            unimplemented!()
        }

        if pow_bits <= BN256_POSEIDON2_ROUNDS_PER_INVOCAITON.trailing_zeros() {
            // serial case
            log::info!("Do serial PoW");
            for challenge in 0u64..(BN256_POSEIDON2_NO_RESULT - 1) {
                // we expect somewhat "good" hash distribution
                let mut new_transcript = base_transcript.clone();

                let (low, high) = (challenge as u32, (challenge >> 32) as u32);
                let low = F::from_u64_unchecked(low as u64);
                let high = F::from_u64_unchecked(high as u64);

                new_transcript.absorb_single_small_field(&low);
                new_transcript.absorb_single_small_field(&high);

                if new_transcript.finalize()[0].into_repr().as_ref()[0].trailing_zeros() >= pow_bits {
                    return challenge;
                }
            }
        }

        use std::sync::atomic::AtomicU64;
        use std::sync::atomic::Ordering;

        let result = std::sync::Arc::new(AtomicU64::new(BN256_POSEIDON2_NO_RESULT));

        log::info!("Do parallel PoW");

        let pow_rounds_per_invocation = BN256_POSEIDON2_ROUNDS_PER_INVOCAITON as u64;
        // it's good to parallelize
        let num_workers = worker.num_cores as u64;
        worker.scope(0, |scope, _| {
            for worker_idx in 0..num_workers {
                let base_transcript = base_transcript.clone();
                let result = std::sync::Arc::clone(&result);
                scope.spawn(move |_| {
                    for i in
                        0..((BN256_POSEIDON2_NO_RESULT - 1) / num_workers / pow_rounds_per_invocation)
                    {
                        let base = (worker_idx + i * num_workers) * pow_rounds_per_invocation;
                        let current_flag = result.load(Ordering::Relaxed);
                        if current_flag == BN256_POSEIDON2_NO_RESULT {
                            for j in 0..pow_rounds_per_invocation {
                                let challenge_u64 = base + j;

                                let mut new_transcript = base_transcript.clone();

                                let (low, high) = (challenge_u64 as u32, (challenge_u64 >> 32) as u32);
                                let low = F::from_u64_unchecked(low as u64);
                                let high = F::from_u64_unchecked(high as u64);

                                new_transcript.absorb_single_small_field(&low);
                                new_transcript.absorb_single_small_field(&high);

                                if new_transcript.finalize()[0].into_repr().as_ref()[0].trailing_zeros() >= pow_bits {
                                    let _ = result.compare_exchange(
                                        BN256_POSEIDON2_NO_RESULT,
                                        challenge_u64,
                                        Ordering::Acquire,
                                        Ordering::Relaxed,
                                    );

                                    break;
                                }
                            }
                        } else {
                            break;
                        }
                    }
                })
            }
        });

        let challenge_u64 = result.load(Ordering::SeqCst);

        assert!(Self::verify_from_field_elements(seed, pow_bits, challenge_u64));

        challenge_u64
    }
    
    fn verify_from_field_elements<FF: SmallField>(
        seed: Vec<FF>,
        pow_bits: u32,
        challenge: u64,
    ) -> bool {
        assert!(pow_bits <= 32);
        let mut base_transcript = Self::new();

        // We expect that F == FF == Goldilocks
        if F::CHAR >= FF::CHAR {
            for el in seed.iter() {
                base_transcript.absorb_single_small_field(
                    &F::from_u64(el.as_u64_reduced()).expect("Should be in range")
                );
            }
        } else {
            unimplemented!()
        }

        let (low, high) = (challenge as u32, (challenge >> 32) as u32);
        let low = F::from_u64_unchecked(low as u64);
        let high = F::from_u64_unchecked(high as u64);

        base_transcript.absorb_single_small_field(&low);
        base_transcript.absorb_single_small_field(&high);
        
        base_transcript.finalize()[0].into_repr().as_ref()[0].trailing_zeros() >= pow_bits
    }
}
