use crate::tests::init_cs;
use franklin_crypto::bellman::pairing::bn256::{Bn256, Fr};
use franklin_crypto::boojum::algebraic_props::round_function::AbsorptionModeTrait;
use franklin_crypto::boojum::cs::implementations::pow::PoWRunner;
use franklin_crypto::boojum::field::goldilocks::GoldilocksField;
use franklin_crypto::boojum::field::SmallField;
use franklin_crypto::boojum::field::U64Representable;
use franklin_crypto::boojum::worker::Worker;
use franklin_crypto::plonk::circuit::{allocated_num::Num, linear_combination::LinearCombination};
use rand::Rand;
use rand::Rng;

use crate::circuit::poseidon2::{circuit_poseidon2_hash, circuit_poseidon2_round_function};
use crate::poseidon::{poseidon_hash, poseidon_round_function};
use crate::poseidon2::{poseidon2_hash, poseidon2_round_function};

use super::Poseidon2Sponge;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct TestingAbsorption;

impl AbsorptionModeTrait<Fr> for TestingAbsorption {
    #[inline(always)]
    fn absorb(dst: &mut Fr, src: &Fr) {
        *dst = *src;
    }
    #[inline(always)]
    fn pad(_dst: &mut Fr) {}
}

#[test]
fn test_different_absorbtions() {
    let num_elements = 1000;
    let mut rng = rand::thread_rng();
    let buffer: Vec<_> = (0..num_elements).map(|_| Fr::rand(&mut rng)).collect();

    // absorb by 1
    let mut hash = Poseidon2Sponge::<Bn256, GoldilocksField, TestingAbsorption, 2, 3>::new();
    for element in buffer.iter() {
        hash.absorb_single(element);
    }
    let commitment1 = hash.finalize();

    // absorb simultaneously
    let mut hash = Poseidon2Sponge::<Bn256, GoldilocksField, TestingAbsorption, 2, 3>::new();
    hash.absorb(&buffer);
    let commitment2 = hash.finalize();

    assert_eq!(commitment1, commitment2);
    dbg!(commitment1);
}

#[ignore]
#[test]
fn test_vs_poseidon() {
    const NUM_ELEMENTS: usize = 10000;
    let mut rng = rand::thread_rng();
    let buffer = [0; NUM_ELEMENTS].map(|_| Fr::rand(&mut rng));

    // hash by poseidon
    let start = std::time::Instant::now();
    poseidon_hash::<Bn256, NUM_ELEMENTS>(&buffer);
    dbg!(start.elapsed());

    // hash by poseidon2
    let start = std::time::Instant::now();
    poseidon2_hash::<Bn256, NUM_ELEMENTS>(&buffer);
    dbg!(start.elapsed());

    let mut buffer = buffer.to_vec();

    // round functions by poseidon
    let params = crate::PoseidonParams::<Bn256, 2, 3>::default();

    let start = std::time::Instant::now();
    for chunk in buffer.chunks_exact_mut(3) {
        poseidon_round_function::<Bn256, _, 2, 3>(&params, chunk.try_into().unwrap());
    }
    dbg!(start.elapsed());

    // round functions by poseidon2
    let params = crate::poseidon2::Poseidon2Params::<Bn256, 2, 3>::default();

    let start = std::time::Instant::now();
    for chunk in buffer.chunks_exact_mut(3) {
        poseidon2_round_function::<Bn256, 2, 3>(chunk.try_into().unwrap(), &params);
    }
    dbg!(start.elapsed());
}

#[test]
fn test_of_sponge_state() {
    let num_elements = 5;
    let mut rng = rand::thread_rng();
    let buffer1: Vec<_> = (0..num_elements).map(|_| Fr::rand(&mut rng)).collect();

    let mut rng = rand::thread_rng();
    let buffer2: Vec<_> = (0..num_elements).map(|_| GoldilocksField::from_u64_unchecked(rng.gen_range(0, GoldilocksField::CHAR))).collect();

    dbg!(&buffer1, &buffer2);

    let mut hash = Poseidon2Sponge::<Bn256, GoldilocksField, TestingAbsorption, 2, 3>::new();

    dbg!(&hash);

    hash.absorb_single(&buffer1[0]);
    dbg!(&hash);

    hash.absorb_single_small_field(&buffer2[0]);
    dbg!(&hash);

    hash.absorb_single_small_field(&buffer2[1]);
    dbg!(&hash);

    hash.absorb_single(&buffer1[1]);
    dbg!(&hash);

    hash.absorb_single(&buffer1[2]);
    dbg!(&hash);

    hash.absorb_single_small_field(&buffer2[2]);
    dbg!(&hash);

    dbg!(&hash.finalize());
    dbg!(&hash);
}

#[test]
fn test_circuit_round_function() {
    let params = crate::poseidon2::Poseidon2Params::<Bn256, 2, 3>::default();

    let cs = &mut init_cs::<Bn256>();

    let mut rng = rand::thread_rng();
    let mut state = [0; 3].map(|_| Fr::rand(&mut rng));
    let mut circuit_state = state.map(|x| Num::alloc(cs, Some(x)).unwrap().into());

    // out of circuit round function
    poseidon2_round_function::<Bn256, 2, 3>(&mut state, &params);

    // circuit round function
    circuit_poseidon2_round_function(cs, &params, &mut circuit_state).unwrap();

    assert_eq!(state, circuit_state.map(|x| x.get_value().unwrap()));
}

#[test]
fn test_circuit_hash() {
    let cs = &mut init_cs::<Bn256>();

    const NUM_ELEMENTS: usize = 10;
    let mut rng = rand::thread_rng();
    let buffer = [0; NUM_ELEMENTS].map(|_| Fr::rand(&mut rng));
    let num_buffer = buffer.map(|x| Num::alloc(cs, Some(x)).unwrap());

    // out of circuit round function
    let hash1 = poseidon2_hash::<Bn256, NUM_ELEMENTS>(&buffer);

    // circuit round function
    let hash2 = circuit_poseidon2_hash(cs, &num_buffer, None).unwrap();

    assert_eq!(hash1, hash2.map(|x| x.get_value().unwrap()));
}

#[test]
fn test_pow_runner() {
    let worker = Worker::new();
    let mut rng = rand::thread_rng();
    let buffer: Vec<_> = (0..4).map(|_| GoldilocksField::from_u64_unchecked(rng.gen_range(0, GoldilocksField::CHAR))).collect();

    let challenge = Poseidon2Sponge::<Bn256, GoldilocksField, TestingAbsorption, 2, 3>::run_from_field_elements(buffer, 10, &worker);

    dbg!(challenge);
}
