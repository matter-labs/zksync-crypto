use franklin_crypto::bellman::Engine;
use franklin_crypto::bellman::Field;
use franklin_crypto::bellman::{
    bn256::{Bn256, Fr},
    SynthesisError,
};
use franklin_crypto::plonk::circuit::allocated_num::{AllocatedNum, Num};
use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::{TrivialAssembly, Width4MainGateWithDNext},
    plonk::circuit::Width4WithCustomGates,
};
use rescue_poseidon::rand::{Rand, SeedableRng, XorShiftRng};
use rescue_poseidon::PoseidonParams;
use rescue_poseidon::RescueParams;
use rescue_poseidon::RescuePrimeParams;
#[allow(dead_code)]
use rescue_poseidon::{generic_hash, CircuitGenericSponge, DomainStrategy, GenericSponge};
use std::convert::TryInto;

pub(crate) fn init_rng() -> XorShiftRng {
    const TEST_SEED: [u32; 4] = [0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654];
    XorShiftRng::from_seed(TEST_SEED)
}

pub(crate) fn init_cs<E: Engine>() -> TrivialAssembly<E, Width4WithCustomGates, Width4MainGateWithDNext> {
    TrivialAssembly::<E, Width4WithCustomGates, Width4MainGateWithDNext>::new()
}

fn test_input<E: Engine, const L: usize>() -> ([E::Fr; L], [Num<E>; L]) {
    let rng = &mut init_rng();
    let mut input = [E::Fr::zero(); L];
    let mut input_as_nums = [Num::Constant(E::Fr::zero()); L];
    for (inp, as_num) in input.iter_mut().zip(input_as_nums.iter_mut()) {
        *inp = E::Fr::rand(rng);
        *as_num = Num::Constant(*inp);
    }

    (input, input_as_nums)
}

fn main() {
    run_simple_fixed_len_rescue_hash();
    run_generic_hash_fixed_length::<Bn256>();
    run_generic_hash_var_length::<Bn256>();
    run_circuit_generic_hash_fixed_length::<Bn256>().ok();
    run_circuit_generic_hash_var_length::<Bn256>().ok();
}

fn run_simple_fixed_len_rescue_hash() {
    use rescue_poseidon::rescue_hash;

    const INPUT_LENGTH: usize = 2;
    let input = [Fr::one(); INPUT_LENGTH]; // dummy input

    let result = rescue_hash::<Bn256, INPUT_LENGTH>(&input);
    assert_eq!(result.len(), 2);
}

fn run_generic_hash_fixed_length<E: Engine>() {
    const RATE: usize = 2;
    const WIDTH: usize = 3;
    const INPUT_LENGTH: usize = 2;
    let rng = &mut init_rng();
    let (input, _) = test_input::<Bn256, INPUT_LENGTH>();
    // we can send all type of params so lets start with rescue
    let rescue_params = RescueParams::<Bn256, RATE, WIDTH>::default();
    let result = generic_hash(&rescue_params, &input, None);
    assert_eq!(result.len(), RATE);

    // now, hash with poseidon params
    let poseidon_params = PoseidonParams::<Bn256, RATE, WIDTH>::default();
    let result = generic_hash(&poseidon_params, &input, None);
    assert_eq!(result.len(), RATE);

    // now, hash with rescue prime params
    // in this case we use original domain strategy used in RescuePrime paper
    let rescue_prime_params = RescuePrimeParams::<Bn256, RATE, WIDTH>::default();
    let result = generic_hash(&rescue_prime_params, &input, Some(DomainStrategy::FixedLength));
    assert_eq!(result.len(), RATE);
}

fn run_generic_hash_var_length<E: Engine>() {
    const RATE: usize = 2;
    const WIDTH: usize = 3;
    const INPUT_LENGTH: usize = 4;
    let rng = &mut init_rng();
    let (input, _) = test_input::<Bn256, INPUT_LENGTH>();

    // we can send all type of params so lets start with rescue
    let rescue_params = RescueParams::<Bn256, RATE, WIDTH>::default();
    let mut rescue_hasher = GenericSponge::new();
    rescue_hasher.absorb_multiple(&input, &rescue_params);
    let _ = rescue_hasher.squeeze(&rescue_params).expect("squeezed eleme");

    // go with poseidon
    let poseidon_params = PoseidonParams::<Bn256, RATE, WIDTH>::default();
    let mut poseidon_hasher = GenericSponge::new();
    poseidon_hasher.absorb_multiple(&input, &poseidon_params);
    let _ = poseidon_hasher.squeeze(&poseidon_params).expect("squeezed eleme");
}

fn run_circuit_generic_hash_fixed_length<E: Engine>() -> Result<(), SynthesisError> {
    const RATE: usize = 2;
    const WIDTH: usize = 3;
    const INPUT_LENGTH: usize = 2;
    let cs = &mut init_cs();

    let (_, input) = test_input::<E, INPUT_LENGTH>();
    // we can send all type of params so lets start with rescue
    let rescue_params = RescueParams::<E, RATE, WIDTH>::default();
    let result = CircuitGenericSponge::hash(cs, &input, &rescue_params, None)?;
    assert_eq!(result.len(), RATE);

    // now, hash with poseidon params
    let poseidon_params = PoseidonParams::<E, RATE, WIDTH>::default();
    let result = CircuitGenericSponge::hash(cs, &input, &poseidon_params, None)?;
    assert_eq!(result.len(), RATE);

    // now, hash with rescue prime params
    let rescue_prime_params = RescuePrimeParams::<E, RATE, WIDTH>::default();
    let result = CircuitGenericSponge::hash(cs, &input, &rescue_prime_params, None)?;
    assert_eq!(result.len(), RATE);

    Ok(())
}

fn run_circuit_generic_hash_var_length<E: Engine>() -> Result<(), SynthesisError> {
    const RATE: usize = 2;
    const WIDTH: usize = 3;
    const INPUT_LENGTH: usize = 2;
    let cs = &mut init_cs();

    let (_, input) = test_input::<E, INPUT_LENGTH>();
    // we can send all type of params so lets start with rescue
    let rescue_params = RescueParams::<E, RATE, WIDTH>::default();
    let mut rescue_hasher = CircuitGenericSponge::new();
    rescue_hasher.absorb_multiple(cs, &input, &rescue_params)?;
    let _ = rescue_hasher.squeeze(cs, &rescue_params)?;

    // now, hash with poseidon params
    let poseidon_params = PoseidonParams::<E, RATE, WIDTH>::default();
    let mut poseidon_hasher = CircuitGenericSponge::new();
    poseidon_hasher.absorb_multiple(cs, &input, &rescue_params)?;
    let _ = poseidon_hasher.squeeze(cs, &poseidon_params)?;

    // now, hash with rescue prime params
    let rescue_prime_params = RescuePrimeParams::<E, RATE, WIDTH>::default();
    let mut rescue_prime_hasher = CircuitGenericSponge::new();
    rescue_prime_hasher.absorb_multiple(cs, &input, &rescue_params)?;
    // lets output some num instead of linear-combination
    let _ = rescue_prime_hasher.squeeze_num(cs, &rescue_prime_params)?;

    Ok(())
}
