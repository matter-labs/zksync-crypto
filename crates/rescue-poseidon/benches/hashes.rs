use criterion::Criterion;

use franklin_crypto::bellman::{
    pairing::bn256::{Bn256, Fr},
    Field,
};
use franklin_crypto::{
    bellman::plonk::better_better_cs::cs::{TrivialAssembly, Width4MainGateWithDNext},
    bellman::Engine,
    plonk::circuit::Width4WithCustomGates,
};

use rand::{Rand, SeedableRng, XorShiftRng};
use rescue_poseidon::generic_round_function;
use rescue_poseidon::{PoseidonParams, RescueParams, RescuePrimeParams};

fn init_rng() -> XorShiftRng {
    const TEST_SEED: [u32; 4] = [0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654];
    XorShiftRng::from_seed(TEST_SEED)
}
fn init_cs<E: Engine>() -> TrivialAssembly<E, Width4WithCustomGates, Width4MainGateWithDNext> {
    TrivialAssembly::<E, Width4WithCustomGates, Width4MainGateWithDNext>::new()
}

fn test_inputs() -> Vec<Fr> {
    let rng = &mut init_rng();
    let mut inputs = vec![];
    for _ in 0..2 {
        inputs.push(Fr::rand(rng));
    }
    inputs
}
fn test_state_inputs() -> [Fr; 3] {
    let rng = &mut init_rng();
    let mut inputs = [Fr::zero(); 3];
    for i in 0..3 {
        inputs[i] = Fr::rand(rng);
    }
    inputs
}

// fn bench_poseidon_round_function_comparison(crit: &mut Criterion) {
//     let params = PoseidonParams::<Bn256, 2, 3>::default();
//     let mut group = crit.benchmark_group("Poseidon Comparison");

//     group.bench_function("New Poseidon Round Function", |b| {
//         b.iter(|| generic_round_function(&params, &mut test_state_inputs(), None));
//     });
//     use poseidon_hash::bn256::Bn256PoseidonParams;
//     use poseidon_hash::StatefulSponge;

//     let params = Bn256PoseidonParams::new_checked_2_into_1();

//     let mut poseidon = StatefulSponge::<Bn256>::new(&params);
//     group.bench_function("Old Poseidon Round Function", |b| {
//         b.iter(|| poseidon.absorb(&test_inputs()))
//     });
// }

fn bench_rescue_round_function_comparison(crit: &mut Criterion) {
    let params = RescueParams::<Bn256, 2, 3>::default();
    let mut group = crit.benchmark_group("Rescue Round Function Comparison");

    group.bench_function("New Rescue", |b| {
        b.iter(|| generic_round_function(&params, &mut test_state_inputs()));
    });

    use franklin_crypto::rescue::bn256::Bn256RescueParams;
    use franklin_crypto::rescue::StatefulRescue;

    let params = Bn256RescueParams::new_checked_2_into_1();

    let mut rescue = StatefulRescue::<Bn256>::new(&params);
    group.bench_function("Old Rescue", |b| b.iter(|| rescue.absorb(&test_inputs())));
}

fn bench_rescue_round_function(crit: &mut Criterion) {
    let params = RescueParams::<Bn256, 2, 3>::default();
    crit.bench_function("Rescue Round Function", |b| {
        b.iter(|| generic_round_function(&params, &mut test_state_inputs()));
    });
}
fn bench_rescue_round_function_via_addition_chain(crit: &mut Criterion) {
    let params = RescueParams::<Bn256, 2, 3>::specialized_for_num_rounds(8, 120);
    crit.bench_function("Rescue Round Function via addition chain", |b| {
        b.iter(|| generic_round_function(&params, &mut test_state_inputs()));
    });
}
fn bench_poseidon_round_function(crit: &mut Criterion) {
    let params = PoseidonParams::<Bn256, 2, 3>::default();
    crit.bench_function("Poseidon Round Function", |b| {
        b.iter(|| generic_round_function(&params, &mut test_state_inputs()));
    });
}

fn bench_rescue_prime_round_function(crit: &mut Criterion) {
    let params = RescuePrimeParams::<Bn256, 2, 3>::default();
    crit.bench_function("RescuePrime Round Function", |b| {
        b.iter(|| generic_round_function(&params, &mut test_state_inputs()));
    });
}

pub fn group(crit: &mut Criterion) {
    bench_rescue_round_function(crit);
    bench_poseidon_round_function(crit);
    bench_rescue_round_function_comparison(crit);
    bench_rescue_round_function_via_addition_chain(crit);
    // bench_poseidon_round_function_comparison(crit);
    bench_rescue_prime_round_function(crit);
}
