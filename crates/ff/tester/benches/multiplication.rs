extern crate ff;
extern crate ff_ce_tester;

use self::ff::*;

mod fr {
    use crate::ff::*;

    #[derive(PrimeField)]
    #[PrimeFieldModulus = "21888242871839275222246405745257275088696311157297823662689037894645226208583"]
    #[PrimeFieldGenerator = "2"]
    pub struct Fr(FrRepr);
}

mod frcios{
    use crate::ff::*;

    #[derive(PrimeField)]
    #[PrimeFieldModulus = "21888242871839275222246405745257275088696311157297823662689037894645226208583"]
    #[PrimeFieldGenerator = "2"]
    #[OptimisticCIOSMultiplication = "true"]
    #[OptimisticCIOSSquaring = "true"]
    pub struct FrCios(FrCiosRepr);
}

use self::fr::Fr;
use self::frcios::FrCios;

use criterion::{black_box, criterion_group, criterion_main, Criterion};

// #[inline(always)]
// fn mul<F: PrimeField>(a: F, b: &F) -> F {
//     let mut c = a;
//     c.mul_assign(b);

//     c
// }

// fn multiplication_benchmark(c: &mut Criterion) {
//     use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

//     let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
//     let a: Fr = rng.gen();
//     let b: Fr = rng.gen();

//     c.bench_function("Mont mul 256", |bencher| bencher.iter(|| mul(black_box(a), &black_box(b))));
// }

fn mul_assing_benchmark(c: &mut Criterion) {
    use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fr = rng.gen();
    let b: Fr = rng.gen();

    c.bench_function("Mont mul assign 256", |bencher| bencher.iter(|| black_box(a).mul_assign(&black_box(b))));
}

fn mul_assing_cios_benchmark(c: &mut Criterion) {
    use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: FrCios = rng.gen();
    let b: FrCios = rng.gen();

    c.bench_function("Mont mul assign 256 CIOS derive", |bencher| bencher.iter(|| black_box(a).mul_assign(&black_box(b))));
}

fn square_benchmark(c: &mut Criterion) {
    use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fr = rng.gen();

    c.bench_function("Mont square assign 256", |bencher| bencher.iter(|| black_box(a).square()));
}

fn square_cios_benchmark(c: &mut Criterion) {
    use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: FrCios = rng.gen();

    c.bench_function("Mont square assign 256 CIOS derive", |bencher| bencher.iter(|| black_box(a).square()));
}
// fn mul_assing_custom_benchmark(c: &mut Criterion) {
//     use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

//     use self::ff_ce_tester::mul_variant0::Fs;
//     // use self::mul_variant0::Fs;

//     let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
//     let a: Fs = rng.gen();
//     let b: Fs = rng.gen();

//     c.bench_function("Mont mul assign 256 custom", |bencher| bencher.iter(|| black_box(a).mul_assign(&black_box(b))));
// }

fn mul_assing_rps_benchmark(c: &mut Criterion) {
    use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

    use ff_ce_tester::mul_variant0::Fs;
    // use self::mul_variant0::Fs;

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fs = rng.gen();
    let b: Fs = rng.gen();

    c.bench_function("Mont mul assign 256 RPS", |bencher| bencher.iter(|| black_box(a).rps_mul_assign(&black_box(b))));
}

// fn mul_assing_optimistic_cios_benchmark(c: &mut Criterion) {
//     use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

//     use self::ff_ce_tester::mul_variant0::Fs;
//     // use self::mul_variant0::Fs;

//     let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
//     let a: Fs = rng.gen();
//     let b: Fs = rng.gen();

//     c.bench_function("Mont mul assign 256 optimistic CIOS", |bencher| bencher.iter(|| black_box(a).optimistic_cios_mul_assign(&black_box(b))));
// }

fn mul_assing_optimistic_cios_by_value_benchmark(c: &mut Criterion) {
    use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

    use self::ff_ce_tester::mul_variant0::Fs;
    // use self::mul_variant0::Fs;

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fs = rng.gen();
    let b: Fs = rng.gen();

    c.bench_function("Mont mul assign 256 optimistic CIOS by value", |bencher| bencher.iter(|| black_box(a).optimistic_cios_by_value(black_box(b))));
}

// fn mul_assing_optimistic_cios_by_value_with_partial_red_benchmark(c: &mut Criterion) {
//     use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

//     use self::ff_ce_tester::mul_variant0::Fs;
//     // use self::mul_variant0::Fs;

//     let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
//     let a: Fs = rng.gen();
//     let b: Fs = rng.gen();

//     c.bench_function("Mont mul assign 256 optimistic CIOS by value with_partial_red", |bencher| bencher.iter(|| black_box(a).optimistic_cios_by_value_with_partial_red(black_box(b))));
// }

// fn mulx_mul_assing_benchmark(c: &mut Criterion) {
//     use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

//     use self::ff_ce_tester::mul_variant0::Fs;
//     // use self::mul_variant0::Fs;

//     let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
//     let a: Fs = rng.gen();
//     let b: Fs = rng.gen();

//     c.bench_function("Mont mul assign 256 with MULX latency", |bencher| bencher.iter(|| black_box(a).mulx_latency_mul_assign(&black_box(b))));
// }

fn llvm_asm_mul_assing_benchmark(c: &mut Criterion) {
    use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

    use self::ff_ce_tester::mul_variant0::{Fs, mont_mul_asm};
    // use self::mul_variant0::Fs;

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fs = rng.gen();
    let b: Fs = rng.gen();

    let a = a.into_raw_repr().0;
    let b = b.into_raw_repr().0;

    c.bench_function("Mont mul assign 256 with LLVM assembly", |bencher| bencher.iter(|| mont_mul_asm(
            black_box(&a),
            black_box(&b)
        )
    ));
}

fn new_asm_mul_assing_benchmark(c: &mut Criterion) {
    use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

    use self::ff_ce_tester::assembly_4::*;
    use ff_ce_tester::mul_variant0::Fs;

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fs = rng.gen();
    let b: Fs = rng.gen();

    let a = a.into_raw_repr().0;
    let b = b.into_raw_repr().0;

    c.bench_function("Mont mul assign 256 with new assembly", |bencher| bencher.iter(|| mont_mul_asm(
            black_box(&a),
            black_box(&b)
        )
    ));
}

fn adx_asm_mul_assing_benchmark(c: &mut Criterion) {
    use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

    use self::ff_ce_tester::assembly_4::*;
    use ff_ce_tester::mul_variant0::Fs;

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fs = rng.gen();
    let b: Fs = rng.gen();

    let a = a.into_raw_repr().0;
    let b = b.into_raw_repr().0;

    c.bench_function("Mont mul assign 256 with new assembly with ADX", |bencher| bencher.iter(|| mont_mul_asm_adx(
            black_box(&a),
            black_box(&b)
        )
    ));
}

fn asm_mul_assing_with_register_abi_benchmark(c: &mut Criterion) {
    use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

    use self::ff_ce_tester::assembly_4::*;
    use ff_ce_tester::mul_variant0::Fs;

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fs = rng.gen();
    let b: Fs = rng.gen();

    let [a0, a1, a2, a3] = a.into_raw_repr().0;
    let b = b.into_raw_repr().0;

    c.bench_function("Mont mul assign 256 with new assembly and ABI through registers", |bencher| bencher.iter(|| mont_mul_asm_through_registers(
            black_box(a0),
            black_box(a1),
            black_box(a2),
            black_box(a3),
            black_box(&b)
        )
    ));
}

fn proth_adx_asm_mul_assing_benchmark(c: &mut Criterion) {
    use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

    use self::ff_ce_tester::assembly_4::*;
    use ff_ce_tester::mul_variant0::Fs;

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fs = rng.gen();
    let b: Fs = rng.gen();

    let a = a.into_raw_repr().0;
    let b = b.into_raw_repr().0;

    c.bench_function("Mont mul assign for Proth prime with new assembly with ADX", |bencher| bencher.iter(|| mont_mul_asm_adx_for_proth_prime(
            black_box(&a),
            black_box(&b)
        )
    ));
}

// fn mul_assing_optimistic_cios_with_different_semantics_benchmark(c: &mut Criterion) {
//     use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

//     use self::ff_ce_tester::mul_variant0::Fs;
//     // use self::mul_variant0::Fs;

//     let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
//     let a: Fs = rng.gen();
//     let b: Fs = rng.gen();

//     c.bench_function("Mont mul assign 256 optimistic CIOS different semantics", |bencher| bencher.iter(|| black_box(a).optimistic_cios_mul_assign_with_different_semantics(&black_box(b))));
// }

// fn mul_assing_vector_benchmark(c: &mut Criterion) {
//     use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

//     let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
//     let mut a = [Fr::zero(); 1024];
//     let mut b = [Fr::zero(); 1024];
//     for (a, b) in a.iter_mut().zip(b.iter_mut()) {
//         *a = rng.gen();
//         *b = rng.gen();
//     }

//     c.bench_function("Mont mul assign vector 256", |bencher| bencher.iter(|| 
//         {
//             let mut a = black_box(a);
//             let b = black_box(b);
//             for (a, b) in a.iter_mut().zip(b.iter()) {
//                 a.mul_assign(b);
//             }
//         }
//     ));
// }

// fn mul_assing_custom_vector_benchmark(c: &mut Criterion) {
//     use ff_ce_tester::rand::{Rng, XorShiftRng, SeedableRng};

//     use self::ff_ce_tester::mul_variant0::Fs;

//     let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
//     let mut a = [Fs::zero(); 1024];
//     let mut b = [Fs::zero(); 1024];
//     for (a, b) in a.iter_mut().zip(b.iter_mut()) {
//         *a = rng.gen();
//         *b = rng.gen();
//     }

//     c.bench_function("Mont mul assign custom vector 256", |bencher| bencher.iter(|| 
//         {
//             let mut a = black_box(a);
//             let b = black_box(b);
//             for (a, b) in a.iter_mut().zip(b.iter()) {
//                 a.mul_assign(b);
//             }
//         }  
//     ));
// }

// criterion_group!(benches, mul_assing_benchmark, mul_assing_rps_benchmark, mul_assing_optimistic_cios_benchmark);
// criterion_group!(benches, mul_assing_vector_benchmark, mul_assing_custom_vector_benchmark);
// criterion_main!(benches);

// criterion_group!(
// 	name = advanced;
//     config = Criterion::default().warm_up_time(std::time::Duration::from_secs(10));
//     targets = mul_assing_benchmark, mul_assing_rps_benchmark, mul_assing_optimistic_cios_benchmark, mul_assing_optimistic_cios_with_different_semantics_benchmark
// );

// criterion_group!(
// 	name = advanced;
//     config = Criterion::default().warm_up_time(std::time::Duration::from_secs(10));
//     targets = mul_assing_optimistic_cios_benchmark, mul_assing_optimistic_cios_by_value_benchmark, asm_mul_assing_benchmark
// );
// criterion_group!(
// 	name = advanced;
//     config = Criterion::default().warm_up_time(std::time::Duration::from_secs(5));
//     targets = mul_assing_benchmark, mul_assing_cios_benchmark, mul_assing_optimistic_cios_by_value_benchmark, llvm_asm_mul_assing_benchmark, new_asm_mul_assing_benchmark, adx_asm_mul_assing_benchmark, asm_mul_assing_with_register_abi_benchmark, proth_adx_asm_mul_assing_benchmark
// );
criterion_group!(
	name = advanced;
    config = Criterion::default().warm_up_time(std::time::Duration::from_secs(5));
    targets = llvm_asm_mul_assing_benchmark, new_asm_mul_assing_benchmark, adx_asm_mul_assing_benchmark, proth_adx_asm_mul_assing_benchmark
);
criterion_main!(advanced);
