extern crate ff;
extern crate rand;
extern crate asm_tester;

use self::ff::*;
use asm_tester::test_large_field::{Fr, FrAsm};

use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn mul_assing_benchmark(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fr = rng.gen();
    let b: Fr = rng.gen();

    c.bench_function("Mont mul assign 256", |bencher| bencher.iter(|| black_box(a).mul_assign(&black_box(b))));
}

fn mul_assing_asm_benchmark(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: FrAsm = rng.gen();
    let b: FrAsm = rng.gen();

    c.bench_function("Mont mul assign 256 ASM", |bencher| bencher.iter(|| black_box(a).mul_assign(&black_box(b))));
}

fn add_assing_benchmark(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fr = rng.gen();
    let b: Fr = rng.gen();

    c.bench_function("Add assign 256", |bencher| bencher.iter(|| black_box(a).add_assign(&black_box(b))));
}

fn add_assing_asm_benchmark(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: FrAsm = rng.gen();
    let b: FrAsm = rng.gen();

    c.bench_function("Add assign 256 ASM", |bencher| bencher.iter(|| black_box(a).add_assign(&black_box(b))));
}

fn sub_assing_benchmark(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fr = rng.gen();
    let b: Fr = rng.gen();

    c.bench_function("Sub assign 256", |bencher| bencher.iter(|| black_box(a).sub_assign(&black_box(b))));
}

fn sub_assing_asm_benchmark(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: FrAsm = rng.gen();
    let b: FrAsm = rng.gen();

    c.bench_function("Sub assign 256 ASM", |bencher| bencher.iter(|| black_box(a).sub_assign(&black_box(b))));
}

fn double_benchmark(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fr = rng.gen();

    c.bench_function("Double 256", |bencher| bencher.iter(|| black_box(a).double()));
}

fn double_asm_benchmark(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: FrAsm = rng.gen();

    c.bench_function("Double 256 ASM", |bencher| bencher.iter(|| black_box(a).double()));
}

fn square_benchmark(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: Fr = rng.gen();

    c.bench_function("Mont square 256", |bencher| bencher.iter(|| black_box(a).square()));
}

fn square_asm_benchmark(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let a: FrAsm = rng.gen();

    c.bench_function("Mont square 256 ASM", |bencher| bencher.iter(|| black_box(a).square()));
}

fn mul_assing_vector_benchmark(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let mut a = [Fr::zero(); 1024];
    let mut b = [Fr::zero(); 1024];
    for (a, b) in a.iter_mut().zip(b.iter_mut()) {
        *a = rng.gen();
        *b = rng.gen();
    }

    c.bench_function("Mont mul assign vector 256", |bencher| bencher.iter(|| 
        {
            let mut a = black_box(a);
            let b = black_box(b);
            for (a, b) in a.iter_mut().zip(b.iter()) {
                a.mul_assign(b);
            }
        }
    ));
}

fn mul_assing_asm_vector_benchmark(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let mut a = [FrAsm::zero(); 1024];
    let mut b = [FrAsm::zero(); 1024];
    for (a, b) in a.iter_mut().zip(b.iter_mut()) {
        *a = rng.gen();
        *b = rng.gen();
    }

    c.bench_function("Mont mul assign vector 256 ASM", |bencher| bencher.iter(|| 
        {
            let mut a = black_box(a);
            let b = black_box(b);
            for (a, b) in a.iter_mut().zip(b.iter()) {
                a.mul_assign(b);
            }
        }
    ));
}

fn square_vector_benchmark(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let mut a = [Fr::zero(); 1024];
    for a in a.iter_mut() {
        *a = rng.gen();
    }

    c.bench_function("Mont square vector 256", |bencher| bencher.iter(|| 
        {
            let mut a = black_box(a);
            for a in a.iter_mut() {
                a.square();
            }
        }
    ));
}

fn square_asm_vector_benchmark(c: &mut Criterion) {
    use rand::{Rng, XorShiftRng, SeedableRng};

    let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);
    let mut a = [FrAsm::zero(); 1024];
    for a in a.iter_mut() {
        *a = rng.gen();
    }

    c.bench_function("Mont square vector 256 ASM", |bencher| bencher.iter(|| 
        {
            let mut a = black_box(a);
            for a in a.iter_mut() {
                a.square();
            }
        }
    ));
}
criterion_group!(
	name = advanced;
    config = Criterion::default().warm_up_time(std::time::Duration::from_secs(5));
    targets = mul_assing_benchmark, mul_assing_asm_benchmark, add_assing_benchmark, add_assing_asm_benchmark, sub_assing_benchmark, sub_assing_asm_benchmark, double_benchmark, double_asm_benchmark, square_benchmark, square_asm_benchmark, mul_assing_vector_benchmark, mul_assing_asm_vector_benchmark, square_vector_benchmark, square_asm_vector_benchmark
);
criterion_main!(advanced);
