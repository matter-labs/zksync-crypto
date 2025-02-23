#![allow(internal_features)]
// for generic const expr
#![allow(incomplete_features)]
#![cfg_attr(not(test), no_std)]
#![feature(generic_const_exprs)]
#![feature(associated_type_defaults)]
#![feature(core_intrinsics)]
#![feature(const_eval_select)]
#![cfg_attr(target_feature = "avx512f", feature(stdarch_x86_avx512))]
#![feature(const_mut_refs)]
#![feature(const_option)]

use core::fmt::Debug;
use core::fmt::Display;
use core::fmt::Formatter;
/// The prime field `F_p` where `p = 2^31 - 1`.
use core::hash::Hash;
use core::hash::Hasher;
use rand::Rng;

pub(crate) mod ops;

pub mod base;
pub mod complex;
pub mod field;
pub mod field_like;
pub mod quartic;

pub use self::base::*;
pub use self::complex::*;
pub use self::quartic::*;

pub use self::field::*;
pub use self::field_like::*;

#[cfg(target_feature = "avx512f")]
pub mod avx_512_impl;
#[cfg(target_feature = "avx512f")]
pub mod ext2_avx_512_impl;
#[cfg(target_feature = "avx512f")]
pub mod ext_avx_512_impl;
#[cfg(target_feature = "avx512f")]
pub mod ext_avx_512_interleaved_impl;

#[cfg(target_feature = "avx512f")]
pub use self::avx_512_impl::*;
#[cfg(target_feature = "avx512f")]
pub use self::ext2_avx_512_impl::*;
#[cfg(target_feature = "avx512f")]
pub use self::ext_avx_512_impl::*;
#[cfg(target_feature = "avx512f")]
pub use self::ext_avx_512_interleaved_impl::*;

#[cfg(not(target_feature = "avx512f"))]
pub mod arm_impl;
#[cfg(not(target_feature = "avx512f"))]
pub mod ext_arm_impl;
#[cfg(not(target_feature = "avx512f"))]
pub mod ext_arm_interleaved_impl;

#[cfg(not(target_feature = "avx512f"))]
pub use self::arm_impl::*;
#[cfg(not(target_feature = "avx512f"))]
pub use self::ext_arm_impl::*;
#[cfg(not(target_feature = "avx512f"))]
pub use self::ext_arm_interleaved_impl::*;

const _: () = const {
    #[cfg(all(any(feature = "use_division", feature = "modular_ops"), not(target_arch = "riscv32")))]
    compile_error!("`use_division` and `modular ops` features are intended for simulated (provable) machines and should not be activated otherwise");

    ()
};

pub fn batch_inverse_checked<F: Field>(input: &mut [F], tmp_buffer: &mut [F]) -> bool {
    assert_eq!(input.len(), tmp_buffer.len());
    if input.is_empty() {
        return true;
    }

    // We do Montgomery batch inversion trick, and reuse a buffer, but record if we encountered any zero in the meantime.
    // We also skip such zero element
    tmp_buffer[0] = F::ONE;
    let mut zero_encountered = false;
    let mut accumulator: F;
    if input[0].is_zero() {
        zero_encountered = true;
        accumulator = F::ONE;
    } else {
        accumulator = input[0];
    }
    for (el, out) in input.iter().zip(tmp_buffer.iter_mut()).skip(1) {
        *out = accumulator;
        if el.is_zero() {
            zero_encountered = true;
        } else {
            accumulator.mul_assign(el);
        }
    }

    // for a set of a, b, c, d we have
    // - input = [1, a, ab, abc],
    // - accumulator = abcd
    let mut grand_inverse = accumulator.inverse().expect("batch inverse must be called on sets without zeroes");

    // grand_inverse = a^-1 b^-1 c^-1 d^-1
    for (tmp, original) in tmp_buffer.iter().rev().zip(input.iter_mut().rev()) {
        let mut tmp = *tmp; // abc
        tmp.mul_assign(&grand_inverse); // d^-1
        if original.is_zero() == false {
            grand_inverse.mul_assign(original); // e.g. it's now a^-1 b^-1 c^-1
            *original = tmp;
        }
    }

    !zero_encountered
}

pub trait Rand {
    fn random_element<R: Rng + ?Sized>(rng: &mut R) -> Self;
}

pub fn rand_fp_from_rng<R: rand::Rng>(rng: &mut R) -> Mersenne31Field {
    Mersenne31Field::from_u64_unchecked(rng.gen_range(0..((1 << 31) - 1)))
}

pub fn rand_fp2_from_rng<R: rand::Rng>(rng: &mut R) -> Mersenne31Complex {
    let a = Mersenne31Field::from_u64_unchecked(rng.gen_range(0..((1 << 31) - 1)));
    let b = Mersenne31Field::from_u64_unchecked(rng.gen_range(0..((1 << 31) - 1)));
    Mersenne31Complex::new(a, b)
}

pub fn rand_fp4_from_rng<R: rand::Rng>(rng: &mut R) -> Mersenne31Quartic {
    let a = Mersenne31Field::from_u64_unchecked(rng.gen_range(0..((1 << 31) - 1)));
    let b = Mersenne31Field::from_u64_unchecked(rng.gen_range(0..((1 << 31) - 1)));
    let c = Mersenne31Field::from_u64_unchecked(rng.gen_range(0..((1 << 31) - 1)));
    let d = Mersenne31Field::from_u64_unchecked(rng.gen_range(0..((1 << 31) - 1)));
    Mersenne31Quartic::from_array_of_base([a, b, c, d])
}
