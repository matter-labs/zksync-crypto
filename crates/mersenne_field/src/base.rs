use super::*;
use crate::ops;
use core::ops::{Add, Sub};

// If we use divisions for reduction (that makes sense on M1 family and in proved environment),
// then representation of the field element is always canonical, other wise it's <= MODULUS (so zero has two representations)

#[derive(Clone, Copy, serde::Serialize, serde::Deserialize)]
#[repr(transparent)]
pub struct Mersenne31Field(pub u32);

const _: () = const {
    assert!(core::mem::size_of::<Mersenne31Field>() == core::mem::size_of::<u32>());
    assert!(core::mem::align_of::<Mersenne31Field>() == core::mem::align_of::<u32>());

    ()
};

impl Mersenne31Field {
    pub const ORDER: u32 = (1 << 31) - 1;
    pub const MSBITMASK: u32 = 1 << 31;

    #[cfg(not(feature = "use_division"))]
    #[inline(always)]
    pub const fn new(value: u32) -> Self {
        debug_assert!((value >> 31) == 0);

        Self(value)
    }

    #[cfg(feature = "use_division")]
    #[inline(always)]
    pub const fn new(value: u32) -> Self {
        debug_assert!(value < Self::ORDER);

        Self(value)
    }

    #[cfg(not(feature = "use_division"))]
    #[inline(always)]
    pub const fn to_reduced_u32(&self) -> u32 {
        // our canonical representation is 0..=modulus (31 bits full range), but not larger
        let mut c = self.0;
        if c >= Self::ORDER {
            c -= Self::ORDER;
        }
        c
    }

    #[cfg(feature = "use_division")]
    #[inline(always)]
    pub const fn to_reduced_u32(&self) -> u32 {
        self.0
    }

    #[cfg(not(feature = "use_division"))]
    pub const fn from_nonreduced_u32(c: u32) -> Self {
        let mut c = c as u32;
        if c >= Self::ORDER {
            c -= Self::ORDER;
        }
        if c >= Self::ORDER {
            c -= Self::ORDER;
        }
        Self::new(c)
    }

    #[cfg(feature = "use_division")]
    #[inline(always)]
    pub const fn from_nonreduced_u32(c: u32) -> Self {
        Self(ops::reduce_with_division(c))
    }

    pub const fn mul_2exp_u64(&self, exp: u64) -> Self {
        // In a Mersenne field, multiplication by 2^k is just a left rotation by k bits.
        let exp = (exp % 31) as u8;
        let left = (self.0 << exp) & ((1 << 31) - 1);
        let right = self.0 >> (31 - exp);
        let rotated = left | right;
        Self::new(rotated)
    }

    #[inline]
    pub const fn div_2exp_u64(&self, exp: u64) -> Self {
        // In a Mersenne field, division by 2^k is just a right rotation by k bits.
        let exp = (exp % 31) as u8;
        let left = self.0 >> exp;
        let right = (self.0 << (31 - exp)) & ((1 << 31) - 1);
        let rotated = left | right;
        Self::new(rotated)
    }

    #[inline(always)]
    pub const fn from_u62(x: u64) -> Self {
        let product_low = (x as u32) & ((1 << 31) - 1);
        let product_high = (x >> 31) as u32;
        let result = ops::add_mod(product_low, product_high);

        Self(result)
    }
}

impl Default for Mersenne31Field {
    fn default() -> Self {
        Self(0u32)
    }
}

impl PartialEq for Mersenne31Field {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.to_reduced_u32() == other.to_reduced_u32()
    }
}
impl Eq for Mersenne31Field {}

impl Hash for Mersenne31Field {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u32(self.to_reduced_u32())
    }
}

impl Ord for Mersenne31Field {
    fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        self.to_reduced_u32().cmp(&other.to_reduced_u32())
    }
}

impl PartialOrd for Mersenne31Field {
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Display for Mersenne31Field {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        Display::fmt(&self.as_u64_reduced(), f)
    }
}

impl Debug for Mersenne31Field {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        Debug::fmt(&self.as_u64_reduced(), f)
    }
}

impl Mersenne31Field {
    #[cfg(not(feature = "use_division"))]
    #[inline(always)]
    pub(crate) const fn is_zero_impl(&self) -> bool {
        self.to_reduced_u32() == 0
    }

    #[cfg(feature = "use_division")]
    #[inline(always)]
    pub(crate) const fn is_zero_impl(&self) -> bool {
        self.0 == 0
    }

    #[inline(always)]
    pub(crate) const fn exp_power_of_2_impl(&mut self, power_log: usize) {
        let mut i = 0;
        while i < power_log {
            self.square_impl();
            i += 1;
        }
    }

    pub(crate) const fn inverse_impl(&self) -> Option<Self> {
        //Since the nonzero elements of GF(pn) form a finite group with respect to multiplication,
        // a^p^n−1 = 1 (for a ≠ 0), thus the inverse of a is a^p^n−2.
        if self.is_zero_impl() {
            return None;
        }

        let mut p101 = *self;
        p101.exp_power_of_2_impl(2);
        p101.mul_assign_impl(&self);

        let mut p1111 = p101;
        p1111.square_impl();
        p1111.mul_assign_impl(&p101);

        let mut p11111111 = p1111;
        p11111111.exp_power_of_2_impl(4);
        p11111111.mul_assign_impl(&p1111);

        let mut p111111110000 = p11111111;
        p111111110000.exp_power_of_2_impl(4);

        let mut p111111111111 = p111111110000;
        p111111111111.mul_assign_impl(&p1111);

        let mut p1111111111111111 = p111111110000;
        p1111111111111111.exp_power_of_2_impl(4);
        p1111111111111111.mul_assign_impl(&p11111111);

        let mut p1111111111111111111111111111 = p1111111111111111;
        p1111111111111111111111111111.exp_power_of_2_impl(12);
        p1111111111111111111111111111.mul_assign_impl(&p111111111111);

        let mut p1111111111111111111111111111101 = p1111111111111111111111111111;
        p1111111111111111111111111111101.exp_power_of_2_impl(3);
        p1111111111111111111111111111101.mul_assign_impl(&p101);
        Some(p1111111111111111111111111111101)

        //Some(self.pow(Mersenne31Field::ORDER - 2))
    }

    pub fn sqrt(&self) -> Option<Self> {
        // p+1 = 2^31, and (p+1)/4 = 2^29
        let mut candidate = *self;
        candidate.exp_power_of_2(29);

        let mut t = candidate;
        t.square();
        if t == *self {
            Some(candidate)
        } else {
            None
        }
    }

    #[inline(always)]
    pub(crate) const fn add_assign_impl(&'_ mut self, other: &Self) -> &'_ mut Self {
        self.0 = ops::add_mod(self.0, other.0);

        self
    }

    #[inline(always)]
    pub(crate) const fn sub_assign_impl(&'_ mut self, other: &Self) -> &'_ mut Self {
        self.0 = ops::sub_mod(self.0, other.0);

        self
    }

    #[inline(always)]
    pub(crate) const fn mul_assign_impl(&'_ mut self, other: &Self) -> &'_ mut Self {
        self.0 = ops::mul_mod(self.0, other.0);
        self
    }

    #[inline(always)]
    pub(crate) const fn square_impl(&'_ mut self) -> &'_ mut Self {
        let t = *self;
        self.mul_assign_impl(&t)
    }

    #[cfg(not(feature = "use_division"))]
    #[inline(always)]
    pub(crate) const fn negate_impl(&'_ mut self) -> &'_ mut Self {
        if self.is_zero_impl() == false {
            *self = Self(Self::ORDER - self.to_reduced_u32());
        }

        self
    }

    #[cfg(all(feature = "use_division", not(feature = "modular_ops")))]
    #[inline(always)]
    pub(crate) const fn negate_impl(&'_ mut self) -> &'_ mut Self {
        *self = Self(ops::reduce_with_division(Self::ORDER.wrapping_sub(self.0)));

        self
    }

    #[cfg(feature = "modular_ops")]
    #[inline(always)]
    pub(crate) const fn negate_impl(&'_ mut self) -> &'_ mut Self {
        self.0 = ops::sub_mod(0, self.0);

        self
    }

    #[cfg(not(feature = "use_division"))]
    #[inline(always)]
    pub(crate) const fn double_impl(&'_ mut self) -> &'_ mut Self {
        let mut sum = self.0 << 1;
        let msb = sum & Self::MSBITMASK;
        sum ^= msb;
        sum += (msb != 0) as u32;
        //if sum >= Self::ORDER { sum -= Self::ORDER };
        self.0 = sum;

        self
    }

    #[cfg(feature = "use_division")]
    #[inline(always)]
    pub(crate) const fn double_impl(&'_ mut self) -> &'_ mut Self {
        let t = *self;
        self.add_assign_impl(&t);

        self
    }

    #[inline(always)]
    pub(crate) const fn mul_by_non_residue_impl(elem: &mut Self) {
        elem.negate_impl();
    }
}

impl Field for Mersenne31Field {
    const ZERO: Self = Self(0);
    const ONE: Self = Self(1);

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.is_zero_impl()
    }

    fn inverse(&self) -> Option<Self> {
        self.inverse_impl()
    }

    #[inline(always)]
    fn add_assign(&'_ mut self, other: &Self) -> &'_ mut Self {
        self.add_assign_impl(other)
    }

    #[inline(always)]
    fn sub_assign(&'_ mut self, other: &Self) -> &'_ mut Self {
        self.sub_assign_impl(other)
    }

    #[inline(always)]
    fn mul_assign(&'_ mut self, other: &Self) -> &'_ mut Self {
        self.mul_assign_impl(other)
    }

    #[inline(always)]
    fn square(&'_ mut self) -> &'_ mut Self {
        self.square_impl()
    }

    #[inline(always)]
    fn negate(&'_ mut self) -> &'_ mut Self {
        self.negate_impl()
    }

    #[inline(always)]
    fn double(&'_ mut self) -> &'_ mut Self {
        self.double_impl()
    }

    #[inline]
    fn exp_power_of_2(&mut self, power_log: usize) {
        self.exp_power_of_2_impl(power_log);
    }

    // TODO: could be optimized a little further?
    #[inline]
    fn mul_by_two(&'_ mut self) -> &'_ mut Self {
        *self = self.mul_2exp_u64(1);
        self
    }

    #[inline]
    fn div_by_two(&'_ mut self) -> &'_ mut Self {
        *self = self.div_2exp_u64(1);
        self
    }

    #[inline(always)]
    fn fused_mul_add_assign(&'_ mut self, a: &Self, b: &Self) -> &'_ mut Self {
        self.0 = ops::fma_mod(self.0, a.0, b.0);
        self
    }
}

impl Add for Mersenne31Field {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        let lhs = self;
        let rhs = rhs;
        let mut res = lhs;
        res.add_assign(&rhs);
        res
    }
}

impl Sub for Mersenne31Field {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        let lhs = self;
        let rhs = rhs;
        let mut res = lhs;
        res.sub_assign(&rhs);
        res
    }
}

impl PrimeField for Mersenne31Field {
    const TWO: Self = Self(2);
    const MINUS_ONE: Self = Self(Self::ORDER - 1);
    const NUM_BYTES_IN_REPR: usize = 4;
    const CHAR_BITS: usize = 31;
    const CHARACTERISTICS: u64 = Self::ORDER as u64;

    #[inline(always)]
    fn as_u64(self) -> u64 {
        self.0 as u64
    }

    #[inline(always)]
    fn from_u64_unchecked(value: u64) -> Self {
        Self::new(value.try_into().expect("Too large"))
    }
    #[inline(always)]
    fn from_u64(value: u64) -> Option<Self> {
        if value as u32 >= Self::ORDER {
            None
        } else {
            Some(Self(value as u32))
        }
    }

    #[inline(always)]
    fn from_u64_with_reduction(value: u64) -> Self {
        // p = 2^31 - 1, so 2^32 = 2p + 2 = 2
        Self((value % Self::ORDER as u64) as u32)
    }

    #[inline(always)]
    fn as_u64_reduced(&self) -> u64 {
        self.to_reduced_u32() as u64
    }

    #[track_caller]
    fn as_boolean(&self) -> bool {
        let as_uint = self.to_reduced_u32();
        assert!(
            as_uint == 0 || as_uint == 1,
            "expected boolean value, got {}",
            as_uint
        );

        as_uint != 0
    }

    fn from_boolean(flag: bool) -> Self {
        if flag {
            Self::ONE
        } else {
            Self::ZERO
        }
    }

    fn to_le_bytes(self) -> [u8; Self::NUM_BYTES_IN_REPR] {
        self.0.to_le_bytes()
    }

    fn increment_unchecked(&'_ mut self) {
        self.0 += 1;
    }
}

impl BaseField for Mersenne31Field {
    const QUADRATIC_NON_RESIDUE: Mersenne31Field = Mersenne31Field::MINUS_ONE;

    #[inline(always)]
    fn mul_by_non_residue(elem: &mut Self) {
        Self::mul_by_non_residue_impl(elem);
    }
}

#[cfg(test)]
mod tests {
    use std::hash::DefaultHasher;

    use super::*;

    #[test]
    // New assumes that u32 is inside the field.
    fn test_new() {
        let a = Mersenne31Field::new(0);
        assert_eq!(a.0, 0);

        let b = Mersenne31Field::new(Mersenne31Field::ORDER - 1);
        assert_eq!(b.0, Mersenne31Field::ORDER - 1);
    }
    #[test]
    fn test_from_nonreduced_u32() {
        let a = Mersenne31Field::from_nonreduced_u32(Mersenne31Field::ORDER);
        assert_eq!(a.0, 0);

        let b = Mersenne31Field::from_nonreduced_u32(Mersenne31Field::ORDER + 1);
        assert_eq!(b.0, 1);

        let c = Mersenne31Field::from_nonreduced_u32(2 * Mersenne31Field::ORDER - 1);
        assert_eq!(c.0, Mersenne31Field::ORDER - 1);

        let d = Mersenne31Field::from_nonreduced_u32(2 * Mersenne31Field::ORDER);
        assert_eq!(d.0, 0);

        let e = Mersenne31Field::from_nonreduced_u32(u32::MAX);
        assert_eq!(e.0, 1);
    }

    #[test]
    fn test_two_zeros() {
        let d = Mersenne31Field::new(Mersenne31Field::ORDER - 1);
        let e = Mersenne31Field::new(1);
        let f = d + e;
        // Zero here can be represented as 0 or 2^31 - 1.
        assert_eq!(f.to_reduced_u32(), 0);
        assert!(f.is_zero());


        let h = Mersenne31Field::default();
        assert!(h.is_zero());

        assert_eq!(h, f);

        let h = Mersenne31Field::default();
        assert_eq!(h, f);

        let mut h_hasher = DefaultHasher::default();
        h.hash(&mut h_hasher);

        let mut f_hasher = DefaultHasher::default();
        f.hash(&mut f_hasher);

        assert_eq!(f_hasher.finish(), h_hasher.finish());

        let one = Mersenne31Field::ONE;

        assert!(f < one);
        assert!(h < one);

    }
    
 
    #[test]
    fn test_add() {
        let a = Mersenne31Field::new(1);
        let b = Mersenne31Field::new(2);
        let c = a + b;
        assert_eq!(c.0, 3);

        let d = Mersenne31Field::new(Mersenne31Field::ORDER - 1);
        let e = Mersenne31Field::new(1);
        let f = d + e;
        // Zero here can be represented as 0 or 2^31 - 1.
        assert_eq!(f.to_reduced_u32(), 0);
        assert!(f.is_zero_impl());
        // But adding +1 to 0 should give 1.
        let g = f + e;
        assert_eq!(g.0, 1);
        assert!(g.is_zero_impl() == false);
    }

    #[test]
    fn test_sub() {
        let a = Mersenne31Field::new(3);
        let b = Mersenne31Field::new(2);
        let c = a - b;
        assert_eq!(c.0, 1);

        let d = Mersenne31Field::new(0);
        let e = Mersenne31Field::new(1);
        let f = d - e;
        assert_eq!(f.0, Mersenne31Field::ORDER - 1);
    }
 
    #[test]
    fn test_mul() {
        let mut a = Mersenne31Field::new(2);
        let b = Mersenne31Field::new(3);
        a.mul_assign(&b);
        assert_eq!(a.0, 6);
        
        let mut d = Mersenne31Field::new(Mersenne31Field::ORDER - 1);
        let e = Mersenne31Field::new(2);
        d.mul_assign(&e);
        assert_eq!(d.0, Mersenne31Field::ORDER - 2);
    }
 
    #[test]
    fn test_inverse() {
        let mut a = Mersenne31Field::new(3);
        let inv_a = a.inverse().unwrap();
        let one = a.mul_assign(&inv_a);
        assert_eq!(one.0, 1);

        let zero = Mersenne31Field::new(0);
        assert!(zero.inverse().is_none());
    }
 
    #[test]
    fn test_sqrt() {
        let a = Mersenne31Field::new(4);
        let sqrt_a = a.sqrt().unwrap();
        assert_eq!(sqrt_a.0, 2);

        let b = Mersenne31Field::new(5);
        assert!(b.sqrt().is_none());
    }
 
    #[test]
    fn test_large_numbers() {
        let mut a = Mersenne31Field::new(Mersenne31Field::ORDER - 1);
        let b = Mersenne31Field::new(Mersenne31Field::ORDER - 2);
        let c = a + b;
        assert_eq!(c.0, Mersenne31Field::ORDER - 3);

        let d = a.mul_assign(&b);
        assert_eq!(d.0, 2);
    }

    #[test]
    fn test_mul_2exp_u64_identity() {
        // Check that rotating by 31 bits (full rotation) returns the original element.
        let a = Mersenne31Field::new(0x1ABCDEF);
        assert_eq!(a.mul_2exp_u64(31).0, a.0);
        assert_eq!(a.mul_2exp_u64(0).0, a.0);
        assert_eq!(a.mul_2exp_u64(31 * 800_000).0, a.0);
    }

    #[test]
    fn test_mul_2exp_u64_rotation() {
        // For a rotation by 1 bit, compute the expected value manually.
        let a = Mersenne31Field::new(0x1234567);
        let expected = 0x2468ace;
        assert_eq!(a.mul_2exp_u64(1).0, expected);
    }

    #[test]
    fn test_div_2exp_u64_identity() {
        // Check that rotating right by 31 bits (full rotation) returns the original element.
        let a = Mersenne31Field::new(0x1ABCDEF);
        assert_eq!(a.div_2exp_u64(31).0, a.0);
        assert_eq!(a.div_2exp_u64(0).0, a.0);
        assert_eq!(a.div_2exp_u64(31 * 800_000).0, a.0);
    }

    #[test]
    fn test_div_2exp_u64_rotation() {
        // For a rotation by 1 bit to the right, compute the expected value manually.
        let a = Mersenne31Field::new(0x2468ace);

        assert_eq!(a.div_2exp_u64(1).0, 0x1234567);
    }

    #[test]
    fn test_from_u62() {
        // For x = 1 << 31, we have:
        // product_low = 0, product_high = 1, so result = add_mod(0, 1) = 1.
        let a = Mersenne31Field::from_u62(1u64 << 31);
        assert_eq!(a.0, 1);

        // For x less than 2^31, from_u62 should act as the identity.
        let b = Mersenne31Field::from_u62(42);
        assert_eq!(b.0, 42);

        // For a composite x = (high << 31) + low.
        // Let high = 10 and low = 100.
        let x = (10u64 << 31) + 100;
        let expected = 110;
        let c = Mersenne31Field::from_u62(x);
        assert_eq!(c.0, expected);
    }

    #[test]
    fn test_mul_by_two() {
        // Test that multiplying by two is equivalent to addition with itself.
        let a = Mersenne31Field::new(7);
        let mut b = a;
        b.mul_by_two();
        let expected = a + a;
        assert_eq!(b.to_reduced_u32(), expected.to_reduced_u32());
        assert_eq!(b.to_reduced_u32(), 14);

    }

    #[test]
    fn test_div_by_two() {
        // Test that dividing by two is the inverse of multiplying by two.
        let a = Mersenne31Field::new(14);
        let mut b = a;
        b.div_by_two();
        // Since 14 / 2 = 7 modulo ORDER
        let expected = Mersenne31Field::new(7);
        assert_eq!(b.to_reduced_u32(), expected.to_reduced_u32());

        b.div_by_two();
        // 1073741827 * 2 % ORDER = 7
        assert_eq!(b.to_reduced_u32(), 1073741827);
    }

    #[test]
    fn test_fused_mul_add_assign() {
        // Test fused multiply-add: result = self + a * b.
        let a = Mersenne31Field::new(2);
        let b = Mersenne31Field::new(3);
        let mut result = Mersenne31Field::new(4);
        result.fused_mul_add_assign(&a, &b);
        // Expected: 4 + (2 * 3) = 10.
        let expected = Mersenne31Field::new(10);
        assert_eq!(result.to_reduced_u32(), expected.to_reduced_u32());
    }

    #[test]
    fn test_negate() {
        // Test negate for a nonzero element.
        let mut a = Mersenne31Field::new(5);
        a.negate();
        let expected = Mersenne31Field::new(Mersenne31Field::ORDER - 5);
        assert_eq!(a.to_reduced_u32(), expected.to_reduced_u32());

        // Test negate on zero remains zero.
        let mut zero = Mersenne31Field::new(0);
        zero.negate();
        assert_eq!(zero.to_reduced_u32(), 0);
    }

    #[test]
    fn test_double() {
        // Test double for an element.
        let mut a = Mersenne31Field::new(6);
        a.double();
        assert_eq!(a.to_reduced_u32(), 12);

        // Test double when wrapping around the modulus.
        let mut b = Mersenne31Field::new(Mersenne31Field::ORDER - 1);
        b.double();
        // (ORDER - 1) * 2 mod ORDER should equal ORDER - 2.
        let expected = Mersenne31Field::new(Mersenne31Field::ORDER - 2);
        assert_eq!(b.to_reduced_u32(), expected.to_reduced_u32());
    }

}
