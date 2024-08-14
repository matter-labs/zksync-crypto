use ff::{adc, sbb, mac_with_carry};
use ff::{Field, PrimeField, SqrtField, PrimeFieldRepr, PrimeFieldDecodingError, LegendreSymbol};
use ff::LegendreSymbol::*;

// s = 6554484396890773809930967563523245729705921265872317281365359162392183254199
const MODULUS: FsRepr = FsRepr([0xd0970e5ed6f72cb7, 0xa6682093ccc81082, 0x6673b0101343b00, 0xe7db4ea6533afa9]);

// The number of bits needed to represent the modulus.
const MODULUS_BITS: u32 = 252;

// The number of bits that must be shaved from the beginning of
// the representation when randomly sampling.
const REPR_SHAVE_BITS: u32 = 4;

// R = 2**256 % s
const R: FsRepr = FsRepr([0x25f80bb3b99607d9, 0xf315d62f66b6e750, 0x932514eeeb8814f4, 0x9a6fc6f479155c6]);

// R2 = R^2 % s
const R2: FsRepr = FsRepr([0x67719aa495e57731, 0x51b0cef09ce3fc26, 0x69dab7fac026e9a5, 0x4f6547b8d127688]);

// INV = -(s^{-1} mod 2^64) mod s
const INV: u64 = 0x1ba3a358ef788ef9;

// GENERATOR = 6 (multiplicative generator of r-1 order, that is also quadratic nonresidue)
const GENERATOR: FsRepr = FsRepr([0x720b1b19d49ea8f1, 0xbf4aa36101f13a58, 0x5fa8cc968193ccbb, 0xe70cbdc7dccf3ac]);

// 2^S * t = MODULUS - 1 with t odd
const S: u32 = 1;

// 2^S root of unity computed by GENERATOR^t
const ROOT_OF_UNITY: FsRepr = FsRepr([0xaa9f02ab1d6124de, 0xb3524a6466112932, 0x7342261215ac260b, 0x4d6b87b1da259e2]);

// -((2**256) mod s) mod s
const NEGATIVE_ONE: Fs = Fs(FsRepr([0xaa9f02ab1d6124de, 0xb3524a6466112932, 0x7342261215ac260b, 0x4d6b87b1da259e2]));

/// This is the underlying representation of an element of `Fs`.
#[derive(Copy, Clone, PartialEq, Eq, Default, Debug, Hash)]
pub struct FsRepr(pub [u64; 4]);

impl ::rand::Rand for FsRepr {
    #[inline(always)]
    fn rand<R: ::rand::Rng>(rng: &mut R) -> Self {
        FsRepr(rng.gen())
    }
}

impl ::std::fmt::Display for FsRepr
{
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "0x")?;
        for i in self.0.iter().rev() {
            write!(f, "{:016x}", *i)?;
        }

        Ok(())
    }
}

impl AsRef<[u64]> for FsRepr {
    #[inline(always)]
    fn as_ref(&self) -> &[u64] {
        &self.0
    }
}

impl AsMut<[u64]> for FsRepr {
    #[inline(always)]
    fn as_mut(&mut self) -> &mut [u64] {
        &mut self.0
    }
}

impl From<u64> for FsRepr {
    #[inline(always)]
    fn from(val: u64) -> FsRepr {
        let mut repr = Self::default();
        repr.0[0] = val;
        repr
    }
}

impl Ord for FsRepr {
    #[inline(always)]
    fn cmp(&self, other: &FsRepr) -> ::std::cmp::Ordering {
        for (a, b) in self.0.iter().rev().zip(other.0.iter().rev()) {
            if a < b {
                return ::std::cmp::Ordering::Less
            } else if a > b {
                return ::std::cmp::Ordering::Greater
            }
        }

        ::std::cmp::Ordering::Equal
    }
}

impl PartialOrd for FsRepr {
    #[inline(always)]
    fn partial_cmp(&self, other: &FsRepr) -> Option<::std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PrimeFieldRepr for FsRepr {
    #[inline(always)]
    fn is_odd(&self) -> bool {
        self.0[0] & 1 == 1
    }

    #[inline(always)]
    fn is_even(&self) -> bool {
        !self.is_odd()
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.0.iter().all(|&e| e == 0)
    }

    #[inline(always)]
    fn shr(&mut self, mut n: u32) {
        if n >= 64 * 4 {
            *self = Self::from(0);
            return;
        }

        while n >= 64 {
            let mut t = 0;
            for i in self.0.iter_mut().rev() {
                ::std::mem::swap(&mut t, i);
            }
            n -= 64;
        }

        if n > 0 {
            let mut t = 0;
            for i in self.0.iter_mut().rev() {
                let t2 = *i << (64 - n);
                *i >>= n;
                *i |= t;
                t = t2;
            }
        }
    }

    #[inline(always)]
    fn div2(&mut self) {
        let mut t = 0;
        for i in self.0.iter_mut().rev() {
            let t2 = *i << 63;
            *i >>= 1;
            *i |= t;
            t = t2;
        }
    }

    #[inline(always)]
    fn mul2(&mut self) {
        let mut last = 0;
        for i in &mut self.0 {
            let tmp = *i >> 63;
            *i <<= 1;
            *i |= last;
            last = tmp;
        }
    }

    #[inline(always)]
    fn shl(&mut self, mut n: u32) {
        if n >= 64 * 4 {
            *self = Self::from(0);
            return;
        }

        while n >= 64 {
            let mut t = 0;
            for i in &mut self.0 {
                ::std::mem::swap(&mut t, i);
            }
            n -= 64;
        }

        if n > 0 {
            let mut t = 0;
            for i in &mut self.0 {
                let t2 = *i >> (64 - n);
                *i <<= n;
                *i |= t;
                t = t2;
            }
        }
    }

    #[inline(always)]
    fn num_bits(&self) -> u32 {
        let mut ret = (4 as u32) * 64;
        for i in self.0.iter().rev() {
            let leading = i.leading_zeros();
            ret -= leading;
            if leading != 64 {
                break;
            }
        }

        ret
    }

    #[inline(always)]
    fn add_nocarry(&mut self, other: &FsRepr) {
        let mut carry = 0;

        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            *a = adc(*a, *b, &mut carry);
        }
    }

    #[inline(always)]
    fn sub_noborrow(&mut self, other: &FsRepr) {
        let mut borrow = 0;

        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            *a = sbb(*a, *b, &mut borrow);
        }
    }
}

/// This is an element of the scalar field of the Jubjub curve.
#[derive(Copy, Clone, PartialEq, Eq, Debug, Hash, Default)]
pub struct Fs(FsRepr);

impl ::std::fmt::Display for Fs
{
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fs({})", self.into_repr())
    }
}

impl ::rand::Rand for Fs {
    fn rand<R: ::rand::Rng>(rng: &mut R) -> Self {
        loop {
            let mut tmp = Fs(FsRepr::rand(rng));

            // Mask away the unused bits at the beginning.
            tmp.0.as_mut()[3] &= 0xffffffffffffffff >> REPR_SHAVE_BITS;

            if tmp.is_valid() {
                return tmp
            }
        }
    }
}

impl From<Fs> for FsRepr {
    fn from(e: Fs) -> FsRepr {
        e.into_repr()
    }
}

impl PrimeField for Fs {
    type Repr = FsRepr;

    fn from_repr(r: FsRepr) -> Result<Fs, PrimeFieldDecodingError> {
        let mut r = Fs(r);
        if r.is_valid() {
            r.mul_assign(&Fs(R2));

            Ok(r)
        } else {
            Err(PrimeFieldDecodingError::NotInField(format!("{}", r.0)))
        }
    }

    fn from_raw_repr(r: FsRepr) -> Result<Fs, PrimeFieldDecodingError> {
        let r = Fs(r);
        if r.is_valid() {
            Ok(r)
        } else {
            Err(PrimeFieldDecodingError::NotInField(format!("{}", r.0)))
        }
    }

    fn into_repr(&self) -> FsRepr {
        let mut r = *self;
        r.mont_reduce((self.0).0[0], (self.0).0[1],
                      (self.0).0[2], (self.0).0[3],
                      0, 0, 0, 0);
        r.0
    }

    fn into_raw_repr(&self) -> FsRepr {
        let r = *self;
        r.0
    }

    fn char() -> FsRepr {
        MODULUS
    }

    const NUM_BITS: u32 = MODULUS_BITS;

    const CAPACITY: u32 = Self::NUM_BITS - 1;

    fn multiplicative_generator() -> Self {
        Fs(GENERATOR)
    }

    const S: u32 = S;

    fn root_of_unity() -> Self {
        Fs(ROOT_OF_UNITY)
    }
}

impl Field for Fs {
    #[inline]
    fn zero() -> Self {
        Fs(FsRepr::from(0))
    }

    #[inline]
    fn one() -> Self {
        Fs(R)
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    #[inline]
    fn add_assign(&mut self, other: &Fs) {
        // This cannot exceed the backing capacity.
        self.0.add_nocarry(&other.0);

        // However, it may need to be reduced.
        self.reduce();
    }

    #[inline]
    fn double(&mut self) {
        // This cannot exceed the backing capacity.
        self.0.mul2();

        // However, it may need to be reduced.
        self.reduce();
    }

    #[inline]
    fn sub_assign(&mut self, other: &Fs) {
        // If `other` is larger than `self`, we'll need to add the modulus to self first.
        if other.0 > self.0 {
            self.0.add_nocarry(&MODULUS);
        }

        self.0.sub_noborrow(&other.0);
    }

    #[inline]
    fn negate(&mut self) {
        if !self.is_zero() {
            let mut tmp = MODULUS;
            tmp.sub_noborrow(&self.0);
            self.0 = tmp;
        }
    }

    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            // Guajardo Kumar Paar Pelzl
            // Efficient Software-Implementation of Finite Fields with Applications to Cryptography
            // Algorithm 16 (BEA for Inversion in Fp)

            let one = FsRepr::from(1);

            let mut u = self.0;
            let mut v = MODULUS;
            let mut b = Fs(R2); // Avoids unnecessary reduction step.
            let mut c = Self::zero();

            while u != one && v != one {
                while u.is_even() {
                    u.div2();

                    if b.0.is_even() {
                        b.0.div2();
                    } else {
                        b.0.add_nocarry(&MODULUS);
                        b.0.div2();
                    }
                }

                while v.is_even() {
                    v.div2();

                    if c.0.is_even() {
                        c.0.div2();
                    } else {
                        c.0.add_nocarry(&MODULUS);
                        c.0.div2();
                    }
                }

                if v < u {
                    u.sub_noborrow(&v);
                    b.sub_assign(&c);
                } else {
                    v.sub_noborrow(&u);
                    c.sub_assign(&b);
                }
            }

            if u == one {
                Some(b)
            } else {
                Some(c)
            }
        }
    }

    #[inline(always)]
    fn frobenius_map(&mut self, _: usize) {
        // This has no effect in a prime field.
    }

    #[inline]
    fn mul_assign(&mut self, other: &Fs)
    {
        let mut carry = 0;
        let r0 = mac_with_carry(0, (self.0).0[0], (other.0).0[0], &mut carry);
        let r1 = mac_with_carry(0, (self.0).0[0], (other.0).0[1], &mut carry);
        let r2 = mac_with_carry(0, (self.0).0[0], (other.0).0[2], &mut carry);
        let r3 = mac_with_carry(0, (self.0).0[0], (other.0).0[3], &mut carry);
        let r4 = carry;
        let mut carry = 0;
        let r1 = mac_with_carry(r1, (self.0).0[1], (other.0).0[0], &mut carry);
        let r2 = mac_with_carry(r2, (self.0).0[1], (other.0).0[1], &mut carry);
        let r3 = mac_with_carry(r3, (self.0).0[1], (other.0).0[2], &mut carry);
        let r4 = mac_with_carry(r4, (self.0).0[1], (other.0).0[3], &mut carry);
        let r5 = carry;
        let mut carry = 0;
        let r2 = mac_with_carry(r2, (self.0).0[2], (other.0).0[0], &mut carry);
        let r3 = mac_with_carry(r3, (self.0).0[2], (other.0).0[1], &mut carry);
        let r4 = mac_with_carry(r4, (self.0).0[2], (other.0).0[2], &mut carry);
        let r5 = mac_with_carry(r5, (self.0).0[2], (other.0).0[3], &mut carry);
        let r6 = carry;
        let mut carry = 0;
        let r3 = mac_with_carry(r3, (self.0).0[3], (other.0).0[0], &mut carry);
        let r4 = mac_with_carry(r4, (self.0).0[3], (other.0).0[1], &mut carry);
        let r5 = mac_with_carry(r5, (self.0).0[3], (other.0).0[2], &mut carry);
        let r6 = mac_with_carry(r6, (self.0).0[3], (other.0).0[3], &mut carry);
        let r7 = carry;
        self.mont_reduce(r0, r1, r2, r3, r4, r5, r6, r7);
    }

    #[inline]
    fn square(&mut self)
    {
        let mut carry = 0;
        let r1 = mac_with_carry(0, (self.0).0[0], (self.0).0[1], &mut carry);
        let r2 = mac_with_carry(0, (self.0).0[0], (self.0).0[2], &mut carry);
        let r3 = mac_with_carry(0, (self.0).0[0], (self.0).0[3], &mut carry);
        let r4 = carry;
        let mut carry = 0;
        let r3 = mac_with_carry(r3, (self.0).0[1], (self.0).0[2], &mut carry);
        let r4 = mac_with_carry(r4, (self.0).0[1], (self.0).0[3], &mut carry);
        let r5 = carry;
        let mut carry = 0;
        let r5 = mac_with_carry(r5, (self.0).0[2], (self.0).0[3], &mut carry);
        let r6 = carry;

        let r7 = r6 >> 63;
        let r6 = (r6 << 1) | (r5 >> 63);
        let r5 = (r5 << 1) | (r4 >> 63);
        let r4 = (r4 << 1) | (r3 >> 63);
        let r3 = (r3 << 1) | (r2 >> 63);
        let r2 = (r2 << 1) | (r1 >> 63);
        let r1 = r1 << 1;

        let mut carry = 0;
        let r0 = mac_with_carry(0, (self.0).0[0], (self.0).0[0], &mut carry);
        let r1 = adc(r1, 0, &mut carry);
        let r2 = mac_with_carry(r2, (self.0).0[1], (self.0).0[1], &mut carry);
        let r3 = adc(r3, 0, &mut carry);
        let r4 = mac_with_carry(r4, (self.0).0[2], (self.0).0[2], &mut carry);
        let r5 = adc(r5, 0, &mut carry);
        let r6 = mac_with_carry(r6, (self.0).0[3], (self.0).0[3], &mut carry);
        let r7 = adc(r7, 0, &mut carry);
        self.mont_reduce(r0, r1, r2, r3, r4, r5, r6, r7);
    }
}

impl Fs {
    /// Determines if the element is really in the field. This is only used
    /// internally.
    #[inline(always)]
    pub fn is_valid(&self) -> bool {
        self.0 < MODULUS
    }

    /// Subtracts the modulus from this element if this element is not in the
    /// field. Only used internally.
    #[inline(always)]
    pub fn reduce(&mut self) {
        if !self.is_valid() {
            self.0.sub_noborrow(&MODULUS);
        }
    }

    #[inline(always)]
    fn mont_reduce(
        &mut self,
        r0: u64,
        mut r1: u64,
        mut r2: u64,
        mut r3: u64,
        mut r4: u64,
        mut r5: u64,
        mut r6: u64,
        mut r7: u64
    )
    {
        // The Montgomery reduction here is based on Algorithm 14.32 in
        // Handbook of Applied Cryptography
        // <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.

        let k = r0.wrapping_mul(INV);
        let mut carry = 0;
        mac_with_carry(r0, k, MODULUS.0[0], &mut carry);
        r1 = mac_with_carry(r1, k, MODULUS.0[1], &mut carry);
        r2 = mac_with_carry(r2, k, MODULUS.0[2], &mut carry);
        r3 = mac_with_carry(r3, k, MODULUS.0[3], &mut carry);
        r4 = adc(r4, 0, &mut carry);
        let carry2 = carry;
        let k = r1.wrapping_mul(INV);
        let mut carry = 0;
        mac_with_carry(r1, k, MODULUS.0[0], &mut carry);
        r2 = mac_with_carry(r2, k, MODULUS.0[1], &mut carry);
        r3 = mac_with_carry(r3, k, MODULUS.0[2], &mut carry);
        r4 = mac_with_carry(r4, k, MODULUS.0[3], &mut carry);
        r5 = adc(r5, carry2, &mut carry);
        let carry2 = carry;
        let k = r2.wrapping_mul(INV);
        let mut carry = 0;
        mac_with_carry(r2, k, MODULUS.0[0], &mut carry);
        r3 = mac_with_carry(r3, k, MODULUS.0[1], &mut carry);
        r4 = mac_with_carry(r4, k, MODULUS.0[2], &mut carry);
        r5 = mac_with_carry(r5, k, MODULUS.0[3], &mut carry);
        r6 = adc(r6, carry2, &mut carry);
        let carry2 = carry;
        let k = r3.wrapping_mul(INV);
        let mut carry = 0;
        mac_with_carry(r3, k, MODULUS.0[0], &mut carry);
        r4 = mac_with_carry(r4, k, MODULUS.0[1], &mut carry);
        r5 = mac_with_carry(r5, k, MODULUS.0[2], &mut carry);
        r6 = mac_with_carry(r6, k, MODULUS.0[3], &mut carry);
        r7 = adc(r7, carry2, &mut carry);
        (self.0).0[0] = r4;
        (self.0).0[1] = r5;
        (self.0).0[2] = r6;
        (self.0).0[3] = r7;
        self.reduce();
    }

    fn branchless_reduction(&mut self) {
        let mut borrow = 0;
        let mut tmp = *self;

        for (a, &b) in (tmp.0).0.iter_mut().zip(MODULUS.0.iter()) {
            *a = sbb(*a, b, &mut borrow);
        }

        for (target, reduced) in (self.0).0.iter_mut().zip((tmp.0).0.iter()) {
            // if there was a borrow then 
            // - it's equal to 2^64 - 1
            // - we take original value
            // otherwise 
            // - bitflip borrow
            // take reduced value

            *target = (*target & borrow) | (*reduced & (!borrow));
        }
    }

    // fn mul_bits<S: AsRef<[u64]>>(&self, bits: BitIterator<S>) -> Self {
    //     let mut res = Self::zero();
    //     for bit in bits {
    //         res.double();

    //         if bit {
    //             res.add_assign(self)
    //         }
    //     }
    //     res
    // }

    pub fn rps_mul_assign(&mut self, other: &Fs) {
        let [b0, b1, b2, b3] = (other.0).0;

        let a0 = (self.0).0[0];
        let a1 = (self.0).0[1];

        // two temporary registers
        let (t1, t0) = full_width_mul(a1, b1);
        // make product a0 * b1, with low part going into z1 (which is empty)
        // and high part is summed with t0 and propagated to t1
        let (z1, t0, t1) = mul_and_add_existing_high(a0, b1, t0, t1);
        // make product a0 * b1, with low part beign added with z1, and then propagated
        // and high part is summed with t0 and propagated to t1
        let (z1, t0, t1) = mul_and_add_existing(a1, b0, z1, t0, t1);
        // make product a0 * b0, and propagate everything
        let (z0, z1, z2, z3, c0) = mul_and_full_block_propagate(a0, b0, 0, z1, 0, 0, t0, t1);

        // round 2

        let (t1, t0) = full_width_mul(a1, b3);
        let (z3, t0, t1) = mul_and_add_existing(a0, b3, z3, t0, t1);
        let (z3, t0, t1) = mul_and_add_existing(a1, b2, z3, t0, t1);
        // we place c0 instead of empty (yet) z4
        let (z2, z3, z4, z5, c1) = mul_and_full_block_propagate(a0, b2, z2, z3, c0, 0, t0, t1);

        // round 3

        drop(a0);
        drop(a1);

        let a2 = (self.0).0[2];
        let a3 = (self.0).0[3];

        let (t1, t0) = full_width_mul(a3, b1);
        let (z3, t0, t1) = mul_and_add_existing(a2, b1, z3, t0, t1);
        let (z3, t0, t1) = mul_and_add_existing(a3, b0, z3, t0, t1);
        let (z2, z3, z4, z5, c1) = mul_and_full_block_propagate_into_existing_carry_catcher(a0, b2, z2, z3, z4, z5, t0, t1, c1);

        // round 4

        let (t1, t0) = full_width_mul(a3, b3);
        let (z5, t0, t1) = mul_and_add_existing(a3, b2, z5, t0, t1);
        let (z5, t0, t1) = mul_and_add_existing(a2, b3, z5, t0, t1);
        let (z4, z5, z6, z7) = mul_and_full_block_propagate_without_carry_catch(a2, b2, z4, z5, c1, 0, t0, t1);

        self.mont_reduce(z0, z1, z2, z3, z4, z5, z6, z7);
    }

    pub fn optimistic_cios_mul_assign(&mut self, other: &Fs) {
        let mut m;

        let [q0, q1, q2, q3] = MODULUS.0;
        let [a0, a1, a2, a3] = (self.0).0;

        // round 0 
        let b0 = (other.0).0[0];
        let (r0, carry) = full_width_mul(a0, b0);
        m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = wrapping_mac_by_value(r0, m, q0);

        // loop over the rest
        let (r1, carry) = mac_by_value(carry, a1, b0);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, q1, red_carry);

        let (r2, carry) = mac_by_value(carry, a2, b0);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, q2, red_carry);

        let (r3, carry) = mac_by_value(carry, a3, b0);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, q3, red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(b0);

        // round 1
        let b1 = (other.0).0[1];
        let (r0, carry) = mac_by_value(r0, a0, b1);
        m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = wrapping_mac_by_value(r0, m, q0);

        // loop over the rest
        let (r1, carry) = mac_with_carry_by_value(r1, carry, a1, b1);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, q1, red_carry);

        let (r2, carry) = mac_with_carry_by_value(r2, carry, a2, b1);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, q2, red_carry);

        let (r3, carry) = mac_with_carry_by_value(r3, carry, a3, b1);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, q3, red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(b1);

        // round 2
        let b2 = (other.0).0[2];
        let (r0, carry) = mac_by_value(r0, a0, b2);
        m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = wrapping_mac_by_value(r0, m, q0);

        // loop over the rest
        let (r1, carry) = mac_with_carry_by_value(r1, carry, a1, b2);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, q1, red_carry);

        let (r2, carry) = mac_with_carry_by_value(r2, carry, a2, b2);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, q2, red_carry);

        let (r3, carry) = mac_with_carry_by_value(r3, carry, a3, b2);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, q3, red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(b2);

        // round 3
        let b3 = (other.0).0[3];
        let (r0, carry) = mac_by_value(r0, a0, b3);
        m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = wrapping_mac_by_value(r0, m, q0);

        // loop over the rest
        let (r1, carry) = mac_with_carry_by_value(r1, carry, a1, b3);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, q1, red_carry);

        let (r2, carry) = mac_with_carry_by_value(r2, carry, a2, b3);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, q2, red_carry);

        let (r3, carry) = mac_with_carry_by_value(r3, carry, a3, b3);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, q3, red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(b3);

        (self.0).0 = [r0, r1, r2, r3];
        self.reduce();
    }

    pub fn optimistic_cios_mul_assign_with_different_semantics(&mut self, other: &Fs) {
        let mut m;

        // round 0 
        let b0 = (other.0).0[0];
        let (r0, carry) = full_width_mul((self.0).0[0], (other.0).0[0]);
        m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = wrapping_mac_by_value(r0, m, MODULUS.0[1]);

        // loop over the rest
        let (r1, carry) = mac_by_value(carry, (self.0).0[1], b0);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, MODULUS.0[1], red_carry);

        let (r2, carry) = mac_by_value(carry, (self.0).0[2], b0);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, MODULUS.0[2], red_carry);

        let (r3, carry) = mac_by_value(carry, (self.0).0[3], b0);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, MODULUS.0[3], red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(b0);

        // round 1
        let b1 = (other.0).0[1];
        let (r0, carry) = mac_by_value(r0, (self.0).0[0], b1);
        m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = wrapping_mac_by_value(r0, m, MODULUS.0[0]);

        // loop over the rest
        let (r1, carry) = mac_with_carry_by_value(r1, carry, (self.0).0[1], b1);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, MODULUS.0[1], red_carry);

        let (r2, carry) = mac_with_carry_by_value(r2, carry, (self.0).0[2], b1);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, MODULUS.0[2], red_carry);

        let (r3, carry) = mac_with_carry_by_value(r3, carry, (self.0).0[3], b1);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, MODULUS.0[3], red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(b1);

        // round 2
        let b2 = (other.0).0[2];
        let (r0, carry) = mac_by_value(r0, (self.0).0[0], b2);
        m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = wrapping_mac_by_value(r0, m, MODULUS.0[0]);

        // loop over the rest
        let (r1, carry) = mac_with_carry_by_value(r1, carry, (self.0).0[1], b2);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, MODULUS.0[1], red_carry);

        let (r2, carry) = mac_with_carry_by_value(r2, carry, (self.0).0[2], b2);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, MODULUS.0[2], red_carry);

        let (r3, carry) = mac_with_carry_by_value(r3, carry, (self.0).0[3], b2);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, MODULUS.0[3], red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(b2);

        // round 3
        let b3 = (other.0).0[3];
        let (r0, carry) = mac_by_value(r0, (self.0).0[0], b3);
        m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = wrapping_mac_by_value(r0, m, MODULUS.0[0]);

        // loop over the rest
        let (r1, carry) = mac_with_carry_by_value(r1, carry, (self.0).0[1], b3);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, MODULUS.0[1], red_carry);

        let (r2, carry) = mac_with_carry_by_value(r2, carry, (self.0).0[2], b3);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, MODULUS.0[2], red_carry);

        let (r3, carry) = mac_with_carry_by_value(r3, carry, (self.0).0[3], b3);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, MODULUS.0[3], red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(b3);

        (self.0).0 = [r0, r1, r2, r3];
        self.reduce();
    }

    pub fn optimistic_cios_by_value(self, other: Fs) -> Self {
        let [q0, q1, q2, q3] = MODULUS.0;
        let [b0, b1, b2, b3] = (other.0).0;

        // round 0 
        let a0 = (self.0).0[0];
        let (r0, carry) = full_width_mul(a0, b0);
        let m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = mac_by_value_return_carry_only(r0, m, q0);

        // loop over the rest
        let (r1, carry) = mac_by_value(carry, a0, b1);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, q1, red_carry);

        let (r2, carry) = mac_by_value(carry, a0, b2);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, q2, red_carry);

        let (r3, carry) = mac_by_value(carry, a0, b3);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, q3, red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(a0);
        drop(m);

        // round 1
        let a1 = (self.0).0[1];
        let (r0, carry) = mac_by_value(r0, a1, b0);
        let m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = mac_by_value_return_carry_only(r0, m, q0);

        // loop over the rest
        let (r1, carry) = mac_with_carry_by_value(r1, a1, b1, carry);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, q1, red_carry);

        let (r2, carry) = mac_with_carry_by_value(r2, a1, b2, carry);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, q2, red_carry);

        let (r3, carry) = mac_with_carry_by_value(r3, a1, b3, carry);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, q3, red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(a1);
        drop(m);

        // round 2
        let a2 = (self.0).0[2];
        let (r0, carry) = mac_by_value(r0, a2, b0);
        let m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = mac_by_value_return_carry_only(r0, m, q0);

        // loop over the rest
        let (r1, carry) = mac_with_carry_by_value(r1, a2, b1, carry);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, q1, red_carry);

        let (r2, carry) = mac_with_carry_by_value(r2, a2, b2, carry);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, q2, red_carry);

        let (r3, carry) = mac_with_carry_by_value(r3, a2, b3, carry);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, q3, red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(a2);
        drop(m);

        // round 3
        let a3 = (self.0).0[3];
        let (r0, carry) = mac_by_value(r0, a3, b0);
        let m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = mac_by_value_return_carry_only(r0, m, q0);

        // loop over the rest
        let (r1, carry) = mac_with_carry_by_value(r1, a3, b1, carry);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, q1, red_carry);

        let (r2, carry) = mac_with_carry_by_value(r2, a3, b2, carry);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, q2, red_carry);

        let (r3, carry) = mac_with_carry_by_value(r3, a3, b3, carry);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, q3, red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(b3);
        drop(m);

        let mut result = Fs(FsRepr([r0, r1, r2, r3]));
        result.reduce();

        result
    }

    pub fn optimistic_cios_by_value_with_partial_red(self, other: Fs) -> Self {
        let mut m;

        let [q0, q1, q2, q3] = MODULUS.0;
        let [a0, a1, a2, a3] = (self.0).0;

        // round 0 
        let b0 = (other.0).0[0];
        let (r0, carry) = full_width_mul(a0, b0);
        m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = mac_by_value_return_carry_only(r0, m, q0);

        // loop over the rest
        let (r1, carry) = mac_by_value(carry, a1, b0);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, q1, red_carry);

        let (r2, carry) = mac_by_value(carry, a2, b0);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, q2, red_carry);

        let (r3, carry) = mac_by_value(carry, a3, b0);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, q3, red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(b0);

        // round 1
        let b1 = (other.0).0[1];
        let (r0, carry) = mac_by_value(r0, a0, b1);
        m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = mac_by_value_return_carry_only(r0, m, q0);

        // loop over the rest
        let (r1, carry) = mac_with_carry_by_value(r1, carry, a1, b1);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, q1, red_carry);

        let (r2, carry) = mac_with_carry_by_value(r2, carry, a2, b1);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, q2, red_carry);

        let (r3, carry) = mac_with_carry_by_value(r3, carry, a3, b1);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, q3, red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(b1);

        // round 2
        let b2 = (other.0).0[2];
        let (r0, carry) = mac_by_value(r0, a0, b2);
        m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = wrapping_mac_by_value(r0, m, q0);

        // loop over the rest
        let (r1, carry) = mac_with_carry_by_value(r1, carry, a1, b2);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, q1, red_carry);

        let (r2, carry) = mac_with_carry_by_value(r2, carry, a2, b2);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, q2, red_carry);

        let (r3, carry) = mac_with_carry_by_value(r3, carry, a3, b2);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, q3, red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(b2);

        // round 3
        let b3 = (other.0).0[3];
        let (r0, carry) = mac_by_value(r0, a0, b3);
        m = r0.wrapping_mul(INV);
        // everywhere semantic is arg0 + (arg1 * arg2)
        let red_carry = wrapping_mac_by_value(r0, m, q0);

        // loop over the rest
        let (r1, carry) = mac_with_carry_by_value(r1, carry, a1, b3);
        let (r0, red_carry) = mac_with_carry_by_value(r1, m, q1, red_carry);

        let (r2, carry) = mac_with_carry_by_value(r2, carry, a2, b3);
        let (r1, red_carry) = mac_with_carry_by_value(r2, m, q2, red_carry);

        let (r3, carry) = mac_with_carry_by_value(r3, carry, a3, b3);
        let (r2, red_carry) = mac_with_carry_by_value(r3, m, q3, red_carry);

        // this will check overflow in debug
        let r3 = red_carry + carry;
        drop(b3);

        Fs(FsRepr([r0, r1, r2, r3]))
    }

    pub fn mulx_latency_mul_assign(&mut self, other: &Fs) {
        let [b0, b1, b2, b3] = (other.0).0;
        let [a0, a1, a2, a3] = (self.0).0;

        // round 0 
        let (r1, r0) = full_width_mul(a0, b0);
        let (r2_0, r1_0) = full_width_mul(a1, b0);
        let (r2_1, r1_1) = full_width_mul(a0, b1);
        let (r3, r2) = full_width_mul(a1, b1);

        // now propagate carries
        let (r1, carry) = add(r1, r1_0);
        // hope for carry chains with add and addx
        let (r1, carry2) = add(r1, r1_1);
        // 
        let (r2, carry) = add_three(r2_0, r2, carry);
        let (r2_1, carry2) = add(r2_1, carry2);
        // 
        let (r2, carry2) = add(r2, r2_1); 
        let r3 = r3 + carry + carry2;

        // round 1
        let (r3_0, r2_0) = full_width_mul(a2, b0);
        let (r4_0, r3_0) = full_width_mul(a3, b0);
        let (r4_1, r3_1) = full_width_mul(a2, b1);
        let (r5, r4) = full_width_mul(a3, b1);

        let (r2, carry) = add(r2, r2_0);

        // now propagate carries
        let (r3, carry) = add_three(r3, r3_0, carry);
        // hope for carry chains with add and addx
        let (r3, carry2) = add(r3, r3_1);
        // 
        let (r4, carry) = add_three(r4_0, r4, carry);
        let (r4_1, carry2) = add(r4_1, carry2);
        // 
        let (r4, carry2) = add(r4, r4_1); 
        let (r5, r6) = add_three(r5, carry, carry2);

        // round 3
        let (r3_0, r2_0) = full_width_mul(a0, b2);
        let (r4_0, r3_0) = full_width_mul(a0, b3);
        let (r4_1, r3_1) = full_width_mul(a1, b2);
        let (r5_0, r4_0) = full_width_mul(a1, b3);

        let (r2, carry) = add(r2, r2_0);

        // now propagate carries
        let (r3, carry) = add_three(r3, r3_0, carry);
        // hope for carry chains with add and addx
        let (r3, carry2) = add(r3, r3_1);
        // 
        let (r4, carry) = add_three(r4_0, r4, carry);
        let (r4_1, carry2) = add(r4_1, carry2);
        // 
        let (r4, carry2) = add(r4, r4_1); 
        let (r5, carry) = add_four(r5, r5_0, carry, carry2);

        let r6 = r6 + carry;

        // round 4
        let (r5_0, r4_0) = full_width_mul(a2, b2);
        let (r6_0, r5_0) = full_width_mul(a2, b3);
        let (r6_1, r5_1) = full_width_mul(a3, b2);
        let (r7, r6_0) = full_width_mul(a3, b3);

        let (r4, carry) = add(r4, r2_0);

        // now propagate carries
        let (r5, carry) = add_three(r5, r5_0, carry);
        // hope for carry chains with add and addx
        let (r5, carry2) = add(r5, r5_1);
        // 
        let (r6, carry) = add_three(r6_0, r6, carry);
        let (r6_1, carry2) = add(r6_1, carry2);
        // 
        let (r6, carry2) = add(r6, r6_1); 
        let r7 = r7 + carry + carry2;

        self.mont_reduce(r0, r1, r2, r3, r4, r5, r6, r7);
    }

    pub fn asm_mul_assign(&mut self, other: &Fs) {
        *self = Fs(FsRepr(mont_mul_asm(&(self.0).0, &(other.0).0)));
    }
}

impl SqrtField for Fs {

    fn legendre(&self) -> LegendreSymbol {
        // s = self^((s - 1) // 2)
        let s = self.pow([0x684b872f6b7b965b, 0x53341049e6640841, 0x83339d80809a1d80, 0x73eda753299d7d4]);
        if s == Self::zero() { Zero }
        else { QuadraticNonResidue }
    }

    fn sqrt(&self) -> Option<Self> {
        // Shank's algorithm for s mod 4 = 3
        // https://eprint.iacr.org/2012/685.pdf (page 9, algorithm 2)

        // a1 = self^((s - 3) // 4)
        let mut a1 = self.pow([0xb425c397b5bdcb2d, 0x299a0824f3320420, 0x4199cec0404d0ec0, 0x39f6d3a994cebea]);
        let mut a0 = a1;
        a0.square();
        a0.mul_assign(self);

        if a0 == NEGATIVE_ONE
        {
            None
        }
        else
        {
            a1.mul_assign(self);
            Some(a1)
        }
    }
}

#[inline(always)]
pub fn full_width_mul(a: u64, b: u64) -> (u64, u64) {
    let tmp = (a as u128) * (b as u128);

    return (tmp as u64, (tmp >> 64) as u64);
}

#[inline(always)]
pub fn mac_with_carry_by_value(a: u64, b: u64, c: u64, carry: u64) -> (u64, u64) {
    let tmp = ((b as u128) * (c as u128)) + (a as u128) + (carry as u128);

    (tmp as u64, (tmp >> 64) as u64)
}

#[inline(always)]
pub fn add(a: u64, b: u64) -> (u64, u64) {
    let tmp = (a as u128) + (b as u128);

    (tmp as u64, (tmp >> 64) as u64)
}

#[inline(always)]
pub fn add_three(a: u64, b: u64, c: u64) -> (u64, u64) {
    let tmp = (a as u128) + (b as u128) + (c as u128);
    
    (tmp as u64, (tmp >> 64) as u64)
}

#[inline(always)]
pub fn add_four(a: u64, b: u64, c: u64, d: u64) -> (u64, u64) {
    let tmp = (a as u128) + (b as u128) + (c as u128) + (d as u128);
    
    (tmp as u64, (tmp >> 64) as u64)
}

#[inline(always)]
pub fn wrapping_mac_by_value(a: u64, b: u64, c: u64) -> u64 {
    b.wrapping_mul(c).wrapping_add(a)
}

#[inline(always)]
pub fn mac_by_value_return_carry_only(a: u64, b: u64, c: u64) -> u64 {
    let tmp = ((b as u128) * (c as u128)) + (a as u128);

    (tmp >> 64) as u64
}

#[inline(always)]
pub fn mac_by_value(a: u64, b: u64, c: u64) -> (u64, u64) {
    let tmp = ((b as u128) * (c as u128)) + (a as u128);

    (tmp as u64, (tmp >> 64) as u64)
}

#[inline(always)]
pub fn mul_with_high_carry_by_value(a: u64, b: u64, c: u64, mut carry: u64) -> (u64, u64, u64) {
    let (hi, lo) = full_width_mul(b, c);
    let (hi, of) = hi.overflowing_add(a);
    if of {
        carry += 1;
    }

    (hi, lo, carry)
}

#[inline(always)]
pub fn mul_and_add_existing_high(a: u64, b: u64, existing_hi: u64, carry: u64) -> (u64, u64, u64) {
    let tmp = (a as u128) * (b as u128);

    let hi = tmp >> 64;
    let lo = tmp as u64;
    
    let tmp = hi + (existing_hi as u128) + ((carry as u128) << 64);

    let carry = (tmp >> 64) as u64; 
    let hi = tmp as u64;

    (lo, hi, carry)
}

#[inline(always)]
pub fn mul_and_add_existing(a: u64, b: u64, existing_lo: u64, existing_hi: u64, carry: u64) -> (u64, u64, u64) {
    let tmp = ((a as u128) * (b as u128)) + (existing_lo as u128);

    let hi = tmp >> 64;
    let lo = tmp as u64;
    
    let tmp = hi + (existing_hi as u128) + ((carry as u128) << 64);

    let carry = (tmp >> 64) as u64; 
    let hi = tmp as u64;

    (lo, hi, carry)
}

#[inline(always)]
pub fn mul_and_full_block_propagate(a: u64, b: u64, z0: u64, z1: u64, z2: u64, z3: u64, t0: u64, t1: u64) -> (u64, u64, u64, u64, u64) {
    let tmp = ((a as u128) * (b as u128)) + (z0 as u128);

    let hi = tmp >> 64;
    let z0 = tmp as u64;
    
    // let tmp = hi + (z1 as u128);

    let tmp = hi + (z1 as u128) + ((t0 as u128) << 64);

    let hi = tmp >> 64; 
    let z1 = tmp as u64;

    let tmp = hi + (z2 as u128) + ((t1 as u128) << 64);

    let hi = tmp >> 64; 
    let z2 = tmp as u64;

    let tmp = hi + (z3 as u128);

    let c0 = (tmp >> 64) as u64; 
    let z3 = tmp as u64;

    (z0, z1, z2, z3, c0)
}

#[inline(always)]
pub fn mul_and_full_block_propagate_into_existing_carry_catcher(a: u64, b: u64, z0: u64, z1: u64, z2: u64, z3: u64, t0: u64, t1: u64, c: u64) -> (u64, u64, u64, u64, u64) {
    let tmp = ((a as u128) * (b as u128)) + (z0 as u128);

    let hi = tmp >> 64;
    let z0 = tmp as u64;
    
    // let tmp = hi + (z1 as u128);

    let tmp = hi + (z1 as u128) + ((t0 as u128) << 64);

    let hi = tmp >> 64; 
    let z1 = tmp as u64;

    let tmp = hi + (z2 as u128) + ((t1 as u128) << 64);

    let hi = tmp >> 64; 
    let z2 = tmp as u64;

    let tmp = hi + (z3 as u128) + ((c as u128) << 64);

    let c0 = (tmp >> 64) as u64; 
    let z3 = tmp as u64;

    (z0, z1, z2, z3, c0)
}

#[inline(always)]
pub fn mul_and_full_block_propagate_without_carry_catch(a: u64, b: u64, z0: u64, z1: u64, z2: u64, z3: u64, t0: u64, t1: u64) -> (u64, u64, u64, u64) {
    let tmp = ((a as u128) * (b as u128)) + (z0 as u128);

    let hi = tmp >> 64;
    let z0 = tmp as u64;
    
    // let tmp = hi + (z1 as u128);

    let tmp = hi + (z1 as u128) + ((t0 as u128) << 64);

    let hi = tmp >> 64; 
    let z1 = tmp as u64;

    let tmp = hi + (z2 as u128) + ((t1 as u128) << 64);

    let hi = (tmp >> 64) as u64; 
    let z2 = tmp as u64;

    let z3 = hi + z3;

    (z0, z1, z2, z3)
}


// // Computes one "row" of multiplication: [a0, a1, a2, a3] * b0
// // Uses MULX
// #[allow(dead_code)]
// #[inline(always)]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
// pub(crate) fn mul_1_asm(a: u64, b0: u64, b1: u64, b2: u64, b3: u64) -> (u64, u64, u64, u64, u64) {
//     let r0: u64;
//     let r1: u64;
//     let r2: u64;
//     let r3: u64;
//     let r4: u64;
//     let _lo: u64;
//     // Binding `_lo` will not be used after assignment.
//     #[allow(clippy::used_underscore_binding)]
//     unsafe {
//         asm!(r"
//         mulx $7, $0, $1      // (r0, r1) = a * b0
//         mulx $8, $5, $2      // (lo, r2) = a * b1
//         add $5, $1           // r1 += lo (carry in CF)
//         mulx $9, $5, $3      // (lo, r3) = a * b2
//         adc $5, $2           // r2 += lo + CF (carry in CF)
//         mulx $10, $5, $4     // (lo, r4) = a * b3
//         adc $5, $3           // r3 += lo + CF (carry in CF)
//         adc $11, $4          // r4 += 0 + CF (no carry, CF to 0)
//         "
//         : // Output constraints
//             "=&r"(r0),   // $0 r0..4 are in registers
//             "=&r"(r1),   // $1
//             "=&r"(r2),   // $2
//             "=&r"(r3),   // $3
//             "=&r"(r4)    // $4
//             "=&r"(_lo)   // $5 Temporary values can be in any register
//         : // Input constraints
//             "{rdx}"(a), // $6 a must be in RDX for MULX to work
//             "rm"(b0),   // $7 b0..b3 can be register or memory
//             "rm"(b1),   // $8
//             "rm"(b2),   // $9
//             "rm"(b3),   // $10
//             "i"(0)      // $11 Immediate zero
//         : // Clobbers
//            "cc"         // Flags
//         )
//     }
//     (r0, r1, r2, r3, r4)
// }

// // Computes r[0..4] += a * b[0..4], returns carry
// // Uses MULX and ADCX/ADOX carry chain
// // Currently unused
// #[allow(dead_code)]
// #[allow(clippy::too_many_arguments)]
// #[inline(always)]
// #[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
// pub(crate) fn mul_add_1_asm(
//     r0: &mut u64,
//     r1: &mut u64,
//     r2: &mut u64,
//     r3: &mut u64,
//     a: u64,
//     b0: u64,
//     b1: u64,
//     b2: u64,
//     b3: u64,
// ) -> u64 {
//     let _lo: u64;
//     let _hi: u64;
//     let r4: u64;
//     // Bindings `_lo` and `_hi` will not be used after assignment.
//     #[allow(clippy::used_underscore_binding)]
//     unsafe {
//         asm!(r"
//         xor $4, $4            // r4 = CF = OF 0
//         mulx $8, $5, $6       // a * b0
//         adcx $5, $0           // r0 += lo + CF (carry in CF)
//         adox $6, $1           // r1 += hi + OF (carry in OF)
//         mulx $9, $5, $6       // a * b1
//         adcx $5, $1           // r1 += lo + CF (carry in CF)
//         adox $6, $2           // r2 += hi + OF (carry in OF)
//         mulx $10, $5, $6      // a * b2
//         adcx $5, $2           // r2 += lo + CF (carry in CF)
//         adox $6, $3           // r3 += hi + OF (carry in OF)
//         mulx $11, $5, $6      // a * b3
//         adcx $5, $3           // r3 += lo + CF (carry in CF)
//         adcx $4, $4           // r4 += CF (no carry, CF = 0)
//         adox $6, $4           // r4 += hi + OF (no carry, OF = 0)
//         "
//         : // Output constraints
//             "+r"(*r0),   // $0 r0..3 are in register and modified in place
//             "+r"(*r1),   // $1
//             "+r"(*r2),   // $2
//             "+r"(*r3),   // $3
//             "=&r"(r4)    // $4 r4 is output to a register
//             "=&r"(_lo),  // $5 Temporary values can be in any register
//             "=&r"(_hi)   // $6
//         : // Input constraints
//             "{rdx}"(a), // $7 a must be in RDX for MULX to work
//             "rm"(b0),   // $8 Second operand can be register or memory
//             "rm"(b1),   // $9 Second operand can be register or memory
//             "rm"(b2),   // $10 Second operand can be register or memory
//             "rm"(b3)    // $11 Second operand can be register or memory
//         : // Clobbers
//            "cc"         // Flags
//         )
//     }
//     r4
// }

// Currently unused
#[allow(dead_code)]
#[allow(clippy::too_many_lines)]
#[inline(always)]
#[cfg(all(target_arch = "x86_64", target_feature = "adx"))]
pub fn mont_mul_asm(a: &[u64; 4], b: &[u64; 4]) -> [u64; 4] {
    use core::mem::MaybeUninit;

    const ZERO: u64 = 0; // $3

    let mut result = MaybeUninit::<[u64; 4]>::uninit();
    // MULX dst_high, dst_low, src_b (src_a = %rdx)
    // src_b can be register or memory, not immediate
    unsafe {
        llvm_asm!(r"
            // Assembly from Aztec's Barretenberg implementation, see 
            // <https://github.com/AztecProtocol/barretenberg/blob/master/src/barretenberg/fields/asm_macros.hpp>
            movq 0($1), %rdx
            xorq %r8, %r8
            mulxq 8($2), %r8, %r9
            mulxq 24($2), %rdi, %r12
            mulxq 0($2), %r13, %r14
            mulxq 16($2), %r15, %r10
            movq %r13, %rdx
            mulxq $8, %rdx, %r11
            adcxq %r8, %r14
            adoxq %rdi, %r10
            adcxq %r9, %r15
            adoxq $3, %r12
            adcxq $3, %r10
            mulxq $4, %r8, %r9
            mulxq $5, %rdi, %r11
            adoxq %r8, %r13
            adcxq %rdi, %r14
            adoxq %r9, %r14
            adcxq %r11, %r15
            mulxq $6, %r8, %r9
            mulxq $7, %rdi, %r11
            adoxq %r8, %r15
            adcxq %rdi, %r10
            adoxq %r9, %r10
            adcxq %r11, %r12
            adoxq $3, %r12
            movq 8($1), %rdx
            mulxq 0($2), %r8, %r9
            mulxq 8($2), %rdi, %r11
            adcxq %r8, %r14
            adoxq %r9, %r15
            adcxq %rdi, %r15
            adoxq %r11, %r10
            mulxq 16($2), %r8, %r9
            mulxq 24($2), %rdi, %r13
            adcxq %r8, %r10
            adoxq %rdi, %r12
            adcxq %r9, %r12
            adoxq $3, %r13
            adcxq $3, %r13
            movq %r14, %rdx
            mulxq $8, %rdx, %r8
            mulxq $4, %r8, %r9
            mulxq $5, %rdi, %r11
            adoxq %r8, %r14
            adcxq %rdi, %r15
            adoxq %r9, %r15
            adcxq %r11, %r10
            mulxq $6, %r8, %r9
            mulxq $7, %rdi, %r11
            adoxq %r8, %r10
            adcxq %r9, %r12
            adoxq %rdi, %r12
            adcxq %r11, %r13
            adoxq $3, %r13
            movq 16($1), %rdx
            mulxq 0($2), %r8, %r9
            mulxq 8($2), %rdi, %r11
            adcxq %r8, %r15
            adoxq %r9, %r10
            adcxq %rdi, %r10
            adoxq %r11, %r12
            mulxq 16($2), %r8, %r9
            mulxq 24($2), %rdi, %r14
            adcxq %r8, %r12
            adoxq %r9, %r13
            adcxq %rdi, %r13
            adoxq $3, %r14
            adcxq $3, %r14
            movq %r15, %rdx
            mulxq $8, %rdx, %r8
            mulxq $4, %r8, %r9
            mulxq $5, %rdi, %r11
            adoxq %r8, %r15
            adcxq %r9, %r10
            adoxq %rdi, %r10
            adcxq %r11, %r12
            mulxq $6, %r8, %r9
            mulxq $7, %rdi, %r11
            adoxq %r8, %r12
            adcxq %r9, %r13
            adoxq %rdi, %r13
            adcxq %r11, %r14
            adoxq $3, %r14
            movq 24($1), %rdx
            mulxq 0($2), %r8, %r9
            mulxq 8($2), %rdi, %r11
            adcxq %r8, %r10
            adoxq %r9, %r12
            adcxq %rdi, %r12
            adoxq %r11, %r13
            mulxq 16($2), %r8, %r9
            mulxq 24($2), %rdi, %r15
            adcxq %r8, %r13
            adoxq %r9, %r14
            adcxq %rdi, %r14
            adoxq $3, %r15
            adcxq $3, %r15
            movq %r10, %rdx
            mulxq $8, %rdx, %r8
            mulxq $4, %r8, %r9
            mulxq $5, %rdi, %r11
            adoxq %r8, %r10
            adcxq %r9, %r12
            adoxq %rdi, %r12
            adcxq %r11, %r13
            mulxq $6, %r8, %r9
            mulxq $7, %rdi, %rdx
            adoxq %r8, %r13
            adcxq %r9, %r14
            adoxq %rdi, %r14
            adcxq %rdx, %r15
            adoxq $3, %r15
            movq %r12, 0($0)
            movq %r13, 8($0)
            movq %r14, 16($0)
            movq %r15, 24($0)
            "
            :
            : "r"(result.as_mut_ptr()),
              "r"(a), "r"(b),
              "m"(ZERO),
              "m"(MODULUS.0[0]),
              "m"(MODULUS.0[1]),
              "m"(MODULUS.0[2]),
              "m"(MODULUS.0[3]),
              "m"(INV)
            : "rdx", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "cc", "memory"
        );
    }
    let result = unsafe { result.assume_init() };

    result
}

#[cfg(test)]
mod test {
    use super::Fs;
    use ff::{Field, PrimeField};

    use rand::{*};

    #[test]
    fn test_optimistic_cios() {
        let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        for _ in 0..10000 {
            let a: Fs = rng.gen(); 
            let b: Fs = rng.gen();

            let mut c = a;
            c.mul_assign(&b);

            let d = a.optimistic_cios_by_value(b);

            assert_eq!(c, d);
        }
    }
}