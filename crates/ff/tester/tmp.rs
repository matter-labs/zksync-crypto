mod test_large_cios_field {
    use ff::*;
    #[PrimeFieldModulus = "21888242871839275222246405745257275088696311157297823662689037894645226208583"]
    #[PrimeFieldGenerator = "2"]
    #[OptimisticCIOSMultiplication = "true"]
    #[OptimisticCIOSSquaring = "true"]
    pub(crate) struct Fr(FrRepr);
    /// This is the modulus m of the prime field
    const MODULUS: FrRepr = FrRepr([
        4332616871279656263u64,
        10917124144477883021u64,
        13281191951274694749u64,
        3486998266802970665u64,
    ]);
    /// The number of bits needed to represent the modulus.
    const MODULUS_BITS: u32 = 254u32;
    /// The number of bits that must be shaved from the beginning of
    /// the representation when randomly sampling.
    const REPR_SHAVE_BITS: u32 = 2u32;
    /// 2^{limbs*64} mod m
    const R: FrRepr = FrRepr([
        15230403791020821917u64,
        754611498739239741u64,
        7381016538464732716u64,
        1011752739694698287u64,
    ]);
    /// 2^{limbs*64*2} mod m
    const R2: FrRepr = FrRepr([
        17522657719365597833u64,
        13107472804851548667u64,
        5164255478447964150u64,
        493319470278259999u64,
    ]);
    /// -(m^{-1} mod m) mod m
    const INV: u64 = 9786893198990664585u64;
    /// Multiplicative generator of `MODULUS` - 1 order, also quadratic
    /// nonresidue.
    const GENERATOR: FrRepr = FrRepr([
        12014063508332092218u64,
        1509222997478479483u64,
        14762033076929465432u64,
        2023505479389396574u64,
    ]);
    /// 2^s * t = MODULUS - 1 with t odd
    const S: u32 = 1u32;
    /// 2^s root of unity computed by GENERATOR^t
    const ROOT_OF_UNITY: FrRepr = FrRepr([
        15230403791020821917u64,
        754611498739239741u64,
        7381016538464732716u64,
        1011752739694698287u64,
    ]);
    pub struct FrRepr(pub [u64; 4usize]);
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::marker::Copy for FrRepr {}
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::clone::Clone for FrRepr {
        #[inline]
        fn clone(&self) -> FrRepr {
            {
                let _: ::core::clone::AssertParamIsClone<[u64; 4usize]>;
                *self
            }
        }
    }
    impl ::core::marker::StructuralPartialEq for FrRepr {}
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::cmp::PartialEq for FrRepr {
        #[inline]
        fn eq(&self, other: &FrRepr) -> bool {
            match *other {
                FrRepr(ref __self_1_0) => match *self {
                    FrRepr(ref __self_0_0) => (*__self_0_0) == (*__self_1_0),
                },
            }
        }
        #[inline]
        fn ne(&self, other: &FrRepr) -> bool {
            match *other {
                FrRepr(ref __self_1_0) => match *self {
                    FrRepr(ref __self_0_0) => (*__self_0_0) != (*__self_1_0),
                },
            }
        }
    }
    impl ::core::marker::StructuralEq for FrRepr {}
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::cmp::Eq for FrRepr {
        #[inline]
        #[doc(hidden)]
        fn assert_receiver_is_total_eq(&self) -> () {
            {
                let _: ::core::cmp::AssertParamIsEq<[u64; 4usize]>;
            }
        }
    }
    #[automatically_derived]
    #[allow(unused_qualifications)]
    impl ::core::default::Default for FrRepr {
        #[inline]
        fn default() -> FrRepr {
            FrRepr(::core::default::Default::default())
        }
    }
    impl ::std::fmt::Debug for FrRepr {
        fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
            f.write_fmt(::core::fmt::Arguments::new_v1(
                &["0x"],
                &match () {
                    () => [],
                },
            ))?;
            for i in self.0.iter().rev() {
                f.write_fmt(::core::fmt::Arguments::new_v1_formatted(
                    &[""],
                    &match (&*i,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::LowerHex::fmt,
                        )],
                    },
                    &[::core::fmt::rt::v1::Argument {
                        position: 0usize,
                        format: ::core::fmt::rt::v1::FormatSpec {
                            fill: ' ',
                            align: ::core::fmt::rt::v1::Alignment::Unknown,
                            flags: 8u32,
                            precision: ::core::fmt::rt::v1::Count::Implied,
                            width: ::core::fmt::rt::v1::Count::Is(16usize),
                        },
                    }],
                ))?;
            }
            Ok(())
        }
    }
    impl ::rand::Rand for FrRepr {
        #[inline(always)]
        fn rand<R: ::rand::Rng>(rng: &mut R) -> Self {
            FrRepr(rng.gen())
        }
    }
    impl ::std::fmt::Display for FrRepr {
        fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
            f.write_fmt(::core::fmt::Arguments::new_v1(
                &["0x"],
                &match () {
                    () => [],
                },
            ))?;
            for i in self.0.iter().rev() {
                f.write_fmt(::core::fmt::Arguments::new_v1_formatted(
                    &[""],
                    &match (&*i,) {
                        (arg0,) => [::core::fmt::ArgumentV1::new(
                            arg0,
                            ::core::fmt::LowerHex::fmt,
                        )],
                    },
                    &[::core::fmt::rt::v1::Argument {
                        position: 0usize,
                        format: ::core::fmt::rt::v1::FormatSpec {
                            fill: ' ',
                            align: ::core::fmt::rt::v1::Alignment::Unknown,
                            flags: 8u32,
                            precision: ::core::fmt::rt::v1::Count::Implied,
                            width: ::core::fmt::rt::v1::Count::Is(16usize),
                        },
                    }],
                ))?;
            }
            Ok(())
        }
    }
    impl std::hash::Hash for FrRepr {
        fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
            for limb in self.0.iter() {
                limb.hash(state);
            }
        }
    }
    impl AsRef<[u64]> for FrRepr {
        #[inline(always)]
        fn as_ref(&self) -> &[u64] {
            &self.0
        }
    }
    impl AsMut<[u64]> for FrRepr {
        #[inline(always)]
        fn as_mut(&mut self) -> &mut [u64] {
            &mut self.0
        }
    }
    impl From<u64> for FrRepr {
        #[inline(always)]
        fn from(val: u64) -> FrRepr {
            use std::default::Default;
            let mut repr = Self::default();
            repr.0[0] = val;
            repr
        }
    }
    impl Ord for FrRepr {
        #[inline(always)]
        fn cmp(&self, other: &FrRepr) -> ::std::cmp::Ordering {
            for (a, b) in self.0.iter().rev().zip(other.0.iter().rev()) {
                if a < b {
                    return ::std::cmp::Ordering::Less;
                } else if a > b {
                    return ::std::cmp::Ordering::Greater;
                }
            }
            ::std::cmp::Ordering::Equal
        }
    }
    impl PartialOrd for FrRepr {
        #[inline(always)]
        fn partial_cmp(&self, other: &FrRepr) -> Option<::std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }
    impl crate::ff::PrimeFieldRepr for FrRepr {
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
            if n as usize >= 64 * 4usize {
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
            if n as usize >= 64 * 4usize {
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
            let mut ret = (4usize as u32) * 64;
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
        fn add_nocarry(&mut self, other: &FrRepr) {
            let mut carry = 0;
            for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
                *a = crate::ff::adc(*a, *b, &mut carry);
            }
        }
        #[inline(always)]
        fn sub_noborrow(&mut self, other: &FrRepr) {
            let mut borrow = 0;
            for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
                *a = crate::ff::sbb(*a, *b, &mut borrow);
            }
        }
    }
    impl ::std::marker::Copy for Fr {}
    impl ::std::clone::Clone for Fr {
        fn clone(&self) -> Fr {
            *self
        }
    }
    impl ::std::cmp::PartialEq for Fr {
        fn eq(&self, other: &Fr) -> bool {
            self.0 == other.0
        }
    }
    impl ::std::cmp::Eq for Fr {}
    impl ::std::fmt::Debug for Fr {
        fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
            f.write_fmt(::core::fmt::Arguments::new_v1(
                &["", "(", ")"],
                &match (&"Fr", &self.into_repr()) {
                    (arg0, arg1) => [
                        ::core::fmt::ArgumentV1::new(arg0, ::core::fmt::Display::fmt),
                        ::core::fmt::ArgumentV1::new(arg1, ::core::fmt::Debug::fmt),
                    ],
                },
            ))
        }
    }
    /// Elements are ordered lexicographically.
    impl Ord for Fr {
        #[inline(always)]
        fn cmp(&self, other: &Fr) -> ::std::cmp::Ordering {
            self.into_repr().cmp(&other.into_repr())
        }
    }
    impl PartialOrd for Fr {
        #[inline(always)]
        fn partial_cmp(&self, other: &Fr) -> Option<::std::cmp::Ordering> {
            Some(self.cmp(other))
        }
    }
    impl ::std::fmt::Display for Fr {
        fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
            f.write_fmt(::core::fmt::Arguments::new_v1(
                &["", "(", ")"],
                &match (&"Fr", &self.into_repr()) {
                    (arg0, arg1) => [
                        ::core::fmt::ArgumentV1::new(arg0, ::core::fmt::Display::fmt),
                        ::core::fmt::ArgumentV1::new(arg1, ::core::fmt::Display::fmt),
                    ],
                },
            ))
        }
    }
    impl ::rand::Rand for Fr {
        /// Computes a uniformly random element using rejection sampling.
        fn rand<R: ::rand::Rng>(rng: &mut R) -> Self {
            loop {
                let mut tmp = Fr(FrRepr::rand(rng));
                tmp.0.as_mut()[3usize] &= 0xffffffffffffffff >> REPR_SHAVE_BITS;
                if tmp.is_valid() {
                    return tmp;
                }
            }
        }
    }
    impl From<Fr> for FrRepr {
        fn from(e: Fr) -> FrRepr {
            e.into_repr()
        }
    }
    impl crate::ff::PrimeField for Fr {
        type Repr = FrRepr;
        fn from_repr(r: FrRepr) -> Result<Fr, crate::ff::PrimeFieldDecodingError> {
            let mut r = Fr(r);
            if r.is_valid() {
                r.mul_assign(&Fr(R2));
                Ok(r)
            } else {
                Err(crate::ff::PrimeFieldDecodingError::NotInField({
                    let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                        &[""],
                        &match (&r.0,) {
                            (arg0,) => [::core::fmt::ArgumentV1::new(
                                arg0,
                                ::core::fmt::Display::fmt,
                            )],
                        },
                    ));
                    res
                }))
            }
        }
        fn from_raw_repr(r: FrRepr) -> Result<Self, crate::ff::PrimeFieldDecodingError> {
            let mut r = Fr(r);
            if r.is_valid() {
                Ok(r)
            } else {
                Err(crate::ff::PrimeFieldDecodingError::NotInField({
                    let res = ::alloc::fmt::format(::core::fmt::Arguments::new_v1(
                        &[""],
                        &match (&r.0,) {
                            (arg0,) => [::core::fmt::ArgumentV1::new(
                                arg0,
                                ::core::fmt::Display::fmt,
                            )],
                        },
                    ));
                    res
                }))
            }
        }
        fn into_repr(&self) -> FrRepr {
            let mut r = *self;
            r.mont_reduce(
                (self.0).0[0usize],
                (self.0).0[1usize],
                (self.0).0[2usize],
                (self.0).0[3usize],
                0,
                0,
                0,
                0,
            );
            r.0
        }
        fn into_raw_repr(&self) -> FrRepr {
            let r = *self;
            r.0
        }
        fn char() -> FrRepr {
            MODULUS
        }
        const NUM_BITS: u32 = MODULUS_BITS;
        const CAPACITY: u32 = Self::NUM_BITS - 1;
        fn multiplicative_generator() -> Self {
            Fr(GENERATOR)
        }
        const S: u32 = S;
        fn root_of_unity() -> Self {
            Fr(ROOT_OF_UNITY)
        }
    }
    impl crate::ff::Field for Fr {
        #[inline]
        fn zero() -> Self {
            Fr(FrRepr::from(0))
        }
        #[inline]
        fn one() -> Self {
            Fr(R)
        }
        #[inline]
        fn is_zero(&self) -> bool {
            self.0.is_zero()
        }
        #[inline]
        fn add_assign(&mut self, other: &Fr) {
            self.0.add_nocarry(&other.0);
            self.reduce();
        }
        #[inline]
        fn double(&mut self) {
            self.0.mul2();
            self.reduce();
        }
        #[inline]
        fn sub_assign(&mut self, other: &Fr) {
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
                let one = FrRepr::from(1);
                let mut u = self.0;
                let mut v = MODULUS;
                let mut b = Fr(R2);
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
        fn frobenius_map(&mut self, _: usize) {}
        #[inline]
        fn mul_assign(&mut self, other: &Fr) {
            let [b0, b1, b2, b3] = (other.0).0;
            let a = (self.0).0[0usize];
            let (r0, carry) = crate::ff::full_width_mul(a, b0);
            let m = r0.wrapping_mul(INV);
            let red_carry = crate::ff::mac_by_value_return_carry_only(r0, m, MODULUS.0[0]);
            let (r1, carry) = crate::ff::mac_by_value(carry, a, b1);
            let (r0, red_carry) =
                crate::ff::mac_with_carry_by_value(r1, m, MODULUS.0[1usize], red_carry);
            let (r2, carry) = crate::ff::mac_by_value(carry, a, b2);
            let (r1, red_carry) =
                crate::ff::mac_with_carry_by_value(r2, m, MODULUS.0[2usize], red_carry);
            let (r3, carry) = crate::ff::mac_by_value(carry, a, b3);
            let (r2, red_carry) =
                crate::ff::mac_with_carry_by_value(r3, m, MODULUS.0[3usize], red_carry);
            let r3 = red_carry + carry;
            let a = (self.0).0[1usize];
            let (r0, carry) = crate::ff::mac_by_value(r0, a, b0);
            let m = r0.wrapping_mul(INV);
            let red_carry = crate::ff::mac_by_value_return_carry_only(r0, m, MODULUS.0[0]);
            let (r1, carry) = crate::ff::mac_with_carry_by_value(r1, a, b1, carry);
            let (r0, red_carry) =
                crate::ff::mac_with_carry_by_value(r1, m, MODULUS.0[1usize], red_carry);
            let (r2, carry) = crate::ff::mac_with_carry_by_value(r2, a, b2, carry);
            let (r1, red_carry) =
                crate::ff::mac_with_carry_by_value(r2, m, MODULUS.0[2usize], red_carry);
            let (r3, carry) = crate::ff::mac_with_carry_by_value(r3, a, b3, carry);
            let (r2, red_carry) =
                crate::ff::mac_with_carry_by_value(r3, m, MODULUS.0[3usize], red_carry);
            let r3 = red_carry + carry;
            let a = (self.0).0[2usize];
            let (r0, carry) = crate::ff::mac_by_value(r0, a, b0);
            let m = r0.wrapping_mul(INV);
            let red_carry = crate::ff::mac_by_value_return_carry_only(r0, m, MODULUS.0[0]);
            let (r1, carry) = crate::ff::mac_with_carry_by_value(r1, a, b1, carry);
            let (r0, red_carry) =
                crate::ff::mac_with_carry_by_value(r1, m, MODULUS.0[1usize], red_carry);
            let (r2, carry) = crate::ff::mac_with_carry_by_value(r2, a, b2, carry);
            let (r1, red_carry) =
                crate::ff::mac_with_carry_by_value(r2, m, MODULUS.0[2usize], red_carry);
            let (r3, carry) = crate::ff::mac_with_carry_by_value(r3, a, b3, carry);
            let (r2, red_carry) =
                crate::ff::mac_with_carry_by_value(r3, m, MODULUS.0[3usize], red_carry);
            let r3 = red_carry + carry;
            let a = (self.0).0[3usize];
            let (r0, carry) = crate::ff::mac_by_value(r0, a, b0);
            let m = r0.wrapping_mul(INV);
            let red_carry = crate::ff::mac_by_value_return_carry_only(r0, m, MODULUS.0[0]);
            let (r1, carry) = crate::ff::mac_with_carry_by_value(r1, a, b1, carry);
            let (r0, red_carry) =
                crate::ff::mac_with_carry_by_value(r1, m, MODULUS.0[1usize], red_carry);
            let (r2, carry) = crate::ff::mac_with_carry_by_value(r2, a, b2, carry);
            let (r1, red_carry) =
                crate::ff::mac_with_carry_by_value(r2, m, MODULUS.0[2usize], red_carry);
            let (r3, carry) = crate::ff::mac_with_carry_by_value(r3, a, b3, carry);
            let (r2, red_carry) =
                crate::ff::mac_with_carry_by_value(r3, m, MODULUS.0[3usize], red_carry);
            let r3 = red_carry + carry;
            *self = Fr(FrRepr([r0, r1, r2, r3]));
            self.reduce();
        }
        #[inline]
        fn square(&mut self) {
            let [a0, a1, a2, a3] = (self.0).0;
            let (r0, carry) = crate::ff::full_width_mul(a0, a0);
            let m = r0.wrapping_mul(INV);
            let red_carry = crate::ff::mac_by_value_return_carry_only(r0, m, MODULUS.0[0]);
            let (r1, carry, superhi) = crate::ff::mul_double_add_by_value(carry, a0, a1);
            let (r0, red_carry) =
                crate::ff::mac_with_carry_by_value(r1, m, MODULUS.0[1usize], red_carry);
            let (r2, carry, superhi) =
                crate::ff::mul_double_add_low_and_high_carry_by_value(a0, a2, carry, superhi);
            let (r1, red_carry) =
                crate::ff::mac_with_carry_by_value(r2, m, MODULUS.0[2usize], red_carry);
            let (r3, carry) = crate::ff::mul_double_add_low_and_high_carry_by_value_ignore_superhi(
                a0, a3, carry, superhi,
            );
            let (r2, r3) = crate::ff::mac_with_low_and_high_carry_by_value(
                red_carry,
                m,
                MODULUS.0[3usize],
                r3,
                carry,
            );
            let m = r0.wrapping_mul(INV);
            let red_carry = crate::ff::mac_by_value_return_carry_only(r0, m, MODULUS.0[0]);
            let (r1, carry) = crate::ff::mac_by_value(r1, a1, a1);
            let (r0, red_carry) =
                crate::ff::mac_with_carry_by_value(r1, m, MODULUS.0[1usize], red_carry);
            let (r2, carry, superhi) =
                crate::ff::mul_double_add_add_carry_by_value(r2, a1, a2, carry);
            let (r1, red_carry) = mac_with_carry_by_value(r2, m, MODULUS.0[2usize], red_carry);
            let (r3, carry) =
                crate::ff::mul_double_add_add_low_and_high_carry_by_value_ignore_superhi(
                    r3, a1, a3, carry, superhi,
                );
            let (r2, r3) = crate::ff::mac_with_low_and_high_carry_by_value(
                red_carry,
                m,
                MODULUS.0[3usize],
                r3,
                carry,
            );
            let m = r0.wrapping_mul(INV);
            let red_carry = crate::ff::mac_by_value_return_carry_only(r0, m, MODULUS.0[0]);
            let (r0, red_carry) =
                crate::ff::mac_with_carry_by_value(r1, m, MODULUS.0[1usize], red_carry);
            let (r2, carry) = crate::ff::mac_by_value(r2, a2, a2);
            let (r1, red_carry) =
                crate::ff::mac_with_carry_by_value(r2, m, MODULUS.0[2usize], red_carry);
            let (r3, carry) =
                crate::ff::mul_double_add_add_carry_by_value_ignore_superhi(r3, a2, a3, carry);
            let (r2, r3) = crate::ff::mac_with_low_and_high_carry_by_value(
                red_carry,
                m,
                MODULUS.0[3usize],
                r3,
                carry,
            );
            let m = r0.wrapping_mul(INV);
            let red_carry = crate::ff::mac_by_value_return_carry_only(r0, m, MODULUS.0[0]);
            let (r0, red_carry) =
                crate::ff::mac_with_carry_by_value(r1, m, MODULUS.0[1usize], red_carry);
            let (r1, red_carry) =
                crate::ff::mac_with_carry_by_value(r2, m, MODULUS.0[2usize], red_carry);
            let (r3, carry) = crate::ff::mac_by_value(r3, a3, a3);
            let (r2, r3) = crate::ff::mac_with_low_and_high_carry_by_value(
                red_carry,
                m,
                MODULUS.0[3usize],
                r3,
                carry,
            );
            *self = Fr(FrRepr([r0, r1, r2, r3]));
            self.reduce();
        }
    }
    impl std::default::Default for Fr {
        fn default() -> Self {
            Self::zero()
        }
    }
    impl std::hash::Hash for Fr {
        fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
            for limb in self.0.as_ref().iter() {
                limb.hash(state);
            }
        }
    }
    impl Fr {
        /// Determines if the element is really in the field. This is only used
        /// internally.
        #[inline(always)]
        fn is_valid(&self) -> bool {
            self.0 < MODULUS
        }
        /// Subtracts the modulus from this element if this element is not in the
        /// field. Only used interally.
        #[inline(always)]
        fn reduce(&mut self) {
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
            mut r7: u64,
        ) {
            let k = r0.wrapping_mul(INV);
            let mut carry = 0;
            crate::ff::mac_with_carry(r0, k, MODULUS.0[0], &mut carry);
            r1 = crate::ff::mac_with_carry(r1, k, MODULUS.0[1usize], &mut carry);
            r2 = crate::ff::mac_with_carry(r2, k, MODULUS.0[2usize], &mut carry);
            r3 = crate::ff::mac_with_carry(r3, k, MODULUS.0[3usize], &mut carry);
            r4 = crate::ff::adc(r4, 0, &mut carry);
            let carry2 = carry;
            let k = r1.wrapping_mul(INV);
            let mut carry = 0;
            crate::ff::mac_with_carry(r1, k, MODULUS.0[0], &mut carry);
            r2 = crate::ff::mac_with_carry(r2, k, MODULUS.0[1usize], &mut carry);
            r3 = crate::ff::mac_with_carry(r3, k, MODULUS.0[2usize], &mut carry);
            r4 = crate::ff::mac_with_carry(r4, k, MODULUS.0[3usize], &mut carry);
            r5 = crate::ff::adc(r5, carry2, &mut carry);
            let carry2 = carry;
            let k = r2.wrapping_mul(INV);
            let mut carry = 0;
            crate::ff::mac_with_carry(r2, k, MODULUS.0[0], &mut carry);
            r3 = crate::ff::mac_with_carry(r3, k, MODULUS.0[1usize], &mut carry);
            r4 = crate::ff::mac_with_carry(r4, k, MODULUS.0[2usize], &mut carry);
            r5 = crate::ff::mac_with_carry(r5, k, MODULUS.0[3usize], &mut carry);
            r6 = crate::ff::adc(r6, carry2, &mut carry);
            let carry2 = carry;
            let k = r3.wrapping_mul(INV);
            let mut carry = 0;
            crate::ff::mac_with_carry(r3, k, MODULUS.0[0], &mut carry);
            r4 = crate::ff::mac_with_carry(r4, k, MODULUS.0[1usize], &mut carry);
            r5 = crate::ff::mac_with_carry(r5, k, MODULUS.0[2usize], &mut carry);
            r6 = crate::ff::mac_with_carry(r6, k, MODULUS.0[3usize], &mut carry);
            r7 = crate::ff::adc(r7, carry2, &mut carry);
            (self.0).0[0usize] = r4;
            (self.0).0[1usize] = r5;
            (self.0).0[2usize] = r6;
            (self.0).0[3usize] = r7;
            self.reduce();
        }
    }
    impl crate::ff::SqrtField for Fr {
        fn legendre(&self) -> crate::ff::LegendreSymbol {
            let s = self.pow([
                11389680472494603939u64,
                14681934109093717318u64,
                15863968012492123182u64,
                1743499133401485332u64,
            ]);
            if s == Self::zero() {
                crate::ff::LegendreSymbol::Zero
            } else if s == Self::one() {
                crate::ff::LegendreSymbol::QuadraticResidue
            } else {
                crate::ff::LegendreSymbol::QuadraticNonResidue
            }
        }
        fn sqrt(&self) -> Option<Self> {
            let mut a1 = self.pow([
                5694840236247301969u64,
                7340967054546858659u64,
                7931984006246061591u64,
                871749566700742666u64,
            ]);
            let mut a0 = a1;
            a0.square();
            a0.mul_assign(self);
            if a0.0
                == FrRepr([
                    7548957153968385962u64,
                    10162512645738643279u64,
                    5900175412809962033u64,
                    2475245527108272378u64,
                ])
            {
                None
            } else {
                a1.mul_assign(self);
                Some(a1)
            }
        }
    }
}
