use crate::const_repr::BigintRepresentation;
use crate::const_repr::FullMultiplication;

use ff::*;

// #[macro_use]
// use crunchy::*;

pub struct PrimeFieldElement<
    P,
    const N: usize
>(pub BigintRepresentation<{N}>, std::marker::PhantomData<P>) 
    where BigintRepresentation<{N}>: FullMultiplication, 
        P: FieldParameters<{N}>, 
        [u64; N]: std::array::LengthAtMost32;

pub trait FieldParameters<const N:usize>: Sized + Copy +  Send + Sync + 'static {
    const NUM_BITS: u32;
    const CAPACITY: u32;
    const REPR_SHAVE_BITS: u32;
    const S: u32;
    const MULTIPLICATIVE_GENERATOR: BigintRepresentation<{N}>;
    const ROOT_OF_UNITY: BigintRepresentation<{N}>;
    const MODULUS: BigintRepresentation<{N}>;
    const R: BigintRepresentation<{N}>;
    const R2: BigintRepresentation<{N}>;
    const INV: u64;
}

impl<P, const N: usize> PrimeFieldElement<P, {N}> 
    where P: FieldParameters<{N}>, 
    [u64; N]: std::array::LengthAtMost32 
{
    #[inline(always)]
    fn is_valid(&self) -> bool {
        self.0 < P::MODULUS
    }

    #[inline(always)]
    fn reduce(&mut self) {
        if !self.is_valid() {
            self.0.sub_noborrow(&P::MODULUS);
        }
    }

    #[inline(always)]
    // fn mont_reduce(&mut self, mut mul_res: BigintRepresentation<{N*2}>) {
    fn mont_reduce(&mut self, mut mul_res:< BigintRepresentation<{N}> as FullMultiplication >::MulResult) {
        let mut carry2 = 0u64;
        let mut carry = 0u64;
        for j in 0..N {
            let k = mul_res.0[j].wrapping_mul(P::INV);
            for i in 0..N {
                mul_res.0[i + j] = ::ff::mac_with_carry(mul_res.0[i + j], k, P::MODULUS.0[i], &mut carry);
            }
            mul_res.0[N + j] = ::ff::adc(mul_res.0[{N} + j], carry2, &mut carry);
            carry2 = carry;
            carry = 0u64;
        }

        for j in 0..N {
            (self.0).0[j] = (mul_res.0)[N + j];
        }

        self.reduce();
    } 
}

impl<P, const N: usize> Copy for PrimeFieldElement<P, {N}> 
    where P: FieldParameters<{N}>, 
    [u64; N]: std::array::LengthAtMost32 {}

impl<P, const N: usize> Clone for PrimeFieldElement<P, {N}> 
    where P: FieldParameters<{N}>, 
    [u64; N]: std::array::LengthAtMost32 {
    fn clone(&self) -> Self {
        *self
    }
}

impl<P, const N: usize> std::cmp::PartialEq for PrimeFieldElement<P, {N}> 
    where P: FieldParameters<{N}>, 
    [u64; N]: std::array::LengthAtMost32 {
    fn eq(&self, other: &Self) -> bool {
        self.0.cmp(&other.0) == std::cmp::Ordering::Equal
    }
}

impl<P, const N: usize> std::cmp::Eq for PrimeFieldElement<P, {N}> 
    where P: FieldParameters<{N}>, 
    [u64; N]: std::array::LengthAtMost32 {}

impl<P, const N: usize> std::fmt::Debug for PrimeFieldElement<P, {N}> 
    where P: FieldParameters<{N}>, 
    [u64; N]: std::array::LengthAtMost32
{
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fp({})", self.into_repr())
    }
}

impl<P, const N: usize> ::rand::Rand for PrimeFieldElement<P, {N}> 
    where P: FieldParameters<{N}>, 
    [u64; N]: std::array::LengthAtMost32 {

    #[inline(always)]
    fn rand<R: ::rand::Rng>(rng: &mut R) -> Self {
        let s = BigintRepresentation::<{N}>::rand(rng);
        // TODO: shave
        Self(s, std::marker::PhantomData)
    }
}

impl<P, const N: usize> std::fmt::Display for PrimeFieldElement<P, {N}> 
    where P: FieldParameters<{N}>, 
    [u64; N]: std::array::LengthAtMost32 {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fp({})", self.into_repr())
    }
}

impl<P, const N: usize> Field for PrimeFieldElement<P, {N}> 
    where P: FieldParameters<{N}>, 
    [u64; N]: std::array::LengthAtMost32 
{
    #[inline(always)]
    fn zero() -> Self {
        Self(BigintRepresentation::<{N}>::default(), std::marker::PhantomData)
    }

    #[inline(always)]
    fn one() -> Self {
        Self(P::R, std::marker::PhantomData)
    }

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    #[inline(always)]
    fn add_assign(&mut self, other: &Self) {
        self.0.add_nocarry(&other.0);
        self.reduce();
    }

    #[inline(always)]
    fn double(&mut self) {
        self.0.mul2();
        self.reduce();
    }

    #[inline(always)]
    fn sub_assign(&mut self, other: &Self) {

        if other.0 > self.0 {
            self.0.add_nocarry(&P::MODULUS);
        }

        self.0.sub_noborrow(&other.0);
    }

    #[inline(always)]
    fn negate(&mut self) {
        if !self.is_zero() {
            let mut tmp = P::MODULUS;
            tmp.sub_noborrow(&self.0);
            self.0 = tmp;
        }
    }

    fn inverse(&self) -> Option<Self> {
        None
        // if self.is_zero() {
        //     None
        // } else {
        //     // Guajardo Kumar Paar Pelzl
        //     // Efficient Software-Implementation of Finite Fields with Applications to Cryptography
        //     // Algorithm 16 (BEA for Inversion in Fp)

        //     let one = #repr::from(1);

        //     let mut u = self.0;
        //     let mut v = MODULUS;
        //     let mut b = #name(R2); // Avoids unnecessary reduction step.
        //     let mut c = Self::zero();

        //     while u != one && v != one {
        //         while u.is_even() {
        //             u.div2();

        //             if b.0.is_even() {
        //                 b.0.div2();
        //             } else {
        //                 b.0.add_nocarry(&MODULUS);
        //                 b.0.div2();
        //             }
        //         }

        //         while v.is_even() {
        //             v.div2();

        //             if c.0.is_even() {
        //                 c.0.div2();
        //             } else {
        //                 c.0.add_nocarry(&MODULUS);
        //                 c.0.div2();
        //             }
        //         }

        //         if v < u {
        //             u.sub_noborrow(&v);
        //             b.sub_assign(&c);
        //         } else {
        //             v.sub_noborrow(&u);
        //             c.sub_assign(&b);
        //         }
        //     }

        //     if u == one {
        //         Some(b)
        //     } else {
        //         Some(c)
        //     }
        // }
    }

    #[inline(always)]
    fn frobenius_map(&mut self, _: usize) {
        // This has no effect in a prime field.
    }

    #[inline(always)]
    fn mul_assign(&mut self, other: &Self) {
        // let mut interm = BigintRepresentation::<{N*2}>::default();
        let mut interm = < BigintRepresentation<{N}> as FullMultiplication >::MulResult::default();
        let mut carry = 0u64;
        for j in 0..N {
            let this_limb = (self.0).0[j];
            for i in 0..N {
                interm.0[i + j] = ::ff::mac_with_carry(interm.0[i + j], this_limb, (other.0).0[i], &mut carry);
            }
            interm.0[N + j] = carry;
            carry = 0u64;
        }

        self.mont_reduce(interm);
    }

    #[inline(always)]
    fn square(&mut self) {
        // let mut interm = BigintRepresentation::<{N*2}>::default();
        let mut interm = < BigintRepresentation<{N}> as FullMultiplication >::MulResult::default();
        let mut carry = 0u64;

        for j in 0..N {
            let this_limb = (self.0).0[j];
            for i in (j+1)..N {
                interm.0[i + j] = ::ff::mac_with_carry(interm.0[i + j], this_limb, (self.0).0[i], &mut carry);
            }
            interm.0[N + j] = carry;
            carry = 0u64;
        }

        interm.0[2*N - 1] = interm.0[2*N - 2] >> 63;

        for j in (2..=(2*N - 2)).rev() {
            interm.0[j] = (interm.0[j] << 1) | (interm.0[j-1] >> 63);
        }

        interm.0[1] = interm.0[1] << 1;
        for j in 0..N {
            let this_limb = (self.0).0[j];
            let idx = 2*j;
            interm.0[idx] = ::ff::mac_with_carry(interm.0[idx], this_limb, this_limb, &mut carry);
            interm.0[idx+1] = ::ff::adc(interm.0[idx+1], 0u64, &mut carry);
        }

        self.mont_reduce(interm);
    }
}

impl<P, const N: usize> From<PrimeFieldElement<P, {N}>> for BigintRepresentation<{N}>
    where P: FieldParameters<{N}>, 
    [u64; N]: std::array::LengthAtMost32 
{
    fn from(e: PrimeFieldElement<P, {N}>) -> Self {
        e.into_repr()
    }
}

impl<P, const N: usize> PrimeField for PrimeFieldElement<P, {N}> 
    where P: FieldParameters<{N}>, 
    [u64; N]: std::array::LengthAtMost32 
{
    const NUM_BITS: u32 = P::NUM_BITS;
    const CAPACITY: u32 = P::CAPACITY;
    const S: u32 = P::S;

    type Repr = BigintRepresentation<{N}>;

    fn from_repr(repr: Self::Repr) -> Result<Self, PrimeFieldDecodingError> {
        let mut r = Self(repr, std::marker::PhantomData);
        if r.is_valid() {
            r.mul_assign(&Self(P::R2, std::marker::PhantomData));

            Ok(r)
        } else {
            Err(PrimeFieldDecodingError::NotInField(format!("{}", r.0)))
        }
    }

    // fn from_raw_repr(repr: Self::Repr) -> Result<Self, PrimeFieldDecodingError> {
    //     let r = Self(repr, std::marker::PhantomData);
    //     if r.is_valid() {
    //         Ok(r)
    //     } else {
    //         Err(PrimeFieldDecodingError::NotInField(format!("{}", r.0)))
    //     }
    // }

    fn into_repr(&self) -> Self::Repr {
        let mut r = *self;
        // let mut interm = BigintRepresentation::<{N*2}>::default();
        let mut interm = < BigintRepresentation<{N}> as FullMultiplication >::MulResult::default();
        for j in 0..N {
            interm.0[j] = (self.0).0[j];
        }
        r.mont_reduce(interm);

        r.0
    }

    // fn into_raw_repr(&self) -> Self::Repr {
    //     self.0
    // }

    fn char() -> Self::Repr {
        P::MODULUS
    }

    fn multiplicative_generator() -> Self {
        Self(P::MULTIPLICATIVE_GENERATOR, std::marker::PhantomData)
    }

    fn root_of_unity() -> Self {
        Self(P::ROOT_OF_UNITY, std::marker::PhantomData)
    }

}