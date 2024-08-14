use ff;
use rand;

pub struct BigintRepresentation<
    const N: usize
>(pub [u64; N]);

impl<const N: usize> Copy for BigintRepresentation<{N}> {}
impl<const N: usize> Clone for BigintRepresentation<{N}> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<const N: usize> std::cmp::PartialEq for BigintRepresentation<{N}> {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(&other) == std::cmp::Ordering::Equal
    }
}

impl<const N: usize> std::cmp::Eq for BigintRepresentation<{N}> {}

impl<const N: usize> std::fmt::Debug for BigintRepresentation<{N}> 
{
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "0x")?;
        for i in self.0.iter().rev() {
            write!(f, "{:016x}", *i)?;
        }

        Ok(())
    }
}

impl<const N: usize> ::rand::Rand for BigintRepresentation<{N}> {
    #[inline(always)]
    fn rand<R: ::rand::Rng>(rng: &mut R) -> Self {
        let mut s = Self::default();
        for el in s.0.iter_mut() {
            *el = rng.gen();
        }

        s
    }
}

impl<const N: usize> std::fmt::Display for BigintRepresentation<{N}> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "0x")?;
        for i in self.0.iter().rev() {
            write!(f, "{:016x}", *i)?;
        }

        Ok(())
    }
}

impl<const N: usize> std::default::Default for BigintRepresentation<{N}> {
    #[inline(always)]
    fn default() -> Self {
        use std::mem::{MaybeUninit};
        let mut s: Self = unsafe {MaybeUninit::uninit().assume_init() };
        for el in s.0.iter_mut() {
            *el = 0u64;
        }

        s
        // let repr: [u64; {N}] = [0u64; {N}];
        // BigintRepresentation::<{N}>(repr)
    }
}

impl<const N: usize> AsRef<[u64]> for BigintRepresentation<{N}> {
    #[inline(always)]
    fn as_ref(&self) -> &[u64] {
        &self.0
    }
}

impl<const N: usize> AsMut<[u64]> for BigintRepresentation<{N}> {
    #[inline(always)]
    fn as_mut(&mut self) -> &mut [u64] {
        &mut self.0
    }
}

impl<const N: usize> From<u64> for BigintRepresentation<{N}> {
    #[inline(always)]
    fn from(val: u64) -> Self {
        let mut repr = Self::default();
        repr.0[0] = val;
        repr
    }
}

impl<const N: usize> Ord for BigintRepresentation<{N}> {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        for (a, b) in self.0.iter().rev().zip(other.0.iter().rev()) {
            if a < b {
                return std::cmp::Ordering::Less
            } else if a > b {
                return std::cmp::Ordering::Greater
            }
        }

        std::cmp::Ordering::Equal
    }
}

impl<const N: usize> PartialOrd for BigintRepresentation<{N}> {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<const N: usize> ::ff::PrimeFieldRepr for BigintRepresentation<{N}> where [u64; N]: std::array::LengthAtMost32 {
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
        if n as usize >= 64 * N {
            *self = Self::from(0);
            return;
        }

        while n >= 64 {
            let mut t = 0;
            for i in self.0.iter_mut().rev() {
                std::mem::swap(&mut t, i);
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
        if n as usize >= 64 * N {
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
        let mut ret = (N as u32) * 64;
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
    fn add_nocarry(&mut self, other: &Self) {
        let mut carry = 0;

        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            *a = ::ff::adc(*a, *b, &mut carry);
        }
    }

    #[inline(always)]
    fn sub_noborrow(&mut self, other: &Self) {
        let mut borrow = 0;

        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            *a = ::ff::sbb(*a, *b, &mut borrow);
        }
    }
}

pub trait FullMultiplication {
    type Multiplicand;
    type MulResult;
}

impl<const N: usize> FullMultiplication for BigintRepresentation<{N}> 
{
    type Multiplicand = BigintRepresentation<{N}>;
    type MulResult = BigintRepresentation<{N*2}>;
}

