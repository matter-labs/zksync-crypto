use pairing::ff::Field;
use pairing::{ff::PrimeField, BitIterator};
use std::sync::Arc;

use super::{fq12::Fq12, fq2::Fq2, fq6::Fq6, params::TorusExtension12Params};
use crate::config::{CSConfig, CSWitnessEvaluationConfig};
use crate::gadgets::non_native_field::implementations::NonNativeFieldOverU16;
use crate::gadgets::tower_extension::params::{Extension2Params, Extension6Params};
use crate::gadgets::traits::witnessable::WitnessHookable;
use crate::{
    cs::traits::cs::ConstraintSystem,
    field::SmallField,
    gadgets::{
        boolean::Boolean,
        non_native_field::traits::NonNativeField,
        traits::{hardexp_compatible::HardexpCompatible, selectable::Selectable},
    },
};

/// [`TorusWrapper`] is an algebraic compression of the `Fq12` element via underlying encoding of `Fq6`.
/// In compressed form operations over Fq12 are less expensive.
///
/// The implementation is based on the following paper:
/// https://eprint.iacr.org/2022/1162.pdf.
#[derive(Clone, Debug, Copy)]
pub struct TorusWrapper<F, T, NN, P>
where
    F: SmallField,
    T: PrimeField,
    NN: NonNativeField<F, T>,
    P: TorusExtension12Params<T>,
{
    pub encoding: Fq6<F, T, NN, P::Ex6>,
}

// TODO: Probably, this could be implemented generally for any two Fqk and Fq(k/2) elements.
impl<F, T, P, const N: usize> TorusWrapper<F, T, NonNativeFieldOverU16<F, T, N>, P>
where
    F: SmallField,
    T: PrimeField,
    P: TorusExtension12Params<T>,
    [(); N + 1]:,
{
    /// Creates a new instance of the [`TorusWrapper`] with the given encoding.
    pub fn new(encoding: Fq6<F, T, NonNativeFieldOverU16<F, T, N>, P::Ex6>) -> Self {
        Self { encoding }
    }

    pub fn one<CS>(
        cs: &mut CS,
        params: &Arc<<NonNativeFieldOverU16<F, T, N> as NonNativeField<F, T>>::Params>,
    ) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        let encoding = Fq6::one(cs, params);
        Self::new(encoding)
    }

    pub fn zero<CS>(
        cs: &mut CS,
        params: &Arc<<NonNativeFieldOverU16<F, T, N> as NonNativeField<F, T>>::Params>,
    ) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        let encoding = Fq6::zero(cs, params);
        Self::new(encoding)
    }

    /// Returns the underlying parameters of the encoded `Fq6` element.
    pub fn get_params(
        &self,
    ) -> &Arc<<NonNativeFieldOverU16<F, T, N> as NonNativeField<F, T>>::Params> {
        self.encoding.get_params()
    }

    /// Normalizes the encoding of the `Fq6` element.
    pub fn normalize<CS>(&mut self, cs: &mut CS)
    where
        CS: ConstraintSystem<F>,
    {
        self.encoding.normalize(cs);
    }

    /// Returns an instance if `flag` is `true`, otherwise returns a zero element.
    pub fn mask<CS>(&mut self, cs: &mut CS, flag: Boolean<F>) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        let zero = Fq6::zero(cs, self.get_params());
        let new_encoding =
            <Fq6<F, T, NonNativeFieldOverU16<F, T, N>, P::Ex6>>::conditionally_select(
                cs,
                flag,
                &self.encoding,
                &zero,
            );

        Self::new(new_encoding)
    }
    pub fn replace_by_constant_if_trivial<CS: ConstraintSystem<F>>(
        &mut self,
        cs: &mut CS,
        cnst: P::Witness,
    ) -> (Self, Boolean<F>) {
        let params = self.get_params();
        let mut elem = self.encoding.clone();
    
        let is_trivial = elem.is_zero(cs);
    
        let (c0, c1) = P::convert_from_structured_witness(cnst);
        let mut compressed_cnst = c1;
        compressed_cnst.mul_assign(&c0);
        compressed_cnst.negate();
    
        let constant_encoding = Fq6::constant(cs, compressed_cnst, params);
    
        // Conditionally select either this constant encoding (if trivial) or the original encoding.
        let new_encoding = <Fq6<F, T, NonNativeFieldOverU16<F, T, N>, P::Ex6> as NonNativeField<F, T>>::conditionally_select(
            cs, 
            is_trivial, 
            &constant_encoding, 
            &self.encoding
        );
 
        (Self::new(new_encoding), is_trivial)
    }
    

    /// Compresses the `Fq12` element `c0 + c1*w` to the Torus (`T2`) element.
    ///
    /// Uses the formula `m <- (1 + c0) / c1` to compress the `Fq12` element with the additional
    /// check for the exceptional case when `c1` is zero.
    ///
    /// If `SAFE=false`, then the function will not check for the exceptional case when `c1` is zero.
    pub fn compress<CS>(
        cs: &mut CS,
        f: &mut Fq12<F, T, NonNativeFieldOverU16<F, T, N>, P>,
        is_safe: bool
    ) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        let params = f.get_params();
        let mut c0 = f.c0.clone();
        let mut c1 = f.c1.clone();

        let mut encoding = if is_safe {
            // Preparing flags for exception cases
            let is_exceptional = Fq6::is_zero(&mut c1, cs);
            let mut c0_is_one = Fq6::one(cs, params);
            let c0_is_one = c0_is_one.equals(cs, &mut c0);
            let mut is_exceptional = Fq6::from_boolean(cs, is_exceptional, params);
            let mut c0_is_one = Fq6::from_boolean(cs, c0_is_one, params);

            // m <- (1 + c0) / c1 if c1 is non-zero. However, to account for the case where
            // c1 is zero, we set numerator to 1 + c0 - 2*c0_is_one and denominator to c1 + is_exceptional.
            let mut numerator = Fq6::one(cs, params);
            let mut numerator = numerator.add(cs, &mut c0);
            let mut c0_is_one_doubled = c0_is_one.double(cs);
            let mut numerator = numerator.sub(cs, &mut c0_is_one_doubled);
            let mut denominator = f.c1.add(cs, &mut is_exceptional);
            denominator.normalize(cs);

            let encoding = numerator.div(cs, &mut denominator);
            encoding
        } else {
            // Verifying that c1 is non-zero
            let boolean_false = Boolean::allocated_constant(cs, false);
            let c1_is_zero = c1.is_zero(cs);
            Boolean::enforce_equal(cs, &c1_is_zero, &boolean_false);

            // m <- (1 + c0) / c1
            let mut encoding = Fq6::one(cs, params);
            let mut encoding = encoding.add(cs, &mut f.c0);
            let encoding = encoding.div(cs, &mut f.c1);

            encoding
        };

        encoding.normalize(cs);
        Self::new(encoding)
    }

    /// Decompresses the Torus (`T2`) element `g` back to the `Fq12` element by using the formula
    ///
    /// `zeta^{-1} = (g + w)/(g - w)`
    pub fn decompress<CS>(&self, cs: &mut CS) -> Fq12<F, T, NonNativeFieldOverU16<F, T, N>, P>
    where
        CS: ConstraintSystem<F>,
    {
        let params = self.get_params();
        let mut one = Fq6::one(cs, params);
        let negative_one = one.negated(cs);

        // Since `g` is a pure `Fq6` element, `g+w` is just an `Fq12` element with `c0 = g` and `c1 = 1`.
        let mut numerator = Fq12::new(self.encoding.clone(), one);
        // Since `g` is a pure `Fq6` element, `g-w` is just an `Fq12` element with `c0 = g` and `c1 = -1`.
        let mut denominator = Fq12::new(self.encoding.clone(), negative_one);

        // zeta^{-1} = (g + w)/(g - w)
        let decompressed = numerator.div(cs, &mut denominator);

        decompressed
    }

    /// Computes the inverse of the Torus element using the formula g -> -g.
    pub fn inverse<CS>(&mut self, cs: &mut CS) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        let encoding = self.encoding.negated(cs);
        Self::new(encoding)
    }

    /// Computes the conjugate of the Torus element using the formula g -> -g.
    /// Note that the conjugate of the Torus element is the same as its inverse.
    pub fn conjugate<CS>(&mut self, cs: &mut CS) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        self.inverse(cs)
    }

    /// Computes the Frobenius map of the Torus element with the given power using the formula
    ///
    /// frob_map(g, i) = g^(p^i) / \gamma^{(p^i-1)/2}
    pub fn frobenius_map<CS>(&mut self, cs: &mut CS, power: usize) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        // We compute frobenius map unconstrained:
        let witness_self = self.encoding_to_witness(cs);
        let witness_frob = P::torus_frobenius_map(witness_self, power);

        // Now, we constraint the frobenius map with a cheaper version:
        // Suppose r = f(g,i) / (f(w,i) * w^{-1}). Then, we require:
        // f(g, i) = f(w, i) * (w^{-1}) * r
        // Notice that `f(w,i)*w^{-1}` must yield an element
        // from Fq6. Thus, we need one frobenius map + mul over Fq6, and
        // one frobenius map + mul over Fq12.
        let params = self.encoding.get_params();
        let mut encoding_new = Fq6::allocate_from_witness(cs, witness_frob, params);

        // rhs = f(w, i) * (w^{-1}) * r
        // First, allocating the w^{-1}
        let w_inverse = P::get_w_inverse_coeffs_c5();
        let mut w_inverse: Fq2<_, _, _, <P::Ex6 as Extension6Params<T>>::Ex2> =
            Fq2::constant(cs, w_inverse, params);

        let mut rhs: Fq12<F, T, NonNativeFieldOverU16<F, T, N>, P> =
            Fq12::one_imaginary(cs, params);
        rhs = rhs.frobenius_map(cs, power);
        rhs = rhs.mul_by_c5(cs, &mut w_inverse);

        // Asserting that c1 is zero since rhs must be a pure Fq6 element at this point.
        let boolean_true = Boolean::allocated_constant(cs, true);
        let c1_is_zero = rhs.c1.is_zero(cs);
        Boolean::enforce_equal(cs, &c1_is_zero, &boolean_true);
        let mut rhs = rhs.c0.clone();

        // Finishing rhs by multiplying by result
        rhs = rhs.mul(cs, &mut encoding_new);

        // lhs = f(g, i)
        let mut lhs = self.encoding.clone();
        lhs = lhs.frobenius_map(cs, power);

        // Asserting that lhs == rhs
        Fq6::enforce_equal(cs, &lhs, &rhs);

        Self::new(encoding_new)
    }

    /// Computes the product of two Torus elements using the formula
    ///
    /// `(g, g') -> (g * g' + \gamma) / (g + g')`
    ///
    /// The formula handles the exceptional case when `g + g'` is zero.
    pub fn mul<CS>(&mut self, cs: &mut CS, other: &mut Self,) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        // We compute multiplication unconstrained:
        let witness_self = self.encoding_to_witness(cs);
        let witness_other = other.encoding_to_witness(cs);
        let witness_mul = P::torus_mul(witness_self, witness_other);

        // Now, we constraint the multiplication with a cheaper version:
        // g'' = (g * g' + \gamma) / (g + g') is equivalent to
        // g'' * (g + g') = (g * g' + \gamma)
        // Here, g'' is the new encoding.
        let params = self.encoding.get_params();
        let encoding_new = Fq6::allocate_from_witness(cs, witness_mul, params);

        // lhs = g'' * (g + g')
        let mut sum = self.encoding.clone().add(cs, &mut other.encoding);
        let lhs = encoding_new.clone().mul(cs, &mut sum);

        // rhs = g * g' + \gamma
        let mut gamma = Fq6::gamma(cs, params);
        let mut rhs = self.encoding.clone().mul(cs, &mut other.encoding);
        let rhs = rhs.add(cs, &mut gamma);

        // Enforce equality
        Fq6::enforce_equal(cs, &lhs, &rhs);

        Self::new(encoding_new)
    }
    pub fn mul_optimal<CS>(&mut self, cs: &mut CS, other: &mut Self, is_safe: bool) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        let params = self.get_params();
        let mut gamma = Fq6::gamma(cs, params);
    
        if is_safe {
            // exceptions in case g2 = - g1
            // modified formula looks like (here flag = exception_flag):
            // x = g1 * g2 + \gamma
            // g = (x - flag * x) / (g1 + g2 + flag)
            let mut lhs = self.encoding.clone();
            let mut rhs = other.encoding.clone().negated(cs);
            let exc_flag = lhs.equals(cs, &mut rhs);
    
            let flag_as_fe = Fq6::from_boolean(cs, exc_flag, params);
    
            // x = g1 * g2 + gamma
            let mut x = self.encoding.clone().mul(cs, &mut other.encoding.clone());
            x = x.add(cs, &mut gamma);
    
            // If exceptional, we want to subtract x from itself, otherwise subtract zero
            let zero = Fq6::zero(cs, params);
            let mut y = <Fq6<F, T, NonNativeFieldOverU16<F, T, N>, P::Ex6>>::conditionally_select(cs, exc_flag, &x, &zero);
    
            let mut numerator = x.sub(cs, &mut y);
    
            let mut denominator = self.encoding.clone().add(cs, &mut other.encoding.clone());
            denominator = denominator.add(cs, &mut flag_as_fe.clone());
    
            let encoding = numerator.div(cs, &mut denominator);
    
            Self::new(encoding)
        } else {
            // g = (g1 * g2 + \gamma) / (g1 + g2)
            // assume that are in the exceptional case: g2 = -g1
            // we are going to enforce relation of the form: g * 0 = g1 * g2 + \gamma
            // unless g1 * g2 + \gamma == 0 g would be never underconstrained
            // if g1 * g2 + \gamma = \gamma - g1^2 = 0 and hence g1 is the root of polynomial X^2 - \gamma = 0,
            // and hence this poly is not irreducible - contradiction with F_q^2 = F_q[w] / (w^2 - \gamma)
            // This means, we are completely safe here and no additional checks are requierd
            let mut numerator = self.encoding.clone().mul(cs, &mut other.encoding.clone());
            numerator = numerator.add(cs, &mut gamma);
    
            let mut denominator = self.encoding.clone().add(cs, &mut other.encoding.clone());
    
            let encoding = numerator.div(cs, &mut denominator);
            Self::new(encoding)
        } 
    }

    pub fn pow_naf_decomposition<CS, S: AsRef<[i8]>>(
        &mut self,
        cs: &mut CS,
        decomposition: S,
        is_safe: bool
    ) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        // Intializing the result with 1
        let mut result = self.clone();

        // Preparing self and self inverse in advance
        let mut self_cloned = self.clone();
        let mut self_inverse = self.conjugate(cs);

        for bit in decomposition.as_ref().iter().skip(1) {
            result = result.square(cs);

            // If bit is 1, multiply by initial torus
            let bit_is_one = Boolean::allocated_constant(cs, *bit == 1);
            let result_times_self = result.mul_optimal(cs, &mut self_cloned, is_safe);
            result = Self::conditionally_select(cs, bit_is_one, &result_times_self, &result);

            // If bit is -1, multiply by inverse initial torus
            let bit_is_minus_one = Boolean::allocated_constant(cs, *bit == -1);
            let result_times_self_inverse = result.mul_optimal(cs, &mut self_inverse, is_safe);
            result = Self::conditionally_select(
                cs,
                bit_is_minus_one,
                &result_times_self_inverse,
                &result,
            );
        }

        result
    }

    pub fn pow_u32<CS, S: AsRef<[u64]>>(&mut self, cs: &mut CS, exponent: S) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        let mut result = Self::one(cs, self.get_params());
        let mut base = self.clone();

        for i in BitIterator::new(exponent) {
            let mut squared = result.square(cs);
            let mut squared_and_multiplied = squared.mul(cs, &mut base);
            let shall_multiply = Boolean::allocated_constant(cs, i);

            result = Self::conditionally_select(
                cs,
                shall_multiply,
                &mut squared_and_multiplied,
                &mut squared,
            );

            result.normalize(cs);
        }

        result
    }

    pub fn square<CS>(&mut self, cs: &mut CS) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        // We compute squaring unconstrained:
        let witness = self.encoding_to_witness(cs);
        let witness_squared = P::torus_square(witness);

        // Now, we constraint squaring with a cheaper version:
        // g' = (1/2)(g + \gamma/g) is equivalent to
        // (2g' - g)*g = gamma
        let params = self.encoding.get_params();
        let encoding_new = Fq6::allocate_from_witness(cs, witness_squared, params);

        // lhs = (2g' - g)*g
        let mut lhs = encoding_new.clone();
        lhs = lhs.double(cs);
        lhs = lhs.sub(cs, &mut self.encoding.clone());
        let lhs = self.encoding.clone().mul(cs, &mut lhs);

        // rhs = g == 0 ? zero : gamma
        let zero = Fq6::zero(cs, params);
        let gamma = Fq6::gamma(cs, params);
        let is_zero_g = self.encoding.is_zero(cs);
        let rhs = <Fq6<F, T, NonNativeFieldOverU16<F, T, N>, P::Ex6>>::conditionally_select(
            cs, is_zero_g, &zero, &gamma,
        );

        // We can just enforce equality without subbing
        Fq6::enforce_equal(cs, &lhs, &rhs);
        Self::new(encoding_new)
    }

    pub fn square_optimal<CS>(&mut self, cs: &mut CS, is_safe: bool) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        let params = self.get_params();
        let mut gamma = Fq6::gamma(cs, params);
        let mut g = self.encoding.clone();
    
        // exception_free formula looks like (here flag := is_exceptional)
        // res = 1/2 (g + [(\gamma * flag!) / (g + flag)])
        // unsafe formula is : res = 1/2 (g + \gamma / g);
        // we are going to do with them simultaneouly, rewriting the formula as: res = 1/2 (g + tmp)
        // where tmp := (\gamma * flag!) / (g + flag) in the first case and tmp := \gamma / g in the second
        let tmp = if is_safe {
            let is_exceptional = g.is_zero(cs);
    
            // denom = g + flag
            let flag_as_fe = Fq6::from_boolean(cs, is_exceptional, params);
            let mut denom = g.clone().add(cs, &mut flag_as_fe.clone());
    
            // numerator = gamma if not exceptional, else zero
            let zero = Fq6::zero(cs, params);
            let tmp = is_exceptional.negated(cs);
            let mut numerator = <Fq6<F, T, NonNativeFieldOverU16<F, T, N>, P::Ex6>>::conditionally_select(cs, tmp, &gamma, &zero);
    
            numerator.div(cs, &mut denom)
        } else {
            // tmp = gamma / g
            let mut g_clone = g.clone();
            gamma.div(cs, &mut g_clone)
        };
    
        // res = (g + tmp)/2
        let mut sum = g.clone().add(cs, &mut tmp.clone());
    
        let mut two = Fq6::one(cs, params);
        two = two.double(cs);
        let mut inv_2 = two.inverse(cs);
    
        let mut res = sum.clone().mul(cs, &mut inv_2);

        let double_res = res.clone().double(cs);
        Fq6::enforce_equal(cs, &sum, &double_res);
    
        Self::new(res)
    }

    /// Converts the encoding of the `Fq6` element to the structured witness.
    pub(super) fn encoding_to_witness<CS>(
        &self,
        cs: &mut CS,
    ) -> <P::Ex6 as Extension6Params<T>>::Witness
    where
        CS: ConstraintSystem<F>,
    {
        match self.encoding.witness_hook(cs)() {
            Some((c0, c1, c2)) => {
                let (c0_c0, c0_c1) = c0;
                let (c1_c0, c1_c1) = c1;
                let (c2_c0, c2_c1) = c2;

                let (c0_c0, c0_c1) = (c0_c0.get(), c0_c1.get());
                let (c1_c0, c1_c1) = (c1_c0.get(), c1_c1.get());
                let (c2_c0, c2_c1) = (c2_c0.get(), c2_c1.get());

                let c0 = <P::Ex6 as Extension6Params<T>>::Ex2::convert_to_structured_witness(
                    c0_c0, c0_c1,
                );
                let c1 = <P::Ex6 as Extension6Params<T>>::Ex2::convert_to_structured_witness(
                    c1_c0, c1_c1,
                );
                let c2 = <P::Ex6 as Extension6Params<T>>::Ex2::convert_to_structured_witness(
                    c2_c0, c2_c1,
                );

                P::Ex6::convert_to_structured_witness(c0, c1, c2)
            }
            None => {
                assert!(!<CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS);
                P::one()
            }
        }
    }
}

impl<F, T, P, const N: usize> Selectable<F>
    for TorusWrapper<F, T, NonNativeFieldOverU16<F, T, N>, P>
where
    F: SmallField,
    T: PrimeField,
    P: TorusExtension12Params<T>,
    [(); N + 1]:,
{
    fn conditionally_select<CS>(cs: &mut CS, flag: Boolean<F>, a: &Self, b: &Self) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        let encoding = <Fq6<F, T, NonNativeFieldOverU16<F, T, N>, P::Ex6>>::conditionally_select(
            cs,
            flag,
            &a.encoding,
            &b.encoding,
        );

        Self::new(encoding)
    }
}

impl<F, T, const N: usize, P> HardexpCompatible<F>
    for TorusWrapper<F, T, NonNativeFieldOverU16<F, T, N>, P>
where
    F: SmallField,
    T: PrimeField,
    P: TorusExtension12Params<T>,
    [(); N + 1]:,
{
    fn mul<CS>(&mut self, cs: &mut CS, other: &mut Self) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        self.mul(cs, other)
    }

    fn square<CS>(&mut self, cs: &mut CS) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        self.square(cs)
    }

    fn conjugate<CS>(&mut self, cs: &mut CS) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        self.conjugate(cs)
    }

    fn inverse<CS>(&mut self, cs: &mut CS) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        self.inverse(cs)
    }

    fn frobenius_map<CS>(&mut self, cs: &mut CS, power: usize) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        self.frobenius_map(cs, power)
    }

    fn pow_u32<CS, S: AsRef<[u64]>>(&mut self, cs: &mut CS, exponent: S) -> Self
    where
        CS: ConstraintSystem<F>,
    {
        self.pow_u32(cs, exponent)
    }

    fn normalize<CS>(&mut self, cs: &mut CS)
    where
        CS: ConstraintSystem<F>,
    {
        self.normalize(cs);
    }
}
