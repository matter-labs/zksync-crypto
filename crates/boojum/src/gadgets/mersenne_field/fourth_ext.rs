use mersenne_field::{Mersenne31Complex, Mersenne31Quartic};

use super::second_ext::*;
use super::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct MersenneQuartic<F: SmallField> {
    pub x: MersenneComplex<F>,
    pub y: MersenneComplex<F>,
}

impl<F: SmallField> MersenneQuartic<F> {
    pub fn zero<CS: ConstraintSystem<F>>(cs: &mut CS) -> Self {
        Self {
            x: MersenneComplex::zero(cs),
            y: MersenneComplex::zero(cs),
        }
    }

    pub fn one<CS: ConstraintSystem<F>>(cs: &mut CS) -> Self {
        Self {
            x: MersenneComplex::one(cs),
            y: MersenneComplex::zero(cs),
        }
    }

    pub fn minus_one<CS: ConstraintSystem<F>>(cs: &mut CS) -> Self {
        Self {
            x: MersenneComplex::minus_one(cs),
            y: MersenneComplex::zero(cs),
        }
    }

    pub fn get_modulus_num<CS: ConstraintSystem<F>>(cs: &mut CS) -> Num<F> {
        Num::from_variable(cs.allocate_constant(F::from_u64_unchecked(M31_MODULUS)))
    }

    pub fn get_variables(&self) -> [Variable; 4] {
        [
            self.x.x.variable,
            self.x.y.variable,
            self.y.x.variable,
            self.y.y.variable,
        ]
    }

    pub fn into_nums(&self) -> [Num<F>; 4] {
        [
            self.x.x.into_num(),
            self.x.y.into_num(),
            self.y.x.into_num(),
            self.y.y.into_num(),
        ]
    }

    pub fn into_uint32s(&self) -> [UInt32<F>; 4] {
        [
            self.x.x.into_uint32(),
            self.x.y.into_uint32(),
            self.y.x.into_uint32(),
            self.y.y.into_uint32(),
        ]
    }

    pub fn into_coeffs(&self) -> [MersenneField<F>; 4] {
        [self.x.x, self.x.y, self.y.x, self.y.y]
    }

    pub fn from_coeffs(coefficients: [MersenneField<F>; 4]) -> Self {
        Self {
            x: MersenneComplex::from_coeffs(coefficients[0..2].try_into().unwrap()),
            y: MersenneComplex::from_coeffs(coefficients[2..4].try_into().unwrap()),
        }
    }

    /// The coordinate values should be in range [0, 2^31 - 2]
    fn from_variables_checked<CS: ConstraintSystem<F>>(
        cs: &mut CS,
        variables: [Variable; 4],
        reduced: bool,
    ) -> Self {
        Self {
            x: MersenneComplex::from_variables_checked(cs, [variables[0], variables[1]], reduced),
            y: MersenneComplex::from_variables_checked(cs, [variables[2], variables[3]], reduced),
        }
    }

    /// The coordinate values should be in range [0, 2^31 - 2]
    fn allocate_checked_without_value<CS: ConstraintSystem<F>>(cs: &mut CS, reduced: bool) -> Self {
        Self {
            x: MersenneComplex::allocate_checked_without_value(cs, reduced),
            y: MersenneComplex::allocate_checked_without_value(cs, reduced),
        }
    }

    pub fn allocate_checked<CS: ConstraintSystem<F>>(
        cs: &mut CS,
        witness: Mersenne31Quartic,
        reduced: bool,
    ) -> Self {
        Self {
            x: MersenneComplex::allocate_checked(cs, witness.c0, reduced),
            y: MersenneComplex::allocate_checked(cs, witness.c1, reduced),
        }
    }

    pub fn enforce_reduced<CS: ConstraintSystem<F>>(&mut self, cs: &mut CS) {
        self.x.enforce_reduced(cs);
        self.y.enforce_reduced(cs);
    }

    pub fn from_base<CS: ConstraintSystem<F>>(cs: &mut CS, value: MersenneField<F>) -> Self {
        Self {
            x: MersenneComplex::from_base(cs, value),
            y: MersenneComplex::zero(cs),
        }
    }

    pub fn from_complex<CS: ConstraintSystem<F>>(cs: &mut CS, value: MersenneComplex<F>) -> Self {
        Self {
            x: value,
            y: MersenneComplex::zero(cs),
        }
    }

    pub fn add<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other: &Self) -> Self {
        Self {
            x: self.x.add(cs, &other.x),
            y: self.y.add(cs, &other.y),
        }
    }

    pub fn double<CS: ConstraintSystem<F>>(&self, cs: &mut CS) -> Self {
        Self {
            x: self.x.double(cs),
            y: self.y.double(cs),
        }
    }

    pub fn sub<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other: &Self) -> Self {
        Self {
            x: self.x.sub(cs, &other.x),
            y: self.y.sub(cs, &other.y),
        }
    }

    pub fn negated<CS: ConstraintSystem<F>>(&self, cs: &mut CS) -> Self {
        Self {
            x: self.x.negated(cs),
            y: self.y.negated(cs),
        }
    }

    pub fn mul<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other: &Self) -> Self {
        // (a + bj)(c + dj) = (ac + kbd) + (ad + bc)j
        let kbd = self.y.mul(cs, &other.y).mul_by_non_residue(cs);
        let bc = self.y.mul(cs, &other.x);

        Self {
            x: self.x.mul_and_add(cs, &other.x, &kbd),
            y: self.x.mul_and_add(cs, &other.y, &bc),
        }
    }

    pub fn mul_and_add<CS: ConstraintSystem<F>>(
        &self,
        cs: &mut CS,
        other_mul: &Self,
        other_add: &Self,
    ) -> Self {
        // (a + bj)(c + dj) + (e + fj) = (ac + kbd + e) + (ad + bc + f)j
        let kbd_plus_e = self
            .y
            .mul_by_non_residue(cs)
            .mul_and_add(cs, &other_mul.y, &other_add.x);
        let bc_plus_f = self.y.mul_and_add(cs, &other_mul.x, &other_add.y);

        Self {
            x: self.x.mul_and_add(cs, &other_mul.x, &kbd_plus_e),
            y: self.x.mul_and_add(cs, &other_mul.y, &bc_plus_f),
        }
    }

    pub fn mul_optimized<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other: &Self) -> Self {
        // (a, b, c, d) = (a1, b1, c1, d1)(a2, b2, c2, d2)
        // (a, b) = (a1, b1)(a2, b2) + (c1, d1)(c2, d2)(2, 1)
        // (c, d) = (a1, b1)(c2, d2) + (c1, d1)(a2, b2)
        // a = a1a2 - b1b2 + 2(c1c2 - d1d2) - (c1d2 + d1c2)
        // b = a1b2 + b1a2 + 2(c1d2 + d1c2) + (c1c2 - d1d2)
        // c = a1c2 - b1d2 + c1a2 - d1b2
        // d = a1d2 + b1c2 + c1b2 + d1a2

        let [a1, b1, c1, d1] = self.into_coeffs();
        let [a2, b2, c2, d2] = other.into_coeffs();

        let c1c2_minus_d1d2 = c1.two_mul_and_sub(cs, &c2, &d1, &d2);
        let c1d2_plus_d1c2 = c1.two_mul_and_add(cs, &d2, &d1, &c2);

        let one = cs.allocate_constant(F::ONE);

        // Computing a
        // (1) a1a2 + modulus^2 + modulus = tmp1
        // (2) tmp1 - b1b2 = tmp2
        // (3) tmp2 + 2(c1c2 - d1d2) = tmp3
        // (4) tmp3 - (c1d2 + d1c2) = tmp4
        // (5) tmp4 - reduce_a * modulus = a
        // reduce_a has 33 bits
        // a has 31 bits
        let tmp1 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (a1.variable, a2.variable),
            F::from_u64_unchecked(M31_MODULUS * (M31_MODULUS + 1)),
            one,
        );
        let tmp2 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::MINUS_ONE,
            (b1.variable, b2.variable),
            F::ONE,
            tmp1,
        );
        let tmp3 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::TWO,
            (c1c2_minus_d1d2.variable, one),
            F::ONE,
            tmp2,
        );
        let tmp4 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::MINUS_ONE,
            (c1d2_plus_d1c2.variable, one),
            F::ONE,
            tmp3,
        );

        let (a, reduce_a) = reduce_mersenne31(cs, tmp4);
        range_check_33_bits(cs, reduce_a);

        // Computing b
        // (1) a1b2 + (c1c2 - d1d2) = tmp6
        // (2) tmp6 + b1a2 = tmp7
        // (3) tmp7 + 2(c1d2 + d1c2) = tmp8
        // (4) tmp8 - reduce_b * modulus = b
        // reduce_b has 33 bits
        // b has 31 bits
        let tmp6 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (a1.variable, b2.variable),
            F::ONE,
            c1c2_minus_d1d2.variable,
        );
        let tmp7 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (b1.variable, a2.variable),
            F::ONE,
            tmp6,
        );
        let tmp8 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::TWO,
            (c1d2_plus_d1c2.variable, one),
            F::ONE,
            tmp7,
        );

        let (b, reduce_b) = reduce_mersenne31(cs, tmp8);
        range_check_33_bits(cs, reduce_b);

        // Computing c
        // (1) a1c2 + 2*modulus^2 = tmp10
        // (2) tmp10 - b1d2 = tmp11
        // (3) tmp11 + c1a2 = tmp12
        // (4) tmp12 - d1b2 = tmp13
        // (6) tmp13 - reduce_c * modulus = c
        // reduce_c has 33 bits
        // c has 31 bits
        let tmp10 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (a1.variable, c2.variable),
            F::from_u64_unchecked(2 * M31_MODULUS * M31_MODULUS),
            one,
        );
        let tmp11 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::MINUS_ONE,
            (b1.variable, d2.variable),
            F::ONE,
            tmp10,
        );
        let tmp12 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (c1.variable, a2.variable),
            F::ONE,
            tmp11,
        );
        let tmp13 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::MINUS_ONE,
            (d1.variable, b2.variable),
            F::ONE,
            tmp12,
        );

        let (c, reduce_c) = reduce_mersenne31(cs, tmp13);
        range_check_33_bits(cs, reduce_c);

        // Computing d
        // (1) a1d2 = tmp15
        // (2) tmp15 + b1c2 = tmp16
        // (3) tmp16 + c1b2 = tmp17
        // (4) tmp17 + d1a2 = tmp18
        // (5) tmp18 - reduce_d * modulus = d
        // reduce_d has 33 bits
        // d has 31 bits
        let tmp15 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (a1.variable, d2.variable),
            F::ZERO,
            one,
        );
        let tmp16 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (b1.variable, c2.variable),
            F::ONE,
            tmp15,
        );
        let tmp17 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (c1.variable, b2.variable),
            F::ONE,
            tmp16,
        );
        let tmp18 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (d1.variable, a2.variable),
            F::ONE,
            tmp17,
        );

        let (d, reduce_d) = reduce_mersenne31(cs, tmp18);
        range_check_33_bits(cs, reduce_d);

        Self {
            x: MersenneComplex { x: a, y: b },
            y: MersenneComplex { x: c, y: d },
        }
    }

    pub fn mul_and_add_optimized<CS: ConstraintSystem<F>>(
        &self,
        cs: &mut CS,
        other_mul: &Self,
        other_add: &Self,
    ) -> Self {
        // (a, b, c, d) = (a1, b1, c1, d1)(a2, b2, c2, d2) + (a3, b3, c3, d3)
        // (a, b) = (a1, b1)(a2, b2) + (c1, d1)(c2, d2)(2, 1) + (a3, b3)
        // (c, d) = (a1, b1)(c2, d2) + (c1, d1)(a2, b2) + (c3, d3)
        // a = a1a2 - b1b2 + 2(c1c2 - d1d2) - (c1d2 + d1c2) + a3
        // b = a1b2 + b1a2 + 2(c1d2 + d1c2) + (c1c2 - d1d2) + b3
        // c = a1c2 - b1d2 + c1a2 - d1b2 + c3
        // d = a1d2 + b1c2 + c1b2 + d1a2 + d3

        let [a1, b1, c1, d1] = self.into_coeffs();
        let [a2, b2, c2, d2] = other_mul.into_coeffs();
        let [a3, b3, c3, d3] = other_add.into_coeffs();

        let c1c2_minus_d1d2 = c1.two_mul_and_sub(cs, &c2, &d1, &d2);
        let c1d2_plus_d1c2 = c1.two_mul_and_add(cs, &d2, &d1, &c2);

        let one = cs.allocate_constant(F::ONE);

        // Computing a
        // (1) a1a2 + a3 = tmp1
        // (2) tmp1 - b1b2 = tmp2
        // (3) tmp2 + 2(c1c2 - d1d2) = tmp3
        // (4) tmp3 - (c1d2 + d1c2) = tmp4
        // (5) tmp4 + modulus^2 + modulus = tmp5
        // (6) tmp5 - reduce_a * modulus = a
        // reduce_a has 33 bits
        // a has 31 bits
        let tmp1 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (a1.variable, a2.variable),
            F::ONE,
            a3.variable,
        );
        let tmp2 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::MINUS_ONE,
            (b1.variable, b2.variable),
            F::ONE,
            tmp1,
        );
        let tmp3 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::TWO,
            (c1c2_minus_d1d2.variable, one),
            F::ONE,
            tmp2,
        );
        let tmp4 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::MINUS_ONE,
            (c1d2_plus_d1c2.variable, one),
            F::ONE,
            tmp3,
        );
        let tmp5 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::from_u64_unchecked(M31_MODULUS * (M31_MODULUS + 1)),
            (one, one),
            F::ONE,
            tmp4,
        );

        let (a, reduce_a) = reduce_mersenne31(cs, tmp5);
        range_check_33_bits(cs, reduce_a);

        // Computing b
        // (1) a1b2 + b3 = tmp6
        // (2) tmp6 + b1a2 = tmp7
        // (3) tmp7 + 2(c1d2 + d1c2) = tmp8
        // (4) tmp8 + (c1c2 - d1d2) = tmp9
        // (5) tmp9 - reduce_b * modulus = b
        // reduce_b has 33 bits
        // b has 31 bits
        let tmp6 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (a1.variable, b2.variable),
            F::ONE,
            b3.variable,
        );
        let tmp7 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (b1.variable, a2.variable),
            F::ONE,
            tmp6,
        );
        let tmp8 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::TWO,
            (c1d2_plus_d1c2.variable, one),
            F::ONE,
            tmp7,
        );
        let tmp9 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (c1c2_minus_d1d2.variable, one),
            F::ONE,
            tmp8,
        );

        let (b, reduce_b) = reduce_mersenne31(cs, tmp9);
        range_check_33_bits(cs, reduce_b);

        // Computing c
        // (1) a1c2 + c3 = tmp10
        // (2) tmp10 - b1d2 = tmp11
        // (3) tmp11 + c1a2 = tmp12
        // (4) tmp12 - d1b2 = tmp13
        // (5) tmp13 + 2*modulus^2 = tmp14
        // (6) tmp14 - reduce_c * modulus = c
        // reduce_c has 33 bits
        // c has 31 bits
        let tmp10 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (a1.variable, c2.variable),
            F::ONE,
            c3.variable,
        );
        let tmp11 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::MINUS_ONE,
            (b1.variable, d2.variable),
            F::ONE,
            tmp10,
        );
        let tmp12 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (c1.variable, a2.variable),
            F::ONE,
            tmp11,
        );
        let tmp13 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::MINUS_ONE,
            (d1.variable, b2.variable),
            F::ONE,
            tmp12,
        );
        let tmp14 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::from_u64_unchecked(2 * M31_MODULUS * M31_MODULUS),
            (one, one),
            F::ONE,
            tmp13,
        );

        let (c, reduce_c) = reduce_mersenne31(cs, tmp14);
        range_check_33_bits(cs, reduce_c);

        // Computing d
        // (1) a1d2 + d3 = tmp15
        // (2) tmp15 + b1c2 = tmp16
        // (3) tmp16 + c1b2 = tmp17
        // (4) tmp17 + d1a2 = tmp18
        // (5) tmp18 - reduce_d * modulus = d
        // reduce_d has 33 bits
        // d has 31 bits
        let tmp15 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (a1.variable, d2.variable),
            F::ONE,
            d3.variable,
        );
        let tmp16 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (b1.variable, c2.variable),
            F::ONE,
            tmp15,
        );
        let tmp17 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (c1.variable, b2.variable),
            F::ONE,
            tmp16,
        );
        let tmp18 = FmaGateInBaseFieldWithoutConstant::compute_fma(
            cs,
            F::ONE,
            (d1.variable, a2.variable),
            F::ONE,
            tmp17,
        );

        let (d, reduce_d) = reduce_mersenne31(cs, tmp18);
        range_check_33_bits(cs, reduce_d);

        Self {
            x: MersenneComplex { x: a, y: b },
            y: MersenneComplex { x: c, y: d },
        }
    }

    // pub fn mul_by_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other: &MersenneField<F>) -> Self {
    //     Self {
    //         x: self.x.mul_by_base(cs, other),
    //         y: self.y.mul_by_base(cs, other),
    //     }
    // }

    // pub fn mul_by_2nd_ext<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other: &MersenneComplex<F>) -> Self {
    //     Self {
    //         x: self.x.mul(cs, other),
    //         y: self.y.mul(cs, other),
    //     }
    // }

    // pub fn mul_by_base_and_add<CS: ConstraintSystem<F>>(&self, cs: &mut CS, coeff: &MersenneField<F>, other: &Self) -> Self {
    //     Self {
    //         x: self.x.mul_by_base_and_add(cs, coeff, &other.x),
    //         y: self.y.mul_by_base_and_add(cs, coeff, &other.y),
    //     }
    // }

    pub fn square<CS: ConstraintSystem<F>>(&self, cs: &mut CS) -> Self {
        // (a, b, c, d) = (a1, b1, c1, d1)^2
        // (a, b) = (a1, b1)^2 + (c1, d1)^2(2, 1)
        // (c, d) = 2(a1, b1)(c1, d1)
        // a = a1^2 - b1^2 + 2(c1^2 - d1^2) - 2c1d1
        // b = 2a1b2 + 4c1d1 + (c1^2 - d1^2)
        // c = 2a1c1 - 2b1d1
        // d = 2a1d1 + 2b1c1
        // TODO: optimize
        self.mul(cs, self)
    }

    pub fn exp_power_of_2<CS: ConstraintSystem<F>>(&self, cs: &mut CS, power_log: usize) -> Self {
        let mut result = self.clone();
        for _ in 0..power_log {
            result = result.square(cs);
        }
        result
    }

    pub fn pow_const<CS: ConstraintSystem<F>>(&self, cs: &mut CS, mut power: usize) -> Self {
        if power == 0 {
            return Self::one(cs);
        }

        let mut bits = vec![];
        while power > 0 {
            bits.push(power & 1);
            power >>= 1;
        }

        let mut result = self.clone();

        for bit in bits.into_iter().rev().skip(1) {
            result = result.square(cs);
            if bit == 1 {
                result = result.mul(cs, self);
            }
        }

        result
    }

    pub fn pow<CS: ConstraintSystem<F>>(&self, cs: &mut CS, power_bits: &[Boolean<F>]) -> Self {
        let one = Self::one(cs);
        let mut result = Self::conditionally_select(cs, *power_bits.last().unwrap(), &self, &one);

        for bit in power_bits.iter().rev().skip(1) {
            result = result.square(cs);

            let res_mul = result.mul(cs, &self);
            result = Self::conditionally_select(cs, *bit, &res_mul, &result);
        }

        result
    }

    /// Computes the division of the value by the other value or zero if the other value is zero
    pub fn div<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other: &Self) -> Self {
        let other_inv = other.inverse_or_zero(cs);
        self.mul(cs, &other_inv)
    }

    /// Computes the inverse of the value or zero if the value is zero
    pub fn inverse_or_zero<CS: ConstraintSystem<F>>(&self, cs: &mut CS) -> Self {
        // TODO: optimize
        let x_square = self.x.square(cs);
        let ky_square = self.y.square(cs).mul_by_non_residue(cs);
        let x_square_minus_ky_square = x_square.sub(cs, &ky_square);
        let x_square_minus_ky_square_inv = x_square_minus_ky_square.inverse_or_zero(cs);
        Self {
            x: self.x.mul(cs, &x_square_minus_ky_square_inv),
            y: self.y.mul(cs, &x_square_minus_ky_square_inv).negated(cs),
        }
    }

    pub fn is_zero<CS: ConstraintSystem<F>>(&mut self, cs: &mut CS) -> Boolean<F> {
        // Could be optimized with concatenation
        let x_is_zero = self.x.is_zero(cs);
        let y_is_zero = self.y.is_zero(cs);
        x_is_zero.and(cs, y_is_zero)
    }

    pub fn equals<CS: ConstraintSystem<F>>(&mut self, cs: &mut CS, other: &mut Self) -> Boolean<F> {
        // Could be optimized with concatenation
        let x_equals = self.x.equals(cs, &mut other.x);
        let y_equals = self.y.equals(cs, &mut other.y);
        x_equals.and(cs, y_equals)
    }

    pub fn enforce_equal<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other: &Self) {
        self.x.enforce_equal(cs, &other.x);
        self.y.enforce_equal(cs, &other.y);
    }

    pub fn mask<CS: ConstraintSystem<F>>(&self, cs: &mut CS, masking_bit: Boolean<F>) -> Self {
        Self {
            x: self.x.mask(cs, masking_bit),
            y: self.y.mask(cs, masking_bit),
        }
    }

    pub fn mask_negated<CS: ConstraintSystem<F>>(
        &self,
        cs: &mut CS,
        masking_bit: Boolean<F>,
    ) -> Self {
        Self {
            x: self.x.mask_negated(cs, masking_bit),
            y: self.y.mask_negated(cs, masking_bit),
        }
    }
}

fn range_check_33_bits<F: SmallField, CS: ConstraintSystem<F>>(cs: &mut CS, variable: Variable) {
    use crate::gadgets::impls::limbs_decompose::decompose_into_limbs;
    use crate::gadgets::non_native_field::implementations::get_16_bits_range_check_table;
    use crate::gadgets::u8::get_8_by_8_range_check_table;

    if let Some(table_id) = get_16_bits_range_check_table(&*cs) {
        let [limb0, limb1, limb2] =
            decompose_into_limbs::<F, CS, 3>(cs, F::from_u64_unchecked(1u64 << 16), variable);

        let zero = cs.allocate_constant(F::ZERO);
        match cs.get_lookup_params().lookup_width() {
            1 => {
                cs.enforce_lookup::<1>(table_id, &[limb0]);
                cs.enforce_lookup::<1>(table_id, &[limb1]);
            }
            3 => {
                cs.enforce_lookup::<3>(table_id, &[limb0, zero, zero]);
                cs.enforce_lookup::<3>(table_id, &[limb1, zero, zero]);
            }
            4 => {
                cs.enforce_lookup::<4>(table_id, &[limb0, zero, zero, zero]);
                cs.enforce_lookup::<4>(table_id, &[limb1, zero, zero, zero]);
            }
            _ => unimplemented!(),
        }
        let _ = Boolean::from_variable_checked(cs, limb2);
    } else if let Some(_table_id) = get_8_by_8_range_check_table(&*cs) {
        let _ = UInt32::from_variable_checked(cs, variable);
    } else {
        unimplemented!()
    }
}

impl<F: SmallField> CSAllocatable<F> for MersenneQuartic<F> {
    type Witness = Mersenne31Quartic;

    fn placeholder_witness() -> Self::Witness {
        Mersenne31Quartic::ZERO
    }
    fn allocate_without_value<CS: ConstraintSystem<F>>(cs: &mut CS) -> Self {
        Self::allocate_checked_without_value(cs, true)
    }
    fn allocate<CS: ConstraintSystem<F>>(cs: &mut CS, witness: Self::Witness) -> Self {
        Self::allocate_checked(cs, witness, false)
    }
    fn allocate_constant<CS: ConstraintSystem<F>>(cs: &mut CS, witness: Self::Witness) -> Self {
        Self {
            x: MersenneComplex::allocate_constant(cs, witness.c0),
            y: MersenneComplex::allocate_constant(cs, witness.c1),
        }
    }
}

impl<F: SmallField> WitnessCastable<F, [F; 4]> for Mersenne31Quartic {
    fn cast_from_source(witness: [F; 4]) -> Self {
        Mersenne31Quartic {
            c0: Mersenne31Complex::cast_from_source(witness[0..2].try_into().unwrap()),
            c1: Mersenne31Complex::cast_from_source(witness[2..4].try_into().unwrap()),
        }
    }

    fn cast_into_source(self) -> [F; 4] {
        [
            F::from_u64_unchecked(self.c0.c0.to_reduced_u32() as u64),
            F::from_u64_unchecked(self.c0.c1.to_reduced_u32() as u64),
            F::from_u64_unchecked(self.c1.c0.to_reduced_u32() as u64),
            F::from_u64_unchecked(self.c1.c1.to_reduced_u32() as u64),
        ]
    }
}

impl<F: SmallField> CSWitnessable<F, 4> for MersenneQuartic<F> {
    type ConversionFunction = Convertor<F, [F; 4], Mersenne31Quartic>;

    fn witness_from_set_of_values(values: [F; 4]) -> Self::Witness {
        <Mersenne31Quartic as WitnessCastable<F, [F; 4]>>::cast_from_source(values)
    }

    fn as_variables_set(&self) -> [Variable; 4] {
        [
            self.x.x.variable,
            self.x.y.variable,
            self.y.x.variable,
            self.y.y.variable,
        ]
    }
}

impl<F: SmallField> WitnessHookable<F> for MersenneQuartic<F> {
    fn witness_hook<CS: ConstraintSystem<F>>(
        &self,
        cs: &CS,
    ) -> Box<dyn FnOnce() -> Option<Self::Witness>> {
        let raw_witness = self.get_witness(cs);
        Box::new(move || raw_witness.wait())
    }
}

impl<F: SmallField> Selectable<F> for MersenneQuartic<F> {
    #[must_use]
    fn conditionally_select<CS: ConstraintSystem<F>>(
        cs: &mut CS,
        flag: Boolean<F>,
        a: &Self,
        b: &Self,
    ) -> Self {
        Self {
            x: MersenneComplex::conditionally_select(cs, flag, &a.x, &b.x),
            y: MersenneComplex::conditionally_select(cs, flag, &a.y, &b.y),
        }
    }
    // const SUPPORTS_PARALLEL_SELECT: bool = true;

    // #[must_use]
    // fn parallel_select<CS: ConstraintSystem<F>, const N: usize>(
    //     cs: &mut CS,
    //     flag: Boolean<F>,
    //     a: &[Self; N],
    //     b: &[Self; N],
    // ) -> [Self; N] {
    //     let ax_nums = a.map(|el| Num::from_variable(el.x.variable));
    //     let ay_nums = a.map(|el| Num::from_variable(el.y.variable));
    //     let bx_nums = b.map(|el| Num::from_variable(el.x.variable));
    //     let by_nums = b.map(|el| Num::from_variable(el.y.variable));

    //     let tmp_x = Num::parallel_select(cs, flag, &ax_nums, &bx_nums);
    //     let tmp_y = Num::parallel_select(cs, flag, &ay_nums, &by_nums);

    //     let mut res = [Self::zero(cs); N];

    //     for i in 0..N {
    //         res[i].x.variable = tmp_x[i].variable;
    //         res[i].y.variable = tmp_y[i].variable;

    //         if a[i].x.reduced && b[i].x.reduced {
    //             res[i].x.reduced = true;
    //         }

    //         if a[i].y.reduced && b[i].y.reduced {
    //             res[i].y.reduced = true;
    //         }
    //     }

    //     res
    // }
}

#[cfg(test)]
mod tests {
    use std::alloc::Global;

    use super::*;
    use crate::cs::*;

    use crate::cs::gates::*;
    use crate::cs::traits::gate::GatePlacementStrategy;
    use crate::dag::CircuitResolverOpts;
    use crate::field::goldilocks::GoldilocksField;
    use crate::gadgets::tables::range_check_16_bits::{
        create_range_check_15_bits_table, create_range_check_16_bits_table, RangeCheck15BitsTable,
        RangeCheck16BitsTable,
    };
    use crate::gadgets::traits::witnessable::WitnessHookable;
    use crate::worker::Worker;
    use mersenne_field::FieldExtension;

    type F = GoldilocksField;

    #[test]
    fn test_mersenne_quartic_field() {
        let geometry = CSGeometry {
            num_columns_under_copy_permutation: 60,
            num_witness_columns: 0,
            num_constant_columns: 4,
            max_allowed_constraint_degree: 4,
        };

        use crate::config::DevCSConfig;
        type RCfg = <DevCSConfig as CSConfig>::ResolverConfig;
        use crate::cs::cs_builder_reference::*;
        let builder_impl =
            CsReferenceImplementationBuilder::<F, F, DevCSConfig>::new(geometry, 1 << 18);
        use crate::cs::cs_builder::new_builder;
        let builder = new_builder::<_, F>(builder_impl);

        let builder = builder.allow_lookup(
            crate::cs::LookupParameters::UseSpecializedColumnsWithTableIdAsConstant {
                width: 1,
                num_repetitions: 10,
                share_table_id: true,
            },
        );

        let builder = ConstantsAllocatorGate::configure_builder(
            builder,
            GatePlacementStrategy::UseGeneralPurposeColumns,
        );
        let builder = FmaGateInBaseFieldWithoutConstant::configure_builder(
            builder,
            GatePlacementStrategy::UseGeneralPurposeColumns,
        );
        let builder = ReductionGate::<F, 4>::configure_builder(
            builder,
            GatePlacementStrategy::UseGeneralPurposeColumns,
        );
        let builder = DotProductGate::<4>::configure_builder(
            builder,
            GatePlacementStrategy::UseGeneralPurposeColumns,
        );
        let builder = UIntXAddGate::<16>::configure_builder(
            builder,
            GatePlacementStrategy::UseGeneralPurposeColumns,
        );
        let builder = SelectionGate::configure_builder(
            builder,
            GatePlacementStrategy::UseGeneralPurposeColumns,
        );
        let builder =
            NopGate::configure_builder(builder, GatePlacementStrategy::UseGeneralPurposeColumns);

        let builder = ReductionGate::<F, 2>::configure_builder(
            builder,
            GatePlacementStrategy::UseGeneralPurposeColumns,
        );
        let builder = ReductionGate::<F, 3>::configure_builder(
            builder,
            GatePlacementStrategy::UseGeneralPurposeColumns,
        );

        let mut owned_cs = builder.build(CircuitResolverOpts::new(1 << 20));

        // add tables
        let table = create_range_check_16_bits_table();
        owned_cs.add_lookup_table::<RangeCheck16BitsTable<1>, 1>(table);

        let table = create_range_check_15_bits_table();
        owned_cs.add_lookup_table::<RangeCheck15BitsTable<1>, 1>(table);

        let cs = &mut owned_cs;

        let rand_base_witness =
            [0; 2].map(|_| Mersenne31Field::new(rand::random::<u32>() % M31_MODULUS as u32));
        let rand_base_vars =
            rand_base_witness.map(|w| MersenneField::<F>::allocate_checked(cs, w, false));

        let rand_witness = [0; 3].map(|_| Mersenne31Quartic {
            c0: Mersenne31Complex {
                c0: Mersenne31Field::new(rand::random::<u32>() % M31_MODULUS as u32),
                c1: Mersenne31Field::new(rand::random::<u32>() % M31_MODULUS as u32),
            },
            c1: Mersenne31Complex {
                c0: Mersenne31Field::new(rand::random::<u32>() % M31_MODULUS as u32),
                c1: Mersenne31Field::new(rand::random::<u32>() % M31_MODULUS as u32),
            },
        });
        let mut rand_vars =
            rand_witness.map(|w| MersenneQuartic::<F>::allocate_checked(cs, w, false));

        // enforce reduced
        for var in rand_vars.iter_mut() {
            var.enforce_reduced(cs);
        }

        // add
        let mut res_witness = rand_witness[0];
        res_witness.add_assign(&rand_witness[1]);
        let res_var = rand_vars[0].add(cs, &rand_vars[1]);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // double
        let mut res_witness = rand_witness[0];
        res_witness.double();
        let res_var = rand_vars[0].double(cs);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // negated
        let mut res_witness = rand_witness[0];
        res_witness.negate();
        let res_var = rand_vars[0].negated(cs);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // sub
        let mut res_witness = rand_witness[0];
        res_witness.sub_assign(&rand_witness[1]);
        let res_var = rand_vars[0].sub(cs, &rand_vars[1]);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // mul
        let mut res_witness = rand_witness[0];
        res_witness.mul_assign(&rand_witness[1]);
        let res_var = rand_vars[0].mul(cs, &rand_vars[1]);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // mul_and_add
        let mut res_witness = rand_witness[0];
        res_witness.mul_assign(&rand_witness[1]);
        res_witness.add_assign(&rand_witness[2]);
        let res_var = rand_vars[0].mul_and_add(cs, &rand_vars[1], &rand_vars[2]);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // square
        let mut res_witness = rand_witness[0];
        res_witness.square();
        let res_var = rand_vars[0].square(cs);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // mul_optimized
        let mut res_witness = rand_witness[0];
        res_witness.mul_assign(&rand_witness[1]);
        let res_var = rand_vars[0].mul_optimized(cs, &rand_vars[1]);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // mul_and_add_optimized
        let mut res_witness = rand_witness[0];
        res_witness.mul_assign(&rand_witness[1]);
        res_witness.add_assign(&rand_witness[2]);
        let res_var = rand_vars[0].mul_and_add_optimized(cs, &rand_vars[1], &rand_vars[2]);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // pow_const
        let rand_power = rand::random::<u32>();
        let res_witness = rand_witness[0].pow(rand_power);
        let res_var = rand_vars[0].pow_const(cs, rand_power as usize);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // pow
        let rand_power_bits: Vec<_> = (0..32)
            .map(|i| Boolean::allocate(cs, (rand_power >> i) & 1 == 1))
            .collect();
        let res_var = rand_vars[0].pow(cs, &rand_power_bits);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // div
        let mut res_witness = rand_witness[0];
        res_witness.mul_assign(&rand_witness[1].inverse().unwrap_or(Mersenne31Quartic::ZERO));
        let res_var = rand_vars[0].div(cs, &rand_vars[1]);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // inverse_or_zero
        let mut res_witness = rand_witness[0];
        res_witness = res_witness.inverse().unwrap_or(Mersenne31Quartic::ZERO);
        let res_var = rand_vars[0].inverse_or_zero(cs);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        let worker = Worker::new_with_num_threads(8);

        drop(cs);
        owned_cs.pad_and_shrink();
        let mut owned_cs = owned_cs.into_assembly::<Global>();
        assert!(owned_cs.check_if_satisfied(&worker));
    }
}
