use mersenne_field::field::BaseField;
use mersenne_field::Mersenne31Complex;

use super::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct MersenneComplex<F: SmallField> {
    pub x: MersenneField<F>,
    pub y: MersenneField<F>,
}

impl<F: SmallField> MersenneComplex<F> {
    pub fn zero<CS: ConstraintSystem<F>>(cs: &mut CS) -> Self {
        Self {
            x: MersenneField::zero(cs),
            y: MersenneField::zero(cs),
        }
    }

    pub fn one<CS: ConstraintSystem<F>>(cs: &mut CS) -> Self {
        Self {
            x: MersenneField::one(cs),
            y: MersenneField::zero(cs),
        }
    }

    pub fn minus_one<CS: ConstraintSystem<F>>(cs: &mut CS) -> Self {
        Self {
            x: MersenneField::minus_one(cs),
            y: MersenneField::zero(cs),
        }
    }

    /// Returns the quadratic non-residue in the field: 2 + i
    pub fn non_residue<CS: ConstraintSystem<F>>(cs: &mut CS) -> Self {
        Self::allocate_constant(cs, Mersenne31Complex::QUADRATIC_NON_RESIDUE)
    }

    pub fn get_modulus_num<CS: ConstraintSystem<F>>(cs: &mut CS) -> Num<F> {
        Num::from_variable(cs.allocate_constant(F::from_u64_unchecked(M31_MODULUS)))
    }

    pub fn get_variables(&self) -> [Variable; 2] {
        [self.x.variable, self.y.variable]
    }

    pub fn into_nums(&self) -> [Num<F>; 2] {
        [self.x.into_num(), self.y.into_num()]
    }

    pub fn into_uint32s(&self) -> [UInt32<F>; 2] {
        [self.x.into_uint32(), self.y.into_uint32()]
    }

    /// The coordinate values should be in range [0, 2^31 - 2]
    pub fn from_variables_checked<CS: ConstraintSystem<F>>(
        cs: &mut CS,
        variables: [Variable; 2],
        reduced: bool,
    ) -> Self {
        Self {
            x: MersenneField::from_variable_checked(cs, variables[0], reduced),
            y: MersenneField::from_variable_checked(cs, variables[1], reduced),
        }
    }

    /// The coordinate values should be in range [0, 2^31 - 2]
    pub fn allocate_checked_without_value<CS: ConstraintSystem<F>>(
        cs: &mut CS,
        reduced: bool,
    ) -> Self {
        Self {
            x: MersenneField::allocate_checked_without_value(cs, reduced),
            y: MersenneField::allocate_checked_without_value(cs, reduced),
        }
    }

    pub fn allocate_checked<CS: ConstraintSystem<F>>(
        cs: &mut CS,
        witness: Mersenne31Complex,
        reduced: bool,
    ) -> Self {
        Self {
            x: MersenneField::allocate_checked(cs, witness.c0, reduced),
            y: MersenneField::allocate_checked(cs, witness.c1, reduced),
        }
    }

    pub fn enforce_reduced<CS: ConstraintSystem<F>>(&mut self, cs: &mut CS) {
        self.x.enforce_reduced(cs);
        self.y.enforce_reduced(cs);
    }

    pub fn from_coeffs(coefficients: [MersenneField<F>; 2]) -> Self {
        Self {
            x: coefficients[0],
            y: coefficients[1],
        }
    }

    pub fn from_base<CS: ConstraintSystem<F>>(cs: &mut CS, value: MersenneField<F>) -> Self {
        Self {
            x: value,
            y: MersenneField::zero(cs),
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
        Self {
            x: self.x.two_mul_and_sub(cs, &other.x, &other.y, &self.y),
            y: self.x.two_mul_and_add(cs, &other.y, &other.x, &self.y),
        }
    }

    pub fn mul_and_add<CS: ConstraintSystem<F>>(
        &self,
        cs: &mut CS,
        other_mul: &Self,
        other_add: &Self,
    ) -> Self {
        Self {
            x: self.x.two_mul_and_sub_and_add(
                cs,
                &other_mul.x,
                &other_mul.y,
                &self.y,
                &other_add.x,
            ),
            y: self
                .x
                .two_mul_and_two_add(cs, &other_mul.y, &other_mul.x, &self.y, &other_add.y),
        }
    }

    pub fn square<CS: ConstraintSystem<F>>(&self, cs: &mut CS) -> Self {
        Self {
            x: self.x.two_mul_and_sub(cs, &self.x, &self.y, &self.y),
            y: self.x.mul_times_2(cs, &self.y),
        }
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

    /// Computes (x + y * i)(2 + i) = (2x - y) + (2y + x)i
    pub fn mul_by_non_residue<CS: ConstraintSystem<F>>(&self, cs: &mut CS) -> Self {
        // We will use the following system of constraints:
        // (1) 2x + modulus = tmp
        // (2) tmp - y = tmp2
        // (3) tmp2 - reduce1 * modulus = result_x
        // (4) 2y + x = tmp3
        // (5) tmp3 - reduce2 * modulus = result_y
        // (6) reduce1 and reduce2 have 8 bits
        // (7) result_x has 31 bits
        // (8) result_y has 31 bits

        let tmp = Num::allocate_without_value(cs);
        let tmp2 = Num::allocate_without_value(cs);
        let tmp3 = Num::allocate_without_value(cs);
        let reduce1 = Num::allocate_without_value(cs);
        let reduce2 = Num::allocate_without_value(cs);
        crate::gadgets::u8::range_check_u8_pair(
            cs,
            &[reduce1.get_variable(), reduce2.get_variable()],
        ); // 6th constraint
        let result_x = MersenneField::allocate_checked_without_value(cs, false); // 7th constraint
        let result_y = MersenneField::allocate_checked_without_value(cs, false); // 8th constraint

        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 2]| {
                let x = inputs[0].as_u64();
                let y = inputs[1].as_u64();

                let tmp = 2 * x + M31_MODULUS;
                let tmp2 = tmp - y;
                let reduce1 = tmp2 / M31_MODULUS;
                let result_x = tmp2 % M31_MODULUS;
                let tmp3 = 2 * y + x;
                let reduce2 = tmp3 / M31_MODULUS;
                let result_y = tmp3 % M31_MODULUS;

                [
                    F::from_u64_unchecked(tmp),
                    F::from_u64_unchecked(tmp2),
                    F::from_u64_unchecked(tmp3),
                    F::from_u64_unchecked(reduce1),
                    F::from_u64_unchecked(reduce2),
                    F::from_u64_unchecked(result_x),
                    F::from_u64_unchecked(result_y),
                ]
            };

            let dependencies = Place::from_variables([self.x.variable, self.y.variable]);
            let outputs = Place::from_variables([
                tmp.get_variable(),
                tmp2.get_variable(),
                tmp3.get_variable(),
                reduce1.variable,
                reduce2.variable,
                result_x.variable,
                result_y.variable,
            ]);

            cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
        }

        if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
            if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
                // (1) 2x + modulus = tmp
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::TWO,
                    linear_term_coeff: F::from_u64_unchecked(M31_MODULUS),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (self.x.variable, cs.allocate_constant(F::ONE)),
                    linear_part: cs.allocate_constant(F::ONE),
                    rhs_part: tmp.variable,
                };
                gate.add_to_cs(cs);

                // (2) tmp - y = tmp2
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: F::MINUS_ONE,
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (tmp.variable, cs.allocate_constant(F::ONE)),
                    linear_part: self.y.variable,
                    rhs_part: tmp2.variable,
                };
                gate.add_to_cs(cs);

                // (3) tmp2 - reduce1 * modulus = result_x
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: *F::from_u64_unchecked(M31_MODULUS).negate(),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (tmp2.variable, cs.allocate_constant(F::ONE)),
                    linear_part: reduce1.variable,
                    rhs_part: result_x.variable,
                };
                gate.add_to_cs(cs);

                // (4) 2y + x = tmp3
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::TWO,
                    linear_term_coeff: F::ONE,
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (self.y.variable, cs.allocate_constant(F::ONE)),
                    linear_part: self.x.variable,
                    rhs_part: tmp3.variable,
                };
                gate.add_to_cs(cs);

                // (5) tmp3 - reduce2 * modulus = result_y
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: *F::from_u64_unchecked(M31_MODULUS).negate(),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (tmp3.variable, cs.allocate_constant(F::ONE)),
                    linear_part: reduce2.variable,
                    rhs_part: result_y.variable,
                };
                gate.add_to_cs(cs);
            }
        } else {
            unimplemented!()
        }

        Self {
            x: result_x,
            y: result_y,
        }
    }

    /// Computes the division of the value by the other value or zero if the other value is zero
    pub fn div<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other: &Self) -> Self {
        let other_inv = other.inverse_or_zero(cs);
        self.mul(cs, &other_inv)
    }

    /// Computes the inverse of the value or zero if the value is zero
    pub fn inverse_or_zero<CS: ConstraintSystem<F>>(&self, cs: &mut CS) -> Self {
        let x_square_plus_y_square = self.x.two_mul_and_add(cs, &self.x, &self.y, &self.y);
        let x_square_plus_y_square_inv = x_square_plus_y_square.inverse_or_zero(cs);
        Self {
            x: self.x.mul(cs, &x_square_plus_y_square_inv),
            y: self.y.mul_negate(cs, &x_square_plus_y_square_inv),
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

impl<F: SmallField> CSAllocatable<F> for MersenneComplex<F> {
    type Witness = Mersenne31Complex;

    fn placeholder_witness() -> Self::Witness {
        Mersenne31Complex::ZERO
    }
    fn allocate_without_value<CS: ConstraintSystem<F>>(cs: &mut CS) -> Self {
        Self::allocate_checked_without_value(cs, true)
    }
    fn allocate<CS: ConstraintSystem<F>>(cs: &mut CS, witness: Self::Witness) -> Self {
        Self::allocate_checked(cs, witness, false)
    }
    fn allocate_constant<CS: ConstraintSystem<F>>(cs: &mut CS, witness: Self::Witness) -> Self {
        Self {
            x: MersenneField::allocate_constant(cs, witness.c0),
            y: MersenneField::allocate_constant(cs, witness.c1),
        }
    }
}

impl<F: SmallField> WitnessCastable<F, [F; 2]> for Mersenne31Complex {
    fn cast_from_source(witness: [F; 2]) -> Self {
        let value_c0 = witness[0].as_u64();
        let value_c1 = witness[1].as_u64();
        assert!(value_c0 < M31_MODULUS);
        assert!(value_c1 < M31_MODULUS);
        Mersenne31Complex {
            c0: Mersenne31Field::new(value_c0 as u32),
            c1: Mersenne31Field::new(value_c1 as u32),
        }
    }

    fn cast_into_source(self) -> [F; 2] {
        [
            F::from_u64_unchecked(self.c0.to_reduced_u32() as u64),
            F::from_u64_unchecked(self.c1.to_reduced_u32() as u64),
        ]
    }
}

impl<F: SmallField> CSWitnessable<F, 2> for MersenneComplex<F> {
    type ConversionFunction = Convertor<F, [F; 2], Mersenne31Complex>;

    fn witness_from_set_of_values(values: [F; 2]) -> Self::Witness {
        <Mersenne31Complex as WitnessCastable<F, [F; 2]>>::cast_from_source(values)
    }

    fn as_variables_set(&self) -> [Variable; 2] {
        [self.x.variable, self.y.variable]
    }
}

impl<F: SmallField> WitnessHookable<F> for MersenneComplex<F> {
    fn witness_hook<CS: ConstraintSystem<F>>(
        &self,
        cs: &CS,
    ) -> Box<dyn FnOnce() -> Option<Self::Witness>> {
        let raw_witness = self.get_witness(cs);
        Box::new(move || raw_witness.wait())
    }
}

impl<F: SmallField> Selectable<F> for MersenneComplex<F> {
    #[must_use]
    fn conditionally_select<CS: ConstraintSystem<F>>(
        cs: &mut CS,
        flag: Boolean<F>,
        a: &Self,
        b: &Self,
    ) -> Self {
        Self {
            x: MersenneField::conditionally_select(cs, flag, &a.x, &b.x),
            y: MersenneField::conditionally_select(cs, flag, &a.y, &b.y),
        }
    }
    const SUPPORTS_PARALLEL_SELECT: bool = true;

    #[must_use]
    fn parallel_select<CS: ConstraintSystem<F>, const N: usize>(
        cs: &mut CS,
        flag: Boolean<F>,
        a: &[Self; N],
        b: &[Self; N],
    ) -> [Self; N] {
        let ax_nums = a.map(|el| Num::from_variable(el.x.variable));
        let ay_nums = a.map(|el| Num::from_variable(el.y.variable));
        let bx_nums = b.map(|el| Num::from_variable(el.x.variable));
        let by_nums = b.map(|el| Num::from_variable(el.y.variable));

        let tmp_x = Num::parallel_select(cs, flag, &ax_nums, &bx_nums);
        let tmp_y = Num::parallel_select(cs, flag, &ay_nums, &by_nums);

        let mut res = [Self::zero(cs); N];

        for i in 0..N {
            res[i].x.variable = tmp_x[i].variable;
            res[i].y.variable = tmp_y[i].variable;

            if a[i].x.reduced && b[i].x.reduced {
                res[i].x.reduced = true;
            }

            if a[i].y.reduced && b[i].y.reduced {
                res[i].y.reduced = true;
            }
        }

        res
    }
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
    fn test_mersenne_complex_field() {
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

        let rand_witness = [0; 2].map(|_| Mersenne31Complex {
            c0: Mersenne31Field::new(rand::random::<u32>() % M31_MODULUS as u32),
            c1: Mersenne31Field::new(rand::random::<u32>() % M31_MODULUS as u32),
        });
        let mut rand_vars =
            rand_witness.map(|w| MersenneComplex::<F>::allocate_checked(cs, w, false));

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

        // square
        let mut res_witness = rand_witness[0];
        res_witness.square();
        let res_var = rand_vars[0].square(cs);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // mul_by_non_residue
        let mut res_witness = rand_witness[0];
        res_witness.mul_assign(&Mersenne31Complex::QUADRATIC_NON_RESIDUE);
        let res_var = rand_vars[0].mul_by_non_residue(cs);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // div
        let mut res_witness = rand_witness[0];
        res_witness.mul_assign(&rand_witness[1].inverse().unwrap_or(Mersenne31Complex::ZERO));
        let res_var = rand_vars[0].div(cs, &rand_vars[1]);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // inverse_or_zero
        let mut res_witness = rand_witness[0];
        res_witness = res_witness.inverse().unwrap_or(Mersenne31Complex::ZERO);
        let res_var = rand_vars[0].inverse_or_zero(cs);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        let worker = Worker::new_with_num_threads(8);

        drop(cs);
        owned_cs.pad_and_shrink();
        let mut owned_cs = owned_cs.into_assembly::<Global>();
        assert!(owned_cs.check_if_satisfied(&worker));
    }
}
