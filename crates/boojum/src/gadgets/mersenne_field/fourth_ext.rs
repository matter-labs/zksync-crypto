use mersenne_field::{Mersenne31Complex, Mersenne31Quartic};

use super::*;
use super::second_ext::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct MersenneQuartic<F: SmallField> {
    pub x: MersenneComplex<F>,
    pub y: MersenneComplex<F>,
}

impl<F: SmallField> MersenneQuartic<F> {
    pub fn allocated_constant<CS: ConstraintSystem<F>>(cs: &mut CS, value: Mersenne31Quartic) -> Self {
        Self {
            x: MersenneComplex::allocated_constant(cs, value.c0),
            y: MersenneComplex::allocated_constant(cs, value.c1),
        }
    }

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
        [self.x.x.variable, self.x.y.variable, self.y.x.variable, self.y.y.variable]
    }

    pub fn into_nums(&self) -> [Num<F>; 4] {
        [self.x.x.into_num(), self.x.y.into_num(), self.y.x.into_num(), self.y.y.into_num()]
    }

    /// The coordinate values should be in range [0, 2^31 - 2]
    fn from_variables_checked<CS: ConstraintSystem<F>>(cs: &mut CS, variables: [Variable; 4], reduced: bool) -> Self {
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

    pub fn from_base<CS: ConstraintSystem<F>>(cs: &mut CS, value: MersenneFiled<F>) -> Self {
        Self {
            x: MersenneComplex::from_base(cs, value),
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
        // (a + bi)(c + di) = (ac + kbd) + (ad + bc)i
        let ac = self.x.mul(cs, &other.x);
        let kbd = self.y.mul(cs, &other.y).mul_by_non_residue(cs);
        let ad = self.x.mul(cs, &other.y);
        let bc = self.y.mul(cs, &other.x);

        Self {
            x: ac.add(cs, &kbd),
            y: ad.add(cs, &bc),
        }
    }

    pub fn mul_by_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other: &MersenneFiled<F>) -> Self {
        Self {
            x: self.x.mul_by_base(cs, other),
            y: self.y.mul_by_base(cs, other),
        }
    }

    pub fn square<CS: ConstraintSystem<F>>(&self, cs: &mut CS) -> Self {
        // TODO: optimize
        self.mul(cs, self)
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

    pub fn mask<CS: ConstraintSystem<F>>(&self, cs: &mut CS, masking_bit: Boolean<F>) -> Self {
        Self {
            x: self.x.mask(cs, masking_bit),
            y: self.y.mask(cs, masking_bit),
        }
    }

    pub fn mask_negated<CS: ConstraintSystem<F>>(&self, cs: &mut CS, masking_bit: Boolean<F>) -> Self {
        Self {
            x: self.x.mask_negated(cs, masking_bit),
            y: self.y.mask_negated(cs, masking_bit),
        }
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
}

impl<F: SmallField> WitnessCastable<F, [F; 4]> for Mersenne31Quartic {
    fn cast_from_source(witness: [F; 4]) -> Self {
        Mersenne31Quartic{
            c0: Mersenne31Complex::cast_from_source(witness[0..2].try_into().unwrap()),
            c1: Mersenne31Complex::cast_from_source(witness[2..4].try_into().unwrap()),
        }
    }

    fn cast_into_source(self) -> [F; 4] {
        [
            F::from_u64_unchecked(self.c0.c0.to_reduced_u32() as u64),
            F::from_u64_unchecked(self.c0.c1.to_reduced_u32() as u64),
            F::from_u64_unchecked(self.c1.c0.to_reduced_u32() as u64),
            F::from_u64_unchecked(self.c1.c1.to_reduced_u32() as u64)
        ]
    }
}

impl<F: SmallField> CSWitnessable<F, 4> for MersenneQuartic<F> {
    type ConversionFunction = Convertor<F, [F; 4], Mersenne31Quartic>;

    fn witness_from_set_of_values(values: [F; 4]) -> Self::Witness {
        <Mersenne31Quartic as WitnessCastable<F, [F; 4]>>::cast_from_source(values)
    }

    fn as_variables_set(&self) -> [Variable; 4] {
        [self.x.x.variable, self.x.y.variable, self.y.x.variable, self.y.y.variable]
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

    use mersenne_field::FieldExtension;
    use crate::cs::gates::*;
    use crate::cs::traits::gate::GatePlacementStrategy;
    use crate::dag::CircuitResolverOpts;
    use crate::field::goldilocks::GoldilocksField;
    use crate::gadgets::tables::range_check_16_bits::{
        create_range_check_16_bits_table, RangeCheck16BitsTable,
        create_range_check_15_bits_table, RangeCheck15BitsTable,
    };
    use crate::gadgets::traits::witnessable::WitnessHookable;
    use crate::worker::Worker;

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

        let mut owned_cs = builder.build(CircuitResolverOpts::new(1 << 20));

        // add tables
        let table = create_range_check_16_bits_table();
        owned_cs.add_lookup_table::<RangeCheck16BitsTable<1>, 1>(table);

        let table = create_range_check_15_bits_table();
        owned_cs.add_lookup_table::<RangeCheck15BitsTable<1>, 1>(table);

        let cs = &mut owned_cs;

        let rand_base_witness = [0; 2].map(|_| Mersenne31Field::new(rand::random::<u32>() % M31_MODULUS as u32));
        let rand_base_vars = rand_base_witness.map(|w| MersenneFiled::<F>::allocate_checked(cs, w, false));

        let rand_witness = [0; 2].map(|_|
            Mersenne31Quartic {
                c0: Mersenne31Complex {
                    c0: Mersenne31Field::new(rand::random::<u32>() % M31_MODULUS as u32),
                    c1: Mersenne31Field::new(rand::random::<u32>() % M31_MODULUS as u32),
                },
                c1: Mersenne31Complex {
                    c0: Mersenne31Field::new(rand::random::<u32>() % M31_MODULUS as u32),
                    c1: Mersenne31Field::new(rand::random::<u32>() % M31_MODULUS as u32),
                },
            }
        );
        let mut rand_vars = rand_witness.map(|w| MersenneQuartic::<F>::allocate_checked(cs, w, false));

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

        // mul_by_base
        let mut res_witness = rand_witness[0];
        res_witness.mul_assign_by_base(&rand_base_witness[0]);
        let res_var = rand_vars[0].mul_by_base(cs, &rand_base_vars[0]);
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
