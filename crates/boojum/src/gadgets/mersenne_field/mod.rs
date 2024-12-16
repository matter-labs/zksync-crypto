use crate::cs::gates::FmaGateInBaseFieldWithoutConstant;
use crate::cs::gates::ReductionGate;
use crate::gadgets::traits::allocatable::CSAllocatable;
use crate::gadgets::SmallField;
use crate::gadgets::num::Num;
use crate::cs::traits::cs::ConstraintSystem;
use crate::cs::Variable;
use crate::cs::gates::ConstantAllocatableCS;
use crate::gadgets::traits::witnessable::CSWitnessable;
use crate::gadgets::traits::castable::WitnessCastable;
use crate::gadgets::traits::castable::Convertor;
use crate::config::CSSetupConfig;
use crate::gadgets::boolean::Boolean;
use crate::gadgets::traits::witnessable::WitnessHookable;
use crate::cs::gates::FmaGateInBaseWithoutConstantParams;
use crate::cs::gates::ReductionGateParams;
use crate::gadgets::u32::UInt32;
use crate::gadgets::traits::selectable::Selectable;
use crate::gadgets::Place;
use crate::config::CSConfig;
use crate::config::CSWitnessEvaluationConfig;

use mersenne_field::Mersenne31Field;
use mersenne_field::field::Field;

pub mod second_ext;
pub mod fourth_ext;

const M31_MODULUS: u64 = (1 << 31) - 1;

// #[derive(Derivative, serde::Serialize, serde::Deserialize)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct MersenneFiled<F: SmallField> {
    // the inner witness value is always reduced with the modulus
    // if reduced is true, then the reduction is proved and max possible value is 2^31 - 2
    // if reduced is false, then only 31-bit check is proved and max possible value is 2^31 - 1
    // NOTE: F should have at least 63 bits capacity for correctness of arithmetical operations
    pub(crate) variable: Variable,
    pub(crate) reduced: bool,
    pub(crate) _marker: std::marker::PhantomData<F>,
}

impl<F: SmallField> MersenneFiled<F> {
    pub fn allocated_constant<CS: ConstraintSystem<F>>(cs: &mut CS, value: Mersenne31Field) -> Self {
        let variable = cs.allocate_constant(F::from_u64_unchecked(value.to_reduced_u32() as u64));

        Self { variable, reduced: true, _marker: std::marker::PhantomData }
    }

    pub fn zero<CS: ConstraintSystem<F>>(cs: &mut CS) -> Self {
        let variable = cs.allocate_constant(F::from_u64_unchecked(0));

        Self { variable, reduced: true, _marker: std::marker::PhantomData }
    }

    pub fn one<CS: ConstraintSystem<F>>(cs: &mut CS) -> Self {
        let variable = cs.allocate_constant(F::from_u64_unchecked(1));

        Self { variable, reduced: true, _marker: std::marker::PhantomData }
    }

    pub fn minus_one<CS: ConstraintSystem<F>>(cs: &mut CS) -> Self {
        let variable = cs.allocate_constant(F::from_u64_unchecked(M31_MODULUS - 1));

        Self { variable, reduced: true, _marker: std::marker::PhantomData }
    }

    pub fn get_modulus_num<CS: ConstraintSystem<F>>(cs: &mut CS) -> Num<F> {
        Num::from_variable(cs.allocate_constant(F::from_u64_unchecked(M31_MODULUS)))
    }

    pub fn get_variable(&self) -> Variable {
        self.variable
    }

    pub fn into_num(&self) -> Num<F> {
        Num::from_variable(self.variable)
    }

    /// The value should be in range [0, 2^31 - 2]
    fn from_variable_checked<CS: ConstraintSystem<F>>(cs: &mut CS, variable: Variable, reduced: bool) -> Self {
        let mut result = Self {
            variable,
            reduced: false,
            _marker: std::marker::PhantomData,
        };

        range_check_31_bits(cs, variable);

        if reduced {
            result.enforce_reduced(cs);
        }

        result
    }

    /// The value should be in range [0, 2^31 - 2]
    fn allocate_checked_without_value<CS: ConstraintSystem<F>>(cs: &mut CS, reduced: bool) -> Self {
        let variable = cs.alloc_variable_without_value();

        Self::from_variable_checked(cs, variable, reduced)
    }

    pub fn allocate_checked<CS: ConstraintSystem<F>>(
        cs: &mut CS,
        witness: Mersenne31Field,
        reduced: bool,
    ) -> Self {
        let variable = cs.alloc_single_variable_from_witness(
            F::from_u64_unchecked(witness.to_reduced_u32() as u64)
        );

        Self::from_variable_checked(cs, variable, reduced)
    }

    pub fn enforce_reduced<CS: ConstraintSystem<F>>(&mut self, cs: &mut CS) {
        // The value is already checked to be in range [0, 2^31 - 1]
        // so we should only check that value != 2^31 - 1
        // We can just use this gate:
        // ((2^31 - 1) - self) * inv = 1

        if self.reduced {
            return;
        }

        let inv = cs.alloc_variable_without_value();
        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 1]| {
                let mut a = F::from_u64_unchecked(M31_MODULUS);
                a.sub_assign(&inputs[0]);

                [a.inverse().unwrap()]
            };

            let dependencies = Place::from_variables([self.variable]);
            let outputs = Place::from_variables([inv]);

            cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
        }

        if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
            if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
                // ((2^31 - 1) - self) * inv = 1
                // -self*inv + (2^31 - 1)*inv = 1
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::MINUS_ONE,
                    linear_term_coeff: F::from_u64_unchecked(M31_MODULUS),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (self.variable, inv),
                    linear_part: inv,
                    rhs_part: cs.allocate_constant(F::ONE),
                };
                gate.add_to_cs(cs);
            }
        } else {
            unimplemented!()
        }

        self.reduced = true;
    }

    pub fn add<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other: &Self) -> Self {
        // We will use the following system of constraints:
        // (1) self + other = tmp
        // (2) tmp - reduce * modulus = result
        // (3) reduce is boolean
        // (4) result has 31 bits

        let tmp = self.into_num().add(cs, &other.into_num()); // 1st constraint
        let reduce = Boolean::allocate_without_value(cs); // 3rd constraint
        let result = Self::allocate_checked_without_value(cs, false); // 4th constraint
        
        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 1]| {
                let sum = inputs[0].as_u64();
                
                let (reduce, result) = if sum >= M31_MODULUS {
                    (F::ONE, sum - M31_MODULUS)
                } else {
                    (F::ZERO, sum)
                };

                [reduce, F::from_u64_unchecked(result)]
            };

            let dependencies = Place::from_variables([tmp.get_variable()]);
            let outputs = Place::from_variables([reduce.get_variable(), result.variable]);

            cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
        }

        if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
            if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
                // (2) tmp + reduce * (-modulus) = result
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: *F::from_u64_unchecked(M31_MODULUS).negate(),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (tmp.variable, cs.allocate_constant(F::ONE)),
                    linear_part: reduce.get_variable(),
                    rhs_part: result.variable,
                };
                gate.add_to_cs(cs);
            }
        } else {
            unimplemented!()
        }

        result
    }

    pub fn double<CS: ConstraintSystem<F>>(&self, cs: &mut CS) -> Self {
        // We will use the following system of constraints:
        // (1) 2 * self - reduce * modulus = result
        // (2) reduce is boolean
        // (3) result has 31 bits

        let reduce = Boolean::allocate_without_value(cs); // 2nd constraint
        let result = Self::allocate_checked_without_value(cs, false); // 3rd constraint
        
        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 1]| {
                let a = inputs[0].as_u64();
                
                let (reduce, result) = if a >= 1<<30 {
                    (F::ONE, 2*a - M31_MODULUS)
                } else {
                    (F::ZERO, 2*a)
                };

                [reduce, F::from_u64_unchecked(result)]
            };

            let dependencies = Place::from_variables([self.variable]);
            let outputs = Place::from_variables([reduce.get_variable(), result.variable]);

            cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
        }

        if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
            if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
                // (1) 2 * self + reduce * (-modulus) = result
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::TWO,
                    linear_term_coeff: *F::from_u64_unchecked(M31_MODULUS).negate(),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (self.variable, cs.allocate_constant(F::ONE)),
                    linear_part: reduce.get_variable(),
                    rhs_part: result.variable,
                };
                gate.add_to_cs(cs);
            }
        } else {
            unimplemented!()
        }

        result
    }

    pub fn negated<CS: ConstraintSystem<F>>(&self, cs: &mut CS) -> Self {
        // We will use the following system of constraints:
        // (1) -self + reduce * modulus = result
        // (2) reduce is boolean
        // (3) result has 31 bits

        let reduce = Boolean::allocate_without_value(cs); // 2nd constraint
        let result = Self::allocate_checked_without_value(cs, false); // 3rd constraint

        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 1]| {
                let a = inputs[0].as_u64();
                
                let (reduce, result) = if a == 0 {
                    (F::ZERO, 0)
                } else {
                    (F::ONE, M31_MODULUS - a)
                };

                [reduce, F::from_u64_unchecked(result)]
            };

            let dependencies = Place::from_variables([self.variable]);
            let outputs = Place::from_variables([reduce.get_variable(), result.variable]);

            cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
        }

        if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
            if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
                // (1) -self + reduce * modulus = result
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::MINUS_ONE,
                    linear_term_coeff: F::from_u64_unchecked(M31_MODULUS),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (self.variable, cs.allocate_constant(F::ONE)),
                    linear_part: reduce.get_variable(),
                    rhs_part: result.variable,
                };
                gate.add_to_cs(cs);
            }
        } else {
            unimplemented!()
        }

        let zero = Self::zero(cs);
        zero.sub(cs, self)
    }

    pub fn sub<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other: &Self) -> Self {
        // We will use the following system of constraints:
        // (1) self - other = tmp
        // (2) tmp + reduce * modulus = result
        // (3) reduce is boolean
        // (4) result has 31 bits

        let tmp = self.into_num().sub(cs, &other.into_num()); // 1st constraint
        let reduce = Boolean::allocate_without_value(cs); // 3rd constraint
        let result = Self::allocate_checked_without_value(cs, false); // 4th constraint
        
        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 2]| {
                let a = inputs[0].as_u64();
                let b = inputs[1].as_u64();
                
                let (reduce, result) = if a < b {
                    (F::ONE, M31_MODULUS + a - b)
                } else {
                    (F::ZERO, a - b)
                };

                [reduce, F::from_u64_unchecked(result)]
            };

            let dependencies = Place::from_variables([self.variable, other.variable]);
            let outputs = Place::from_variables([reduce.get_variable(), result.variable]);

            cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
        }

        if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
            if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
                // (2) tmp + reduce * modulus = result
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: F::from_u64_unchecked(M31_MODULUS),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (tmp.variable, cs.allocate_constant(F::ONE)),
                    linear_part: reduce.get_variable(),
                    rhs_part: result.variable,
                };
                gate.add_to_cs(cs);
            }
        } else {
            unimplemented!()
        }

        result
    }

    pub fn mul<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other: &Self) -> Self {
        // We will use the following system of constraints:
        // (1) self * other = result + reduce * modulus
        // (2) reduce has 32 bits
        // (3) result has 31 bits
        
        let reduce = Num::allocate_without_value(cs);
        range_check_32_bits(cs, reduce.get_variable()); // 2nd constraint
        let result = Self::allocate_checked_without_value(cs, false); // 3rd constraint
        
        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 2]| {
                let a = inputs[0].as_u64();
                let b = inputs[1].as_u64();
                
                let prod = a * b;
                let reduce = prod / M31_MODULUS;
                let result = prod % M31_MODULUS;
                assert!(reduce < 1<<32);

                [F::from_u64_unchecked(reduce), F::from_u64_unchecked(result)]
            };

            let dependencies = Place::from_variables([self.variable, other.variable]);
            let outputs = Place::from_variables([reduce.get_variable(), result.variable]);

            cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
        }

        if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
            if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
                // (1) self * other + reduce * (-modulus) = result
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: *F::from_u64_unchecked(M31_MODULUS).negate(),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (self.variable, other.variable),
                    linear_part: reduce.get_variable(),
                    rhs_part: result.variable,
                };
                gate.add_to_cs(cs);
            }
        } else {
            unimplemented!()
        }

        result
    }

    pub fn square<CS: ConstraintSystem<F>>(&self, cs: &mut CS) -> Self {
        self.mul(cs, self)
    }

    /// Computes - self * other
    pub fn mul_negate<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other: &Self) -> Self {
        // We will use the following system of constraints:
        // (1) reduce * modulus - self * other = result
        // (2) reduce has 32 bits
        // (3) result has 31 bits
        
        let reduce = Num::allocate_without_value(cs);
        range_check_32_bits(cs, reduce.get_variable()); // 2nd constraint
        let result = Self::allocate_checked_without_value(cs, false); // 3rd constraint
        
        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 2]| {
                let a = inputs[0].as_u64();
                let b = inputs[1].as_u64();
                
                let prod = a * b;
                let mut reduce = prod / M31_MODULUS;
                let mut result = prod % M31_MODULUS;
                if result != 0 {
                    reduce += 1;
                    result = M31_MODULUS - result;
                }

                assert!(reduce < 1<<32);

                [F::from_u64_unchecked(reduce), F::from_u64_unchecked(result)]
            };

            let dependencies = Place::from_variables([self.variable, other.variable]);
            let outputs = Place::from_variables([reduce.get_variable(), result.variable]);

            cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
        }

        if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
            if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
                // (1) reduce * modulus - self * other = result
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::MINUS_ONE,
                    linear_term_coeff: F::from_u64_unchecked(M31_MODULUS),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (self.variable, other.variable),
                    linear_part: reduce.get_variable(),
                    rhs_part: result.variable,
                };
                gate.add_to_cs(cs);
            }
        } else {
            unimplemented!()
        }

        result
    }

    /// Computes 2 * self * other
    pub fn mul_times_2<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other: &Self) -> Self {
        // We will use the following system of constraints:
        // (1) 2 * self * other = result + reduce * modulus
        // (2) reduce has 32 bits
        // (3) result has 31 bits
        // NOTE: max value of number is 2^31 - 1, so the max value of 2*a*b is 2*(2^31-1)^2
        // 2*(2^31-1)^2 = 2*(2^62 - 2*2^31 + 1) = 2^63 - 4*2^31 + 2 fits to 63 bits, so no Goldilocks overflow
        // also the max reduce value we can have is 2*(2^31 - 1) = 2^32 - 2 fits to 32 bits
        
        let reduce = Num::allocate_without_value(cs);
        range_check_32_bits(cs, reduce.get_variable()); // 2nd constraint
        let result = Self::allocate_checked_without_value(cs, false); // 3rd constraint
        
        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 2]| {
                let a = inputs[0].as_u64();
                let b = inputs[1].as_u64();
                
                let prod = 2 * a * b;
                let reduce = prod / M31_MODULUS;
                let result = prod % M31_MODULUS;
                assert!(reduce < 1<<32);

                [F::from_u64_unchecked(reduce), F::from_u64_unchecked(result)]
            };

            let dependencies = Place::from_variables([self.variable, other.variable]);
            let outputs = Place::from_variables([reduce.get_variable(), result.variable]);

            cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
        }

        if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
            if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
                // (1) 2 * self * other + reduce * (-modulus) = result
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::TWO,
                    linear_term_coeff: *F::from_u64_unchecked(M31_MODULUS).negate(),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (self.variable, other.variable),
                    linear_part: reduce.get_variable(),
                    rhs_part: result.variable,
                };
                gate.add_to_cs(cs);
            }
        } else {
            unimplemented!()
        }

        result
    }

    /// Computes self * other_mul + other_add
    pub fn mul_and_add<CS: ConstraintSystem<F>>(&self, cs: &mut CS, other_mul: &Self, other_add: &Self) -> Self {
        // We will use the following system of constraints:
        // (1) self * other_mul + other_add = tmp
        // (2) tmp - reduce * modulus = result
        // (3) reduce has 32 bits
        // (4) result has 31 bits
        
        let tmp = Num::allocate_without_value(cs);
        let reduce = Num::allocate_without_value(cs);
        range_check_32_bits(cs, reduce.get_variable()); // 3rd constraint
        let result = Self::allocate_checked_without_value(cs, false); // 4th constraint
        
        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 3]| {
                let a = inputs[0].as_u64();
                let b = inputs[1].as_u64();
                let c = inputs[2].as_u64();
                
                let tmp = a * b + c;
                let reduce = tmp / M31_MODULUS;
                let result = tmp % M31_MODULUS;
                assert!(reduce < 1<<32);

                [F::from_u64_unchecked(tmp), F::from_u64_unchecked(reduce), F::from_u64_unchecked(result)]
            };

            let dependencies = Place::from_variables([self.variable, other_mul.variable, other_add.variable]);
            let outputs = Place::from_variables([tmp.get_variable(), reduce.get_variable(), result.variable]);

            cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
        }

        if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
            if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
                // (1) self * other_mul + other_add = tmp
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: F::ONE,
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (self.variable, other_mul.variable),
                    linear_part: other_add.get_variable(),
                    rhs_part: tmp.get_variable(),
                };
                gate.add_to_cs(cs);

                // (2) tmp - reduce * modulus = result
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: *F::from_u64_unchecked(M31_MODULUS).negate(),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (tmp.get_variable(), cs.allocate_constant(F::ONE)),
                    linear_part: reduce.get_variable(),
                    rhs_part: result.variable,
                };
                gate.add_to_cs(cs);
            }
        } else {
            unimplemented!()
        }

        result
    }

    /// Computes self * other_mul + other_add_1 * other_add_2
    pub fn two_mul_and_add<CS: ConstraintSystem<F>>(
        &self, 
        cs: &mut CS, 
        other_mul: &Self,
        other_add_1: &Self,
        other_add_2: &Self,
    ) -> Self {
        // We will use the following system of constraints:
        // (1) self * other_mul - reduce * modulus = tmp
        // (2) tmp + other_add_1 * other_add_2 = result
        // (3) reduce has 32 bits
        // (4) result has 31 bits
        // Note that max value of number is 2^31 - 1, so the max value of a*b + c*d is 2*(2^31-1)^2
        // 2*(2^31-1)^2 = 2*(2^62 - 2*2^31 + 1) = 2^63 - 4*2^31 + 2 fits to 63 bits, so no Goldilocks overflow
        // also the max reduce value we can have is 2*(2^31 - 1) = 2^32 - 2 fits to 32 bits
        
        let tmp = Num::allocate_without_value(cs);
        let reduce = Num::allocate_without_value(cs);
        range_check_32_bits(cs, reduce.get_variable()); // 3rd constraint
        let result = Self::allocate_checked_without_value(cs, false); // 4th constraint
        
        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 4]| {
                let a = inputs[0].as_u64();
                let b = inputs[1].as_u64();
                let c = inputs[2].as_u64();
                let d = inputs[3].as_u64();
                
                let unredused_result = a * b + c * d;
                let reduce = unredused_result / M31_MODULUS;
                let result = unredused_result % M31_MODULUS;
                assert!(reduce < 1<<32);

                let mut tmp = F::from_u64_unchecked(a * b);
                tmp.sub_assign(&F::from_u64_unchecked(reduce * M31_MODULUS));

                [tmp, F::from_u64_unchecked(reduce), F::from_u64_unchecked(result)]
            };

            let dependencies = Place::from_variables([
                self.variable, other_mul.variable, other_add_1.variable, other_add_2.variable]);
            let outputs = Place::from_variables([tmp.get_variable(), reduce.get_variable(), result.variable]);

            cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
        }

        if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
            if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
                // (1) self * other_mul - reduce * modulus = tmp
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: *F::from_u64_unchecked(M31_MODULUS).negate(),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (self.variable, other_mul.variable),
                    linear_part: reduce.get_variable(),
                    rhs_part: tmp.get_variable(),
                };
                gate.add_to_cs(cs);

                // (2) tmp + other_add_1 * other_add_2 = result
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: F::ONE,
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (other_add_1.variable, other_add_2.variable),
                    linear_part: tmp.get_variable(),
                    rhs_part: result.variable,
                };
                gate.add_to_cs(cs);
            }
        } else {
            unimplemented!()
        }

        result
    }

    /// Computes self * other_mul + other_add_1 * other_add_2 + third_add
    pub fn two_mul_and_two_add<CS: ConstraintSystem<F>>(
        &self, 
        cs: &mut CS, 
        other_mul: &Self,
        other_add_1: &Self,
        other_add_2: &Self,
        third_add: &Self,
    ) -> Self {
        // We will use the following system of constraints:
        // (1) self * other_mul + third_add = tmp
        // (2) tmp + other_add_1 * other_add_2 = tmp2
        // (3) tmp2 - reduce * modulus= result
        // (4) reduce has 32 bits
        // (5) result has 31 bits
        // Note that max value of number is 2^31 - 1, so the max value of a*b + c*d + e is 2*(2^31-1)^2 + (2^31-1)
        // 2*(2^31-1)^2 + (2^31-1) = 2*(2^62 - 2*2^31 + 1) + (2^31-1) = 2^63 - 3*2^31 + 1
        // This fits to 63 bits, so no Goldilocks overflow
        // Also the max reduce value we can have is 2*(2^31 - 1) + 1 = 2^32 - 1 fits to 32 bits
        
        let tmp = Num::allocate_without_value(cs);
        let tmp2 = Num::allocate_without_value(cs);
        let reduce = Num::allocate_without_value(cs);
        range_check_32_bits(cs, reduce.get_variable()); // 4th constraint
        let result = Self::allocate_checked_without_value(cs, false); // 5th constraint
        
        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 5]| {
                let a = inputs[0].as_u64();
                let b = inputs[1].as_u64();
                let c = inputs[2].as_u64();
                let d = inputs[3].as_u64();
                let e = inputs[4].as_u64();
                
                let tmp = a * b + e;
                let tmp2 = tmp + c * d;
                let reduce = tmp2 / M31_MODULUS;
                let result = tmp2 % M31_MODULUS;
                assert!(reduce < 1<<32);

                [
                    F::from_u64_unchecked(tmp),
                    F::from_u64_unchecked(tmp2),
                    F::from_u64_unchecked(reduce),
                    F::from_u64_unchecked(result)
                ]
            };

            let dependencies = Place::from_variables(
                [self.variable, other_mul.variable, other_add_1.variable, other_add_2.variable, third_add.variable]);
            let outputs = Place::from_variables(
                [tmp.get_variable(), tmp2.get_variable(), reduce.get_variable(), result.variable]);

            cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
        }

        if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
            if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
                // (1) self * other_mul + third_add = tmp
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: F::ONE,
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (self.variable, other_mul.variable),
                    linear_part: third_add.get_variable(),
                    rhs_part: tmp.get_variable(),
                };
                gate.add_to_cs(cs);

                // (2) tmp + other_add_1 * other_add_2 = tmp2
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: F::ONE,
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (other_add_1.variable, other_add_2.variable),
                    linear_part: tmp.get_variable(),
                    rhs_part: tmp2.get_variable(),
                };
                gate.add_to_cs(cs);

                // (3) tmp2 - reduce * modulus= result
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: *F::from_u64_unchecked(M31_MODULUS).negate(),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (tmp2.get_variable(), cs.allocate_constant(F::ONE)),
                    linear_part: reduce.get_variable(),
                    rhs_part: result.get_variable(),
                };
                gate.add_to_cs(cs);
            }
        } else {
            unimplemented!()
        }

        result
    }

    /// Computes self * other_mul - other_sub_1 * other_sub_2
    pub fn two_mul_and_sub<CS: ConstraintSystem<F>>(
        &self, 
        cs: &mut CS, 
        other_mul: &Self,
        other_sub_1: &Self,
        other_sub_2: &Self,
    ) -> Self {
        // We will use the following system of constraints:
        // (1) self * other_mul + modulus^2 = tmp
        // (2) tmp - other_sub_1 * other_sub_2 = tmp2
        // (3) tmp2 - reduce * modulus = result
        // (4) reduce has 32 bits
        // (5) result has 31 bits
        // Note that max value of number is 2^31 - 1, so the value of a*b - c*d + modulus^2 is between 0 and 2*(2^31-1)^2
        // 2*(2^31-1)^2 = 2*(2^62 - 2*2^31 + 1) = 2^63 - 4*2^31 + 2 fits to 63 bits, so no Goldilocks overflow
        // also the max reduce value we can have is 2*(2^31 - 1) = 2^32 - 2 fits to 32 bits
        
        let tmp = Num::allocate_without_value(cs);
        let tmp2 = Num::allocate_without_value(cs);
        let reduce = Num::allocate_without_value(cs);
        range_check_32_bits(cs, reduce.get_variable()); // 4rd constraint
        let result = Self::allocate_checked_without_value(cs, false); // 5th constraint
        
        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 4]| {
                let a = inputs[0].as_u64();
                let b = inputs[1].as_u64();
                let c = inputs[2].as_u64();
                let d = inputs[3].as_u64();
                
                let tmp = a * b + M31_MODULUS * M31_MODULUS;
                let tmp2 = tmp - c * d;
                
                let reduce = tmp2 / M31_MODULUS;
                let result = tmp2 % M31_MODULUS;


                [F::from_u64_unchecked(tmp), F::from_u64_unchecked(tmp2), F::from_u64_unchecked(reduce), F::from_u64_unchecked(result)]
            };

            let dependencies = Place::from_variables([
                self.variable, other_mul.variable, other_sub_1.variable, other_sub_2.variable]);
            let outputs = Place::from_variables(
                [tmp.get_variable(), tmp2.get_variable(), reduce.get_variable(), result.variable]);

            cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
        }

        if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
            if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
                // (1) self * other_mul + modulus^2 = tmp
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: F::from_u64_unchecked(M31_MODULUS * M31_MODULUS),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (self.variable, other_mul.variable),
                    linear_part: cs.allocate_constant(F::ONE),
                    rhs_part: tmp.get_variable(),
                };
                gate.add_to_cs(cs);

                // (2) tmp - other_sub_1 * other_sub_2 = tmp2
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::MINUS_ONE,
                    linear_term_coeff: F::ONE,
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (other_sub_1.variable, other_sub_2.variable),
                    linear_part: tmp.get_variable(),
                    rhs_part: tmp2.variable,
                };
                gate.add_to_cs(cs);

                // (3) tmp2 - reduce * modulus = result
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: *F::from_u64_unchecked(M31_MODULUS).negate(),
                    linear_term_coeff: F::ONE,
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (reduce.variable, cs.allocate_constant(F::ONE)),
                    linear_part: tmp2.get_variable(),
                    rhs_part: result.variable,
                };
                gate.add_to_cs(cs);
            }
        } else {
            unimplemented!()
        }

        result
    }    


    /// Computes self * other_mul - other_sub_1 * other_sub_2 + third_add
    pub fn two_mul_and_sub_and_add<CS: ConstraintSystem<F>>(
        &self, 
        cs: &mut CS, 
        other_mul: &Self,
        other_sub_1: &Self,
        other_sub_2: &Self,
        third_add: &Self,
    ) -> Self {
        // We will use the following system of constraints:
        // (1) self * other_mul + modulus^2 = tmp
        // (2) tmp - other_sub_1 * other_sub_2 = tmp2
        // (3) tmp2 + third_add = tmp3
        // (4) tmp3 - reduce * modulus = result
        // (5) reduce has 32 bits
        // (6) result has 31 bits
        // Note that max value of number is 2^31 - 1, so the value of a*b - c*d + e modulus^2 is between 0 and 2*(2^31-1)^2 + (2^31-1)
        // 2*(2^31-1)^2 + (2^31-1) = 2*(2^62 - 2*2^31 + 1) + (2^31-1) = 2^63 - 3*2^31 + 1 fits to 63 bits, so no Goldilocks overflow
        // also the max reduce value we can have is 2*(2^31 - 1) + 1 = 2^32 - 1 fits to 32 bits
        
        let tmp = Num::allocate_without_value(cs);
        let tmp2 = Num::allocate_without_value(cs);
        let tmp3 = Num::allocate_without_value(cs);
        let reduce = Num::allocate_without_value(cs);
        range_check_32_bits(cs, reduce.get_variable()); // 5th constraint
        let result = Self::allocate_checked_without_value(cs, false); // 6th constraint
        
        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 5]| {
                let a = inputs[0].as_u64();
                let b = inputs[1].as_u64();
                let c = inputs[2].as_u64();
                let d = inputs[3].as_u64();
                let e = inputs[4].as_u64();
                
                let tmp = a * b + M31_MODULUS * M31_MODULUS;
                let tmp2 = tmp - c * d;
                let tmp3 = tmp2 + e;
                
                let reduce = tmp3 / M31_MODULUS;
                let result = tmp3 % M31_MODULUS;
                assert!(reduce < 1<<32);

                [
                    F::from_u64_unchecked(tmp),
                    F::from_u64_unchecked(tmp2),
                    F::from_u64_unchecked(tmp3),
                    F::from_u64_unchecked(reduce),
                    F::from_u64_unchecked(result)
                ]
            };

            let dependencies = Place::from_variables(
                [self.variable, other_mul.variable, other_sub_1.variable, other_sub_2.variable, third_add.variable]);
            let outputs = Place::from_variables(
                [tmp.get_variable(), tmp2.get_variable(), tmp3.get_variable(), reduce.get_variable(), result.variable]);

            cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
        }

        if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
            if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
                // (1) self * other_mul + modulus^2 = tmp
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: F::from_u64_unchecked(M31_MODULUS * M31_MODULUS),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (self.variable, other_mul.variable),
                    linear_part: cs.allocate_constant(F::ONE),
                    rhs_part: tmp.get_variable(),
                };
                gate.add_to_cs(cs);

                // (2) tmp - other_sub_1 * other_sub_2 = tmp2
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::MINUS_ONE,
                    linear_term_coeff: F::ONE,
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (other_sub_1.variable, other_sub_2.variable),
                    linear_part: tmp.get_variable(),
                    rhs_part: tmp2.variable,
                };
                gate.add_to_cs(cs);

                // (3) tmp2 + third_add = tmp3
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: F::ONE,
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (third_add.variable, cs.allocate_constant(F::ONE)),
                    linear_part: tmp2.get_variable(),
                    rhs_part: tmp3.variable,
                };
                gate.add_to_cs(cs);

                // (4) tmp3 - reduce * modulus = result
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: *F::from_u64_unchecked(M31_MODULUS).negate(),
                    linear_term_coeff: F::ONE,
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (reduce.variable, cs.allocate_constant(F::ONE)),
                    linear_part: tmp3.get_variable(),
                    rhs_part: result.variable,
                };
                gate.add_to_cs(cs);
            }
        } else {
            unimplemented!()
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
        // We will use the following system of constraints:
        // (1) result * self = not_zero + reduce * modulus
        // (2) self * (1 - not_zero) = 0
        // (3) result * (1 - not_zero) = 0
        // (4) reduce has 32 bits
        // (5) result has 31 bits
        //
        // so if self == 0 -(1, 4)-> not_zero == 0 -(3)-> result == 0
        // else if self != 0 -(2)-> not_zero == 1 -(1, 4, 5)-> result * self = 1 (mod modulus)
        // NOTE: if self == 2^31 - 1 -(2)-> not_zero == 1 -(1, 4, 5)-> contradiction (such result doesn't exist)

        let result = Self::allocate_checked_without_value(cs, false); // 5th constraint
        let not_zero = Num::allocate_without_value(cs);
        let reduce = Num::allocate_without_value(cs);
        range_check_32_bits(cs, reduce.get_variable()); // 4th constraint

        if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
            let value_fn = move |inputs: [F; 1]| {
                let a = inputs[0].as_u64();
                let a_inv = Mersenne31Field::new(a as u32)
                    .inverse()
                    .unwrap_or(Mersenne31Field::ZERO)
                    .to_reduced_u32() as u64;

                let (not_zero, reduce) = if a == 0 {
                    (0, 0)
                } else {
                    (1, (a * a_inv) / M31_MODULUS)
                };

                [F::from_u64_unchecked(a_inv), F::from_u64_unchecked(not_zero), F::from_u64_unchecked(reduce)]
            };

            let dependencies = Place::from_variables([self.variable]);
            let outputs = Place::from_variables([result.variable, not_zero.get_variable(), reduce.get_variable()]);

            cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
        }

        if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
            if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
                // (1) result * self + reduce * (-modulus) = not_zero
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::ONE,
                    linear_term_coeff: *F::from_u64_unchecked(M31_MODULUS).negate(),
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (self.variable, result.variable),
                    linear_part: reduce.get_variable(),
                    rhs_part: not_zero.get_variable(),
                };
                gate.add_to_cs(cs);

                // (2) self * (1 - not_zero) = 0
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::MINUS_ONE,
                    linear_term_coeff: F::ONE,
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (self.variable, not_zero.get_variable()),
                    linear_part: self.variable,
                    rhs_part: cs.allocate_constant(F::ZERO),
                };
                gate.add_to_cs(cs);

                // (3) result * (1 - not_zero) = 0
                let params = FmaGateInBaseWithoutConstantParams {
                    coeff_for_quadtaric_part: F::MINUS_ONE,
                    linear_term_coeff: F::ONE,
                };

                let gate = FmaGateInBaseFieldWithoutConstant {
                    params,
                    quadratic_part: (result.variable, not_zero.get_variable()),
                    linear_part: result.variable,
                    rhs_part: cs.allocate_constant(F::ZERO),
                };
                gate.add_to_cs(cs);
            }
        } else {
            unimplemented!()
        }

        result
    }

    pub fn is_zero<CS: ConstraintSystem<F>>(&mut self, cs: &mut CS) -> Boolean<F> {
        self.enforce_reduced(cs);
        self.into_num().is_zero(cs)
    }

    pub fn equals<CS: ConstraintSystem<F>>(&mut self, cs: &mut CS, other: &mut Self) -> Boolean<F> {
        self.enforce_reduced(cs);
        other.enforce_reduced(cs);
        Num::equals(cs, &self.into_num(), &other.into_num())
    }

    pub fn mask<CS: ConstraintSystem<F>>(&self, cs: &mut CS, masking_bit: Boolean<F>) -> Self {
        let inner = self.into_num().mask(cs, masking_bit);
        Self {
            variable: inner.get_variable(),
            reduced: self.reduced,
            _marker: std::marker::PhantomData,
        }
    }

    pub fn mask_negated<CS: ConstraintSystem<F>>(&self, cs: &mut CS, masking_bit: Boolean<F>) -> Self {
        let inner = self.into_num().mask_negated(cs, masking_bit);
        Self {
            variable: inner.get_variable(),
            reduced: self.reduced,
            _marker: std::marker::PhantomData,
        }
    }
}

fn range_check_32_bits<F: SmallField, CS: ConstraintSystem<F>>(
    cs: &mut CS,
    variable: Variable,
) {
    use crate::gadgets::non_native_field::implementations::get_16_bits_range_check_table;
    use crate::gadgets::u8::get_8_by_8_range_check_table;
    use crate::gadgets::impls::limbs_decompose::decompose_into_limbs;

    if let Some(table_id) = get_16_bits_range_check_table(&*cs) {
        let limbs = decompose_into_limbs::<F, CS, 2>(
            cs,
            F::from_u64_unchecked(1u64 << 16),
            variable,
        );

        cs.enforce_lookup::<1>(table_id, &[limbs[0]]);
        cs.enforce_lookup::<1>(table_id, &[limbs[1]]);
    } else if let Some(_table_id) = get_8_by_8_range_check_table(&*cs) {
        let _ = UInt32::from_variable_checked(cs, variable);
    } else {
        unimplemented!()
    }
}

fn range_check_31_bits<F: SmallField, CS: ConstraintSystem<F>>(
    cs: &mut CS,
    variable: Variable,
) {
    use crate::gadgets::non_native_field::implementations::get_16_bits_range_check_table;
    use crate::gadgets::impls::limbs_decompose::decompose_into_limbs;

    if let (Some(table_id_16), Some(table_id_15)) 
        = (get_16_bits_range_check_table(&*cs), get_15_bits_range_check_table(&*cs)) 
    {
        let limbs = decompose_into_limbs::<F, CS, 2>(
            cs,
            F::from_u64_unchecked(1u64 << 16),
            variable,
        );

        cs.enforce_lookup::<1>(table_id_16, &[limbs[0]]);
        cs.enforce_lookup::<1>(table_id_15, &[limbs[1]]);
    } else {
        unimplemented!()
    }
}
pub fn get_16_bits_range_check_table<F: SmallField, CS: ConstraintSystem<F>>(
    cs: &CS,
) -> Option<u32> {
    use crate::gadgets::tables::range_check_16_bits::RangeCheck16BitsTable;
    cs.get_table_id_for_marker::<RangeCheck16BitsTable<1>>()
}
pub fn get_15_bits_range_check_table<F: SmallField, CS: ConstraintSystem<F>>(
    cs: &CS,
) -> Option<u32> {
    use crate::gadgets::tables::range_check_16_bits::RangeCheck15BitsTable;
    cs.get_table_id_for_marker::<RangeCheck15BitsTable<1>>()
}

/// Returns a and reduce_a such that unreduced_a = a + reduce_a * modulus
/// a is constrainted to be in [0, 2^31 - 1]
/// reduce_a has no additional constraints
pub fn reduce_mersenne31<F: SmallField, CS: ConstraintSystem<F>>(
    cs: &mut CS,
    unreduced_a: Variable,
) -> (MersenneFiled<F>, Variable) {
    let a = MersenneFiled::allocate_checked_without_value(cs, false);
    let reduce_a = cs.alloc_variable_without_value();

    if <CS::Config as CSConfig>::WitnessConfig::EVALUATE_WITNESS {
        let value_fn = move |inputs: [F; 1]| {
            let unreduced_a = inputs[0].as_u64();

            [F::from_u64_unchecked(unreduced_a % M31_MODULUS), F::from_u64_unchecked(unreduced_a / M31_MODULUS)]
        };

        let dependencies = Place::from_variables([unreduced_a]);
        let outputs = Place::from_variables([a.get_variable(), reduce_a]);

        cs.set_values_with_dependencies(&dependencies, &outputs, value_fn);
    }
    if cs.gate_is_allowed::<ReductionGate<F, 2>>() {
        let params = ReductionGateParams {
            reduction_constants: [F::ONE, F::from_u64_unchecked(M31_MODULUS)],
        };

        let gate = ReductionGate {
            params,
            terms: [a.get_variable(), reduce_a],
            reduction_result: unreduced_a,
        };
        gate.add_to_cs(cs);
    } else if cs.gate_is_allowed::<FmaGateInBaseFieldWithoutConstant<F>>() {
        if <CS::Config as CSConfig>::SetupConfig::KEEP_SETUP == true {
            // (1) self * other_mul + modulus^2 = tmp
            let params = FmaGateInBaseWithoutConstantParams {
                coeff_for_quadtaric_part: F::ONE,
                linear_term_coeff: *F::from_u64_unchecked(M31_MODULUS).negate(),
            };

            let gate = FmaGateInBaseFieldWithoutConstant {
                params,
                quadratic_part: (unreduced_a, cs.allocate_constant(F::ONE)),
                linear_part: reduce_a,
                rhs_part: a.get_variable(),
            };
            gate.add_to_cs(cs);
        }
    } else {
        unimplemented!()
    }

    (a, reduce_a)
}

impl<F: SmallField> CSAllocatable<F> for MersenneFiled<F> {
    type Witness = Mersenne31Field;

    fn placeholder_witness() -> Self::Witness {
        Mersenne31Field::ZERO
    }
    fn allocate_without_value<CS: ConstraintSystem<F>>(cs: &mut CS) -> Self {
        Self::allocate_checked_without_value(cs, true)
    }
    fn allocate<CS: ConstraintSystem<F>>(cs: &mut CS, witness: Self::Witness) -> Self {
        Self::allocate_checked(cs, witness, false)
    }
}

impl<F: SmallField> WitnessCastable<F, [F; 1]> for Mersenne31Field {
    fn cast_from_source(witness: [F; 1]) -> Self {
        let value = witness[0].as_u64();
        assert!(value < M31_MODULUS);
        Mersenne31Field::new(value as u32)
    }

    fn cast_into_source(self) -> [F; 1] {
        [F::from_u64_unchecked(self.to_reduced_u32() as u64)]
    }
}

impl<F: SmallField> CSWitnessable<F, 1> for MersenneFiled<F> {
    type ConversionFunction = Convertor<F, [F; 1], Mersenne31Field>;

    fn witness_from_set_of_values(values: [F; 1]) -> Self::Witness {
        <Mersenne31Field as WitnessCastable<F, [F; 1]>>::cast_from_source(values)
    }

    fn as_variables_set(&self) -> [Variable; 1] {
        [self.variable]
    }
}

impl<F: SmallField> WitnessHookable<F> for MersenneFiled<F> {
    fn witness_hook<CS: ConstraintSystem<F>>(
        &self,
        cs: &CS,
    ) -> Box<dyn FnOnce() -> Option<Self::Witness>> {
        let raw_witness = self.get_witness(cs);
        Box::new(move || raw_witness.wait())
    }
}

impl<F: SmallField> Selectable<F> for MersenneFiled<F> {
    #[must_use]
    fn conditionally_select<CS: ConstraintSystem<F>>(
        cs: &mut CS,
        flag: Boolean<F>,
        a: &Self,
        b: &Self,
    ) -> Self {
        let res_variable = {
            let a_num = Num::from_variable(a.variable);
            let b_num = Num::from_variable(b.variable);
            let res_num = Num::conditionally_select(cs, flag, &a_num, &b_num);
            res_num.get_variable()
        };

        Self {
            variable: res_variable,
            reduced: a.reduced && b.reduced,
            _marker: std::marker::PhantomData,
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
        let a_nums = a.map(|el| Num::from_variable(el.variable));
        let b_nums = b.map(|el| Num::from_variable(el.variable));

        let tmp = Num::parallel_select(cs, flag, &a_nums, &b_nums);
        let mut res = unsafe { tmp.map(|el| 
            Self {
                variable: el.variable,
                reduced: false,
                _marker: std::marker::PhantomData,
            }
        )};

        for i in 0..N {
            if a[i].reduced && b[i].reduced {
                res[i].reduced = true;
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
        create_range_check_16_bits_table, RangeCheck16BitsTable,
        create_range_check_15_bits_table, RangeCheck15BitsTable,
    };
    use crate::gadgets::traits::witnessable::WitnessHookable;
    use crate::worker::Worker;

    type F = GoldilocksField;

    #[test]
    fn test_mersenne_field() {
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

        let rand_witness = [0; 4].map(|_| Mersenne31Field::new(rand::random::<u32>() % M31_MODULUS as u32));
        let mut rand_vars = rand_witness.map(|w| MersenneFiled::<F>::allocate_checked(cs, w, false));

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

        // mul_negate
        let mut res_witness = rand_witness[0];
        res_witness.mul_assign(&rand_witness[1]);
        res_witness.negate();
        let res_var = rand_vars[0].mul_negate(cs, &rand_vars[1]);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // mul_times_2
        let mut res_witness = rand_witness[0];
        res_witness.mul_assign(&rand_witness[1]);   
        res_witness.double();
        let res_var = rand_vars[0].mul_times_2(cs, &rand_vars[1]);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // mul_and_add
        let mut res_witness = rand_witness[0];
        res_witness.mul_assign(&rand_witness[1]);
        res_witness.add_assign(&rand_witness[2]);
        let res_var = rand_vars[0].mul_and_add(cs, &rand_vars[1], &rand_vars[2]);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // two_mul_and_add
        let mut res_witness = rand_witness[0];
        res_witness.mul_assign(&rand_witness[1]);
        let mut tmp = rand_witness[2];
        tmp.mul_assign(&rand_witness[3]);
        res_witness.add_assign(&tmp);
        let res_var = rand_vars[0].two_mul_and_add(cs, &rand_vars[1], &rand_vars[2], &rand_vars[3]);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // two_mul_and_sub
        let mut res_witness = rand_witness[0];
        res_witness.mul_assign(&rand_witness[1]);
        let mut tmp = rand_witness[2];
        tmp.mul_assign(&rand_witness[3]);
        res_witness.sub_assign(&tmp);
        let res_var = rand_vars[0].two_mul_and_sub(cs, &rand_vars[1], &rand_vars[2], &rand_vars[3]);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // div
        let mut res_witness = rand_witness[0];
        res_witness.mul_assign(&rand_witness[1].inverse().unwrap_or(Mersenne31Field::ZERO));
        let res_var = rand_vars[0].div(cs, &rand_vars[1]);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());

        // inverse_or_zero
        let mut res_witness = rand_witness[0];
        res_witness = res_witness.inverse().unwrap_or(Mersenne31Field::ZERO);
        let res_var = rand_vars[0].inverse_or_zero(cs);
        assert_eq!(res_witness, res_var.witness_hook(&*cs)().unwrap());


        let worker = Worker::new_with_num_threads(8);

        drop(cs);
        owned_cs.pad_and_shrink();
        let mut owned_cs = owned_cs.into_assembly::<Global>();
        assert!(owned_cs.check_if_satisfied(&worker));
    }
}
