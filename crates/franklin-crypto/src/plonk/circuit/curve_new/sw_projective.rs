use super::*;
use crate::plonk::circuit::bigint_new::*;

use crate::bellman::pairing::{Engine, GenericCurveAffine, GenericCurveProjective};

use crate::bellman::pairing::ff::{BitIterator, Field, PrimeField, PrimeFieldRepr, ScalarEngine};

use crate::bellman::SynthesisError;

use crate::bellman::plonk::better_better_cs::cs::{
    ArithmeticTerm, Coefficient, ConstraintSystem, Gate, GateInternal, LinearCombinationOfTerms, MainGate, MainGateTerm, PlonkConstraintSystemParams, PlonkCsWidth4WithNextStepParams,
    PolynomialInConstraint, PolynomialMultiplicativeTerm, TimeDilation, TrivialAssembly, Variable, Width4MainGateWithDNext,
};

use super::super::allocated_num::{AllocatedNum, Num};
use super::super::boolean::{AllocatedBit, Boolean};
use super::super::linear_combination::LinearCombination;
use super::super::simple_term::Term;
use crate::plonk::circuit::Assignment;

#[derive(Clone, Debug)]
pub struct ProjectivePoint<'a, E: Engine, G: GenericCurveAffine>
where
    <G as GenericCurveAffine>::Base: PrimeField,
{
    pub x: FieldElement<'a, E, G::Base>,
    pub y: FieldElement<'a, E, G::Base>,
    pub z: FieldElement<'a, E, G::Base>,
    pub value: Option<G::Projective>,
}

impl<'a, E: Engine, G: GenericCurveAffine> From<AffinePoint<'a, E, G>> for ProjectivePoint<'a, E, G>
where
    <G as GenericCurveAffine>::Base: PrimeField,
{
    fn from(affine_pt: AffinePoint<'a, E, G>) -> Self {
        let params = affine_pt.x.representation_params;
        let AffinePoint { x, y, value } = affine_pt;

        ProjectivePoint::<E, G> {
            x,
            y,
            z: FieldElement::one(&params),
            value: value.map(|x| x.into_projective()),
        }
    }
}

impl<'a, E: Engine, G: GenericCurveAffine> ProjectivePoint<'a, E, G>
where
    <G as GenericCurveAffine>::Base: PrimeField,
{
    pub fn get_x(&self) -> FieldElement<'a, E, G::Base> {
        self.x.clone()
    }

    pub fn get_y(&self) -> FieldElement<'a, E, G::Base> {
        self.y.clone()
    }

    pub fn get_z(&self) -> FieldElement<'a, E, G::Base> {
        self.z.clone()
    }

    pub fn zero(params: &'a RnsParameters<E, G::Base>) -> Self {
        let x = FieldElement::zero(params);
        let y = FieldElement::one(params);
        let z = FieldElement::zero(params);
        let value = Some(G::Projective::zero());

        Self { x, y, z, value }
    }

    pub fn is_constant(&self) -> bool {
        self.x.is_constant() & self.y.is_constant() & self.z.is_constant()
    }

    pub fn get_value(&self) -> Option<G> {
        self.value.map(|el| el.into_affine())
    }

    pub fn negate<CS: ConstraintSystem<E>>(&self, cs: &mut CS) -> Result<Self, SynthesisError> {
        let y_negated = self.y.negate(cs)?;
        let new_value = self.value.map(|el| {
            let mut t = el;
            t.negate();
            t
        });

        let new = Self {
            x: self.x.clone(),
            y: y_negated,
            z: self.z.clone(),
            value: new_value,
        };

        Ok(new)
    }

    pub fn sub<CS: ConstraintSystem<E>>(&self, cs: &mut CS, other: &Self) -> Result<Self, SynthesisError> {
        let other_negated = other.negate(cs)?;
        self.add(cs, &other_negated)
    }

    pub unsafe fn convert_to_affine<CS>(&self, cs: &mut CS) -> Result<AffinePoint<'a, E, G>, SynthesisError>
    where
        CS: ConstraintSystem<E>,
    {
        let x = self.x.div(cs, &self.z)?;
        let y = self.y.div(cs, &self.z)?;
        let value = self.get_value();

        Ok(AffinePoint { x, y, value })
    }

    pub fn convert_to_affine_or_default<CS: ConstraintSystem<E>>(&mut self, cs: &mut CS, default: &AffinePoint<'a, E, G>) -> Result<(AffinePoint<'a, E, G>, Boolean), SynthesisError> {
        let params = self.x.representation_params;
        let is_point_at_infty = self.z.is_zero(cs)?;
        let safe_z = FieldElement::conditionally_select(cs, &is_point_at_infty, &FieldElement::one(params), &self.z)?;
        let x_for_safe_z = self.x.div(cs, &safe_z)?;
        let y_for_safe_z = self.y.div(cs, &safe_z)?;
        let x = FieldElement::conditionally_select(cs, &is_point_at_infty, &default.x, &x_for_safe_z)?;
        let y = FieldElement::conditionally_select(cs, &is_point_at_infty, &default.y, &y_for_safe_z)?;

        let value = match (is_point_at_infty.get_value(), self.get_value(), default.get_value()) {
            (Some(true), _, Some(val)) | (Some(false), Some(val), _) => Some(val),
            _ => None,
        };

        let new = AffinePoint { x, y, value };
        Ok((new, is_point_at_infty))
    }

    #[track_caller]
    pub fn add<CS: ConstraintSystem<E>>(&self, cs: &mut CS, other: &Self) -> Result<Self, SynthesisError> {
        // this formula is only valid for curve with zero j-ivariant
        assert!(G::a_coeff().is_zero());

        let params = self.x.representation_params;
        let curve_b = G::b_coeff();
        let mut curve_b3 = curve_b;
        curve_b3.double();
        curve_b3.add_assign(&curve_b);
        let b3 = FieldElement::constant(curve_b3, params);

        let this_value = self.value;
        let other_value = other.value;
        let x1 = self.x.clone();
        let y1 = self.y.clone();
        let z1 = self.z.clone();
        let x2 = other.x.clone();
        let y2 = other.y.clone();
        let z2 = other.z.clone();

        // exception free addition in projective coordiantes
        // 1. t0 ← X1 · X2
        let t0 = x1.mul(cs, &x2)?;
        // 2. t1 ← Y1 · Y2
        let t1 = y1.mul(cs, &y2)?;
        // 3. t2 ← Z1 · Z2
        let t2 = z1.mul(cs, &z2)?;
        // 4. t3 ← X1 + Y1
        let t3 = x1.add(cs, &y1)?;
        // 5. t4 ← X2 + Y2
        let t4 = x2.add(cs, &y2)?;
        // 6. t3 ← t3 · t4
        let t3 = t3.mul(cs, &t4)?;
        // 7. t4 ← t0 + t1
        let t4 = t0.add(cs, &t1)?;
        // 8. t3 ← t3 − t4
        let t3 = t3.sub(cs, &t4)?;
        // 9. t4 ← Y1 + Z1
        let t4 = y1.add(cs, &z1)?;
        // 10. X3 ← Y2 + Z2
        let x3 = y2.add(cs, &z2)?;
        // 11. t4 ← t4 · X3
        let t4 = t4.mul(cs, &x3)?;
        // 12. X3 ← t1 + t2
        let x3 = t1.add(cs, &t2)?;
        // 13. t4 ← t4 − X3
        let t4 = t4.sub(cs, &x3)?;
        // 14. X3 ← X1 + Z1
        let x3 = x1.add(cs, &z1)?;
        // 15. Y3 ← X2 + Z2
        let y3 = x2.add(cs, &z2)?;
        // 16. X3 ← X3 · Y3
        let x3 = x3.mul(cs, &y3)?;
        // 17. Y3 ← t0 + t2
        let y3 = t0.add(cs, &t2)?;
        // 18. Y3 ← X3 − Y3
        let y3 = x3.sub(cs, &y3)?;
        // 19. X3 ← t0 + t0
        let x3 = t0.double(cs)?;
        // 20. t0 ← X3 + t0
        let t0 = x3.add(cs, &t0)?;
        // 21. t2 ← b3 · t2
        let t2 = b3.mul(cs, &t2)?;
        // 22. Z3 ← t1 + t2
        let z3 = t1.add(cs, &t2)?;
        // 23. t1 ← t1 − t2
        let t1 = t1.sub(cs, &t2)?;
        // 24. Y3 ← b3 · Y3
        let y3 = b3.mul(cs, &y3)?;
        // 25. X3 ← t4 · Y3
        let x3 = t4.mul(cs, &y3)?;
        // 26. t2 ← t3 · t1
        let t2 = t3.mul(cs, &t1)?;
        // 27. X3 ← t2 − X3
        let x3 = t2.sub(cs, &x3)?;
        // 28. Y3 ← Y3 · t0
        let y3 = y3.mul(cs, &t0)?;
        // 29. t1 ← t1 · Z3
        let t1 = t1.mul(cs, &z3)?;
        // 30. Y3 ← t1 + Y3
        let y3 = t1.add(cs, &y3)?;
        // 31. t0 ← t0 · t3
        let t0 = t0.mul(cs, &t3)?;
        // 32. Z3 ← Z3 · t4
        let z3 = z3.mul(cs, &t4)?;
        // 33. Z3 ← Z3 + t0
        let z3 = z3.add(cs, &t0)?;

        let new_value = match (this_value, other_value) {
            (Some(this), Some(other)) => {
                let mut tmp = this;
                tmp.add_assign(&other);

                Some(tmp)
            }
            _ => None,
        };

        let new = Self {
            x: x3,
            y: y3,
            z: z3,
            value: new_value,
        };
        Ok(new)
    }

    #[track_caller]
    pub fn double<CS: ConstraintSystem<E>>(&self, cs: &mut CS) -> Result<Self, SynthesisError> {
        // this formula is only valid for curve with zero j-ivariant
        assert!(G::a_coeff().is_zero());

        let params = self.x.representation_params;
        let curve_b = G::b_coeff();
        let mut curve_b3 = curve_b;
        curve_b3.double();
        curve_b3.add_assign(&curve_b);
        let b3 = FieldElement::constant(curve_b3, params);
        let this_value = self.value;

        let x = self.x.clone();
        let y = self.y.clone();
        let z = self.z.clone();

        // 1. t0 ← Y · Y
        let t0 = y.square(cs)?;
        // 2. Z3 ← t0 + t0
        let z3 = t0.double(cs)?;
        // 3. Z3 ← Z3 + Z3
        let z3 = z3.double(cs)?;
        // 4. Z3 ← Z3 + Z3
        let z3 = z3.double(cs)?;
        // 5. t1 ← Y · Z
        let t1 = y.mul(cs, &z)?;
        // 6. t2 ← Z · Z
        let t2 = z.square(cs)?;
        // 7. t2 ← b3 · t2
        let t2 = b3.mul(cs, &t2)?;
        // 8. X3 ← t2 · Z3
        let x3 = t2.mul(cs, &z3)?;
        // 9. Y3 ← t0 + t2
        let y3 = t0.add(cs, &t2)?;
        // 10. Z3 ← t1 · Z3
        let z3 = t1.mul(cs, &z3)?;
        // 11. t1 ← t2 + t2
        let t1 = t2.double(cs)?;
        // 12. t2 ← t1 + t2
        let t2 = t1.add(cs, &t2)?;
        // 13. t0 ← t0 − t2
        let t0 = t0.sub(cs, &t2)?;
        // 14. Y3 ← t0 · Y3
        let y3 = t0.mul(cs, &y3)?;
        // 15. Y3 ← X3 + Y3
        let y3 = x3.add(cs, &y3)?;
        // 16. t1 ← X · Y
        let t1 = x.mul(cs, &y)?;
        // 17. X3 ← t0 · t1
        let x3 = t0.mul(cs, &t1)?;
        // 18. X3 ← X3 + X3
        let x3 = x3.double(cs)?;

        let new_value = this_value.map(|el| {
            let mut tmp = el;
            tmp.double();
            tmp
        });

        let new = Self {
            x: x3,
            y: y3,
            z: z3,
            value: new_value,
        };
        Ok(new)
    }

    #[track_caller]
    pub fn add_mixed<CS>(&self, cs: &mut CS, other: &AffinePoint<'a, E, G>) -> Result<Self, SynthesisError>
    where
        CS: ConstraintSystem<E>,
    {
        // this formula is only valid for curve with zero j-ivariant
        assert!(G::a_coeff().is_zero());

        let params = self.x.representation_params;
        let curve_b = G::b_coeff();
        let mut curve_b3 = curve_b;
        curve_b3.double();
        curve_b3.add_assign(&curve_b);
        let b3 = FieldElement::constant(curve_b3, params);

        let x1 = self.x.clone();
        let y1 = self.y.clone();
        let z1 = self.z.clone();
        let x2 = other.x.clone();
        let y2 = other.y.clone();

        // 1. t0 ← X1 · X2
        let t0 = x1.mul(cs, &x2)?;
        // 2. t1 ← Y1 · Y2
        let t1 = y1.mul(cs, &y2)?;
        // 3. t3 ← X2 + Y2
        let t3 = x2.add(cs, &y2)?;
        // 4. t4 ← X1 + Y1
        let t4 = x1.add(cs, &y1)?;
        // 5. t3 ← t3 · t4
        let t3 = t3.mul(cs, &t4)?;
        // 6. t4 ← t0 + t1
        let t4 = t0.add(cs, &t1)?;
        // 7. t3 ← t3 − t4
        let t3 = t3.sub(cs, &t4)?;
        // 8. t4 ← Y2 · Z1
        let t4 = y2.mul(cs, &z1)?;
        // 9. t4 ← t4 + Y1
        let t4 = t4.add(cs, &y1)?;
        // 10. Y3 ← X2 · Z1
        let y3 = x2.mul(cs, &z1)?;
        // 11. Y3 ← Y3 + X1
        let y3 = y3.add(cs, &x1)?;
        // 12. X3 ← t0 + t0
        let x3 = t0.double(cs)?;
        // 13. t0 ← X3 + t0
        let t0 = x3.add(cs, &t0)?;
        // 14. t2 ← b3 · Z1
        let t2 = b3.mul(cs, &z1)?;
        // 15. Z3 ← t1 + t2
        let z3 = t1.add(cs, &t2)?;
        // 16. t1 ← t1 − t2
        let t1 = t1.sub(cs, &t2)?;
        // 17. Y3 ← b3 · Y3
        let y3 = b3.mul(cs, &y3)?;
        // 18. X3 ← t4 · Y3
        let x3 = t4.mul(cs, &y3)?;
        // 19. t2 ← t3 · t1
        let t2 = t3.mul(cs, &t1)?;
        // 20. X3 ← t2 − X3
        let x3 = t2.sub(cs, &x3)?;
        // 21. Y3 ← Y3 · t0
        let y3 = y3.mul(cs, &t0)?;
        // 22. t1 ← t1 · Z3
        let t1 = t1.mul(cs, &z3)?;
        // 23. Y3 ← t1 + Y3
        let y3 = t1.add(cs, &y3)?;
        // 24. t0 ← t0 · t3
        let t0 = t0.mul(cs, &t3)?;
        // 25. Z3 ← Z3 · t4
        let z3 = z3.mul(cs, &t4)?;
        // 26. Z3 ← Z3 + t0
        let z3 = z3.add(cs, &t0)?;

        let new_value = match (self.value, other.get_value()) {
            (Some(this), Some(other)) => {
                let mut tmp = this;
                tmp.add_assign_mixed(&other);
                Some(tmp)
            }
            _ => None,
        };

        let new = Self {
            x: x3,
            y: y3,
            z: z3,
            value: new_value,
        };
        Ok(new)
    }

    pub fn conditionally_select<CS: ConstraintSystem<E>>(cs: &mut CS, flag: &Boolean, first: &Self, second: &Self) -> Result<Self, SynthesisError> {
        let first_value = first.value;
        let second_value = second.value;
        let x = FieldElement::conditionally_select(cs, flag, &first.x, &second.x)?;
        let y = FieldElement::conditionally_select(cs, flag, &first.y, &second.y)?;
        let z = FieldElement::conditionally_select(cs, flag, &first.z, &second.z)?;

        let value = match (flag.get_value(), first_value, second_value) {
            (Some(true), Some(p), _) => Some(p),
            (Some(false), _, Some(p)) => Some(p),
            (_, _, _) => None,
        };

        let selected = Self { x, y, z, value };
        Ok(selected)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::bellman::pairing::bn256::{Bn256, Fq, Fr, G1Affine};
    use bellman::plonk::better_better_cs::cs::*;
    use bellman::plonk::better_better_cs::gates::{self, selector_optimized_with_d_next::SelectorOptimizedWidth4MainGateWithDNext};
    use plonk::circuit::Width4WithCustomGates;
    use crate::rand::{Rng, SeedableRng, XorShiftRng};

    #[test]
    fn test_arithmetic_for_projective_bn256_curve() {
        let mut cs = TrivialAssembly::<Bn256, Width4WithCustomGates, SelectorOptimizedWidth4MainGateWithDNext>::new();
        inscribe_default_bitop_range_table(&mut cs).unwrap();
        let params = RnsParameters::<Bn256, Fq>::new_optimal(&mut cs, 80usize);
        let scalar_params = RnsParameters::<Bn256, Fr>::new_optimal(&mut cs, 80usize);
        let mut rng = crate::rand::thread_rng();

        let a: G1Affine = rng.gen();
        let b: G1Affine = rng.gen();
        let mut tmp = a.into_projective();
        tmp.add_assign_mixed(&b);
        let result = tmp.into_affine();

        let mut a_affine = AffinePoint::alloc(&mut cs, Some(a), &params).unwrap();
        let mut a_proj = ProjectivePoint::from(a_affine);
        let mut b_affine = AffinePoint::alloc(&mut cs, Some(b), &params).unwrap();
        let mut actual_result = AffinePoint::alloc(&mut cs, Some(result), &params).unwrap();
        let naive_mul_start = cs.get_current_step_number();
        let mut result = a_proj.add_mixed(&mut cs, &b_affine).unwrap();
        let naive_mul_end = cs.get_current_step_number();
        println!("num of gates: {}", naive_mul_end - naive_mul_start);

        let mut result = unsafe { result.convert_to_affine(&mut cs).unwrap() };
        println!("WOW");
        // //result.x.normalize(&mut cs).unwrap();
        // result.y.normalize(&mut cs).unwrap();
        AffinePoint::enforce_equal(&mut cs, &mut result, &mut actual_result).unwrap();
        assert!(cs.is_satisfied());
        println!("PROJ MIXED ADD 2");
    }
}
