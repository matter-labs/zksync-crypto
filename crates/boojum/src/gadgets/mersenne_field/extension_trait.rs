use super::*;
use fourth_ext::MersenneQuartic;
use second_ext::MersenneComplex;

pub trait CircuitFieldExpression<F: SmallField, BaseField> {
    fn add_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &BaseField) -> Self;
    fn sub_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &BaseField) -> Self;
    fn mul_by_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &BaseField) -> Self;
    fn mul_by_base_and_add<CS: ConstraintSystem<F>>(
        &self,
        cs: &mut CS,
        mul_base: &BaseField,
        add: &Self,
    ) -> Self;
}

impl<F: SmallField> CircuitFieldExpression<F, MersenneField<F>> for MersenneField<F> {
    fn add_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneField<F>) -> Self {
        self.add(cs, base)
    }

    fn sub_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneField<F>) -> Self {
        self.sub(cs, base)
    }

    fn mul_by_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneField<F>) -> Self {
        self.mul(cs, base)
    }

    fn mul_by_base_and_add<CS: ConstraintSystem<F>>(
        &self,
        cs: &mut CS,
        mul_base: &MersenneField<F>,
        add: &Self,
    ) -> Self {
        self.mul_and_add(cs, mul_base, add)
    }
}

impl<F: SmallField> CircuitFieldExpression<F, MersenneField<F>> for MersenneComplex<F> {
    fn add_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneField<F>) -> Self {
        Self {
            x: self.x.add(cs, base),
            y: self.y,
        }
    }

    fn sub_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneField<F>) -> Self {
        Self {
            x: self.x.sub(cs, base),
            y: self.y,
        }
    }

    fn mul_by_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneField<F>) -> Self {
        Self {
            x: self.x.mul(cs, base),
            y: self.y.mul(cs, base),
        }
    }

    fn mul_by_base_and_add<CS: ConstraintSystem<F>>(
        &self,
        cs: &mut CS,
        mul_base: &MersenneField<F>,
        add: &Self,
    ) -> Self {
        Self {
            x: self.x.mul_and_add(cs, mul_base, &add.x),
            y: self.y.mul_and_add(cs, mul_base, &add.y),
        }
    }
}

impl<F: SmallField> CircuitFieldExpression<F, MersenneField<F>> for MersenneQuartic<F> {
    fn add_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneField<F>) -> Self {
        Self {
            x: self.x.add_base(cs, base),
            y: self.y,
        }
    }

    fn sub_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneField<F>) -> Self {
        Self {
            x: self.x.sub_base(cs, base),
            y: self.y,
        }
    }

    fn mul_by_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneField<F>) -> Self {
        Self {
            x: self.x.mul_by_base(cs, base),
            y: self.y.mul_by_base(cs, base),
        }
    }

    fn mul_by_base_and_add<CS: ConstraintSystem<F>>(
        &self,
        cs: &mut CS,
        mul_base: &MersenneField<F>,
        add: &Self,
    ) -> Self {
        Self {
            x: self.x.mul_by_base_and_add(cs, mul_base, &add.x),
            y: self.y.mul_by_base_and_add(cs, mul_base, &add.y),
        }
    }
}

impl<F: SmallField> CircuitFieldExpression<F, MersenneComplex<F>> for MersenneComplex<F> {
    fn add_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneComplex<F>) -> Self {
        self.add(cs, base)
    }

    fn sub_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneComplex<F>) -> Self {
        self.sub(cs, base)
    }

    fn mul_by_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneComplex<F>) -> Self {
        self.mul(cs, base)
    }

    fn mul_by_base_and_add<CS: ConstraintSystem<F>>(
        &self,
        cs: &mut CS,
        mul_base: &MersenneComplex<F>,
        add: &Self,
    ) -> Self {
        self.mul_and_add(cs, mul_base, add)
    }
}

impl<F: SmallField> CircuitFieldExpression<F, MersenneComplex<F>> for MersenneQuartic<F> {
    fn add_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneComplex<F>) -> Self {
        Self {
            x: self.x.add(cs, base),
            y: self.y,
        }
    }

    fn sub_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneComplex<F>) -> Self {
        Self {
            x: self.x.sub(cs, base),
            y: self.y,
        }
    }

    fn mul_by_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneComplex<F>) -> Self {
        Self {
            x: self.x.mul(cs, base),
            y: self.y.mul(cs, base),
        }
    }

    fn mul_by_base_and_add<CS: ConstraintSystem<F>>(
        &self,
        cs: &mut CS,
        mul_base: &MersenneComplex<F>,
        add: &Self,
    ) -> Self {
        Self {
            x: self.x.mul_and_add(cs, mul_base, &add.x),
            y: self.y.mul_and_add(cs, mul_base, &add.y),
        }
    }
}

impl<F: SmallField> CircuitFieldExpression<F, MersenneQuartic<F>> for MersenneQuartic<F> {
    fn add_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneQuartic<F>) -> Self {
        self.add(cs, base)
    }

    fn sub_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneQuartic<F>) -> Self {
        self.sub(cs, base)
    }

    fn mul_by_base<CS: ConstraintSystem<F>>(&self, cs: &mut CS, base: &MersenneQuartic<F>) -> Self {
        self.mul(cs, base)
    }

    fn mul_by_base_and_add<CS: ConstraintSystem<F>>(
        &self,
        cs: &mut CS,
        mul_base: &MersenneQuartic<F>,
        add: &Self,
    ) -> Self {
        self.mul_and_add(cs, mul_base, add)
    }
}
