use super::*;
use crate::bellman::plonk::better_better_cs::cs::{
    Variable, Index, ConstraintSystem, ArithmeticTerm, MainGateTerm, Gate, MainGate
};
use super::super::simple_term::Term;
use crate::bellman::plonk::better_better_cs::cs::PlonkConstraintSystemParams;
use crate::indexmap::{IndexMap, IndexSet};
use crate::plonk::circuit::utils::is_selector_specialized_gate;

const REQUIRED_STATE_WIDTH : usize = 4;


#[derive(Clone, Debug)]
struct GateConstructorHelper<E: Engine> {
    a: Variable, b: Variable, c: Variable, d: Variable,
    a_mul_b_coef: E::Fr, a_mul_c_coef: E::Fr,
    a_linear_coef: E::Fr, b_linear_coef: E::Fr, c_linear_coef: E::Fr, d_linear_coef: E::Fr,
    cnst: E::Fr, free_vars_start_idx: usize, free_vars_end_idx: usize, is_final: bool 
}

impl<E: Engine> GateConstructorHelper<E> {
    pub fn new<CS: ConstraintSystem<E>>(cs: &mut CS) -> Self {
        let dummy = AllocatedNum::zero(cs).get_variable();
        GateConstructorHelper::<E> {
            a: dummy, b: dummy, c: dummy, d: dummy, a_mul_b_coef: E::Fr::zero(), a_mul_c_coef: E::Fr::zero(),
            a_linear_coef: E::Fr::zero(), b_linear_coef: E::Fr::zero(), c_linear_coef: E::Fr::zero(), 
            d_linear_coef: E::Fr::zero(), cnst: E::Fr::zero(), 
            free_vars_start_idx: 0, free_vars_end_idx: REQUIRED_STATE_WIDTH, is_final: false 
        }
    }

    pub fn new_for_mul<CS: ConstraintSystem<E>>(cs: &mut CS, vars: OrderedVariablePair, coef: E::Fr) -> Self {
        let dummy = AllocatedNum::zero(cs).get_variable();
        GateConstructorHelper::<E> {
            a: vars.first, b: vars.second, c: dummy, d: dummy, a_mul_b_coef: coef, a_mul_c_coef: E::Fr::zero(),
            a_linear_coef: E::Fr::zero(), b_linear_coef: E::Fr::zero(), c_linear_coef: E::Fr::zero(), 
            d_linear_coef: E::Fr::zero(), cnst: E::Fr::zero(), 
            free_vars_start_idx: 2, free_vars_end_idx: REQUIRED_STATE_WIDTH, is_final: false 
        }
    }

    pub fn new_for_pair_of_muls<CS: ConstraintSystem<E>>(
        cs: &mut CS, a: Variable, b: Variable, c: Variable, a_mul_b_coef: E::Fr, a_mul_c_coef: E::Fr
    ) -> Self {
        let dummy = AllocatedNum::zero(cs).get_variable();
        let zero = E::Fr::zero();
        GateConstructorHelper::<E> {
            a, b, c, d: dummy, a_mul_b_coef, a_mul_c_coef, 
            a_linear_coef: zero, b_linear_coef: zero, c_linear_coef: zero, d_linear_coef: zero, 
            cnst: zero, free_vars_start_idx: 3, free_vars_end_idx: REQUIRED_STATE_WIDTH, is_final: false 
        }
    }

    pub fn set_finality_flag(&mut self) {
        self.is_final = true;
    }

    pub fn add_constant_term(&mut self, coef: E::Fr) {
        self.cnst = coef;
    }

    pub fn add_next_trace_step_term(&mut self, var: Variable) {
        self.d = var;
        self.d_linear_coef = E::Fr::one();
        self.free_vars_end_idx -= 1;
    }

    pub fn add_linear_coefficients_for_bound_variables(&mut self, var_map: &mut IndexMap<Variable, E::Fr>) {
        let iter = std::array::IntoIter::new([
            (self.a, &mut self.a_linear_coef), (self.b, &mut self.b_linear_coef), 
            (self.c, &mut self.c_linear_coef), (self.d, &mut self.d_linear_coef)
        ]).take(self.free_vars_start_idx);
        for (var, pos) in iter {
            let fr = var_map.remove(&var).unwrap_or(E::Fr::zero());
            *pos = fr;
        }
    }

    pub fn add_linear_coefficients_for_free_variables(
        &mut self, var_map: &mut IndexMap<Variable, E::Fr>, free_vars_set: &mut IndexSet<Variable>
    ) {
        // Rust is not an expressive language at all, that's why we unroll the loop manually
        // How much better it would look in C++.. Mmm...
        let mut vars = [self.a, self.b, self.c, self.d];
        let mut coefs = [self.a_linear_coef, self.b_linear_coef, self.c_linear_coef, self.d_linear_coef];
        
        let num_elems = self.free_vars_end_idx - self.free_vars_start_idx;
        for (var_ptr, coef_ptr) in vars.iter_mut().zip(coefs.iter_mut()).skip(self.free_vars_start_idx).take(num_elems) {
            if free_vars_set.is_empty() {
                break;
            }

            let elt = free_vars_set.iter().next().cloned().unwrap();
            let fr = var_map.remove(&elt).unwrap_or(E::Fr::zero());
            *var_ptr = elt;
            *coef_ptr = fr;
            free_vars_set.take(&elt).unwrap();
            self.free_vars_start_idx += 1;
        }

        // boilerplate code, because it is rusty rust
        self.a = vars[0];
        self.b = vars[1];
        self.c = vars[2];
        self.d = vars[3];

        self.a_linear_coef = coefs[0];
        self.b_linear_coef = coefs[1];
        self.c_linear_coef = coefs[2];
        self.d_linear_coef = coefs[3];
    }

    pub fn materialize<CS: ConstraintSystem<E>>(&self, cs: &mut CS) -> Result<(), SynthesisError> {
        let range_of_multiplicative_terms = CS::MainGate::range_of_multiplicative_term();
        let range_of_linear_terms = CS::MainGate::range_of_linear_terms();
        let index_for_constant_term = CS::MainGate::index_for_constant_term();
        let range_of_next_step_linear_terms =  CS::MainGate::range_of_next_step_linear_terms();

        let dummy = AllocatedNum::zero(cs).get_variable();
        let gate_term = MainGateTerm::new();
        let (_vars, mut coefs) = CS::MainGate::format_term(gate_term, dummy)?;
        let vars = [self.a, self.b, self.c, self.d];

        for (pos, fr) in range_of_multiplicative_terms.zip(&[self.a_mul_b_coef, self.a_mul_c_coef]) {
            coefs[pos] = *fr; 
        }

        let linear_coefs = [self.a_linear_coef, self.b_linear_coef, self.c_linear_coef, self.d_linear_coef];
        for (pos, fr) in range_of_linear_terms.zip(&linear_coefs[..]) {
            coefs[pos] = *fr; 
        }

        coefs[index_for_constant_term] = self.cnst;
        if !self.is_final {
            let mut minus_one = E::Fr::one();
            minus_one.negate();
            coefs[range_of_next_step_linear_terms.last().unwrap()] = minus_one;
        }

        let mg = CS::MainGate::default();
        cs.new_single_gate_for_trace_step(&mg, &coefs, &vars, &[])
    }
        
    pub fn evaluate_next_trace_step_value<CS: ConstraintSystem<E>>(&self, cs: &mut CS) -> Result<E::Fr, SynthesisError> {
        let a_val = cs.get_value(self.a)?;
        let b_val = cs.get_value(self.b)?;
        let c_val = cs.get_value(self.c)?;
        let d_val = cs.get_value(self.d)?;
        
        let mut result = self.cnst;
        let products = [
            &[a_val, b_val, self.a_mul_b_coef][..], &[a_val, c_val, self.a_mul_c_coef][..], 
            &[a_val, self.a_linear_coef][..], &[b_val, self.b_linear_coef][..], 
            &[c_val, self.c_linear_coef][..], &[d_val, self.d_linear_coef][..]
        ];
        for chunk in std::array::IntoIter::new(products) {
            let mut tmp = chunk[0].clone();
            for multiplier in chunk.iter().skip(1) {
                tmp.mul_assign(&multiplier);
            }
            result.add_assign(&tmp);
        }

        Ok(result)
    }
}


#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct OrderedVariablePair {
    first: Variable,
    second: Variable
}

impl OrderedVariablePair {
    pub fn new(a: Variable, b: Variable) -> Self {
        let lt_flag = match (a.get_unchecked(), b.get_unchecked()) {
            (Index::Aux(a_idx), Index::Aux(b_idx)) => a_idx < b_idx,
            _ => unimplemented!(),
        };
        let (first, second) = if lt_flag { (a, b) } else { (b, a) };
        OrderedVariablePair { first, second }
    }

    pub fn get_associate(&self, elem: Variable) -> Variable {
        if elem == self.first {
            self.second
        } else if elem == self.second {
            self.first
        } else {
            unreachable!("element should be equal to at least one of components of Variable pair")
        }
    }
}


// module containing amplified version of linear combination that supports both linear and multiplicative terms
pub struct AmplifiedLinearCombination<E: Engine> {
    value: Option<E::Fr>,
    linear_terms: IndexMap<Variable, E::Fr>,
    quadratic_terms: IndexMap<OrderedVariablePair, E::Fr>,
    constant: E::Fr,
}

impl<E: Engine> std::fmt::Debug for AmplifiedLinearCombination<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AmplifiedLinearCombination")
            .field("value", &self.value)
            .field("number of terms", &(self.linear_terms.len() + self.quadratic_terms.len()))
            .field("linear_terms", &self.linear_terms)
            .field("quadratic_terms", &self.quadratic_terms)
            .field("constant", &self.constant)
            .finish()
    }
}

impl<E: Engine> From<AllocatedNum<E>> for AmplifiedLinearCombination<E> {
    fn from(num: AllocatedNum<E>) -> AmplifiedLinearCombination<E> {
        Self {
            value: num.value,
            linear_terms: IndexMap::from([(num.variable, E::Fr::one())]),
            quadratic_terms: IndexMap::new(),
            constant: E::Fr::zero()
        }
    }
}

impl<E: Engine> From<Num<E>> for AmplifiedLinearCombination<E> {
    fn from(num: Num<E>) -> AmplifiedLinearCombination<E> {
        match num {
            Num::Variable(allocated) => {
                Self::from(allocated)
            },
            Num::Constant(constant) => {
                Self {
                    value: Some(constant),
                    linear_terms: IndexMap::new(),
                    quadratic_terms: IndexMap::new(),
                    constant: constant
                }
            }
        }
    }
}

impl<E: Engine> Clone for AmplifiedLinearCombination<E> {
    fn clone(&self) -> Self {
        Self {
            value: self.value,
            linear_terms: self.linear_terms.clone(),
            quadratic_terms: self.quadratic_terms.clone(),
            constant: self.constant,
        }
    }
}

impl<E: Engine> AmplifiedLinearCombination<E> {
    pub fn zero() -> Self {
        Self {
            value: Some(E::Fr::zero()),
            linear_terms: IndexMap::new(),
            quadratic_terms: IndexMap::new(),
            constant: E::Fr::zero(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn len(&self) -> usize {
        self.linear_terms.len() + self.quadratic_terms.len()
    }

    pub fn get_value(&self) -> Option<E::Fr> {
        self.value
    }

    pub fn scale(&mut self, coeff: &E::Fr) {
        if coeff.is_zero() {
            self.linear_terms.clear();
            self.quadratic_terms.clear();
            self.constant = E::Fr::zero();
            self.value = Some(E::Fr::zero());
            return;
        }

        if let Some(ref mut val) = self.value.as_mut() {
            val.mul_assign(&coeff);
        }

        for (_k, v) in self.linear_terms.iter_mut() {
            v.mul_assign(&coeff);
        }
        for (_k, v) in self.quadratic_terms.iter_mut() {
            v.mul_assign(&coeff);
        }

        self.constant.mul_assign(&coeff);
    }

    pub fn add_assign_number_with_coeff(&mut self, number: &Num<E>, coeff: E::Fr) {
        match number {
            Num::Variable(ref allocated_number) => {
                self.add_assign_variable_with_coeff(allocated_number, coeff);
            },
            Num::Constant(constant) => {
                let mut res = coeff;
                res.mul_assign(&constant);
                self.add_assign_constant(res);
            } 
        }
    }

    pub fn add_assign_variable_with_coeff(&mut self, variable: &AllocatedNum<E>, coeff: E::Fr) {
        if coeff.is_zero() {
            return;
        }

        let newval = match (self.value, variable.get_value()) {
            (Some(mut curval), Some(val)) => {
                let mut tmp = val;
                tmp.mul_assign(&coeff);

                curval.add_assign(&tmp);

                Some(curval)
            },
            _ => None
        };

        self.value = newval;
        self.linear_terms.entry(variable.get_variable()).and_modify(|fr| fr.add_assign(&coeff)).or_insert(coeff);
    }

    pub fn add_assign_term(&mut self, term: &Term<E>)
    {
        if term.is_constant() {
            self.add_assign_constant(term.get_constant_value());
            return;
        }

        // otherwise add constant and scaled num separately
        self.add_assign_constant(term.constant_term);
        self.add_assign_number_with_coeff(&term.num, term.coeff);
    }

    pub fn sub_assign_term(&mut self, term: &Term<E>)
    {
        let mut negated_term = term.clone();
        negated_term.negate();
        self.add_assign_term(&negated_term);
    }

    pub fn add_assign_term_with_coeff(&mut self, term: &Term<E>, coeff: E::Fr)
    {
        let mut scaled_term = term.clone();
        scaled_term.scale(&coeff);
        self.add_assign_term(&scaled_term)
    }

    pub fn sub_assign_term_with_coeff(&mut self, term: &Term<E>, coeff: E::Fr)
    {
        let mut scaled_term = term.clone();
        scaled_term.scale(&coeff);
        scaled_term.negate();
        self.add_assign_term(&scaled_term)
    }
   
    pub fn add_assign_constant(&mut self, coeff: E::Fr){
        if coeff.is_zero() {
            return;
        }

        let newval = match self.value {
            Some(mut curval) => {
                curval.add_assign(&coeff);

                Some(curval)
            },
            None => {
                None
            }
        };

        self.value = newval;
        self.constant.add_assign(&coeff);
    }

    pub fn sub_assign_constant(&mut self, coeff: E::Fr) {
        let mut c = coeff;
        c.negate();
        self.add_assign_constant(c);
    }

    pub fn add_assign_product_of_terms(&mut self, a: &Term<E>, b: &Term<E>) {
        self.add_assign_product_of_terms_with_coeff(a, b, E::Fr::one())
    }

    pub fn sub_assign_product_of_terms(&mut self, a: &Term<E>, b: &Term<E>) {
        let mut minus_one = E::Fr::one();
        minus_one.negate();
        self.add_assign_product_of_terms_with_coeff(a, b, minus_one)
    }

    pub fn add_assign_product_of_terms_with_coeff(&mut self, a: &Term<E>, b: &Term<E>, coeff: E::Fr) {
        let mut a_scaled = a.clone();
        a_scaled.scale(&coeff);

        if let Some(fr) = a_scaled.try_into_constant_value() {
            self.add_assign_term_with_coeff(b, fr);    
        } else if let Some(fr) = b.try_into_constant_value() {
            self.add_assign_term_with_coeff(&a_scaled, fr);
        } else {
            self.value = match (self.value, a_scaled.get_value(), b.get_value()) {
                (Some(lc_val), Some(a_val), Some(b_val)) => {
                    let mut tmp = a_val;
                    tmp.mul_assign(&b_val);
                    tmp.add_assign(&lc_val);
                    Some(tmp)
                },
                _ => None,
            };

            // let a = p_1 * x + c_1, b = p_2 * y + c_2, then:
            // a * b = (p_1 * x + c_1) * (p_2 * y + c_2) = 
            // = p_1 * p_2 * x * y + p_1 * c_2 * x + p_2 * c_1 * y + c_1 * c_2
            let (x, p_1, c_1) = a_scaled.unpack(); 
            let (y, p_2, c_2) = b.unpack();
        
            // add p_1 * p_2 * x * y
            let var_pair = OrderedVariablePair::new(x, y);
            let mut coeff = p_1.clone();
            coeff.mul_assign(&p_2); 
            self.quadratic_terms.entry(var_pair).and_modify(|fr| fr.add_assign(&coeff)).or_insert(coeff);

            // add p_1 * c_2 * x
            let mut coeff = p_1.clone();
            coeff.mul_assign(&c_2); 
            self.linear_terms.entry(x).and_modify(|fr| fr.add_assign(&coeff)).or_insert(coeff);

            // add p_2 * c_1 * y
            let mut coeff = c_1.clone();
            coeff.mul_assign(&p_2); 
            self.linear_terms.entry(y).and_modify(|fr| fr.add_assign(&coeff)).or_insert(coeff);

            // add c_1 * c_2
            let mut coeff = c_1.clone();
            coeff.mul_assign(&c_2);
            self.constant.add_assign(&coeff); 
        }
    }

    pub fn normalize(&mut self) {
        self.linear_terms.retain(|&_, v| !v.is_zero());
        self.quadratic_terms.retain(|&_, v| !v.is_zero());
    }

    fn get_linear_terms_only_variables(&self) -> IndexSet<Variable> {
        let mut quad = IndexSet::new();
        for var_pair in self.quadratic_terms.keys() {
            quad.insert(var_pair.first);
            quad.insert(var_pair.second);
        }
        
        let mut lin = IndexSet::new();
        for var in self.linear_terms.keys() {
            lin.insert(var.clone());
        }

        lin.difference(&quad).cloned().collect::<IndexSet<_>>()
    }

    #[track_caller]
    pub fn enforce_zero<CS: ConstraintSystem<E>>(mut self, cs: &mut CS) -> Result<usize, SynthesisError> {
        if let Some(value) = self.get_value() {
            assert!(value.is_zero(), "ALC is not actually zero with value {}", value);
        }
        assert!(CS::Params::CAN_ACCESS_NEXT_TRACE_STEP);
        assert!(CS::Params::STATE_WIDTH == REQUIRED_STATE_WIDTH);
        assert!(is_selector_specialized_gate::<E, CS>());
        self.normalize();
        if self.len() == 0 { return Ok(0); }
 
        let mut linear_terms_only_vars = self.get_linear_terms_only_variables();
        let flattened_quad_releations : Vec<(OrderedVariablePair, E::Fr)> = self.quadratic_terms.into_iter().collect();
        let flattened_arr_len = flattened_quad_releations.len();
        let mut arr_indexer : IndexMap<Variable, usize> = IndexMap::new();
        let mut gate_templates : Vec<GateConstructorHelper<E>> = vec![];
        let mut num_gates_allocated : usize = 0;

        for i in 0..flattened_arr_len {
            let (var_pair_i, fr_i) = flattened_quad_releations[i].clone();
            let mut insertion_flag = true;
            for var in [var_pair_i.first, var_pair_i.second].iter() {
                if let Some(j) = arr_indexer.get(var) {
                    let (var_pair_j, fr_j) = flattened_quad_releations[*j].clone();
                    
                    // construct gate from two multiplications:
                    let a = *var;
                    let b = var_pair_i.get_associate(a);
                    let c = var_pair_j.get_associate(a);
                    let gate = GateConstructorHelper::new_for_pair_of_muls(cs, a, b, c, fr_i, fr_j);
                    gate_templates.push(gate);
                    arr_indexer.remove(&var_pair_j.first);
                    arr_indexer.remove(&var_pair_j.second);
                    insertion_flag = false;
                    break;
                }
            }
            for var in [var_pair_i.first, var_pair_i.second].iter() {
                if insertion_flag {
                    arr_indexer.insert(*var, i);
                }
            }
        }
        
        let unconsumed_idxes = IndexSet::<usize>::from_iter(arr_indexer.into_values());
        for i in unconsumed_idxes {
            let (var_pair, fr) = flattened_quad_releations[i];
            let gate_template = GateConstructorHelper::new_for_mul(cs, var_pair, fr);
            gate_templates.push(gate_template);
        }

        let mut is_first = true;
        let mut idx = 0;
        let mut next_trace_step_var = AllocatedNum::zero(cs).get_variable();
        loop {
            let mut gate_template = gate_templates.get(idx).cloned().unwrap_or_else(|| GateConstructorHelper::new(cs));
            idx += 1;
            if is_first {
                gate_template.add_constant_term(self.constant);
                is_first = false;
            }
            else {
                gate_template.add_next_trace_step_term(next_trace_step_var);
            }
            gate_template.add_linear_coefficients_for_bound_variables(&mut self.linear_terms);
            gate_template.add_linear_coefficients_for_free_variables(
                &mut self.linear_terms, &mut linear_terms_only_vars
            );

            let is_final = (idx >= gate_templates.len()) && (linear_terms_only_vars.is_empty());
            if is_final {
                gate_template.set_finality_flag();
            }
            else {
                // we manually allocate the new variable
                let may_be_new_intermediate_value = gate_template.evaluate_next_trace_step_value(cs);
                next_trace_step_var = cs.alloc(|| { may_be_new_intermediate_value })?;
            }
            gate_template.materialize(cs)?;
            num_gates_allocated += 1;

            if is_final {
                break
            }
        } 
        Ok(num_gates_allocated)
    }

    pub fn into_num<CS: ConstraintSystem<E>>(self, cs: &mut CS) -> Result<Num<E>, SynthesisError> {
        let (res, _num_gates) = self.into_num_ext(cs)?;
        Ok(res)
    }

    pub fn into_num_ext<CS: ConstraintSystem<E>>(mut self, cs: &mut CS) -> Result<(Num<E>, usize), SynthesisError> {
        self.normalize();
        let value = self.get_value();
        let num = if self.len() == 0 {
            Num::Constant(value.unwrap())
        } else {
            let var = AllocatedNum::alloc(cs, || value.grab())?;
            Num::Variable(var)
        };

        let mut minus_one = E::Fr::one();
        minus_one.negate();
        self.add_assign_number_with_coeff(&num, minus_one);
        let num_gates = self.enforce_zero(cs)?;
        
        Ok((num, num_gates))
    }
}
