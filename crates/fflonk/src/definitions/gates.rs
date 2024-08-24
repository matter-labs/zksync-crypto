use super::*;

pub const STD_MAIN_GATE_NAME: &'static str = "standard main gate";
pub const STD_MAIN_GATE_NAME_WITH_DNEXT: &'static str = "main gate of width 4 with D_next";
pub const SELECTOR_OPTIMIZED_MAIN_GATE_NAME: &'static str = "main gate of width 4 with D_next and selector optimization";
pub const SBOX_CUSTOM_GATE_NAME: &'static str = "Alpha=5 custom gate for Rescue/Poseidon";
pub const TEST_BIT_CUSTOM_GATE_NAME: &'static str = "Test bit gate on A";

pub const fn supported_main_gates() -> [&'static str; 3] {
    [STD_MAIN_GATE_NAME, STD_MAIN_GATE_NAME_WITH_DNEXT, SELECTOR_OPTIMIZED_MAIN_GATE_NAME]
}
pub const fn supported_custom_gates() -> [&'static str; 2] {
    [SBOX_CUSTOM_GATE_NAME, TEST_BIT_CUSTOM_GATE_NAME]
}

pub fn sorted_gates_from_circuit_definitions<E: Engine, C: Circuit<E>>() -> Vec<Box<dyn GateInternal<E>>> {
    let mut sorted_gates = vec![];
    // main gate first
    // custom gate is next if there is one
    let gates = C::declare_used_gates().unwrap();
    if supported_main_gates().contains(&gates[0].name()) {
        sorted_gates.push(gates[0].clone());
    }
    if gates.len() > 1 {
        if supported_main_gates().contains(&gates[1].name()) {
            sorted_gates.push(gates[1].clone());
        }
    }

    if supported_custom_gates().contains(&gates[0].name()) {
        sorted_gates.push(gates[0].clone());
    }
    if gates.len() > 1 {
        if supported_custom_gates().contains(&gates[1].name()) {
            sorted_gates.push(gates[1].clone());
        }
    }

    sorted_gates
}

#[inline(always)]
pub fn has_custom_gate<E: Engine, C: Circuit<E>>() -> bool {
    C::declare_used_gates().unwrap().len() > 1
}

pub struct TraceAndGateMonomials<'a, E: Engine> {
    pub trace_monomials: Vec<&'a Polynomial<E::Fr, Coefficients>>,
    pub main_gate_quotient_monomial: Polynomial<E::Fr, Coefficients>,
    pub custom_gate_quotient_monomial: Option<Polynomial<E::Fr, Coefficients>>,
    pub num_state_polys: usize,
    pub num_witness_polys: usize,
    pub main_gate: Box<dyn GateInternal<E>>,
    pub custom_gate: Option<Box<dyn GateInternal<E>>>,
    pub num_polys: usize,
}

impl<'a, E: Engine> TraceAndGateMonomials<'a, E> {
    pub fn flatten(&self) -> Vec<&Polynomial<E::Fr, Coefficients>> {
        let mut flattened = vec![];

        flattened.extend(self.trace_monomials.iter().map(|p| p));
        flattened.push(&self.main_gate_quotient_monomial);
        self.custom_gate_quotient_monomial.as_ref().map(|m| flattened.push(m));

        flattened
    }

    pub fn interpolation_size(&self) -> usize {
        self.num_polys.next_power_of_two()
    }
}

#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize, PartialEq, Eq)]
#[serde(bound = "")]
pub struct TraceAndGateEvaluations<F: PrimeField> {
    pub trace_evaluations_at_z: Vec<F>,
    pub trace_evaluations_at_z_omega: Option<Vec<F>>,
    pub main_gate_quotient_at_z: F,
    pub main_gate_quotient_at_z_omega: Option<F>,
    pub custom_gate_quotient_at_z: Option<F>,
    pub custom_gate_quotient_at_z_omega: Option<F>,
    pub num_polys: usize,
    has_custom_gate: bool,
}

#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct TraceAndGateEvaluationOffsets {
    pub trace_evaluations_at_z: usize,
    pub trace_evaluations_at_z_omega: usize,
    pub main_gate_quotient_at_z_omega: usize,
    pub custom_gate_quotient_at_z_omega: usize,

    pub a_term: usize,
    pub b_term: usize,
    pub c_term: usize,
    pub d_term: usize,
}

impl<F: PrimeField> Flatten<F> for TraceAndGateEvaluations<F> {
    fn flatten(&self) -> (Vec<F>, Vec<F>) {
        let mut evals = vec![];
        let mut shifted_evals = vec![];

        evals.extend_from_slice(&self.trace_evaluations_at_z);

        if self.requires_opening_at_shifted_point() {
            self.trace_evaluations_at_z_omega.as_ref().map(|values| shifted_evals.extend_from_slice(values));
            self.main_gate_quotient_at_z_omega.map(|value| shifted_evals.push(value));
            self.custom_gate_quotient_at_z_omega.map(|value| shifted_evals.push(value));

            assert_eq!(evals.len(), shifted_evals.len());
        }

        (evals, shifted_evals)
    }
}

impl<F: PrimeField> TraceAndGateEvaluations<F> {
    pub fn requires_opening_at_shifted_point(&self) -> bool {
        self.main_gate_quotient_at_z_omega.is_some()
    }

    pub fn interpolation_size(&self) -> usize {
        self.num_polys.next_power_of_two()
    }
}

pub(crate) fn compute_quotient_of_main_gate_at_z<F: PrimeField>(
    main_gate_name: &'static str,
    trace_evals_at_z: &TraceAndGateEvaluations<F>,
    setup_evals_at_z: &SetupEvaluations<F>,
    public_inputs_at_z: F,
) -> F {
    if main_gate_name == STD_MAIN_GATE_NAME {
        compute_quotient_of_std_main_gate_at_z(main_gate_name, trace_evals_at_z, setup_evals_at_z, public_inputs_at_z)
    } else if main_gate_name == STD_MAIN_GATE_NAME_WITH_DNEXT {
        compute_quotient_of_std_main_gate_with_dnext_at_z(main_gate_name, trace_evals_at_z, setup_evals_at_z, public_inputs_at_z)
    } else if main_gate_name == SELECTOR_OPTIMIZED_MAIN_GATE_NAME {
        compute_quotient_of_selector_optimized_main_gate_at_z(main_gate_name, trace_evals_at_z, setup_evals_at_z, public_inputs_at_z)
    } else {
        unreachable!("unknown main gate type");
    }
}

pub(crate) fn compute_quotient_of_main_gate_at_z_flattened<F: PrimeField>(main_gate_name: &'static str, evaluations: &[F], public_inputs_at_z: F, offsets: &EvaluationOffsets) -> F {
    if main_gate_name == STD_MAIN_GATE_NAME {
        compute_quotient_of_std_main_gate_at_z_flattend(main_gate_name, evaluations, public_inputs_at_z, offsets)
    } else if main_gate_name == STD_MAIN_GATE_NAME_WITH_DNEXT {
        compute_quotient_of_std_main_gate_with_dnext_at_z_flattend(main_gate_name, evaluations, public_inputs_at_z, offsets)
    } else if main_gate_name == SELECTOR_OPTIMIZED_MAIN_GATE_NAME {
        compute_quotient_of_selector_optimized_main_gate_at_z_flattened(main_gate_name, evaluations, public_inputs_at_z, offsets)
    } else {
        unreachable!("unknown main gate type");
    }
}

fn compute_quotient_of_std_main_gate_at_z<F: PrimeField>(name: &'static str, trace_evals_at_z: &TraceAndGateEvaluations<F>, setup_evals_at_z: &SetupEvaluations<F>, public_inputs_at_z: F) -> F {
    assert_eq!(name, STD_MAIN_GATE_NAME);
    // a, b, c
    let mut trace_evals_iter = trace_evals_at_z.trace_evaluations_at_z.iter().cloned();
    let a = trace_evals_iter.next().unwrap();
    let b = trace_evals_iter.next().unwrap();
    let c = trace_evals_iter.next().unwrap();

    // q_a, q_b, q_c
    let mut setup_evals_iter = setup_evals_at_z.gate_setups_at_z.iter().cloned();
    let q_a = setup_evals_iter.next().unwrap();
    let q_b = setup_evals_iter.next().unwrap();
    let q_c = setup_evals_iter.next().unwrap();
    let q_m = setup_evals_iter.next().unwrap();
    let q_const = setup_evals_iter.next().unwrap();

    let mut sum = q_const;
    sum.add_assign(&public_inputs_at_z);
    for (selector, value) in [q_a, q_b, q_c].into_iter().zip([a, b, c]) {
        let mut tmp = selector;
        tmp.mul_assign(&value);
        sum.add_assign(&tmp);
    }

    let mut tmp = q_m;
    tmp.mul_assign(&a);
    tmp.mul_assign(&b);
    sum.add_assign(&tmp);

    if trace_evals_at_z.has_custom_gate {
        assert!(setup_evals_at_z.gate_selectors_at_z.len() > 1);
        sum.mul_assign(&setup_evals_at_z.gate_selectors_at_z[0]);
    }

    sum
}

fn compute_quotient_of_std_main_gate_with_dnext_at_z<F: PrimeField>(
    name: &'static str,
    trace_evals_at_z: &TraceAndGateEvaluations<F>,
    setup_evals_at_z: &SetupEvaluations<F>,
    public_inputs_at_z: F,
) -> F {
    assert_eq!(name, STD_MAIN_GATE_NAME_WITH_DNEXT);
    // a, b, c, d
    let mut trace_evals_iter = trace_evals_at_z.trace_evaluations_at_z.iter().cloned();
    let a = trace_evals_iter.next().unwrap();
    let b = trace_evals_iter.next().unwrap();
    let c = trace_evals_iter.next().unwrap();
    let d = trace_evals_iter.next().unwrap();

    // q_a, q_b, q_c, q_d
    let mut setup_evals_iter = setup_evals_at_z.gate_setups_at_z.iter().cloned();
    let q_a = setup_evals_iter.next().unwrap();
    let q_b = setup_evals_iter.next().unwrap();
    let q_c = setup_evals_iter.next().unwrap();
    let q_d = setup_evals_iter.next().unwrap();
    let q_m = setup_evals_iter.next().unwrap();
    let q_const = setup_evals_iter.next().unwrap();
    let q_d_next = setup_evals_iter.next().unwrap();

    let mut sum = q_const;
    sum.add_assign(&public_inputs_at_z);
    for (selector, value) in [q_a, q_b, q_c, q_d].into_iter().zip([a, b, c, d]) {
        let mut tmp = selector;
        tmp.mul_assign(&value);
        sum.add_assign(&tmp);
    }
    // D_next is d*w, prover already shows all trace polys at z*w
    let mut tmp = trace_evals_at_z.trace_evaluations_at_z_omega.as_ref().expect("trace evaluations at z oemga")[3];
    tmp.mul_assign(&q_d_next);
    sum.add_assign(&tmp);

    let mut tmp = q_m;
    tmp.mul_assign(&a);
    tmp.mul_assign(&b);
    sum.add_assign(&tmp);

    if trace_evals_at_z.has_custom_gate {
        assert!(setup_evals_at_z.gate_selectors_at_z.len() > 1);
        sum.mul_assign(&setup_evals_at_z.gate_selectors_at_z[0]);
    }

    sum
}

fn compute_quotient_of_std_main_gate_at_z_flattend<F: PrimeField>(name: &'static str, evaluations: &[F], public_inputs_at_z: F, offsets: &EvaluationOffsets) -> F {
    assert_eq!(name, STD_MAIN_GATE_NAME);
    // a, b, c
    let mut trace_evals_iter = evaluations[offsets.trace.trace_evaluations_at_z..offsets.trace.trace_evaluations_at_z + 3].iter().cloned();
    let a = trace_evals_iter.next().unwrap();
    let b = trace_evals_iter.next().unwrap();
    let c = trace_evals_iter.next().unwrap();
    assert!(trace_evals_iter.next().is_none());

    // q_a, q_b, q_c, qAB, qConst
    let mut setup_evals_iter = evaluations[offsets.setup.gate_setups_at_z..offsets.setup.gate_setups_at_z + 5].iter().cloned();
    let q_a = setup_evals_iter.next().unwrap();
    let q_b = setup_evals_iter.next().unwrap();
    let q_c = setup_evals_iter.next().unwrap();
    let q_m = setup_evals_iter.next().unwrap();
    let q_const = setup_evals_iter.next().unwrap();
    assert!(setup_evals_iter.next().is_none());

    let mut sum = q_const;
    sum.add_assign(&public_inputs_at_z);
    for (selector, value) in [q_a, q_b, q_c].into_iter().zip([a, b, c]) {
        let mut tmp = selector;
        tmp.mul_assign(&value);
        sum.add_assign(&tmp);
    }

    let mut tmp = q_m;
    tmp.mul_assign(&a);
    tmp.mul_assign(&b);
    sum.add_assign(&tmp);

    // main gate 3 variables doesn't need selector as there is no custom gate in this approach
    sum
}

fn compute_quotient_of_std_main_gate_with_dnext_at_z_flattend<F: PrimeField>(name: &'static str, evaluations: &[F], public_inputs_at_z: F, offsets: &EvaluationOffsets) -> F {
    assert_eq!(name, STD_MAIN_GATE_NAME_WITH_DNEXT);
    // a, b, c, d
    let mut trace_evals_iter = evaluations[offsets.trace.trace_evaluations_at_z..offsets.trace.trace_evaluations_at_z + 4].iter().cloned();
    let a = trace_evals_iter.next().unwrap();
    let b = trace_evals_iter.next().unwrap();
    let c = trace_evals_iter.next().unwrap();
    let d = trace_evals_iter.next().unwrap();
    assert!(trace_evals_iter.next().is_none());

    // q_a, q_b, q_c, q_d
    let mut setup_evals_iter = evaluations[offsets.setup.gate_setups_at_z..offsets.setup.gate_setups_at_z + 7].iter().cloned();
    let q_a = setup_evals_iter.next().unwrap();
    let q_b = setup_evals_iter.next().unwrap();
    let q_c = setup_evals_iter.next().unwrap();
    let q_d = setup_evals_iter.next().unwrap();
    let q_m = setup_evals_iter.next().unwrap();
    let q_const = setup_evals_iter.next().unwrap();
    let q_d_next = setup_evals_iter.next().unwrap();
    assert!(setup_evals_iter.next().is_none());

    let mut sum = q_const;
    sum.add_assign(&public_inputs_at_z);
    for (selector, value) in [q_a, q_b, q_c, q_d].into_iter().zip([a, b, c, d]) {
        let mut tmp = selector;
        tmp.mul_assign(&value);
        sum.add_assign(&tmp);
    }
    // D_next is d*w, prover already shows all trace polys at z*w
    let mut tmp = d;
    tmp.mul_assign(&q_d_next);
    sum.add_assign(&tmp);

    let mut tmp = q_m;
    tmp.mul_assign(&a);
    tmp.mul_assign(&b);
    sum.add_assign(&tmp);

    if offsets.has_custom_gate {
        sum.mul_assign(&evaluations[offsets.setup.gate_selectors_at_z]);
    }

    sum
}

fn compute_quotient_of_selector_optimized_main_gate_at_z<F: PrimeField>(
    name: &'static str,
    trace_evals_at_z: &TraceAndGateEvaluations<F>,
    setup_evals_at_z: &SetupEvaluations<F>,
    public_inputs_at_z: F,
) -> F {
    assert_eq!(name, SELECTOR_OPTIMIZED_MAIN_GATE_NAME);
    // a, b, c, d
    let mut trace_evals_iter = trace_evals_at_z.trace_evaluations_at_z.iter().cloned();
    let a = trace_evals_iter.next().unwrap();
    let b = trace_evals_iter.next().unwrap();
    let c = trace_evals_iter.next().unwrap();
    let d = trace_evals_iter.next().unwrap();

    // q_a, q_b, q_c, q_d
    let mut setup_evals_iter = setup_evals_at_z.gate_setups_at_z.iter().cloned();
    let q_a = setup_evals_iter.next().unwrap();
    let q_b = setup_evals_iter.next().unwrap();
    let q_c = setup_evals_iter.next().unwrap();
    let q_d = setup_evals_iter.next().unwrap();
    let q_m_ab = setup_evals_iter.next().unwrap();
    let q_m_ac = setup_evals_iter.next().unwrap();
    let q_const = setup_evals_iter.next().unwrap();
    let q_d_next = setup_evals_iter.next().unwrap();

    let mut sum = q_const;
    sum.add_assign(&public_inputs_at_z);
    for (selector, value) in [q_a, q_b, q_c, q_d].into_iter().zip([a, b, c, d]) {
        let mut tmp = selector;
        tmp.mul_assign(&value);
        sum.add_assign(&tmp);
    }

    // D_next is d*w, prover already shows all trace polys at z*w
    let mut tmp = trace_evals_at_z.trace_evaluations_at_z_omega.as_ref().expect("trace evaluations at z omega")[3];
    tmp.mul_assign(&q_d_next);
    sum.add_assign(&tmp);

    let mut tmp = q_m_ab;
    tmp.mul_assign(&a);
    tmp.mul_assign(&b);
    sum.add_assign(&tmp);

    let mut tmp = q_m_ac;
    tmp.mul_assign(&a);
    tmp.mul_assign(&c);
    sum.add_assign(&tmp);

    if trace_evals_at_z.has_custom_gate {
        assert!(setup_evals_at_z.gate_selectors_at_z.len() > 1);
        sum.mul_assign(&setup_evals_at_z.gate_selectors_at_z[0]);
    }
    sum
}

fn compute_quotient_of_selector_optimized_main_gate_at_z_flattened<F: PrimeField>(name: &'static str, evaluations: &[F], public_inputs_at_z: F, offsets: &EvaluationOffsets) -> F {
    assert_eq!(name, SELECTOR_OPTIMIZED_MAIN_GATE_NAME);
    // a, b, c, d
    let mut trace_evals_iter = evaluations[offsets.trace.trace_evaluations_at_z..offsets.trace.trace_evaluations_at_z + 4].iter().cloned();
    let a = trace_evals_iter.next().unwrap();
    let b = trace_evals_iter.next().unwrap();
    let c = trace_evals_iter.next().unwrap();
    let d = trace_evals_iter.next().unwrap();
    assert!(trace_evals_iter.next().is_none());

    // q_a, q_b, q_c, q_d
    let mut setup_evals_iter = evaluations[offsets.setup.gate_setups_at_z..offsets.setup.gate_setups_at_z + 8].iter().cloned();
    let q_a = setup_evals_iter.next().unwrap();
    let q_b = setup_evals_iter.next().unwrap();
    let q_c = setup_evals_iter.next().unwrap();
    let q_d = setup_evals_iter.next().unwrap();
    let q_m_ab = setup_evals_iter.next().unwrap();
    let q_m_ac = setup_evals_iter.next().unwrap();
    let q_const = setup_evals_iter.next().unwrap();
    let q_d_next = setup_evals_iter.next().unwrap();
    assert!(setup_evals_iter.next().is_none());

    let mut sum = q_const;
    sum.add_assign(&public_inputs_at_z);
    for (selector, value) in [q_a, q_b, q_c, q_d].into_iter().zip([a, b, c, d]) {
        let mut tmp = selector;
        tmp.mul_assign(&value);
        sum.add_assign(&tmp);
    }

    // D_next is d*w, prover already shows all trace polys at z*w
    let mut tmp = evaluations[offsets.trace.trace_evaluations_at_z_omega + 3];
    tmp.mul_assign(&q_d_next);
    sum.add_assign(&tmp);

    let mut tmp = q_m_ab;
    tmp.mul_assign(&a);
    tmp.mul_assign(&b);
    sum.add_assign(&tmp);

    let mut tmp = q_m_ac;
    tmp.mul_assign(&a);
    tmp.mul_assign(&c);
    sum.add_assign(&tmp);

    if offsets.has_custom_gate {
        sum.mul_assign(&evaluations[offsets.setup.gate_selectors_at_z]);
    }
    sum
}

pub(crate) fn compute_quotient_of_custom_gate_at_z<F: PrimeField>(custom_gate_name: &'static str, trace_evals_at_z: &TraceAndGateEvaluations<F>, setup_evals_at_z: &SetupEvaluations<F>) -> F {
    if custom_gate_name == SBOX_CUSTOM_GATE_NAME {
        compute_quotient_of_sbox_custom_gate_at_z(custom_gate_name, trace_evals_at_z, setup_evals_at_z)
    } else if custom_gate_name == TEST_BIT_CUSTOM_GATE_NAME {
        compute_quotient_of_test_bit_custom_gate_at_z(custom_gate_name, trace_evals_at_z, setup_evals_at_z)
    } else {
        unreachable!("unknown custom gate");
    }
}

pub(crate) fn compute_quotient_of_custom_gate_at_z_flattened<F: PrimeField>(custom_gate_name: &'static str, evaluations: &[F], offsets: &EvaluationOffsets) -> F {
    if custom_gate_name == SBOX_CUSTOM_GATE_NAME {
        compute_quotient_of_sbox_custom_gate_at_z_flattened(custom_gate_name, evaluations, offsets)
    } else if custom_gate_name == TEST_BIT_CUSTOM_GATE_NAME {
        compute_quotient_of_test_bit_custom_gate_at_z_flattened(custom_gate_name, evaluations, offsets)
    } else {
        unreachable!("unknown custom gate");
    }
}

fn compute_quotient_of_sbox_custom_gate_at_z<F: PrimeField>(custom_gate_name: &'static str, trace_evals_at_z: &TraceAndGateEvaluations<F>, setup_evals_at_z: &SetupEvaluations<F>) -> F {
    assert_eq!(custom_gate_name, SBOX_CUSTOM_GATE_NAME);
    // x, x^2, x^4, x^5
    // a, b, c, d
    let mut trace_evals_iter = trace_evals_at_z.trace_evaluations_at_z.iter().cloned();
    let a = trace_evals_iter.next().unwrap();
    let b = trace_evals_iter.next().unwrap();
    let c = trace_evals_iter.next().unwrap();
    let d = trace_evals_iter.next().unwrap();

    let mut sum = F::zero();

    let mut a_square_minus_b = a.clone();
    a_square_minus_b.square();
    a_square_minus_b.sub_assign(&b);
    sum.add_assign(&a_square_minus_b);

    // b^2 - c = 0
    let mut b_square_minus_c = b.clone();
    b_square_minus_c.square();
    b_square_minus_c.sub_assign(&c);
    sum.add_assign(&b_square_minus_c);

    let mut a_c_d = c;
    a_c_d.mul_assign(&a);
    a_c_d.sub_assign(&d);
    sum.add_assign(&a_c_d);
    assert!(setup_evals_at_z.gate_selectors_at_z.len() > 1);
    sum.mul_assign(&setup_evals_at_z.gate_selectors_at_z[1]);

    sum
}

fn compute_quotient_of_sbox_custom_gate_at_z_flattened<F: PrimeField>(custom_gate_name: &'static str, evaluations: &[F], offsets: &EvaluationOffsets) -> F {
    assert_eq!(custom_gate_name, SBOX_CUSTOM_GATE_NAME);
    // x, x^2, x^4, x^5
    // a, b, c, d
    let mut trace_evals_iter = evaluations[offsets.trace.trace_evaluations_at_z..offsets.trace.trace_evaluations_at_z + 4].iter().cloned();
    let a = trace_evals_iter.next().unwrap();
    let b = trace_evals_iter.next().unwrap();
    let c = trace_evals_iter.next().unwrap();
    let d = trace_evals_iter.next().unwrap();
    assert!(trace_evals_iter.next().is_none());

    let mut sum = F::zero();

    let mut a_square_minus_b = a.clone();
    a_square_minus_b.square();
    a_square_minus_b.sub_assign(&b);
    sum.add_assign(&a_square_minus_b);

    // b^2 - c = 0
    let mut b_square_minus_c = b.clone();
    b_square_minus_c.square();
    b_square_minus_c.sub_assign(&c);
    sum.add_assign(&b_square_minus_c);

    let mut a_c_d = c;
    a_c_d.mul_assign(&a);
    a_c_d.sub_assign(&d);
    sum.add_assign(&a_c_d);
    sum.mul_assign(&evaluations[offsets.setup.gate_selectors_at_z + 1]);

    sum
}

fn compute_quotient_of_test_bit_custom_gate_at_z<F: PrimeField>(custom_gate_name: &'static str, trace_evals_at_z: &TraceAndGateEvaluations<F>, setup_evals_at_z: &SetupEvaluations<F>) -> F {
    assert_eq!(custom_gate_name, TEST_BIT_CUSTOM_GATE_NAME);

    // a*(a -1) = 0
    let a = trace_evals_at_z.trace_evaluations_at_z[0];
    let mut sum = a;
    sum.sub_assign(&F::one());
    sum.mul_assign(&a);

    assert!(setup_evals_at_z.gate_selectors_at_z.len() > 1);
    sum.mul_assign(&setup_evals_at_z.gate_selectors_at_z[1]);

    sum
}

fn compute_quotient_of_test_bit_custom_gate_at_z_flattened<F: PrimeField>(custom_gate_name: &'static str, evaluations: &[F], offsets: &EvaluationOffsets) -> F {
    assert_eq!(custom_gate_name, TEST_BIT_CUSTOM_GATE_NAME);

    // a*(a -1) = 0
    let a = evaluations[offsets.trace.trace_evaluations_at_z];
    let mut sum = a;
    sum.sub_assign(&F::one());
    sum.mul_assign(&a);

    sum.mul_assign(&evaluations[offsets.setup.gate_selectors_at_z + 1]);

    sum
}

pub(crate) fn compute_gate_quotients<E: Engine, MG: MainGate<E>>(
    all_gates: Vec<Box<dyn GateInternal<E>>>,
    public_inputs: &[E::Fr],
    domain_size: usize,
    ldes_storage: &mut AssembledPolynomialStorage<E>,
    monomials_storage: &AssembledPolynomialStorageForMonomialForms<E>,
    omegas_bitreversed: &BitReversedOmegas<E::Fr>,
    _omegas_inv_bitreversed: &OmegasInvBitreversed<E::Fr>,
    inverse_divisor_on_coset_lde_natural_ordering: &Polynomial<E::Fr, Values>,
    coset_factor: E::Fr,
    worker: &Worker,
) -> Result<(Polynomial<E::Fr, Coefficients>, Option<Polynomial<E::Fr, Coefficients>>), SynthesisError> {
    assert!(all_gates.len() <= 2);
    let main_gate_quotient_degree = main_gate_quotient_degree(&all_gates);
    let mut gates_iter = all_gates.into_iter();
    let main_gate_as_internal = gates_iter.next().expect("main gate");
    let main_gate = MG::default();
    assert!(&main_gate.clone().into_internal() == &main_gate_as_internal);
    let lde_factor = ldes_storage.lde_factor;

    assert_eq!(inverse_divisor_on_coset_lde_natural_ordering.size(), domain_size * lde_factor);

    let main_gate_challenges = vec![E::Fr::one(); main_gate_as_internal.num_quotient_terms()];
    let mut main_gate_contrib_lde = MainGate::contribute_into_quotient_for_public_inputs(
        &main_gate,
        domain_size,
        public_inputs,
        ldes_storage,
        monomials_storage,
        &main_gate_challenges,
        omegas_bitreversed,
        _omegas_inv_bitreversed,
        worker,
    )?;

    let custom_gate_contrib_monomial = if let Some(custom_gate_as_internal) = gates_iter.next() {
        assert!(gates_iter.next().is_none());
        let key = PolyIdentifier::GateSelector(main_gate_as_internal.name());
        let main_gate_selector_monomial = monomials_storage.gate_selectors.get(&key).unwrap().as_ref();
        let main_gate_selector_lde = main_gate_selector_monomial
            .clone_padded_to_domain()?
            .bitreversed_lde_using_bitreversed_ntt(&worker, lde_factor, omegas_bitreversed, &coset_factor)?;
        main_gate_contrib_lde.mul_assign(worker, &main_gate_selector_lde);

        let custom_gate_challenges = vec![E::Fr::one(); custom_gate_as_internal.num_quotient_terms()];
        let mut custom_gate_contrib_lde = custom_gate_as_internal.contribute_into_quotient(
            domain_size,
            ldes_storage,
            monomials_storage,
            &custom_gate_challenges,
            omegas_bitreversed,
            _omegas_inv_bitreversed,
            worker,
        )?;

        let key = PolyIdentifier::GateSelector(custom_gate_as_internal.name());
        let custom_gate_selector_monomial = monomials_storage.gate_selectors.get(&key).unwrap().as_ref();
        let custom_gate_selector_lde = custom_gate_selector_monomial
            .clone_padded_to_domain()?
            .bitreversed_lde_using_bitreversed_ntt(&worker, lde_factor, omegas_bitreversed, &coset_factor)?;
        custom_gate_contrib_lde.mul_assign(worker, &custom_gate_selector_lde);
        let custom_gate_quotient_degree = custom_gate_as_internal.degree();
        let custom_gate_contrib_monomial = compute_quotient_monomial(
            worker,
            custom_gate_contrib_lde,
            inverse_divisor_on_coset_lde_natural_ordering,
            coset_factor,
            lde_factor,
            domain_size,
            custom_gate_quotient_degree,
        )?;
        assert_eq!(custom_gate_contrib_monomial.size(), custom_gate_quotient_degree * domain_size);

        Some(custom_gate_contrib_monomial)
    } else {
        None
    };

    let main_gate_contrib_monomial = compute_quotient_monomial(
        worker,
        main_gate_contrib_lde,
        inverse_divisor_on_coset_lde_natural_ordering,
        coset_factor,
        lde_factor,
        domain_size,
        main_gate_quotient_degree,
    )?;
    assert_eq!(main_gate_contrib_monomial.size(), main_gate_quotient_degree * domain_size);

    Ok((main_gate_contrib_monomial, custom_gate_contrib_monomial))
}

pub(crate) fn evaluate_trace_and_gate_monomials<E: Engine>(worker: &Worker, monomial_storage: &TraceAndGateMonomials<E>, z: E::Fr, z_omega: E::Fr) -> TraceAndGateEvaluations<E::Fr> {
    let mut evals = TraceAndGateEvaluations::default();
    evals.num_polys = monomial_storage.num_polys;

    evals.trace_evaluations_at_z = monomial_storage.trace_monomials.iter().map(|p| p.evaluate_at(worker, z)).collect();

    evals.main_gate_quotient_at_z = monomial_storage.main_gate_quotient_monomial.evaluate_at(worker, z);

    evals.custom_gate_quotient_at_z = monomial_storage.custom_gate_quotient_monomial.as_ref().map(|p| p.evaluate_at(worker, z));

    if requires_trace_polys_opening_at_shifted_point(monomial_storage.main_gate.clone()) {
        evals.trace_evaluations_at_z_omega = Some(monomial_storage.trace_monomials.iter().map(|p| p.evaluate_at(worker, z_omega)).collect());

        evals.main_gate_quotient_at_z_omega = Some(monomial_storage.main_gate_quotient_monomial.evaluate_at(worker, z_omega));
        evals.custom_gate_quotient_at_z_omega = monomial_storage.custom_gate_quotient_monomial.as_ref().map(|p| p.evaluate_at(worker, z_omega));
    }

    evals
}

pub fn requires_trace_polys_opening_at_shifted_point<E: Engine>(main_gate: Box<dyn GateInternal<E>>) -> bool {
    for el in main_gate.all_queried_polynomials() {
        let (id, dilation) = el.into_id_and_raw_dilation();
        match id {
            PolyIdentifier::VariablesPolynomial(_) => {
                if dilation > 0 {
                    return true;
                }
            }
            _ => (),
        }
    }

    return false;
}
