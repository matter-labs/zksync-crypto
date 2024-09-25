use bellman::plonk::better_better_cs::utils::FieldBinop;
use franklin_crypto::plonk::circuit::bigint::{biguint_to_fe, fe_to_biguint, repr_to_biguint};
use num_bigint::BigUint;
use num_traits::{One, Zero};

use super::*;

pub fn lcm(numbers: &[usize]) -> usize {
    fn gcd(mut a: usize, mut b: usize) -> usize {
        while b != 0 {
            let temp = b;
            b = a % b;
            a = temp;
        }
        a
    }
    let lcm = numbers.iter().fold(1, |acc, &num| acc / gcd(acc, num) * num);
    numbers.iter().for_each(|value| assert_eq!(lcm % value, 0));

    lcm
}

pub(crate) fn compute_max_combined_degree_from_assembly<E: Engine, P: PlonkConstraintSystemParams<E>, MG: MainGate<E>, S: SynthesisMode, C: Circuit<E>>(assembly: &Assembly<E, P, MG, S>) -> usize {
    let has_custom_gate = assembly.sorted_gates.len() > 1;
    let has_lookup = assembly.num_table_lookups > 0 && assembly.tables.len() > 0;
    let main_gate_quotient_degree = main_gate_quotient_degree(&assembly.sorted_gates);
    let custom_gate_quotient_degree = if has_custom_gate { custom_gate_quotient_degree(&assembly.sorted_gates) } else { 0 };
    let copy_permutation_quotient_degree = P::STATE_WIDTH;
    let lookup_quotient_degree = if has_lookup { 2 } else { 0 };

    let (num_setup_polys, num_first_round_polys, num_second_round_polys, _) = num_system_polys_from_assembly::<_, _, _, _, C>(&assembly);

    let domain_size = assembly.n() + 1;
    assert!(domain_size.is_power_of_two());

    let max_combined_degree = [
        (num_setup_polys, 1),
        (num_first_round_polys, main_gate_quotient_degree),
        (num_first_round_polys, custom_gate_quotient_degree),
        (num_second_round_polys, copy_permutation_quotient_degree),
        (num_second_round_polys, lookup_quotient_degree),
    ]
    .into_iter()
    .map(|(num_polys, degree)| num_polys * degree * domain_size + num_polys - 1)
    .max()
    .unwrap();

    max_combined_degree
}

pub(crate) fn combine_monomials<F: PrimeField>(monomials: &[&Polynomial<F, Coefficients>], domain_size: usize) -> Result<Polynomial<F, Coefficients>, SynthesisError> {
    let num_polys = monomials.len();
    assert!(num_polys.is_power_of_two());
    assert!(domain_size.is_power_of_two());

    let full_size = domain_size * num_polys + num_polys - 1;
    let mut combined = vec![F::zero(); full_size];

    // fe(x) = c0 + c1x
    // fo(x) = c2 + c3x
    // f(x) = fe(x^2) + xfo(x^2)
    // f(x) = c0 + c1x^2 + x*(c2 + c3x^2)
    // f(x) = c0 + c2x + c1x^2 + c3x^3

    for (poly_idx, poly) in monomials.iter().enumerate() {
        assert_eq!(
            poly.as_ref().len(),
            domain_size,
            "poly idx {} has mistmatchin size {} whereas domain size {}",
            poly_idx,
            poly.as_ref().len(),
            domain_size,
        );
        for (el_idx, &coeff) in poly.as_ref().iter().enumerate() {
            let idx = el_idx * num_polys + poly_idx;
            combined[idx] = coeff;
        }
    }

    Polynomial::from_coeffs_unpadded(combined)
}

pub(crate) fn combine_mixed_degree_monomials<F: PrimeField>(monomials: &[&Polynomial<F, Coefficients>], domain_size: usize) -> Result<Polynomial<F, Coefficients>, SynthesisError> {
    let num_polys = monomials.len();
    let max_degree = monomials.iter().map(|m| m.size()).max().unwrap();
    assert!(domain_size.is_power_of_two());
    // C(x) = f0(x^k) + x*f1(x^k) + .. + x^{k-1}*f{k-1}(x^k) where
    // deg(C) = k*n+k-1 and
    let full_degree = max_degree * num_polys + num_polys - 1;
    let mut combined = vec![F::zero(); full_degree];

    for (poly_idx, poly) in monomials.iter().enumerate() {
        assert!(
            poly.as_ref().len() <= max_degree,
            "poly idx {} has degree {} higher than allowed degree {}",
            poly_idx,
            poly.size(),
            max_degree,
        );
        for (el_idx, &coeff) in poly.as_ref().iter().enumerate() {
            let idx = el_idx * num_polys + poly_idx;
            combined[idx] = coeff;
        }
    }
    Polynomial::from_coeffs_unpadded(combined)
}

pub(crate) fn multiply_monomial_with_multiple_sparse_polys<F: PrimeField>(
    worker: &Worker,
    poly: &Polynomial<F, Coefficients>,
    pairs: &[(usize, F)],
) -> Result<Polynomial<F, Coefficients>, SynthesisError> {
    // P(x) * (X^n - c0)(X^m-c1)
    // first multipliy all sparse polynomials together
    let (first_degree, first_c) = pairs[0].clone();
    let mut previous = vec![F::zero(); first_degree + 1];
    previous[0] = first_c;
    previous[0].negate();
    previous[first_degree] = F::one();
    for (degree, c) in pairs.iter().cloned().skip(1) {
        let mut current = vec![F::zero(); degree + 1];
        current[0] = c;
        current[0].negate();
        current[degree] = F::one();

        previous = multiply_monomials(&previous, &current);
    }
    let total_degree: usize = pairs.iter().cloned().map(|(d, _)| d).sum();
    assert_eq!(previous.len(), total_degree + 1);

    // then compute partial multiplication and sum into the result
    let full_degree = poly.size() + total_degree;
    let mut result = Polynomial::from_coeffs_unpadded(vec![F::zero(); full_degree])?;
    for (pos, coeff) in previous.iter().cloned().rev().enumerate() {
        let current_degree = total_degree - pos;
        let mut tmp = poly.clone();
        tmp.scale(worker, coeff);
        let mut tmp2 = Polynomial::from_coeffs_unpadded(vec![F::zero(); current_degree + poly.size()])?;
        tmp2.as_mut()[current_degree..].copy_from_slice(tmp.as_ref());
        result.add_assign(worker, &tmp2);
    }

    Ok(result)
}

pub(crate) fn multiply_monomial_with_sparse_poly<F: PrimeField>(
    worker: &Worker,
    poly: &Polynomial<F, Coefficients>,
    degree_of_sparse: usize,
    c: F,
) -> Result<Polynomial<F, Coefficients>, SynthesisError> {
    // P(x) * (X^n - c)
    let mut coeffs = vec![F::zero(); degree_of_sparse + poly.size()];
    coeffs[degree_of_sparse..].copy_from_slice(poly.as_ref());
    let mut result = Polynomial::from_coeffs_unpadded(coeffs)?;

    let mut minus_constant = c;
    minus_constant.negate();
    let mut tmp = poly.clone();
    tmp.scale(&worker, minus_constant);
    result.add_assign(&worker, &tmp);

    Ok(result)
}

pub(crate) fn multiply_monomials<F: PrimeField>(poly1: &[F], poly2: &[F]) -> Vec<F> {
    assert!(!poly1.is_empty());
    assert!(!poly2.is_empty());

    let coeffs1 = poly1.as_ref();
    let coeffs2 = poly2.as_ref();
    let total_degree = (coeffs1.len() - 1) + (coeffs2.len() - 1) + 1;
    let mut result = vec![F::zero(); total_degree];

    for (i, p1) in coeffs1.iter().enumerate() {
        for (j, p2) in coeffs2.iter().enumerate() {
            let mut tmp = p1.clone();
            tmp.mul_assign(p2);
            result[i + j].add_assign(&tmp);
        }
    }
    result
}

pub(crate) fn compute_lagrange_basis_inverses<F: PrimeField>(num_polys: usize, h: F, y: F) -> Vec<F> {
    assert!(num_polys.is_power_of_two());
    let degree = num_polys as u64;
    // L_i(x) = (w_i/(N*h^{N-1})) * (X^N-h^N)/(X-w_i*h)
    // so compute inverses
    // 1/(N*h^{N-1}*(X-w_i*h));

    let domain = Domain::new_for_size(degree).unwrap();
    let omega = domain.generator;

    let degree_as_fe = F::from_str(&degree.to_string()).unwrap();
    let mut constant_part = h.pow(&[degree - 1]);
    constant_part.mul_assign(&degree_as_fe);

    let mut current_omega = F::one();
    let mut inverses = vec![];
    for _ in 0..degree {
        let mut denum = current_omega;
        denum.mul_assign(&h);
        denum.negate();
        denum.add_assign(&y);
        denum.mul_assign(&constant_part);
        inverses.push(denum);
        current_omega.mul_assign(&omega);
    }

    batch_inversion(&mut inverses);

    inverses
}

pub(crate) fn compute_lagrange_basis_inverses_for_union_set<F: PrimeField>(num_polys: usize, h: F, h_shifted: F, y: F, omega: F) -> (Vec<F>, Vec<F>) {
    let degree = num_polys as usize;
    let degree_as_fe = F::from_str(&degree.to_string()).unwrap();

    // precompute constant parts
    let mut h_pows = vec![F::one(); 2 * degree];
    let mut h_shifted_pows = vec![F::one(); 2 * degree];
    let mut y_pows = vec![F::one(); 2 * degree + 1];
    let mut acc1 = F::one();
    let mut acc2 = F::one();
    let mut acc3 = F::one();
    for ((a, b), c) in h_pows.iter_mut().zip(h_shifted_pows.iter_mut()).zip(y_pows.iter_mut()) {
        *a = acc1;
        acc1.mul_assign(&h);

        *b = acc2;
        acc2.mul_assign(&h_shifted);

        *c = acc3;
        acc3.mul_assign(&y);
    }
    y_pows[2 * degree] = acc3;
    assert_eq!(acc3, y.pow(&[2 * degree as u64]));

    assert_eq!(h_pows.len(), 2 * degree);
    assert_eq!(h_shifted_pows.len(), 2 * degree);
    assert_eq!(y_pows.len(), 2 * degree + 1);

    assert_eq!(omega.pow(&[degree as u64]), F::one());

    // L_i(x) = w_i/N*(h^{2N-1}-h_s^N*h^{N-1})*(X^{2N}-X^N*(h^N+h_s^N)+(h*h_s)^{2N})/(X-h*w_i)

    // constant part of first lagrange poly
    // N*(h1^{2N-1} - h1^{N-1}*h2^N)
    let mut denum_constant_part_first = h_pows[degree - 1];
    denum_constant_part_first.mul_assign(&h_shifted_pows[degree]);
    denum_constant_part_first.negate();
    denum_constant_part_first.add_assign(&h_pows[2 * degree - 1]);
    denum_constant_part_first.mul_assign(&degree_as_fe);

    // constant part of second lagrange poly
    // N*(h_s^{2N-1} - h_s^{N-1}*h^N)
    let mut denum_constant_part_second = h_shifted_pows[degree - 1];
    denum_constant_part_second.mul_assign(&h_pows[degree]);
    denum_constant_part_second.negate();
    denum_constant_part_second.add_assign(&h_shifted_pows[2 * degree - 1]);
    denum_constant_part_second.mul_assign(&degree_as_fe);

    let mut first_lagrange_basis_evals = Vec::with_capacity(degree);
    let mut second_lagrange_basis_evals = Vec::with_capacity(degree);

    let mut first = y;
    first.sub_assign(&h);
    first.mul_assign(&denum_constant_part_first);
    first_lagrange_basis_evals.push(first);

    let mut second = y;
    second.sub_assign(&h_shifted);
    second.mul_assign(&denum_constant_part_second);
    second_lagrange_basis_evals.push(second);

    let mut current_omega_h = omega;
    current_omega_h.mul_assign(&h);

    let mut current_omega_h_shifted = omega;
    current_omega_h_shifted.mul_assign(&h_shifted);
    for _ in 1..num_polys {
        let mut first = y;
        first.sub_assign(&current_omega_h);
        first.mul_assign(&denum_constant_part_first);
        first_lagrange_basis_evals.push(first);

        let mut second = y;
        second.sub_assign(&current_omega_h_shifted);
        second.mul_assign(&denum_constant_part_second);
        second_lagrange_basis_evals.push(second);

        current_omega_h.mul_assign(&omega);
        current_omega_h_shifted.mul_assign(&omega);
    }

    let mut evals = first_lagrange_basis_evals;
    evals.extend(second_lagrange_basis_evals);

    batch_inversion(&mut evals);

    let inverses = evals.clone();

    // numerator is the same for both lagrange polys
    // (X^{2N}-X^N*(h^N + h_s^N)+(h*h_s)^N)
    let mut tmp = h_pows[degree];
    tmp.mul_assign(&h_shifted_pows[degree]);

    let mut num_at_y = h_pows[degree];
    num_at_y.add_assign(&h_shifted_pows[degree]);
    num_at_y.mul_assign(&y_pows[degree]);
    num_at_y.negate();
    num_at_y.add_assign(&tmp);
    num_at_y.add_assign(&y_pows[2 * degree]);

    let mut current_omega = F::one();
    for idx in 0..degree {
        evals[idx].mul_assign(&num_at_y);
        evals[idx].mul_assign(&current_omega);

        evals[idx + degree].mul_assign(&num_at_y);
        evals[idx + degree].mul_assign(&current_omega);

        current_omega.mul_assign(&omega);
    }

    (inverses, evals)
}

pub(crate) fn precompute_lagrange_basis_evaluations_from_inverses<F: PrimeField>(num_polys: usize, inverses_of_denums: &[F], h: F, y: F, omega: F) -> Vec<F> {
    assert!(num_polys.is_power_of_two());
    assert_eq!(inverses_of_denums.len(), num_polys);

    // L_i(x) = (w_i/(N*h^{N-1})) * (X^N-h^N)/(X-w_i*h)
    // we have 1/(N*h^{N-1}*(X-w_i*h))
    let degree = num_polys as u64;
    let mut h_pow = h.pow(&[degree - 1]);
    let degree_as_fe = F::from_str(&degree.to_string()).unwrap();
    let mut constant_denum = degree_as_fe;
    constant_denum.mul_assign(&h_pow);
    h_pow.mul_assign(&h);
    let y_pow = y.pow(&[degree]);
    let mut num = h_pow;
    num.negate();
    num.add_assign(&y_pow);

    let mut current_omega = F::one();

    let mut result = vec![];
    for inv in inverses_of_denums.iter() {
        let mut denum = current_omega;
        denum.mul_assign(&h);
        denum.negate();
        denum.add_assign(&y);
        denum.mul_assign(&constant_denum);

        denum.mul_assign(inv);
        assert_eq!(denum, F::one());

        let mut value = num;
        value.mul_assign(&inv);
        value.mul_assign(&current_omega);
        result.push(value);

        current_omega.mul_assign(&omega);
    }
    assert_eq!(result.len(), inverses_of_denums.len());

    result
}

pub(crate) fn precompute_lagrange_basis_evaluations_from_inverses_for_union_set<F: PrimeField>(num_polys: usize, inverses_of_denums: &[F], h: F, h_shifted: F, y: F, omega: F) -> Vec<F> {
    let degree = num_polys as usize;
    assert_eq!(inverses_of_denums.len(), num_polys as usize * 2);
    let degree_as_fe = F::from_str(&degree.to_string()).unwrap();

    // precompute constant parts
    let mut h_pows = vec![F::one(); 2 * degree];
    let mut h_shifted_pows = vec![F::one(); 2 * degree];
    let mut y_pows = vec![F::one(); 2 * degree + 1];
    let mut acc1 = F::one();
    let mut acc2 = F::one();
    let mut acc3 = F::one();
    for ((a, b), c) in h_pows.iter_mut().zip(h_shifted_pows.iter_mut()).zip(y_pows.iter_mut()) {
        *a = acc1;
        acc1.mul_assign(&h);

        *b = acc2;
        acc2.mul_assign(&h_shifted);

        *c = acc3;
        acc3.mul_assign(&y);
    }
    y_pows[2 * degree] = acc3;
    assert_eq!(acc3, y.pow(&[2 * degree as u64]));

    assert_eq!(h_pows.len(), 2 * degree);
    assert_eq!(h_shifted_pows.len(), 2 * degree);
    assert_eq!(y_pows.len(), 2 * degree + 1);

    // L_i(x) = w_i/N*(h^{2N-1}-h_s^N*h^{N-1})*(X^{2N}-X^N*(h^N+h_s^N)+(h*h_s)^{2N})/(X-h*w_i)

    // constant part of first lagrange poly
    // N*(h1^{2N-1} - h1^{N-1}*h2^N)
    let mut denum_constant_part_first = h_pows[degree - 1];
    denum_constant_part_first.mul_assign(&h_shifted_pows[degree]);
    denum_constant_part_first.negate();
    denum_constant_part_first.add_assign(&h_pows[2 * degree - 1]);
    denum_constant_part_first.mul_assign(&degree_as_fe);

    // constant part of second lagrange poly
    // N*(h_s^{2N-1} - h_s^{N-1}*h^N)
    let mut denum_constant_part_second = h_shifted_pows[degree - 1];
    denum_constant_part_second.mul_assign(&h_pows[degree]);
    denum_constant_part_second.negate();
    denum_constant_part_second.add_assign(&h_shifted_pows[2 * degree - 1]);
    denum_constant_part_second.mul_assign(&degree_as_fe);

    let mut first_lagrange_basis_evals = Vec::with_capacity(degree as usize);
    let mut second_lagrange_basis_evals = Vec::with_capacity(degree as usize);

    let mut first = y;
    first.sub_assign(&h);
    first.mul_assign(&denum_constant_part_first);
    first_lagrange_basis_evals.push(first);

    let mut second = y;
    second.sub_assign(&h_shifted);
    second.mul_assign(&denum_constant_part_second);
    second_lagrange_basis_evals.push(second);

    let mut current_omega_h = omega;
    current_omega_h.mul_assign(&h);

    let mut current_omega_h_shifted = omega;
    current_omega_h_shifted.mul_assign(&h_shifted);
    for _ in 1..degree {
        let mut first = y;
        first.sub_assign(&current_omega_h);
        first.mul_assign(&denum_constant_part_first);
        first_lagrange_basis_evals.push(first);

        let mut second = y;
        second.sub_assign(&current_omega_h_shifted);
        second.mul_assign(&denum_constant_part_second);
        second_lagrange_basis_evals.push(second);

        current_omega_h.mul_assign(&omega);
        current_omega_h_shifted.mul_assign(&omega);
    }

    let mut evals = first_lagrange_basis_evals;
    evals.extend(second_lagrange_basis_evals);

    for (inverse, actual) in inverses_of_denums.iter().zip(evals.iter()) {
        let mut tmp = inverse.clone();
        tmp.mul_assign(actual);
        assert_eq!(tmp, F::one());
    }
    let mut evals = inverses_of_denums.to_vec();
    // numerator is the same for both lagrange polys
    // (X^{2N}-X^N*(h^N + h_s^N)+(h*h_s)^N)
    let mut tmp = h_pows[degree];
    tmp.mul_assign(&h_shifted_pows[degree]);

    let mut num_at_y = h_pows[degree];
    num_at_y.add_assign(&h_shifted_pows[degree]);
    num_at_y.mul_assign(&y_pows[degree]);
    num_at_y.negate();
    num_at_y.add_assign(&tmp);
    num_at_y.add_assign(&y_pows[2 * degree]);

    let mut current_omega = F::one();
    for idx in 0..degree {
        evals[idx].mul_assign(&num_at_y);
        evals[idx].mul_assign(&current_omega);

        evals[idx + degree].mul_assign(&num_at_y);
        evals[idx + degree].mul_assign(&current_omega);

        current_omega.mul_assign(&omega);
    }

    evals
}

pub(crate) fn precompute_all_lagrange_basis_evaluations<F: PrimeField>(
    interpolation_size_of_setup: usize,
    interpolation_size_of_first_round: usize,
    interpolation_size_of_second_round: usize,
    h0: (F, Option<F>),
    h1: (F, Option<F>),
    h2: (F, F),
    y: F,
    setup_requires_opening_at_shifted_point: bool,
    first_round_requires_opening_at_shifted_point: bool,
) -> (Vec<F>, [Vec<F>; 3]) {
    assert!(interpolation_size_of_setup.is_power_of_two());
    assert!(interpolation_size_of_first_round.is_power_of_two());
    assert_eq!(interpolation_size_of_second_round, 3);

    let [omega_setup, omega_first_round, omega_second_round] = compute_generators(interpolation_size_of_setup, interpolation_size_of_first_round, interpolation_size_of_second_round);

    let (h0, h0_shifted) = h0;
    let (h1, h1_shifted) = h1;

    let (lagrange_basis_inverses_of_setup_polys, lagrange_basis_evals_of_setup_polys) = if setup_requires_opening_at_shifted_point {
        let h0_shifted = h0_shifted.expect("h0 shifted");
        let (lagrange_basis_inverses_of_setup, _lagrange_basis_evals_of_setup) = compute_lagrange_basis_inverses_for_union_set(interpolation_size_of_setup, h0, h0_shifted, y, omega_setup);
        let lagrange_basis_evals_of_setup =
            precompute_lagrange_basis_evaluations_from_inverses_for_union_set(interpolation_size_of_setup, &lagrange_basis_inverses_of_setup, h0, h0_shifted, y, omega_setup);
        if SANITY_CHECK {
            assert_eq!(lagrange_basis_evals_of_setup, _lagrange_basis_evals_of_setup);
        }

        (lagrange_basis_inverses_of_setup, lagrange_basis_evals_of_setup)
    } else {
        let lagrange_basis_inverses_of_setup = compute_lagrange_basis_inverses(interpolation_size_of_setup, h0, y);

        let lagrange_basis_evals_of_setup = precompute_lagrange_basis_evaluations_from_inverses(interpolation_size_of_setup, &lagrange_basis_inverses_of_setup, h0, y, omega_setup);

        (lagrange_basis_inverses_of_setup, lagrange_basis_evals_of_setup)
    };

    let (lagrange_basis_inverses_of_first_round, lagrange_basis_evals_of_first_round) = if first_round_requires_opening_at_shifted_point {
        let h1_shifted = h1_shifted.expect("h1 shifted");
        let (lagrange_basis_inverses_of_first_round, _lagrange_basis_evals_of_first_round) =
            compute_lagrange_basis_inverses_for_union_set(interpolation_size_of_first_round, h1, h1_shifted, y, omega_first_round);
        let lagrange_basis_evals_of_first_round =
            precompute_lagrange_basis_evaluations_from_inverses_for_union_set(interpolation_size_of_first_round, &lagrange_basis_inverses_of_first_round, h1, h1_shifted, y, omega_first_round);
        if SANITY_CHECK {
            assert_eq!(lagrange_basis_evals_of_first_round, _lagrange_basis_evals_of_first_round);
        }

        (lagrange_basis_inverses_of_first_round, lagrange_basis_evals_of_first_round)
    } else {
        let lagrange_basis_inverses_of_first_round = compute_lagrange_basis_inverses(interpolation_size_of_first_round, h1, y);

        let lagrange_basis_evals_of_first_round =
            precompute_lagrange_basis_evaluations_from_inverses(interpolation_size_of_first_round, &lagrange_basis_inverses_of_first_round, h1, y, omega_first_round);

        (lagrange_basis_inverses_of_first_round, lagrange_basis_evals_of_first_round)
    };

    let (h2, h2_shifted) = h2;
    let (lagrange_basis_inverses_of_second_round, _lagrange_basis_evals_of_second_round) =
        compute_lagrange_basis_inverses_for_union_set(interpolation_size_of_second_round, h2, h2_shifted, y, omega_second_round);
    let lagrange_basis_evals_of_second_round =
        precompute_lagrange_basis_evaluations_from_inverses_for_union_set(interpolation_size_of_second_round, &lagrange_basis_inverses_of_second_round, h2, h2_shifted, y, omega_second_round);
    if SANITY_CHECK {
        assert_eq!(lagrange_basis_evals_of_second_round, _lagrange_basis_evals_of_second_round);
    }
    let mut flattened_inverses = lagrange_basis_inverses_of_setup_polys;
    flattened_inverses.extend(lagrange_basis_inverses_of_first_round);
    flattened_inverses.extend(lagrange_basis_inverses_of_second_round);
    (
        flattened_inverses,
        [lagrange_basis_evals_of_setup_polys, lagrange_basis_evals_of_first_round, lagrange_basis_evals_of_second_round],
    )
}

pub fn precompute_all_lagrange_basis_evaluations_from_inverses<F: PrimeField>(
    inverses_of_denums: &[F],
    interpolation_size_of_setup: usize,
    interpolation_size_of_first_round: usize,
    interpolation_size_of_second_round: usize,
    h0: (F, Option<F>),
    h1: (F, Option<F>),
    h2: (F, F),
    y: F,
    setup_requires_opening_at_shifted_point: bool,
    first_round_requires_opening_at_shifted_point: bool,
) -> [Vec<F>; 3] {
    assert!(interpolation_size_of_setup.is_power_of_two());
    assert!(interpolation_size_of_first_round.is_power_of_two());
    assert_eq!(interpolation_size_of_second_round, 3);

    let [omega_setup, omega_first_round, omega_second_round] = compute_generators(interpolation_size_of_setup, interpolation_size_of_first_round, interpolation_size_of_second_round);

    let h2 = (h2.0, Some(h2.1));

    let mut inverses_iter = inverses_of_denums;

    let mut result = [vec![], vec![], vec![]];

    for (idx, (interpolation_set_size, omega, h, requires_opening_at_shifted_point)) in [
        (interpolation_size_of_setup, omega_setup, h0, setup_requires_opening_at_shifted_point),
        (interpolation_size_of_first_round, omega_first_round, h1, first_round_requires_opening_at_shifted_point),
        (interpolation_size_of_second_round, omega_second_round, h2, true),
    ]
    .into_iter()
    .enumerate()
    {
        let (h, h_shifted) = h;
        let basis_evals = if requires_opening_at_shifted_point {
            let (current_inverses, rest) = inverses_iter.split_at(2 * interpolation_set_size);
            inverses_iter = rest;

            precompute_lagrange_basis_evaluations_from_inverses_for_union_set(interpolation_set_size, &current_inverses, h, h_shifted.expect("h shifted"), y, omega)
        } else {
            let (current_inverses, rest) = inverses_iter.split_at(interpolation_set_size);
            inverses_iter = rest;
            precompute_lagrange_basis_evaluations_from_inverses(interpolation_set_size, &current_inverses, h, y, omega)
        };
        result[idx] = basis_evals;
    }

    result
}

pub(crate) fn batch_inversion<F: PrimeField>(values: &mut [F]) {
    // a    1       a^-1        a^-1    -
    // b    a       (ab)^-1     b^-1    a^-1
    // c    ab      (abc)^-1    c^-1    (ab)^-1
    // d    abc     (abcd)^-1   d^-1    (abc)^-1

    let mut products = vec![];
    let mut acc = F::one();

    for el in values.iter() {
        products.push(acc);
        acc.mul_assign(&el);
    }
    assert_eq!(values.len(), products.len());

    let mut inv = acc.inverse().expect("inverse");

    for (el, src) in products.iter_mut().rev().zip(values.iter().rev()) {
        el.mul_assign(&inv);
        inv.mul_assign(src);
    }

    for (a, b) in values.iter().zip(products.iter()) {
        let mut tmp = a.clone();
        tmp.mul_assign(b);
        assert_eq!(tmp, F::one());
    }

    values.copy_from_slice(&products);
}

pub(crate) fn horner_evaluation<F: PrimeField>(coeffs: &[F], x: F) -> F {
    // c0 + c1*x + c2*x^2 + c3*x^3
    // c0 + x*(c1 + x*(c2 + x*c3))
    let mut sum = coeffs.last().unwrap().clone();
    for coeff in coeffs.iter().rev().skip(1) {
        sum.mul_assign(&x);
        sum.add_assign(&coeff);
    }

    sum
}

pub fn compute_generators<F: PrimeField>(interpolation_size_of_setup: usize, interpolation_size_of_first_round: usize, interpolation_size_of_second_round: usize) -> [F; 3] {
    let omega_setup = Domain::new_for_size(interpolation_size_of_setup as u64).unwrap().generator;
    let omega_first_round = Domain::new_for_size(interpolation_size_of_first_round as u64).unwrap().generator;
    let omega_second_round: F = hardcoded_generator_of_3rd_roots_of_unity();
    assert_eq!(omega_second_round.pow(&[interpolation_size_of_second_round as u64]), F::one());

    [omega_setup, omega_first_round, omega_second_round]
}

pub fn compute_opening_points<F: PrimeField>(
    r: F,
    z: F,
    z_omega: F,
    power: usize,
    interpolation_size_of_setup: usize,
    interpolation_size_of_first_round: usize,
    interpolation_size_of_second_round: usize,
    domain_size: usize,
    setup_requires_opening_at_shifted_point: bool,
    first_round_requires_opening_at_shifted_point: bool,
) -> ((F, Option<F>), (F, Option<F>), (F, F)) {
    let mut points = vec![];
    for (interpolation_set_size, requires_opening_at_shifted_point) in [
        (interpolation_size_of_second_round, true),
        (interpolation_size_of_first_round, first_round_requires_opening_at_shifted_point),
        (interpolation_size_of_setup, setup_requires_opening_at_shifted_point),
    ] {
        assert_eq!(power % interpolation_set_size, 0);
        let exp = (power / interpolation_set_size) as u64;
        let h = r.pow(&[exp]);
        assert_eq!(h.pow(&[interpolation_set_size as u64]), z);

        let h_shifted = if requires_opening_at_shifted_point {
            assert_eq!(interpolation_set_size, 3);
            // h_s = (w)^{1/3}*h
            let mut h_shifted: F = compute_cubic_root_of_domain(domain_size as u64);
            h_shifted.mul_assign(&h);

            assert_eq!(h_shifted.pow(&[interpolation_set_size as u64]), z_omega);

            Some(h_shifted)
        } else {
            None
        };

        points.push((h, h_shifted));
    }

    let h0 = points.pop().unwrap();
    let h1 = points.pop().unwrap();
    let (h2, h2_shifted) = points.pop().unwrap();

    (h0, h1, (h2, h2_shifted.expect("h2 shifted")))
}

#[derive(Default, Debug)]
pub struct RecomptuedQuotientEvaluations<F: PrimeField> {
    pub main_gate_quotient_at_z: F,
    pub custom_gate_quotient_at_z: Option<F>,
    pub copy_permutation_first_quotient_at_z: F,
    pub copy_permutation_second_quotient_at_z: F,
    pub lookup_first_quotient_at_z: Option<F>,
    pub lookup_second_quotient_at_z: Option<F>,
    pub lookup_third_quotient_at_z: Option<F>,
}

impl<F: PrimeField> RecomptuedQuotientEvaluations<F> {
    pub fn flatten(&self) -> Vec<F> {
        let mut flattened = vec![self.main_gate_quotient_at_z];
        self.custom_gate_quotient_at_z.map(|v| flattened.push(v));
        flattened.push(self.copy_permutation_first_quotient_at_z);
        flattened.push(self.copy_permutation_second_quotient_at_z);
        self.lookup_first_quotient_at_z.map(|v| flattened.push(v));
        self.lookup_second_quotient_at_z.map(|v| flattened.push(v));
        self.lookup_third_quotient_at_z.map(|v| flattened.push(v));
        debug_assert_eq!(flattened.len(), 3);

        flattened
    }
}

fn vanishing_inv_at_z<F: PrimeField>(z: F, domain_size: usize) -> F {
    let vanishing_at_z = evaluate_vanishing_for_size(&z, domain_size as u64);

    vanishing_at_z.inverse().unwrap()
}

pub fn recompute_quotients_from_evaluations<F: PrimeField>(
    all_evaluations: &[F],
    evaluation_offsets: &EvaluationOffsets,
    public_inputs: &[F],
    z: F,
    domain_size: usize,
    main_gate: &'static str,
    beta_for_copy_permutation: F,
    gamma_for_copy_permutation: F,
    non_residues: &[F],
    num_state_polys: usize,
    custom_gate_name: Option<&'static str>,
    eta_for_lookup: Option<F>,
    beta_for_lookup: Option<F>,
    gamma_for_lookup: Option<F>,
) -> RecomptuedQuotientEvaluations<F> {
    let mut recomputed_quotients = RecomptuedQuotientEvaluations::default();

    recomputed_quotients.main_gate_quotient_at_z = recompute_main_gate_quotient_from_evaluations(all_evaluations, evaluation_offsets, public_inputs, z, domain_size, main_gate);
    if let Some(custom_gate_name) = custom_gate_name {
        recomputed_quotients.custom_gate_quotient_at_z = Some(recompute_custom_gate_quotient_from_evaluations(&all_evaluations, evaluation_offsets, z, domain_size, custom_gate_name));
    }

    let (copy_perm_first_quotient, copy_perm_second_quotient) = recompute_copy_perm_quotients_from_evaluations(
        all_evaluations,
        evaluation_offsets,
        beta_for_copy_permutation,
        gamma_for_copy_permutation,
        non_residues,
        z,
        num_state_polys,
        domain_size,
    );
    recomputed_quotients.copy_permutation_first_quotient_at_z = copy_perm_first_quotient;
    recomputed_quotients.copy_permutation_second_quotient_at_z = copy_perm_second_quotient;

    if let (Some(beta_for_lookup), Some(gamma_for_lookup), Some(eta_for_lookup)) = (beta_for_lookup, gamma_for_lookup, eta_for_lookup) {
        let (lookup_first_quotient, lookup_second_quotient, lookup_third_quotient) =
            recompute_lookup_quotients_from_evaluations(all_evaluations, evaluation_offsets, eta_for_lookup, beta_for_lookup, gamma_for_lookup, z, domain_size, num_state_polys);
        recomputed_quotients.lookup_first_quotient_at_z = Some(lookup_first_quotient);
        recomputed_quotients.lookup_second_quotient_at_z = Some(lookup_second_quotient);
        recomputed_quotients.lookup_third_quotient_at_z = Some(lookup_third_quotient);
    }

    recomputed_quotients
}

fn recompute_main_gate_quotient_from_evaluations<F: PrimeField>(
    evaluations: &[F],
    evaluations_offsets: &EvaluationOffsets,
    public_inputs: &[F],
    z: F,
    domain_size: usize,
    main_gate: &'static str,
) -> F {
    let vanishing_at_z_inv = vanishing_inv_at_z(z, domain_size);

    let mut public_inputs_at_z = F::zero();
    let domain = Domain::<F>::new_for_size(domain_size as u64).unwrap();
    for (poly_idx, input) in public_inputs.iter().enumerate() {
        let mut lagrange_at_z = evaluate_lagrange_poly_at_point(poly_idx, &domain, z).unwrap();
        lagrange_at_z.mul_assign(input);
        public_inputs_at_z.add_assign(&lagrange_at_z);
    }

    let mut main_gate_rhs = compute_quotient_of_main_gate_at_z_flattened(main_gate, evaluations, public_inputs_at_z, evaluations_offsets);
    main_gate_rhs.mul_assign(&vanishing_at_z_inv);

    main_gate_rhs
}

fn recompute_custom_gate_quotient_from_evaluations<F: PrimeField>(evaluations: &[F], evaluation_offsets: &EvaluationOffsets, z: F, domain_size: usize, custom_gate_name: &'static str) -> F {
    let vanishing_at_z_inv = vanishing_inv_at_z(z, domain_size);

    let mut custom_gate_rhs = compute_quotient_of_custom_gate_at_z_flattened(custom_gate_name, evaluations, &evaluation_offsets);
    custom_gate_rhs.mul_assign(&vanishing_at_z_inv);

    custom_gate_rhs
}

fn recompute_copy_perm_quotients_from_evaluations<F: PrimeField>(
    evaluations: &[F],
    evaluations_offsets: &EvaluationOffsets,
    beta_for_copy_permutation: F,
    gamma_for_copy_permutation: F,
    non_residues: &[F],
    z: F,
    num_state_polys: usize,
    domain_size: usize,
) -> (F, F) {
    let l_0_at_z = evaluate_l0_at_point(domain_size as u64, z).unwrap();
    let vanishing_at_z_inv = vanishing_inv_at_z(z, domain_size);
    // we have only 2 main gate types where both has the same number of variables
    // z(X)(A + beta*X + gamma)(B + beta*k_1*X + gamma)(C + beta*K_2*X + gamma)(D + beta*K_3*X + gamma) -
    // - (A + beta*perm_a(X) + gamma)(B + beta*perm_b(X) + gamma)(C + beta*perm_c(X) + gamma)*(D + beta*perm_d(X) + gamma)*Z(X*Omega)== 0
    let mut copy_permutation_first_quotient_rhs_num_part = z;
    copy_permutation_first_quotient_rhs_num_part.mul_assign(&beta_for_copy_permutation);
    copy_permutation_first_quotient_rhs_num_part.add_assign(&gamma_for_copy_permutation);
    copy_permutation_first_quotient_rhs_num_part.add_assign(&evaluations[evaluations_offsets.trace.trace_evaluations_at_z]);

    assert_eq!(non_residues.len() + 1, num_state_polys);
    for (non_residue, state_poly) in non_residues.iter().zip(
        evaluations[evaluations_offsets.trace.trace_evaluations_at_z..evaluations_offsets.trace.trace_evaluations_at_z + num_state_polys]
            .iter()
            .skip(1),
    ) {
        let mut tmp = z;
        tmp.mul_assign(&non_residue);
        tmp.mul_assign(&beta_for_copy_permutation);
        tmp.add_assign(&gamma_for_copy_permutation);
        tmp.add_assign(state_poly);
        copy_permutation_first_quotient_rhs_num_part.mul_assign(&tmp);
    }
    copy_permutation_first_quotient_rhs_num_part.mul_assign(&evaluations[evaluations_offsets.copy_permutation.grand_product_at_z]);

    let mut copy_permutation_first_quotient_rhs_denum_part = evaluations[evaluations_offsets.copy_permutation.grand_product_at_z_omega];
    for (permutation, state_poly) in evaluations[evaluations_offsets.setup.permutations_at_z..evaluations_offsets.setup.permutations_at_z + num_state_polys]
        .iter()
        .zip(evaluations[evaluations_offsets.trace.trace_evaluations_at_z..evaluations_offsets.trace.trace_evaluations_at_z + num_state_polys].iter())
    {
        let mut tmp = beta_for_copy_permutation;
        tmp.mul_assign(&permutation);
        tmp.add_assign(&gamma_for_copy_permutation);
        tmp.add_assign(state_poly);
        copy_permutation_first_quotient_rhs_denum_part.mul_assign(&tmp);
    }

    let mut copy_permutation_first_quotient_rhs = copy_permutation_first_quotient_rhs_num_part;
    copy_permutation_first_quotient_rhs.sub_assign(&copy_permutation_first_quotient_rhs_denum_part);
    copy_permutation_first_quotient_rhs.mul_assign(&vanishing_at_z_inv);
    // recomputed_quotients.copy_permutation_firt_quotient_at_z = copy_permutation_first_quotient_rhs;

    // (Z(x) - 1) * L_{0} == 0
    let mut copy_permutation_second_quotient_rhs = evaluations[evaluations_offsets.copy_permutation.grand_product_at_z];
    copy_permutation_second_quotient_rhs.sub_assign(&F::one());
    copy_permutation_second_quotient_rhs.mul_assign(&l_0_at_z);
    copy_permutation_second_quotient_rhs.mul_assign(&vanishing_at_z_inv);
    // recomputed_quotients.copy_permutation_second_quotient_at_z = copy_permutation_second_quotient_rhs;

    (copy_permutation_first_quotient_rhs, copy_permutation_second_quotient_rhs)
}

fn recompute_lookup_quotients_from_evaluations<F: PrimeField>(
    evaluations: &[F],
    evaluation_offsets: &EvaluationOffsets,
    eta_for_lookup: F,
    beta_for_lookup: F,
    gamma_for_lookup: F,
    z: F,
    domain_size: usize,
    num_state_polys: usize,
) -> (F, F, F) {
    let vanishing_at_z_inv = vanishing_inv_at_z(z, domain_size);
    let mut beta_gamma = beta_for_lookup;
    beta_gamma.add_assign(&F::one());
    beta_gamma.mul_assign(&gamma_for_lookup);
    // lookup identities
    // ( Z(x*omega)*(\gamma*(1 + \beta) + s(x) + \beta * s(x*omega))) -
    // - Z(x) * (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega)) )*(X - omega^{n-1})
    // LHS
    // f(z) = f0(z) + z*f1(z)
    let lookup_offsets = evaluation_offsets.lookup.as_ref().expect("lookup offsets");

    // RHS
    let mut lookup_first_quotient_rhs_denum_part = evaluations[lookup_offsets.s_poly_at_z_omega];
    lookup_first_quotient_rhs_denum_part.mul_assign(&beta_for_lookup);
    lookup_first_quotient_rhs_denum_part.add_assign(&evaluations[lookup_offsets.s_poly_at_z]);
    lookup_first_quotient_rhs_denum_part.add_assign(&beta_gamma);
    lookup_first_quotient_rhs_denum_part.mul_assign(&evaluations[lookup_offsets.grand_product_at_z_omega]);
    // Prover doesn't open aggregated columns of table rather it opens each of them
    // seperately because they are committed in the first round
    // and there is no reandomness.

    // aggregate witnesses a + eta*b + eta^2*c + eta^3*table_type
    // expands into (((table_type*eta + c)*eta  + b)*eta + a)
    let mut aggregated_lookup_f_at_z = evaluations[evaluation_offsets.setup.lookup_table_type_at_z];
    for col in evaluations[evaluation_offsets.trace.trace_evaluations_at_z..evaluation_offsets.trace.trace_evaluations_at_z + num_state_polys]
        .iter()
        .take(num_state_polys - 1)
        .rev()
    {
        aggregated_lookup_f_at_z.mul_assign(&eta_for_lookup);
        aggregated_lookup_f_at_z.add_assign(col);
    }
    aggregated_lookup_f_at_z.mul_assign(&evaluations[evaluation_offsets.setup.lookup_selector_at_z]);
    // col0 + eta * col1 + eta^2*col2 + eta^3*table_type
    let mut aggregated_lookup_table_cols_at_z = evaluations[evaluation_offsets.setup.lookup_tables_at_z + 3];
    let mut aggregated_lookup_table_cols_at_z_omega = evaluations[evaluation_offsets.setup.lookup_tables_at_z_omega + 3];
    for (at_z, at_z_omega) in evaluations[evaluation_offsets.setup.lookup_tables_at_z..evaluation_offsets.setup.lookup_tables_at_z + 3]
        .iter()
        .take(num_state_polys - 1)
        .rev()
        .zip(
            evaluations[evaluation_offsets.setup.lookup_tables_at_z_omega..evaluation_offsets.setup.lookup_tables_at_z_omega + 3]
                .iter()
                .rev(),
        )
    {
        aggregated_lookup_table_cols_at_z.mul_assign(&eta_for_lookup);
        aggregated_lookup_table_cols_at_z.add_assign(at_z);

        aggregated_lookup_table_cols_at_z_omega.mul_assign(&eta_for_lookup);
        aggregated_lookup_table_cols_at_z_omega.add_assign(at_z_omega);
    }
    aggregated_lookup_f_at_z.add_assign(&gamma_for_lookup);
    // We also need to aggregate shifted table columns to construct t(z*w)
    // First identity is for multiset-equality
    let mut lookup_first_quotient_rhs_num_part = aggregated_lookup_table_cols_at_z_omega;
    lookup_first_quotient_rhs_num_part.mul_assign(&beta_for_lookup);
    lookup_first_quotient_rhs_num_part.add_assign(&aggregated_lookup_table_cols_at_z);
    lookup_first_quotient_rhs_num_part.add_assign(&beta_gamma);

    lookup_first_quotient_rhs_num_part.mul_assign(&aggregated_lookup_f_at_z);
    let mut beta_one = beta_for_lookup;
    beta_one.add_assign(&F::one());
    lookup_first_quotient_rhs_num_part.mul_assign(&beta_one);
    lookup_first_quotient_rhs_num_part.mul_assign(&evaluations[lookup_offsets.grand_product_at_z]);

    let mut lookup_first_quotient_rhs = lookup_first_quotient_rhs_denum_part;
    lookup_first_quotient_rhs.sub_assign(&lookup_first_quotient_rhs_num_part);

    let domain = Domain::<F>::new_for_size(domain_size as u64).unwrap();
    let last_omega = domain.generator.pow(&[domain_size as u64 - 1]);
    let mut tmp = z;
    tmp.sub_assign(&last_omega);
    lookup_first_quotient_rhs.mul_assign(&tmp);
    lookup_first_quotient_rhs.mul_assign(&vanishing_at_z_inv);

    let l_0_at_z = evaluate_l0_at_point(domain_size as u64, z).unwrap();
    // Then verify that first element of the grand product poly equals to 1
    // (Z(x) - 1) * L_{0} == 0
    let mut lookup_second_quotient_rhs = evaluations[lookup_offsets.grand_product_at_z];
    lookup_second_quotient_rhs.sub_assign(&F::one());
    lookup_second_quotient_rhs.mul_assign(&l_0_at_z);
    lookup_second_quotient_rhs.mul_assign(&vanishing_at_z_inv);

    // Also verify that last element is equals to expected value
    // (Z(x) - expected) * L_{n-1} == 0
    let expected = beta_gamma.pow([(domain_size - 1) as u64]);
    let l_last = evaluate_lagrange_poly_at_point(domain_size - 1, &domain, z).unwrap();
    let mut lookup_third_quotient_rhs = evaluations[lookup_offsets.grand_product_at_z];
    lookup_third_quotient_rhs.sub_assign(&expected);
    lookup_third_quotient_rhs.mul_assign(&l_last);
    lookup_third_quotient_rhs.mul_assign(&vanishing_at_z_inv);

    (lookup_first_quotient_rhs, lookup_second_quotient_rhs, lookup_third_quotient_rhs)
}

pub fn evaluate_r_polys_at_point_with_flattened_evals_and_precomputed_basis<F: PrimeField>(
    all_evaluations: Vec<F>,
    recomputed_quotient_evaluations: &RecomptuedQuotientEvaluations<F>,
    num_setup_polys: usize,
    num_first_round_polys: usize,
    num_second_round_polys: usize,
    h0: (F, Option<F>),
    h1: (F, Option<F>),
    h2: (F, F),
    precomputed_basis_evals: [Vec<F>; 3],
    setup_requires_opening_at_shifted_point: bool,
    first_round_requires_opening_at_shifted_point: bool,
) -> [F; 3] {
    let [setup_omega, first_round_omega, second_round_omega] = compute_generators(num_setup_polys, num_first_round_polys, num_second_round_polys);
    // r_i(X) polys needs to be evaluated at random point.
    // We are going to construct C_i(w_i*h) from existing
    // evaluations of the system polynomials. Then calculate r_i(y)
    // by evaluating each basis poly at the given point

    let mut r_evals = [F::zero(); 3];

    let mut all_evaluations_iter = all_evaluations.into_iter();
    let setup_evals: Vec<_> = all_evaluations_iter.by_ref().take(num_setup_polys).collect();
    let mut first_round_evals: Vec<_> = all_evaluations_iter.by_ref().take(num_first_round_polys - 1).collect();
    first_round_evals.push(recomputed_quotient_evaluations.main_gate_quotient_at_z);
    let mut second_round_evals: Vec<_> = all_evaluations_iter.by_ref().take(num_second_round_polys - 2).collect();
    second_round_evals.push(recomputed_quotient_evaluations.copy_permutation_first_quotient_at_z);
    second_round_evals.push(recomputed_quotient_evaluations.copy_permutation_second_quotient_at_z);

    let setup_evals_shifted = if setup_requires_opening_at_shifted_point {
        all_evaluations_iter.by_ref().take(num_setup_polys).collect()
    } else {
        vec![]
    };

    let first_round_evals_shifted = if first_round_requires_opening_at_shifted_point {
        all_evaluations_iter.by_ref().take(num_first_round_polys).collect()
    } else {
        vec![]
    };

    let second_round_evals_shifted: Vec<_> = all_evaluations_iter.by_ref().take(num_second_round_polys).collect();
    assert_eq!(second_round_evals.len(), second_round_evals_shifted.len());

    assert!(all_evaluations_iter.next().is_none());

    for ((r, (interpolation_set_size, h, omega, requires_opening_at_shifted_point, evals, evals_shifted)), precomputed_basis_evals) in r_evals
        .iter_mut()
        .zip([
            (num_setup_polys, h0, setup_omega, setup_requires_opening_at_shifted_point, setup_evals, setup_evals_shifted),
            (
                num_first_round_polys,
                h1,
                first_round_omega,
                first_round_requires_opening_at_shifted_point,
                first_round_evals,
                first_round_evals_shifted,
            ),
            (num_second_round_polys, (h2.0, Some(h2.1)), second_round_omega, true, second_round_evals, second_round_evals_shifted),
        ])
        .zip(precomputed_basis_evals.into_iter())
    {
        let (h, h_shifted) = h;
        let mut current_omega = h;
        if requires_opening_at_shifted_point {
            assert_eq!(precomputed_basis_evals.len(), 2 * interpolation_set_size);
            let (basis_evals, basis_evals_shifted) = precomputed_basis_evals.split_at(interpolation_set_size);
            let mut current_omega_shifted = h_shifted.expect("h shifted");
            for (basis, basis_shifted) in basis_evals.iter().zip(basis_evals_shifted.iter()) {
                let mut sum = horner_evaluation(&evals, current_omega);
                sum.mul_assign(basis);
                r.add_assign(&sum);

                let mut sum_shifted = horner_evaluation(&evals_shifted, current_omega_shifted);
                sum_shifted.mul_assign(basis_shifted);
                r.add_assign(&sum_shifted);

                current_omega.mul_assign(&omega);
                current_omega_shifted.mul_assign(&omega);
            }
        } else {
            assert_eq!(precomputed_basis_evals.len(), interpolation_set_size);
            for basis in precomputed_basis_evals.iter() {
                let mut sum = horner_evaluation(&evals, current_omega);
                sum.mul_assign(basis);
                r.add_assign(&sum);

                current_omega.mul_assign(&omega);
            }
        }
    }

    r_evals
}

pub(crate) fn interpolate_r_monomial<F: PrimeField>(evaluations: &[F], degree: usize, point: F, omega: F) -> Vec<F> {
    assert!(degree.is_power_of_two());
    assert_eq!(evaluations.len(), degree);

    // L_i(x) = (w_i/h^{N-1})*(X^N-h^N)/(X-w_i*h)
    let point_pow = point.pow(&[degree as u64 - 1]);
    let mut num = vec![F::zero(); degree + 1];
    num[degree] = F::one();
    num[0] = point_pow;
    num[0].mul_assign(&point);
    num[0].negate();
    let inv_degree_as_fe = F::from_str(&degree.to_string()).unwrap();
    let mut constant_denum = inv_degree_as_fe;
    constant_denum.mul_assign(&point_pow);
    let constant_denum = constant_denum.inverse().unwrap();

    let mut current_omega = F::one();
    let mut result = vec![F::zero(); degree];
    for value in evaluations.iter() {
        let mut constant_part = current_omega;
        constant_part.mul_assign(&constant_denum);
        constant_part.mul_assign(value);

        let mut point = point;
        point.mul_assign(&current_omega);
        let mut quotient = divide_by_linear_term(&num, point);
        assert_eq!(quotient.len(), result.len());

        for (coeff, r) in quotient.iter_mut().zip(result.iter_mut()) {
            coeff.mul_assign(&constant_part);
            r.add_assign(coeff);
        }
        if SANITY_CHECK {
            assert_eq!(horner_evaluation(&quotient, point), *value);
        }
        current_omega.mul_assign(&omega);
    }

    result
}

pub(crate) fn interpolate_r_monomial_from_union_set<F: PrimeField>(values_over_set: &[F], values_over_shifted_set: &[F], h: F, h_shifted: F, omega: F, num_polys: usize) -> Vec<F> {
    assert_eq!(values_over_set.len(), num_polys);
    assert_eq!(values_over_shifted_set.len(), num_polys);
    let degree = num_polys;
    let degree_as_fe = F::from_str(&degree.to_string()).unwrap();
    // we have to show opening for two sets
    // L_i(x) = w_i/N*(h^{2N-1}-h_s^N*h^{N-1})*(X^{2N}-X^N*(h^N+h_s^N)+(h*h_s)^{N})/(X-h*w_i)
    // and
    // L_j(x) = w_j/N*(h_s^{2N-1}-h^N*h_s^{N-1})*(X^{2N}-X^N*(h^N+h_s^N)+(h*h_s)^{N})/(X-h_s*w_j))

    let mut h_pows = vec![F::one(); 2 * degree];
    let mut h_shifted_pows = vec![F::one(); 2 * degree];
    let mut acc1 = F::one();
    let mut acc2 = F::one();
    for (a, b) in h_pows.iter_mut().zip(h_shifted_pows.iter_mut()) {
        *a = acc1;
        acc1.mul_assign(&h);

        *b = acc2;
        acc2.mul_assign(&h_shifted);
    }

    // constant part of first lagrange poly
    // N*(h1^{2N-1} - h1^{N-1}*h2^N)
    let mut denum_constant_part_first = h_pows[degree - 1];
    denum_constant_part_first.mul_assign(&h_shifted_pows[degree]);
    denum_constant_part_first.negate();
    denum_constant_part_first.add_assign(&h_pows[2 * degree - 1]);
    denum_constant_part_first.mul_assign(&degree_as_fe);
    let denum_constant_part_first = denum_constant_part_first.inverse().unwrap();

    // constant part of second lagrange poly
    // N*(h_s^{2N-1} - h_s^{N-1}*h^N)
    let mut denum_constant_part_second = h_shifted_pows[degree - 1];
    denum_constant_part_second.mul_assign(&h_pows[degree]);
    denum_constant_part_second.negate();
    denum_constant_part_second.add_assign(&h_shifted_pows[2 * degree - 1]);
    denum_constant_part_second.mul_assign(&degree_as_fe);
    let denum_constant_part_second = denum_constant_part_second.inverse().unwrap();

    // numerator is the same for both lagrange polys
    // (X^{2N}-X^N*(h^N + h_s^N)+(h*h_s)^N)
    let mut num = vec![F::zero(); 2 * degree + 1];
    num[0] = h_pows[degree];
    num[0].mul_assign(&h_shifted_pows[degree]);
    num[degree] = h_pows[degree];
    num[degree].add_assign(&h_shifted_pows[degree]);
    num[degree].negate();
    num[2 * degree] = F::one();

    let mut powers_of_omega = vec![];
    let mut acc = F::one();
    for _ in 0..degree {
        powers_of_omega.push(acc);
        acc.mul_assign(&omega);
    }

    let mut all_lagrange_monomials = vec![];
    for (omega_i, (value, value_shifted)) in powers_of_omega.iter().cloned().zip(values_over_set.iter().zip(values_over_shifted_set.iter())) {
        // (X-h*w_i)
        let mut divide_by = omega_i;
        divide_by.mul_assign(&h);
        let mut q1 = divide_by_linear_term(&num, divide_by);
        if SANITY_CHECK {
            assert_eq!(horner_evaluation(&num, divide_by), F::zero());

            let q = Polynomial::from_coeffs_unpadded(q1.clone()).unwrap();
            let actual = multiply_monomial_with_sparse_poly(&Worker::new(), &q, 1, divide_by).unwrap();
            assert_eq!(&num, &actual.as_ref()[..num.len()]);
        }

        let mut divide_by = omega_i;
        divide_by.mul_assign(&h_shifted);
        let mut q2 = divide_by_linear_term(&num, divide_by);
        if SANITY_CHECK {
            assert_eq!(horner_evaluation(&num, divide_by), F::zero());

            let q = Polynomial::from_coeffs_unpadded(q2.clone()).unwrap();
            let actual = multiply_monomial_with_sparse_poly(&Worker::new(), &q, 1, divide_by).unwrap();
            assert_eq!(&num, &actual.as_ref()[..num.len()]);
        }

        let mut tmp = denum_constant_part_first;
        tmp.mul_assign(&omega_i);
        tmp.mul_assign(&value);
        for el in q1.iter_mut() {
            el.mul_assign(&tmp);
        }
        let mut tmp = denum_constant_part_second;
        tmp.mul_assign(&omega_i);
        tmp.mul_assign(&value_shifted);
        for el in q2.iter_mut() {
            el.mul_assign(&tmp);
        }
        if SANITY_CHECK {
            // L_i(h*w_i) = 1 and 0 otherwise
            let mut point = h;
            point.mul_assign(&omega_i);
            assert_eq!(horner_evaluation(&q1, point), *value);

            let mut point = h_shifted;
            point.mul_assign(&omega_i);
            assert_eq!(horner_evaluation(&q2, point), *value_shifted);

            let mut next_point = h;
            next_point.mul_assign(&omega_i);
            next_point.mul_assign(&omega);
            assert_eq!(horner_evaluation(&q1, next_point), F::zero());

            let mut next_point = h_shifted;
            next_point.mul_assign(&omega_i);
            next_point.mul_assign(&omega);
            assert_eq!(horner_evaluation(&q2, next_point), F::zero());
        }
        all_lagrange_monomials.push(q1);
        all_lagrange_monomials.push(q2);
    }

    let degree = all_lagrange_monomials[0].len();
    all_lagrange_monomials.iter().for_each(|m| assert_eq!(m.len(), degree));

    let mut result = all_lagrange_monomials.pop().unwrap();
    for m in all_lagrange_monomials.iter() {
        for (a, b) in result.iter_mut().zip(m.iter()) {
            a.add_assign(b)
        }
    }
    if SANITY_CHECK {
        for (omega_i, (value, value_shifted)) in powers_of_omega.into_iter().zip(values_over_set.iter().zip(values_over_shifted_set.iter())) {
            let mut point = h;
            point.mul_assign(&omega_i);
            assert_eq!(horner_evaluation(&result, point), *value);

            let mut point = h_shifted;
            point.mul_assign(&omega_i);
            assert_eq!(horner_evaluation(&result, point), *value_shifted);
        }
    }

    result
}

pub(crate) fn divide_by_linear_term<F: PrimeField>(poly: &[F], opening_point: F) -> Vec<F> {
    // we are only interested in quotient without a reminder, so we actually don't need opening value
    let mut b = opening_point;
    b.negate();

    let mut q = vec![F::zero(); poly.len()];

    let mut tmp = F::zero();
    let mut found_one = false;
    for (q, r) in q.iter_mut().rev().skip(1).zip(poly.iter().rev()) {
        if !found_one {
            if r.is_zero() {
                continue;
            } else {
                found_one = true;
            }
        }

        let mut lead_coeff = *r;
        lead_coeff.sub_assign(&tmp);
        *q = lead_coeff;
        tmp = lead_coeff;
        tmp.mul_assign(&b);
    }
    assert!(q.last().unwrap().is_zero());
    let _ = q.pop().unwrap();
    q
}

pub(crate) fn divide_by_higher_degree_term<F: PrimeField>(poly: &[F], divisor_degree: usize, point: F) -> Vec<F> {
    let dividend_degree = poly.len() - 1;
    let quotient_degree = dividend_degree - divisor_degree;
    let mut tmp = vec![F::zero(); poly.len()]; // TODO
    let mut quotient = vec![F::zero(); quotient_degree + 1];

    let mut minus_point = point;
    minus_point.negate();
    for (pos, coeff) in poly.iter().rev().enumerate() {
        let cur_degree = dividend_degree - pos;
        let mut lead_coeff = coeff.clone();
        lead_coeff.sub_assign(&tmp[cur_degree]);
        tmp[cur_degree] = F::zero();

        if lead_coeff.is_zero() {
            continue;
        }
        let cur_quotient_degree = cur_degree - divisor_degree;
        quotient[cur_quotient_degree] = lead_coeff;
        lead_coeff.mul_assign(&minus_point);
        tmp[cur_quotient_degree] = lead_coeff;
        if cur_degree == divisor_degree {
            break;
        }
    }

    quotient
}

pub(crate) fn divide_by_multiple_higher_degree_sparse_polys<F: PrimeField>(poly: &[F], mut sparse_polys: Vec<(usize, F)>) -> Vec<F> {
    let total_degree: usize = sparse_polys.iter().map(|(degree, _)| *degree).sum();
    let (degree, at) = sparse_polys.pop().unwrap();
    let mut result = divide_by_higher_degree_term(poly, degree, at);
    for (degree, at) in sparse_polys.iter().cloned() {
        result = divide_by_higher_degree_term(&result, degree, at);
    }
    if SANITY_CHECK {
        assert_eq!(result.len() + total_degree, poly.len());
    }
    result
}

pub(crate) fn ensure_single_poly_in_map_or_create<'a, 'b, E: Engine>(
    worker: &Worker,
    domain_size: usize,
    omegas_bitreversed: &BitReversedOmegas<E::Fr>,
    lde_factor: usize,
    coset_factor: E::Fr,
    monomials_map: &AssembledPolynomialStorageForMonomialForms<'a, E>,
    ldes_map: &mut AssembledPolynomialStorage<'b, E>,
    idx: PolyIdentifier,
) -> Result<(), SynthesisError> {
    let idx = PolynomialInConstraint::from_id(idx);
    ensure_in_map_or_create(worker, idx, domain_size, omegas_bitreversed, lde_factor, coset_factor, monomials_map, ldes_map)
}

pub fn compute_power_of_two_root_of_generator<F: PrimeField>(domain_size: u64, k_th_root: usize) -> F {
    assert!(k_th_root.is_power_of_two());
    let k_th_root = k_th_root as usize;

    let mut power_of_two = 0;
    let mut k = domain_size;
    while k != 1 {
        k >>= 1;
        power_of_two += 1;
    }
    let max_power_of_two = F::S as u64;
    if power_of_two > max_power_of_two {
        panic!("exceeded allowed power")
    }
    // (w^(2^n))^(2^(S-n)) = 1
    // g=(w^(2^n))^1/k)=(w^(2^n*2^-2^log(k)))
    // g=w^(2^(n-log(k)))
    let mut generator = F::root_of_unity();
    for _ in power_of_two..(max_power_of_two - k_th_root.trailing_zeros() as u64) {
        generator.square()
    }
    if SANITY_CHECK {
        // wj^32 =wj^32 = wi where wi is the generator of the main domain
        // and wj is the generator of the evaluation set
        let gen_of_main_domain = Domain::new_for_size(domain_size).unwrap().generator;
        assert_eq!(generator.pow(&[k_th_root as u64]), gen_of_main_domain);
    }

    generator
}

pub(crate) fn compute_quotient_monomial<F: PrimeField>(
    worker: &Worker,
    mut lde: Polynomial<F, Values>,
    divisor: &Polynomial<F, Values>,
    coset_factor: F,
    lde_factor: usize,
    domain_size: usize,
    quotient_degree: usize,
) -> Result<Polynomial<F, Coefficients>, SynthesisError> {
    assert_eq!(lde.size(), lde_factor * domain_size);
    assert_eq!(divisor.size(), lde_factor * domain_size);
    assert!(quotient_degree <= lde_factor);
    lde.bitreverse_enumeration(worker);
    lde.mul_assign(&worker, &divisor);
    let monomial = lde.icoset_fft_for_generator(worker, &coset_factor);
    // Our trace has n+1 points that interpolates to a degree n monomial
    // divison cancels n roots so that new_degree + n = original_degree
    let new_degree = quotient_degree * domain_size;
    // assert!(new_degree + domain_size <= monomial.size());
    if SANITY_CHECK {
        for el in monomial.as_ref()[new_degree..].iter() {
            assert!(el.is_zero())
        }
    }
    let mut coeffs = monomial.into_coeffs();
    coeffs.truncate(new_degree);
    Polynomial::from_coeffs_unpadded(coeffs)
}

pub(crate) fn evaluate_vanishing_for_size<F: PrimeField>(point: &F, vanishing_domain_size: u64) -> F {
    let mut result = point.pow(&[vanishing_domain_size]);
    result.sub_assign(&F::one());

    result
}

pub(crate) fn evaluate_l0_at_point<F: PrimeField>(domain_size: u64, point: F) -> Result<F, SynthesisError> {
    let size_as_fe = F::from_str(&format!("{}", domain_size)).unwrap();

    let mut den = point;
    den.sub_assign(&F::one());
    den.mul_assign(&size_as_fe);

    let den = den.inverse().ok_or(SynthesisError::DivisionByZero)?;

    let mut num = point.pow(&[domain_size]);
    num.sub_assign(&F::one());
    num.mul_assign(&den);

    Ok(num)
}

pub(crate) fn evaluate_lagrange_poly_at_point<F: PrimeField>(poly_number: usize, domain: &Domain<F>, point: F) -> Result<F, SynthesisError> {
    // lagrange polynomials have a form
    // (omega^i / N) / (X - omega^i) * (X^N - 1)

    let mut num = evaluate_vanishing_for_size(&point, domain.size);
    let omega_power = domain.generator.pow(&[poly_number as u64]);
    num.mul_assign(&omega_power);

    let size_as_fe = F::from_str(&format!("{}", domain.size)).unwrap();

    let mut den = point;
    den.sub_assign(&omega_power);
    den.mul_assign(&size_as_fe);

    let den = den.inverse().ok_or(SynthesisError::DivisionByZero)?;

    num.mul_assign(&den);

    Ok(num)
}

pub(crate) fn evaluate_multiple_sparse_polys<F: PrimeField>(pair: Vec<(usize, F)>, at_point: F) -> F {
    let mut acc = F::one();
    for (degree, c) in pair.into_iter() {
        // f(x) = X^n-c
        let mut tmp = at_point.pow(&[degree as u64]);
        tmp.sub_assign(&c);
        acc.mul_assign(&tmp);
    }

    acc
}

pub(crate) fn num_system_polys_from_assembly<E: Engine, MG: MainGate<E>, P: PlonkConstraintSystemParams<E>, S: SynthesisMode, C: Circuit<E>>(
    assembly: &Assembly<E, P, MG, S>,
) -> (usize, usize, usize, usize) {
    let gates = assembly.sorted_gates.clone();

    let has_custom_gate = gates.len() > 1;
    let has_lookup = assembly.num_table_lookups > 0;

    let num_state_polys = P::STATE_WIDTH;
    assert_eq!(gates[0].variable_polynomials().len(), num_state_polys);
    let num_witness_polys = P::WITNESS_WIDTH;

    let num_setup_polys = {
        let num_gate_setups = gates[0].setup_polynomials().len();
        let num_permutation_polys = num_state_polys;
        let num_gate_selectors = if has_custom_gate { 2 } else { 0 };
        let num_lookup_table_setup_polys = if has_lookup {
            // selector + table type + num cols
            1 + 1 + 4
        } else {
            0
        };

        num_gate_setups + num_permutation_polys + num_gate_selectors + num_lookup_table_setup_polys
    };

    num_system_polys(has_lookup, has_custom_gate, num_setup_polys, num_state_polys, num_witness_polys)
}

pub fn num_system_polys_from_vk<E: Engine, C: Circuit<E>>(vk: &FflonkVerificationKey<E, C>) -> (usize, usize, usize, usize) {
    let gates = sorted_gates_from_circuit_definitions::<_, C>();

    let has_custom_gate = gates.len() > 1;
    let has_lookup = vk.total_lookup_entries_length > 0;

    let num_state_polys = vk.num_state_polys;
    let num_witness_polys = vk.num_witness_polys;
    let mut num_setup_polys = gates[0].setup_polynomials().len();
    num_setup_polys += num_state_polys; // permutations
    if has_custom_gate {
        num_setup_polys += 2;
    }
    if has_lookup {
        num_setup_polys += 6;
    }
    num_system_polys(has_lookup, has_custom_gate, num_setup_polys, num_state_polys, num_witness_polys)
}

pub(crate) fn num_system_polys(has_lookup: bool, has_custom_gate: bool, num_setup_polys: usize, num_state_polys: usize, num_witness_polys: usize) -> (usize, usize, usize, usize) {
    let mut num_first_round_polys = num_state_polys + num_witness_polys + 1; // main gate quotient
    if has_custom_gate {
        num_first_round_polys += 1; // custom gate quotient
    }

    let num_lookup_polys = if has_lookup {
        // s(x), z(x), T0(x), T1(x), T2(x)
        1 + 1 + 1 + 1 + 1
    } else {
        0
    };
    let num_copy_permutation_polys = 1 + 1 + 1;
    let num_second_round_polys = num_lookup_polys + num_copy_permutation_polys;

    let max_num_polys = [num_setup_polys, num_first_round_polys, num_second_round_polys].into_iter().max().unwrap();

    (num_setup_polys, num_first_round_polys, num_second_round_polys, max_num_polys)
}
pub fn interpolate_union_set<F: PrimeField>(evals: Vec<F>, evals_shifted: Vec<F>, interpolation_set_size: usize, h: (F, Option<F>), omega: F, requires_shifted_opening: bool) -> Vec<F> {
    if requires_shifted_opening {
        assert_eq!(evals.len(), evals_shifted.len());
    }

    let mut c_evals = vec![];
    let mut c_evals_shifted = vec![];
    let (h, h_shifted) = h;

    let mut current_omega = h;
    let mut current_omega_shifted = h_shifted.unwrap_or(F::zero());

    for _ in 0..interpolation_set_size {
        let sum = horner_evaluation(&evals, current_omega);
        c_evals.push(sum);
        current_omega.mul_assign(&omega);

        if requires_shifted_opening {
            let sum_shifted = horner_evaluation(&evals_shifted, current_omega_shifted);
            c_evals_shifted.push(sum_shifted);
            current_omega_shifted.mul_assign(&omega);
        }
    }

    let monomial = if requires_shifted_opening {
        assert_eq!(c_evals.len(), c_evals_shifted.len());

        let monomial = interpolate_r_monomial_from_union_set(&c_evals, &c_evals_shifted, h, h_shifted.expect("h shifted"), omega, interpolation_set_size);
        assert_eq!(monomial.len(), 2 * interpolation_set_size);
        monomial
    } else {
        let monomial = interpolate_r_monomial(&c_evals, interpolation_set_size, h, omega);
        assert_eq!(monomial.len(), interpolation_set_size);
        monomial
    };

    monomial
}

pub(crate) fn construct_r_monomials<F: PrimeField>(
    setup_evaluations: &SetupEvaluations<F>,
    first_round_evaluations: &FirstRoundEvaluations<F>,
    second_round_evaluations: &SecondRoundEvaluations<F>,
    recomputed_quotient_evaluations: &RecomptuedQuotientEvaluations<F>,
    h0: (F, Option<F>),
    h1: (F, Option<F>),
    h2: (F, F),
) -> (Vec<F>, Vec<F>, Vec<F>) {
    let interpolation_size_of_setup = setup_evaluations.interpolation_size();
    assert!(interpolation_size_of_setup.is_power_of_two());
    let interpolation_size_of_first_round = first_round_evaluations.interpolation_size();
    assert!(interpolation_size_of_first_round.is_power_of_two());
    let interpolation_size_of_second_round = second_round_evaluations.interpolation_size();
    assert_eq!(interpolation_size_of_second_round, 3);

    let [setup_omega, first_round_omega, second_round_omega] = compute_generators(interpolation_size_of_setup, interpolation_size_of_first_round, interpolation_size_of_second_round);

    let (evals, evals_shifted) = setup_evaluations.flatten();
    let setup_r_monomial = interpolate_union_set(
        evals,
        evals_shifted,
        interpolation_size_of_setup,
        h0,
        setup_omega,
        setup_evaluations.requires_opening_at_shifted_point(),
    );

    let (mut evals, evals_shifted) = first_round_evaluations.flatten();
    evals.push(recomputed_quotient_evaluations.main_gate_quotient_at_z);
    assert!(recomputed_quotient_evaluations.custom_gate_quotient_at_z.is_none());

    let first_round_r_monomial = interpolate_union_set(
        evals,
        evals_shifted,
        interpolation_size_of_first_round,
        h1,
        first_round_omega,
        first_round_evaluations.requires_opening_at_shifted_point(),
    );
    let h2 = (h2.0, Some(h2.1));
    let (mut evals, evals_shifted) = second_round_evaluations.flatten();
    evals.push(recomputed_quotient_evaluations.copy_permutation_first_quotient_at_z);
    evals.push(recomputed_quotient_evaluations.copy_permutation_second_quotient_at_z);
    assert!(recomputed_quotient_evaluations.lookup_first_quotient_at_z.is_none());
    assert!(recomputed_quotient_evaluations.lookup_second_quotient_at_z.is_none());
    assert!(recomputed_quotient_evaluations.lookup_third_quotient_at_z.is_none());
    assert_eq!(second_round_evaluations.interpolation_size(), 3);
    let second_round_r_monomial = interpolate_union_set(evals, evals_shifted, interpolation_size_of_second_round, h2, second_round_omega, true);

    (setup_r_monomial, first_round_r_monomial, second_round_r_monomial)
}

pub fn main_gate_quotient_degree<E: Engine>(sorted_gates: &[Box<dyn GateInternal<E>>]) -> usize {
    let has_custom_gate = sorted_gates.len() > 1;

    if has_custom_gate {
        sorted_gates[0].degree()
    } else {
        sorted_gates[0].degree() - 1
    }
}

pub fn custom_gate_quotient_degree<E: Engine>(sorted_gates: &[Box<dyn GateInternal<E>>]) -> usize {
    assert!(sorted_gates.len() > 1);
    sorted_gates[0].degree()
}

pub(crate) fn construct_set_difference_monomials<F: PrimeField>(
    z: F,
    z_omega: F,
    interpolation_size_of_setup: usize,
    interpolation_size_of_first_round: usize,
    interpolation_size_of_second_round: usize,
    first_round_requires_opening_at_shifted_point: bool,
) -> [Vec<(usize, F)>; 4] {
    let first_round_part = if first_round_requires_opening_at_shifted_point {
        vec![(interpolation_size_of_first_round, z), (interpolation_size_of_first_round, z_omega)]
    } else {
        vec![(interpolation_size_of_first_round, z)]
    };
    // Z_{T\S0}(x)
    let mut sparse_polys_for_setup = vec![(interpolation_size_of_second_round, z), (interpolation_size_of_second_round, z_omega)];
    sparse_polys_for_setup.extend_from_slice(&first_round_part);

    // Z_{T\S1}(x)
    let sparse_polys_for_first_round = vec![(interpolation_size_of_setup, z), (interpolation_size_of_second_round, z), (interpolation_size_of_second_round, z_omega)];

    // Z_{T\S2}(x)
    let mut sparse_polys_for_second_round = vec![(interpolation_size_of_setup, z)];
    sparse_polys_for_second_round.extend_from_slice(&first_round_part);

    // Z_T(x)
    let mut sparse_polys = vec![(interpolation_size_of_setup, z), (interpolation_size_of_second_round, z), (interpolation_size_of_second_round, z_omega)];

    sparse_polys.extend_from_slice(&first_round_part);

    [sparse_polys_for_setup, sparse_polys_for_first_round, sparse_polys_for_second_round, sparse_polys]
}

pub fn binop_over_slices<F: PrimeField, B: FieldBinop<F>>(worker: &Worker, binop: &B, dest: &mut [F], source: &[F]) {
    assert_eq!(dest.len(), source.len());
    worker.scope(dest.len(), |scope, chunk| {
        for (dest, source) in dest.chunks_mut(chunk).zip(source.chunks(chunk)) {
            scope.spawn(move |_| {
                for (dest, source) in dest.iter_mut().zip(source.iter()) {
                    binop.apply(dest, source);
                }
            });
        }
    });
}

pub fn calculate_lagrange_poly<F: PrimeField>(worker: &Worker, poly_size: usize, poly_number: usize) -> Result<Polynomial<F, Coefficients>, SynthesisError> {
    assert!(poly_size.is_power_of_two());
    assert!(poly_number < poly_size);

    let mut poly = Polynomial::<F, Values>::from_values(vec![F::zero(); poly_size])?;

    poly.as_mut()[poly_number] = F::one();

    Ok(poly.ifft(&worker))
}

pub fn evaluate_vanishing_polynomial_of_degree_on_domain_size<F: PrimeField>(
    vanishing_degree: u64,
    coset_factor: &F,
    domain_size: u64,
    worker: &Worker,
) -> Result<Polynomial<F, Values>, SynthesisError> {
    let domain = Domain::<F>::new_for_size(domain_size)?;
    let domain_generator = domain.generator;

    let coset_factor = coset_factor.pow(&[vanishing_degree]);

    let domain_generator_in_vanishing_power = domain_generator.pow(&[vanishing_degree]);

    let mut minus_one = F::one();
    minus_one.negate();

    let mut result = vec![minus_one; domain.size as usize];

    worker.scope(result.len(), |scope, chunk_size| {
        for (chunk_id, chunk) in result.chunks_mut(chunk_size).enumerate() {
            scope.spawn(move |_| {
                let start = chunk_id * chunk_size;
                let mut pow = domain_generator_in_vanishing_power.pow(&[start as u64]);
                pow.mul_assign(&coset_factor);
                for el in chunk.iter_mut() {
                    el.add_assign(&pow);
                    pow.mul_assign(&domain_generator_in_vanishing_power);
                }
            });
        }
    });

    Polynomial::from_values(result)
}

pub fn materialize_domain_elements_with_natural_enumeration<F: PrimeField>(domain: &Domain<F>, worker: &Worker) -> Vec<F> {
    let mut values = vec![F::zero(); domain.size as usize];
    let generator = domain.generator;

    worker.scope(values.len(), |scope, chunk| {
        for (i, values) in values.chunks_mut(chunk).enumerate() {
            scope.spawn(move |_| {
                let mut current_power = generator.pow(&[(i * chunk) as u64]);

                for p in values {
                    *p = current_power;
                    current_power.mul_assign(&generator);
                }
            });
        }
    });

    values
}

pub fn commit_point_as_xy<E: Engine, T: Transcript<E::Fr>>(transcript: &mut T, point: &E::G1Affine) {
    if point.is_zero() {
        transcript.commit_fe(&E::Fq::zero());
        transcript.commit_fe(&E::Fq::zero());
    } else {
        let (x, y) = point.into_xy_unchecked();
        transcript.commit_fe(&x);
        transcript.commit_fe(&y);
    }
}

pub fn compute_cubic_root_of_domain<F: PrimeField>(domain_size: u64) -> F {
    let modulus = repr_to_biguint::<F>(&F::char());
    let domain = Domain::<F>::new_for_size(domain_size).unwrap();
    let omega = fe_to_biguint(&domain.generator);
    let root = compute_cube_root_via_tonelli_shanks(omega, modulus).expect("cube root for {domain_size}");

    biguint_to_fe(root)
}

pub(crate) fn compute_cube_root_via_tonelli_shanks(n: BigUint, p: BigUint) -> Option<BigUint> {
    // Check if n is a cubic residue
    let exp = (&p - 1u32) / 3u32;
    if n.modpow(&exp, &p) != BigUint::one() {
        return None; // n is not a cubic residue
    }

    // Find Q and S, where p - 1 = 3^S * Q
    let mut q = &p - 1u32;
    let mut s = 0u32;
    while &q % 3u32 == BigUint::zero() {
        s += 1;
        q /= 3u32;
    }

    // Find a cubic non-residue z
    let mut z = BigUint::from(2u32);
    while z.modpow(&exp, &p) == BigUint::one() {
        z += 1u32;
    }

    let mut m = s;
    let mut c = z.modpow(&q, &p);
    let mut t = n.modpow(&q, &p);
    let mut r = n.modpow(&((q + 1u32) / 3u32), &p);

    loop {
        if t == BigUint::one() {
            assert_eq!(n, (&r * &r * &r) % &p);
            return Some(r);
        }

        let mut i = 0u32;
        let mut t2 = t.clone();
        while t2 != BigUint::one() && i < m {
            t2 = (&t2 * &t2 * &t2) % &p;
            i += 1;
        }

        if i == m {
            return None; // n is not a cubic residue
        }

        let b = c.modpow(&BigUint::from(3u32.pow(m - i - 1)), &p);
        m = i;
        c = (&b * &b * &b) % &p;
        t = (&t * &c) % &p;
        r = (&r * &b) % &p;
    }
}

pub(crate) fn compute_3rd_roots_of_unity<F: PrimeField>() -> F {
    let g = Fr::multiplicative_generator();
    let modulus = repr_to_biguint::<Fr>(&Fr::char());
    let p_minus_one = &modulus - &BigUint::one();
    let three = BigUint::from(3u8);
    assert_eq!(&p_minus_one % &three, BigUint::zero());
    let three_inv = three.modinv(&modulus).expect("inverse of 3");
    let exp = (&p_minus_one * &three_inv) % &modulus;

    let g = fe_to_biguint(&g);
    let omega = g.modpow(&exp, &modulus);
    assert_eq!(omega.modpow(&three, &modulus), BigUint::one());

    let omega_squared = &omega * &omega;
    assert_eq!(omega_squared.modpow(&three, &modulus), BigUint::one());

    biguint_to_fe(omega)
}

pub(crate) fn hardcoded_cube_root_of_generator_of_trace_domain<F: PrimeField>() -> F {
    let root: F = biguint_to_fe(BigUint::parse_bytes(b"925f0bd364638ec3084b45fc27895f8f3f6f079096600fe946c8e9db9a47124", 16).unwrap());

    let omega: F = Domain::new_for_size(1 << L1_VERIFIER_DOMAIN_SIZE_LOG).unwrap().generator;
    if SANITY_CHECK {
        let mut actual = root;
        actual.mul_assign(&root);
        actual.mul_assign(&root);
        assert_eq!(omega, actual);
    }

    root
}
pub(crate) fn hardcoded_generator_of_3rd_roots_of_unity<F: PrimeField>() -> F {
    biguint_to_fe(BigUint::parse_bytes(b"0000000000000000b3c4d79d41a917585bfc41088d8daaa78b17ea66b99c90dd", 16).unwrap())
}

pub struct FflonkTestCircuit;

impl Circuit<Bn256> for FflonkTestCircuit {
    type MainGate = NaiveMainGate;

    fn synthesize<CS: circuit_definitions::snark_wrapper::franklin_crypto::bellman::plonk::better_better_cs::cs::ConstraintSystem<Bn256> + 'static>(
        &self,
        cs: &mut CS,
    ) -> Result<(), circuit_definitions::snark_wrapper::franklin_crypto::bellman::SynthesisError> {
        use circuit_definitions::snark_wrapper::franklin_crypto::bellman::Field;
        use circuit_definitions::snark_wrapper::franklin_crypto::plonk::circuit::allocated_num::Num;
        let a = Fr::from_str(&65.to_string()).unwrap();
        let b = Fr::from_str(&66.to_string()).unwrap();
        let mut c = a;
        c.add_assign(&b);

        let a_var = Num::alloc(cs, Some(a))?;
        let b_var = Num::alloc(cs, Some(b))?;
        let c_var = Num::alloc(cs, Some(c))?;

        for _ in 0..1 << 5 {
            let mut lc = LinearCombination::zero();
            lc.add_assign_number_with_coeff(&a_var, Fr::one());
            lc.add_assign_number_with_coeff(&b_var, Fr::one());
            let mut minus_one = Fr::one();
            minus_one.negate();
            lc.add_assign_number_with_coeff(&c_var, minus_one);

            let _ = lc.into_num(cs)?;
        }

        let _input = cs.alloc_input(|| Ok(Fr::one()))?;

        Ok(())
    }
}

#[test]
fn test_poly_combination() {
    use rand::{Rng, SeedableRng, XorShiftRng};
    fn init_rng() -> XorShiftRng {
        XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654])
    }

    let worker = Worker::new();
    let rng = &mut init_rng();

    let domain_size = 1 << 10;
    let num_polys = 2u64;
    let padded_num_polys = num_polys.next_power_of_two();

    let mut all_monomials = vec![];

    for _ in 0..num_polys {
        let coeffs: Vec<Fr> = (0..domain_size).map(|_| rng.gen()).collect();
        let monomial = Polynomial::from_coeffs(coeffs).expect("monomial");
        all_monomials.push(monomial);
    }
    let domain = Domain::new_for_size(num_polys).expect("domain");
    let omega: Fr = domain.generator;
    assert_eq!(omega.pow(&[padded_num_polys]), Fr::one());

    let point = rng.gen::<Fr>();
    let size = domain.size;
    let omega: Fr = domain.generator;
    let mut evaluation_set = vec![point; num_polys as usize];
    let mut cur = omega;
    for p in evaluation_set.iter_mut().skip(1) {
        p.mul_assign(&cur);
        cur.mul_assign(&omega);
    }

    evaluation_set.iter().for_each(|el| assert_eq!(el.pow(&[size]), point.pow(&[size])));

    let z = point.pow(&[padded_num_polys]);

    let all_evaluations: Vec<Fr> = all_monomials.iter().map(|p| p.evaluate_at(&worker, z)).collect();

    let combined_monomial = combine_monomials(&all_monomials.iter().map(|p| p).collect::<Vec<_>>(), domain_size).expect("combined monomial");
    assert_eq!(combined_monomial.size(), padded_num_polys as usize * domain_size);
    // f(x) = f0(x^2) + x*f1(x^2)
    // we have f(x^2) and f1(x^2)
    // and need to compute
    // f(r) = f0(r^2) + r*f1(r^2)
    // let mut previous_result = F::zero();
    let mut values_over_set = vec![];
    for (_, eval_point) in evaluation_set.iter().enumerate() {
        let expected = combined_monomial.evaluate_at(&worker, eval_point.clone());
        assert_eq!(eval_point.pow(&[padded_num_polys as u64]), point.pow(&[padded_num_polys as u64]));
        let mut actual = Fr::zero();
        let mut current_factor = Fr::one();
        for p in all_evaluations.iter() {
            let mut contrib = p.clone();
            contrib.mul_assign(&current_factor);
            actual.add_assign(&contrib);
            current_factor.mul_assign(&eval_point);
        }
        assert_eq!(actual, expected);
        values_over_set.push(expected);
    }
}

#[test]
fn test_divide_by_higher_degree() {
    use rand::{Rng, SeedableRng, XorShiftRng};
    fn init_rng() -> XorShiftRng {
        XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654])
    }

    let worker = Worker::new();
    let rng = &mut init_rng();

    let poly_degree = (1 << 20) - 1;
    let first_degree = 32;
    let sparse_degree = 8;

    let total_degree = poly_degree + first_degree + sparse_degree;
    println!(
        "poly degree {}\nfirst degree {}\nsecond degree {}\n total degree {}",
        poly_degree, first_degree, sparse_degree, total_degree
    );

    let z: Fr = rng.gen();
    let mut minus_z = z.clone();
    minus_z.negate();

    let quotient_coeffs: Vec<_> = (0..total_degree).map(|_| rng.gen()).collect();
    let quotient_poly = Polynomial::from_coeffs_unpadded(quotient_coeffs).unwrap();

    let mut sparse_coeffs = vec![Fr::zero(); sparse_degree + 1];
    sparse_coeffs[0] = minus_z;
    sparse_coeffs[sparse_degree] = Fr::one();

    let dividend_poly = multiply_monomial_with_sparse_poly(&worker, &quotient_poly, sparse_degree, z).unwrap();

    let actual_quotient = divide_by_higher_degree_term(dividend_poly.as_ref(), sparse_degree, z);
    let actual_quotient_poly = Polynomial::from_coeffs_unpadded(actual_quotient).unwrap();

    assert_eq!(quotient_poly.size(), actual_quotient_poly.size());
    assert_eq!(quotient_poly.as_ref(), actual_quotient_poly.as_ref());
}

#[test]
fn test_alternative_monomial_multiplication() -> Result<(), SynthesisError> {
    fn multiply_polys_in_monomial<F: PrimeField>(poly1: &Polynomial<F, Coefficients>, poly2: &Polynomial<F, Coefficients>) -> Result<Polynomial<F, Coefficients>, SynthesisError> {
        let coeffs = multiply_monomials(poly1.as_ref(), poly2.as_ref());

        Polynomial::from_coeffs_unpadded(coeffs)
    }

    use rand::{Rng, SeedableRng, XorShiftRng};
    fn init_rng() -> XorShiftRng {
        XorShiftRng::from_seed([0x3dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654])
    }

    let worker = Worker::new();
    let rng = &mut init_rng();

    let poly_degree = (1 << 20) - 1;
    let first_degree = 32;
    let second_degree = 8;

    let total_degree = poly_degree + first_degree + second_degree;
    println!(
        "poly degree {}\nfirst degree {}\nsecond degree {}\n total degree {}",
        poly_degree, first_degree, second_degree, total_degree
    );

    let z: Fr = rng.gen();
    let mut minus_z = z.clone();
    minus_z.negate();

    let z2: Fr = rng.gen();
    let mut minus_z2 = z2.clone();
    minus_z2.negate();

    let coeffs: Vec<Fr> = (0..(poly_degree + 1)).map(|_| rng.gen()).collect();
    let poly = Polynomial::from_coeffs_unpadded(coeffs)?;
    let mut first_coeffs: Vec<Fr> = vec![Fr::zero(); first_degree + 1];
    first_coeffs[0] = minus_z;
    first_coeffs[first_degree] = Fr::one();
    let mut second_coeffs: Vec<Fr> = vec![Fr::zero(); second_degree + 1];
    second_coeffs[0] = minus_z2;
    second_coeffs[second_degree] = Fr::one();

    let first = Polynomial::from_coeffs_unpadded(first_coeffs)?;
    let second = Polynomial::from_coeffs_unpadded(second_coeffs)?;
    let now = std::time::Instant::now();
    let tmp = multiply_polys_in_monomial(&poly, &first)?;
    println!("first mul degree {}", tmp.size() - 1);
    let expected = multiply_polys_in_monomial(&tmp, &second)?;
    println!("naive elapsed {}", now.elapsed().as_millis());
    println!("second mul degree {}", expected.size() - 1);

    let pairs = vec![(first_degree, z), (second_degree, z2)];
    let now = std::time::Instant::now();
    let actual = multiply_monomial_with_multiple_sparse_polys(&worker, &poly, &pairs)?;
    println!("alternative elapsed {}", now.elapsed().as_millis());
    assert_eq!(expected.size(), actual.size());
    assert_eq!(expected.as_ref(), actual.as_ref());

    Ok(())
}

#[test]
fn test_cube_root() {
    use super::*;
    use franklin_crypto::bellman::plonk::domains::Domain;
    use franklin_crypto::plonk::circuit::bigint::{fe_to_biguint, repr_to_biguint};
    let modulus = repr_to_biguint::<Fr>(&Fr::char());
    let omega = Domain::<Fr>::new_for_size(1 << L1_VERIFIER_DOMAIN_SIZE_LOG).unwrap().generator;
    let omega_big = fe_to_biguint(&omega);
    let cube_root = compute_cube_root_via_tonelli_shanks(omega_big.clone(), modulus.clone()).expect("cube root");
    println!("cube_root of size 2^{} domain generator {:x}", L1_VERIFIER_DOMAIN_SIZE_LOG, cube_root);
    assert_eq!((&cube_root * &cube_root * &cube_root) % &modulus, omega_big);

    let cube_root = hardcoded_cube_root_of_generator_of_trace_domain::<Fr>();
    let mut actual = cube_root;
    actual.mul_assign(&cube_root);
    actual.mul_assign(&cube_root);
    assert_eq!(omega, actual);

    let g = compute_3rd_roots_of_unity::<Fr>();
    println!("3rd roots of unity {g}");
    use franklin_crypto::bellman::Field;
    let mut actual = g;
    actual.mul_assign(&g);
    actual.mul_assign(&g);
    assert_eq!(actual, Fr::one());

    let mut current = cube_root;
    assert_eq!(current.pow(&[3]), omega);
    current.mul_assign(&g);
    assert_eq!(current.pow(&[3]), omega);
    current.mul_assign(&g);
    assert_eq!(current.pow(&[3]), omega);

    let r = Fr::from_str(&65.to_string()).unwrap();

    let z = r.pow(&[24]);
    let mut z_omega = z;
    z_omega.mul_assign(&omega);

    let h = r.pow(&[8]);
    assert_eq!(h.pow(&[3]), z);
    let mut h_s: Fr = hardcoded_cube_root_of_generator_of_trace_domain();
    h_s.mul_assign(&h);
    assert_eq!(h_s.pow(&[3]), z_omega);
}
