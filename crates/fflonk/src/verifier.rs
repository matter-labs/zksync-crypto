use super::*;
use bellman::plonk::better_cs::keys::read_curve_affine;
use bellman::plonk::better_cs::keys::read_curve_affine_vector;
use bellman::plonk::better_cs::keys::read_fr_vec;
use bellman::plonk::better_cs::keys::write_curve_affine;
use bellman::plonk::better_cs::keys::write_curve_affine_vec;
use bellman::plonk::better_cs::keys::write_fr_vec;
use byteorder::BigEndian;
use byteorder::ReadBytesExt;
use byteorder::WriteBytesExt;
use std::io::{Read, Write};

#[derive(Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct FflonkVerificationKey<E: Engine, C: Circuit<E>> {
    pub n: usize,
    pub c0: E::G1Affine,
    pub num_inputs: usize,
    // TODO
    pub num_state_polys: usize,
    pub num_witness_polys: usize,
    pub total_lookup_entries_length: usize,
    pub non_residues: Vec<E::Fr>,
    pub g2_elements: [E::G2Affine; 2],

    #[serde(skip_serializing, skip_deserializing, default)]
    #[serde(bound(serialize = ""))]
    #[serde(bound(deserialize = ""))]
    _marker: std::marker::PhantomData<C>,
}

impl<E: Engine, C: Circuit<E>> FflonkVerificationKey<E, C> {
    pub fn new(n: usize, c0: E::G1Affine, num_inputs: usize, num_state_polys: usize, num_witness_polys: usize, total_lookup_entries_length: usize, g2_elements: [E::G2Affine; 2]) -> Self {
        let non_residues = make_non_residues(num_state_polys - 1);
        FflonkVerificationKey {
            n,
            c0,
            num_inputs,
            num_state_polys,
            num_witness_polys,
            non_residues,
            g2_elements,
            total_lookup_entries_length,
            _marker: std::marker::PhantomData,
        }
    }
    pub fn from_setup(setup: &FflonkSetup<E, C>, crs: &Crs<E, CrsForMonomialForm>) -> Result<Self, SynthesisError> {
        let FflonkSetup { original_setup, c0_commitment: c0 } = setup;
        Ok(Self {
            n: original_setup.n,
            num_inputs: original_setup.num_inputs,
            c0: c0.clone(),
            num_state_polys: original_setup.state_width,
            num_witness_polys: original_setup.num_witness_polys,
            total_lookup_entries_length: original_setup.total_lookup_entries_length,
            non_residues: original_setup.non_residues.clone(),
            g2_elements: [crs.g2_monomial_bases[0], crs.g2_monomial_bases[1]],
            _marker: std::marker::PhantomData,
        })
    }

    pub fn read<R: Read>(mut src: R) -> Result<Self, std::io::Error> {
        let n = src.read_u64::<BigEndian>()? as usize;
        let c0 = read_curve_affine(&mut src)?;
        let num_inputs = src.read_u64::<BigEndian>()? as usize;
        let num_state_polys = src.read_u64::<BigEndian>()? as usize;
        let num_witness_polys = src.read_u64::<BigEndian>()? as usize;
        let total_lookup_entries_length = src.read_u64::<BigEndian>()? as usize;
        let non_residues = read_fr_vec(&mut src)?;
        let g2_elements = read_curve_affine_vector(&mut src)?;

        Ok(Self {
            n,
            c0,
            num_inputs,
            num_state_polys,
            num_witness_polys,
            total_lookup_entries_length,
            non_residues,
            g2_elements: g2_elements.try_into().unwrap(),
            _marker: std::marker::PhantomData,
        })
    }

    pub fn write<W: Write>(&self, mut dst: W) -> Result<(), std::io::Error> {
        dst.write_u64::<BigEndian>(self.n as u64)?;
        write_curve_affine(&self.c0, &mut dst)?;
        dst.write_u64::<BigEndian>(self.num_inputs as u64)?;
        dst.write_u64::<BigEndian>(self.num_state_polys as u64)?;
        dst.write_u64::<BigEndian>(self.num_witness_polys as u64)?;
        dst.write_u64::<BigEndian>(self.total_lookup_entries_length as u64)?;
        write_fr_vec(&self.non_residues, &mut dst)?;
        write_curve_affine_vec(&self.g2_elements, &mut dst)?;

        Ok(())
    }
}

pub fn verify<E: Engine, C: Circuit<E>, T: Transcript<E::Fr>>(
    vk: &FflonkVerificationKey<E, C>,
    proof: &FflonkProof<E, C>,
    transcript_params: Option<T::InitializationParameters>,
) -> Result<bool, SynthesisError> {
    let mut transcript = if let Some(params) = transcript_params { T::new_from_params(params) } else { T::new() };

    let sorted_gates = sorted_gates_from_circuit_definitions::<_, C>();
    assert!(sorted_gates.len() > 0);

    let FflonkProof { evaluations, inputs, .. } = proof;

    let n = vk.n;
    let domain_size = n + 1;
    assert!(domain_size.is_power_of_two());
    assert!(domain_size.trailing_zeros() <= 23);

    let (num_setup_polys, num_first_round_polys, num_second_round_polys, _) = num_system_polys_from_vk(vk);
    let num_state_polys = vk.num_state_polys;
    println!("num state polys: {}", vk.num_state_polys);
    println!("num witness polys: {}", vk.num_witness_polys);
    println!("num first round polys: {}", num_setup_polys);
    println!("num second round polys: {}", num_first_round_polys);
    println!("num third round polys: {}", num_second_round_polys);
    let non_residues = vk.non_residues.clone();

    let has_lookup = vk.total_lookup_entries_length > 0;
    let has_custom_gate = sorted_gates.len() > 1;

    // Commit data in the transcript  then get challenges
    // Note that at this point transcript only has public inputs
    // But luckily prover doesn't need any randomness in the first round
    // so that prover has no control over the values because quotients are
    // seperated(there is no quotient aggregation neither in this round nor all rounds)
    if proof.inputs.len() != vk.num_inputs {
        return Ok(false);
    }
    

    for inp in proof.inputs.iter() {
        transcript.commit_field_element(inp);
    }
    // commit first round commitment: setup
    commit_point_as_xy::<E, T>(&mut transcript, &vk.c0);

    // commit second round commitment: witnesses
    commit_point_as_xy::<E, T>(&mut transcript, &proof.commitments[0]);

    // copy-permutation challenges
    let beta_for_copy_permutation = transcript.get_challenge();
    let gamma_for_copy_permutation = transcript.get_challenge();
    let (eta_for_lookup, beta_for_lookup, gamma_for_lookup) = if has_lookup {
        // lookup challenges
        let eta = transcript.get_challenge();
        let beta = transcript.get_challenge();
        let gamma = transcript.get_challenge();

        (Some(eta), Some(beta), Some(gamma))
    } else {
        (None, None, None)
    };
    commit_point_as_xy::<E, T>(&mut transcript, &proof.commitments[1]);
    // evaluation challenge
    let r = transcript.get_challenge();
    // commit evaluations
    for eval in proof.evaluations.iter() {
        transcript.commit_field_element(eval);
    }
    let c1 = proof.commitments[0];
    let c2 = proof.commitments[1];
    let w = proof.commitments[2];
    let w_prime = proof.commitments[3];
    // opening linearization challenge
    let alpha = transcript.get_challenge();
    commit_point_as_xy::<E, T>(&mut transcript, &w);

    // last opening challenge
    let y = transcript.get_challenge();

    // all system polynomials will be evaluated at z
    // then combined polynomials will be opened at h_i = r^power_i
    // then it becomes e.g C_i(x) = f_0(x^2) + x*f(x^2) in case of two polynomials
    assert_eq!(num_second_round_polys, 3);
    let power = lcm(&[num_setup_polys.next_power_of_two(), num_first_round_polys.next_power_of_two(), num_second_round_polys]);
    println!("LCM {power}");
    let z = r.pow(&[power as u64]);
    let evaluation_offsets = EvaluationOffsets::from_vk(vk);

    let mut all_gates_iter = sorted_gates.clone().into_iter();
    let main_gate_internal = all_gates_iter.next().unwrap();
    assert!(&C::MainGate::default().into_internal() == &main_gate_internal);
    let main_gate = C::MainGate::default();

    let custom_gate_name = if has_custom_gate {
        let custom_gate = all_gates_iter.next().unwrap();
        Some(custom_gate.name())
    } else {
        None
    };

    let recomputed_quotients = recompute_quotients_from_evaluations(
        evaluations,
        &evaluation_offsets,
        inputs,
        z,
        domain_size,
        main_gate.name(),
        beta_for_copy_permutation,
        gamma_for_copy_permutation,
        &non_residues,
        num_state_polys,
        custom_gate_name,
        eta_for_lookup,
        beta_for_lookup,
        gamma_for_lookup,
    );

    let setup_requires_opening_at_shifted_point = has_lookup; // TODO
    let first_round_requires_opening_at_shifted_point = requires_trace_polys_opening_at_shifted_point(main_gate_internal);
    aggregate_points_and_check_pairing(
        vk,
        r,
        z,
        alpha,
        y,
        num_setup_polys,
        num_first_round_polys,
        num_second_round_polys,
        proof.evaluations.to_vec(),
        recomputed_quotients,
        &proof.lagrange_basis_inverses,
        c1,
        c2,
        w,
        w_prime,
        setup_requires_opening_at_shifted_point,
        first_round_requires_opening_at_shifted_point,
    )
}

fn aggregate_points_and_check_pairing<E: Engine, C: Circuit<E>>(
    vk: &FflonkVerificationKey<E, C>,
    r: E::Fr,
    z: E::Fr,
    alpha: E::Fr,
    y: E::Fr,
    num_setup_polys: usize,
    num_first_round_polys: usize,
    num_second_round_polys: usize,
    all_evaluations: Vec<E::Fr>,
    recomputed_quotient_evaluations: RecomptuedQuotientEvaluations<E::Fr>,
    lagrange_basis_inverses: &[E::Fr],
    c1: E::G1Affine,
    c2: E::G1Affine,
    w: E::G1Affine,
    w_prime: E::G1Affine,
    setup_requires_opening_at_shifted_point: bool,
    first_round_requires_opening_at_shifted_point: bool,
) -> Result<bool, SynthesisError> {
    let domain_size = vk.n + 1;
    assert!(domain_size.is_power_of_two());

    // Now it is time to combine all rounds in a combined poly
    // C0(x) = f0(X^k0) + X*f1(X^k0) + X*f_{k0-1}(X^k0) where k0 is the total number of

    // Since fft-style combination requires roots of unity, we need total number of polyn
    // be power of two, in our case none of them are power of two so that we apply a padding
    // with zeroes here.
    // Compute h0 = r^(power/k0) here, it will be plugged in the C0(X)
    let interpolation_size_of_setup = num_setup_polys.next_power_of_two();
    let interpolation_size_of_first_round = num_first_round_polys.next_power_of_two();
    let interpolation_size_of_second_round = num_second_round_polys;

    let power = lcm(&[interpolation_size_of_setup, interpolation_size_of_first_round, interpolation_size_of_second_round]);
    let mut z_omega = z;
    let omega = Domain::new_for_size(domain_size as u64).unwrap().generator;
    z_omega.mul_assign(&omega);

    let (h0, h1, h2) = compute_opening_points(
        r,
        z,
        z_omega,
        power,
        interpolation_size_of_setup,
        interpolation_size_of_first_round,
        interpolation_size_of_second_round,
        domain_size,
        setup_requires_opening_at_shifted_point,
        first_round_requires_opening_at_shifted_point,
    );

    // In order to verify openings, we should construct r(x) such that r(h_i*w_i) = C_i(h_i*w_i)
    // We have all the necessary evaluations of the system polynomial and can reconstruct r_i(x)
    // from those; e.g in case of copy-permuation combined poly we have
    // - C_i(h0*w0) = Z(h0_w0) + (h0*w0)*T1(h0*w0) + ((h0*w0)^2)*(T2(h0*w0))
    // - C_i(h0*w1) = Z(h0_w1) + (h0*w1)*T1(h0*w1) + ((h0*w1)^2)*(T2(h0*w1))
    // ...
    // - C_i(h0*w^{k-1}) = Z(h0_w1^{k-1}) + (h0*w1^{k-1})*T1(h0*w1^{k-1}) + ((h0*w1^{k-1})^2)*(T2(h0*w1^{k-1}))

    // Now openings
    // f(x) = Z_{T\S0}(x)(C0(x)- r0(x)) + alpha*(Z_{T\S1}(x)*(C1(x)- r1(x))) + alpha^2*(Z_{T\S2}*(C2(x)- r2(x)))
    // Note that, in our case set differences(Z_T\{S_i}) are:
    // - Z_{T\S0}(x): (X^k1-z)*(X^k2-z)*(X^k2-z*w)
    // - Z_{T\S1}(x): (X^k0-z)*(X^k2-z)*(X^k2-z*w)
    // - Z_{T\S2}(x): (X^k0-z)*(X^k1-z) where
    // k0, k1, and k2 are number of the polynomials for setup, first and second
    // round respectively

    // W(x) = f(x) / Z_T(x) where Z_T(x) = (X^k0-z)(X^k1-z)*(X^k2-z)*(X^k2-z*w)
    // we need to check that
    // f(x) - W(x) * Z_T(x) = 0

    // L(x) = Z_{T\S0}(y)(C0(x)- r0(y)) + alpha*Z_{T\S1}(y)*(C1(x)- r1(y)) + alpha^2*Z_{T\S2}(y)*(C2(x)- r2(y)) - Z_T(x)*W(x)
    // W'(x) = L(x) / (Z_{T\S0}(y)*(x-y))
    // the identity check is reduced into following
    // L(x) - W'(x)*Z_{T\S0}(y)(x-y) == 0
    // verifier has commitments to the C_i(x) polynomials
    // verifer also recomputed r_i(y)
    // group constant and commitment parts
    // first prepare L(x)/Z_{T\S0}
    // C(x) = C0(x) + (alpha*Z_{T\S1}/Z_{T\S0})*C1(x) + (alpha^2*Z_{T\S2}/Z_{T\S0})*C2(x)
    // r(y) = r0(y) + (alpha*Z_{T\S1}/Z_{T\S0})*r1(y) + (alpha^2*Z_{T\S2}/Z_{T\S0})*r2(y)
    // now construct
    // L(x)/Z_{T\S0} = C(x) - r(y) - (Z_T(y)/Z_{T\S0})*W(x)
    // now check following identity
    // C(x) - r(y) - (Z_t(y)/Z_{T\S0}(y))*W(x) - W'(x)*(x-y)) = 0
    // [C(x)] - [r(y)*G1] - (Z_T(y)/Z_{T\S0}(y))*[W] - [(x-y)*W'] = 0
    // [C(x)] - [r(y)*G1] - (Z_T(y)/Z_{T\S0}(y))*[W] - [x*W'] + [y*W]' = 0
    // [C(x)] - [r(y)*G1] - (Z_T(y)/Z_{T\S0}(y))*[W] + [y*W'] - [x*W'] = 0
    // points with x will be multiplied in the exponent via pairing
    // so final pairing would ne
    // e([C(x)] - [r(y)*G1] - [Z_T(y)/Z_{T\S0}(y)*W] + [y*W'], G2)*e(-W', x*G2) = 1
    // F = [C]
    // E = [r(y)*G1]
    // J = [Z_T(y)*W]
    // e(F- E -J + [y*W'], G2[0]) * e(-W', x*G2[0]) = 1

    let mut alpha_squared = alpha;
    alpha_squared.mul_assign(&alpha);

    // Construct evaluations of C_i(x) polynomials using existing evaluations
    // of the system polynomials
    // Since evaluation sets are not constant, rather than lagrange interpolation
    // barycentric interpolation is utilized here.
    let precomputed_basis_evals = precompute_all_lagrange_basis_evaluations_from_inverses(
        lagrange_basis_inverses,
        interpolation_size_of_setup,
        interpolation_size_of_first_round,
        interpolation_size_of_second_round,
        h0,
        h1,
        h2,
        y,
        setup_requires_opening_at_shifted_point,
        first_round_requires_opening_at_shifted_point,
    );

    let [setup_r_at_y, mut first_round_r_at_y, mut second_round_r_at_y] = evaluate_r_polys_at_point_with_flattened_evals_and_precomputed_basis(
        all_evaluations,
        &recomputed_quotient_evaluations,
        num_setup_polys,
        num_first_round_polys,
        num_second_round_polys,
        h0,
        h1,
        h2,
        precomputed_basis_evals,
        setup_requires_opening_at_shifted_point,
        first_round_requires_opening_at_shifted_point,
    );

    let [
            sparse_polys_for_setup, // Z_{T\S0}(x)
            sparse_polys_for_first_round, // Z_{T\S1}(x)
            sparse_polys_for_second_round,// Z_{T\S2}(x)
            sparse_polys,// Z_T(x)
        ] = construct_set_difference_monomials(
            z,
            z_omega,
            interpolation_size_of_setup,
            interpolation_size_of_first_round,
            interpolation_size_of_second_round,
            first_round_requires_opening_at_shifted_point,
        );
    let sparse_polys_for_setup_at_y = evaluate_multiple_sparse_polys(sparse_polys_for_setup, y);
    let inv_sparse_polys_for_setup_at_y = sparse_polys_for_setup_at_y.inverse().unwrap();
    let sparse_polys_for_first_round_at_y = evaluate_multiple_sparse_polys(sparse_polys_for_first_round, y);
    let sparse_polys_for_second_round_at_y = evaluate_multiple_sparse_polys(sparse_polys_for_second_round, y);

    // r0(y)
    let mut aggregated_r_at_y = setup_r_at_y;

    // + (alpha*Z_{T\S1}(y)/Z_{T\S0}(y))*r1(y)
    first_round_r_at_y.mul_assign(&alpha);
    first_round_r_at_y.mul_assign(&sparse_polys_for_first_round_at_y);
    first_round_r_at_y.mul_assign(&inv_sparse_polys_for_setup_at_y);
    aggregated_r_at_y.add_assign(&first_round_r_at_y);

    // + (alpha^2*Z_{T\S2}(y)/Z_{T\S0}(y))*r2(y)
    second_round_r_at_y.mul_assign(&alpha_squared);
    second_round_r_at_y.mul_assign(&sparse_polys_for_second_round_at_y);
    second_round_r_at_y.mul_assign(&inv_sparse_polys_for_setup_at_y);
    aggregated_r_at_y.add_assign(&second_round_r_at_y);

    // C0
    let mut aggregated_commitment = vk.c0.into_projective();

    // + (alpha*Z_{T\S0}(y)/Z_{T\S0}(y))*C1
    let mut factor = alpha;
    factor.mul_assign(&sparse_polys_for_first_round_at_y);
    factor.mul_assign(&inv_sparse_polys_for_setup_at_y);
    let tmp = c1.mul(factor.into_repr());
    aggregated_commitment.add_assign(&tmp);

    // + (alpha^2*Z_{T\S0}(y)/Z_{T\S0}(y))*C2
    let mut factor = alpha_squared;
    factor.mul_assign(&sparse_polys_for_second_round_at_y);
    factor.mul_assign(&inv_sparse_polys_for_setup_at_y);
    let tmp = c2.mul(factor.into_repr());
    aggregated_commitment.add_assign(&tmp);

    let one = E::G1Affine::one();
    let e = one.mul(aggregated_r_at_y.into_repr());

    // (Z_T(y)/Z_{T\S0}(y))*W
    let mut z_t_at_y = evaluate_multiple_sparse_polys(sparse_polys, y);
    z_t_at_y.mul_assign(&inv_sparse_polys_for_setup_at_y);
    let j = w.mul(z_t_at_y.into_repr());

    // y*W'
    let w_prime_by_y = w_prime.mul(y.into_repr());

    aggregated_commitment.sub_assign(&e);
    aggregated_commitment.sub_assign(&j);
    aggregated_commitment.add_assign(&w_prime_by_y);
    let pair_with_generator = aggregated_commitment.into_affine();

    let mut pair_with_x = w_prime;
    pair_with_x.negate();

    let valid = E::final_exponentiation(&E::miller_loop(&[
        (&pair_with_generator.prepare(), &vk.g2_elements[0].prepare()),
        (&pair_with_x.prepare(), &vk.g2_elements[1].prepare()),
    ]))
    .ok_or(SynthesisError::Unsatisfiable)?
        == E::Fqk::one();

    Ok(valid)
}
