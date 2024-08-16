use std::io::{Read, Write};
use bellman::plonk::better_cs::keys::read_curve_affine;
use bellman::plonk::better_cs::keys::read_curve_affine_vector;
use bellman::plonk::better_cs::keys::read_fr_vec;
use bellman::plonk::better_cs::keys::write_curve_affine;
use bellman::plonk::better_cs::keys::write_curve_affine_vec;
use bellman::plonk::better_cs::keys::write_fr_vec;
use byteorder::BigEndian;
use byteorder::ReadBytesExt;
use byteorder::WriteBytesExt;
use super::*;

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
    pub fn from_setup(
        setup: &FflonkSetup<E, C>,
        crs: &Crs<E, CrsForMonomialForm>,
    ) -> Result<Self, SynthesisError> {
        let FflonkSetup {
            original_setup,
            c0_commitment: c0,
        } = setup;
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
        
        Ok(Self { n, c0, num_inputs, num_state_polys, num_witness_polys, 
            total_lookup_entries_length, non_residues, 
            g2_elements: g2_elements.try_into().unwrap(), 
            _marker: std::marker::PhantomData })
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
    let mut transcript = if let Some(params) = transcript_params {
        T::new_from_params(params)
    } else {
        T::new()
    };

    let sorted_gates = sorted_gates_from_circuit_definitions::<_, C>();
    assert!(sorted_gates.len() > 0);

    assert!((vk.n + 1).is_power_of_two());
    let required_domain_size = vk.n.next_power_of_two();

    let domain = Domain::<E::Fr>::new_for_size(required_domain_size as u64)?;

    let FflonkProof {
        setup_evaluations,
        first_round_evaluations,
        second_round_evaluations,
        inputs: public_inputs,
        ..
    } = proof;

    let FirstRoundEvaluations {
        trace_and_gate_evaluations,
    } = first_round_evaluations;

    let n = vk.n;
    let required_domain_size = n + 1;
    assert!(required_domain_size.is_power_of_two());
    assert!(required_domain_size.trailing_zeros() <= 23);

    let num_state_polys = vk.num_state_polys;
    let num_witness_polys = vk.num_witness_polys;
    let (num_setup_polys, num_first_round_polys, num_second_round_polys, _) =
        num_system_polys_from_vk(vk);
    println!("num state polys: {}", num_state_polys);
    println!("num witness polys: {}", num_witness_polys);
    println!("num şetup polys: {}", num_setup_polys);
    println!("num first round polys: {}", num_first_round_polys);
    println!("num second round polys: {}", num_second_round_polys);    

    let non_residues = vk.non_residues.clone();

    let has_lookup = vk.total_lookup_entries_length > 0;
    let has_custom_gate = sorted_gates.len() > 1;

    // Commit data in the transcript  then get challenges
    // Note that at this point transcript only has public inputs
    // But luckily prover doesn't need any randomness in the first round
    // so that prover has no control over the values because quotients are
    // seperated(there is no quotient aggregation neither in this round nor all rounds)
    for inp in proof.inputs.iter() {
        transcript.commit_field_element(inp);
    }
    // commit first round commitment: setup
    commit_point_as_xy::<E, T>(&mut transcript, &vk.c0);

    // commit second round commitment: witnesses
    commit_point_as_xy::<E, T>(&mut transcript, &proof.c1);

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
    commit_point_as_xy::<E, T>(&mut transcript, &proof.c2);
    // evaluation challenge
    let r = transcript.get_challenge();
    // commit evaluations
    let all_evaluations = flatten_all_evaluations(
        &setup_evaluations,
        &first_round_evaluations,
        &second_round_evaluations,
    );
    for eval in all_evaluations.iter() {
        transcript.commit_field_element(eval);
    }
    // opening linearization challenge
    let alpha = transcript.get_challenge();
    commit_point_as_xy::<E, T>(&mut transcript, &proof.w);

    // last opening challenge
    let y = transcript.get_challenge();

    // all system polynomials will be evaluated at z
    // then combined polynomials will be opened at h_i = r^power_i
    // then it becomes e.g C_i(x) = f_0(x^2) + x*f(x^2) in case of two polynomials
    assert_eq!(num_second_round_polys, 3);
    let power = lcm(&[
        num_setup_polys.next_power_of_two(),
        num_first_round_polys.next_power_of_two(),
        num_second_round_polys,
    ]);
    println!("LCM {power}");
    let z = r.pow(&[power as u64]);

    // all identity testing requires vanishing poly Z(X)
    // to be opened at z such that Q(z)=P(z)/V(z)
    let vanishing_at_z = evaluate_vanishing_for_size(&z, required_domain_size as u64);

    // Each degree k*n round polynomials broken into k sub polynomias each degree n
    // f(X) = f0(X) + X^n*f1(x) + .. x^{k-1}*f_{n-1}(X)
    // So LHS of each quotient needs to be constructed at z
    // f(z) = f0(z) + z^n*f1(z) + .. z^{k-1}*f_{n-1}(z)
    // So that compute powers of z = {1, z^n, z^2n .. z^{(k-1)*n}} at once
    // and later on reuse it in each quotient identity testing f(z) = (..) / V(z)
    // and test f(z)*V(z)- (...) = 0
    // Start to calculate LHS of main gate quotients in the same way explained above
    let mut main_gate_lhs: E::Fr = trace_and_gate_evaluations.main_gate_quotient_at_z;
    main_gate_lhs.mul_assign(&vanishing_at_z);

    let mut all_gates_iter = sorted_gates.clone().into_iter();
    let main_gate_internal = all_gates_iter.next().unwrap();
    assert!(&C::MainGate::default().into_internal() == &main_gate_internal);
    let main_gate = C::MainGate::default();

    let mut public_inputs_at_z = E::Fr::zero();
    for (idx, inp) in public_inputs.iter().enumerate() {
        let mut tmp = evaluate_lagrange_poly_at_point(idx, &domain, z)?;
        tmp.mul_assign(&inp);
        public_inputs_at_z.add_assign(&tmp);
    }

    let main_gate_rhs = compute_quotient_of_main_gate_at_z(
        main_gate.name(),
        &trace_and_gate_evaluations,
        &setup_evaluations,
        public_inputs_at_z,
    );

    // check main gate identity
    assert_eq!(main_gate_lhs, main_gate_rhs);

    // Do the same for custom gate
    if has_custom_gate {
        let mut custom_gate_lhs = trace_and_gate_evaluations
            .custom_gate_quotient_at_z
            .expect("custom gate quotients at z");
        custom_gate_lhs.mul_assign(&vanishing_at_z);

        let custom_gate = all_gates_iter.next().expect("custom gate");
        let custom_gate_rhs = compute_quotient_of_custom_gate_at_z(
            custom_gate.name(),
            &trace_and_gate_evaluations,
            &setup_evaluations,
        );
        assert_eq!(custom_gate_lhs, custom_gate_rhs);
    }

    let l_0_at_z = evaluate_l0_at_point(required_domain_size as u64, z)?;

    let SecondRoundEvaluations {
        copy_permutation_evaluations,
        lookup_evaluations,
    } = second_round_evaluations;

    // copy permutation identities

    // we have only 2 main gate types where both has the same number of variables
    // z(X)(A + beta*X + gamma)(B + beta*k_1*X + gamma)(C + beta*K_2*X + gamma)(D + beta*K_3*X + gamma) -
    // - (A + beta*perm_a(X) + gamma)(B + beta*perm_b(X) + gamma)(C + beta*perm_c(X) + gamma)*(D + beta*perm_d(X) + gamma)*Z(X*Omega)== 0
    let mut copy_permutation_first_quotient_lhs = copy_permutation_evaluations.first_quotient_at_z;
    copy_permutation_first_quotient_lhs.mul_assign(&vanishing_at_z);
    let mut copy_permutation_first_quotient_rhs_num_part = z;
    copy_permutation_first_quotient_rhs_num_part.mul_assign(&beta_for_copy_permutation);
    copy_permutation_first_quotient_rhs_num_part.add_assign(&gamma_for_copy_permutation);
    copy_permutation_first_quotient_rhs_num_part
        .add_assign(&trace_and_gate_evaluations.trace_evaluations_at_z[0]);

    assert_eq!(non_residues.len() + 1, num_state_polys);
    for (non_residue, state_poly) in non_residues.iter().zip(
        trace_and_gate_evaluations
            .trace_evaluations_at_z
            .iter()
            .take(num_state_polys)
            .skip(1),
    ) {
        let mut tmp = z;
        tmp.mul_assign(&non_residue);
        tmp.mul_assign(&beta_for_copy_permutation);
        tmp.add_assign(&gamma_for_copy_permutation);
        tmp.add_assign(state_poly);
        copy_permutation_first_quotient_rhs_num_part.mul_assign(&tmp);
    }
    copy_permutation_first_quotient_rhs_num_part
        .mul_assign(&copy_permutation_evaluations.grand_product_at_z);

    let mut copy_permutation_first_quotient_rhs_denum_part =
        copy_permutation_evaluations.grand_product_at_z_omega;
    assert_eq!(
        setup_evaluations.permutations_at_z.len(),
        trace_and_gate_evaluations.trace_evaluations_at_z.len()
    );
    for (permutation, state_poly) in setup_evaluations.permutations_at_z.iter().zip(
        trace_and_gate_evaluations
            .trace_evaluations_at_z
            .iter()
            .take(num_state_polys),
    ) {
        let mut tmp = beta_for_copy_permutation;
        tmp.mul_assign(&permutation);
        tmp.add_assign(&gamma_for_copy_permutation);
        tmp.add_assign(state_poly);
        copy_permutation_first_quotient_rhs_denum_part.mul_assign(&tmp);
    }

    let mut copy_permutation_first_quotient_rhs = copy_permutation_first_quotient_rhs_num_part;
    copy_permutation_first_quotient_rhs.sub_assign(&copy_permutation_first_quotient_rhs_denum_part);
    assert_eq!(
        copy_permutation_first_quotient_lhs,
        copy_permutation_first_quotient_rhs
    );

    // (Z(x) - 1) * L_{0} == 0
    let mut copy_permutation_second_quotient_rhs = copy_permutation_evaluations.grand_product_at_z;
    copy_permutation_second_quotient_rhs.sub_assign(&E::Fr::one());
    copy_permutation_second_quotient_rhs.mul_assign(&l_0_at_z);
    let mut second_quotient_at_z = copy_permutation_evaluations.second_quotient_at_z;
    second_quotient_at_z.mul_assign(&vanishing_at_z);
    assert_eq!(second_quotient_at_z, copy_permutation_second_quotient_rhs);

    if has_lookup {
        let eta_for_lookup = eta_for_lookup.expect("eta for lookup");
        let beta_for_lookup = beta_for_lookup.expect("beta for lookup");
        let gamma_for_lookup = gamma_for_lookup.expect("gamma for lookup");
        let mut beta_gamma = beta_for_lookup;
        beta_gamma.add_assign(&E::Fr::one());
        beta_gamma.mul_assign(&gamma_for_lookup);
        // lookup identities
        // ( Z(x*omega)*(\gamma*(1 + \beta) + s(x) + \beta * s(x*omega))) -
        // - Z(x) * (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega)) )*(X - omega^{n-1})
        let LookupEvaluations {
            s_poly_at_z,
            s_poly_at_z_omega,
            grand_product_at_z,
            grand_product_at_z_omega,
            first_quotient_at_z,
            second_quotient_at_z,
            third_quotient_at_z,
            ..
        } = lookup_evaluations.as_ref().expect("lookup evaluations");
        // LHS
        // f(z) = f0(z) + z*f1(z)
        let mut lookup_first_quotient_lhs = *first_quotient_at_z;
        lookup_first_quotient_lhs.mul_assign(&vanishing_at_z);

        // RHS
        let mut lookup_first_quotient_rhs_denum_part = *s_poly_at_z_omega;
        lookup_first_quotient_rhs_denum_part.mul_assign(&beta_for_lookup);
        lookup_first_quotient_rhs_denum_part.add_assign(s_poly_at_z);
        lookup_first_quotient_rhs_denum_part.add_assign(&beta_gamma);
        lookup_first_quotient_rhs_denum_part.mul_assign(grand_product_at_z_omega);

        // Prover doesn't open aggregated columns of table rather it opens each of them
        // seperately because they are committed in the first round
        // and there is no reandomness.

        // aggregate witnesses a + eta*b + eta^2*c + eta^3*table_type
        // expands into (((table_type*eta + c)*eta  + b)*eta + a)
        let lookup_table_type = setup_evaluations
            .lookup_table_type_at_z
            .as_ref()
            .expect("lookup table tyoe")
            .clone();
        let mut aggregated_lookup_f_at_z = lookup_table_type;
        for col in trace_and_gate_evaluations
            .trace_evaluations_at_z
            .iter()
            .take(num_state_polys - 1)
            .rev()
        {
            aggregated_lookup_f_at_z.mul_assign(&eta_for_lookup);
            aggregated_lookup_f_at_z.add_assign(col);
        }
        aggregated_lookup_f_at_z.mul_assign(
            &setup_evaluations
                .lookup_selector_at_z
                .expect("lookup selector"),
        );
        // col0 + eta * col1 + eta^2*col2 + eta^3*table_type
        let table_evals_at_z = setup_evaluations
            .lookup_tables_at_z
            .as_ref()
            .expect("table col evaluations at z");
        let table_evals_at_z_omega = setup_evaluations
            .lookup_tables_at_z_omega
            .as_ref()
            .expect("table col evaluations at z*w");
        let mut aggregated_lookup_table_cols_at_z = table_evals_at_z.last().unwrap().clone();
        let mut aggregated_lookup_table_cols_at_z_omega =
            table_evals_at_z_omega.last().unwrap().clone();
        for (at_z, at_z_omega) in table_evals_at_z.iter().take(num_state_polys - 1).rev().zip(
            table_evals_at_z_omega
                .iter()
                .take(num_state_polys - 1)
                .rev(),
        ) {
            aggregated_lookup_table_cols_at_z.mul_assign(&eta_for_lookup);
            aggregated_lookup_table_cols_at_z.add_assign(at_z);

            aggregated_lookup_table_cols_at_z_omega.mul_assign(&eta_for_lookup);
            aggregated_lookup_table_cols_at_z_omega.add_assign(at_z_omega);
        }
        // We also need to aggregate shifted table columns to construct t(z*w)
        // First identity is for multiset-equality
        let mut lookup_first_quotient_rhs_num_part = aggregated_lookup_table_cols_at_z_omega;
        lookup_first_quotient_rhs_num_part.mul_assign(&beta_for_lookup);
        lookup_first_quotient_rhs_num_part.add_assign(&aggregated_lookup_table_cols_at_z);
        lookup_first_quotient_rhs_num_part.add_assign(&beta_gamma);
        aggregated_lookup_f_at_z.add_assign(&gamma_for_lookup);
        lookup_first_quotient_rhs_num_part.mul_assign(&aggregated_lookup_f_at_z);
        let mut beta_one = beta_for_lookup;
        beta_one.add_assign(&E::Fr::one());
        lookup_first_quotient_rhs_num_part.mul_assign(&beta_one);
        lookup_first_quotient_rhs_num_part.mul_assign(&grand_product_at_z);

        let mut lookup_first_quotient_rhs = lookup_first_quotient_rhs_denum_part;
        lookup_first_quotient_rhs.sub_assign(&lookup_first_quotient_rhs_num_part);

        let last_omega = domain.generator.pow(&[required_domain_size as u64 - 1]);
        let mut tmp = z;
        tmp.sub_assign(&last_omega);
        lookup_first_quotient_rhs.mul_assign(&tmp);
        assert_eq!(lookup_first_quotient_lhs, lookup_first_quotient_rhs);

        // Then verify that first element of the grand product poly equals to 1
        // (Z(x) - 1) * L_{0} == 0
        let mut lookup_second_quotient_rhs = *grand_product_at_z;
        lookup_second_quotient_rhs.sub_assign(&E::Fr::one());
        lookup_second_quotient_rhs.mul_assign(&l_0_at_z);
        let mut second_quotient_at_z = *second_quotient_at_z;
        second_quotient_at_z.mul_assign(&vanishing_at_z);
        assert_eq!(second_quotient_at_z, lookup_second_quotient_rhs);

        // Also verify that last element is equals to expected value
        // (Z(x) - expected) * L_{n-1} == 0
        let expected = beta_gamma.pow([(required_domain_size - 1) as u64]);
        let l_last = evaluate_lagrange_poly_at_point(required_domain_size - 1, &domain, z)?;
        let mut lookup_second_quotient_rhs = *grand_product_at_z;
        lookup_second_quotient_rhs.sub_assign(&expected);
        lookup_second_quotient_rhs.mul_assign(&l_last);
        let mut third_quotient_at_z = *third_quotient_at_z;
        third_quotient_at_z.mul_assign(&vanishing_at_z);
        assert_eq!(third_quotient_at_z, lookup_second_quotient_rhs);
    }

    let FflonkProof {
        c1,
        c2,
        w,
        w_prime,
        lagrange_basis_inverses,
        ..
    } = proof;
    aggregate_points_and_check_pairing(
        vk,
        r,
        z,
        alpha,
        y,
        num_setup_polys,
        num_first_round_polys,
        num_second_round_polys,
        &all_evaluations,
        lagrange_basis_inverses,
        c1.clone(),
        c2.clone(),
        w.clone(),
        w_prime.clone(),
        setup_evaluations.requires_opening_at_shifted_point(),
        first_round_evaluations.requires_opening_at_shifted_point(),
    )
}

pub fn verify_flattened_proof<E: Engine, C: Circuit<E>, T: Transcript<E::Fr>>(
    vk: &FflonkVerificationKey<E, C>,
    proof: &FlattenedFflonkProof<E, C>,
    transcript_params: Option<T::InitializationParameters>,
) -> Result<bool, SynthesisError> {
    let mut transcript = if let Some(params) = transcript_params {
        T::new_from_params(params)
    } else {
        T::new()
    };

    let sorted_gates = sorted_gates_from_circuit_definitions::<_, C>();
    assert!(sorted_gates.len() > 0);

    assert!((vk.n + 1).is_power_of_two());
    let required_domain_size = vk.n.next_power_of_two();

    let domain = Domain::<E::Fr>::new_for_size(required_domain_size as u64)?;

    let FlattenedFflonkProof {
        evaluations,
        commitments,
        inputs,
        ..
    } = proof;

    let n = vk.n;
    let required_domain_size = n + 1;
    assert!(required_domain_size.is_power_of_two());
    assert!(required_domain_size.trailing_zeros() <= 23);

    let (num_setup_polys, num_first_round_polys, num_second_round_polys, _) =
        num_system_polys_from_vk(vk);
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
    assert!(proof.inputs.is_empty() == false);
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
    let power = lcm(&[
        num_setup_polys.next_power_of_two(),
        num_first_round_polys.next_power_of_two(),
        num_second_round_polys,
    ]);
    println!("LCM {power}");
    let z = r.pow(&[power as u64]);

    // all identity testing requires vanishing poly Z(X)
    // to be opened at z such that Q(z)=P(z)/V(z)
    let vanishing_at_z = evaluate_vanishing_for_size(&z, required_domain_size as u64);

    let offsets = EvaluationOffsets::from_vk(vk);

    // Each degree k*n round polynomials broken into k sub polynomias each degree n
    // f(X) = f0(X) + X^n*f1(x) + .. x^{k-1}*f_{n-1}(X)
    // So LHS of each quotient needs to be constructed at z
    // f(z) = f0(z) + z^n*f1(z) + .. z^{k-1}*f_{n-1}(z)
    // So that compute powers of z = {1, z^n, z^2n .. z^{(k-1)*n}} at once
    // and later on reuse it in each quotient identity testing f(z) = (..) / V(z)
    // and test f(z)*V(z)- (...) = 0

    // Start to calculate LHS of main gate quotients in the same way explained above
    let mut main_gate_lhs = evaluations[offsets.trace.main_gate_quotient_at_z];
    main_gate_lhs.mul_assign(&vanishing_at_z);

    let mut all_gates_iter = sorted_gates.clone().into_iter();
    let main_gate_internal = all_gates_iter.next().unwrap();
    assert!(&C::MainGate::default().into_internal() == &main_gate_internal);
    let main_gate = C::MainGate::default();

    let mut public_inputs_at_z = E::Fr::zero();
    assert!(inputs.is_empty() == false);
    for (idx, inp) in inputs.iter().enumerate() {
        let mut tmp = evaluate_lagrange_poly_at_point(idx, &domain, z)?;
        tmp.mul_assign(&inp);
        public_inputs_at_z.add_assign(&tmp);
    }

    let main_gate_rhs = compute_quotient_of_main_gate_at_z_flattened(
        main_gate.name(),
        &evaluations,
        public_inputs_at_z,
        &offsets,
    );
    // check main gate identity
    assert_eq!(main_gate_lhs, main_gate_rhs);

    // Do the same for custom gate
    if has_custom_gate {
        let mut custom_gate_lhs = evaluations[offsets.trace.custom_gate_quotient_at_z];
        custom_gate_lhs.mul_assign(&vanishing_at_z);
        let custom_gate = all_gates_iter.next().expect("custom gate");
        let custom_gate_rhs = compute_quotient_of_custom_gate_at_z_flattened(
            custom_gate.name(),
            &evaluations,
            &offsets,
        );

        assert_eq!(custom_gate_lhs, custom_gate_rhs);
    }

    let l_0_at_z = evaluate_l0_at_point(required_domain_size as u64, z)?;

    // copy permutation identities

    // we have only 2 main gate types where both has the same number of variables
    // z(X)(A + beta*X + gamma)(B + beta*k_1*X + gamma)(C + beta*K_2*X + gamma)(D + beta*K_3*X + gamma) -
    // - (A + beta*perm_a(X) + gamma)(B + beta*perm_b(X) + gamma)(C + beta*perm_c(X) + gamma)*(D + beta*perm_d(X) + gamma)*Z(X*Omega)== 0
    let mut copy_permutation_first_quotient_lhs =
        evaluations[offsets.copy_permutation.first_quotient_at_z];
    copy_permutation_first_quotient_lhs.mul_assign(&vanishing_at_z);

    let mut copy_permutation_first_quotient_rhs_num_part = z;
    copy_permutation_first_quotient_rhs_num_part.mul_assign(&beta_for_copy_permutation);
    copy_permutation_first_quotient_rhs_num_part.add_assign(&gamma_for_copy_permutation);
    copy_permutation_first_quotient_rhs_num_part
        .add_assign(&evaluations[offsets.trace.trace_evaluations_at_z]);

    assert_eq!(non_residues.len() + 1, num_state_polys);
    for (non_residue, state_poly) in non_residues.iter().zip(
        evaluations[offsets.trace.trace_evaluations_at_z
            ..offsets.trace.trace_evaluations_at_z + num_state_polys]
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
    copy_permutation_first_quotient_rhs_num_part
        .mul_assign(&evaluations[offsets.copy_permutation.grand_product_at_z]);

    let mut copy_permutation_first_quotient_rhs_denum_part =
        evaluations[offsets.copy_permutation.grand_product_at_z_omega];
    for (permutation, state_poly) in evaluations
        [offsets.setup.permutations_at_z..offsets.setup.permutations_at_z + num_state_polys]
        .iter()
        .zip(
            evaluations[offsets.trace.trace_evaluations_at_z
                ..offsets.trace.trace_evaluations_at_z + num_state_polys]
                .iter(),
        )
    {
        let mut tmp = beta_for_copy_permutation;
        tmp.mul_assign(&permutation);
        tmp.add_assign(&gamma_for_copy_permutation);
        tmp.add_assign(state_poly);
        copy_permutation_first_quotient_rhs_denum_part.mul_assign(&tmp);
    }

    let mut copy_permutation_first_quotient_rhs = copy_permutation_first_quotient_rhs_num_part;
    copy_permutation_first_quotient_rhs.sub_assign(&copy_permutation_first_quotient_rhs_denum_part);
    assert_eq!(
        copy_permutation_first_quotient_lhs,
        copy_permutation_first_quotient_rhs
    );

    // (Z(x) - 1) * L_{0} == 0
    let mut copy_permutation_second_quotient_rhs =
        evaluations[offsets.copy_permutation.grand_product_at_z];
    copy_permutation_second_quotient_rhs.sub_assign(&E::Fr::one());
    copy_permutation_second_quotient_rhs.mul_assign(&l_0_at_z);
    let mut second_quotient_at_z = evaluations[offsets.copy_permutation.second_quotient_at_z];
    second_quotient_at_z.mul_assign(&vanishing_at_z);
    assert_eq!(second_quotient_at_z, copy_permutation_second_quotient_rhs);

    if has_lookup {
        let eta_for_lookup = eta_for_lookup.expect("eta for lookup");
        let beta_for_lookup = beta_for_lookup.expect("beta for lookup");
        let gamma_for_lookup = gamma_for_lookup.expect("gamma for lookup");
        let mut beta_gamma = beta_for_lookup;
        beta_gamma.add_assign(&E::Fr::one());
        beta_gamma.mul_assign(&gamma_for_lookup);
        // lookup identities
        // ( Z(x*omega)*(\gamma*(1 + \beta) + s(x) + \beta * s(x*omega))) -
        // - Z(x) * (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega)) )*(X - omega^{n-1})
        // LHS
        // f(z) = f0(z) + z*f1(z)
        let lookup_offsets = offsets.lookup.expect("lookup offsets");
        let mut lookup_first_quotient_lhs = evaluations[lookup_offsets.first_quotient_at_z];
        lookup_first_quotient_lhs.mul_assign(&vanishing_at_z);

        // RHS
        let mut lookup_first_quotient_rhs_denum_part =
            evaluations[lookup_offsets.s_poly_at_z_omega];
        lookup_first_quotient_rhs_denum_part.mul_assign(&beta_for_lookup);
        lookup_first_quotient_rhs_denum_part.add_assign(&evaluations[lookup_offsets.s_poly_at_z]);
        lookup_first_quotient_rhs_denum_part.add_assign(&beta_gamma);
        lookup_first_quotient_rhs_denum_part
            .mul_assign(&evaluations[lookup_offsets.grand_product_at_z_omega]);
        // Prover doesn't open aggregated columns of table rather it opens each of them
        // seperately because they are committed in the first round
        // and there is no reandomness.

        // aggregate witnesses a + eta*b + eta^2*c + eta^3*table_type
        // expands into (((table_type*eta + c)*eta  + b)*eta + a)
        let mut aggregated_lookup_f_at_z = evaluations[offsets.setup.lookup_table_type_at_z];
        for col in evaluations[offsets.trace.trace_evaluations_at_z
            ..offsets.trace.trace_evaluations_at_z + num_state_polys]
            .iter()
            .take(num_state_polys - 1)
            .rev()
        {
            aggregated_lookup_f_at_z.mul_assign(&eta_for_lookup);
            aggregated_lookup_f_at_z.add_assign(col);
        }
        aggregated_lookup_f_at_z.mul_assign(&evaluations[offsets.setup.lookup_selector_at_z]);
        // col0 + eta * col1 + eta^2*col2 + eta^3*table_type
        let mut aggregated_lookup_table_cols_at_z =
            evaluations[offsets.setup.lookup_tables_at_z + 3];
        let mut aggregated_lookup_table_cols_at_z_omega =
            evaluations[offsets.setup.lookup_tables_at_z_omega + 3];
        for (at_z, at_z_omega) in evaluations
            [offsets.setup.lookup_tables_at_z..offsets.setup.lookup_tables_at_z + 3]
            .iter()
            .take(num_state_polys - 1)
            .rev()
            .zip(
                evaluations[offsets.setup.lookup_tables_at_z_omega
                    ..offsets.setup.lookup_tables_at_z_omega + 3]
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
        beta_one.add_assign(&E::Fr::one());
        lookup_first_quotient_rhs_num_part.mul_assign(&beta_one);
        lookup_first_quotient_rhs_num_part
            .mul_assign(&evaluations[lookup_offsets.grand_product_at_z]);

        let mut lookup_first_quotient_rhs = lookup_first_quotient_rhs_denum_part;
        lookup_first_quotient_rhs.sub_assign(&lookup_first_quotient_rhs_num_part);

        let last_omega = domain.generator.pow(&[required_domain_size as u64 - 1]);
        let mut tmp = z;
        tmp.sub_assign(&last_omega);
        lookup_first_quotient_rhs.mul_assign(&tmp);
        assert_eq!(lookup_first_quotient_lhs, lookup_first_quotient_rhs);

        // Then verify that first element of the grand product poly equals to 1
        // (Z(x) - 1) * L_{0} == 0
        let mut lookup_second_quotient_rhs = evaluations[lookup_offsets.grand_product_at_z];
        lookup_second_quotient_rhs.sub_assign(&E::Fr::one());
        lookup_second_quotient_rhs.mul_assign(&l_0_at_z);
        let mut second_quotient_at_z = evaluations[lookup_offsets.second_quotient_at_z];
        second_quotient_at_z.mul_assign(&vanishing_at_z);
        assert_eq!(second_quotient_at_z, lookup_second_quotient_rhs);

        // Also verify that last element is equals to expected value
        // (Z(x) - expected) * L_{n-1} == 0
        let expected = beta_gamma.pow([(required_domain_size - 1) as u64]);
        let l_last = evaluate_lagrange_poly_at_point(required_domain_size - 1, &domain, z)?;
        let mut lookup_second_quotient_rhs = evaluations[lookup_offsets.grand_product_at_z];
        lookup_second_quotient_rhs.sub_assign(&expected);

        lookup_second_quotient_rhs.mul_assign(&l_last);
        let mut third_quotient_at_z = evaluations[lookup_offsets.third_quotient_at_z];
        third_quotient_at_z.mul_assign(&vanishing_at_z);
        assert_eq!(third_quotient_at_z, lookup_second_quotient_rhs);
    }

    // Now it is time to combine all rounds in a combined poly
    // C0(x) = f0(X^k0) + X*f1(X^k0) + X*f_{k0-1}(X^k0) where k0 is the total number of

    // Since fft-style combination requires roots of unity, we need total number of polyn
    // be power of two, in our case none of them are power of two so that we apply a padding
    // with zeroes here.
    // Compute h0 = r^(power/k0) here, it will be plugged in the C0(X)
    let setup_requires_opening_at_shifted_point = has_lookup; // TODO
    let first_round_requires_opening_at_shifted_point =
        requires_trace_polys_opening_at_shifted_point(main_gate_internal);
    aggregate_points_and_check_pairing(
        vk,
        r,
        z,
        alpha,
        y,
        num_setup_polys,
        num_first_round_polys,
        num_second_round_polys,
        &proof.evaluations,
        &proof.lagrange_basis_inverses,
        c1,
        c2,
        w,
        w_prime,
        setup_requires_opening_at_shifted_point,
        first_round_requires_opening_at_shifted_point,
    )
}

pub fn aggregate_points_and_check_pairing<E: Engine, C: Circuit<E>>(
    vk: &FflonkVerificationKey<E, C>,
    r: E::Fr,
    z: E::Fr,
    alpha: E::Fr,
    y: E::Fr,
    num_setup_polys: usize,
    num_first_round_polys: usize,
    num_second_round_polys: usize,
    all_evaluations: &[E::Fr],
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

    let power = lcm(&[
        interpolation_size_of_setup,
        interpolation_size_of_first_round,
        interpolation_size_of_second_round,
    ]);
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
    

    let [setup_r_at_y, mut first_round_r_at_y, mut second_round_r_at_y] =
        evaluate_r_polys_at_point_with_flattened_evals_and_precomputed_basis(
            &all_evaluations,
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
    let sparse_polys_for_first_round_at_y =
        evaluate_multiple_sparse_polys(sparse_polys_for_first_round, y);
    let sparse_polys_for_second_round_at_y =
        evaluate_multiple_sparse_polys(sparse_polys_for_second_round, y);

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
    assert!(valid, "pairing check failed");

    Ok(valid)
}