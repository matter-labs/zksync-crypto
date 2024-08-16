use super::*;

pub fn create_proof<
    E: Engine,
    C: Circuit<E>,
    P: PlonkConstraintSystemParams<E>,
    MG: MainGate<E>,
    S: SynthesisMode,
    T: Transcript<<E as ScalarEngine>::Fr>,
>(
    assembly: Assembly<E, P, MG, S>,
    worker: &Worker,
    setup: &FflonkSetup<E, C>,
    mon_crs: &Crs<E, CrsForMonomialForm>,
    transcript_params: Option<T::InitializationParameters>,
) -> Result<FflonkProof<E, C>, SynthesisError> {
    assert!(<SynthesisModeTesting as SynthesisMode>::PRODUCE_WITNESS);
    assert!(assembly.is_finalized);

    let mut transcript = if let Some(params) = transcript_params {
        T::new_from_params(params)
    } else {
        T::new()
    };

    let input_values = assembly.input_assingments.clone();
    assert!(input_values.is_empty() == false);
    println!("num public inputs {}", input_values.len());
    for inp in input_values.iter() {
        transcript.commit_field_element(inp);
    }
    commit_point_as_xy::<E, _>(&mut transcript, &setup.c0_commitment);

    let required_domain_size = assembly.n() + 1;
    assert!(required_domain_size.is_power_of_two());
    assert!(required_domain_size.trailing_zeros() <= 23, "Only trace length up to 2^23 is allowed");

    let has_lookup = assembly.tables.len() > 0;
    assert_eq!(has_lookup, false, "lookup support is disabled");
    if has_lookup {
        assert!(assembly.num_table_lookups > 0);
    }
    assert_eq!(
        &sorted_gates_from_circuit_definitions::<_, C>(),
        &assembly.sorted_gates
    );
    let main_gate = assembly.sorted_gates[0].clone();
    assert_eq!(main_gate.variable_polynomials().len(), P::STATE_WIDTH);

    let has_custom_gate = assembly.sorted_gates.len() > 1;
    let custom_gate = if has_custom_gate {
        assert!(P::HAS_CUSTOM_GATES);
        Some(assembly.sorted_gates[1].clone())
    } else {
        None
    };
    assert_eq!(has_custom_gate, false, "custom gate support is disabled");

    let num_state_polys = P::STATE_WIDTH;
    assert_eq!(num_state_polys, main_gate.variable_polynomials().len());
    let num_witness_polys = P::WITNESS_WIDTH;
    let main_gate = assembly.sorted_gates[0].clone();
    assert_eq!(num_state_polys, main_gate.variable_polynomials().len());
    let (num_setup_polys, num_first_round_polys, num_second_round_polys, _) =
        num_system_polys_from_assembly::<_, _, _, _, C>(&assembly);

    let max_combined_degree = compute_max_combined_degree_from_assembly::<_, _, _, _, C>(&assembly);
    assert!(max_combined_degree <= mon_crs.g1_bases.len());
    dbg!(num_setup_polys);
    dbg!(num_first_round_polys);
    dbg!(num_second_round_polys);
    let coset_factor = E::Fr::multiplicative_generator();

    let mut lde_factor = num_state_polys.next_power_of_two();
    for g in assembly.sorted_gates.iter() {
        let degree = g.degree();
        if degree > lde_factor {
            lde_factor = degree;
        }
    }
    assert!(lde_factor <= 4);
    assert_eq!(
        assembly.max_constraint_degree.next_power_of_two(),
        lde_factor
    );
    let domain: Domain<E::Fr> = Domain::new_for_size(required_domain_size as u64)?;
    let lde_domain_size = required_domain_size * lde_factor;
    let lde_domain: Domain<E::Fr> = Domain::new_for_size(lde_domain_size as u64)?;

    let mut values_storage = assembly.make_assembled_poly_storage(worker, true)?;

    let omegas_bitreversed = BitReversedOmegas::<E::Fr>::new_for_domain_size(required_domain_size);
    let omegas_inv_bitreversed =
        <OmegasInvBitreversed<E::Fr> as CTPrecomputations<E::Fr>>::new_for_domain_size(
            required_domain_size,
        );
    let setup = &setup.original_setup;
    if S::PRODUCE_SETUP {
        let permutation_polys = assembly.make_permutations(&worker)?;
        assert_eq!(permutation_polys.len(), num_state_polys);

        for (idx, poly) in permutation_polys.into_iter().enumerate() {
            let key = PolyIdentifier::PermutationPolynomial(idx);
            let poly = PolynomialProxy::from_owned(poly);
            values_storage.setup_map.insert(key, poly);
        }
    } else {
        // compute from setup
        for idx in 0..num_state_polys {
            let key = PolyIdentifier::PermutationPolynomial(idx);
            // let vals = setup.permutation_monomials[idx].clone().fft(&worker).into_coeffs();
            let vals = setup.permutation_monomials[idx]
                .clone()
                .fft_using_bitreversed_ntt(&worker, &omegas_bitreversed, &E::Fr::one())?
                .into_coeffs();
            let poly = Polynomial::from_values_unpadded(vals)?;
            let poly = PolynomialProxy::from_owned(poly);
            values_storage.setup_map.insert(key, poly);
        }
    }

    let mut ldes_storage = AssembledPolynomialStorage::<E>::new(
        true,
        assembly.max_constraint_degree.next_power_of_two(),
    );

    let mut monomials_storage = Assembly::<E, P, MG, S>::create_monomial_storage(
        &worker,
        &omegas_inv_bitreversed,
        &values_storage,
        true,
    )?;

    monomials_storage.extend_from_setup(setup)?;

    let inverse_divisor_on_coset_lde_natural_ordering = {
        let mut vanishing_poly_inverse_bitreversed =
            evaluate_vanishing_polynomial_of_degree_on_domain_size::<E::Fr>(
                required_domain_size as u64,
                &E::Fr::multiplicative_generator(),
                (required_domain_size * lde_factor) as u64,
                &worker,
            )?;
        vanishing_poly_inverse_bitreversed.batch_inversion(&worker)?;
        // vanishing_poly_inverse_bitreversed.bitreverse_enumeration(&worker)?;

        vanishing_poly_inverse_bitreversed
    };

    let mut trace_monomials = vec![];
    for i in 0..num_state_polys {
        let key = PolyIdentifier::VariablesPolynomial(i);
        let poly_ref = monomials_storage.get_poly(key);
        trace_monomials.push(poly_ref);
    }

    for i in 0..num_witness_polys {
        let key = PolyIdentifier::WitnessPolynomial(i);
        let poly_ref = monomials_storage.get_poly(key);
        trace_monomials.push(poly_ref);
    }

    let (main_gate_quotient_monomial, custom_gate_quotient_monomial) =
        compute_gate_quotients::<_, MG>(
            assembly.sorted_gates.clone(),
            &input_values,
            required_domain_size,
            &mut ldes_storage,
            &monomials_storage,
            &omegas_bitreversed,
            &omegas_inv_bitreversed,
            &inverse_divisor_on_coset_lde_natural_ordering,
            E::Fr::multiplicative_generator(),
            worker,
        )?;
    let main_gate = assembly.sorted_gates[0].clone();
    let trace_and_gate_monomials = TraceAndGateMonomials {
        trace_monomials,
        main_gate_quotient_monomial,
        custom_gate_quotient_monomial,
        num_state_polys,
        num_witness_polys,
        main_gate: main_gate.clone(),
        custom_gate,
        num_polys: num_first_round_polys,
    };

    let first_round_monomials = FirstRoundMonomials {
        trace_and_gate_monomials,
    };

    let flattened_first_round_monomials = first_round_monomials.flatten();
    assert_eq!(flattened_first_round_monomials.len(), num_first_round_polys);
    let c1 =
        combine_mixed_degree_monomials(&flattened_first_round_monomials, required_domain_size)?;
    let expected_degree = num_first_round_polys
        * main_gate_quotient_degree(&assembly.sorted_gates)
        * required_domain_size + num_first_round_polys -1;
    assert_eq!(c1.size(), expected_degree);
    assert!(c1.size() <= max_combined_degree, "left:{}\nright:{}", c1.size(), max_combined_degree);
    let c1_commitment = commit_using_monomials(&c1, mon_crs, worker)?;
    commit_point_as_xy::<E, T>(&mut transcript, &c1_commitment);
    let mut x_poly_lde_bitreversed =
        Polynomial::from_values(vec![coset_factor; required_domain_size * lde_factor])?;
    x_poly_lde_bitreversed.distribute_powers(&worker, lde_domain.generator);
    x_poly_lde_bitreversed.bitreverse_enumeration(&worker);

    assert_eq!(x_poly_lde_bitreversed.size(), lde_domain_size);

    let l_0 =
        calculate_lagrange_poly::<E::Fr>(&worker, required_domain_size.next_power_of_two(), 0)?;

    let l_0_coset_lde_bitreversed = l_0.bitreversed_lde_using_bitreversed_ntt(
        &worker,
        lde_factor,
        &omegas_bitreversed,
        &coset_factor,
    )?;

    // step 2 - grand product arguments

    let beta_for_copy_permutation = transcript.get_challenge();
    let gamma_for_copy_permutation = transcript.get_challenge();

    // copy permutation grand product argument

    let mut grand_products_protos_with_gamma = vec![];

    for i in 0..num_state_polys {
        let id = PolyIdentifier::VariablesPolynomial(i);

        let mut p = values_storage.state_map.get(&id).unwrap().as_ref().clone();
        p.add_constant(&worker, &gamma_for_copy_permutation);

        grand_products_protos_with_gamma.push(p);
    }

    let mut domain_elements =
        materialize_domain_elements_with_natural_enumeration(&domain, &worker);

    domain_elements
        .pop()
        .expect("must pop last element for omega^i");

    let non_residues = make_non_residues::<E::Fr>(num_state_polys - 1);

    let mut domain_elements_poly_by_beta = Polynomial::from_values_unpadded(domain_elements)?;
    domain_elements_poly_by_beta.scale(&worker, beta_for_copy_permutation);

    // we take A, B, C, ... values and form (A + beta * X * non_residue + gamma), etc and calculate their grand product

    let mut z_num = {
        let mut grand_products_proto_it = grand_products_protos_with_gamma.iter().cloned();

        let mut z_1 = grand_products_proto_it.next().unwrap();
        z_1.add_assign(&worker, &domain_elements_poly_by_beta);

        for (mut p, non_res) in grand_products_proto_it.zip(non_residues.iter()) {
            p.add_assign_scaled(&worker, &domain_elements_poly_by_beta, non_res);
            z_1.mul_assign(&worker, &p);
        }

        z_1
    };

    // we take A, B, C, ... values and form (A + beta * perm_a + gamma), etc and calculate their grand product

    let mut permutation_polynomials_values_of_size_n_minus_one = vec![];

    for idx in 0..num_state_polys {
        let key = PolyIdentifier::PermutationPolynomial(idx);

        let mut coeffs = values_storage.get_poly(key).clone().into_coeffs();
        coeffs.pop().unwrap();

        let p = Polynomial::from_values_unpadded(coeffs)?;
        permutation_polynomials_values_of_size_n_minus_one.push(p);
    }

    let z_den = {
        assert_eq!(
            permutation_polynomials_values_of_size_n_minus_one.len(),
            grand_products_protos_with_gamma.len()
        );
        let mut grand_products_proto_it = grand_products_protos_with_gamma.into_iter();
        let mut permutation_polys_it = permutation_polynomials_values_of_size_n_minus_one.iter();

        let mut z_2 = grand_products_proto_it.next().unwrap();
        z_2.add_assign_scaled(
            &worker,
            permutation_polys_it.next().unwrap(),
            &beta_for_copy_permutation,
        );

        for (mut p, perm) in grand_products_proto_it.zip(permutation_polys_it) {
            // permutation polynomials
            p.add_assign_scaled(&worker, &perm, &beta_for_copy_permutation);
            z_2.mul_assign(&worker, &p);
        }

        z_2.batch_inversion(&worker)?;

        z_2
    };

    z_num.mul_assign(&worker, &z_den);
    drop(z_den);

    let z = z_num.calculate_shifted_grand_product(&worker)?;
    drop(z_num);

    assert!(z.size().is_power_of_two());

    assert!(z.as_ref()[0] == E::Fr::one());

    let copy_permutation_z_monomial =
        z.ifft_using_bitreversed_ntt(&worker, &omegas_inv_bitreversed, &E::Fr::one())?;

    // Prove copy-permutation identities

    // For both Z_1 and Z_2 we first check for grand products
    // z*(X)(A + beta*X + gamma)(B + beta*k_1*X + gamma)(C + beta*K_2*X + gamma) -
    // - (A + beta*perm_a(X) + gamma)(B + beta*perm_b(X) + gamma)(C + beta*perm_c(X) + gamma)*Z(X*Omega)== 0

    // Prepare lde of state polynomials

    let state_poly_ids: Vec<_> = (0..num_state_polys)
        .map(|idx| PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(idx)))
        .collect();

    for poly_idx in state_poly_ids.into_iter() {
        ensure_in_map_or_create(
            &worker,
            poly_idx,
            required_domain_size,
            &omegas_bitreversed,
            lde_factor,
            coset_factor,
            &monomials_storage,
            &mut ldes_storage,
        )?;
    }

    let z_coset_lde_bitreversed = copy_permutation_z_monomial
        .clone()
        .bitreversed_lde_using_bitreversed_ntt(
            &worker,
            lde_factor,
            &omegas_bitreversed,
            &coset_factor,
        )?;

    assert!(z_coset_lde_bitreversed.size() == lde_domain_size);

    let z_shifted_coset_lde_bitreversed =
        z_coset_lde_bitreversed.clone_shifted_assuming_bitreversed(lde_factor, &worker)?;

    assert!(z_shifted_coset_lde_bitreversed.size() == lde_domain_size);

    let mut copy_permutation_first_quotient_rhs_num = z_coset_lde_bitreversed.clone();
    let mut tmp = ldes_storage
        .state_map
        .get(&PolyIdentifier::VariablesPolynomial(0))
        .unwrap()
        .as_ref()
        .clone();
    tmp.add_constant(&worker, &gamma_for_copy_permutation);
    tmp.add_assign_scaled(&worker, &x_poly_lde_bitreversed, &beta_for_copy_permutation);
    copy_permutation_first_quotient_rhs_num.mul_assign(&worker, &tmp);

    assert_eq!(non_residues.len() + 1, num_state_polys);

    for (poly_idx, non_res) in (1..num_state_polys).zip(non_residues.iter()) {
        let mut factor = beta_for_copy_permutation;
        factor.mul_assign(&non_res);

        let key = PolyIdentifier::VariablesPolynomial(poly_idx);
        tmp.reuse_allocation(&ldes_storage.state_map.get(&key).unwrap().as_ref());
        tmp.add_constant(&worker, &gamma_for_copy_permutation);
        tmp.add_assign_scaled(&worker, &x_poly_lde_bitreversed, &factor);
        copy_permutation_first_quotient_rhs_num.mul_assign(&worker, &tmp);
    }

    let mut copy_permutation_first_quotient_rhs_denum = z_shifted_coset_lde_bitreversed;

    // A + beta*perm_a + gamma
    for idx in 0..num_state_polys {
        let key = PolyIdentifier::VariablesPolynomial(idx);

        tmp.reuse_allocation(&ldes_storage.state_map.get(&key).unwrap().as_ref());
        tmp.add_constant(&worker, &gamma_for_copy_permutation);

        let key = PolyIdentifier::PermutationPolynomial(idx);
        let perm = monomials_storage
            .get_poly(key)
            .clone()
            .bitreversed_lde_using_bitreversed_ntt(
                &worker,
                lde_factor,
                &omegas_bitreversed,
                &coset_factor,
            )?;
        tmp.add_assign_scaled(&worker, &perm, &beta_for_copy_permutation);
        copy_permutation_first_quotient_rhs_denum.mul_assign(&worker, &tmp);
        drop(perm);
    }
    copy_permutation_first_quotient_rhs_num
        .sub_assign(&worker, &copy_permutation_first_quotient_rhs_denum);
    let copy_perm_quotient_degree = num_state_polys;
    let first_copy_perm_quotient_monomial = compute_quotient_monomial(
        worker,
        copy_permutation_first_quotient_rhs_num,
        &inverse_divisor_on_coset_lde_natural_ordering,
        coset_factor,
        lde_factor,
        required_domain_size,
        copy_perm_quotient_degree,
    )?;
    assert_eq!(
        first_copy_perm_quotient_monomial.size(),
        copy_perm_quotient_degree * required_domain_size
    );
    drop(tmp);

    // (Z(x) - 1) * L_{0} == 0
    let mut z_minus_one_by_l_0 = z_coset_lde_bitreversed;
    z_minus_one_by_l_0.sub_constant(&worker, &E::Fr::one());
    z_minus_one_by_l_0.mul_assign(&worker, &l_0_coset_lde_bitreversed);
    let second_copy_perm_quotient_monomial = compute_quotient_monomial(
        worker,
        z_minus_one_by_l_0,
        &inverse_divisor_on_coset_lde_natural_ordering,
        coset_factor,
        lde_factor,
        required_domain_size,
        1,
    )?;
    assert_eq!(
        second_copy_perm_quotient_monomial.size(),
        required_domain_size
    );
    let copy_permutation_monomials = CopyPermutationMonomials::<E::Fr> {
        grand_product_monomial: copy_permutation_z_monomial,
        first_quotient: first_copy_perm_quotient_monomial,
        second_quotient: second_copy_perm_quotient_monomial,
    };

    // step 1.5 - if there are lookup tables then draw random "eta" to linearlize over tables
    let lookup_monomials = if has_lookup {
        let eta = transcript.get_challenge();
        let beta = transcript.get_challenge();
        let gamma = transcript.get_challenge();

        // these are selected rows from witness (where lookup applies)

        let (selector_poly, table_type_mononial, table_type_values) = if S::PRODUCE_SETUP {
            let selector_for_lookup_values = assembly.calculate_lookup_selector_values()?;
            assert!((selector_for_lookup_values.len() + 1).is_power_of_two());
            let table_type_values = assembly.calculate_table_type_values()?;

            assert_eq!(selector_for_lookup_values.len(), table_type_values.len());

            let table_type_poly_monomial = {
                let mon = Polynomial::from_values(table_type_values.clone())?;
                let mon = mon.ifft_using_bitreversed_ntt(
                    &worker,
                    &omegas_inv_bitreversed,
                    &E::Fr::one(),
                )?;

                mon
            };

            let selector_poly =
                Polynomial::<E::Fr, Values>::from_values(selector_for_lookup_values)?
                    .ifft_using_bitreversed_ntt(&worker, &omegas_inv_bitreversed, &E::Fr::one())?;

            let selector_poly = PolynomialProxy::from_owned(selector_poly);
            let table_type_poly = PolynomialProxy::from_owned(table_type_poly_monomial);

            (selector_poly, table_type_poly, table_type_values)
        } else {
            let selector_poly_ref = setup
                .lookup_selector_monomial
                .as_ref()
                .expect("setup must contain lookup selector poly");
            let selector_poly = PolynomialProxy::from_borrowed(selector_poly_ref);

            let table_type_poly_ref = setup
                .lookup_table_type_monomial
                .as_ref()
                .expect("setup must contain lookup table type poly");
            let table_type_poly = PolynomialProxy::from_borrowed(table_type_poly_ref);

            // let mut table_type_values = table_type_poly_ref.clone().fft(&worker).into_coeffs();
            let mut table_type_values = table_type_poly_ref
                .clone()
                .fft_using_bitreversed_ntt(&worker, &omegas_bitreversed, &E::Fr::one())?
                .into_coeffs();

            table_type_values.pop().unwrap();

            (selector_poly, table_type_poly, table_type_values)
        };

        assert!((table_type_values.len() + 1).is_power_of_two());
        let witness_len = required_domain_size - 1;
        assert!((witness_len + 1).is_power_of_two());
        assert_eq!(table_type_values.len(), witness_len);

        let f_poly_values_aggregated = {
            let mut table_contributions_values = if S::PRODUCE_SETUP && S::PRODUCE_WITNESS {
                let masked_entries_using_bookkept_bitmasks =
                    assembly.calculate_masked_lookup_entries(&values_storage)?;

                let typical_len = masked_entries_using_bookkept_bitmasks[0].len();
                assert!((typical_len + 1).is_power_of_two());

                masked_entries_using_bookkept_bitmasks
            } else {
                assert!(S::PRODUCE_WITNESS);
                // let selector_values = PolynomialProxy::from_owned(selector_poly.as_ref().clone().fft(&worker));
                let selector_values = selector_poly.as_ref().clone().fft_using_bitreversed_ntt(
                    &worker,
                    &omegas_bitreversed,
                    &E::Fr::one(),
                )?;

                let selector_values = PolynomialProxy::from_owned(selector_values);

                assembly.calculate_masked_lookup_entries_using_selector(
                    &values_storage,
                    &selector_values,
                )?
            };

            assert_eq!(table_type_values.len(), table_contributions_values[0].len());

            assert_eq!(table_contributions_values.len(), 3);

            assert_eq!(witness_len, table_contributions_values[0].len());

            let mut f_poly_values_aggregated = table_contributions_values
                .drain(0..1)
                .collect::<Vec<_>>()
                .pop()
                .unwrap();

            let mut current = eta;
            for t in table_contributions_values.into_iter() {
                let op = BinopAddAssignScaled::new(current);
                binop_over_slices(&worker, &op, &mut f_poly_values_aggregated, &t);

                current.mul_assign(&eta);
            }

            // add table type marker
            let op = BinopAddAssignScaled::new(current);
            binop_over_slices(
                &worker,
                &op,
                &mut f_poly_values_aggregated,
                &table_type_values,
            );

            Polynomial::from_values_unpadded(f_poly_values_aggregated)?
        };

        let (t_poly_values, t_poly_values_shifted, t_poly_monomial) = if S::PRODUCE_SETUP {
            // these are unsorted rows of lookup tables
            let mut t_poly_ends =
                assembly.calculate_t_polynomial_values_for_single_application_tables()?;
            assert_eq!(t_poly_ends.len(), 4);

            let mut t_poly_values_aggregated =
                t_poly_ends.drain(0..1).collect::<Vec<_>>().pop().unwrap();
            let mut current = eta;
            for t in t_poly_ends.into_iter() {
                let op = BinopAddAssignScaled::new(current);
                binop_over_slices(&worker, &op, &mut t_poly_values_aggregated, &t);

                current.mul_assign(&eta);
            }

            let copy_start = witness_len - t_poly_values_aggregated.len();
            let mut full_t_poly_values = vec![E::Fr::zero(); witness_len];
            let mut full_t_poly_values_shifted = full_t_poly_values.clone();

            full_t_poly_values[copy_start..].copy_from_slice(&t_poly_values_aggregated);
            full_t_poly_values_shifted[(copy_start - 1)..(witness_len - 1)]
                .copy_from_slice(&t_poly_values_aggregated);

            assert!(full_t_poly_values[0].is_zero());

            let t_poly_monomial = {
                let mon = Polynomial::from_values(full_t_poly_values.clone())?;
                let mon = mon.ifft_using_bitreversed_ntt(
                    &worker,
                    &omegas_inv_bitreversed,
                    &E::Fr::one(),
                )?;

                mon
            };

            (
                PolynomialProxy::from_owned(Polynomial::from_values_unpadded(full_t_poly_values)?),
                PolynomialProxy::from_owned(Polynomial::from_values_unpadded(
                    full_t_poly_values_shifted,
                )?),
                PolynomialProxy::from_owned(t_poly_monomial),
            )
        } else {
            let mut t_poly_values_monomial_aggregated = setup.lookup_tables_monomials[0].clone();
            let mut current = eta;
            for idx in 1..4 {
                let to_aggregate_ref = &setup.lookup_tables_monomials[idx];
                t_poly_values_monomial_aggregated.add_assign_scaled(
                    &worker,
                    to_aggregate_ref,
                    &current,
                );

                current.mul_assign(&eta);
            }

            assert!(t_poly_values_monomial_aggregated.size().is_power_of_two());

            let mut t_poly_values = t_poly_values_monomial_aggregated
                .clone()
                .fft_using_bitreversed_ntt(&worker, &omegas_bitreversed, &E::Fr::one())?;
            assert!(t_poly_values.as_ref().last().unwrap().is_zero());
            assert!(t_poly_values.size().is_power_of_two());

            // let mut t_values_shifted_coeffs = vec![E::Fr::zero(); t_poly_values.size()];
            // // manually shift by 1
            // t_values_shifted_coeffs[1..].copy_from_slice(&t_poly_values.as_ref()[0..(t_poly_values.size()-1)]);
            // t_values_shifted_coeffs[0] = t_poly_values.as_ref()[(t_poly_values.size()-1)];

            let mut t_values_shifted_coeffs = t_poly_values.clone().into_coeffs();
            let _last = t_poly_values.pop_last()?;
            assert!(_last.is_zero());
            let _: Vec<_> = t_values_shifted_coeffs.drain(0..1).collect();

            let t_poly_values_shifted = Polynomial::from_values_unpadded(t_values_shifted_coeffs)?;

            assert_eq!(witness_len, t_poly_values.size());
            assert_eq!(witness_len, t_poly_values_shifted.size());

            (
                PolynomialProxy::from_owned(t_poly_values),
                PolynomialProxy::from_owned(t_poly_values_shifted),
                PolynomialProxy::from_owned(t_poly_values_monomial_aggregated),
            )
        };

        let (s_poly_monomial, s_poly_unpadded_values, s_shifted_unpadded_values) = {
            let s_poly_values_aggregated =
                assembly.calculate_s_poly_contributions_from_witness(eta)?;

            let sorted_copy_start = witness_len - s_poly_values_aggregated.len();

            let mut full_s_poly_values = vec![E::Fr::zero(); witness_len];
            let mut full_s_poly_values_shifted = full_s_poly_values.clone();

            full_s_poly_values[sorted_copy_start..].copy_from_slice(&s_poly_values_aggregated);
            full_s_poly_values_shifted[(sorted_copy_start - 1)..(witness_len - 1)]
                .copy_from_slice(&s_poly_values_aggregated);

            assert!(full_s_poly_values[0].is_zero());

            let s_poly_monomial = {
                let mon = Polynomial::from_values(full_s_poly_values.clone())?;
                let mon = mon.ifft_using_bitreversed_ntt(
                    &worker,
                    &omegas_inv_bitreversed,
                    &E::Fr::one(),
                )?;

                mon
            };

            (
                s_poly_monomial,
                Polynomial::from_values_unpadded(full_s_poly_values)?,
                Polynomial::from_values_unpadded(full_s_poly_values_shifted)?,
            )
        };

        // compute grand product
        let mut beta_plus_one = beta.clone();
        beta_plus_one.add_assign(&E::Fr::one());
        let mut gamma_beta = gamma;
        gamma_beta.mul_assign(&beta_plus_one);

        let expected = gamma_beta.pow([(required_domain_size - 1) as u64]);

        let mut z_num = {
            // (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega))
            let mut t = t_poly_values.as_ref().clone();
            t.add_assign_scaled(&worker, t_poly_values_shifted.as_ref(), &beta);
            t.add_constant(&worker, &gamma_beta);

            let mut tmp = f_poly_values_aggregated.clone();
            tmp.add_constant(&worker, &gamma);
            tmp.scale(&worker, beta_plus_one);

            t.mul_assign(&worker, &tmp);
            drop(tmp);

            t
        };

        let z_den = {
            // (\gamma*(1 + \beta) + s(x) + \beta * s(x*omega)))
            let mut t = s_poly_unpadded_values.clone();
            t.add_assign_scaled(&worker, &s_shifted_unpadded_values, &beta);
            t.add_constant(&worker, &gamma_beta);

            t.batch_inversion(&worker)?;

            t
        };

        z_num.mul_assign(&worker, &z_den);
        drop(z_den);

        let z = z_num.calculate_shifted_grand_product(&worker)?;
        drop(z_num);

        assert!(z.size().is_power_of_two());

        assert_eq!(z.as_ref()[0], E::Fr::one());
        assert_eq!(*z.as_ref().last().unwrap(), expected);

        let z_poly_in_monomial_form =
            z.ifft_using_bitreversed_ntt(&worker, &omegas_inv_bitreversed, &E::Fr::one())?;

        // Prove lookup identities

        // Numerator degree is at max 4n, so it's < 4n after division

        // ( Z(x*omega)*(\gamma*(1 + \beta) + s(x) + \beta * s(x*omega))) -
        // - Z(x) * (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega)) )*(X - omega^{n-1})

        // Q(x) = (Z(x*W) - Z(x) * (..)/(...) )/ X^n-1
        // degree of grand product is n+1
        // this is a relation that should be valid over domain except last point
        // vanishing poly is like X^n-1/X-w^{n-1}
        // Q(x) = (Z(x*W) - Z(x) * (..)/(...) )/ ((X^n-1)*X-w^{n-1})
        // Q(x) = (X-w^{n-1})*(Z(x*W) - Z(x) * (..)/(...) )/ ((X^n-1))

        // Z(x*omega)*(\gamma*(1 + \beta) + s(x) + \beta * s(x*omega)))

        let s_lde = s_poly_monomial
            .clone()
            .bitreversed_lde_using_bitreversed_ntt(
                &worker,
                lde_factor,
                &omegas_bitreversed,
                &coset_factor,
            )?;

        let s_lde_shifted = s_lde.clone_shifted_assuming_bitreversed(lde_factor, &worker)?;

        let z_lde = z_poly_in_monomial_form
            .clone()
            .bitreversed_lde_using_bitreversed_ntt(
                &worker,
                lde_factor,
                &omegas_bitreversed,
                &coset_factor,
            )?;

        let z_lde_shifted = z_lde.clone_shifted_assuming_bitreversed(lde_factor, &worker)?;

        let mut lookup_first_quotient_rhs_denum_part = s_lde;
        lookup_first_quotient_rhs_denum_part.add_assign_scaled(&worker, &s_lde_shifted, &beta);
        lookup_first_quotient_rhs_denum_part.add_constant(&worker, &gamma_beta);
        lookup_first_quotient_rhs_denum_part.mul_assign(&worker, &z_lde_shifted);

        drop(s_lde_shifted);
        drop(z_lde_shifted);

        let t_lde = t_poly_monomial
            .as_ref()
            .clone()
            .bitreversed_lde_using_bitreversed_ntt(
                &worker,
                lde_factor,
                &omegas_bitreversed,
                &coset_factor,
            )?;

        let t_lde_shifted = t_lde.clone_shifted_assuming_bitreversed(lde_factor, &worker)?;

        let f_lde = {
            // add up ldes of a,b,c and table_type poly and multiply by selector

            // let a_ref = get_from_map_unchecked(
            //     PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(0)),
            //     &ldes_storage,
            // );
            for p in 0..3 {
                let idx = PolyIdentifier::VariablesPolynomial(p);
                ensure_single_poly_in_map_or_create(
                    &worker,
                    required_domain_size,
                    &omegas_bitreversed,
                    lde_factor,
                    coset_factor,
                    &monomials_storage,
                    &mut ldes_storage,
                    idx,
                )?;
            }

            let a_ref = get_from_map_unchecked(
                PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(0)),
                &ldes_storage,
            );
            let mut tmp = a_ref.clone();
            drop(a_ref);

            let mut current = eta;

            let b_ref = get_from_map_unchecked(
                PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(1)),
                &ldes_storage,
            );

            tmp.add_assign_scaled(&worker, b_ref, &current);

            drop(b_ref);
            current.mul_assign(&eta);

            let c_ref = get_from_map_unchecked(
                PolynomialInConstraint::from_id(PolyIdentifier::VariablesPolynomial(2)),
                &ldes_storage,
            );

            tmp.add_assign_scaled(&worker, c_ref, &current);

            drop(c_ref);
            current.mul_assign(&eta);

            let table_type_lde = table_type_mononial
                .as_ref()
                .clone()
                .bitreversed_lde_using_bitreversed_ntt(
                    &worker,
                    lde_factor,
                    &omegas_bitreversed,
                    &coset_factor,
                )?;

            tmp.add_assign_scaled(&worker, &table_type_lde, &current);

            drop(table_type_lde);

            let lookup_selector_lde = selector_poly
                .as_ref()
                .clone()
                .bitreversed_lde_using_bitreversed_ntt(
                    &worker,
                    lde_factor,
                    &omegas_bitreversed,
                    &coset_factor,
                )?;

            tmp.mul_assign(&worker, &lookup_selector_lde);

            drop(lookup_selector_lde);

            tmp
        };

        //  - Z(x) * (\beta + 1) * (\gamma + f(x)) * (\gamma(1 + \beta) + t(x) + \beta * t(x*omega))
        let mut lookup_first_quotient_rhs_num_part = f_lde;
        lookup_first_quotient_rhs_num_part.add_constant(&worker, &gamma);
        lookup_first_quotient_rhs_num_part.mul_assign(&worker, &z_lde);
        lookup_first_quotient_rhs_num_part.scale(&worker, beta_plus_one);

        let mut tmp = t_lde;
        tmp.add_assign_scaled(&worker, &t_lde_shifted, &beta);
        tmp.add_constant(&worker, &gamma_beta);
        lookup_first_quotient_rhs_num_part.mul_assign(&worker, &tmp);
        drop(tmp);
        drop(t_lde_shifted);

        lookup_first_quotient_rhs_denum_part
            .sub_assign(&worker, &lookup_first_quotient_rhs_num_part);

        // multiply by (X - omega^{n-1})
        let last_omega = domain.generator.pow(&[(required_domain_size - 1) as u64]);
        let mut x_minus_last_omega = x_poly_lde_bitreversed.clone();
        x_minus_last_omega.sub_constant(&worker, &last_omega);
        lookup_first_quotient_rhs_denum_part.mul_assign(&worker, &x_minus_last_omega);
        drop(x_minus_last_omega);
        let lookup_first_quotient_monomial = compute_quotient_monomial(
            worker,
            lookup_first_quotient_rhs_denum_part,
            &inverse_divisor_on_coset_lde_natural_ordering,
            coset_factor,
            lde_factor,
            required_domain_size,
            3,
        )?;
        // TODO
        assert_eq!(
            lookup_first_quotient_monomial.size(),
            3 * required_domain_size
        );
        // check that (Z(x) - 1) * L_{0} == 0
        // TODO copy-perm can reuse L0
        let l_0_coset_lde_bitreversed = l_0_coset_lde_bitreversed.clone();
        let mut tmp = z_lde.clone();
        tmp.sub_constant(&worker, &E::Fr::one());
        tmp.mul_assign(&worker, &l_0_coset_lde_bitreversed);
        drop(l_0_coset_lde_bitreversed);
        let lookup_second_quotient_monomial = compute_quotient_monomial(
            worker,
            tmp,
            &inverse_divisor_on_coset_lde_natural_ordering,
            coset_factor,
            lde_factor,
            required_domain_size,
            1,
        )?;
        assert_eq!(lookup_second_quotient_monomial.size(), required_domain_size);

        // check that (Z(x) - expected) * L_{n-1}  == 0
        let l_last = calculate_lagrange_poly::<E::Fr>(
            &worker,
            required_domain_size.next_power_of_two(),
            required_domain_size - 1,
        )?;

        let l_last_coset_lde_bitreversed = l_last.bitreversed_lde_using_bitreversed_ntt(
            &worker,
            lde_factor,
            &omegas_bitreversed,
            &coset_factor,
        )?;

        let mut tmp = z_lde.clone();
        tmp.sub_constant(&worker, &expected);
        tmp.mul_assign(&worker, &l_last_coset_lde_bitreversed);
        let lookup_third_quotient_monomial = compute_quotient_monomial(
            worker,
            tmp,
            &inverse_divisor_on_coset_lde_natural_ordering,
            coset_factor,
            lde_factor,
            required_domain_size,
            1,
        )?;
        assert_eq!(lookup_third_quotient_monomial.size(), required_domain_size);

        let lookup_monomials: LookupMonomials<E::Fr> = LookupMonomials {
            s_poly_monomial,
            grand_product_monomial: z_poly_in_monomial_form,
            first_quotient: lookup_first_quotient_monomial,
            second_quotient: lookup_second_quotient_monomial,
            third_quotient: lookup_third_quotient_monomial,
        };

        Some(lookup_monomials)
    } else {
        None
    };
    let second_round_monomials = SecondRoundMonomials {
        copy_permutation: copy_permutation_monomials,
        lookup: lookup_monomials,
    };

    let flattened_second_round_monomials = second_round_monomials.flatten();

    assert_eq!(
        flattened_second_round_monomials.len(),
        num_second_round_polys
    );
    let c2 =
        combine_mixed_degree_monomials(&flattened_second_round_monomials, required_domain_size)?;
    let expected_degree =
    num_second_round_polys * copy_perm_quotient_degree * required_domain_size + num_second_round_polys - 1;
    assert_eq!(c2.size(), expected_degree);
    assert!(c2.size() <= max_combined_degree);
    let c2_commitment = commit_using_monomials(&c2, mon_crs, &worker)?;

    commit_point_as_xy::<E, T>(&mut transcript, &c2_commitment);

    // draw opening point
    let r = transcript.get_challenge();

    let interpolation_size_of_setup = num_setup_polys.next_power_of_two();
    let interpolation_size_of_first_round = num_first_round_polys.next_power_of_two();
    let interpolation_size_of_second_round = num_second_round_polys;
    let power = lcm(&[
        interpolation_size_of_setup,
        interpolation_size_of_first_round,
        interpolation_size_of_second_round,
    ]);
    println!("LCM {power}");
    let z = r.pow(&[power as u64]);
    let mut z_omega = z;
    z_omega.mul_assign(&domain.generator);

    // setup  evaluations
    let setup_evaluations = evaluate_setup_monomials(&worker, &setup, z, z_omega);
    assert_eq!(
        setup_evaluations.interpolation_size(),
        interpolation_size_of_setup
    );
    // first round evaluations
    let first_round_evaluations =
        evaluate_first_round_polynomials(&worker, &first_round_monomials, z, z_omega);
    assert_eq!(
        first_round_evaluations.interpolation_size(),
        interpolation_size_of_first_round
    );
    // second round evaluations
    let second_round_evaluations =
        evaluate_second_round_polynomials(&worker, &second_round_monomials, z, z_omega);
    assert_eq!(
        second_round_evaluations.interpolation_size(),
        interpolation_size_of_second_round
    );

    let all_evaluations = flatten_all_evaluations(
        &setup_evaluations,
        &first_round_evaluations,
        &second_round_evaluations,
    );
    let c0 = compute_combined_setup_monomial(&setup, required_domain_size)?;
    let expected_degree = num_setup_polys * required_domain_size + num_setup_polys - 1;
    assert_eq!(c0.size(), expected_degree);
    assert!(c0.size() <= max_combined_degree);
    for eval in all_evaluations.iter() {
        transcript.commit_field_element(eval);
    }
    // openings
    let (h0, h1, h2) = compute_opening_points(
        r, 
        z, 
        z_omega, 
        power, 
        interpolation_size_of_setup, 
        interpolation_size_of_first_round, 
        interpolation_size_of_second_round, 
        required_domain_size,
        setup_evaluations.requires_opening_at_shifted_point(),
        first_round_evaluations.requires_opening_at_shifted_point()
    );
    // draw challenge for f(x)
    let alpha = transcript.get_challenge();
    let mut alpha_squared = alpha.clone();
    alpha_squared.mul_assign(&alpha);

    let (setup_r, first_round_r, second_round_r) = construct_r_monomials(
        &setup_evaluations,
        &first_round_evaluations,
        &second_round_evaluations,
        h0,
        h1, 
        h2
    );

    let setup_r_monomial = Polynomial::from_coeffs_unpadded(setup_r)?;
    let first_round_r_monomial = Polynomial::from_coeffs_unpadded(first_round_r)?;
    let second_round_r_monomial = Polynomial::from_coeffs_unpadded(second_round_r)?;
    let [setup_omega, first_round_omega, second_round_omega] = compute_generators(
        interpolation_size_of_setup, interpolation_size_of_first_round, interpolation_size_of_second_round
    );
    if SANITY_CHECK {
        for (round_idx, (c, r, h, interpolation_set_size, omega, requirest_opening_at_shifted_point)) in [
            (&c0, &setup_r_monomial, h0, interpolation_size_of_setup, setup_omega, setup_evaluations.requires_opening_at_shifted_point()),
            (&c1, &first_round_r_monomial, h1, interpolation_size_of_first_round, first_round_omega, first_round_evaluations.requires_opening_at_shifted_point()),
            (&c2, &second_round_r_monomial, (h2.0, Some(h2.1)), interpolation_size_of_second_round, second_round_omega, true),
            ].into_iter().enumerate(){
                        
            let (h, h_shifted) = h;
            let mut current_omega = h;
            let mut current_omega_shifted = if requirest_opening_at_shifted_point{
                h_shifted.expect("h shifted")
            }else{
                E::Fr::zero()
            };

            for _ in 0..interpolation_set_size {
                assert_eq!(
                    c.evaluate_at(&worker, current_omega),
                    horner_evaluation(r.as_ref(), current_omega),
                    "{round_idx} round evaluations mismatches"
                );
                current_omega.mul_assign(&omega);

                if requirest_opening_at_shifted_point{
                    assert_eq!(
                        c.evaluate_at(&worker, current_omega_shifted),
                        horner_evaluation(r.as_ref(), current_omega_shifted),
                        "{round_idx} round shifted evaluations mismatches"
                    );  
                    current_omega_shifted.mul_assign(&omega);
                }
            }
        }
    }
    // f(x) = Z_{T\S0}(x)*(C0(x)- r0(x)) + alpha*Z_{T\S1}(x)*(C1(x)- r1(x)) + alpha^2*Z_{T\S2}(x)*(C2(x)- r2(x))

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
        first_round_evaluations.requires_opening_at_shifted_point()
    );

    let sparse_poly_degree: usize = sparse_polys_for_first_round
        .iter()
        .map(|(degree, _)| *degree)
        .sum();

    let f_poly_degree = [
        (c0.size() + sparse_poly_degree),
        (c1.size() + sparse_poly_degree),
        (c2.size() + sparse_poly_degree),
    ]
    .into_iter()
    .max()
    .map(|d| d)
    .unwrap();

    // we will compute f(x) in monomial form

    let mut f_poly = Polynomial::from_coeffs_unpadded(vec![E::Fr::zero(); f_poly_degree])?;

    // Z_{T\S0}(x)*C0(x)- r0(x)
    let mut tmp = c0.clone();
    tmp.sub_assign(&worker, &setup_r_monomial);
    let tmp = multiply_monomial_with_multiple_sparse_polys(worker, &tmp, &sparse_polys_for_setup)?;
    f_poly.add_assign(&worker, &tmp);

    // + alpha*Z_{T\S1}(x)*(C1(x)- r1(x))
    let mut tmp = c1.clone();
    tmp.sub_assign(&worker, &first_round_r_monomial);
    let tmp =
        multiply_monomial_with_multiple_sparse_polys(worker, &tmp, &sparse_polys_for_first_round)?;
    f_poly.add_assign_scaled(&worker, &tmp, &alpha);

    // + alpha^2*Z_{T\S2}(x)*(C2(x)- r2(x))
    let mut tmp = c2.clone();
    tmp.sub_assign(&worker, &second_round_r_monomial);
    let tmp =
        multiply_monomial_with_multiple_sparse_polys(worker, &tmp, &sparse_polys_for_second_round)?;
    f_poly.add_assign_scaled(&worker, &tmp, &alpha_squared);

    // we will do W(x) = f(x) / Z_T in monomial form
    let w = divide_by_multiple_higher_degree_sparse_polys(f_poly.as_ref(), sparse_polys.clone());
    let w = Polynomial::from_coeffs_unpadded(w)?;
    if SANITY_CHECK {
        let actual =
            multiply_monomial_with_multiple_sparse_polys(worker, &w, &sparse_polys).unwrap();
        assert_eq!(f_poly.size(), actual.size());
        assert_eq!(f_poly.as_ref(), actual.as_ref());
    }
    let w_commitment = commit_using_monomials(&w, mon_crs, worker)?;
    commit_point_as_xy::<E, T>(&mut transcript, &w_commitment);

    // draw challenge for L(x)
    let y = transcript.get_challenge();

    let (flattened_lagrange_basis_inverses, lagrange_basis_evaluations) =
        precompute_all_lagrange_basis_evaluations(
            interpolation_size_of_setup,
            interpolation_size_of_first_round,
            interpolation_size_of_second_round,
            h0,
            h1,
            h2,
            y,
            setup_evaluations.requires_opening_at_shifted_point(),
            first_round_evaluations.requires_opening_at_shifted_point()
        );

    // L(x) = Z_{T\S0}(y)*(C0(x)- r0(y)) + alpha*Z_{T\S1}(y)*(C1(x)- r1(y)) + alpha^2*Z_{T\S2}(y)*(C2(x)- r2(y)) - Z_T(y)*W(x)
    let [setup_r_at_y, first_round_r_at_y, second_round_r_at_y] =
    evaluate_r_polys_at_point_with_flattened_evals_and_precomputed_basis(
            &all_evaluations,
            num_setup_polys,
            num_first_round_polys,
            num_second_round_polys,
            h0,
            h1,
            h2,
            lagrange_basis_evaluations,
            setup_evaluations.requires_opening_at_shifted_point(),
            first_round_evaluations.requires_opening_at_shifted_point()
        );

    if SANITY_CHECK {
        assert_eq!(
            horner_evaluation(setup_r_monomial.as_ref(), y),
            setup_r_at_y
        );
        assert_eq!(
            horner_evaluation(first_round_r_monomial.as_ref(), y),
            first_round_r_at_y
        );
        assert_eq!(
            horner_evaluation(second_round_r_monomial.as_ref(), y),
            second_round_r_at_y
        );
    }

    let sparse_polys_for_setup_at_y = evaluate_multiple_sparse_polys(sparse_polys_for_setup, y);
    let sparse_polys_for_first_round_at_y =
        evaluate_multiple_sparse_polys(sparse_polys_for_first_round, y);
    let sparse_polys_for_second_round_at_y =
        evaluate_multiple_sparse_polys(sparse_polys_for_second_round, y);

    let l_poly_degree = [c0.size(), c1.size(), c2.size()]
        .into_iter()
        .max()
        .map(|d| d)
        .unwrap();

    let mut l_poly = Polynomial::from_coeffs_unpadded(vec![E::Fr::zero(); l_poly_degree])?;
    // Z_{T\S0}(y)*(C0(x)- r0(y))
    let mut tmp = c0;
    tmp.as_mut()[0].sub_assign(&setup_r_at_y);
    tmp.scale(&worker, sparse_polys_for_setup_at_y);
    l_poly.add_assign(&worker, &tmp);

    // + (alpha*Z_{T\S1}(y))*(C1(x)- r1(y))
    let mut tmp = c1;
    tmp.as_mut()[0].sub_assign(&first_round_r_at_y);
    let mut factor = alpha;
    factor.mul_assign(&sparse_polys_for_first_round_at_y);
    l_poly.add_assign_scaled(&worker, &tmp, &factor);

    // + (alpha^2*Z_{T\S2}(y))*(C2(x)- r2(y))
    let mut tmp = c2;
    tmp.as_mut()[0].sub_assign(&second_round_r_at_y);
    let mut factor = alpha_squared;
    factor.mul_assign(&sparse_polys_for_second_round_at_y);
    l_poly.add_assign_scaled(&worker, &tmp, &factor);

    // -Z_T(y)*W(X)
    let mut w_scaled = w.clone();
    let z_t_at_y = evaluate_multiple_sparse_polys(sparse_polys, y);
    w_scaled.scale(&worker, z_t_at_y);
    l_poly.sub_assign(&worker, &w_scaled);

    // W'(x) = L(x) / (Z_{T\S0}(y)*(x-y))
    // dividing by Z_{T\S0}(y) reduces verifier G1 scalar mul by 1
    // since alpha^0=1 for C0.
    let l_divided_by_y = divide_by_linear_term(l_poly.as_ref(), y);
    let mut w_prime = Polynomial::from_coeffs_unpadded(l_divided_by_y)?;
    let inv_sparse_polys_for_setup_at_y = sparse_polys_for_setup_at_y.inverse().unwrap();
    w_prime.scale(&worker, inv_sparse_polys_for_setup_at_y);
    if SANITY_CHECK {
        let mut actual = multiply_monomial_with_sparse_poly(worker, &w_prime, 1, y)?;
        actual.scale(&worker, sparse_polys_for_setup_at_y);
        assert_eq!(l_poly.as_ref(), &actual.as_ref()[..l_poly.as_ref().len()]);
    }

    let w_prime_commitment = commit_using_monomials(&w_prime, mon_crs, worker)?;

    let mut proof = FflonkProof::<E, C>::empty();
    proof.n = assembly.n();
    proof.inputs = input_values.clone();

    proof.setup_evaluations = setup_evaluations;
    proof.first_round_evaluations = first_round_evaluations;
    proof.second_round_evaluations = second_round_evaluations;

    proof.c1 = c1_commitment;
    proof.c2 = c2_commitment;
    proof.w = w_commitment;
    proof.w_prime = w_prime_commitment;
    proof.lagrange_basis_inverses = flattened_lagrange_basis_inverses;

    Ok(proof)
}
