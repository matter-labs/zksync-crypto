use std::io::{Read, Write};

use super::*;
use bellman::plonk::better_better_cs::cs::SetupAssembly;
use bellman::plonk::better_cs::keys::read_curve_affine;
use bellman::plonk::better_cs::keys::write_curve_affine;

pub struct FflonkSetup<E: Engine, C: Circuit<E>> {
    pub original_setup: Setup<E, C>,
    pub c0_commitment: E::G1Affine,
}

impl<E: Engine, C: Circuit<E>> FflonkSetup<E, C> {
    pub fn create_setup<P: PlonkConstraintSystemParams<E>, MG: MainGate<E>, S: SynthesisMode>(
        setup_assembly: &Assembly<E, P, MG, S>,
        worker: &Worker,
        crs: &Crs<E, CrsForMonomialForm>,
    ) -> Result<Self, SynthesisError> {
        assert!(S::PRODUCE_SETUP, "Assembly should hold setup values");
        let setup = setup_assembly.create_setup(worker)?;
        let domain_size = setup_assembly.n() + 1;
        assert!(domain_size.is_power_of_two());
        assert!(domain_size <= crs.g1_bases.len());

        let combined_setup_monomial = compute_combined_setup_monomial(&setup, domain_size)?;
        let c0 = commit_using_monomials(&combined_setup_monomial, crs, worker)?;

        Ok(Self {
            original_setup: setup,
            c0_commitment: c0,
        })
    }

    pub fn read<R: Read>(mut reader: R) -> std::io::Result<Self> {
        let original_setup: Setup<E, C> = Setup::read(&mut reader)?;
        let c0_commitment = read_curve_affine(&mut reader)?;
        Ok(Self {
            original_setup,
            c0_commitment,
        })
    }

    pub fn write<W: Write>(&self, mut writer: W) -> std::io::Result<()> {
        self.original_setup.write(&mut writer)?;
        write_curve_affine(&self.c0_commitment, &mut writer)?;
        Ok(())
    }
}

pub fn compute_combined_setup_monomial<E: Engine, C: Circuit<E>>(
    setup: &Setup<E, C>,
    domain_size: usize,
) -> Result<Polynomial<E::Fr, Coefficients>, SynthesisError> {
    assert!(domain_size.is_power_of_two());
    assert_eq!(domain_size, setup.n + 1);

    let flattened_setup_monomials = flatten_setup_monomials(&setup);
    combine_monomials(&flattened_setup_monomials, domain_size)
}
pub fn num_setup_polys<E: Engine, C: Circuit<E>>(setup: &Setup<E, C>) -> usize {
    let gates = sorted_gates_from_circuit_definitions::<_, C>();
    let has_custom_gate = gates.len() > 1;
    let has_lookup = setup.total_lookup_entries_length > 0;
    let main_gate = gates[0].clone();
    let num_gate_setups = main_gate.setup_polynomials().len();
    let num_gate_selectors = if has_custom_gate { 2 } else { 0 };
    let num_permutations = main_gate.variable_polynomials().len();
    assert_eq!(num_permutations, setup.state_width);
    let num_lookup_polys = if has_lookup { 1 + 1 + 4 } else { 0 };
    let num_setup_polys = flatten_setup_monomials(setup).len();
    assert_eq!(
        num_setup_polys,
        num_gate_setups + num_gate_selectors + num_permutations + num_lookup_polys
    );

    num_setup_polys
}

pub fn flatten_setup_monomials<E: Engine, C: Circuit<E>>(
    setup: &Setup<E, C>,
) -> Vec<&Polynomial<E::Fr, Coefficients>> {
    let mut setup_monomials = vec![];
    for poly in setup.gate_setup_monomials.iter() {
        setup_monomials.push(poly)
    }
    for poly in setup.gate_selectors_monomials.iter() {
        setup_monomials.push(poly)
    }
    for poly in setup.permutation_monomials.iter() {
        setup_monomials.push(poly)
    }

    setup
        .lookup_selector_monomial
        .as_ref()
        .map(|m| setup_monomials.push(m));

    let has_lookup = setup.lookup_selector_monomial.is_some();
    if has_lookup {
        assert!(setup.total_lookup_entries_length > 0);
        for poly in setup.lookup_tables_monomials.iter() {
            setup_monomials.push(poly)
        }
    }
    setup
        .lookup_table_type_monomial
        .as_ref()
        .map(|m| setup_monomials.push(m));

    setup_monomials
}

#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize, PartialEq, Eq)]
#[serde(bound = "")]
pub struct SetupEvaluations<F: PrimeField> {
    pub gate_setups_at_z: Vec<F>,
    pub gate_setups_at_z_omega: Option<Vec<F>>,
    pub gate_selectors_at_z: Vec<F>,
    pub gate_selectors_at_z_omega: Option<Vec<F>>,
    pub permutations_at_z: Vec<F>,
    pub permutations_at_z_omega: Option<Vec<F>>,
    pub lookup_selector_at_z: Option<F>,
    pub lookup_selector_at_z_omega: Option<F>,
    pub lookup_tables_at_z: Option<Vec<F>>,
    pub lookup_tables_at_z_omega: Option<Vec<F>>,
    pub lookup_table_type_at_z: Option<F>,
    pub lookup_table_type_at_z_omega: Option<F>,
    pub num_polys: usize,
}

impl<F: PrimeField> Flatten<F> for SetupEvaluations<F> {
    fn flatten(&self) -> (Vec<F>, Vec<F>) {
        let mut evals = vec![];
        let mut evals_shifted = vec![];

        evals.extend_from_slice(&self.gate_setups_at_z);
        evals.extend_from_slice(&self.gate_selectors_at_z);
        evals.extend_from_slice(&self.permutations_at_z);

        self.lookup_selector_at_z.map(|value| evals.push(value));
        self.lookup_tables_at_z
            .as_ref()
            .map(|values| evals.extend_from_slice(values));
        self.lookup_table_type_at_z.map(|value| evals.push(value));

        assert_eq!(evals.len(), self.num_polys);

        self.gate_setups_at_z_omega
            .as_ref()
            .map(|values| evals_shifted.extend_from_slice(&values));
        self.gate_selectors_at_z_omega
            .as_ref()
            .map(|values| evals_shifted.extend_from_slice(&values));
        self.permutations_at_z_omega
            .as_ref()
            .map(|values| evals_shifted.extend_from_slice(&values));
        self.lookup_selector_at_z_omega
            .map(|value| evals_shifted.push(value));
        self.lookup_tables_at_z_omega
            .as_ref()
            .map(|values| evals_shifted.extend_from_slice(values));
        self.lookup_table_type_at_z_omega
            .map(|value| evals_shifted.push(value));

        if self.requires_opening_at_shifted_point() {
            assert_eq!(
                self.gate_setups_at_z.len(),
                self.gate_setups_at_z_omega.as_ref().unwrap().len()
            );
            assert_eq!(
                self.gate_selectors_at_z.len(),
                self.gate_selectors_at_z_omega.as_ref().unwrap().len()
            );

            assert_eq!(
                self.permutations_at_z.len(),
                self.permutations_at_z_omega.as_ref().unwrap().len()
            );
            assert_eq!(
                self.lookup_tables_at_z.as_ref().unwrap().len(),
                self.lookup_tables_at_z_omega.as_ref().unwrap().len()
            );
            assert!(self.lookup_selector_at_z_omega.is_some());
            assert!(self.lookup_table_type_at_z_omega.is_some());

            assert_eq!(evals.len(), evals_shifted.len());
        }

        (evals, evals_shifted)
    }
}

impl<F: PrimeField> SetupEvaluations<F> {
    pub fn interpolation_size(&self) -> usize {
        self.num_polys.next_power_of_two()
    }

    pub fn requires_opening_at_shifted_point(&self) -> bool {
        self.lookup_tables_at_z_omega.is_some()
    }
}

pub fn requires_setup_polys_opening_at_shifted_point<E: Engine, C: Circuit<E>>(
    vk: &FflonkVerificationKey<E, C>,
) -> bool {
    vk.total_lookup_entries_length > 0
}
pub fn evaluate_setup_monomials<'a, E: Engine, C: Circuit<E>>(
    worker: &Worker,
    setup: &Setup<E, C>,
    z: E::Fr,
    z_omega: E::Fr,
) -> SetupEvaluations<E::Fr> {
    let mut evals = SetupEvaluations::default();
    evals.num_polys = num_setup_polys(setup);

    let has_lookup = setup.total_lookup_entries_length > 0;

    evals.gate_setups_at_z = setup
        .gate_setup_monomials
        .iter()
        .map(|m| m.evaluate_at(worker, z))
        .collect();

    evals.gate_selectors_at_z = setup
        .gate_selectors_monomials
        .iter()
        .map(|m| m.evaluate_at(worker, z))
        .collect();

    if has_custom_gate::<_, C>() {
        assert_eq!(evals.gate_selectors_at_z.len(), 2);
    } else {
        assert_eq!(evals.gate_selectors_at_z.len(), 0);
    }

    evals.permutations_at_z = setup
        .permutation_monomials
        .iter()
        .map(|m| m.evaluate_at(worker, z))
        .collect();

    if has_lookup {
        // if circuit has lookup support, then all round polynomials
        // should be opened at shifted point.

        evals.gate_setups_at_z_omega = Some(
            setup
                .gate_setup_monomials
                .iter()
                .map(|m| m.evaluate_at(worker, z_omega))
                .collect(),
        );

        evals.gate_selectors_at_z_omega = Some(
            setup
                .gate_selectors_monomials
                .iter()
                .map(|m| m.evaluate_at(worker, z_omega))
                .collect(),
        );

        evals.permutations_at_z_omega = Some(
            setup
                .permutation_monomials
                .iter()
                .map(|m| m.evaluate_at(worker, z_omega))
                .collect(),
        );
        assert!(setup.lookup_selector_monomial.is_some());
        assert_eq!(setup.lookup_tables_monomials.len(), 4);
        assert!(setup.lookup_table_type_monomial.is_some());

        evals.lookup_selector_at_z = setup
            .lookup_selector_monomial
            .as_ref()
            .map(|m| m.evaluate_at(worker, z));

        evals.lookup_selector_at_z_omega = setup
            .lookup_selector_monomial
            .as_ref()
            .map(|m| m.evaluate_at(worker, z_omega));

        evals.lookup_tables_at_z = Some(
            setup
                .lookup_tables_monomials
                .iter()
                .map(|m| m.evaluate_at(worker, z))
                .collect(),
        );
        evals.lookup_tables_at_z_omega = Some(
            setup
                .lookup_tables_monomials
                .iter()
                .map(|m| m.evaluate_at(worker, z_omega))
                .collect(),
        );

        evals.lookup_table_type_at_z = setup
            .lookup_table_type_monomial
            .as_ref()
            .map(|m| m.evaluate_at(worker, z));

        evals.lookup_table_type_at_z_omega = setup
            .lookup_table_type_monomial
            .as_ref()
            .map(|m| m.evaluate_at(worker, z_omega));
    }

    evals
}

#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct SetupEvaluationOffsets {
    pub num_gate_setups: usize,
    pub gate_setups_at_z: usize,
    pub gate_setups_at_z_omega: usize,
    pub num_gate_selectors: usize,
    pub gate_selectors_at_z: usize,
    pub gate_selectors_at_z_omega: usize,
    pub num_permutations: usize,
    pub permutations_at_z: usize,
    pub permutations_at_z_omega: usize,
    pub lookup_selector_at_z: usize,
    pub lookup_selector_at_z_omega: usize,
    pub num_lookup_table_cols: usize,
    pub lookup_tables_at_z: usize,
    pub lookup_tables_at_z_omega: usize,
    pub lookup_table_type_at_z: usize,
    pub lookup_table_type_at_z_omega: usize,
    pub q_const_term: usize,
    pub q_dnext_term: usize,
    pub q_ab_term: usize,
    pub q_ac_term: usize,
}
