use super::*;
mod copy_permutation;
mod gates;
mod lookup;
pub mod proof;
pub mod setup;

pub(crate) use copy_permutation::*;
pub use gates::*;
pub(crate) use lookup::*;
pub use proof::*;
pub use setup::*;

pub trait Flatten<T> {
    fn flatten(&self) -> (Vec<T>, Vec<T>);
}

pub struct FirstRoundMonomials<'a, E: Engine> {
    pub trace_and_gate_monomials: TraceAndGateMonomials<'a, E>,
}

impl<'a, E: Engine> FirstRoundMonomials<'a, E> {
    pub fn flatten(&self) -> Vec<&Polynomial<E::Fr, Coefficients>> {
        self.trace_and_gate_monomials.flatten()
    }
}

#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize, PartialEq, Eq)]
#[serde(bound = "")]
pub struct FirstRoundEvaluations<F: PrimeField> {
    pub trace_and_gate_evaluations: TraceAndGateEvaluations<F>,
}

impl<F: PrimeField> FirstRoundEvaluations<F> {
    pub fn flatten(&self) -> (Vec<F>, Vec<F>) {
        self.trace_and_gate_evaluations.flatten()
    }

    pub fn requires_opening_at_shifted_point(&self) -> bool {
        self.trace_and_gate_evaluations
            .requires_opening_at_shifted_point()
    }

    pub fn interpolation_size(&self) -> usize {
        self.trace_and_gate_evaluations.interpolation_size()
    }
}

pub fn evaluate_first_round_polynomials<E: Engine>(
    worker: &Worker,
    monomials_storage: &FirstRoundMonomials<E>,
    z: E::Fr,
    z_omega: E::Fr,
) -> FirstRoundEvaluations<E::Fr> {
    let trace_and_gates = evaluate_trace_and_gate_monomials(
        worker,
        &monomials_storage.trace_and_gate_monomials,
        z,
        z_omega,
    );

    FirstRoundEvaluations {
        trace_and_gate_evaluations: trace_and_gates,
    }
}

pub struct SecondRoundMonomials<F: PrimeField> {
    pub copy_permutation: CopyPermutationMonomials<F>,
    pub lookup: Option<LookupMonomials<F>>,
}

impl<F: PrimeField> SecondRoundMonomials<F> {
    pub fn flatten(&self) -> Vec<&Polynomial<F, Coefficients>> {
        let mut flattened = self.copy_permutation.flatten();
        self.lookup
            .as_ref()
            .map(|lookup| flattened.extend(lookup.flatten()));

        flattened
    }
}

#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize, PartialEq, Eq)]
#[serde(bound = "")]
pub struct SecondRoundEvaluations<F: PrimeField> {
    pub copy_permutation_evaluations: CopyPermutationEvaluations<F>,
    pub lookup_evaluations: Option<LookupEvaluations<F>>,
}

impl<F: PrimeField> SecondRoundEvaluations<F> {
    pub fn flatten(&self) -> (Vec<F>, Vec<F>) {
        let (mut evals, mut evals_shifted) = self.copy_permutation_evaluations.flatten();
        if let Some(ref lookup_evaluations) = self.lookup_evaluations {
            let (lookup_evals, lookup_evals_shifted) = lookup_evaluations.flatten();
            evals.extend(lookup_evals);
            evals_shifted.extend(lookup_evals_shifted);
        }

        (evals, evals_shifted)
    }

    pub fn interpolation_size(&self) -> usize {
        let num_copy_perm_polys = 3usize;
        let num_lookup_polys = self.lookup_evaluations.as_ref().map(|_| 5).unwrap_or(0);

        num_copy_perm_polys + num_lookup_polys
    }
}

pub fn evaluate_second_round_polynomials<F: PrimeField>(
    worker: &Worker,
    monomials_storage: &SecondRoundMonomials<F>,
    z: F,
    z_omega: F,
) -> SecondRoundEvaluations<F> {
    let copy_permutation = evaluate_copy_permutation_monomials(
        worker,
        &monomials_storage.copy_permutation,
        z,
        z_omega,
    );
    let lookup = monomials_storage.lookup.as_ref().map(|lookup_evaluations| {
        evaluate_lookup_monomials(worker, lookup_evaluations, z, z_omega)
    });

    SecondRoundEvaluations {
        copy_permutation_evaluations: copy_permutation,
        lookup_evaluations: lookup,
    }
}
