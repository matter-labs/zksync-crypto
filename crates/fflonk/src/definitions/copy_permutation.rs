use super::*;

pub struct CopyPermutationMonomials<F: PrimeField> {
    pub grand_product_monomial: Polynomial<F, Coefficients>,
    pub first_quotient: Polynomial<F, Coefficients>,
    pub second_quotient: Polynomial<F, Coefficients>,
}

impl<F: PrimeField> CopyPermutationMonomials<F> {
    pub fn flatten(&self) -> Vec<&Polynomial<F, Coefficients>> {
        let mut flattened = vec![];

        flattened.push(&self.grand_product_monomial);
        flattened.push(&self.first_quotient);
        flattened.push(&self.second_quotient);

        flattened
    }
}

#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize, PartialEq, Eq)]
#[serde(bound = "")]
pub struct CopyPermutationEvaluations<F: PrimeField> {
    pub grand_product_at_z: F,
    pub grand_product_at_z_omega: F,
    pub first_quotient_at_z_omega: F,
    pub second_quotient_at_z_omega: F,
}
#[derive(Clone, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct CopyPermutationEvaluationOffsets {
    pub grand_product_at_z: usize,
    pub grand_product_at_z_omega: usize,
    pub first_quotient_at_z_omega: usize,
    pub second_quotient_at_z_omega: usize,
}
impl<F: PrimeField> Flatten<F> for CopyPermutationEvaluations<F> {
    fn flatten(&self) -> (Vec<F>, Vec<F>) {
        let mut evals = vec![];
        let mut shifted_evals = vec![];

        evals.push(self.grand_product_at_z);
        shifted_evals.push(self.grand_product_at_z_omega);
        shifted_evals.push(self.first_quotient_at_z_omega);
        shifted_evals.push(self.second_quotient_at_z_omega);
        assert_eq!(evals.len() + 2, shifted_evals.len());

        (evals, shifted_evals)
    }
}

pub fn evaluate_copy_permutation_monomials<'a, F: PrimeField>(worker: &Worker, monomial_storage: &CopyPermutationMonomials<F>, z: F, z_omega: F) -> CopyPermutationEvaluations<F> {
    let grand_product_at_z = monomial_storage.grand_product_monomial.evaluate_at(worker, z);
    let grand_product_at_z_omega = monomial_storage.grand_product_monomial.evaluate_at(worker, z_omega);
    let first_quotient_at_z_omega = monomial_storage.first_quotient.evaluate_at(worker, z_omega);
    let second_quotient_at_z_omega = monomial_storage.second_quotient.evaluate_at(worker, z_omega);

    CopyPermutationEvaluations {
        grand_product_at_z,
        grand_product_at_z_omega,
        first_quotient_at_z_omega,
        second_quotient_at_z_omega,
    }
}
