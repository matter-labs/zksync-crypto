use super::*;

#[derive(serde::Serialize)]
pub struct LookupMonomials<F: PrimeField> {
    pub s_poly_monomial: Polynomial<F, Coefficients>,
    pub grand_product_monomial: Polynomial<F, Coefficients>,
    pub first_quotient: Polynomial<F, Coefficients>,
    pub second_quotient: Polynomial<F, Coefficients>,
    pub third_quotient: Polynomial<F, Coefficients>,
}

impl<F: PrimeField> LookupMonomials<F> {
    pub fn flatten(&self) -> Vec<&Polynomial<F, Coefficients>> {
        let mut flattened = vec![];

        flattened.push(&self.s_poly_monomial);
        flattened.push(&self.grand_product_monomial);
        flattened.push(&self.first_quotient);
        flattened.push(&self.second_quotient);
        flattened.push(&self.third_quotient);

        flattened
    }
}

#[derive(Default, Debug, Clone, serde::Serialize, serde::Deserialize, PartialEq, Eq)]
#[serde(bound = "")]
pub struct LookupEvaluations<F: PrimeField> {
    pub s_poly_at_z: F,
    pub s_poly_at_z_omega: F,
    pub grand_product_at_z: F,
    pub grand_product_at_z_omega: F,
    pub first_quotient_at_z: F,
    pub first_quotient_at_z_omega: F,
    pub second_quotient_at_z: F,
    pub second_quotient_at_z_omega: F,
    pub third_quotient_at_z: F,
    pub third_quotient_at_z_omega: F,
}

#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct LookupEvaluationOffsets {
    pub s_poly_at_z: usize,
    pub s_poly_at_z_omega: usize,
    pub grand_product_at_z: usize,
    pub grand_product_at_z_omega: usize,
    pub first_quotient_at_z: usize,
    pub first_quotient_at_z_omega: usize,
    pub second_quotient_at_z: usize,
    pub second_quotient_at_z_omega: usize,
    pub third_quotient_at_z: usize,
    pub third_quotient_at_z_omega: usize,
}
impl<F: PrimeField> Flatten<F> for LookupEvaluations<F> {
    fn flatten(&self) -> (Vec<F>, Vec<F>) {
        let mut evals = vec![];
        let mut shifted_evals = vec![];
        evals.push(self.s_poly_at_z);
        shifted_evals.push(self.s_poly_at_z_omega);
        evals.push(self.grand_product_at_z);
        shifted_evals.push(self.grand_product_at_z_omega);
        evals.push(self.first_quotient_at_z);
        shifted_evals.push(self.first_quotient_at_z_omega);
        evals.push(self.second_quotient_at_z);
        shifted_evals.push(self.second_quotient_at_z_omega);
        evals.push(self.third_quotient_at_z);
        shifted_evals.push(self.third_quotient_at_z_omega);
        assert_eq!(evals.len(), shifted_evals.len());

        (evals, shifted_evals)
    }
}

pub fn evaluate_lookup_monomials<'a, F: PrimeField>(
    worker: &Worker,
    monomials_storage: &LookupMonomials<F>,
    z: F,
    z_omega: F,
) -> LookupEvaluations<F> {
    let s_poly_at_z = monomials_storage.s_poly_monomial.evaluate_at(worker, z);
    let s_poly_at_z_omega = monomials_storage
        .s_poly_monomial
        .evaluate_at(worker, z_omega);

    let grand_product_at_z = monomials_storage
        .grand_product_monomial
        .evaluate_at(worker, z);
    let grand_product_at_z_omega = monomials_storage
        .grand_product_monomial
        .evaluate_at(worker, z_omega);

    let first_quotient_at_z = monomials_storage.first_quotient.evaluate_at(worker, z);
    let first_quotient_at_z_omega = monomials_storage
        .first_quotient
        .evaluate_at(worker, z_omega);
    let second_quotient_at_z = monomials_storage.second_quotient.evaluate_at(worker, z);
    let second_quotient_at_z_omega = monomials_storage
        .second_quotient
        .evaluate_at(worker, z_omega);
    let third_quotient_at_z = monomials_storage.third_quotient.evaluate_at(worker, z);
    let third_quotient_at_z_omega = monomials_storage
        .third_quotient
        .evaluate_at(worker, z_omega);

    LookupEvaluations {
        s_poly_at_z,
        s_poly_at_z_omega,
        grand_product_at_z,
        grand_product_at_z_omega,
        first_quotient_at_z,
        first_quotient_at_z_omega,
        third_quotient_at_z,
        third_quotient_at_z_omega,
        second_quotient_at_z,
        second_quotient_at_z_omega,
    }
}
