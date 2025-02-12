use lazy_static::lazy_static;

use self::{
    algebraic_torus::TorusTestCases,
    field_extensions::{Fq12TestCases, Fq2TestCases, Fq6TestCases},
};

pub mod algebraic_torus;
pub mod field_extensions;
pub mod types;

// All tests gathered in one place
lazy_static! {
    /// Test cases for `Fq2` operations
    pub static ref FQ2_TEST_CASES: Fq2TestCases = field_extensions::load_fq2_test_cases();
    /// Test cases for `Fq6` operations
    pub static ref FQ6_TEST_CASES: Fq6TestCases = field_extensions::load_fq6_test_cases();
    /// Test cases for `Fq12` operations
    pub static ref FQ12_TEST_CASES: Fq12TestCases = field_extensions::load_fq12_test_cases();
    /// Test cases for algebraic torus operations
    pub static ref TORUS_TEST_CASES: TorusTestCases = algebraic_torus::load_torus_test_cases();
}
