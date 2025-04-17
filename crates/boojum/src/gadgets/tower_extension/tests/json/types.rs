// Helper utils for testing

use std::sync::Arc;

use crate::{
    field::goldilocks::GoldilocksField,
    gadgets::{
        curves::sw_projective::SWProjectivePoint,
        non_native_field::implementations::{NonNativeFieldOverU16, NonNativeFieldOverU16Params},
        tower_extension::{
            algebraic_torus::TorusWrapper,
            params::bn256::{BN256Extension12Params, BN256Extension2Params, BN256Extension6Params},
        },
    },
    pairing::{
        bn256::{Fq12, Fq2, Fq6},
        ff::PrimeField,
    },
};
use serde::{Deserialize, Serialize};

use crate::cs::traits::cs::ConstraintSystem;

type F = GoldilocksField;
type P = GoldilocksField;
use crate::pairing::bn256::fq::Fq as BN256Fq;
// Order of group of points for bn256 curve
use crate::pairing::bn256::fr::Fr as BN256Fr;

pub use crate::pairing::bn256::G1Affine as BN256Affine;

use crate::gadgets::tower_extension::{
    fq12::Fq12 as NonNativeFq12, fq2::Fq2 as NonNativeFq2, fq6::Fq6 as NonNativeFq6,
};

/// Params of BN256 base field
pub type BN256BaseNNFieldParams = NonNativeFieldOverU16Params<BN256Fq, 17>;
/// Params of BN256 scalar field
pub type BN256ScalarNNFieldParams = NonNativeFieldOverU16Params<BN256Fr, 17>;
/// Non-native field over u16 for BN256 base field
pub type BN256BaseNNField<F> = NonNativeFieldOverU16<F, BN256Fq, 17>;
/// Non-native field over u16 for BN256 scalar field
pub type BN256ScalarNNField<F> = NonNativeFieldOverU16<F, BN256Fr, 17>;

// P.S. we used 17 bits since 17 bits * 16 bits in u16 = 272 bits > 254 bits
// used in BN254 (so we have some extra space to deal with)

// --- Field extensions for BN256 curve ---
/// Non-native field extension Fq2 for BN256 curve
pub type BN256Fq2NNField<F> = NonNativeFq2<F, BN256Fq, BN256BaseNNField<F>, BN256Extension2Params>;
/// Non-native field extension Fq6 for BN256 curve
pub type BN256Fq6NNField<F> = NonNativeFq6<F, BN256Fq, BN256BaseNNField<F>, BN256Extension6Params>;
/// Non-native field extension Fq12 for BN256 curve
pub type BN256Fq12NNField<F> =
    NonNativeFq12<F, BN256Fq, BN256BaseNNField<F>, BN256Extension12Params>;

// --- Torus compression types for BN256 curve ---
pub type BN256TorusWrapper<F> =
    TorusWrapper<F, BN256Fq, BN256BaseNNField<F>, BN256Extension12Params>;

// --- SW Projective points for BN256 curves: regular and twisted ---
/// SW Projective point for BN256 curve over non-extended base field
pub type BN256SWProjectivePoint<F> = SWProjectivePoint<F, BN256Affine, BN256BaseNNField<F>>;

/// Representation of an elliptic curve point in raw form (as strings)
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RawG1Point {
    pub x: String,
    pub y: String,
}

pub fn bn254_base_field_params() -> BN256BaseNNFieldParams {
    NonNativeFieldOverU16Params::create()
}

impl RawG1Point {
    /// Converts a raw point to a projective point
    pub fn to_projective_point<CS: ConstraintSystem<F>>(
        &self,
        cs: &mut CS,
    ) -> BN256SWProjectivePoint<F> {
        let base_params = Arc::new(bn254_base_field_params());

        let x = BN256Fq::from_str(self.x.as_str()).unwrap();
        let y = BN256Fq::from_str(self.y.as_str()).unwrap();

        let x_nn = BN256BaseNNField::allocate_checked(cs, x, &base_params);
        let y_nn = BN256BaseNNField::allocate_checked(cs, y, &base_params);

        BN256SWProjectivePoint::<F>::from_xy_unchecked(cs, x_nn, y_nn)
    }

    /// Converts a raw point to a the tuple of allocated coordinates `(x, y)`
    pub fn to_coordinates<CS: ConstraintSystem<F>>(
        &self,
        cs: &mut CS,
    ) -> (BN256BaseNNField<F>, BN256BaseNNField<F>) {
        let base_params = Arc::new(bn254_base_field_params());

        let x = BN256Fq::from_str(self.x.as_str()).unwrap();
        let y = BN256Fq::from_str(self.y.as_str()).unwrap();

        let x_nn = BN256BaseNNField::allocate_checked(cs, x, &base_params);
        let y_nn = BN256BaseNNField::allocate_checked(cs, y, &base_params);

        (x_nn, y_nn)
    }
}

/// Representation of a G2 elliptic curve point in raw form (as strings)
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RawG2Point {
    pub x: RawFq2,
    pub y: RawFq2,
}

/// Representation of an `Fq2` element in a raw form (as strings)
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RawFq2 {
    pub c0: String,
    pub c1: String,
}

impl RawFq2 {
    /// Converts a raw point to a non-native fq2 element
    pub fn to_fq2<CS: ConstraintSystem<F>>(&self, cs: &mut CS) -> BN256Fq2NNField<F> {
        let base_params = Arc::new(bn254_base_field_params());

        let c0 = BN256Fq::from_str(self.c0.as_str()).unwrap();
        let c0 = BN256BaseNNField::allocate_checked(cs, c0, &base_params);

        let c1 = BN256Fq::from_str(self.c1.as_str()).unwrap();
        let c1 = BN256BaseNNField::allocate_checked(cs, c1, &base_params);

        BN256Fq2NNField::new(c0, c1)
    }

    pub fn to_native_fq2(&self) -> Fq2 {
        let c0 = BN256Fq::from_str(self.c0.as_str()).unwrap();
        let c1 = BN256Fq::from_str(self.c1.as_str()).unwrap();

        Fq2 { c0, c1 }
    }
}

/// Representation of an `Fq6` element in a raw form (as strings)
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RawFq6 {
    pub c0: RawFq2,
    pub c1: RawFq2,
    pub c2: RawFq2,
}

impl RawFq6 {
    /// Converts a raw point to a non-native `Fq6` element
    pub fn to_fq6<CS: ConstraintSystem<F>>(&self, cs: &mut CS) -> BN256Fq6NNField<F> {
        let c0 = self.c0.to_fq2(cs);
        let c1 = self.c1.to_fq2(cs);
        let c2 = self.c2.to_fq2(cs);

        BN256Fq6NNField::new(c0, c1, c2)
    }

    pub fn to_native_fq6(&self) -> Fq6 {
        let c0 = self.c0.to_native_fq2();
        let c1 = self.c1.to_native_fq2();
        let c2 = self.c2.to_native_fq2();

        Fq6 { c0, c1, c2 }
    }
}

/// Representation of an `Fq12` element in a raw form (as strings)
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RawFq12 {
    pub c0: RawFq6,
    pub c1: RawFq6,
}

impl RawFq12 {
    /// Converts a raw point to a non-native `Fq12` element
    pub fn to_fq12<CS: ConstraintSystem<F>>(&self, cs: &mut CS) -> BN256Fq12NNField<F> {
        let c0 = self.c0.to_fq6(cs);
        let c1 = self.c1.to_fq6(cs);

        BN256Fq12NNField::new(c0, c1)
    }

    pub fn to_native_fq12(&self) -> Fq12 {
        let c0 = self.c0.to_native_fq6();
        let c1 = self.c1.to_native_fq6();

        Fq12 { c0, c1 }
    }
}
