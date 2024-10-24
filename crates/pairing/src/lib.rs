#![allow(clippy::too_many_arguments, clippy::needless_borrows_for_generic_args)]
// Force public structures to implement Debug
#![deny(missing_debug_implementations)]
// Asm is only available on nightly, with this unstable feature
#![cfg_attr(feature = "asm", feature(asm_const))]

extern crate byteorder;
extern crate rand;

#[cfg(test)]
pub mod tests;

pub extern crate ff;

pub use ff::*;

pub mod bls12_381;
pub mod bn256;
pub mod compact_bn256;

mod wnaf;
pub use self::wnaf::Wnaf;

mod base;
pub use self::base::*;

use ff::{Field, PrimeField, PrimeFieldDecodingError, PrimeFieldRepr, ScalarEngine, SqrtField};
use std::error::Error;
use std::fmt;

/// An "engine" is a collection of types (fields, elliptic curve groups, etc.)
/// with well-defined relationships. In particular, the G1/G2 curve groups are
/// of prime order `r`, and are equipped with a bilinear pairing function.
pub trait Engine: ScalarEngine {
    /// The projective representation of an element in G1.
    type G1: CurveProjective<Engine = Self, Base = Self::Fq, Scalar = Self::Fr, Affine = Self::G1Affine> + From<Self::G1Affine>;

    /// The affine representation of an element in G1.
    type G1Affine: CurveAffine<Engine = Self, Base = Self::Fq, Scalar = Self::Fr, Projective = Self::G1, Pair = Self::G2Affine, PairingResult = Self::Fqk> + From<Self::G1> + RawEncodable;

    /// The projective representation of an element in G2.
    type G2: CurveProjective<Engine = Self, Base = Self::Fqe, Scalar = Self::Fr, Affine = Self::G2Affine> + From<Self::G2Affine>;

    /// The affine representation of an element in G2.
    type G2Affine: CurveAffine<Engine = Self, Base = Self::Fqe, Scalar = Self::Fr, Projective = Self::G2, Pair = Self::G1Affine, PairingResult = Self::Fqk> + From<Self::G2>;

    /// The base field that hosts G1.
    type Fq: PrimeField + SqrtField;

    /// The extension field that hosts G2.
    type Fqe: SqrtField;

    /// The extension field that hosts the target group of the pairing.
    type Fqk: Field;

    /// Perform a miller loop with some number of (G1, G2) pairs.
    fn miller_loop<'a, I>(i: I) -> Self::Fqk
    where
        I: IntoIterator<Item = &'a (&'a <Self::G1Affine as CurveAffine>::Prepared, &'a <Self::G2Affine as CurveAffine>::Prepared)>;

    /// Perform final exponentiation of the result of a miller loop.
    fn final_exponentiation(r: &Self::Fqk) -> Option<Self::Fqk>;

    /// Performs a complete pairing operation `(p, q)`.
    fn pairing<G1, G2>(p: G1, q: G2) -> Self::Fqk
    where
        G1: Into<Self::G1Affine>,
        G2: Into<Self::G2Affine>,
    {
        Self::final_exponentiation(&Self::miller_loop([(&(p.into().prepare()), &(q.into().prepare()))].iter())).unwrap()
    }
}

/// Projective representation of an elliptic curve point guaranteed to be
/// in the correct prime order subgroup.
pub trait CurveProjective: PartialEq + Eq + Sized + Copy + Clone + Send + Sync + fmt::Debug + fmt::Display + rand::Rand + 'static {
    type Engine: Engine<Fr = Self::Scalar>;
    type Scalar: PrimeField + SqrtField;
    type Base: SqrtField;
    type Affine: CurveAffine<Projective = Self, Scalar = Self::Scalar, Base = Self::Base>;

    /// Returns the additive identity.
    fn zero() -> Self;

    /// Returns a fixed generator of unknown exponent.
    fn one() -> Self;

    /// Determines if this point is the point at infinity.
    fn is_zero(&self) -> bool;

    /// Normalizes a slice of projective elements so that
    /// conversion to affine is cheap.
    fn batch_normalization(v: &mut [Self]);

    /// Checks if the point is already "normalized" so that
    /// cheap affine conversion is possible.
    fn is_normalized(&self) -> bool;

    /// Doubles this element.
    fn double(&mut self);

    /// Adds another element to this element.
    fn add_assign(&mut self, other: &Self);

    /// Subtracts another element from this element.
    fn sub_assign(&mut self, other: &Self) {
        let mut tmp = *other;
        tmp.negate();
        self.add_assign(&tmp);
    }

    /// Adds an affine element to this element.
    fn add_assign_mixed(&mut self, other: &Self::Affine);

    /// Negates this element.
    fn negate(&mut self);

    /// Performs scalar multiplication of this element.
    fn mul_assign<S: Into<<Self::Scalar as PrimeField>::Repr>>(&mut self, other: S);

    /// Converts this element into its affine representation.
    fn into_affine(&self) -> Self::Affine;

    /// Recommends a wNAF window table size given a scalar. Always returns a number
    /// between 2 and 22, inclusive.
    fn recommended_wnaf_for_scalar(scalar: <Self::Scalar as PrimeField>::Repr) -> usize;

    /// Recommends a wNAF window size given the number of scalars you intend to multiply
    /// a base by. Always returns a number between 2 and 22, inclusive.
    fn recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize;

    /// Returns references to underlying X, Y and Z coordinates. Users should check for infinity
    /// outside of this call
    fn as_xyz(&self) -> (&Self::Base, &Self::Base, &Self::Base) {
        unimplemented!("default implementation does not exist for this function")
    }

    /// Returns underlying X, Y and Z coordinates. Users should check for infinity
    /// outside of this call
    fn into_xyz_unchecked(self) -> (Self::Base, Self::Base, Self::Base) {
        unimplemented!("default implementation does not exist for this function")
    }

    /// Creates a point from raw X, Y and Z coordinates. Point of infinity is encoded as (0,1,0) by default.
    /// On-curve check is NOT performed
    fn from_xyz_unchecked(_x: Self::Base, _y: Self::Base, _z: Self::Base) -> Self {
        unimplemented!("default implementation does not exist for this function")
    }

    /// Creates a point from raw X, Y and Z coordinates. Point of infinity is encoded as (0,1,0) by default.
    /// On-curve check is performed
    fn from_xyz_checked(_x: Self::Base, _y: Self::Base, _z: Self::Base) -> Result<Self, GroupDecodingError> {
        unimplemented!("default implementation does not exist for this function")
    }
}

/// Affine representation of an elliptic curve point guaranteed to be
/// in the correct prime order subgroup.
pub trait CurveAffine: Copy + Clone + Sized + Send + Sync + fmt::Debug + fmt::Display + PartialEq + Eq + 'static + serde::Serialize + serde::de::DeserializeOwned {
    type Engine: Engine<Fr = Self::Scalar>;
    type Scalar: PrimeField + SqrtField;
    type Base: SqrtField;
    type Projective: CurveProjective<Affine = Self, Scalar = Self::Scalar, Base = Self::Base>;
    type Prepared: Clone + Send + Sync + 'static;
    type Uncompressed: EncodedPoint<Affine = Self>;
    type Compressed: EncodedPoint<Affine = Self>;
    type Pair: CurveAffine<Pair = Self>;
    type PairingResult: Field;

    /// Returns the additive identity.
    fn zero() -> Self;

    /// Returns a fixed generator of unknown exponent.
    fn one() -> Self;

    /// Determines if this point represents the point at infinity; the
    /// additive identity.
    fn is_zero(&self) -> bool;

    /// Negates this element.
    fn negate(&mut self);

    /// Performs scalar multiplication of this element with mixed addition.
    fn mul<S: Into<<Self::Scalar as PrimeField>::Repr>>(&self, other: S) -> Self::Projective;

    /// Prepares this element for pairing purposes.
    fn prepare(&self) -> Self::Prepared;

    /// Perform a pairing
    fn pairing_with(&self, other: &Self::Pair) -> Self::PairingResult;

    /// Converts this element into its affine representation.
    fn into_projective(&self) -> Self::Projective;

    /// Converts this element into its compressed encoding, so long as it's not
    /// the point at infinity.
    fn into_compressed(&self) -> Self::Compressed {
        <Self::Compressed as EncodedPoint>::from_affine(*self)
    }

    /// Converts this element into its uncompressed encoding, so long as it's not
    /// the point at infinity.
    fn into_uncompressed(&self) -> Self::Uncompressed {
        <Self::Uncompressed as EncodedPoint>::from_affine(*self)
    }

    /// Returns references to underlying X and Y coordinates. Users should check for infinity
    /// outside of this call
    fn as_xy(&self) -> (&Self::Base, &Self::Base);

    /// Returns underlying X and Y coordinates. Users should check for infinity
    /// outside of this call
    fn into_xy_unchecked(self) -> (Self::Base, Self::Base);

    /// Creates a point from raw X and Y coordinates. Point of infinity is encoded as (0,0) by default.
    /// On-curve check is NOT performed
    fn from_xy_unchecked(x: Self::Base, y: Self::Base) -> Self;

    /// Creates a point from raw X and Y coordinates. Point of infinity is encoded as (0,0) by default.
    /// On-curve check is performed
    fn from_xy_checked(x: Self::Base, y: Self::Base) -> Result<Self, GroupDecodingError>;

    /// returns A coefficient for a short Weierstrass form
    fn a_coeff() -> Self::Base;

    /// returns B coefficient for a short Weierstrass form
    fn b_coeff() -> Self::Base;
}

pub trait RawEncodable: CurveAffine {
    /// Converts this element into its uncompressed encoding, so long as it's not
    /// the point at infinity. Leaves coordinates in Montgommery form
    fn into_raw_uncompressed_le(&self) -> Self::Uncompressed;

    /// Creates a point from raw encoded coordinates without checking on curve
    fn from_raw_uncompressed_le_unchecked(encoded: &Self::Uncompressed, infinity: bool) -> Result<Self, GroupDecodingError>;

    /// Creates a point from raw encoded coordinates
    fn from_raw_uncompressed_le(encoded: &Self::Uncompressed, infinity: bool) -> Result<Self, GroupDecodingError>;
}

/// An encoded elliptic curve point, which should essentially wrap a `[u8; N]`.
pub trait EncodedPoint: Sized + Send + Sync + AsRef<[u8]> + AsMut<[u8]> + Clone + Copy + 'static {
    type Affine: CurveAffine;

    /// Creates an empty representation.
    fn empty() -> Self;

    /// Returns the number of bytes consumed by this representation.
    fn size() -> usize;

    /// Converts an `EncodedPoint` into a `CurveAffine` element,
    /// if the encoding represents a valid element.
    fn into_affine(&self) -> Result<Self::Affine, GroupDecodingError>;

    /// Converts an `EncodedPoint` into a `CurveAffine` element,
    /// without guaranteeing that the encoding represents a valid
    /// element. This is useful when the caller knows the encoding is
    /// valid already.
    ///
    /// If the encoding is invalid, this can break API invariants,
    /// so caution is strongly encouraged.
    fn into_affine_unchecked(&self) -> Result<Self::Affine, GroupDecodingError>;

    /// Creates an `EncodedPoint` from an affine point, as long as the
    /// point is not the point at infinity.
    fn from_affine(affine: Self::Affine) -> Self;
}

/// An error that may occur when trying to decode an `EncodedPoint`.
#[derive(Debug)]
pub enum GroupDecodingError {
    /// The coordinate(s) do not lie on the curve.
    NotOnCurve,
    /// The element is not part of the r-order subgroup.
    NotInSubgroup,
    /// One of the coordinates could not be decoded
    CoordinateDecodingError(&'static str, PrimeFieldDecodingError),
    /// The compression mode of the encoded element was not as expected
    UnexpectedCompressionMode,
    /// The encoding contained bits that should not have been set
    UnexpectedInformation,
}

impl GroupDecodingError {
    fn self_description(&self) -> &str {
        match *self {
            GroupDecodingError::NotOnCurve => "coordinate(s) do not lie on the curve",
            GroupDecodingError::NotInSubgroup => "the element is not part of an r-order subgroup",
            GroupDecodingError::CoordinateDecodingError(..) => "coordinate(s) could not be decoded",
            GroupDecodingError::UnexpectedCompressionMode => "encoding has unexpected compression mode",
            GroupDecodingError::UnexpectedInformation => "encoding has unexpected information",
        }
    }
}

impl Error for GroupDecodingError {
    fn description(&self) -> &str {
        self.self_description()
    }
}

impl fmt::Display for GroupDecodingError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        match *self {
            GroupDecodingError::CoordinateDecodingError(description, ref err) => {
                write!(f, "{} decoding error: {}", description, err)
            }
            _ => write!(f, "{}", self.self_description()),
        }
    }
}
