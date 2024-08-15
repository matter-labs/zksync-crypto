mod circuits;
mod generate;
mod serialize;

use std::str::FromStr;

use circuits::MockCircuit;
pub use generate::{generate, Encoding, MainGateType};

use ethereum_types::U256;
use rescue_poseidon::franklin_crypto::bellman::pairing::bn256::{Bn256, Fr};
use rescue_poseidon::franklin_crypto::bellman::pairing::ff::*;
use rescue_poseidon::franklin_crypto::bellman::pairing::*;
use rescue_poseidon::franklin_crypto::bellman::plonk::better_better_cs::cs::{Circuit, VerificationKey};
use rescue_poseidon::franklin_crypto::bellman::plonk::better_better_cs::proof::Proof;

fn render_scalar_to_hex<F: PrimeField>(el: &F) -> String {
    let mut buff = vec![];
    let repr = el.into_repr();
    repr.write_be(&mut buff).unwrap();

    format!("0x{}", hex::encode(buff))
}

fn render_g1_affine_to_hex<E: Engine>(point: &E::G1Affine) -> [String; 2] {
    if <E::G1Affine as CurveAffine>::is_zero(point) {
        return ["0x0".to_owned(), "0x0".to_owned()];
    }

    let (x, y) = <E::G1Affine as CurveAffine>::into_xy_unchecked(*point);
    [render_scalar_to_hex(&x), render_scalar_to_hex(&y)]
}

fn render_g2_affine_to_hex(point: &<bn256::Bn256 as Engine>::G2Affine) -> [String; 4] {
    if <<bn256::Bn256 as Engine>::G2Affine as CurveAffine>::is_zero(point) {
        return ["0x0".to_owned(), "0x0".to_owned(), "0x0".to_owned(), "0x0".to_owned()];
    }

    let (x, y) = <<bn256::Bn256 as Engine>::G2Affine as CurveAffine>::into_xy_unchecked(*point);

    [render_scalar_to_hex(&x.c0), render_scalar_to_hex(&x.c1), render_scalar_to_hex(&y.c0), render_scalar_to_hex(&y.c1)]
}

fn serialize_g1_for_ethereum(point: &<bn256::Bn256 as Engine>::G1Affine) -> (U256, U256) {
    if <<bn256::Bn256 as Engine>::G1Affine as CurveAffine>::is_zero(point) {
        return (U256::zero(), U256::zero());
    }
    let uncompressed = point.into_uncompressed();

    let (x, y) = <<bn256::Bn256 as Engine>::G1Affine as CurveAffine>::into_xy_unchecked(*point);
    let _ = <<bn256::Bn256 as Engine>::G1Affine as CurveAffine>::from_xy_checked(x, y).unwrap();

    let mut buffer = [0u8; 32];
    x.into_repr().write_be(&mut buffer[..]).unwrap();
    let x = U256::from_big_endian(&buffer);

    let mut buffer = [0u8; 32];
    y.into_repr().write_be(&mut buffer[..]).unwrap();
    let y = U256::from_big_endian(&buffer);

    (x, y)
}

fn serialize_fe_for_ethereum(field_element: &Fr) -> U256 {
    let mut be_bytes = [0u8; 32];
    field_element.into_repr().write_be(&mut be_bytes[..]).expect("get new root BE bytes");
    U256::from_big_endian(&be_bytes[..])
}

pub fn serialize_proof<T: Circuit<Bn256>>(proof: &Proof<Bn256, T>) -> (Vec<U256>, Vec<U256>) {
    let mut inputs = vec![];
    for input in proof.inputs.iter() {
        inputs.push(serialize_fe_for_ethereum(&input));
    }
    let mut serialized_proof = vec![];

    for c in proof.state_polys_commitments.iter() {
        let (x, y) = serialize_g1_for_ethereum(&c);
        serialized_proof.push(x);
        serialized_proof.push(y);
    }

    let (x, y) = serialize_g1_for_ethereum(&proof.copy_permutation_grand_product_commitment);
    serialized_proof.push(x);
    serialized_proof.push(y);

    let (x, y) = serialize_g1_for_ethereum(&proof.lookup_s_poly_commitment.unwrap());
    serialized_proof.push(x);
    serialized_proof.push(y);

    let (x, y) = serialize_g1_for_ethereum(&proof.lookup_grand_product_commitment.unwrap());
    serialized_proof.push(x);
    serialized_proof.push(y);

    for c in proof.quotient_poly_parts_commitments.iter() {
        let (x, y) = serialize_g1_for_ethereum(&c);
        serialized_proof.push(x);
        serialized_proof.push(y);
    }

    assert_eq!(proof.state_polys_openings_at_z.len(), 4);
    for c in proof.state_polys_openings_at_z.iter() {
        serialized_proof.push(serialize_fe_for_ethereum(&c));
    }

    assert_eq!(proof.state_polys_openings_at_dilations.len(), 1);
    for (_, _, c) in proof.state_polys_openings_at_dilations.iter() {
        serialized_proof.push(serialize_fe_for_ethereum(&c));
    }

    assert_eq!(proof.gate_setup_openings_at_z.len(), 0);
    for (_, _, c) in proof.gate_setup_openings_at_z.iter() {
        serialized_proof.push(serialize_fe_for_ethereum(&c));
    }

    assert_eq!(proof.gate_selectors_openings_at_z.len(), 1);
    for (_, c) in proof.gate_selectors_openings_at_z.iter() {
        serialized_proof.push(serialize_fe_for_ethereum(&c));
    }

    for c in proof.copy_permutation_polys_openings_at_z.iter() {
        serialized_proof.push(serialize_fe_for_ethereum(&c));
    }

    serialized_proof.push(serialize_fe_for_ethereum(&proof.copy_permutation_grand_product_opening_at_z_omega));
    serialized_proof.push(serialize_fe_for_ethereum(&proof.lookup_s_poly_opening_at_z_omega.unwrap()));
    serialized_proof.push(serialize_fe_for_ethereum(&proof.lookup_grand_product_opening_at_z_omega.unwrap()));
    serialized_proof.push(serialize_fe_for_ethereum(&proof.lookup_t_poly_opening_at_z.unwrap()));
    serialized_proof.push(serialize_fe_for_ethereum(&proof.lookup_t_poly_opening_at_z_omega.unwrap()));
    serialized_proof.push(serialize_fe_for_ethereum(&proof.lookup_selector_poly_opening_at_z.unwrap()));
    serialized_proof.push(serialize_fe_for_ethereum(&proof.lookup_table_type_poly_opening_at_z.unwrap()));
    serialized_proof.push(serialize_fe_for_ethereum(&proof.quotient_poly_opening_at_z));
    serialized_proof.push(serialize_fe_for_ethereum(&proof.linearization_poly_opening_at_z));

    let (x, y) = serialize_g1_for_ethereum(&proof.opening_proof_at_z);
    serialized_proof.push(x);
    serialized_proof.push(y);

    let (x, y) = serialize_g1_for_ethereum(&proof.opening_proof_at_z_omega);
    serialized_proof.push(x);
    serialized_proof.push(y);

    dbg!(&serialized_proof.len());

    (inputs, serialized_proof)
}
