use franklin_crypto::bellman::{CurveAffine, Engine, PrimeField, PrimeFieldRepr, pairing::bn256::{Bn256, Fq, Fq2, Fr}};


use serde::ser::{Serialize, SerializeStruct, Serializer};

pub struct FieldElement(Fr);

impl From<Fr> for FieldElement {
    fn from(el: Fr) -> Self {
        Self(el)
    }
}

impl Serialize for FieldElement {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("FieldElement", 1)?;

        let mut init_str = String::from("PairingsBn254.new_fr");
        init_str.push_str("(");
        init_str.push_str(&field_element_to_hex(&self.0));
        init_str.push_str(")");
        state.serialize_field("el", &init_str)?;
        state.end()
    }
}

fn field_element_to_hex<F: PrimeField>(el: &F) -> String {
    let mut buf = [0u8; 32];
    el.into_repr()
        .write_be(&mut buf[..])
        .expect("consume buffer");

    let mut result = String::from("0x");
    result.push_str(&hex::encode(buf));

    result
}

pub struct G1Point {
    x: Fq,
    y: Fq,
}

impl G1Point {
    pub(crate) fn from_affine_point(point: <Bn256 as Engine>::G1Affine) -> Self {
        let (x, y) = point.into_xy_unchecked();
        Self { x, y }
    }
}

impl Serialize for G1Point {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("G1Point", 2)?;

        let mut init_str = String::from("PairingsBn254.new_g1");
        init_str.push_str("(");
        init_str.push_str(&field_element_to_hex(&self.x));
        init_str.push_str(",");
        init_str.push_str(&field_element_to_hex(&self.y));
        init_str.push_str(")");
        state.serialize_field("g1", &init_str)?;
        state.end()
    }
}

pub struct G2Point {
    x: Fq2,
    y: Fq2,
}

impl G2Point {
    pub(crate) fn from_affine_point(point: <Bn256 as Engine>::G2Affine) -> Self {
        let (x, y) = point.into_xy_unchecked();
        Self { x, y }
    }
}

impl Serialize for G2Point {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("G2Point", 2)?;
        let mut init_str = String::from("PairingsBn254.new_g2");
        // Encoding of field elements is: X[0] * z + X[1]
        init_str.push_str("([");
        init_str.push_str(&field_element_to_hex(&self.x.c1));
        init_str.push_str(",");
        init_str.push_str(&field_element_to_hex(&self.x.c0));
        init_str.push_str("],[");
        init_str.push_str(&field_element_to_hex(&self.y.c1));
        init_str.push_str(",");
        init_str.push_str(&field_element_to_hex(&self.y.c0));
        init_str.push_str("])");
        state.serialize_field("g2", &init_str)?;
        state.end()
    }
}
