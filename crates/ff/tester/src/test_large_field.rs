use ff::*;

#[derive(PrimeField)]
#[PrimeFieldModulus = "21888242871839275222246405745257275088696311157297823662689037894645226208583"]
#[PrimeFieldGenerator = "2"]
pub(crate) struct Fr(FrRepr);

#[test]
fn test_to_hex() {
    let a = Fr::from_repr(FrRepr::from(2)).unwrap();
    assert_eq!("0000000000000000000000000000000000000000000000000000000000000002", to_hex(&a));
    println!("`2` into hex = {}", to_hex(&a));
}

#[test]
fn test_hash_impl() {
    let mut hashset = std::collections::HashSet::new();

    hashset.insert(Fr::from_repr(FrRepr::from(2)).unwrap());
    hashset.insert(Fr::from_repr(FrRepr::from(4)).unwrap());
    hashset.insert(Fr::from_repr(FrRepr::from(2)).unwrap());

    assert_eq!(hashset.len(), 2);
}
