use crate::cs::traits::cs::ConstraintSystem;
use crate::field::goldilocks::GoldilocksField;
use crate::gadgets::tower_extension::tests::json::types::{
    BN256Fq12NNField, BN256Fq2NNField, BN256Fq6NNField,
};
use crate::gadgets::traits::witnessable::WitnessHookable;

type F = GoldilocksField;

fn equal_fq2<CS: ConstraintSystem<F>>(
    cs: &mut CS,
    a: &BN256Fq2NNField<F>,
    b: &BN256Fq2NNField<F>,
) -> bool {
    let a_c0 = a.c0.witness_hook(cs)().unwrap().get();
    let a_c1 = a.c1.witness_hook(cs)().unwrap().get();

    let b_c0 = b.c0.witness_hook(cs)().unwrap().get();
    let b_c1 = b.c1.witness_hook(cs)().unwrap().get();

    a_c0.eq(&b_c0) && a_c1.eq(&b_c1)
}

fn equal_fq6<CS: ConstraintSystem<F>>(
    cs: &mut CS,
    a: &BN256Fq6NNField<F>,
    b: &BN256Fq6NNField<F>,
) -> bool {
    equal_fq2(cs, &a.c0, &b.c0) && equal_fq2(cs, &a.c1, &b.c1) && equal_fq2(cs, &a.c2, &b.c2)
}

fn equal_fq12<CS: ConstraintSystem<F>>(
    cs: &mut CS,
    a: &BN256Fq12NNField<F>,
    b: &BN256Fq12NNField<F>,
) -> bool {
    equal_fq6(cs, &a.c0, &b.c0) && equal_fq6(cs, &a.c1, &b.c1)
}

pub(in super::super) fn assert_equal_fq2<CS: ConstraintSystem<F>>(
    cs: &mut CS,
    a: &BN256Fq2NNField<F>,
    b: &BN256Fq2NNField<F>,
) {
    assert!(equal_fq2(cs, a, b));
}

pub(in super::super) fn assert_equal_fq6<CS: ConstraintSystem<F>>(
    cs: &mut CS,
    a: &BN256Fq6NNField<F>,
    b: &BN256Fq6NNField<F>,
) {
    assert!(equal_fq6(cs, a, b));
}

pub(in super::super) fn assert_equal_fq12<CS: ConstraintSystem<F>>(
    cs: &mut CS,
    a: &BN256Fq12NNField<F>,
    b: &BN256Fq12NNField<F>,
) {
    assert!(equal_fq12(cs, a, b));
}

pub(in super::super) fn assert_not_equal_fq12<CS: ConstraintSystem<F>>(
    cs: &mut CS,
    a: &BN256Fq12NNField<F>,
    b: &BN256Fq12NNField<F>,
) {
    assert!(!equal_fq12(cs, a, b));
}
