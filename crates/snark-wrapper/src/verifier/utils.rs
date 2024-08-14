use super::*;

pub fn smart_and<E: Engine, CS: ConstraintSystem<E>>(cs: &mut CS, bools: &[Boolean]) -> Result<Boolean, SynthesisError> {
    const LIMIT: usize = 4;
    assert!(bools.len() > 0);
    if bools.len() == 1 {
        return Ok(bools[0]);
    }

    if bools.len() == 2 {
        // 1 gate
        let result = Boolean::and(cs, &bools[0], &bools[1])?;
        return Ok(result);
    }

    // 1 gate for 2,
    // 2 gates for 3, etc
    if bools.len() < LIMIT {
        // 1 gate
        let mut result = Boolean::and(cs, &bools[0], &bools[1])?;
        // len - 2 gates
        for b in bools[2..].iter() {
            result = Boolean::and(cs, &result, &b)?;
        }
        return Ok(result);
    }

    // 1 gate for 3
    // 2 gates for 6
    // 3 gates for 9, etc
    let mut lc = LinearCombination::zero();
    let num_elements_as_fr = E::Fr::from_str(&bools.len().to_string()).unwrap();
    lc.sub_assign_constant(num_elements_as_fr);
    for b in bools.iter() {
        lc.add_assign_boolean_with_coeff(b, E::Fr::one());
    }
    let as_num = lc.into_num(cs)?;

    // 2 gates here
    let all_true = as_num.is_zero(cs)?;

    // so 2 gates for 3
    // 4 gates for 6
    // 5 gates for 9
    // so we win at 4+

    Ok(all_true)
}

pub(crate) fn binary_select<E: Engine, CS: ConstraintSystem<E>>(cs: &mut CS, elements: &[GoldilocksField<E>], bits: &[Boolean]) -> Result<GoldilocksField<E>, SynthesisError> {
    assert_eq!(elements.len(), 1 << bits.len());
    assert!(bits.len() > 0);

    let mut input_space = Vec::with_capacity(elements.len() / 2);
    let mut dst_space = Vec::with_capacity(elements.len() / 2);

    for (idx, bit) in bits.iter().enumerate() {
        let src = if idx == 0 { elements } else { &input_space };

        debug_assert_eq!(elements.len() % 2, 0);
        dst_space.clear();

        for src in src.array_chunks::<2>() {
            let [a, b] = src;
            // NOTE order here
            let selected = GoldilocksField::conditionally_select(cs, *bit, b, a)?;
            dst_space.push(selected);
        }

        std::mem::swap(&mut dst_space, &mut input_space);
    }

    assert_eq!(input_space.len(), 1);

    Ok(input_space.pop().unwrap())
}

pub(crate) fn materialize_powers_serial<E: Engine, CS: ConstraintSystem<E> + 'static>(cs: &mut CS, base: GoldilocksExtAsFieldWrapper<E, CS>, size: usize) -> Vec<GoldilocksExtAsFieldWrapper<E, CS>> {
    if size == 0 {
        return Vec::new();
    }
    let mut storage = Vec::with_capacity(size);
    let mut current = GoldilocksExtAsFieldWrapper::one(cs);
    storage.push(current);
    for idx in 1..size {
        if idx == 1 {
            current = base;
        } else {
            current.mul_assign(&base, cs);
        }
        storage.push(current);
    }

    storage
}
