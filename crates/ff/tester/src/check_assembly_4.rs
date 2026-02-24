#[cfg(test)]
mod test {
    use super::super::mul_variant0::Fs;
    // use super::super::mul_variant0::mont_mul_asm;
    use super::super::assembly_4::*;

    use crate::rand::*;
    use ff::*;

    #[test]
    fn check_mul_asm() {
        let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        for i in 0..10000 {
            let a: Fs = rng.gen();
            let b: Fs = rng.gen();

            let a_asm = unsafe { std::mem::transmute::<_, [u64; 4]>(a) };
            let b_asm = unsafe { std::mem::transmute::<_, [u64; 4]>(b) };

            let mut c = a;
            c.mul_assign(&b);

            // let c_asm = mont_mul_asm(&a_asm, &b_asm);
            // let c_asm = mont_mul_asm_adx(&a_asm, &b_asm);
            let c_asm = mont_mul_asm_adx_with_reduction(&a_asm, &b_asm);

            let mut c_back = unsafe { std::mem::transmute::<_, Fs>(c_asm) };
            // if !c_back.is_valid() {
            //     c_back.reduce();
            // }

            assert_eq!(c, c_back, "failed at iteration {}: a = {:?}, b = {:?}", i, a, b);
        }
    }
}
