mod normal {
    use ff::*;

    #[derive(PrimeField)]
    #[PrimeFieldModulus = "21888242871839275222246405745257275088696311157297823662689037894645226208583"]
    #[PrimeFieldGenerator = "2"]
    pub struct Fr(FrRepr);
}

mod asm {
    use ff::*;

    #[derive(PrimeFieldAsm)]
    #[PrimeFieldModulus = "21888242871839275222246405745257275088696311157297823662689037894645226208583"]
    #[PrimeFieldGenerator = "2"]
    #[UseADX = "true"]
    pub struct FrAsm(FrReprAsm);
}

pub use self::asm::FrAsm;
pub use self::normal::Fr;

#[cfg(test)]
mod test {
    use super::*;

    use crate::rand::*;
    use ff::*;

    #[test]
    fn check_mul_asm() {
        let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        for i in 0..10000 {
            let a: Fr = rng.gen();
            let b: Fr = rng.gen();

            let a_asm = unsafe { std::mem::transmute::<_, FrAsm>(a) };
            let b_asm = unsafe { std::mem::transmute::<_, FrAsm>(b) };

            let mut c = a;
            c.mul_assign(&b);

            let mut c_asm = a_asm;
            c_asm.mul_assign(&b_asm);

            let c_back = unsafe { std::mem::transmute::<_, Fr>(c_asm) };

            assert_eq!(c, c_back, "failed at iteration {}: a = {:?}, b = {:?}", i, a, b);
        }
    }

    #[test]
    fn check_add_asm() {
        let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        for i in 0..10000 {
            let a: Fr = rng.gen();
            let b: Fr = rng.gen();

            let a_asm = unsafe { std::mem::transmute::<_, FrAsm>(a) };
            let b_asm = unsafe { std::mem::transmute::<_, FrAsm>(b) };

            let mut c = a;
            c.add_assign(&b);

            let mut c_asm = a_asm;
            c_asm.add_assign(&b_asm);

            let c_back = unsafe { std::mem::transmute::<_, Fr>(c_asm) };

            assert_eq!(c, c_back, "failed at iteration {}: a = {:?}, b = {:?}", i, a, b);
        }
    }

    #[test]
    fn check_sub_asm() {
        let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        for i in 0..10000 {
            let a: Fr = rng.gen();
            let b: Fr = rng.gen();

            let a_asm = unsafe { std::mem::transmute::<_, FrAsm>(a) };
            let b_asm = unsafe { std::mem::transmute::<_, FrAsm>(b) };

            let mut c = a;
            c.sub_assign(&b);

            let mut c_asm = a_asm;
            c_asm.sub_assign(&b_asm);

            let c_back = unsafe { std::mem::transmute::<_, Fr>(c_asm) };

            assert_eq!(c, c_back, "failed at iteration {}: a = {:?}, b = {:?}", i, a, b);
        }
    }

    #[test]
    fn check_double_asm() {
        let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        for i in 0..10000 {
            let a: Fr = rng.gen();

            let a_asm = unsafe { std::mem::transmute::<_, FrAsm>(a) };

            let mut c = a;
            c.double();

            let mut c_asm = a_asm;
            c_asm.double();

            let c_back = unsafe { std::mem::transmute::<_, Fr>(c_asm) };

            assert_eq!(c, c_back, "failed at iteration {}: a = {:?}", i, a);
        }
    }

    #[test]
    fn check_square_asm() {
        let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        for i in 0..10000 {
            let a: Fr = rng.gen();

            let a_asm = unsafe { std::mem::transmute::<_, FrAsm>(a) };

            let mut c = a;
            c.square();

            let mut c_asm = a_asm;
            c_asm.square();

            let c_back = unsafe { std::mem::transmute::<_, Fr>(c_asm) };

            assert_eq!(c, c_back, "failed at iteration {}: a = {:?}", i, a);
        }
    }
}
