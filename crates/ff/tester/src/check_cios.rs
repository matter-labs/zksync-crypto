#[cfg(test)]
mod test {
    use super::super::test_large_field::Fr as Fr;
    use super::super::test_large_cios_field::Fr as FrCios;

    use rand::*;
    use ff::*;

    #[test]
    fn check_mul() {
        let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        for _ in 0..10000 {
            let a: Fr = rng.gen(); 
            let b: Fr = rng.gen();

            let a_cios = unsafe { std::mem::transmute::<_, FrCios>(a) };
            let b_cios = unsafe { std::mem::transmute::<_, FrCios>(b) };

            let mut c = a;
            c.mul_assign(&b);

            let mut c_cios = a_cios;
            c_cios.mul_assign(&b_cios);

            let c_back = unsafe { std::mem::transmute::<_, Fr>(c_cios) };

            assert_eq!(c, c_back);
        }
    }

    #[test]
    fn check_sqr() {
        let rng = &mut XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

        for _ in 0..10000 {
            let a: Fr = rng.gen(); 

            let a_cios = unsafe { std::mem::transmute::<_, FrCios>(a) };

            let mut c = a;
            c.square();

            let mut c_cios = a_cios;
            c_cios.square();

            let c_back = unsafe { std::mem::transmute::<_, Fr>(c_cios) };

            assert_eq!(c, c_back);
        }
    }
}