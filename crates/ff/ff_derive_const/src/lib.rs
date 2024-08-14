
#![feature(const_generics)]
#![feature(const_generic_impls_guard)]

#![recursion_limit = "1024"]

extern crate ff;
extern crate rand;

pub mod const_repr;
pub mod const_field_element;
// mod alt;

#[cfg(test)]
mod tests {
    use super::const_repr::*;
    use super::const_field_element::*;

    const MODULUS: BigintRepresentation::<4> = BigintRepresentation::<4>([
        0x677297dc392126f1,
        0xab3eedb83920ee0a,
        0x370a08b6d0302b0b,
        0x060c89ce5c263405,
    ]);

    const MODULUS_BITS: u32 = 251;
    const REPR_SHAVE_BITS: u32 = 5;
    const R: BigintRepresentation::<4> = BigintRepresentation::<4>([
        0x073315dea08f9c76,
        0xe7acffc6a098f24b,
        0xf85a9201d818f015,
        0x01f16424e1bb7724,
    ]);

    const R2: BigintRepresentation::<4> = BigintRepresentation::<4>([
        0x35e44abee7ecb21e,
        0x74646cacf5f84ec4,
        0xe472df203faa158f,
        0x0445b524f1ba50a8,
    ]);

    const INV: u64 = 0x532ce5aebc48f5ef;

    const GENERATOR: BigintRepresentation::<4> = BigintRepresentation::<4>([
        0x6380695df1aaf958,
        0xff3d22fdf1ecc3f8,
        0x5c65ec9f484e3a81,
        0x0180a96573d3d9f8,
    ]);

    const S: u32 = 4;

    const ROOT_OF_UNITY: BigintRepresentation::<4> = BigintRepresentation::<4>([
        0xa13885692e7afcb0,
        0xb789766cd18573ca,
        0xd5468c0174efc3b9,
        0x03534b612b0b6f7a,
    ]);

    #[derive(Copy, Clone)]
    struct FsParams;

    impl FieldParameters<4> for FsParams {
        const NUM_BITS: u32 = MODULUS_BITS;
        const CAPACITY: u32 = 250;
        const REPR_SHAVE_BITS: u32 = REPR_SHAVE_BITS;
        const S: u32 = S;
        const MULTIPLICATIVE_GENERATOR: BigintRepresentation::<4> = GENERATOR;
        const ROOT_OF_UNITY: BigintRepresentation::<4> = ROOT_OF_UNITY;
        const MODULUS: BigintRepresentation::<4> = MODULUS;
        const R: BigintRepresentation::<4> = R;
        const R2: BigintRepresentation::<4> = R2;
        const INV: u64 = INV;
    }

    type Fs = PrimeFieldElement<FsParams, 4>;

    #[test]
    fn make_naive() {
        use crate::ff::PrimeField;

        let repr = BigintRepresentation::<4>::from(3u64);
        
        let fe = Fs::from_repr(repr).unwrap();
        println!("{:?}", fe);
    }
}