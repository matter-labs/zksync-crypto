use super::fq2::Fq2;
use super::fq6::Fq6;
use super::fq12::Fq12;
use ff::{Field, PrimeField, PrimeFieldRepr};

cfg_if::cfg_if! {
    if #[cfg(feature = "asm")] {
        use core::arch::asm;
        use ff::PrimeFieldAsm;

        #[derive(PrimeFieldAsm)]
        #[PrimeFieldModulus = "21888242871839275222246405745257275088696311157297823662689037894645226208583"]
        #[PrimeFieldGenerator = "2"]
        #[UseADX = "true"]
        pub struct Fq(FqRepr);
    } else {
        #[derive(PrimeField)]
        #[PrimeFieldModulus = "21888242871839275222246405745257275088696311157297823662689037894645226208583"]
        #[PrimeFieldGenerator = "2"]
        pub struct Fq(FqRepr);
    }
}

// B coefficient of BN256 curve, B = 3
// In Montgommery form with R = 2^256
pub const B_COEFF: Fq = Fq(FqRepr([0x7a17caa950ad28d7, 0x1f6ac17ae15521b9, 0x334bea4e696bd284, 0x2a1f6744ce179d8e]));

pub const B_COEFF_FQ2: Fq2 = Fq2 {
    c0: Fq(FqRepr([0x3bf938e377b802a8, 0x020b1b273633535d, 0x26b7edf049755260, 0x2514c6324384a86d])),
    c1: Fq(FqRepr([0x38e7ecccd1dcff67, 0x65f0b37d93ce0d3e, 0xd749d0dd22ac00aa, 0x0141b9ce4a688d4d])),
};

// The generators of G1/G2

// Generator of G1
// x = 1
// y = 2
pub const G1_GENERATOR_X: Fq = Fq(FqRepr([0xd35d438dc58f0d9d, 0x0a78eb28f5c70b3d, 0x666ea36f7879462c, 0x0e0a77c19a07df2f]));
pub const G1_GENERATOR_Y: Fq = Fq(FqRepr([0xa6ba871b8b1e1b3a, 0x14f1d651eb8e167b, 0xccdd46def0f28c58, 0x1c14ef83340fbe5e]));

// Generator of G2
//
// x = 11559732032986387107991004021392285783925812861821192530917403151452391805634*u
//     + 10857046999023057135944570762232829481370756359578518086990519993285655852781
//
// y = 4082367875863433681332203403145435568316851327593401208105741076214120093531*u
//     + 8495653923123431417604973247489272438418190587263600148770280649306958101930

pub const G2_GENERATOR_X_C0: Fq = Fq(FqRepr([0x8e83b5d102bc2026, 0xdceb1935497b0172, 0xfbb8264797811adf, 0x19573841af96503b]));
pub const G2_GENERATOR_X_C1: Fq = Fq(FqRepr([0xafb4737da84c6140, 0x6043dd5a5802d8c4, 0x09e950fc52a02f86, 0x14fef0833aea7b6b]));
pub const G2_GENERATOR_Y_C0: Fq = Fq(FqRepr([0x619dfa9d886be9f6, 0xfe7fd297f59e9b78, 0xff9e1a62231b7dfe, 0x28fd7eebae9e4206]));
pub const G2_GENERATOR_Y_C1: Fq = Fq(FqRepr([0x64095b56c71856ee, 0xdc57f922327d3cbb, 0x55f935be33351076, 0x0da4a0e693fd6482]));

// Coefficients for the Frobenius automorphism.
pub const FROBENIUS_COEFF_FQ2_C1: [Fq; 2] = [
    // Fq(-1)**(((q^0) - 1) / 2)
    // it's 1 in Montgommery form
    Fq(FqRepr([0xd35d438dc58f0d9d, 0x0a78eb28f5c70b3d, 0x666ea36f7879462c, 0x0e0a77c19a07df2f])),
    // Fq(-1)**(((q^1) - 1) / 2)
    Fq(FqRepr([0x68c3488912edefaa, 0x8d087f6872aabf4f, 0x51e1a24709081231, 0x2259d6b14729c0fa])),
];

// Fq2(u + 9)**(((q^1) - 1) / 2)
pub const XI_TO_Q_MINUS_1_OVER_2: Fq2 = Fq2 {
    c0: Fq(FqRepr([0xe4bbdd0c2936b629, 0xbb30f162e133bacb, 0x31a9d1b6f9645366, 0x253570bea500f8dd])),
    c1: Fq(FqRepr([0xa1d77ce45ffe77c7, 0x07affd117826d1db, 0x6d16bd27bb7edc6b, 0x2c87200285defecc])),
};

pub const FROBENIUS_COEFF_FQ6_C1: [Fq2; 6] = [
    // Fq2(u + 9)**(((q^0) - 1) / 3)
    Fq2 {
        c0: Fq(FqRepr([0xd35d438dc58f0d9d, 0x0a78eb28f5c70b3d, 0x666ea36f7879462c, 0x0e0a77c19a07df2f])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0])),
    },
    // Fq2(u + 9)**(((q^1) - 1) / 3)
    // taken from go-ethereum and also re-calculated manually
    Fq2 {
        c0: Fq(FqRepr([0xb5773b104563ab30, 0x347f91c8a9aa6454, 0x7a007127242e0991, 0x1956bcd8118214ec])),
        c1: Fq(FqRepr([0x6e849f1ea0aa4757, 0xaa1c7b6d89f89141, 0xb6e713cdfae0ca3a, 0x26694fbb4e82ebc3])),
    },
    // Fq2(u + 9)**(((q^2) - 1) / 3)
    // this one and other below are recalculated manually
    Fq2 {
        c0: Fq(FqRepr([0x3350c88e13e80b9c, 0x7dce557cdb5e56b9, 0x6001b4b8b615564a, 0x2682e617020217e0])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0])),
    },
    // Fq2(u + 9)**(((q^3) - 1) / 3)
    Fq2 {
        c0: Fq(FqRepr([0xc9af22f716ad6bad, 0xb311782a4aa662b2, 0x19eeaf64e248c7f4, 0x20273e77e3439f82])),
        c1: Fq(FqRepr([0xacc02860f7ce93ac, 0x3933d5817ba76b4c, 0x69e6188b446c8467, 0x0a46036d4417cc55])),
    },
    // Fq2(u + 9)**(((q^4) - 1) / 3)
    Fq2 {
        c0: Fq(FqRepr([0x71930c11d782e155, 0xa6bb947cffbe3323, 0xaa303344d4741444, 0x2c3b3f0d26594943])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0])),
    },
    // Fq2(u + 9)**(((q^5) - 1) / 3)
    Fq2 {
        c0: Fq(FqRepr([0xf91aba2654e8e3b1, 0x4771cb2fdc92ce12, 0xdcb16ae0fc8bdf35, 0x274aa195cd9d8be4])),
        c1: Fq(FqRepr([0x5cfc50ae18811f8b, 0x4bb28433cb43988c, 0x4fd35f13c3b56219, 0x301949bd2fc8883a])),
    },
];

pub const FROBENIUS_COEFF_FQ6_C2: [Fq2; 6] = [
    // Fq2(u + 1)**(((2q^0) - 2) / 3)
    Fq2 {
        c0: Fq(FqRepr([0xd35d438dc58f0d9d, 0x0a78eb28f5c70b3d, 0x666ea36f7879462c, 0x0e0a77c19a07df2f])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0])),
    },
    // Fq2(u + 1)**(((2q^1) - 2) / 3)
    Fq2 {
        c0: Fq(FqRepr([0x7361d77f843abe92, 0xa5bb2bd3273411fb, 0x9c941f314b3e2399, 0x15df9cddbb9fd3ec])),
        c1: Fq(FqRepr([0x5dddfd154bd8c949, 0x62cb29a5a4445b60, 0x37bc870a0c7dd2b9, 0x24830a9d3171f0fd])),
    },
    // Fq2(u + 1)**(((2q^2) - 2) / 3)
    Fq2 {
        c0: Fq(FqRepr([0x71930c11d782e155, 0xa6bb947cffbe3323, 0xaa303344d4741444, 0x2c3b3f0d26594943])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0])),
    },
    // Fq2(u + 1)**(((2q^3) - 2) / 3)
    Fq2 {
        c0: Fq(FqRepr([0x448a93a57b6762df, 0xbfd62df528fdeadf, 0xd858f5d00e9bd47a, 0x06b03d4d3476ec58])),
        c1: Fq(FqRepr([0x2b19daf4bcc936d1, 0xa1a54e7a56f4299f, 0xb533eee05adeaef1, 0x170c812b84dda0b2])),
    },
    // Fq2(u + 1)**(((2q^4) - 2) / 3)
    Fq2 {
        c0: Fq(FqRepr([0x3350c88e13e80b9c, 0x7dce557cdb5e56b9, 0x6001b4b8b615564a, 0x2682e617020217e0])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0])),
    },
    // Fq2(u + 1)**(((2q^5) - 2) / 3)
    Fq2 {
        c0: Fq(FqRepr([0x843420f1d8dadbd6, 0x31f010c9183fcdb2, 0x436330b527a76049, 0x13d47447f11adfe4])),
        c1: Fq(FqRepr([0xef494023a857fa74, 0x2a925d02d5ab101a, 0x83b015829ba62f10, 0x2539111d0c13aea3])),
    },
];

// non_residue^((modulus^i-1)/6) for i=0,...,11
pub const FROBENIUS_COEFF_FQ12_C1: [Fq2; 12] = [
    // Fq2(u + 1)**(((q^0) - 1) / 6)
    Fq2 {
        c0: Fq(FqRepr([0xd35d438dc58f0d9d, 0x0a78eb28f5c70b3d, 0x666ea36f7879462c, 0x0e0a77c19a07df2f])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0])),
    },
    // Fq2(u + 1)**(((q^1) - 1) / 6)
    Fq2 {
        c0: Fq(FqRepr([0xaf9ba69633144907, 0xca6b1d7387afb78a, 0x11bded5ef08a2087, 0x02f34d751a1f3a7c])),
        c1: Fq(FqRepr([0xa222ae234c492d72, 0xd00f02a4565de15b, 0xdc2ff3a253dfc926, 0x10a75716b3899551])),
    },
    // Fq2(u + 1)**(((q^2) - 1) / 6)
    Fq2 {
        c0: Fq(FqRepr([0xca8d800500fa1bf2, 0xf0c5d61468b39769, 0x0e201271ad0d4418, 0x04290f65bad856e6])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0])),
    },
    // Fq2(u + 1)**(((q^3) - 1) / 6)
    Fq2 {
        c0: Fq(FqRepr([0x365316184e46d97d, 0x0af7129ed4c96d9f, 0x659da72fca1009b5, 0x08116d8983a20d23])),
        c1: Fq(FqRepr([0xb1df4af7c39c1939, 0x3d9f02878a73bf7f, 0x9b2220928caf0ae0, 0x26684515eff054a6])),
    },
    // Fq2(u + 1)**(((q^4) - 1) / 6)
    Fq2 {
        c0: Fq(FqRepr([0x3350c88e13e80b9c, 0x7dce557cdb5e56b9, 0x6001b4b8b615564a, 0x2682e617020217e0])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0])),
    },
    // Fq2(u + 1)**(((q^5) - 1) / 6)
    Fq2 {
        c0: Fq(FqRepr([0x86b76f821b329076, 0x408bf52b4d19b614, 0x53dfb9d0d985e92d, 0x051e20146982d2a7])),
        c1: Fq(FqRepr([0x0fbc9cd47752ebc7, 0x6d8fffe33415de24, 0xbef22cf038cf41b9, 0x15c0edff3c66bf54])),
    },
    // Fq2(u + 1)**(((q^6) - 1) / 6)
    Fq2 {
        c0: Fq(FqRepr([0x68c3488912edefaa, 0x8d087f6872aabf4f, 0x51e1a24709081231, 0x2259d6b14729c0fa])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0])),
    },
    // Fq2(u + 1)**(((q^7) - 1) / 6)
    Fq2 {
        c0: Fq(FqRepr([0x8c84e580a568b440, 0xcd164d1de0c21302, 0xa692585790f737d5, 0x2d7100fdc71265ad])),
        c1: Fq(FqRepr([0x99fdddf38c33cfd5, 0xc77267ed1213e931, 0xdc2052142da18f36, 0x1fbcf75c2da80ad7])),
    },
    // Fq2(u + 1)**(((q^8) - 1) / 6)
    Fq2 {
        c0: Fq(FqRepr([0x71930c11d782e155, 0xa6bb947cffbe3323, 0xaa303344d4741444, 0x2c3b3f0d26594943])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0])),
    },
    // Fq2(u + 1)**(((q^9) - 1) / 6)
    Fq2 {
        c0: Fq(FqRepr([0x05cd75fe8a3623ca, 0x8c8a57f293a85cee, 0x52b29e86b7714ea8, 0x2852e0e95d8f9306])),
        c1: Fq(FqRepr([0x8a41411f14e0e40e, 0x59e26809ddfe0b0d, 0x1d2e2523f4d24d7d, 0x09fc095cf1414b83])),
    },
    // Fq2(u + 1)**(((q^10) - 1) / 6)
    Fq2 {
        c0: Fq(FqRepr([0x08cfc388c494f1ab, 0x19b315148d1373d4, 0x584e90fdcb6c0213, 0x09e1685bdf2f8849])),
        c1: Fq(FqRepr([0x0, 0x0, 0x0, 0x0])),
    },
    // Fq2(u + 1)**(((q^11) - 1) / 6)
    Fq2 {
        c0: Fq(FqRepr([0xb5691c94bd4a6cd1, 0x56f575661b581478, 0x64708be5a7fb6f30, 0x2b462e5e77aecd82])),
        c1: Fq(FqRepr([0x2c63ef42612a1180, 0x29f16aae345bec69, 0xf95e18c648b216a4, 0x1aa36073a4cae0d4])),
    },
];

// -((2**256) mod q) mod q
pub const NEGATIVE_ONE: Fq = Fq(FqRepr([0x974bc177a0000006, 0xf13771b2da58a367, 0x51e1a2470908122e, 0x2259d6b14729c0fa]));


// used in pairing ceritificate
pub const FQ12_CUBIC_NON_RESIDUE : Fq12 = Fq12 { 
    c0: Fq6 { 
        c0: Fq2 { 
            c0: Fq(FqRepr([10142183612076115840, 10336180170933656228, 470517678848978077, 1225265646762263217])), 
            c1: Fq(FqRepr([11348019713242347020, 2754915599916440970, 15254070011878377052, 582093106373754232])), 
        }, 
        c1: Fq2 { 
            c0: Fq(FqRepr([397496596320767129, 15618912924840201202, 4512184083626728156, 117261539653304508])),
            c1: Fq(FqRepr([10874872625388672849, 8697447527460198735, 17201932474568901787, 3239181218336645708])),  
        }, 
        c2: Fq2 { 
            c0: Fq(FqRepr([3319746353403080977, 12107531098513439570, 7977768322294476906, 2777360315899794408])),
            c1: Fq(FqRepr([4282545994011025486, 13446103940465242757, 621706771151745465, 121235365032730123])),
        }
    }, 
    c1: Fq6 { 
        c0: Fq2 { 
            c0: Fq(FqRepr([10835993777695610176, 7471612555186838325, 5669943730362164093, 1589478061029939273])),
            c1: Fq(FqRepr([13798745335319118500, 9933747591465430734, 17986311991554267884, 2884267725610847413])), 
        }, 
        c1: Fq2 { 
            c0: Fq(FqRepr([9010191163654582374, 16916065212229622062, 13444508417525226022, 3305365841183382359])),
            c1: Fq(FqRepr([16399630071410982746, 11060057513443281354, 11904104725095752194, 2793879429500272640])),
        }, 
        c2: Fq2 { 
            c0: Fq(FqRepr([12230496698620887775, 2879834267624199743, 1836535646417438336, 2949238563047019140])),
            c1: Fq(FqRepr([2887317464157129146, 9549941403471574669, 13305129467053585697, 2048982498456539113])),
        } 
    } 
};

pub const ROOT_27_OF_UNITY : Fq6 = Fq6 { 
    c0: Fq2 { 
        c0: Fq(FqRepr([0, 0, 0, 0])), 
        c1: Fq(FqRepr([0, 0, 0, 0])) 
    }, 
    c1: Fq2 { 
        c0: Fq(FqRepr([0, 0, 0, 0])), 
        c1: Fq(FqRepr([0, 0, 0, 0])) 
    }, 
    c2: Fq2 { 
        c0: Fq(FqRepr([194717170545481332, 10713050357451161581, 34962824991068036, 114905412056219140])),
        c1: Fq(FqRepr([4228108808443533428, 10892079493576126185, 17015400802854381179, 1820057530234360532])),
    } 
};


#[cfg(test)]
use rand::{Rand, SeedableRng, XorShiftRng};

#[test]
fn test_fq_repr_from() {
    assert_eq!(FqRepr::from(100), FqRepr([100, 0, 0, 0]));
    assert_eq!(FqRepr::from(3), FqRepr([3, 0, 0, 0]));
}

#[test]
fn test_fq_repr_is_odd() {
    assert!(!FqRepr::from(0).is_odd());
    assert!(FqRepr::from(0).is_even());
    assert!(FqRepr::from(1).is_odd());
    assert!(!FqRepr::from(1).is_even());
    assert!(!FqRepr::from(324834872).is_odd());
    assert!(FqRepr::from(324834872).is_even());
    assert!(FqRepr::from(324834873).is_odd());
    assert!(!FqRepr::from(324834873).is_even());
}

#[test]
fn test_fq_repr_num_bits() {
    let mut a = FqRepr::from(0);
    assert_eq!(0, a.num_bits());
    a = FqRepr::from(1);
    for i in 1..257 {
        assert_eq!(i, a.num_bits());
        a.mul2();
    }
    assert_eq!(0, a.num_bits());
}

#[test]
fn test_fq_is_valid() {
    print!("modulus = {}\n", MODULUS);
    print!("R = {}\n", R);
    let mut a = Fq(MODULUS);
    assert!(!a.is_valid());
    a.0.sub_noborrow(&FqRepr::from(1));
    assert!(a.is_valid());
    assert!(Fq(FqRepr::from(0)).is_valid());
    assert!(Fq(FqRepr([0xdf4671abd14dab3e, 0xe2dc0c9f534fbd33, 0x31ca6c880cc444a6, 0x257a67e70ef33359])).is_valid());
    assert!(!Fq(FqRepr([0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff,])).is_valid());

    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    for _ in 0..1000 {
        let a = Fq::rand(&mut rng);
        assert!(a.is_valid());
    }
}

#[test]
fn test_fq_repr_display() {
    assert_eq!(
        format!("{}", Fq::into_repr(&Fq::one())),
        "0x0000000000000000000000000000000000000000000000000000000000000001".to_string()
    );
    assert_eq!(format!("{}", FqRepr([0, 0, 0, 0])), "0x0000000000000000000000000000000000000000000000000000000000000000".to_string());
}

#[test]
fn test_fq_num_bits() {
    assert_eq!(Fq::NUM_BITS, 254);
    assert_eq!(Fq::CAPACITY, 253);
}

#[test]
fn test_fq_sqrt() {
    use ff::SqrtField;

    let mut rng = XorShiftRng::from_seed([0x5dbe6259, 0x8d313d76, 0x3237db17, 0xe5bc0654]);

    assert_eq!(Fq::zero().sqrt().unwrap(), Fq::zero());

    for _ in 0..1000 {
        // Ensure sqrt(a^2) = a or -a
        let a = Fq::rand(&mut rng);
        let mut nega = a;
        nega.negate();
        let mut b = a;
        b.square();

        let b = b.sqrt().unwrap();

        assert!(a == b || nega == b);
    }

    for _ in 0..1000 {
        // Ensure sqrt(a)^2 = a for random a
        let a = Fq::rand(&mut rng);

        if let Some(mut tmp) = a.sqrt() {
            tmp.square();

            assert_eq!(a, tmp);
        }
    }
}

#[test]
fn test_fq_sqrt_2() {
    use ff::SqrtField;

    let x = Fq::from_str("4").unwrap();
    print!("x = {}\n", x);
    if let Some(y) = x.sqrt() {
        print!("y = {}\n", y);
        let mut y_other = y;
        y_other.negate();
        print!("y' = {}\n", y_other);
    }
}

#[test]
fn fq_field_tests() {
    crate::tests::field::random_field_tests::<Fq>();
    crate::tests::field::random_sqrt_tests::<Fq>();
    crate::tests::field::random_frobenius_tests::<Fq, _>(Fq::char(), 13);
    crate::tests::field::from_str_tests::<Fq>();
}
