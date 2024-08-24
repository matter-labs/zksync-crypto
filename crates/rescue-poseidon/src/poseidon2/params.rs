use franklin_crypto::bellman::{Engine, Field};

use crate::common::params::InnerHashParameters;
use crate::traits::{CustomGate, HashFamily, HashParams, Sbox};
use franklin_crypto::bellman::PrimeField;

#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Poseidon2Params<E: Engine, const RATE: usize, const WIDTH: usize> {
    #[serde(serialize_with = "crate::serialize_array_of_arrays")]
    #[serde(deserialize_with = "crate::deserialize_array_of_arrays")]
    pub(crate) mds_external_matrix: [[E::Fr; WIDTH]; WIDTH],
    #[serde(with = "crate::BigArraySerde")]
    pub(crate) diag_internal_matrix: [E::Fr; WIDTH],
    #[serde(serialize_with = "crate::serialize_vec_of_arrays")]
    #[serde(deserialize_with = "crate::deserialize_vec_of_arrays")]
    pub(crate) round_constants: Vec<[E::Fr; WIDTH]>,
    pub(crate) alpha: Sbox,
    pub(crate) full_rounds: usize,
    pub(crate) partial_rounds: usize,
    pub(crate) custom_gate: CustomGate,
}

impl<E: Engine, const RATE: usize, const WIDTH: usize> PartialEq for Poseidon2Params<E, RATE, WIDTH> {
    fn eq(&self, other: &Self) -> bool {
        self.hash_family() == other.hash_family()
    }
}

impl<E: Engine, const RATE: usize, const WIDTH: usize> Default for Poseidon2Params<E, RATE, WIDTH> {
    fn default() -> Self {
        let security_level = 80; // TODO: check, but we actually don't use it anywhere

        // Number of rounds from the original Poseidon2 implementation
        // https://github.com/HorizenLabs/poseidon2
        let full_rounds = 8;
        let partial_rounds = 56;

        let mut params = InnerHashParameters::<E, RATE, WIDTH>::new(security_level, full_rounds, partial_rounds);

        // Same constants as in the Poseidon
        let number_of_rounds = full_rounds + partial_rounds;
        let rounds_tag = b"Rescue_f";
        params.compute_round_constants(number_of_rounds, rounds_tag);

        let mds_external_matrix = poseidon2_external_matrix::<E, WIDTH>();
        let diag_internal_matrix = poseidon2_internal_matrix::<E, WIDTH>();

        let mut round_constants = params.round_constants().to_owned();
        for i in (full_rounds / 2)..(full_rounds / 2) + partial_rounds {
            for j in 1..WIDTH {
                round_constants[i][j] = E::Fr::zero();
            }
        }

        let alpha = 5u64;

        Self {
            alpha: Sbox::Alpha(alpha),
            full_rounds: params.full_rounds,
            partial_rounds: params.partial_rounds,
            custom_gate: CustomGate::QuinticWidth4,

            mds_external_matrix,
            diag_internal_matrix,
            round_constants,
        }
    }
}

impl<E: Engine, const RATE: usize, const WIDTH: usize> HashParams<E, RATE, WIDTH> for Poseidon2Params<E, RATE, WIDTH> {
    fn hash_family(&self) -> HashFamily {
        HashFamily::Poseidon2
    }

    fn constants_of_round(&self, round: usize) -> &[E::Fr; WIDTH] {
        &self.round_constants[round]
    }

    fn mds_matrix(&self) -> &[[E::Fr; WIDTH]; WIDTH] {
        panic!("Poseidon2 uses 2 types of matrixes")
    }

    fn number_of_full_rounds(&self) -> usize {
        self.full_rounds
    }

    fn number_of_partial_rounds(&self) -> usize {
        self.partial_rounds
    }

    fn alpha(&self) -> &Sbox {
        &self.alpha
    }

    fn alpha_inv(&self) -> &Sbox {
        unimplemented!("Poseidon doesn't have inverse direction")
    }

    fn optimized_round_constants(&self) -> &[[E::Fr; WIDTH]] {
        unimplemented!("Poseidon doesn't use optimized constants")
    }

    fn optimized_mds_matrixes(&self) -> (&[[E::Fr; WIDTH]; WIDTH], &[[[E::Fr; WIDTH]; WIDTH]]) {
        unimplemented!("Poseidon doesn't use optimized matrixes")
    }

    fn custom_gate(&self) -> CustomGate {
        self.custom_gate
    }

    fn use_custom_gate(&mut self, custom_gate: CustomGate) {
        self.custom_gate = custom_gate;
    }

    fn try_to_poseidon2_params(&self) -> Option<&crate::poseidon2::Poseidon2Params<E, RATE, WIDTH>> {
        Some(self)
    }
}

fn poseidon2_external_matrix<E: Engine, const WIDTH: usize>() -> [[E::Fr; WIDTH]; WIDTH] {
    let one = E::Fr::one();
    let two = E::Fr::from_str("2").unwrap();

    let mut result = [[E::Fr::zero(); WIDTH]; WIDTH];
    match WIDTH {
        2 => {
            // circ(2, 1)
            result[0][0] = two;
            result[0][1] = one;
            result[1][0] = one;
            result[1][1] = two;
        }
        3 => {
            // circ(2, 1, 1)
            result[0][0] = two;
            result[0][1] = one;
            result[0][2] = one;
            result[1][0] = one;
            result[1][1] = two;
            result[1][2] = one;
            result[2][0] = one;
            result[2][1] = one;
            result[2][2] = two;
        }
        _ => {
            assert!(WIDTH > 0 && WIDTH % 4 == 0);

            // [5, 7, 1, 3]
            // [4, 6, 1, 1]
            // [1, 3, 5, 7]
            // [1, 1, 4, 6]
            let three = E::Fr::from_str("3").unwrap();
            let four = E::Fr::from_str("4").unwrap();
            let five = E::Fr::from_str("5").unwrap();
            let six = E::Fr::from_str("6").unwrap();
            let seven = E::Fr::from_str("7").unwrap();

            let m_4_mds_matrix = [[five, seven, one, three], [four, six, one, one], [one, three, five, seven], [one, one, four, six]];

            // circ(2*M4, M4, ..., M4)
            for i in 0..WIDTH {
                for j in 0..WIDTH {
                    result[i][j] = m_4_mds_matrix[i % 4][j % 4];
                    if i / 4 == j / 4 {
                        result[i][j].mul_assign(&two);
                    }
                }
            }
        }
    };

    result
}

fn poseidon2_internal_matrix<E: Engine, const WIDTH: usize>() -> [E::Fr; WIDTH] {
    let two = E::Fr::from_str("2").unwrap();
    let three = E::Fr::from_str("3").unwrap();

    let mut result = [E::Fr::zero(); WIDTH];
    match WIDTH {
        3 => {
            result[0] = two;
            result[1] = two;
            result[2] = three;
        }
        _ => todo!("poseidon_2_internal_matrix for WIDTH == {}", WIDTH),
    };

    result
}
