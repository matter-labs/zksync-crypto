use crate::traits::HashParams;
use franklin_crypto::bellman::{Engine, Field, PrimeField};
use crate::common::domain_strategy::DomainStrategy;
use super::params::Poseidon2Params;
use crate::traits::Sbox;

/// Receives inputs whose length `known` prior(fixed-length).
/// Also uses custom domain strategy which basically sets value of capacity element to
/// length of input and applies a padding rule which makes input size equals to multiple of
/// rate parameter.
pub fn poseidon2_hash<
    E: Engine,
    const L: usize
>(input: &[E::Fr; L]) -> [E::Fr; 2] {
    const WIDTH: usize = 3;
    const RATE: usize = 2;

    let params = Poseidon2Params::<E, RATE, WIDTH>::default();
    crate::generic_hash(&params, input, None)
}

pub(crate) fn poseidon2_round_function<
    E: Engine,
    const RATE: usize,
    const WIDTH: usize,
>(
    state: &mut [E::Fr; WIDTH],
    params: &Poseidon2Params<E, RATE, WIDTH>,
) {
    debug_assert!(params.full_rounds & 1 == 0);
    let half_of_full_rounds = params.number_of_full_rounds() / 2;

    // Linear layer at beginning
    poseidon2_matmul_external::<E, WIDTH>(state);

    for r in 0..half_of_full_rounds {
        add_rc::<E, WIDTH>(state, &params.round_constants[r]);
        apply_sbox::<E>(state, &params.alpha);
        poseidon2_matmul_external::<E, WIDTH>(state);
    }

    for r in half_of_full_rounds..(half_of_full_rounds + params.partial_rounds) {
        state[0].add_assign(&params.round_constants[r][0]);
        apply_sbox::<E>(&mut state[..1], &params.alpha);
        poseidon2_matmul_internal::<E, WIDTH>(state, &params.diag_internal_matrix);
    }
    
    for r in (half_of_full_rounds + params.partial_rounds)..(2*half_of_full_rounds + params.partial_rounds) {
        add_rc::<E, WIDTH>(state, &params.round_constants[r]);
        apply_sbox::<E>(state, &params.alpha);
        poseidon2_matmul_external::<E, WIDTH>(state);
    }
}

pub(crate) fn poseidon2_matmul_external<
    E: Engine,
    const WIDTH: usize,
>(
    state: &mut [E::Fr; WIDTH]
) {
    match WIDTH {
        2 => {
            // Matrix circ(2, 1)
            let mut sum = state[0];
            sum.add_assign(&state[1]);
            state[0].add_assign(&sum);
            state[1].add_assign(&sum);
        }
        3 => {
            // Matrix circ(2, 1, 1)
            let mut sum = state[0];
            sum.add_assign(&state[1]);
            sum.add_assign(&state[2]);
            state[0].add_assign(&sum);
            state[1].add_assign(&sum);
            state[2].add_assign(&sum);
        }
        4 => {
            // Applying cheap 4x4 MDS matrix to each 4-element part of the state
            matmul_m4::<E, WIDTH>(state);
        }
        8 | 12 | 16 | 20 | 24 => {
            // Applying cheap 4x4 MDS matrix to each 4-element part of the state
            matmul_m4::<E, WIDTH>(state);

            // Applying second cheap matrix for t > 4
            let t4 = WIDTH / 4;
            let mut stored = [E::Fr::zero(); 4];
            for l in 0..4 {
                stored[l] = state[l];
                for j in 1..t4 {
                    stored[l].add_assign(&state[4 * j + l]);
                }
            }
            for i in 0..WIDTH {
                state[i].add_assign(&stored[i % 4]);
            }
        }
        _ => {
            panic!()
        }
    }
}

fn matmul_m4<
    E: Engine,
    const WIDTH: usize,
>(
    state: &mut [E::Fr; WIDTH]
) {
    // Mul each 4-element chunk by
    // [5, 7, 1, 3]
    // [4, 6, 1, 1]
    // [1, 3, 5, 7]
    // [1, 1, 4, 6]

    let t4 = WIDTH / 4;
    for i in 0..t4 {
        let start_index = i * 4;
        let mut t_0 = state[start_index];
        t_0.add_assign(&state[start_index + 1]);
        let mut t_1 = state[start_index + 2];
        t_1.add_assign(&state[start_index + 3]);
        let mut t_2 = state[start_index + 1];
        t_2.double();
        t_2.add_assign(&t_1);
        let mut t_3 = state[start_index + 3];
        t_3.double();
        t_3.add_assign(&t_0);
        let mut t_4 = t_1;
        t_4.double();
        t_4.double();
        t_4.add_assign(&t_3);
        let mut t_5 = t_0;
        t_5.double();
        t_5.double();
        t_5.add_assign(&t_2);
        let mut t_6 = t_3;
        t_6.add_assign(&t_5);
        let mut t_7 = t_2;
        t_7.add_assign(&t_4);
        state[start_index] = t_6;
        state[start_index + 1] = t_5;
        state[start_index + 2] = t_7;
        state[start_index + 3] = t_4;
    }
}


pub(crate) fn poseidon2_matmul_internal<
    E: Engine,
    const WIDTH: usize,
>(
    state: &mut [E::Fr; WIDTH],
    diag_internal_matrix: &[E::Fr; WIDTH]
) {
    match WIDTH {
        2 => {
            // [2, 1]
            // [1, 3]
            debug_assert_eq!(diag_internal_matrix[0], E::Fr::from_str("2").unwrap());
            debug_assert_eq!(diag_internal_matrix[1], E::Fr::from_str("3").unwrap());

            let mut sum = state[0];
            sum.add_assign(&state[1]);
            state[0].add_assign(&sum);
            state[1].double();
            state[1].add_assign(&sum);
        }
        3 => {
            // [2, 1, 1]
            // [1, 2, 1]
            // [1, 1, 3]
            debug_assert_eq!(diag_internal_matrix[0], E::Fr::from_str("2").unwrap());
            debug_assert_eq!(diag_internal_matrix[1], E::Fr::from_str("2").unwrap());
            debug_assert_eq!(diag_internal_matrix[2], E::Fr::from_str("3").unwrap());

            let mut sum = state[0];
            sum.add_assign(&state[1]);
            sum.add_assign(&state[2]);
            state[0].add_assign(&sum);
            state[1].add_assign(&sum);
            state[2].double();
            state[2].add_assign(&sum);
        }
        4 | 8 | 12 | 16 | 20 | 24 => {
            // Compute state sum
            let mut sum = state[0];
            state
                .iter()
                .skip(1)
                .take(WIDTH-1)
                .for_each(|el| sum.add_assign(el));
            // Add sum + (diag entry - 1) * element to each element
            for i in 0..WIDTH {
                let mut coeff = diag_internal_matrix[i];
                coeff.sub_assign(&E::Fr::one());
                state[i].mul_assign(&coeff);
                state[i].add_assign(&sum);
            }
        }
        _ => {
            panic!()
        }
    }
}

pub(crate) fn apply_sbox<
    E: Engine,
>(
    elements: &mut [E::Fr],
    sbox: &Sbox
) {
    debug_assert!(sbox == &Sbox::Alpha(5));

    for element in elements.iter_mut() {
        let mut res = *element;
        res.square();
        res.square();
        res.mul_assign(&element);
        *element = res;
    }
}

pub(crate) fn add_rc<
    E: Engine,
    const WIDTH: usize,
>(
    elements: &mut [E::Fr; WIDTH],
    constants: &[E::Fr; WIDTH]
) {
    for (element, constant) in elements.iter_mut().zip(constants.iter()) {
        element.add_assign(constant);
    }
}
