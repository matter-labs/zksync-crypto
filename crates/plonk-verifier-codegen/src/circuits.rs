use rescue_poseidon::franklin_crypto::{
    bellman::{
        bn256::{Bn256, Fr},
        kate_commitment::{Crs, CrsForMonomialForm},
        plonk::{
            better_better_cs::{
                cs::{
                    ArithmeticTerm, Circuit, ConstraintSystem, Gate, GateInternal, LookupTableApplication, MainGate, MainGateTerm, PlonkConstraintSystemParams,
                    PlonkCsWidth4WithNextStepAndCustomGatesParams, PlonkCsWidth4WithNextStepParams, PolyIdentifier, TrivialAssembly, VerificationKey, Width4MainGateWithDNext,
                },
                gates::selector_optimized_with_d_next::SelectorOptimizedWidth4MainGateWithDNext,
                proof::Proof,
            },
            commitments::transcript::{keccak_transcript::RollingKeccakTranscript, Transcript},
        },
        worker::Worker,
        Engine, Field, PrimeField, SynthesisError,
    },
    plonk::circuit::{
        allocated_num::{AllocatedNum, Num},
        boolean::Boolean,
        custom_rescue_gate::Rescue5CustomGate,
        linear_combination::LinearCombination,
        Width4WithCustomGates,
    },
};

use rescue_poseidon::{circuit_generic_hash, CustomGate, HashParams, RescueParams};

use paste::paste;

macro_rules! circuit_inner {
    ($id:ident, $main_gate:ty, $declare_rescue:expr, $( $synth:tt ),*) => {
        #[derive(Clone, Debug)]
        pub struct $id;

        impl<E: Engine> Circuit<E> for $id {
            type MainGate = $main_gate;

            fn synthesize<CS: ConstraintSystem<E>>(
                &self,
                cs: &mut CS,
            ) -> Result<(), SynthesisError> {
                inner_circuit_main_gate_part(cs)?;
                $(
                    $synth(cs)?;
                )*

                Ok(())
            }

            fn declare_used_gates() -> Result<Vec<Box<dyn GateInternal<E>>>, SynthesisError> {
                let has_rescue = $declare_rescue;
                if has_rescue{
                    Ok(vec![
                        Self::MainGate::default().into_internal(),
                        Rescue5CustomGate::default().into_internal(),
                    ])
                }else{
                    Ok(vec![
                        Self::MainGate::default().into_internal(),
                    ])
                }
            }
        }
    };
}

#[macro_export]
macro_rules! circuit {
    ($id:ident, $main_gate:ty) => {
        circuit_inner!($id, $main_gate, false, inner_circuit_main_gate_part);
        paste! {
            circuit_inner!([<$id WithLookup>], $main_gate, false,  inner_circuit_lookup_part);
        }
        paste! {
            circuit_inner!([<$id WithRescue>], $main_gate, true, inner_circuit_rescue_part);
        }
        paste! {
            circuit_inner!([<$id WithLookupAndRescue>], $main_gate, true, inner_circuit_lookup_part, inner_circuit_rescue_part);
        }
    };
}

circuit!(MockCircuit, Width4MainGateWithDNext);
circuit!(MockCircuitSelectorOptimized, SelectorOptimizedWidth4MainGateWithDNext);

fn inner_circuit_main_gate_part<E: Engine, CS: ConstraintSystem<E>>(cs: &mut CS) -> Result<(), SynthesisError> {
    for _ in 0..32 {
        let a = Num::alloc(cs, Some(E::Fr::one()))?;
        let b = Num::alloc(cs, Some(E::Fr::zero()))?;
        let flag = Boolean::alloc(cs, Some(true))?;
        let c = Num::conditionally_select(cs, &flag, &a, &b)?;
        let is_equal = Num::equals(cs, &a, &c)?;
        Boolean::enforce_equal(cs, &is_equal, &Boolean::Constant(true))?;

        let mut lc = LinearCombination::zero();
        for idx in 0..6 {
            let el = E::Fr::from_str(&idx.to_string()).unwrap();
            let allocated = Num::alloc(cs, Some(el))?;
            lc.add_assign_number_with_coeff(&allocated, E::Fr::one());
        }
        let sum = Num::alloc(cs, Some(E::Fr::from_str(&20.to_string()).unwrap()))?;
        let mut minus_one = E::Fr::one();
        minus_one.negate();
        lc.add_assign_number_with_coeff(&sum, minus_one);
    }
    let a = Num::alloc(cs, Some(E::Fr::one()))?;
    let b = Num::alloc(cs, Some(E::Fr::one()))?;
    let expected = cs.alloc_input(|| Ok(E::Fr::from_str(&2.to_string()).unwrap()))?;
    let actual = a.add(cs, &b)?;

    let mut term = MainGateTerm::<E>::new();
    term.add_assign(ArithmeticTerm::from_variable(expected));
    term.sub_assign(ArithmeticTerm::from_variable(actual.get_variable().get_variable()));
    cs.allocate_main_gate(term)?;

    Ok(())
}

fn inner_circuit_lookup_part<E: Engine, CS: ConstraintSystem<E>>(cs: &mut CS) -> Result<(), SynthesisError> {
    // add dummy lookup table queries
    let dummy = CS::get_dummy_variable();
    // need to create a table (any)
    let columns = vec![PolyIdentifier::VariablesPolynomial(0), PolyIdentifier::VariablesPolynomial(1), PolyIdentifier::VariablesPolynomial(2)];
    let range_table = LookupTableApplication::new_range_table_of_width_3(2, columns.clone())?;
    let _range_table_name = range_table.functional_name();

    let xor_table = LookupTableApplication::new_xor_table(2, columns.clone())?;
    let _xor_table_name = xor_table.functional_name();

    let and_table = LookupTableApplication::new_and_table(2, columns)?;
    let and_table_name = and_table.functional_name();

    cs.add_table(range_table)?;
    cs.add_table(xor_table)?;
    cs.add_table(and_table)?;

    let binary_x_value = E::Fr::from_str("3").unwrap();
    let binary_y_value = E::Fr::from_str("1").unwrap();

    let t = AllocatedNum::zero(cs);
    let tt = AllocatedNum::one(cs);
    let ttt = t.mul(cs, &tt)?;
    ttt.inputize(cs)?;

    let binary_x = cs.alloc(|| Ok(binary_x_value))?;

    let binary_y = cs.alloc(|| Ok(binary_y_value))?;

    let table = cs.get_table(&and_table_name)?;
    let num_keys_and_values = table.width();

    let and_result_value = table.query(&[binary_x_value, binary_y_value])?[0];

    let binary_z = cs.alloc(|| Ok(and_result_value))?;

    cs.begin_gates_batch_for_step()?;

    let vars = [binary_x, binary_y, binary_z, dummy];
    cs.allocate_variables_without_gate(&vars, &[])?;

    cs.apply_single_lookup_gate(&vars[..num_keys_and_values], table)?;

    cs.end_gates_batch_for_step()?;

    Ok(())
}

fn inner_circuit_rescue_part<E: Engine, CS: ConstraintSystem<E>>(cs: &mut CS) -> Result<(), SynthesisError> {
    // make single rescue hash to satisfy gate requirements of declaration
    let mut params = RescueParams::default();
    params.use_custom_gate(CustomGate::QuinticWidth4);

    let elem = Num::alloc(cs, Some(E::Fr::from_str("42").unwrap()))?;
    let _ = circuit_generic_hash::<_, _, _, 2, 3, 2>(cs, &[elem, elem], &params, None)?;

    Ok(())
}

#[test]
fn test_create_proof_for_all_circuits() {
    type T = RollingKeccakTranscript<Fr>;

    let base_dir = std::env::var("PLONK_VERIFIER_DATA_DIR").expect("output dir for output files");
    let out_dir = format!("{}/std", base_dir,);
    println!("out dir {}", out_dir);
    println!("Std main gate");
    create_proof_for_circuit::<_, PlonkCsWidth4WithNextStepParams, Width4MainGateWithDNext, T>(MockCircuit {}, &format!("{}/std", &out_dir));
    println!("Std main gate with lookup");
    create_proof_for_circuit::<_, PlonkCsWidth4WithNextStepParams, Width4MainGateWithDNext, T>(MockCircuitWithLookup {}, &format!("{}/std_with_lookup", &out_dir));
    println!("Std main gate with sbox custom gate");
    create_proof_for_circuit::<_, PlonkCsWidth4WithNextStepAndCustomGatesParams, Width4MainGateWithDNext, T>(MockCircuitWithRescue {}, &format!("{}/std_with_rescue", &out_dir));
    println!("Std main gate with lookup and sbox custom gate");
    create_proof_for_circuit::<_, PlonkCsWidth4WithNextStepAndCustomGatesParams, Width4MainGateWithDNext, T>(MockCircuitWithLookupAndRescue {}, &format!("{}/std_with_lookup_and_rescue", &out_dir));

    let out_dir = format!("{}/optimized", base_dir,);
    println!("SelectorOptimized main gate");
    create_proof_for_circuit::<_, PlonkCsWidth4WithNextStepAndCustomGatesParams, SelectorOptimizedWidth4MainGateWithDNext, T>(MockCircuitSelectorOptimized {}, &format!("{}/optimized", &out_dir));
    println!("SelectorOptimized main gate with lookup");
    create_proof_for_circuit::<_, PlonkCsWidth4WithNextStepAndCustomGatesParams, SelectorOptimizedWidth4MainGateWithDNext, T>(
        MockCircuitSelectorOptimizedWithLookup {},
        &format!("{}/optimized_with_lookup", &out_dir),
    );
    println!("SelectorOptimized main gate with sbox custom gate");
    create_proof_for_circuit::<_, PlonkCsWidth4WithNextStepAndCustomGatesParams, SelectorOptimizedWidth4MainGateWithDNext, T>(
        MockCircuitSelectorOptimizedWithRescue {},
        &format!("{}/optimized_with_rescue", &out_dir),
    );
    println!("SelectorOptimized main gate with lookup and sbox custom gate");
    create_proof_for_circuit::<_, PlonkCsWidth4WithNextStepAndCustomGatesParams, SelectorOptimizedWidth4MainGateWithDNext, T>(
        MockCircuitSelectorOptimizedWithLookupAndRescue {},
        &format!("{}/optimized_with_lookup_and_rescue", &out_dir),
    );
}

fn init_crs(worker: &Worker, domain_size: usize) -> Crs<Bn256, CrsForMonomialForm> {
    let mon_crs = if let Ok(crs_file_path) = std::env::var("CRS_FILE") {
        println!("using crs file at {crs_file_path}");
        let crs_file = std::fs::File::open(&crs_file_path).expect(&format!("crs file at {}", crs_file_path));
        let mon_crs = Crs::<Bn256, CrsForMonomialForm>::read(crs_file).expect(&format!("read crs file at {}", crs_file_path));
        assert!(domain_size <= mon_crs.g1_bases.len());

        mon_crs
    } else {
        Crs::<Bn256, CrsForMonomialForm>::crs_42(domain_size, &worker)
    };

    mon_crs
}

fn create_proof_for_circuit<C: Circuit<Bn256>, P: PlonkConstraintSystemParams<Bn256> + 'static, MG: MainGate<Bn256>, T: Transcript<Fr>>(
    circuit: C,
    out_dir: &str,
) -> (VerificationKey<Bn256, C>, Proof<Bn256, C>) {
    let worker = Worker::new();
    println!("Synthesizing circuit");
    let mut assembly = TrivialAssembly::<Bn256, P, MG>::new();
    circuit.synthesize(&mut assembly).unwrap();
    assert!(assembly.is_satisfied());
    assembly.finalize();

    let domain_size = assembly.n() + 1;
    assert!(domain_size.is_power_of_two());
    println!("Generating setup");
    let setup = assembly.create_setup::<C>(&worker).unwrap();

    let crs = init_crs(&worker, domain_size);
    println!("Generating Vk");
    let vk = VerificationKey::from_setup(&setup, &worker, &crs).expect("vk from setup");

    println!("Generating proof");
    let proof = assembly.create_proof::<C, T>(&worker, &setup, &crs, None).expect("proof");

    save_proof_and_vk_into_file(&proof, &vk, &out_dir);

    (vk, proof)
}

pub fn save_proof_and_vk_into_file<C: Circuit<Bn256>>(proof: &Proof<Bn256, C>, vk: &VerificationKey<Bn256, C>, output_path: &str) {
    let proof_file_path = format!("{}_proof.json", output_path);
    let proof_file = std::fs::File::create(&proof_file_path).unwrap();
    serde_json::to_writer(proof_file, &proof).unwrap();
    println!("proof saved at {proof_file_path}");

    let vk_file_path = format!("{}_vk.json", output_path);
    let vk_file = std::fs::File::create(&vk_file_path).unwrap();
    serde_json::to_writer(vk_file, &vk).unwrap();
    println!("vk saved at {vk_file_path}");
}
