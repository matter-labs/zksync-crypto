use std::path::PathBuf;

use franklin_crypto::bellman::{
    bn256::{Bn256, Fr},
    kate_commitment::{Crs, CrsForMonomialForm},
    plonk::{
        better_better_cs::{
            cs::{
                Assembly, Circuit, ConstraintSystem, Gate, GateInternal, LookupTableApplication,
                PlonkCsWidth4WithNextStepAndCustomGatesParams, PolyIdentifier, Setup,
                TrivialAssembly,
            },
            gates::{
                main_gate_with_d_next::Width4MainGateWithDNext,
                selector_optimized_with_d_next::SelectorOptimizedWidth4MainGateWithDNext,
            },
            proof::Proof,
            setup::VerificationKey,
            verifier,
        },
        commitments::transcript::{keccak_transcript::RollingKeccakTranscript, Transcript},
    },
    worker::Worker,
    Engine, Field, PrimeField, ScalarEngine, SynthesisError,
};

use crate::{circuits::*, generate};

fn generate_setup_vk_and_proof_for_std_main_gate<E: Engine, C: Circuit<E>, T: Transcript<E::Fr>>(
    circuit: &C,
    transcript_params: Option<T::InitializationParameters>,
    prefix: &str,
) -> Result<(), SynthesisError> {
    let worker = Worker::new();

    let mut cs = TrivialAssembly::<
        E,
        PlonkCsWidth4WithNextStepAndCustomGatesParams,
        Width4MainGateWithDNext,
    >::new();

    circuit.synthesize(&mut cs)?;
    cs.finalize();

    assert!(cs.is_satisfied());

    let setup = cs.create_setup::<C>(&worker)?;

    println!("domain size {}", setup.n);

    let domain_size = setup.n.clone();

    let crs = Crs::<E, CrsForMonomialForm>::crs_42(domain_size.next_power_of_two(), &worker);

    let vk = VerificationKey::from_setup(&setup, &worker, &crs)?;
    dbg!(vk.total_lookup_entries_length);

    let vk_file_name = format!("/tmp/{}_{}", prefix, "vk_keccak.key");
    let proof_file_name = format!("/tmp/{}_{}", prefix, "proof_keccak.proof");

    let mut vk_writer = std::fs::File::create(vk_file_name).expect("create vk file");
    vk.write(&mut vk_writer).expect("write vk into file");
    let proof = cs.create_proof::<_, T>(&worker, &setup, &crs, transcript_params.clone())?;

    let mut proof_writer = std::fs::File::create(proof_file_name).expect("create proof file");
    proof
        .write(&mut proof_writer)
        .expect("write proof into file");
    let verified = verifier::verify::<E, _, T>(&vk, &proof, transcript_params)?;

    assert!(verified, "proof verification failed");

    Ok(())
}
fn generate_setup_vk_and_proof_for_selector_optimized_main_gate<
    E: Engine,
    C: Circuit<E>,
    T: Transcript<E::Fr>,
>(
    circuit: &C,
    transcript_params: Option<T::InitializationParameters>,
    prefix: &str,
) -> Result<(), SynthesisError> {
    let worker = Worker::new();
    let mut cs = TrivialAssembly::<
        E,
        PlonkCsWidth4WithNextStepAndCustomGatesParams,
        SelectorOptimizedWidth4MainGateWithDNext,
    >::new();

    circuit.synthesize(&mut cs)?;
    cs.finalize();

    assert!(cs.is_satisfied());

    let setup = cs.create_setup::<C>(&worker)?;

    println!("domain size {}", setup.n);

    let domain_size = setup.n.clone();

    let crs = Crs::<E, CrsForMonomialForm>::crs_42(domain_size.next_power_of_two(), &worker);

    let vk = VerificationKey::from_setup(&setup, &worker, &crs)?;
    dbg!(vk.total_lookup_entries_length);
    // println!("{:#?}", vk);

    let vk_file_name = format!("/tmp/{}_{}", prefix, "vk_keccak.key");
    let proof_file_name = format!("/tmp/{}_{}", prefix, "proof_keccak.proof");

    let mut vk_writer = std::fs::File::create(vk_file_name).expect("create vk file");
    vk.write(&mut vk_writer).expect("write vk into file");
    let proof = cs.create_proof::<_, T>(&worker, &setup, &crs, transcript_params.clone())?;

    let mut proof_writer = std::fs::File::create(proof_file_name).expect("create proof file");
    proof
        .write(&mut proof_writer)
        .expect("write proof into file");
    let verified = verifier::verify::<E, _, T>(&vk, &proof, transcript_params)?;

    assert!(verified, "proof verification failed");

    Ok(())
}

#[test]
fn test_create_vk_and_proof_file_for_std_main_gate() {
    {
        let circuit = DummyCircuit;
        generate_setup_vk_and_proof_for_std_main_gate::<
            Bn256,
            _,
            RollingKeccakTranscript<<Bn256 as ScalarEngine>::Fr>,
        >(&circuit, None, "width4")
        .unwrap()
    }
    {
        let circuit = DummyCircuitWithLookup;
        generate_setup_vk_and_proof_for_std_main_gate::<
            Bn256,
            _,
            RollingKeccakTranscript<<Bn256 as ScalarEngine>::Fr>,
        >(&circuit, None, "width4_with_lookup")
        .unwrap()
    }
    {
        let circuit = DummyCircuitWithRescue;
        generate_setup_vk_and_proof_for_std_main_gate::<
            Bn256,
            _,
            RollingKeccakTranscript<<Bn256 as ScalarEngine>::Fr>,
        >(&circuit, None, "width4_with_rescue")
        .unwrap()
    }
    {
        let circuit = DummyCircuitWithLookupAndRescue;
        generate_setup_vk_and_proof_for_std_main_gate::<
            Bn256,
            _,
            RollingKeccakTranscript<<Bn256 as ScalarEngine>::Fr>,
        >(&circuit, None, "width4_with_lookup_and_rescue")
        .unwrap()
    }
}

#[test]
fn test_create_vk_and_proof_file_for_selector_optimized_main_gate() {
    {
        let circuit = SelectorOptimizedDummyCircuit;
        generate_setup_vk_and_proof_for_selector_optimized_main_gate::<
            Bn256,
            _,
            RollingKeccakTranscript<<Bn256 as ScalarEngine>::Fr>,
        >(&circuit, None, "selector_optimized")
        .unwrap()
    }
    {
        let circuit = SelectorOptimizedDummyCircuitWithLookup;
        generate_setup_vk_and_proof_for_selector_optimized_main_gate::<
            Bn256,
            _,
            RollingKeccakTranscript<<Bn256 as ScalarEngine>::Fr>,
        >(&circuit, None, "selector_optimized_with_lookup")
        .unwrap()
    }
    {
        let circuit = SelectorOptimizedDummyCircuitWithRescue;
        generate_setup_vk_and_proof_for_selector_optimized_main_gate::<
            Bn256,
            _,
            RollingKeccakTranscript<<Bn256 as ScalarEngine>::Fr>,
        >(&circuit, None, "selector_optimized_with_rescue")
        .unwrap()
    }
    {
        let circuit = SelectorOptimizedDummyCircuitWithLookupAndRescue;
        generate_setup_vk_and_proof_for_selector_optimized_main_gate::<
            Bn256,
            _,
            RollingKeccakTranscript<<Bn256 as ScalarEngine>::Fr>,
        >(
            &circuit,
            None,
            "selector_optimized_with_lookup_and_rescue",
        )
        .unwrap()
    }
}

#[test]
fn test_render_vk_with_default_main_gate() {
    generate(
        PathBuf::from("./../block_vk_20_keccak.key"),
        PathBuf::from("./../hardhat/contracts"),
        None,
    );
}

#[test]
fn test_generate_vk_and_proof_for_dummy_circuit_with_selector_optimized_main_gate() {
    generate(
        PathBuf::from("./../block_vk_20_keccak.key"),
        PathBuf::from("./../hardhat/contracts"),
        None,
    );
}

#[test]
fn test_verification_keys() {
    {
        let reader = std::fs::File::open("/tmp/width4_vk_keccak.key").expect("open file");
        let vk = VerificationKey::<Bn256, DummyCircuit>::read(&reader).unwrap();
        assert_eq!(vk.gate_selectors_commitments.len(), 0);
        assert_eq!(
            <DummyCircuit as Circuit<Bn256>>::declare_used_gates()
                .unwrap()
                .len(),
            1
        );
    }
    {
        let reader =
            std::fs::File::open("/tmp/width4_with_lookup_vk_keccak.key").expect("open file");
        let vk = VerificationKey::<Bn256, DummyCircuitWithLookup>::read(&reader).unwrap();
        assert_eq!(vk.gate_selectors_commitments.len(), 0);
        assert_eq!(
            <DummyCircuitWithLookup as Circuit<Bn256>>::declare_used_gates()
                .unwrap()
                .len(),
            1
        );
    }
    {
        let reader =
            std::fs::File::open("/tmp/width4_with_rescue_vk_keccak.key").expect("open file");
        let vk = VerificationKey::<Bn256, DummyCircuitWithRescue>::read(&reader).unwrap();
        assert_eq!(vk.gate_selectors_commitments.len(), 2);
        assert_eq!(
            <DummyCircuitWithRescue as Circuit<Bn256>>::declare_used_gates()
                .unwrap()
                .len(),
            2
        );
    }
    {
        let reader = std::fs::File::open("/tmp/width4_with_lookup_and_rescue_vk_keccak.key")
            .expect("open file");
        let vk = VerificationKey::<Bn256, DummyCircuitWithLookupAndRescue>::read(&reader).unwrap();
        assert_eq!(vk.gate_selectors_commitments.len(), 2);
        assert_eq!(
            <DummyCircuitWithLookupAndRescue as Circuit<Bn256>>::declare_used_gates()
                .unwrap()
                .len(),
            2
        );
    }
}
