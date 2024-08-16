use circuit_definitions::circuit_definitions::{
    aux_layer::{
        wrapper::ZkSyncCompressionWrapper, ZkSyncCompressionForWrapperCircuit, ZkSyncCompressionProofForWrapper, ZkSyncCompressionVerificationKey, ZkSyncCompressionVerificationKeyForWrapper,
    },
    recursion_layer::{ZkSyncRecursionLayerProof, ZkSyncRecursionLayerVerificationKey},
};
use circuit_definitions::snark_wrapper::franklin_crypto::bellman::{
    bn256::{Bn256, Fr},
    plonk::{
        better_better_cs::{
            cs::{Circuit, PlonkCsWidth3Params, SetupAssembly},
            gates::naive_main_gate::NaiveMainGate,
        },
        commitments::transcript::keccak_transcript::RollingKeccakTranscript,
    },
    worker::Worker,
};

use super::convenience::*;

#[test]
#[ignore]
fn test_fflonk_proof_verification() {
    let vk_file_path = std::env::var("FFLONK_VK_FILE").expect("fflonk vk file path");
    let vk_file = std::fs::File::open(&vk_file_path).expect(&format!("vk file at {}", vk_file_path));
    let vk: FflonkSnarkVerifierCircuitVK = serde_json::from_reader(&vk_file).expect("deserialize vk");

    let proof_file_path = std::env::var("FFLONK_PROOF_FILE").expect("fflonk proof file path");
    let proof_file = std::fs::File::open(&proof_file_path).expect(&format!("proof file at {}", proof_file_path));
    let proof: FflonkSnarkVerifierCircuitProof = serde_json::from_reader(&proof_file).expect("deserialize proof");

    let is_valid = crate::verifier::verify_flattened_proof::<_, _, RollingKeccakTranscript<Fr>>(&vk, &proof, None).expect("verify proof");
    assert!(is_valid, "fflonk proof is not corrects");
}

#[test]
#[ignore]
fn run_proof_compression_by_schedule() {
    let scheduler_vk_file = std::fs::File::open("./data/scheduler_recursive_vk.json").unwrap();
    let scheduler_vk: ZkSyncRecursionLayerVerificationKey = serde_json::from_reader(&scheduler_vk_file).unwrap();
    let scheduler_proof_file = std::fs::File::open("./data/scheduler_recursive_proof.json").unwrap();
    let scheduler_proof: ZkSyncRecursionLayerProof = serde_json::from_reader(&scheduler_proof_file).unwrap();

    for compression_schedule in [CompressionSchedule::simplest(), CompressionSchedule::simple(), CompressionSchedule::mid(), CompressionSchedule::hard()] {
        compress_proof(scheduler_proof.clone().into_inner(), scheduler_vk.clone().into_inner(), compression_schedule);
    }
}

#[test]
#[ignore]
fn dump_trace_size_of_compressed_proofs() {
    for compression_schedule in [CompressionSchedule::simplest(), CompressionSchedule::simple(), CompressionSchedule::mid(), CompressionSchedule::hard()] {
        let compression_schedule_name = compression_schedule.name();
        let compression_mode = compression_schedule.compression_steps.last().unwrap().clone() as u8;
        let compression_wrapper_mode = compression_mode as u8 + 1;
        println!("Compression schedule {compression_schedule_name}");
        println!("Compression mode {compression_mode}");
        println!("Compression wrapper mode {}", compression_wrapper_mode);
        let proof_file_path = format!(
            "./data/compression_schedule/{compression_schedule_name}/compression_wrapper_{}_proof.json",
            compression_wrapper_mode as u8
        );
        let proof_file = std::fs::File::open(&proof_file_path).unwrap();
        let proof: ZkSyncCompressionProofForWrapper = serde_json::from_reader(&proof_file).unwrap();
        println!("Proof config for compression wrapper mode {}: {:?}", compression_wrapper_mode, &proof.proof_config);

        let compression_wrapper_vk_file_path = format!("./data/compression_schedule/{compression_schedule_name}/compression_wrapper_{}_vk.json", compression_wrapper_mode as u8);
        let compression_wrapper_vk_file = std::fs::File::open(&compression_wrapper_vk_file_path).unwrap();
        let compression_wrapper_vk: ZkSyncCompressionVerificationKeyForWrapper = serde_json::from_reader(&compression_wrapper_vk_file).unwrap();
        let fixed_parameters = compression_wrapper_vk.fixed_parameters.clone();
        println!("CS Geometry for compression wrapper mode {}: {:?}", compression_wrapper_mode, &fixed_parameters.parameters,);

        let compression_wrapper_function = ZkSyncCompressionWrapper::from_numeric_circuit_type(compression_wrapper_mode);

        let circuit = FflonkSnarkVerifierCircuit {
            witness: Some(proof),
            vk: compression_wrapper_vk,
            fixed_parameters,
            transcript_params: (),
            wrapper_function: compression_wrapper_function,
        };
        let mut setup_assembly = SetupAssembly::<Bn256, PlonkCsWidth3Params, NaiveMainGate>::new();
        circuit.synthesize(&mut setup_assembly).unwrap();
        let raw_size = setup_assembly.n();
        setup_assembly.finalize();
        let domain_size = setup_assembly.n() + 1;
        assert!(domain_size.is_power_of_two());

        println!(
            "SNARK wrapper for compression mode {} has {} constraints {} padded  {} log",
            compression_wrapper_mode,
            raw_size,
            setup_assembly.n(),
            domain_size.trailing_zeros()
        );
    }
}

#[test]
#[ignore]
fn verify_compression_wrapper_proof() {
    let compression_wrapper_mode: u8 = std::env::var("COMPRESSION_WRAPPER_MODE").unwrap().parse::<u8>().unwrap();
    let vk_file_path = std::env::var("INPUT_VK_FILE").expect("input vk file path");
    let vk_file = std::fs::File::open(&vk_file_path).expect(&format!("vk file at {}", vk_file_path));
    let input_vk: ZkSyncCompressionVerificationKey = serde_json::from_reader(&vk_file).expect("deserialize vk");

    let proof_file_path = std::env::var("COMPRESSION_WRAPPER_PROOF_FILE").expect("compression wrapper proof file path");
    let proof_file = std::fs::File::open(&proof_file_path).expect(&format!("proof file at {}", proof_file_path));
    let proof: ZkSyncCompressionProofForWrapper = serde_json::from_reader(&proof_file).expect("deserialize proof");

    let vk_file_path = std::env::var("COMPRESSION_WRAPPER_VK_FILE").expect("input vk file path");
    let vk_file = std::fs::File::open(&vk_file_path).expect(&format!("vk file at {}", vk_file_path));
    let vk: ZkSyncCompressionVerificationKeyForWrapper = serde_json::from_reader(&vk_file).expect("deserialize vk");

    let compression_wrapper_circuit = ZkSyncCompressionForWrapperCircuit::from_witness_and_vk(None, input_vk, compression_wrapper_mode);
    let builder = compression_wrapper_circuit.into_dyn_verifier_builder();
    let verifier = builder.create_verifier();
    let is_proof_valid = match compression_wrapper_circuit {
        ZkSyncCompressionForWrapperCircuit::CompressionMode1Circuit(inner) => verify_compression_circuit(inner, &proof, &vk, verifier),
        ZkSyncCompressionForWrapperCircuit::CompressionMode2Circuit(inner) => verify_compression_circuit(inner, &proof, &vk, verifier),
        ZkSyncCompressionForWrapperCircuit::CompressionMode3Circuit(inner) => verify_compression_circuit(inner, &proof, &vk, verifier),
        ZkSyncCompressionForWrapperCircuit::CompressionMode4Circuit(inner) => verify_compression_circuit(inner, &proof, &vk, verifier),
        ZkSyncCompressionForWrapperCircuit::CompressionMode5Circuit(inner) => verify_compression_circuit(inner, &proof, &vk, verifier),
    };
    if is_proof_valid == false {
        println!("Proof is not correct")
    }
}

#[test]
#[ignore]
fn test_snark_circuit_with_naive_main_gate() {
    let path = format!("./data/compression_schedule/hard");
    let compression_wrapper_mode = if let Ok(compression_wrapper_mode) = std::env::var("COMPRESSION_WRAPPER_MODE") {
        compression_wrapper_mode.parse::<u8>().unwrap()
    } else {
        5u8
    };
    println!("Compression mode {}", compression_wrapper_mode);
    let compression_proof_file_path = if let Ok(file_path) = std::env::var("COMPRESSION_PROOF_FILE") {
        file_path
    } else {
        format!("{}/compression_wrapper_{compression_wrapper_mode}_proof.json", path)
    };
    println!("Reading proof file at {compression_proof_file_path}");
    let compression_vk_file_path = if let Ok(file_path) = std::env::var("COMPRESSION_VK_FILE") {
        file_path
    } else {
        format!("{}/compression_wrapper_{compression_wrapper_mode}_vk.json", path)
    };
    println!("Reading vk file at {compression_vk_file_path}");

    let compression_proof_file = std::fs::File::open(compression_proof_file_path).unwrap();
    let compression_proof: ZkSyncCompressionProofForWrapper = serde_json::from_reader(&compression_proof_file).unwrap();

    let compression_vk_file = std::fs::File::open(compression_vk_file_path).unwrap();
    let compression_vk: ZkSyncCompressionVerificationKeyForWrapper = serde_json::from_reader(&compression_vk_file).unwrap();

    let wrapper_function = ZkSyncCompressionWrapper::from_numeric_circuit_type(compression_wrapper_mode);
    let fixed_parameters = compression_vk.fixed_parameters.clone();

    let circuit = FflonkSnarkVerifierCircuit {
        witness: Some(compression_proof),
        vk: compression_vk,
        fixed_parameters,
        transcript_params: (),
        wrapper_function,
    };
    let worker = Worker::new();
    let (proof, vk) = prove_fflonk_snark_verifier_circuit_single_shot(&circuit, &worker);

    save_fflonk_proof_and_vk_into_file(&proof, &vk, &path);
}

#[test]
#[ignore]
fn full_e2e_test() {
    let start = std::time::Instant::now();
    // prepare circuit input of the first compression step:
    // - scheduler proof
    // - scheduler vk
    let scheduler_vk_file = std::fs::File::open("./data/scheduler_recursive_vk.json").unwrap();
    let scheduler_vk: ZkSyncRecursionLayerVerificationKey = serde_json::from_reader(&scheduler_vk_file).unwrap();
    let scheduler_proof_file = std::fs::File::open("./data/scheduler_recursive_proof.json").unwrap();
    let scheduler_proof: ZkSyncRecursionLayerProof = serde_json::from_reader(&scheduler_proof_file).unwrap();

    // set compression schedule:
    // - "hard" is the strategy that gives smallest final circuit
    let compression_schedule = CompressionSchedule::hard();
    let compression_schedule_name = compression_schedule.name();
    let compression_wrapper_mode = compression_schedule.compression_steps.last().unwrap().clone() as u8 + 1;
    // compress proof step by step: 1 -> 2 -> 3 -> 4 -> 5(wrapper)
    compress_proof(scheduler_proof.into_inner(), scheduler_vk.into_inner(), compression_schedule);
    let path = format!("./data/compression_schedule/{compression_schedule_name}");
    // compression done, take proof and vk of last step as input of the next step
    // - 5 (wrapper) -> fflonk proof
    let compression_wrapper_proof_file_path = format!("{}/compression_wrapper_{}_proof.json", &path, compression_wrapper_mode);
    let compression_wrapper_proof_file = std::fs::File::open(compression_wrapper_proof_file_path).unwrap();
    let compression_wrapper_proof = serde_json::from_reader(&compression_wrapper_proof_file).unwrap();

    let compression_wrapper_vk_file_path = format!("{}/compression_wrapper_{}_vk.json", &path, compression_wrapper_mode);
    let compression_wrapper_vk_file = std::fs::File::open(compression_wrapper_vk_file_path).unwrap();
    let compression_wrapper_vk: ZkSyncCompressionVerificationKeyForWrapper = serde_json::from_reader(&compression_wrapper_vk_file).unwrap();

    // construct fflonk snark verifier circuit
    let wrapper_function = ZkSyncCompressionWrapper::from_numeric_circuit_type(compression_wrapper_mode);
    let fixed_parameters = compression_wrapper_vk.fixed_parameters.clone();
    let circuit = FflonkSnarkVerifierCircuit {
        witness: Some(compression_wrapper_proof),
        vk: compression_wrapper_vk,
        fixed_parameters: fixed_parameters,
        transcript_params: (),
        wrapper_function,
    };
    // create fflonk proof in single shot - without precomputation
    let (proof, vk) = prove_fflonk_snark_verifier_circuit_single_shot(&circuit, &Worker::new());
    println!("Full e2e pipeline takes {} seconds", start.elapsed().as_secs());
    save_fflonk_proof_and_vk_into_file(&proof, &vk, &path);
}

#[test]
#[ignore]
fn test_generate_compression_precomputations() {
    // prepare circuit input of the first compression step:
    // - scheduler vk
    let scheduler_vk_file = std::fs::File::open("./data/scheduler_recursive_vk.json").unwrap();
    let scheduler_vk: ZkSyncRecursionLayerVerificationKey = serde_json::from_reader(&scheduler_vk_file).unwrap();
    let compression_schedule = CompressionSchedule::hard();
    let compression_schedule_name = compression_schedule.name();
    let compression_wrapper_mode = compression_schedule.compression_steps.last().unwrap().clone() as u8 + 1;
    let output_blob_path = format!("./data/compression_schedule/{compression_schedule_name}");
    let output_blob_path = std::path::Path::new(&output_blob_path);
    precompute_and_save_setup_for_compression_steps(scheduler_vk.into_inner(), compression_schedule, &output_blob_path);
    println!("Compression step precomputations are generated");

    let compression_wrapper_proof_file_path = format!("{}/compression_wrapper_{}_proof.json", output_blob_path.to_str().unwrap(), compression_wrapper_mode);
    let compression_wrapper_proof_file = std::fs::File::open(compression_wrapper_proof_file_path).unwrap();
    let compression_wrapper_proof = serde_json::from_reader(&compression_wrapper_proof_file).unwrap();

    let compression_wrapper_vk_file_path = format!("{}/compression_wrapper_{}_vk.json", output_blob_path.to_str().unwrap(), compression_wrapper_mode);
    let compression_wrapper_vk_file = std::fs::File::open(compression_wrapper_vk_file_path).unwrap();
    let compression_wrapper_vk: ZkSyncCompressionVerificationKeyForWrapper = serde_json::from_reader(&compression_wrapper_vk_file).unwrap();

    // construct fflonk snark verifier circuit
    let wrapper_function = ZkSyncCompressionWrapper::from_numeric_circuit_type(compression_wrapper_mode);
    let fixed_parameters = compression_wrapper_vk.fixed_parameters.clone();
    let circuit = FflonkSnarkVerifierCircuit {
        witness: Some(compression_wrapper_proof),
        vk: compression_wrapper_vk,
        fixed_parameters: fixed_parameters,
        transcript_params: (),
        wrapper_function,
    };

    precompute_and_save_setup_for_fflonk_snark_circuit(&circuit, &Worker::new(), &output_blob_path.to_string_lossy());
    println!("fflonk precomputation is generated");
}

#[test]
#[ignore]
fn test_full_precomputed_e2e() {
    let start = std::time::Instant::now();
    // prepare circuit input of the first compression step:
    // - scheduler proof
    // - scheduler vk
    let scheduler_vk_file = std::fs::File::open("./data/scheduler_recursive_vk.json").unwrap();
    let scheduler_vk: ZkSyncRecursionLayerVerificationKey = serde_json::from_reader(&scheduler_vk_file).unwrap();
    let scheduler_proof_file = std::fs::File::open("./data/scheduler_recursive_proof.json").unwrap();
    let scheduler_proof: ZkSyncRecursionLayerProof = serde_json::from_reader(&scheduler_proof_file).unwrap();

    // set compression schedule:
    // - "hard" is the strategy that gives smallest final circuit
    let compression_schedule = CompressionSchedule::hard();
    let compression_schedule_name = compression_schedule.name();
    let compression_wrapper_mode = compression_schedule.compression_steps.last().unwrap().clone() as u8 + 1;
    // compress proof step by step: 1 -> 2 -> 3 -> 4 -> 5(wrapper)
    let output_blob_path = format!("./data/compression_schedule/{compression_schedule_name}");
    let output_blob_path = std::path::Path::new(&output_blob_path);
    compress_proof_with_precomputations(scheduler_proof.into_inner(), scheduler_vk.into_inner(), compression_schedule, &output_blob_path);
    println!("Compression steps with precomputations takes {} seconds", start.elapsed().as_secs());

    // compression done, take proof and vk of last step as input of the next step
    // - 5 (wrapper) -> fflonk proof
    let compression_wrapper_proof_file_path = format!("{}/compression_wrapper_{}_proof.json", output_blob_path.to_str().unwrap(), compression_wrapper_mode);
    let compression_wrapper_proof_file = std::fs::File::open(compression_wrapper_proof_file_path).unwrap();
    let compression_wrapper_proof = serde_json::from_reader(&compression_wrapper_proof_file).unwrap();

    let compression_wrapper_vk_file_path = format!("{}/compression_wrapper_{}_vk.json", output_blob_path.to_str().unwrap(), compression_wrapper_mode);
    let compression_wrapper_vk_file = std::fs::File::open(compression_wrapper_vk_file_path).unwrap();
    let compression_wrapper_vk: ZkSyncCompressionVerificationKeyForWrapper = serde_json::from_reader(&compression_wrapper_vk_file).unwrap();

    // construct fflonk snark verifier circuit
    let wrapper_function = ZkSyncCompressionWrapper::from_numeric_circuit_type(compression_wrapper_mode);
    let fixed_parameters = compression_wrapper_vk.fixed_parameters.clone();
    let circuit = FflonkSnarkVerifierCircuit {
        witness: Some(compression_wrapper_proof),
        vk: compression_wrapper_vk,
        fixed_parameters: fixed_parameters,
        transcript_params: (),
        wrapper_function,
    };
    // create fflonk proof with precomputations
    let (setup, vk) = load_fflonk_setup_and_vk_from_file(&output_blob_path.to_string_lossy());
    let proof = prove_fflonk_snark_verifier_circuit_with_precomputation(&circuit, &setup, &vk, &Worker::new());
    println!("Full e2e with precomputations pipeline takes {} seconds", start.elapsed().as_secs());
    save_fflonk_proof_and_vk_into_file(&proof, &vk, output_blob_path.to_str().unwrap());
}

#[test]
#[ignore]
fn download_and_transform_crs() {
    // 0 run circuit
    // 1 get domain size(required number of points) by running circuit
    // 2 then download CRS over network and
    // 3 then transform into our format
    // 4 then save into file for future use

    // step 0
    let compression_wrapper_mode = if let Ok(compression_wrapper_mode) = std::env::var("COMPRESSION_WRAPPER_MODE") {
        compression_wrapper_mode.parse::<u8>().unwrap()
    } else {
        5u8
    };
    println!("Compression mode {}", compression_wrapper_mode);
    let compression_proof_file_path = if let Ok(file_path) = std::env::var("COMPRESSION_PROOF_FILE") {
        file_path
    } else {
        format!("./data/compression_schedule/hard/compression_wrapper_{compression_wrapper_mode}_proof.json")
    };
    println!("Reading proof file at {compression_proof_file_path}");
    let compression_vk_file_path = if let Ok(file_path) = std::env::var("COMPRESSION_VK_FILE") {
        file_path
    } else {
        format!("./data/compression_schedule/hard/compression_wrapper_{compression_wrapper_mode}_vk.json")
    };
    println!("Reading vk file at {compression_vk_file_path}");

    let compression_proof_file = std::fs::File::open(compression_proof_file_path).unwrap();
    let compression_proof: ZkSyncCompressionProofForWrapper = serde_json::from_reader(&compression_proof_file).unwrap();

    let compression_vk_file = std::fs::File::open(compression_vk_file_path).unwrap();
    let compression_vk: ZkSyncCompressionVerificationKeyForWrapper = serde_json::from_reader(&compression_vk_file).unwrap();

    let wrapper_function = ZkSyncCompressionWrapper::from_numeric_circuit_type(compression_wrapper_mode);
    let fixed_parameters = compression_vk.fixed_parameters.clone();

    let circuit = FflonkSnarkVerifierCircuit {
        witness: Some(compression_proof),
        vk: compression_vk,
        fixed_parameters,
        transcript_params: (),
        wrapper_function,
    };
    let mut assembly = SetupAssembly::<Bn256, PlonkCsWidth3Params, NaiveMainGate>::new();
    circuit.synthesize(&mut assembly).expect("must work");
    assert!(assembly.is_satisfied());
    println!("Trace length {} for compression mode {}", assembly.n(), compression_wrapper_mode);
    assembly.finalize();
    let domain_size = assembly.n() + 1;
    assert!(domain_size.is_power_of_two());
    println!(
        "Trace log length {} and length {} for compression mode {}",
        domain_size.trailing_zeros(),
        domain_size,
        compression_wrapper_mode
    );
    // Step 1
    let max_degree = crate::compute_max_combined_degree_from_assembly::<_, _, _, _, FflonkSnarkVerifierCircuit>(&assembly);
    // Step 2-4
    download_and_transform_ignition_transcripts(max_degree);
}
