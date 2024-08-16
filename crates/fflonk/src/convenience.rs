use super::*;
use bellman::plonk::better_better_cs::cs::ProvingAssembly;
use bellman::plonk::better_better_cs::{
    cs::{Circuit, PlonkCsWidth3Params, SetupAssembly, TrivialAssembly},
    gates::naive_main_gate::NaiveMainGate,
};
use byteorder::{BigEndian, ReadBytesExt};
use circuit_definitions::circuit_definitions::aux_layer::compression_modes::{CompressionTranscriptForWrapper, CompressionTreeHasherForWrapper};
use circuit_definitions::circuit_definitions::recursion_layer::{ZkSyncRecursionProof, ZkSyncRecursionVerificationKey};
use circuit_definitions::{
    circuit_definitions::aux_layer::{
        compression::{CompressionLayerCircuit, ProofCompressionFunction},
        ZkSyncCompressionForWrapperCircuit, ZkSyncCompressionLayerCircuit, ZkSyncCompressionProof, ZkSyncCompressionProofForWrapper, ZkSyncCompressionVerificationKey,
        ZkSyncCompressionVerificationKeyForWrapper, ZkSyncSnarkWrapperCircuitNoLookupCustomGate,
    },
    snark_wrapper::franklin_crypto::bellman::plonk::commitments::transcript::keccak_transcript::RollingKeccakTranscript,
};
use franklin_crypto::boojum::algebraic_props::round_function::AbsorptionModeOverwrite;
use franklin_crypto::boojum::algebraic_props::sponge::GoldilocksPoseidon2Sponge;
use franklin_crypto::boojum::config::{ProvingCSConfig, SetupCSConfig};
use franklin_crypto::boojum::cs::implementations::fast_serialization::MemcopySerializable;
use franklin_crypto::boojum::cs::implementations::hints::{DenseVariablesCopyHint, DenseWitnessCopyHint};
use franklin_crypto::boojum::cs::implementations::polynomial_storage::{SetupBaseStorage, SetupStorage};
use franklin_crypto::boojum::cs::implementations::transcript::GoldilocksPoisedon2Transcript;
use franklin_crypto::boojum::cs::oracle::merkle_tree::MerkleTreeWithCap;
use franklin_crypto::boojum::cs::oracle::TreeHasher;
use franklin_crypto::boojum::field::SmallField;
use franklin_crypto::boojum::{
    config::{CSConfig, DevCSConfig},
    cs::{
        cs_builder::new_builder,
        cs_builder_reference::CsReferenceImplementationBuilder,
        implementations::{prover::ProofConfig, reference_cs::CSReferenceAssembly, setup::FinalizationHintsForProver, verifier::Verifier},
    },
    field::goldilocks::{GoldilocksExt2, GoldilocksField},
};

use franklin_crypto::boojum::cs::implementations::proof::Proof as BoojumProof;
use franklin_crypto::boojum::cs::implementations::verifier::VerificationKey as BoojumVK;
use franklin_crypto::boojum::worker::Worker as BoojumWorker;

pub type FflonkSnarkVerifierCircuit = ZkSyncSnarkWrapperCircuitNoLookupCustomGate;
pub type FflonkSnarkVerifierCircuitVK = FflonkVerificationKey<Bn256, FflonkSnarkVerifierCircuit>;
pub type FflonkSnarkVerifierCircuitProof = FlattenedFflonkProof<Bn256, FflonkSnarkVerifierCircuit>;
pub type FflonkSnarkVerifierCircuitSetup = FflonkSetup<Bn256, FflonkSnarkVerifierCircuit>;

type CompressionTranscript = GoldilocksPoisedon2Transcript;
type CompressionTreeHasher = GoldilocksPoseidon2Sponge<AbsorptionModeOverwrite>;

pub const L1_VERIFIER_DOMAIN_SIZE_LOG: usize = 23;
type F = GoldilocksField;
type EXT = GoldilocksExt2;

pub fn init_crs(worker: &Worker, domain_size: usize, max_combined_degree: usize) -> Crs<Bn256, CrsForMonomialForm> {
    assert!(domain_size <= 1 << L1_VERIFIER_DOMAIN_SIZE_LOG);
    let mon_crs = if let Ok(crs_file_path) = std::env::var("CRS_FILE") {
        println!("using crs file at {crs_file_path}");
        let crs_file = std::fs::File::open(&crs_file_path).expect(&format!("crs file at {}", crs_file_path));
        let mon_crs = Crs::<Bn256, CrsForMonomialForm>::read(crs_file).expect(&format!("read crs file at {}", crs_file_path));
        assert!(max_combined_degree <= mon_crs.g1_bases.len());

        mon_crs
    } else {
        Crs::<Bn256, CrsForMonomialForm>::non_power_of_two_crs_42(max_combined_degree, &worker)
    };

    mon_crs
}

pub fn init_crs_for_vk() -> Crs<Bn256, CrsForMonomialForm> {
    let mon_crs = if let Ok(crs_file_path) = std::env::var("CRS_FILE") {
        println!("using crs file at {crs_file_path}");
        let crs_file = std::fs::File::open(&crs_file_path).expect(&format!("crs file at {}", crs_file_path));
        let mon_crs = Crs::<Bn256, CrsForMonomialForm>::read(crs_file).expect(&format!("read crs file at {}", crs_file_path));
        mon_crs
    } else {
        Crs::<Bn256, CrsForMonomialForm>::non_power_of_two_crs_42(1, &Worker::new())
    };

    mon_crs
}

pub fn precompute_and_save_setup_for_fflonk_snark_circuit(circuit: &FflonkSnarkVerifierCircuit, worker: &Worker, output_blob_path: &str) {
    let compression_wrapper_mode = circuit.wrapper_function.numeric_circuit_type();
    println!("Compression mode: {compression_wrapper_mode}");
    let mut setup_assembly = SetupAssembly::<Bn256, PlonkCsWidth3Params, NaiveMainGate>::new();
    circuit.synthesize(&mut setup_assembly).expect("must work");
    assert!(setup_assembly.is_satisfied());
    setup_assembly.finalize();
    println!("Finalized assembly contains {} gates", setup_assembly.n());
    let max_combined_degree = compute_max_combined_degree_from_assembly::<_, _, _, _, FflonkSnarkVerifierCircuit>(&setup_assembly);

    let domain_size = setup_assembly.n() + 1;
    assert!(domain_size.is_power_of_two());
    assert!(domain_size <= 1 << L1_VERIFIER_DOMAIN_SIZE_LOG);

    let mon_crs = init_crs(&worker, domain_size, max_combined_degree);
    let setup: FflonkSnarkVerifierCircuitSetup = FflonkSetup::create_setup(&setup_assembly, &worker, &mon_crs).expect("fflonk setup");
    let vk = FflonkVerificationKey::from_setup(&setup, &mon_crs).unwrap();

    save_fflonk_setup_and_vk_into_file(&setup, &vk, output_blob_path);
}

pub fn prove_fflonk_snark_verifier_circuit_with_precomputation(
    circuit: &FflonkSnarkVerifierCircuit,
    precomputed_setup: &FflonkSnarkVerifierCircuitSetup,
    vk: &FflonkSnarkVerifierCircuitVK,
    worker: &Worker,
) -> FflonkSnarkVerifierCircuitProof {
    let compression_wrapper_mode = circuit.wrapper_function.numeric_circuit_type();
    let mut assembly = ProvingAssembly::<Bn256, PlonkCsWidth3Params, NaiveMainGate>::new();
    circuit.synthesize(&mut assembly).expect("must work");
    assert!(assembly.is_satisfied());
    assembly.finalize();

    let domain_size = assembly.n() + 1;
    println!("Trace log length {} for compression mode {}", domain_size.trailing_zeros(), compression_wrapper_mode);
    assert!(domain_size.is_power_of_two());

    let max_combined_degree = compute_max_combined_degree_from_assembly::<_, _, _, _, FflonkSnarkVerifierCircuit>(&assembly);

    assert!(domain_size <= 1 << L1_VERIFIER_DOMAIN_SIZE_LOG);

    let mon_crs = init_crs(&worker, domain_size, max_combined_degree);

    let proof = crate::prover::create_proof::<_, FflonkSnarkVerifierCircuit, _, _, _, RollingKeccakTranscript<Fr>>(assembly, &worker, &precomputed_setup, &mon_crs, None).expect("proof");
    let valid = crate::verify::<_, _, RollingKeccakTranscript<Fr>>(&vk, &proof, None).unwrap();
    assert!(valid, "proof verification fails");
    let flattened_proof = proof.clone().flatten();
    compare_proof_vs_flattened_proof(&proof, &flattened_proof, &vk);
    let valid = crate::verify_flattened_proof::<_, _, RollingKeccakTranscript<Fr>>(&vk, &flattened_proof, None).unwrap();
    assert!(valid, "flattened proof verification fails");

    flattened_proof
}

pub fn prove_fflonk_snark_verifier_circuit_single_shot(circuit: &FflonkSnarkVerifierCircuit, worker: &Worker) -> (FflonkSnarkVerifierCircuitProof, FflonkSnarkVerifierCircuitVK) {
    let compression_wrapper_mode = circuit.wrapper_function.numeric_circuit_type();
    let mut assembly = TrivialAssembly::<Bn256, PlonkCsWidth3Params, NaiveMainGate>::new();
    circuit.synthesize(&mut assembly).expect("must work");
    assert!(assembly.is_satisfied());
    assembly.finalize();
    let domain_size = assembly.n() + 1;
    assert!(domain_size.is_power_of_two());
    assert!(domain_size <= 1 << L1_VERIFIER_DOMAIN_SIZE_LOG);
    println!("Trace log length {} for compression mode {}", domain_size.trailing_zeros(), compression_wrapper_mode);

    let max_combined_degree = compute_max_combined_degree_from_assembly::<_, _, _, _, FflonkSnarkVerifierCircuit>(&assembly);
    println!("Max degree is {}", max_combined_degree);
    let mon_crs = init_crs(&worker, domain_size, max_combined_degree);
    let setup = FflonkSetup::create_setup(&assembly, &worker, &mon_crs).expect("setup");
    let vk = FflonkVerificationKey::from_setup(&setup, &mon_crs).unwrap();

    let proof = crate::prover::create_proof::<_, FflonkSnarkVerifierCircuit, _, _, _, RollingKeccakTranscript<Fr>>(assembly, &worker, &setup, &mon_crs, None).expect("proof");
    let valid = crate::verify::<_, _, RollingKeccakTranscript<Fr>>(&vk, &proof, None).unwrap();
    assert!(valid, "proof verification fails");
    let flattened_proof = proof.flatten();
    let valid = crate::verify_flattened_proof::<_, _, RollingKeccakTranscript<Fr>>(&vk, &flattened_proof, None).unwrap();
    assert!(valid, "flattened proof verification fails");

    (flattened_proof, vk)
}

#[derive(Copy, Clone, Debug)]
pub enum CompressionMode {
    One = 1,
    Two = 2,
    Three = 3,
    Four = 4,
    Five = 5,
}

impl CompressionMode {
    pub fn from_compression_mode(compression_mode: u8) -> Self {
        match compression_mode {
            1 => CompressionMode::One,
            2 => CompressionMode::Two,
            3 => CompressionMode::Three,
            4 => CompressionMode::Four,
            5 => CompressionMode::Five,
            _ => unreachable!(),
        }
    }
}

#[derive(Debug)]
pub struct CompressionSchedule {
    name: &'static str,
    pub compression_steps: Vec<CompressionMode>,
}

impl CompressionSchedule {
    pub fn name(&self) -> &'static str {
        self.name
    }

    pub fn simplest() -> Self {
        CompressionSchedule {
            name: "simplest",
            compression_steps: vec![CompressionMode::One],
        }
    }
    pub fn simple() -> Self {
        CompressionSchedule {
            name: "simple",
            compression_steps: vec![CompressionMode::One, CompressionMode::Two],
        }
    }

    pub fn mid() -> Self {
        CompressionSchedule {
            name: "mid",
            compression_steps: vec![CompressionMode::One, CompressionMode::Two, CompressionMode::Three],
        }
    }
    pub fn hard() -> Self {
        CompressionSchedule {
            name: "hard",
            compression_steps: vec![CompressionMode::One, CompressionMode::Two, CompressionMode::Three, CompressionMode::Four],
        }
    }
}

#[derive(Clone)]
pub enum CompressionInput {
    RecursionLayer(Option<ZkSyncRecursionProof>, ZkSyncRecursionVerificationKey, CompressionMode),
    CompressionLayer(Option<ZkSyncCompressionProof>, ZkSyncCompressionVerificationKey, CompressionMode),
    CompressionWrapperLayer(Option<ZkSyncCompressionProof>, ZkSyncCompressionVerificationKey, CompressionMode),
}

impl CompressionInput {
    pub fn into_compression_circuit(self) -> ZkSyncCompressionLayerCircuit {
        match self {
            CompressionInput::RecursionLayer(proof, vk, compression_mode) => {
                assert_eq!(compression_mode as u8, 1);
                ZkSyncCompressionLayerCircuit::from_witness_and_vk(proof, vk, 1)
            }
            CompressionInput::CompressionLayer(proof, vk, compression_mode) => ZkSyncCompressionLayerCircuit::from_witness_and_vk(proof, vk, compression_mode as u8),
            CompressionInput::CompressionWrapperLayer(_, _, _) => {
                unreachable!()
            }
        }
    }

    pub fn into_compression_wrapper_circuit(self) -> ZkSyncCompressionForWrapperCircuit {
        match self {
            CompressionInput::RecursionLayer(_, _, _) => {
                unreachable!()
            }
            CompressionInput::CompressionLayer(_, _, _) => {
                unreachable!()
            }
            CompressionInput::CompressionWrapperLayer(proof, vk, compression_mode) => ZkSyncCompressionForWrapperCircuit::from_witness_and_vk(proof, vk, compression_mode as u8),
        }
    }
}

pub fn compress_proof(proof: ZkSyncRecursionProof, vk: ZkSyncRecursionVerificationKey, schedule: CompressionSchedule) {
    let worker = BoojumWorker::new();
    let mut input = CompressionInput::RecursionLayer(Some(proof), vk, CompressionMode::One);

    dbg!(&schedule);
    let CompressionSchedule {
        name: compression_schedule_name,
        compression_steps,
    } = schedule;

    let last_compression_wrapping_mode = CompressionMode::from_compression_mode(compression_steps.last().unwrap().clone() as u8 + 1);
    dbg!(&last_compression_wrapping_mode);

    /*
        This illustrates how compression enforced for the "hardest" strategy

           input                       compression     verifier          output        compression wrapper
       _____________________________   ____________    ___________     __________      ___________________
       scheduler       proof   vk          1           scheduler   ->  compressed1         compressed2
       compressed1     proof   vk          2           compressed1 ->  compressed2         compressed3
       compressed2     proof   vk          3           compressed2 ->  compressed3         compressed4
       compressed3     proof   vk          4           compressed3 ->  compressed4         compressed5


       compressed5     proof   vk          -       compression wrapper5       ->  fflonk proof
    */

    let num_compression_steps = compression_steps.len();
    let mut compression_modes_iter = compression_steps.into_iter();
    for step_idx in 0..num_compression_steps {
        let compression_mode = compression_modes_iter.next().unwrap();
        let proof_file_path = format!("./data/compression_schedule/{compression_schedule_name}/compression_{}_proof.json", compression_mode as u8);
        let proof_file_path = std::path::Path::new(&proof_file_path);
        let vk_file_path = format!("./data/compression_schedule/{compression_schedule_name}/compression_{}_vk.json", compression_mode as u8);
        let vk_file_path = std::path::Path::new(&vk_file_path);
        if proof_file_path.exists() && vk_file_path.exists() {
            println!("Compression {compression_schedule_name}/{} proof and vk already exist ignoring", compression_mode as u8);
            let proof_file = std::fs::File::open(proof_file_path).unwrap();
            let proof = serde_json::from_reader(&proof_file).unwrap();
            let vk_file = std::fs::File::open(vk_file_path).unwrap();
            let vk = serde_json::from_reader(&vk_file).unwrap();
            if step_idx + 1 == num_compression_steps {
                input = CompressionInput::CompressionWrapperLayer(proof, vk, last_compression_wrapping_mode)
            } else {
                input = CompressionInput::CompressionLayer(proof, vk, CompressionMode::from_compression_mode(compression_mode as u8 + 1))
            }

            continue;
        }
        let compression_circuit = input.into_compression_circuit();
        let circuit_type = compression_circuit.numeric_circuit_type();
        println!("Proving compression {compression_schedule_name}/{}", compression_mode as u8);
        let (proof, vk) = inner_prove_compression_layer_circuit(compression_circuit, &worker);
        println!("Proof for compression {compression_schedule_name}/{} is generated!", compression_mode as u8);

        save_compression_proof_and_vk_into_file(&proof, &vk, compression_schedule_name, circuit_type as u8);
        if step_idx + 1 == num_compression_steps {
            input = CompressionInput::CompressionWrapperLayer(Some(proof), vk, last_compression_wrapping_mode);
        } else {
            input = CompressionInput::CompressionLayer(Some(proof), vk, CompressionMode::from_compression_mode(compression_mode as u8 + 1));
        }
    }

    // last wrapping step
    let proof_file_path = format!(
        "./data/compression_schedule/{compression_schedule_name}/compression_wrapper_{}_proof.json",
        last_compression_wrapping_mode as u8
    );
    let proof_file_path = std::path::Path::new(&proof_file_path);
    let vk_file_path = format!(
        "./data/compression_schedule/{compression_schedule_name}/compression_wrapper_{}_vk.json",
        last_compression_wrapping_mode as u8
    );
    let vk_file_path = std::path::Path::new(&vk_file_path);
    println!("Compression for wrapper level {}", last_compression_wrapping_mode as u8);
    if proof_file_path.exists() && vk_file_path.exists() {
        println!(
            "Compression {compression_schedule_name}/{} for wrapper proof and vk already exist ignoring",
            last_compression_wrapping_mode as u8
        );
    } else {
        println!("Proving compression {compression_schedule_name}/{} for wrapper", last_compression_wrapping_mode as u8);
        let compression_circuit = input.into_compression_wrapper_circuit();
        let (proof, vk) = inner_prove_compression_wrapper_circuit(compression_circuit, &worker);
        println!("Proof for compression wrapper {compression_schedule_name}/{} is generated!", last_compression_wrapping_mode as u8);
        save_compression_wrapper_proof_and_vk_into_file(&proof, &vk, compression_schedule_name, last_compression_wrapping_mode as u8);
        println!("Compression wrapper proof and vk for {compression_schedule_name}/{} saved", last_compression_wrapping_mode as u8);
    }
}

pub fn compress_proof_with_precomputations(proof: ZkSyncRecursionProof, vk: ZkSyncRecursionVerificationKey, schedule: CompressionSchedule, blob_path: &std::path::Path) {
    let worker = BoojumWorker::new();
    let mut input = CompressionInput::RecursionLayer(Some(proof), vk, CompressionMode::One);

    dbg!(&schedule);
    let CompressionSchedule {
        name: compression_schedule_name,
        compression_steps,
    } = schedule;

    let last_compression_wrapping_mode = CompressionMode::from_compression_mode(compression_steps.last().unwrap().clone() as u8 + 1);
    dbg!(&last_compression_wrapping_mode);

    /*
        This illustrates how compression enforced for the "hardest" strategy

           input                       compression     verifier          output        compression wrapper
       _____________________________   ____________    ___________     __________      ___________________
       scheduler       proof   vk          1           scheduler   ->  compressed1         compressed2
       compressed1     proof   vk          2           compressed1 ->  compressed2         compressed3
       compressed2     proof   vk          3           compressed2 ->  compressed3         compressed4
       compressed3     proof   vk          4           compressed3 ->  compressed4         compressed5


       compressed5     proof   vk          -       compression wrapper5       ->  fflonk proof
    */

    let num_compression_steps = compression_steps.len();
    let mut compression_modes_iter = compression_steps.into_iter();
    for step_idx in 0..num_compression_steps {
        let compression_mode = compression_modes_iter.next().unwrap();
        let proof_file_path = format!("./data/compression_schedule/{compression_schedule_name}/compression_{}_proof.json", compression_mode as u8);
        let proof_file_path = std::path::Path::new(&proof_file_path);
        let vk_file_path = format!("./data/compression_schedule/{compression_schedule_name}/compression_{}_vk.json", compression_mode as u8);
        let vk_file_path = std::path::Path::new(&vk_file_path);
        if proof_file_path.exists() && vk_file_path.exists() {
            println!("Compression {compression_schedule_name}/{} proof and vk already exist ignoring", compression_mode as u8);
            let proof_file = std::fs::File::open(proof_file_path).unwrap();
            let proof = serde_json::from_reader(&proof_file).unwrap();
            let vk_file = std::fs::File::open(vk_file_path).unwrap();
            let vk = serde_json::from_reader(&vk_file).unwrap();
            if step_idx + 1 == num_compression_steps {
                input = CompressionInput::CompressionWrapperLayer(proof, vk, last_compression_wrapping_mode)
            } else {
                input = CompressionInput::CompressionLayer(proof, vk, CompressionMode::from_compression_mode(compression_mode as u8 + 1))
            }
            continue;
        }
        let compression_circuit = input.into_compression_circuit();
        assert_eq!(compression_circuit.numeric_circuit_type(), compression_mode as u8);
        println!("Proving compression {compression_schedule_name}/{}", compression_mode as u8);
        let proof_config = compression_circuit.proof_config_for_compression_step();

        let (proof, vk) = match compression_circuit {
            ZkSyncCompressionLayerCircuit::CompressionMode1Circuit(inner) => {
                prove_compression_circuit_with_precomputations(inner, proof_config, &worker, blob_path.to_str().unwrap(), compression_mode as u8)
            }
            ZkSyncCompressionLayerCircuit::CompressionMode2Circuit(inner) => {
                prove_compression_circuit_with_precomputations(inner, proof_config, &worker, blob_path.to_str().unwrap(), compression_mode as u8)
            }
            ZkSyncCompressionLayerCircuit::CompressionMode3Circuit(inner) => {
                prove_compression_circuit_with_precomputations(inner, proof_config, &worker, blob_path.to_str().unwrap(), compression_mode as u8)
            }
            ZkSyncCompressionLayerCircuit::CompressionMode4Circuit(inner) => {
                prove_compression_circuit_with_precomputations(inner, proof_config, &worker, blob_path.to_str().unwrap(), compression_mode as u8)
            }
            ZkSyncCompressionLayerCircuit::CompressionMode5Circuit(inner) => {
                prove_compression_circuit_with_precomputations(inner, proof_config, &worker, blob_path.to_str().unwrap(), compression_mode as u8)
            }
        };
        println!("Proof for compression {compression_schedule_name}/{} is generated!", compression_mode as u8);

        save_compression_proof_and_vk_into_file(&proof, &vk, compression_schedule_name, compression_mode as u8);
        if step_idx + 1 == num_compression_steps {
            input = CompressionInput::CompressionWrapperLayer(Some(proof), vk, last_compression_wrapping_mode);
        } else {
            input = CompressionInput::CompressionLayer(Some(proof), vk, CompressionMode::from_compression_mode(compression_mode as u8 + 1));
        }
    }

    // last wrapping step
    let proof_file_path = format!(
        "./data/compression_schedule/{compression_schedule_name}/compression_wrapper_{}_proof.json",
        last_compression_wrapping_mode as u8
    );
    let proof_file_path = std::path::Path::new(&proof_file_path);
    let vk_file_path = format!(
        "./data/compression_schedule/{compression_schedule_name}/compression_wrapper_{}_vk.json",
        last_compression_wrapping_mode as u8
    );
    let vk_file_path = std::path::Path::new(&vk_file_path);
    println!("Compression for wrapper level {}", last_compression_wrapping_mode as u8);
    if proof_file_path.exists() && vk_file_path.exists() {
        println!(
            "Compression {compression_schedule_name}/{} for wrapper proof and vk already exist ignoring",
            last_compression_wrapping_mode as u8
        );
    } else {
        println!("Proving compression {compression_schedule_name}/{} for wrapper", last_compression_wrapping_mode as u8);
        let compression_wrapper_circuit = input.into_compression_wrapper_circuit();
        assert_eq!(compression_wrapper_circuit.numeric_circuit_type(), last_compression_wrapping_mode as u8);
        let proof_config = compression_wrapper_circuit.proof_config_for_compression_step();
        let (proof, vk) = match compression_wrapper_circuit {
            ZkSyncCompressionForWrapperCircuit::CompressionMode1Circuit(inner) => {
                prove_compression_wrapper_circuit_with_precomputations(inner, proof_config, &worker, blob_path.to_str().unwrap(), last_compression_wrapping_mode as u8)
            }
            ZkSyncCompressionForWrapperCircuit::CompressionMode2Circuit(inner) => {
                prove_compression_wrapper_circuit_with_precomputations(inner, proof_config, &worker, blob_path.to_str().unwrap(), last_compression_wrapping_mode as u8)
            }
            ZkSyncCompressionForWrapperCircuit::CompressionMode3Circuit(inner) => {
                prove_compression_wrapper_circuit_with_precomputations(inner, proof_config, &worker, blob_path.to_str().unwrap(), last_compression_wrapping_mode as u8)
            }
            ZkSyncCompressionForWrapperCircuit::CompressionMode4Circuit(inner) => {
                prove_compression_wrapper_circuit_with_precomputations(inner, proof_config, &worker, blob_path.to_str().unwrap(), last_compression_wrapping_mode as u8)
            }
            ZkSyncCompressionForWrapperCircuit::CompressionMode5Circuit(inner) => {
                prove_compression_wrapper_circuit_with_precomputations(inner, proof_config, &worker, blob_path.to_str().unwrap(), last_compression_wrapping_mode as u8)
            }
        };
        println!("Proof for compression wrapper {compression_schedule_name}/{} is generated!", last_compression_wrapping_mode as u8);
        save_compression_wrapper_proof_and_vk_into_file(&proof, &vk, compression_schedule_name, last_compression_wrapping_mode as u8);
        println!("Compression wrapper proof and vk for {compression_schedule_name}/{} saved", last_compression_wrapping_mode as u8);
    }
}

pub fn precompute_and_save_setup_for_compression_steps(vk: ZkSyncRecursionVerificationKey, schedule: CompressionSchedule, output_blob_path: &std::path::Path) {
    let worker = BoojumWorker::new();
    let mut input = CompressionInput::RecursionLayer(None, vk, CompressionMode::One);

    let CompressionSchedule {
        name: compression_schedule_name,
        compression_steps,
    } = schedule;

    let last_compression_wrapping_mode = CompressionMode::from_compression_mode(compression_steps.last().unwrap().clone() as u8 + 1);

    let num_compression_steps = compression_steps.len();
    let mut compression_modes_iter = compression_steps.into_iter();
    for step_idx in 0..num_compression_steps {
        let compression_mode = compression_modes_iter.next().unwrap();
        let compression_circuit = input.into_compression_circuit();
        println!("Generating precomputed setup for compression layer {compression_schedule_name}/{}", compression_mode as u8);
        let (setup_base, setup, setup_tree, vars_hint, wits_hint, vk) = inner_precompute_setup_for_compression_layer(compression_circuit, &worker);
        save_precomputations_into_file(
            &setup_base,
            &setup,
            &setup_tree,
            &vk,
            &vars_hint,
            &wits_hint,
            output_blob_path.to_str().unwrap(),
            compression_mode as u8,
        );
        if step_idx + 1 == num_compression_steps {
            input = CompressionInput::CompressionWrapperLayer(None, vk, CompressionMode::from_compression_mode(last_compression_wrapping_mode as u8));
        } else {
            input = CompressionInput::CompressionLayer(None, vk, CompressionMode::from_compression_mode(step_idx as u8 + 1));
        }
    }

    // last wrapping step
    let compression_wrapper_circuit = input.into_compression_wrapper_circuit();
    assert_eq!(last_compression_wrapping_mode as u8, compression_wrapper_circuit.numeric_circuit_type());
    println!(
        "Generating precomputed setup for compression wrapper {compression_schedule_name}/{}",
        last_compression_wrapping_mode as u8
    );
    let (setup_base, setup, setup_tree, vars_hint, wits_hint, vk) = inner_precompute_setup_for_compression_wrapper(compression_wrapper_circuit, &worker);
    save_precomputations_into_file(
        &setup_base,
        &setup,
        &setup_tree,
        &vk,
        &vars_hint,
        &wits_hint,
        output_blob_path.to_str().unwrap(),
        last_compression_wrapping_mode as u8,
    );
}

pub fn inner_prove_compression_layer_circuit(circuit: ZkSyncCompressionLayerCircuit, worker: &BoojumWorker) -> (ZkSyncCompressionProof, ZkSyncCompressionVerificationKey) {
    let proof_config = circuit.proof_config_for_compression_step();
    let verifier_builder = circuit.into_dyn_verifier_builder();
    let verifier = verifier_builder.create_verifier();

    let (proof, vk, is_proof_valid) = match circuit {
        ZkSyncCompressionLayerCircuit::CompressionMode1Circuit(inner) => {
            let (proof, vk) = prove_compression_circuit(inner.clone(), proof_config, worker);
            let is_proof_valid = verify_compression_circuit(inner, &proof, &vk, verifier);
            (proof, vk, is_proof_valid)
        }
        ZkSyncCompressionLayerCircuit::CompressionMode2Circuit(inner) => {
            let (proof, vk) = prove_compression_circuit(inner.clone(), proof_config, worker);
            let is_proof_valid = verify_compression_circuit(inner, &proof, &vk, verifier);
            (proof, vk, is_proof_valid)
        }
        ZkSyncCompressionLayerCircuit::CompressionMode3Circuit(inner) => {
            let (proof, vk) = prove_compression_circuit(inner.clone(), proof_config, worker);
            let is_proof_valid = verify_compression_circuit(inner, &proof, &vk, verifier);
            (proof, vk, is_proof_valid)
        }
        ZkSyncCompressionLayerCircuit::CompressionMode4Circuit(inner) => {
            let (proof, vk) = prove_compression_circuit(inner.clone(), proof_config, worker);
            let is_proof_valid = verify_compression_circuit(inner, &proof, &vk, verifier);
            (proof, vk, is_proof_valid)
        }
        ZkSyncCompressionLayerCircuit::CompressionMode5Circuit(_inner) => {
            unreachable!("Only 4 modes of compression is allowed")
        }
    };
    if is_proof_valid == false {
        println!("Proof is invalid");
    }
    (proof, vk)
}

pub fn inner_prove_compression_wrapper_circuit(circuit: ZkSyncCompressionForWrapperCircuit, worker: &BoojumWorker) -> (ZkSyncCompressionProofForWrapper, ZkSyncCompressionVerificationKeyForWrapper) {
    let proof_config = circuit.proof_config_for_compression_step();
    let verifier_builder = circuit.into_dyn_verifier_builder();
    let verifier = verifier_builder.create_verifier();

    let (proof, vk, is_proof_valid) = match circuit {
        ZkSyncCompressionForWrapperCircuit::CompressionMode1Circuit(inner) => {
            let (proof, vk) = prove_compression_circuit(inner.clone(), proof_config, worker);
            let is_proof_valid = verify_compression_circuit(inner, &proof, &vk, verifier);
            (proof, vk, is_proof_valid)
        }
        ZkSyncCompressionForWrapperCircuit::CompressionMode2Circuit(inner) => {
            let (proof, vk) = prove_compression_circuit(inner.clone(), proof_config, worker);
            let is_proof_valid = verify_compression_circuit(inner, &proof, &vk, verifier);
            (proof, vk, is_proof_valid)
        }
        ZkSyncCompressionForWrapperCircuit::CompressionMode3Circuit(inner) => {
            let (proof, vk) = prove_compression_circuit(inner.clone(), proof_config, worker);
            let is_proof_valid = verify_compression_circuit(inner, &proof, &vk, verifier);
            (proof, vk, is_proof_valid)
        }
        ZkSyncCompressionForWrapperCircuit::CompressionMode4Circuit(inner) => {
            let (proof, vk) = prove_compression_circuit(inner.clone(), proof_config, worker);
            let is_proof_valid = verify_compression_circuit(inner, &proof, &vk, verifier);
            (proof, vk, is_proof_valid)
        }
        ZkSyncCompressionForWrapperCircuit::CompressionMode5Circuit(inner) => {
            let (proof, vk) = prove_compression_circuit(inner.clone(), proof_config, worker);
            let is_proof_valid = verify_compression_circuit(inner, &proof, &vk, verifier);
            (proof, vk, is_proof_valid)
        }
    };
    if is_proof_valid == false {
        println!("Proof is invalid");
    }

    (proof, vk)
}

pub fn synthesize_circuit<CF: ProofCompressionFunction, CS: CSConfig>(circuit: CompressionLayerCircuit<CF>) -> (FinalizationHintsForProver, CSReferenceAssembly<F, F, CS>) {
    let geometry = circuit.geometry();
    let (max_trace_len, num_vars) = circuit.size_hint();

    let builder_impl = CsReferenceImplementationBuilder::<GoldilocksField, F, CS>::new(geometry, max_trace_len.unwrap());
    let builder = new_builder::<_, GoldilocksField>(builder_impl);

    let builder = circuit.configure_builder_proxy(builder);
    let mut cs = builder.build(num_vars.unwrap());
    circuit.add_tables(&mut cs);
    circuit.synthesize_into_cs(&mut cs);
    let (domain_size, finalization_hint) = cs.pad_and_shrink();
    let cs = cs.into_assembly::<std::alloc::Global>();

    (finalization_hint, cs)
}

pub fn prove_compression_circuit<CF: ProofCompressionFunction>(
    circuit: CompressionLayerCircuit<CF>,
    proof_config: ProofConfig,
    worker: &BoojumWorker,
) -> (BoojumProof<F, CF::ThisLayerHasher, EXT>, BoojumVK<GoldilocksField, <CF as ProofCompressionFunction>::ThisLayerHasher>) {
    let (_, cs) = synthesize_circuit::<_, DevCSConfig>(circuit);
    cs.prove_one_shot::<_, CF::ThisLayerTranscript, CF::ThisLayerHasher, CF::ThisLayerPoW>(worker, proof_config, CF::this_layer_transcript_parameters())
}

fn inner_precompute_setup_for_compression_layer(
    circuit: ZkSyncCompressionLayerCircuit,
    worker: &BoojumWorker,
) -> (
    SetupBaseStorage<F, F>,
    SetupStorage<F, F>,
    MerkleTreeWithCap<F, CompressionTreeHasher>,
    DenseVariablesCopyHint,
    DenseWitnessCopyHint,
    BoojumVK<F, CompressionTreeHasher>,
) {
    let proof_config = circuit.proof_config_for_compression_step();
    match circuit {
        ZkSyncCompressionLayerCircuit::CompressionMode1Circuit(inner) => {
            let (_, cs) = synthesize_circuit::<_, SetupCSConfig>(inner.clone());
            cs.prepare_base_setup_with_precomputations_and_vk::<CompressionTranscript, CompressionTreeHasher>(proof_config, worker)
        }
        ZkSyncCompressionLayerCircuit::CompressionMode2Circuit(inner) => {
            let (_, cs) = synthesize_circuit::<_, SetupCSConfig>(inner.clone());
            cs.prepare_base_setup_with_precomputations_and_vk::<CompressionTranscript, CompressionTreeHasher>(proof_config, worker)
        }
        ZkSyncCompressionLayerCircuit::CompressionMode3Circuit(inner) => {
            let (_, cs) = synthesize_circuit::<_, SetupCSConfig>(inner.clone());
            cs.prepare_base_setup_with_precomputations_and_vk::<CompressionTranscript, CompressionTreeHasher>(proof_config, worker)
        }
        ZkSyncCompressionLayerCircuit::CompressionMode4Circuit(inner) => {
            let (_, cs) = synthesize_circuit::<_, SetupCSConfig>(inner.clone());
            cs.prepare_base_setup_with_precomputations_and_vk::<CompressionTranscript, CompressionTreeHasher>(proof_config, worker)
        }
        ZkSyncCompressionLayerCircuit::CompressionMode5Circuit(_inner) => {
            unreachable!("Only 4 modes of compression is allowed")
        }
    }
}

fn inner_precompute_setup_for_compression_wrapper(
    circuit: ZkSyncCompressionForWrapperCircuit,
    worker: &BoojumWorker,
) -> (
    SetupBaseStorage<F, F>,
    SetupStorage<F, F>,
    MerkleTreeWithCap<F, CompressionTreeHasherForWrapper>,
    DenseVariablesCopyHint,
    DenseWitnessCopyHint,
    BoojumVK<F, CompressionTreeHasherForWrapper>,
) {
    let proof_config = circuit.proof_config_for_compression_step();
    match circuit {
        ZkSyncCompressionForWrapperCircuit::CompressionMode1Circuit(inner) => {
            let (_, cs) = synthesize_circuit::<_, SetupCSConfig>(inner.clone());
            cs.prepare_base_setup_with_precomputations_and_vk::<CompressionTranscriptForWrapper, CompressionTreeHasherForWrapper>(proof_config, worker)
        }
        ZkSyncCompressionForWrapperCircuit::CompressionMode2Circuit(inner) => {
            let (_, cs) = synthesize_circuit::<_, SetupCSConfig>(inner.clone());
            cs.prepare_base_setup_with_precomputations_and_vk::<CompressionTranscriptForWrapper, CompressionTreeHasherForWrapper>(proof_config, worker)
        }
        ZkSyncCompressionForWrapperCircuit::CompressionMode3Circuit(inner) => {
            let (_, cs) = synthesize_circuit::<_, SetupCSConfig>(inner.clone());
            cs.prepare_base_setup_with_precomputations_and_vk::<CompressionTranscriptForWrapper, CompressionTreeHasherForWrapper>(proof_config, worker)
        }
        ZkSyncCompressionForWrapperCircuit::CompressionMode4Circuit(inner) => {
            let (_, cs) = synthesize_circuit::<_, SetupCSConfig>(inner.clone());
            cs.prepare_base_setup_with_precomputations_and_vk::<CompressionTranscriptForWrapper, CompressionTreeHasherForWrapper>(proof_config, worker)
        }
        ZkSyncCompressionForWrapperCircuit::CompressionMode5Circuit(_inner) => {
            unreachable!("Only 4 modes of compression is allowed")
        }
    }
}

pub fn load_precomputations_from_file<H: TreeHasher<F>>(
    blob_path: &str,
    compression_mode: u8,
) -> (
    SetupBaseStorage<F, F>,
    SetupStorage<F, F>,
    MerkleTreeWithCap<F, H>,
    DenseVariablesCopyHint,
    DenseWitnessCopyHint,
    BoojumVK<F, H>,
)
where
    H: TreeHasher<F>,
    H::Output: serde::de::DeserializeOwned,
{
    let setup_base_file_path = format!("{}/setup_base_{}.setup", blob_path, compression_mode as u8);
    println!("Loading setup base from file at {setup_base_file_path}");
    let setup_base_file = std::fs::File::open(setup_base_file_path).unwrap();
    let setup_base = SetupBaseStorage::read_from_buffer(&setup_base_file).unwrap();

    let setup_file_path = format!("{}/setup_{}.setup", blob_path, compression_mode as u8);
    println!("Loading setup from file at {setup_file_path}");
    let setup_file = std::fs::File::open(setup_file_path).unwrap();
    let setup = SetupStorage::read_from_buffer(&setup_file).unwrap();

    let setup_tree_file_path = format!("{}/setup_tree_{}.setup", blob_path, compression_mode as u8);
    println!("Loading setup tree from file at {setup_tree_file_path}");
    let setup_tree_file = std::fs::File::open(setup_tree_file_path).unwrap();
    // setup_tree.write_into_buffer(&setup_base_file).unwrap();
    let setup_tree = serde_json::from_reader(setup_tree_file).unwrap();

    let vars_hint_file_path = format!("{}/vars_hint{}.setup", blob_path, compression_mode as u8);
    println!("Loading variables hint from file at {vars_hint_file_path}");
    let vars_hint_file = std::fs::File::open(vars_hint_file_path).unwrap();
    let vars_hint = DenseVariablesCopyHint::read_from_buffer(&vars_hint_file).unwrap();

    let wits_hint_file_path = format!("{}/wits_hint{}.setup", blob_path, compression_mode as u8);
    println!("Loading witnesses hint from file at {wits_hint_file_path}");
    let wits_hint_file = std::fs::File::open(wits_hint_file_path).unwrap();
    let wits_hint = DenseWitnessCopyHint::read_from_buffer(&wits_hint_file).unwrap();

    let vk_file_path = format!("{}/setup_compression_{}.setup", blob_path, compression_mode as u8);
    println!("Loading vk from file at {vk_file_path}");
    let vk_file = std::fs::File::open(vk_file_path).unwrap();
    let vk = serde_json::from_reader(&vk_file).unwrap();

    (setup_base, setup, setup_tree, vars_hint, wits_hint, vk)
}

pub fn save_precomputations_into_file<F: SmallField, H>(
    setup_base: &SetupBaseStorage<F, F>,
    setup: &SetupStorage<F, F>,
    setup_tree: &MerkleTreeWithCap<F, H>,
    vk: &BoojumVK<F, H>,
    vars_hint: &DenseVariablesCopyHint,
    wits_hint: &DenseWitnessCopyHint,
    output_blob_path: &str,
    compression_mode: u8,
) where
    H: TreeHasher<F>,
    // TODO H::Output: serde::Serialize + MemcopySerializable,
    H::Output: serde::Serialize,
{
    let setup_base_file_path = format!("{}/setup_base_{}.setup", output_blob_path, compression_mode as u8);
    println!("Saving setup base into file at {setup_base_file_path}");
    let setup_base_file = std::fs::File::create(setup_base_file_path).unwrap();
    setup_base.write_into_buffer(&setup_base_file).unwrap();

    let setup_file_path = format!("{}/setup_{}.setup", output_blob_path, compression_mode as u8);
    println!("Saving setup into file at {setup_file_path}");
    let setup_file = std::fs::File::create(setup_file_path).unwrap();
    setup.write_into_buffer(&setup_file).unwrap();

    let setup_tree_file_path = format!("{}/setup_tree_{}.setup", output_blob_path, compression_mode as u8);
    println!("Saving setup tree into file at {setup_tree_file_path}");
    let setup_tree_file = std::fs::File::create(setup_tree_file_path).unwrap();
    // setup_tree.write_into_buffer(&setup_base_file).unwrap();
    serde_json::to_writer(setup_tree_file, setup_tree).unwrap();

    let vars_hint_file_path = format!("{}/vars_hint{}.setup", output_blob_path, compression_mode as u8);
    println!("Saving variables hint into file at {vars_hint_file_path}");
    let vars_hint_file = std::fs::File::create(vars_hint_file_path).unwrap();
    vars_hint.write_into_buffer(&vars_hint_file).unwrap();

    let wits_hint_file_path = format!("{}/wits_hint{}.setup", output_blob_path, compression_mode as u8);
    println!("Saving witnesses hint into file at {wits_hint_file_path}");
    let wits_hint_file = std::fs::File::create(wits_hint_file_path).unwrap();
    wits_hint.write_into_buffer(&wits_hint_file).unwrap();

    let vk_file_path = format!("{}/setup_compression_{}.setup", output_blob_path, compression_mode as u8);
    println!("Saving vk from into at {vk_file_path}");
    let vk_file = std::fs::File::create(vk_file_path).unwrap();
    serde_json::to_writer(&vk_file, vk).unwrap();
}

// this function is generic over the compression function for the PoW
pub fn prove_compression_circuit_with_precomputations<CF: ProofCompressionFunction>(
    circuit: CompressionLayerCircuit<CF>,
    proof_config: ProofConfig,
    worker: &BoojumWorker,
    blob_path: &str,
    compression_mode: u8,
) -> (BoojumProof<F, CompressionTreeHasher, EXT>, BoojumVK<GoldilocksField, CompressionTreeHasher>) {
    let (setup_base, setup, setup_tree, vars_hint, wits_hint, vk) = load_precomputations_from_file::<CompressionTreeHasher>(blob_path, compression_mode);
    let (_, cs) = synthesize_circuit::<_, ProvingCSConfig>(circuit);
    let proof =
        cs.prove_from_precomputations::<_, CompressionTranscript, CompressionTreeHasher, CF::ThisLayerPoW>(proof_config, &setup_base, &setup, &setup_tree, &vk, &vars_hint, &wits_hint, (), worker);

    (proof, vk)
}

pub fn prove_compression_wrapper_circuit_with_precomputations<CF: ProofCompressionFunction>(
    circuit: CompressionLayerCircuit<CF>,
    proof_config: ProofConfig,
    worker: &BoojumWorker,
    blob_path: &str,
    compression_mode: u8,
) -> (BoojumProof<F, CompressionTreeHasherForWrapper, EXT>, BoojumVK<GoldilocksField, CompressionTreeHasherForWrapper>) {
    let (setup_base, setup, setup_tree, vars_hint, wits_hint, vk) = load_precomputations_from_file::<CompressionTreeHasherForWrapper>(blob_path, compression_mode);
    let (_, cs) = synthesize_circuit::<_, ProvingCSConfig>(circuit);
    let proof = cs.prove_from_precomputations::<_, CompressionTranscriptForWrapper, CompressionTreeHasherForWrapper, CF::ThisLayerPoW>(
        proof_config,
        &setup_base,
        &setup,
        &setup_tree,
        &vk,
        &vars_hint,
        &wits_hint,
        (),
        worker,
    );

    (proof, vk)
}

pub fn verify_compression_circuit<CF: ProofCompressionFunction>(
    _circuit: CompressionLayerCircuit<CF>,
    proof: &BoojumProof<F, CF::ThisLayerHasher, EXT>,
    vk: &BoojumVK<GoldilocksField, <CF as ProofCompressionFunction>::ThisLayerHasher>,
    verifier: Verifier<GoldilocksField, EXT>,
) -> bool {
    verifier.verify::<CF::ThisLayerHasher, CF::ThisLayerTranscript, CF::ThisLayerPoW>(CF::this_layer_transcript_parameters(), vk, proof)
}

pub fn save_compression_proof_and_vk_into_file(proof: &ZkSyncCompressionProof, vk: &ZkSyncCompressionVerificationKey, compression_schedule_name: &'static str, compression_mode: u8) {
    let proof_file = std::fs::File::create(&format!("./data/compression_schedule/{compression_schedule_name}/compression_{}_proof.json", compression_mode)).unwrap();
    serde_json::to_writer(proof_file, &proof).unwrap();
    let vk_file = std::fs::File::create(&format!("./data/compression_schedule/{compression_schedule_name}/compression_{}_vk.json", compression_mode)).unwrap();
    serde_json::to_writer(vk_file, &vk).unwrap();
}

pub fn save_compression_wrapper_proof_and_vk_into_file(
    proof: &ZkSyncCompressionProofForWrapper,
    vk: &ZkSyncCompressionVerificationKeyForWrapper,
    compression_schedule_name: &'static str,
    compression_mode: u8,
) {
    let proof_file = std::fs::File::create(&format!("./data/compression_schedule/{compression_schedule_name}/compression_wrapper_{}_proof.json", compression_mode)).unwrap();
    serde_json::to_writer(proof_file, &proof).unwrap();
    let vk_file = std::fs::File::create(&format!("./data/compression_schedule/{compression_schedule_name}/compression_wrapper_{}_vk.json", compression_mode)).unwrap();
    serde_json::to_writer(vk_file, &vk).unwrap();
}

pub fn save_fflonk_proof_and_vk_into_file(proof: &FflonkSnarkVerifierCircuitProof, vk: &FflonkSnarkVerifierCircuitVK, output_blob_path: &str) {
    let proof_file_path = format!("{}/final_proof.json", output_blob_path);
    let proof_file = std::fs::File::create(&proof_file_path).unwrap();
    serde_json::to_writer(proof_file, &proof).unwrap();
    println!("proof saved at {proof_file_path}");
    let proof_file_path = format!("{}/final_proof_hex.json", output_blob_path);
    let hex_proof_file = std::fs::File::create(&proof_file_path).unwrap();
    proof.serialize_into_evm_format(hex_proof_file).unwrap();
    println!("evm proof saved at {proof_file_path}");
    let vk_file_path = format!("{}/final_vk.json", output_blob_path);
    let vk_file = std::fs::File::create(&vk_file_path).unwrap();
    serde_json::to_writer(vk_file, &vk).unwrap();
    println!("vk saved at {vk_file_path}");
}

pub fn save_fflonk_setup_and_vk_into_file(setup: &FflonkSnarkVerifierCircuitSetup, vk: &FflonkSnarkVerifierCircuitVK, output_blob_path: &str) {
    let setup_file_path = format!("{}/fflonk_snark_setup.json", output_blob_path);
    let setup_file = std::fs::File::create(&setup_file_path).unwrap();
    setup.write(&setup_file).unwrap();
    println!("fflonk precomputed setup saved into file at {setup_file_path}");
    let vk_file_path = format!("{}/fflonk_snark_setup.json", output_blob_path);
    let vk_file = std::fs::File::create(&vk_file_path).unwrap();
    vk.write(&vk_file).unwrap();
    println!("fflonk VK saved into file at {vk_file_path}");
}

pub fn load_fflonk_setup_and_vk_from_file(output_blob_path: &str) -> (FflonkSnarkVerifierCircuitSetup, FflonkSnarkVerifierCircuitVK) {
    let setup_file_path = format!("{}/fflonk_snark_setup.json", output_blob_path);
    println!("reading fflonk precomputed setup from file at {setup_file_path}");
    let setup_file = std::fs::File::open(&setup_file_path).unwrap();
    let setup = FflonkSnarkVerifierCircuitSetup::read(&setup_file).unwrap();
    let vk_file_path = format!("{}/fflonk_snark_setup.json", output_blob_path);
    println!("reading fflonk VK saved from file at {vk_file_path}");
    let vk_file = std::fs::File::open(&vk_file_path).unwrap();
    let vk = FflonkSnarkVerifierCircuitVK::read(&vk_file).unwrap();

    (setup, vk)
}
fn compare_proof_vs_flattened_proof<E: Engine, C: Circuit<E>>(expected_proof: &FflonkProof<E, C>, actual_proof: &FlattenedFflonkProof<E, C>, vk: &FflonkVerificationKey<E, C>) {
    let actual_proof = actual_proof.into_original_proof(vk);
    assert_eq!(expected_proof.n, actual_proof.n);
    assert_eq!(expected_proof.inputs, actual_proof.inputs);
    assert_eq!(expected_proof.c1, actual_proof.c1);
    assert_eq!(expected_proof.c2, actual_proof.c2);
    assert_eq!(expected_proof.w, actual_proof.w);
    assert_eq!(expected_proof.w_prime, actual_proof.w_prime);
    assert_eq!(expected_proof.setup_evaluations.gate_setups_at_z, actual_proof.setup_evaluations.gate_setups_at_z);
    assert_eq!(expected_proof.setup_evaluations.gate_selectors_at_z, actual_proof.setup_evaluations.gate_selectors_at_z);
    assert_eq!(expected_proof.setup_evaluations.permutations_at_z, actual_proof.setup_evaluations.permutations_at_z);
    assert_eq!(expected_proof.setup_evaluations.lookup_selector_at_z, actual_proof.setup_evaluations.lookup_selector_at_z);
    assert_eq!(expected_proof.setup_evaluations.lookup_tables_at_z, actual_proof.setup_evaluations.lookup_tables_at_z);
    assert_eq!(expected_proof.setup_evaluations.lookup_tables_at_z_omega, actual_proof.setup_evaluations.lookup_tables_at_z_omega);
    assert_eq!(expected_proof.setup_evaluations.lookup_table_type_at_z, actual_proof.setup_evaluations.lookup_table_type_at_z);

    assert_eq!(
        expected_proof.first_round_evaluations.trace_and_gate_evaluations.trace_evaluations_at_z,
        actual_proof.first_round_evaluations.trace_and_gate_evaluations.trace_evaluations_at_z
    );
    assert_eq!(
        expected_proof.first_round_evaluations.trace_and_gate_evaluations.trace_evaluations_at_z_omega,
        actual_proof.first_round_evaluations.trace_and_gate_evaluations.trace_evaluations_at_z_omega
    );
    assert_eq!(
        expected_proof.first_round_evaluations.trace_and_gate_evaluations.main_gate_quotient_at_z,
        actual_proof.first_round_evaluations.trace_and_gate_evaluations.main_gate_quotient_at_z
    );
    assert_eq!(
        expected_proof.first_round_evaluations.trace_and_gate_evaluations.main_gate_quotient_at_z_omega,
        actual_proof.first_round_evaluations.trace_and_gate_evaluations.main_gate_quotient_at_z_omega
    );
    assert_eq!(
        expected_proof.first_round_evaluations.trace_and_gate_evaluations.custom_gate_quotient_at_z,
        actual_proof.first_round_evaluations.trace_and_gate_evaluations.custom_gate_quotient_at_z
    );
    assert_eq!(
        expected_proof.first_round_evaluations.trace_and_gate_evaluations.custom_gate_quotient_at_z_omega,
        actual_proof.first_round_evaluations.trace_and_gate_evaluations.custom_gate_quotient_at_z_omega
    );

    assert_eq!(
        expected_proof.second_round_evaluations.copy_permutation_evaluations.grand_product_at_z,
        actual_proof.second_round_evaluations.copy_permutation_evaluations.grand_product_at_z
    );
    assert_eq!(
        expected_proof.second_round_evaluations.copy_permutation_evaluations.grand_product_at_z_omega,
        actual_proof.second_round_evaluations.copy_permutation_evaluations.grand_product_at_z_omega
    );
    assert_eq!(
        expected_proof.second_round_evaluations.copy_permutation_evaluations.first_quotient_at_z,
        actual_proof.second_round_evaluations.copy_permutation_evaluations.first_quotient_at_z
    );
    assert_eq!(
        expected_proof.second_round_evaluations.copy_permutation_evaluations.first_quotient_at_z_omega,
        actual_proof.second_round_evaluations.copy_permutation_evaluations.first_quotient_at_z_omega
    );
    assert_eq!(
        expected_proof.second_round_evaluations.copy_permutation_evaluations.second_quotient_at_z,
        actual_proof.second_round_evaluations.copy_permutation_evaluations.second_quotient_at_z
    );
    assert_eq!(
        expected_proof.second_round_evaluations.copy_permutation_evaluations.second_quotient_at_z_omega,
        actual_proof.second_round_evaluations.copy_permutation_evaluations.second_quotient_at_z_omega
    );

    if let (Some(actual_lookup), Some(expected_lookup)) = (
        expected_proof.second_round_evaluations.lookup_evaluations.as_ref(),
        actual_proof.second_round_evaluations.lookup_evaluations.as_ref(),
    ) {
        assert_eq!(expected_lookup.s_poly_at_z, actual_lookup.s_poly_at_z);
        assert_eq!(expected_lookup.s_poly_at_z_omega, actual_lookup.s_poly_at_z_omega);
        assert_eq!(expected_lookup.grand_product_at_z, actual_lookup.grand_product_at_z);
        assert_eq!(expected_lookup.grand_product_at_z_omega, actual_lookup.grand_product_at_z_omega);
        assert_eq!(expected_lookup.first_quotient_at_z, actual_lookup.first_quotient_at_z);
        assert_eq!(expected_lookup.first_quotient_at_z_omega, actual_lookup.first_quotient_at_z_omega);
        assert_eq!(expected_lookup.second_quotient_at_z, actual_lookup.second_quotient_at_z);
        assert_eq!(expected_lookup.second_quotient_at_z_omega, actual_lookup.second_quotient_at_z_omega);
        assert_eq!(expected_lookup.third_quotient_at_z, actual_lookup.third_quotient_at_z);
        assert_eq!(expected_lookup.third_quotient_at_z_omega, actual_lookup.third_quotient_at_z_omega);
    }
}

fn make_crs_from_ignition_transcript<S: AsRef<std::ffi::OsStr> + ?Sized>(path: &S, num_chunks: usize) -> Result<Crs<bellman::pairing::bn256::Bn256, CrsForMonomialForm>, SynthesisError> {
    use bellman::pairing::bn256::{Fq, Fq2};
    use std::io::BufRead;
    use std::io::Read;

    let chunk_size = 5_040_000;

    let base_path = std::path::Path::new(&path);

    let mut g1_bases = Vec::with_capacity(100800000 + 1);
    g1_bases.push(<Bn256 as Engine>::G1Affine::one());
    let mut g2_bases = vec![<Bn256 as Engine>::G2Affine::one()];

    for i in 0..num_chunks {
        let full_path = base_path.join(&format!("transcript{:02}.dat", i));
        println!("Opening {}", full_path.to_string_lossy());
        let file = std::fs::File::open(full_path).map_err(|e| SynthesisError::IoError(e))?;
        let mut reader = std::io::BufReader::with_capacity(1 << 24, file);

        // skip 28 bytes
        let mut tmp = [0u8; 28];
        reader.read_exact(&mut tmp).expect("must skip 28 bytes");

        let mut fq_repr = <Fq as PrimeField>::Repr::default();
        let b_coeff = Fq::from_str("3").unwrap();

        fq_repr.as_mut()[0] = 0x3bf938e377b802a8;
        fq_repr.as_mut()[1] = 0x020b1b273633535d;
        fq_repr.as_mut()[2] = 0x26b7edf049755260;
        fq_repr.as_mut()[3] = 0x2514c6324384a86d;

        let c0 = Fq::from_raw_repr(fq_repr).expect("c0 for B coeff for G2");

        fq_repr.as_mut()[0] = 0x38e7ecccd1dcff67;
        fq_repr.as_mut()[1] = 0x65f0b37d93ce0d3e;
        fq_repr.as_mut()[2] = 0xd749d0dd22ac00aa;
        fq_repr.as_mut()[3] = 0x0141b9ce4a688d4d;

        let c1 = Fq::from_raw_repr(fq_repr).expect("c0 for B coeff for G2");

        let b_coeff_fq2 = Fq2 { c0: c0, c1: c1 };

        for _ in 0..chunk_size {
            // we have to manually read X and Y coordinates
            for k in 0..4 {
                fq_repr.as_mut()[k] = reader.read_u64::<BigEndian>().expect("must read u64");
            }

            let x = Fq::from_repr(fq_repr).expect("must be valid field element encoding");

            for k in 0..4 {
                fq_repr.as_mut()[k] = reader.read_u64::<BigEndian>().expect("must read u64");
            }

            let y = Fq::from_repr(fq_repr).expect("must be valid field element encoding");

            // manual on-curve check
            {
                let mut lhs = y;
                lhs.square();

                let mut rhs = x;
                rhs.square();
                rhs.mul_assign(&x);
                rhs.add_assign(&b_coeff);

                assert!(lhs == rhs);
            }

            let p = <Bn256 as Engine>::G1Affine::from_xy_unchecked(x, y);

            g1_bases.push(p);
        }

        if i == 0 {
            // read G2
            {
                for k in 0..4 {
                    fq_repr.as_mut()[k] = reader.read_u64::<BigEndian>().expect("must read u64");
                }

                let x_c0 = Fq::from_repr(fq_repr).expect("must be valid field element encoding");

                for k in 0..4 {
                    fq_repr.as_mut()[k] = reader.read_u64::<BigEndian>().expect("must read u64");
                }

                let x_c1 = Fq::from_repr(fq_repr).expect("must be valid field element encoding");

                for k in 0..4 {
                    fq_repr.as_mut()[k] = reader.read_u64::<BigEndian>().expect("must read u64");
                }

                let y_c0 = Fq::from_repr(fq_repr).expect("must be valid field element encoding");

                for k in 0..4 {
                    fq_repr.as_mut()[k] = reader.read_u64::<BigEndian>().expect("must read u64");
                }

                let y_c1 = Fq::from_repr(fq_repr).expect("must be valid field element encoding");

                let x = Fq2 { c0: x_c0, c1: x_c1 };

                let y = Fq2 { c0: y_c0, c1: y_c1 };

                {
                    let mut lhs = y;
                    lhs.square();

                    let mut rhs = x;
                    rhs.square();
                    rhs.mul_assign(&x);
                    rhs.add_assign(&b_coeff_fq2);

                    assert!(lhs == rhs);
                }

                let g2 = <Bn256 as Engine>::G2Affine::from_xy_unchecked(x, y);

                g2_bases.push(g2);

                // sanity check by using pairing
                {
                    // check e(g1, g2^x) == e(g1^{x}, g2)
                    let valid = Bn256::final_exponentiation(&Bn256::miller_loop(&[(&g1_bases[0].prepare(), &g2.prepare())])).unwrap()
                        == Bn256::final_exponentiation(&Bn256::miller_loop(&[(&g1_bases[1].prepare(), &g2_bases[0].prepare())])).unwrap();

                    assert!(valid);
                }
            }
            // read G2
            let mut tmp = [0u8; 128];
            reader.read_exact(&mut tmp).expect("must skip 128 bytes of irrelevant G2 point");
        }

        // read to end
        reader.consume(64);

        assert_eq!(reader.fill_buf().unwrap().len(), 0);
    }

    assert_eq!(g1_bases.len(), chunk_size * num_chunks + 1);
    assert_eq!(g2_bases.len(), 2);

    let new = Crs::new(g1_bases, g2_bases);

    Ok(new)
}

fn download_file(url: &str, output_file: &str) -> Result<(), ureq::Error> {
    use std::io::BufWriter;
    let response = ureq::get(url).call()?;
    if response.status() == 200 {
        let mut dest = BufWriter::new(std::fs::File::create(output_file)?);
        let mut reader = response.into_reader();
        std::io::copy(&mut reader, &mut dest)?;
        println!("Transcript file downloaded successfully.");
    } else {
        println!("Failed to download file: HTTP {}", response.status());
    }

    Ok(())
}

pub(crate) fn download_and_transform_ignition_transcripts(domain_size: usize) {
    let transcripts_dir = std::env::var("IGNITION_TRANSCRIPT_PATH").unwrap_or("./".to_string());
    let chunk_size = 5_040_000usize;
    let num_chunks = domain_size.div_ceil(chunk_size);

    let base_url = "https://aztec-ignition.s3.eu-west-2.amazonaws.com/MAIN+IGNITION/sealed";
    for idx in 0..num_chunks {
        let file_url = format!("{}/transcript{:02}.dat", base_url, idx);
        println!("Downloading file at {file_url}");
        let output_file = format!("{}/transcript{:02}.dat", transcripts_dir, idx);
        download_file(&file_url, &output_file).unwrap();
    }

    // transform
    let crs = make_crs_from_ignition_transcript(&transcripts_dir, num_chunks).unwrap();
    let out_path = format!("{}/full_ignition.key", &transcripts_dir);
    let out_file = std::fs::File::create(&out_path).unwrap();

    let Crs { g1_bases, g2_monomial_bases, .. } = crs;
    assert!(g1_bases.len() >= domain_size);
    let mut g1_bases = std::sync::Arc::try_unwrap(g1_bases).unwrap();
    let g2_monomial_bases = std::sync::Arc::try_unwrap(g2_monomial_bases).unwrap();
    g1_bases.truncate(domain_size);

    let crs: Crs<Bn256, CrsForMonomialForm> = Crs::new(g1_bases, g2_monomial_bases);
    crs.write(&out_file).unwrap();
    println!("full ignition ceremony saved into {out_path}");
}
