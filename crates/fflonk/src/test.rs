use super::*;
use franklin_crypto::bellman::{
    bn256::{Bn256, Fr},
    plonk::{
        better_better_cs::{
            cs::{Circuit, PlonkCsWidth3Params},
            gates::naive_main_gate::NaiveMainGate,
        },
        commitments::transcript::keccak_transcript::RollingKeccakTranscript,
    },
    worker::Worker,
};

#[test]
#[ignore]
fn test_test_circuit_with_naive_main_gate() {
    use crate::definitions::setup::FflonkSetup;
    use crate::utils::compute_max_combined_degree_from_assembly;
    use crate::verifier::FflonkVerificationKey;
    use franklin_crypto::bellman::plonk::better_better_cs::cs::TrivialAssembly;
    let worker = Worker::new();
    let circuit = crate::FflonkTestCircuit {};

    let mut assembly = TrivialAssembly::<Bn256, PlonkCsWidth3Params, NaiveMainGate>::new();
    circuit.synthesize(&mut assembly).expect("must work");
    assert!(assembly.is_satisfied());
    assembly.finalize();
    let domain_size = assembly.n() + 1;
    assert!(domain_size.is_power_of_two());
    assert!(domain_size <= 1 << L1_VERIFIER_DOMAIN_SIZE_LOG);
    println!("Trace log length {}", domain_size.trailing_zeros());

    let max_combined_degree = compute_max_combined_degree_from_assembly::<_, _, _, _, FflonkTestCircuit>(&assembly);
    println!("Max degree is {}", max_combined_degree);
    let mon_crs = init_crs(&worker, domain_size);
    let setup = FflonkSetup::create_setup(&assembly, &worker, &mon_crs).expect("setup");
    let vk = FflonkVerificationKey::from_setup(&setup, &mon_crs).unwrap();
    let vk_file = std::fs::File::create("/tmp/test_vk.json").unwrap();
    serde_json::to_writer(&vk_file, &vk).unwrap();
    println!("vk file saved");

    let proof = crate::prover::create_proof::<_, FflonkTestCircuit, _, _, _, RollingKeccakTranscript<Fr>>(&assembly, &worker, &setup, &mon_crs, None).expect("proof");
    dbg!(&proof.commitments);
    let valid = crate::verify::<_, _, RollingKeccakTranscript<Fr>>(&vk, &proof, None).unwrap();
    assert!(valid, "proof verification fails");
}

pub fn init_crs(worker: &Worker, domain_size: usize) -> Crs<Bn256, CrsForMonomialForm> {
    assert!(domain_size <= 1 << L1_VERIFIER_DOMAIN_SIZE_LOG);
    let num_points = MAX_COMBINED_DEGREE_FACTOR * domain_size;
    let mon_crs = if let Ok(crs_file_path) = std::env::var("CRS_FILE") {
        println!("using crs file at {crs_file_path}");
        let crs_file = std::fs::File::open(&crs_file_path).expect(&format!("crs file at {}", crs_file_path));
        let mon_crs = Crs::<Bn256, CrsForMonomialForm>::read(crs_file).expect(&format!("read crs file at {}", crs_file_path));
        assert!(num_points <= mon_crs.g1_bases.len());

        mon_crs
    } else {
        Crs::<Bn256, CrsForMonomialForm>::non_power_of_two_crs_42(num_points, &worker)
    };

    mon_crs
}
