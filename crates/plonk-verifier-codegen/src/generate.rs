use rescue_poseidon::franklin_crypto::bellman::bn256::G1Affine;
use rescue_poseidon::franklin_crypto::bellman::plonk::better_better_cs::gates::selector_optimized_with_d_next::SelectorOptimizedWidth4MainGateWithDNext;
use rescue_poseidon::franklin_crypto::bellman::plonk::better_better_cs::proof::{self, Proof};

use handlebars::*;
use rescue_poseidon::franklin_crypto::bellman::plonk::commitments::transcript::keccak_transcript::RollingKeccakTranscript;
use rescue_poseidon::franklin_crypto::bellman::PrimeField;
use rescue_poseidon::franklin_crypto::bellman::{
    pairing::bn256::{Bn256, Fr},
    plonk::{
        better_better_cs::{cs::Width4MainGateWithDNext, setup::VerificationKey},
        domains::Domain,
    },
    CurveAffine, Engine,
};
use serde_json::Map;
use std::io::Write;
use std::path::PathBuf;

use serde::ser::Serialize;

use crate::circuits::{
    MockCircuitSelectorOptimized, MockCircuitSelectorOptimizedWithLookup, MockCircuitSelectorOptimizedWithLookupAndRescue, MockCircuitSelectorOptimizedWithRescue, MockCircuitWithLookup,
    MockCircuitWithLookupAndRescue, MockCircuitWithRescue,
};
use crate::{
    circuits::MockCircuit,
    serialize::{FieldElement, G1Point, G2Point},
};

pub enum MainGateType {
    Standard,
    SelectorOptimized,
}

struct TemplateVars {
    has_rescue_custom_gate: bool,
    has_lookup: bool,
    is_selector_optimized_main_gate: bool,
    num_main_gate_selectors: usize,
    ab_coeff_idx: usize,
    ac_coeff_idx: usize,
    constant_coeff_idx: usize,
    d_next_coeff_idx: usize,
}

pub enum Encoding {
    Json,
    Default,
}

// Generally it is okey to use MockCircuit in type definitions
// and then transmute it to the corresponding type of circuit
// by inspecting some values
pub fn generate(vk_path: PathBuf, proof_path: Option<PathBuf>, output_dir: PathBuf, encoding_type: Encoding, mut template_files_path: Vec<&str>) {
    let mut reader = std::fs::File::open(vk_path).expect("vk file");

    let vk = match encoding_type {
        Encoding::Json => serde_json::from_reader(reader).expect("read vk from json encoded data"),
        Encoding::Default => VerificationKey::<Bn256, MockCircuit>::read(&mut reader).expect("read vk from default encoded data"),
    };

    // we know from the fact that vk belongs to a
    // - standart main gate when there are 7 selectors
    // - selector optimized main gate when there are 8 selectors
    let num_selectors_of_main_gate = vk.gate_setup_commitments.len();
    let main_gate = if num_selectors_of_main_gate == 7 {
        MainGateType::Standard
    } else if num_selectors_of_main_gate == 8 {
        MainGateType::SelectorOptimized
    } else {
        unimplemented!()
    };

    let has_rescue_custom_gate = if vk.gate_selectors_commitments.len() > 1 {
        assert_eq!(vk.gate_selectors_commitments.len(), 2, "only sbox custom gate is supported");
        true
    } else {
        assert!(vk.gate_selectors_commitments.is_empty());
        false
    };

    let has_lookup = if vk.total_lookup_entries_length > 0 {
        assert!(vk.lookup_selector_commitment.is_some());
        assert!(vk.lookup_tables_commitments.len() > 0);
        assert!(vk.lookup_table_type_commitment.is_some());
        true
    } else {
        assert!(vk.lookup_selector_commitment.is_none());
        assert!(vk.lookup_tables_commitments.len() == 0);
        assert!(vk.lookup_table_type_commitment.is_none());
        false
    };

    let (num_main_gate_selectors, ab_coeff_idx, constant_coeff_idx, d_next_coeff_idx, ac_coeff_idx) = match main_gate {
        MainGateType::Standard => (
            7,
            Width4MainGateWithDNext::AB_MULTIPLICATION_TERM_COEFF_INDEX,
            Width4MainGateWithDNext::CONSTANT_TERM_COEFF_INDEX,
            Width4MainGateWithDNext::D_NEXT_TERM_COEFF_INDEX,
            None,
        ),
        MainGateType::SelectorOptimized => (
            8,
            SelectorOptimizedWidth4MainGateWithDNext::AB_MULTIPLICATION_TERM_COEFF_INDEX,
            SelectorOptimizedWidth4MainGateWithDNext::CONSTANT_TERM_COEFF_INDEX,
            SelectorOptimizedWidth4MainGateWithDNext::D_NEXT_TERM_COEFF_INDEX,
            Some(SelectorOptimizedWidth4MainGateWithDNext::AC_MULTIPLICATION_TERM_COEFF_INDEX),
        ),
    };

    let is_selector_optimized_main_gate = ac_coeff_idx.is_some();
    let ac_coeff_idx = if let Some(coeff) = ac_coeff_idx { coeff } else { 0 };

    let vars = TemplateVars {
        has_rescue_custom_gate,
        has_lookup,
        is_selector_optimized_main_gate,
        num_main_gate_selectors,
        ab_coeff_idx,
        ac_coeff_idx,
        constant_coeff_idx,
        d_next_coeff_idx,
    };

    let proof_template_dir = template_files_path.pop().unwrap();

    render_verifier(vars, &vk, &output_dir, &template_files_path);

    if let Some(proof_path) = proof_path {
        let mut reader = std::fs::File::open(proof_path).expect("proof file");
        let proof = match encoding_type {
            Encoding::Json => serde_json::from_reader(reader).expect("read proof from json encoded data"),
            Encoding::Default => Proof::<Bn256, MockCircuit>::read(&mut reader).expect("read proof from default encoded data"),
        };
        unsafe { verify_proof(&vk, &proof, main_gate) };
        render_expected_proofs(proof, proof_template_dir, &output_dir, has_rescue_custom_gate, has_lookup);
    }
}

unsafe fn verify_proof(vk: &VerificationKey<Bn256, MockCircuit>, proof: &Proof<Bn256, MockCircuit>, main_gate_type: MainGateType) {
    use rescue_poseidon::franklin_crypto::bellman::plonk::better_better_cs::verifier::verify;

    assert_eq!(vk.n, proof.n);
    assert_eq!(vk.num_inputs, proof.inputs.len());

    let has_lookup = vk.total_lookup_entries_length > 0;
    assert_eq!(has_lookup, proof.lookup_grand_product_commitment.is_some());
    assert_eq!(has_lookup, proof.lookup_s_poly_commitment.is_some());
    assert_eq!(has_lookup, proof.lookup_grand_product_opening_at_z_omega.is_some());
    assert_eq!(has_lookup, proof.lookup_s_poly_opening_at_z_omega.is_some());
    assert_eq!(has_lookup, proof.lookup_selector_poly_opening_at_z.is_some());
    assert_eq!(has_lookup, proof.lookup_t_poly_opening_at_z.is_some());
    assert_eq!(has_lookup, proof.lookup_t_poly_opening_at_z_omega.is_some());
    assert_eq!(has_lookup, proof.lookup_table_type_poly_opening_at_z.is_some());

    let has_custom_gate = vk.gate_selectors_commitments.len() > 0;

    let is_valid = match main_gate_type {
        MainGateType::Standard => {
            if !has_lookup && !has_custom_gate {
                verify::<Bn256, MockCircuit, RollingKeccakTranscript<Fr>>(std::mem::transmute(vk), std::mem::transmute(proof), None).unwrap()
            } else if has_lookup && !has_custom_gate {
                verify::<Bn256, MockCircuitWithLookup, RollingKeccakTranscript<Fr>>(std::mem::transmute(vk), std::mem::transmute(proof), None).unwrap()
            } else if has_custom_gate && !has_lookup {
                verify::<Bn256, MockCircuitWithRescue, RollingKeccakTranscript<Fr>>(std::mem::transmute(vk), std::mem::transmute(proof), None).unwrap()
            } else {
                assert!(has_lookup);
                assert!(has_custom_gate);
                verify::<Bn256, MockCircuitWithLookupAndRescue, RollingKeccakTranscript<Fr>>(std::mem::transmute(vk), std::mem::transmute(proof), None).unwrap()
            }
        }
        MainGateType::SelectorOptimized => {
            if !has_lookup && !has_custom_gate {
                verify::<Bn256, MockCircuitSelectorOptimized, RollingKeccakTranscript<Fr>>(std::mem::transmute(vk), std::mem::transmute(proof), None).unwrap()
            } else if has_lookup && !has_custom_gate {
                verify::<Bn256, MockCircuitSelectorOptimizedWithLookup, RollingKeccakTranscript<Fr>>(std::mem::transmute(vk), std::mem::transmute(proof), None).unwrap()
            } else if has_custom_gate && !has_lookup {
                verify::<Bn256, MockCircuitSelectorOptimizedWithRescue, RollingKeccakTranscript<Fr>>(std::mem::transmute(vk), std::mem::transmute(proof), None).unwrap()
            } else {
                assert!(has_lookup);
                assert!(has_custom_gate);
                verify::<Bn256, MockCircuitSelectorOptimizedWithLookupAndRescue, RollingKeccakTranscript<Fr>>(std::mem::transmute(vk), std::mem::transmute(proof), None).unwrap()
            }
        }
    };

    assert!(is_valid, "proof verification failed at codegen");
}

fn length_of_serialized_proof_from_vk(vk: &VerificationKey<Bn256, MockCircuit>, has_custom_gate: bool, has_lookup: bool, has_dnext: bool, quotient_degree: usize) -> usize {
    assert_eq!(vk.state_width, quotient_degree);

    let mut num_commitments = vk.state_width; // trace
    num_commitments += 1; // copy-perm z(x)
    num_commitments += quotient_degree;
    num_commitments += 2; // opening proofs

    if has_lookup {
        num_commitments += 2; // lookup z(x) + lookup s(x)
    }

    // openings
    let mut num_openings = vk.state_width; // trace
    if has_dnext {
        num_openings += 1;
    }

    if has_custom_gate {
        num_openings += 1; // main gate selector
    }

    num_openings += vk.state_width - 1; // sigmas
                                        // copy-perm z(z) is part of linearizaton
    num_openings += 1; // copy-perm z(z*w)

    if has_lookup {
        // - s(z*w)
        // - z(z*w)
        // - t(z)
        // - t(z*w)
        // - selector(z)
        // - type(z)
        num_openings += 6;
    }

    num_openings += 1; // quotient
    num_openings += 1; // linearization

    2 * num_commitments + num_openings
}

fn compute_quotient_degree(state_width: usize, has_custom_gate: bool, has_lookup: bool) -> usize {
    let mut main_gate_quotient_degree = 2;
    if has_custom_gate {
        main_gate_quotient_degree += 1;
    }

    let copy_perm_quotient_degree = state_width;

    let lookup_quotient_degree = if has_lookup { 2 } else { 0 };

    [main_gate_quotient_degree, copy_perm_quotient_degree, lookup_quotient_degree].iter().cloned().max().unwrap()
}

fn render_verifier(vars: TemplateVars, vk: &VerificationKey<Bn256, MockCircuit>, output_dir: &PathBuf, template_files_path: &[&str]) {
    let mut map = MapWrapper::new();
    let mut handlebars = Handlebars::new();

    map.insert("is_selector_optimized_main_gate", vars.is_selector_optimized_main_gate);
    // main gate + custom rescue
    map.insert("has_lookup", vars.has_lookup);
    map.insert("has_rescue_custom_gate", vars.has_rescue_custom_gate);
    map.insert("MAIN_GATE_AB_COEFF_IDX", vars.ab_coeff_idx);
    map.insert("CONSTANT_TERM_COEFF_INDEX", vars.constant_coeff_idx);
    map.insert("D_NEXT_TERM_COEFF_INDEX", vars.d_next_coeff_idx);
    assert_eq!(vars.ab_coeff_idx, 4);
    map.insert("MAIN_GATE_AC_COEFF_IDX", vars.ac_coeff_idx);
    assert_eq!(vk.gate_setup_commitments.len(), vars.num_main_gate_selectors);
    map.insert("NUM_MAIN_GATE_SELECTORS", vars.num_main_gate_selectors);
    // a, b, c, d
    println!("VK STATE WIDTH {}", vk.state_width);
    map.insert("STATE_WIDTH", vk.state_width);
    map.insert("DNEXT_INDEX", vk.state_width - 1);
    map.insert("NUM_G2_ELS", vk.g2_elements.len());
    map.insert("NUM_LOOKUP_TABLES", vk.lookup_tables_commitments.len());
    let quotient_degree = compute_quotient_degree(vk.state_width, vars.has_rescue_custom_gate, vars.has_lookup);
    let serialized_proof_length = length_of_serialized_proof_from_vk(vk, vars.has_rescue_custom_gate, vars.has_lookup, vars.d_next_coeff_idx > 0, quotient_degree);
    map.insert("SERIALIZED_PROOF_LENGTH", serialized_proof_length);
    let mut num_commitments_at_z = 2 + 4 + 3;
    let mut num_commitments_at_z_omega = 1 + 2;

    let mut num_alpha_challenges = 1 + 2;
    if vars.has_rescue_custom_gate {
        num_commitments_at_z += 1;
        num_alpha_challenges += 3;
    }
    if vars.has_lookup {
        num_commitments_at_z += 3;
        num_alpha_challenges += 3;
        num_commitments_at_z_omega += 3;
    }

    map.insert("rescue_alpha_idx", 1);
    map.insert("num_commitments_at_z", num_commitments_at_z);
    map.insert("num_commitments_at_z_omega", num_commitments_at_z_omega);
    map.insert("NUM_ALPHA_CHALLENGES", num_alpha_challenges);
    if vars.has_rescue_custom_gate {
        map.insert("copy_permutation_alpha_idx", 4);
        map.insert("lookup_alpha_idx", 6);
    } else {
        map.insert("copy_permutation_alpha_idx", 1);
        map.insert("lookup_alpha_idx", 3);
    }

    // domain
    map.insert("num_inputs".into(), vk.num_inputs);
    // assert!(vk.num_inputs > 0);
    let domain: Domain<Fr> = Domain::new_for_size(vk.n as u64).expect("a domain");
    map.insert("domain_size".into(), domain.size);
    map.insert("domain_generator".into(), FieldElement::from(domain.generator));

    // G1Points
    let mut gate_setup_commitments = vec![];
    for cmt in vk.gate_setup_commitments.iter() {
        gate_setup_commitments.push(G1Point::from_affine_point(cmt.clone()))
    }
    map.insert("gate_setup_commitments", gate_setup_commitments);

    let mut gate_selectors_commitments = vec![];
    for cmt in vk.gate_selectors_commitments.iter() {
        gate_selectors_commitments.push(G1Point::from_affine_point(cmt.clone()))
    }
    map.insert("gate_selectors_commitments", gate_selectors_commitments);

    let mut permutation_commitments = vec![];
    for cmt in vk.permutation_commitments.iter() {
        permutation_commitments.push(G1Point::from_affine_point(cmt.clone()))
    }
    map.insert("permutation_commitments", permutation_commitments);

    if vk.total_lookup_entries_length > 0 {
        assert!(vk.lookup_selector_commitment.is_some());
        assert!(vk.lookup_tables_commitments.len() > 0);
        assert!(vk.lookup_table_type_commitment.is_some());

        map.insert("has_lookup", true);
        map.insert(
            "lookup_selector_commitment",
            G1Point::from_affine_point(vk.lookup_selector_commitment.unwrap_or(<Bn256 as Engine>::G1Affine::zero())),
        );
        if vk.total_lookup_entries_length > 0 {
            assert!(vk.lookup_selector_commitment.is_some());
        }
        let mut lookup_tables_commitments = vec![];
        for cmt in vk.lookup_tables_commitments.iter() {
            lookup_tables_commitments.push(G1Point::from_affine_point(cmt.clone()))
        }
        map.insert("lookup_tables_commitments", lookup_tables_commitments);
        map.insert(
            "lookup_table_type_commitment",
            G1Point::from_affine_point(vk.lookup_table_type_commitment.unwrap_or(<Bn256 as Engine>::G1Affine::zero())),
        );
    }

    // non residues
    let mut non_residues = vec![];
    for el in vk.non_residues.iter() {
        non_residues.push(FieldElement::from(el.clone()));
    }
    map.insert("non_residues", non_residues);

    // pairing g2 elements
    let mut g2_elements = vec![];
    for point in vk.g2_elements.iter() {
        g2_elements.push(G2Point::from_affine_point(point.clone()));
    }
    map.insert("g2_elements", g2_elements);

    let mut src_dir = output_dir.clone();
    src_dir.push("src");
    for template_file_path in template_files_path {
        let mut output_path = src_dir.clone();
        output_path.push(template_file_path.split('/').last().unwrap());
        let mut writer = std::fs::File::create(output_path).expect("output file");
        // register template from a file and assign a name to it
        handlebars
            .register_template_file("contract", template_file_path)
            .expect(&format!("must read the template at path {}", template_file_path));

        let rendered = handlebars.render("contract", &map.inner).unwrap();

        writer.write(rendered.as_bytes()).expect("must write to file");
    }
}

fn render_expected_proofs(proof: Proof<Bn256, MockCircuit>, proof_template_dir: &str, output_dir: &PathBuf, has_custom_gate: bool, has_lookup: bool) {
    let output_file = format!("{}/test/HardcodedValues.sol", output_dir.to_string_lossy());
    let mut writer = std::fs::File::create(output_file).expect("output file");
    let mut handlebars = Handlebars::new();
    handlebars
        .register_template_file("contract", proof_template_dir)
        .expect(&format!("must read the template at path {}", proof_template_dir));
    let json_proof = transform_proof_into_json(proof);
    let serialized_inputs = json_proof.inputs.clone();
    let serialized_proof = serialize_proof(json_proof.clone(), has_custom_gate, has_lookup);

    let num_inputs = serialized_inputs.len();
    let serialized_proof_length = serialized_proof.len();
    let hardcoded_values = HardcodedProofValues {
        proof: json_proof,
        serialized_inputs,
        serialized_proof,
        num_inputs,
        serialized_proof_length,
        has_custom_gate,
        has_lookup,
    };
    let rendered = handlebars.render("contract", &hardcoded_values).unwrap();

    writer.write(rendered.as_bytes()).expect("must write to file");
}

#[derive(Clone, Default, serde::Serialize)]
pub struct HardcodedProofValues {
    proof: JsonProof,
    serialized_inputs: Vec<String>,
    serialized_proof: Vec<String>,
    num_inputs: usize,
    serialized_proof_length: usize,
    has_lookup: bool,
    has_custom_gate: bool,
}

#[derive(Clone, Default, serde::Serialize)]
pub struct JsonProof {
    pub n: usize,
    pub inputs: Vec<String>,
    pub state_polys_commitments: Vec<[String; 2]>,
    pub copy_permutation_grand_product_commitment: [String; 2],

    pub lookup_s_poly_commitment: [String; 2],
    pub lookup_grand_product_commitment: [String; 2],

    pub quotient_poly_parts_commitments: Vec<[String; 2]>,

    pub state_polys_openings_at_z: Vec<String>,
    pub state_polys_openings_at_z_omega: Vec<String>,

    pub gate_setup_openings_at_z: Vec<String>,
    pub gate_selectors_openings_at_z: Vec<String>,

    pub copy_permutation_polys_openings_at_z: Vec<String>,
    pub copy_permutation_grand_product_opening_at_z_omega: String,

    pub lookup_s_poly_opening_at_z_omega: String,
    pub lookup_grand_product_opening_at_z_omega: String,

    pub lookup_t_poly_opening_at_z: String,
    pub lookup_t_poly_opening_at_z_omega: String,

    pub lookup_selector_poly_opening_at_z: String,
    pub lookup_table_type_poly_opening_at_z: String,

    pub quotient_poly_opening_at_z: String,

    pub linearization_poly_opening_at_z: String,

    pub opening_proof_at_z: [String; 2],
    pub opening_proof_at_z_omega: [String; 2],
}

pub fn to_hex<F: PrimeField>(el: &F) -> String {
    format!("0x{}", rescue_poseidon::franklin_crypto::bellman::to_hex(el))
}

pub fn point_into_xy(point: &G1Affine) -> [String; 2] {
    let (x, y) = point.as_xy();

    [to_hex(x), to_hex(y)]
}

pub fn transform_proof_into_json(proof: Proof<Bn256, MockCircuit>) -> JsonProof {
    let mut json_proof = JsonProof::default();
    json_proof.n = proof.n;
    json_proof.inputs = proof.inputs.iter().map(|v| to_hex(v)).collect();
    json_proof.state_polys_commitments = proof.state_polys_commitments.iter().map(|p| point_into_xy(p)).collect();

    json_proof.copy_permutation_grand_product_commitment = point_into_xy(&proof.copy_permutation_grand_product_commitment);

    json_proof.lookup_s_poly_commitment = proof.lookup_s_poly_commitment.as_ref().map(|p| point_into_xy(p)).unwrap_or_default();
    json_proof.lookup_grand_product_commitment = proof.lookup_grand_product_commitment.as_ref().map(|p| point_into_xy(p)).unwrap_or_default();

    json_proof.quotient_poly_parts_commitments = proof.quotient_poly_parts_commitments.iter().map(|p| point_into_xy(p)).collect();

    json_proof.opening_proof_at_z = point_into_xy(&proof.opening_proof_at_z);
    json_proof.opening_proof_at_z_omega = point_into_xy(&proof.opening_proof_at_z_omega);

    json_proof.state_polys_openings_at_z = proof.state_polys_openings_at_z.iter().map(|v| to_hex(v)).collect();
    assert_eq!(proof.state_polys_openings_at_dilations.len(), 1, "only one dilation is allowed");
    json_proof.state_polys_openings_at_z_omega = proof.state_polys_openings_at_dilations.iter().map(|(_, _, v)| to_hex(v)).collect();
    json_proof.gate_setup_openings_at_z = proof.gate_setup_openings_at_z.iter().map(|(_, _, v)| to_hex(v)).collect();
    json_proof.gate_selectors_openings_at_z = proof.gate_selectors_openings_at_z.iter().map(|(_, v)| to_hex(v)).collect();
    json_proof.copy_permutation_polys_openings_at_z = proof.copy_permutation_polys_openings_at_z.iter().map(|v| to_hex(v)).collect();
    json_proof.copy_permutation_grand_product_opening_at_z_omega = to_hex(&proof.copy_permutation_grand_product_opening_at_z_omega);

    json_proof.lookup_s_poly_opening_at_z_omega = proof.lookup_s_poly_opening_at_z_omega.as_ref().map(|v| to_hex(v)).unwrap_or_default();

    json_proof.lookup_grand_product_opening_at_z_omega = proof.lookup_grand_product_opening_at_z_omega.as_ref().map(|v| to_hex(v)).unwrap_or_default();

    json_proof.lookup_t_poly_opening_at_z = proof.lookup_t_poly_opening_at_z.as_ref().map(|v| to_hex(v)).unwrap_or_default();
    json_proof.lookup_t_poly_opening_at_z_omega = proof.lookup_t_poly_opening_at_z_omega.as_ref().map(|v| to_hex(v)).unwrap_or_default();
    json_proof.lookup_selector_poly_opening_at_z = proof.lookup_selector_poly_opening_at_z.as_ref().map(|v| to_hex(v)).unwrap_or_default();
    json_proof.lookup_table_type_poly_opening_at_z = proof.lookup_table_type_poly_opening_at_z.as_ref().map(|v| to_hex(v)).unwrap_or_default();

    json_proof.quotient_poly_opening_at_z = to_hex(&proof.quotient_poly_opening_at_z);
    json_proof.linearization_poly_opening_at_z = to_hex(&proof.linearization_poly_opening_at_z);

    json_proof
}

pub fn serialize_proof(proof: JsonProof, has_custom_gate: bool, has_lookup: bool) -> Vec<String> {
    let JsonProof {
        state_polys_commitments,
        copy_permutation_grand_product_commitment,
        lookup_s_poly_commitment,
        lookup_grand_product_commitment,
        quotient_poly_parts_commitments,
        state_polys_openings_at_z,
        state_polys_openings_at_z_omega,
        gate_setup_openings_at_z,
        gate_selectors_openings_at_z,
        copy_permutation_polys_openings_at_z,
        copy_permutation_grand_product_opening_at_z_omega,
        lookup_s_poly_opening_at_z_omega,
        lookup_grand_product_opening_at_z_omega,
        lookup_t_poly_opening_at_z,
        lookup_t_poly_opening_at_z_omega,
        lookup_selector_poly_opening_at_z,
        lookup_table_type_poly_opening_at_z,
        quotient_poly_opening_at_z,
        linearization_poly_opening_at_z,
        opening_proof_at_z,
        opening_proof_at_z_omega,
        ..
    } = proof;

    let mut serialized_proof = vec![];
    serialized_proof.extend(state_polys_commitments.iter().flat_map(|inner| inner.iter().cloned()));

    serialized_proof.extend(copy_permutation_grand_product_commitment);
    if has_lookup {
        serialized_proof.extend(lookup_s_poly_commitment);
        serialized_proof.extend(lookup_grand_product_commitment);
    }

    serialized_proof.extend(quotient_poly_parts_commitments.iter().flat_map(|inner| inner.iter().cloned()));

    serialized_proof.extend(state_polys_openings_at_z);
    serialized_proof.extend(state_polys_openings_at_z_omega);
    serialized_proof.extend(gate_setup_openings_at_z);
    if has_custom_gate {
        // linearization includes gate selector of the sbox gate
        assert_eq!(gate_selectors_openings_at_z.len(), 1);
    } else {
        assert!(gate_selectors_openings_at_z.is_empty());
    }
    serialized_proof.extend(gate_selectors_openings_at_z);
    serialized_proof.extend(copy_permutation_polys_openings_at_z);
    serialized_proof.push(copy_permutation_grand_product_opening_at_z_omega);
    if has_lookup {
        serialized_proof.push(lookup_s_poly_opening_at_z_omega);
        serialized_proof.push(lookup_grand_product_opening_at_z_omega);
        serialized_proof.push(lookup_t_poly_opening_at_z);
        serialized_proof.push(lookup_t_poly_opening_at_z_omega);
        serialized_proof.push(lookup_selector_poly_opening_at_z);
        serialized_proof.push(lookup_table_type_poly_opening_at_z);
    }

    serialized_proof.push(quotient_poly_opening_at_z);
    serialized_proof.push(linearization_poly_opening_at_z);

    serialized_proof.extend(opening_proof_at_z);
    serialized_proof.extend(opening_proof_at_z_omega);

    serialized_proof
}
struct MapWrapper {
    inner: Map<String, JsonValue>,
}
impl MapWrapper {
    fn new() -> Self {
        Self { inner: Map::new() }
    }

    fn insert<T: Serialize>(&mut self, key: &str, value: T) -> Option<JsonValue> {
        self.inner.insert(key.into(), to_json(value))
    }
}
