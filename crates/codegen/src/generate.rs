use franklin_crypto::bellman::plonk::better_better_cs::gates::selector_optimized_with_d_next::SelectorOptimizedWidth4MainGateWithDNext;
use franklin_crypto::bellman::{
    pairing::bn256::{Bn256, Fr},
    plonk::{
        better_better_cs::{cs::Width4MainGateWithDNext, setup::VerificationKey},
        domains::Domain,
    },
    CurveAffine, Engine,
};
use handlebars::*;
use serde_json::Map;
use std::io::Write;
use std::path::PathBuf;

use serde::ser::Serialize;

use crate::{
    circuits::DummyCircuit,
    serialize::{FieldElement, G1Point, G2Point},
};

pub enum MainGateType {
    Standard,
    SelectorOptimized,
}

struct TemplateVars {
    num_gates: usize,
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

pub fn generate(vk_path: PathBuf, output_dir: PathBuf, encoding_type: Encoding, template_files_path: Vec<&str>) {
    let mut reader = std::fs::File::open(vk_path).expect("vk file");

    let vk = match encoding_type {
        Encoding::Json => serde_json::from_reader(reader).expect("read vk from json encoded data"),
        Encoding::Default => VerificationKey::<Bn256, DummyCircuit>::read(&mut reader).expect("read vk from default encoded data"),
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

    let num_gates = if vk.gate_selectors_commitments.len() == 0 { 1 } else { vk.gate_selectors_commitments.len() };

    let has_rescue_custom_gate = if num_gates > 1 { true } else { false };

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
        num_gates,
        has_rescue_custom_gate,
        has_lookup,
        is_selector_optimized_main_gate,
        num_main_gate_selectors,
        ab_coeff_idx,
        ac_coeff_idx,
        constant_coeff_idx,
        d_next_coeff_idx,
    };

    render(vars, vk, output_dir, template_files_path)
}

fn render(vars: TemplateVars, vk: VerificationKey<Bn256, DummyCircuit>, output_dir: PathBuf, template_files_path: Vec<&str>) {
    let mut map = MapWrapper::new();
    let mut handlebars = Handlebars::new();

    map.insert("is_selector_optimized_main_gate", vars.is_selector_optimized_main_gate);
    // main gate + custom rescue
    map.insert("NUM_GATES", vars.num_gates);
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
    map.insert("SERIALIZED_PROOF_LENGTH", 44);
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
    assert!(vk.num_inputs > 0);
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

    for template_file_path in template_files_path {
        let mut output_path = output_dir.clone();
        output_path.push(template_file_path.split('/').last().unwrap());
        let mut writer = std::fs::File::create(output_path).expect("output file");
        // register template from a file and assign a name to it
        handlebars
            .register_template_file("contract", template_file_path.clone())
            .expect(&format!("must read the template at path {}", template_file_path));

        let rendered = handlebars.render("contract", &map.inner).unwrap();

        writer.write(rendered.as_bytes()).expect("must write to file");
    }
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
