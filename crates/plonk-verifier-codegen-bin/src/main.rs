use plonk_verifier_codegen::{generate, Encoding};
use std::path::PathBuf;
use structopt::StructOpt;

const DEFAULT_OUTPUT_FILE: &str = "./crates/plonk-verifier-foundry";

const PAIRING_BN_254_FILE_PATH: &str = "./crates/plonk-verifier-codegen/template/PairingsBn254.sol";
const TRANSCRIPT_LIB_FILE_PATH: &str = "./crates/plonk-verifier-codegen/template/TranscriptLib.sol";
const UNCHECKED_MATH_FILE_PATH: &str = "./crates/plonk-verifier-codegen/template/UncheckedMath.sol";
const PLONK_4_VERIFIER_FILE_PATH: &str = "./crates/plonk-verifier-codegen/template/Plonk4VerifierWithAccessToDNext.sol";
const VEERIFIER_TEMPLATE_FILE_PATH: &str = "./crates/plonk-verifier-codegen/template/Verifier.sol";
const PROOF_TEST_TEMPLATE_FILE_PATH: &str = "./crates/plonk-verifier-codegen/template/HardcodedValues.sol";

#[derive(StructOpt, Debug)]
pub struct Opts {
    /// Path to verification key(required)
    #[structopt(long, parse(from_os_str))]
    verification_key: PathBuf,
    /// Path to proof key(optional)
    #[structopt(long, parse(from_os_str))]
    proof: Option<PathBuf>,
    /// Output directory
    #[structopt(long, parse(from_os_str), default_value = DEFAULT_OUTPUT_FILE)]
    output: PathBuf,

    #[structopt(long)]
    encoding: Option<String>,
}

fn main() {
    let opts = Opts::from_args();
    println!("{:#?}", opts);

    let Opts {
        verification_key,
        proof,
        output,
        encoding,
    } = opts;

    let encoding = match encoding {
        Some(encoding) => match encoding.as_str() {
            "json" => Encoding::Json,
            _ => Encoding::Default,
        },
        None => Encoding::Default,
    };

    generate(
        verification_key,
        proof,
        output.clone(),
        encoding,
        vec![
            VEERIFIER_TEMPLATE_FILE_PATH,
            PLONK_4_VERIFIER_FILE_PATH,
            TRANSCRIPT_LIB_FILE_PATH,
            PAIRING_BN_254_FILE_PATH,
            UNCHECKED_MATH_FILE_PATH,
            PROOF_TEST_TEMPLATE_FILE_PATH,
        ],
    );

    eprintln!("Success!");
}
