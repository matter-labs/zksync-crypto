use std::path::PathBuf;
use structopt::StructOpt;
use zksync_solidity_vk_codegen::{generate, Encoding};

const DEFAULT_OUTPUT_FILE: &str = "./hardhat/contracts";

const PAIRING_BN_254_FILE_PATH: &str = "PairingsBn254.sol";
const TRANSCRIPT_LIB_FILE_PATH: &str = "TranscriptLib.sol";
const UNCHECKED_MATH_FILE_PATH: &str = "UncheckedMath.sol";
const PLONK_4_VERIFIER_FILE_PATH: &str = "Plonk4VerifierWithAccessToDNext.sol";
const VEERIFIER_TEMPLATE_FILE_PATH: &str = "Verifier.sol";

#[derive(StructOpt, Debug)]
pub struct Opts {
    /// Path to verification key (required)
    #[structopt(long, parse(from_os_str))]
    verification_key: PathBuf,
    /// Path to the folder with templates (required)
    /// Should point to the `codegen/template` directory
    #[structopt(long, parse(from_os_str))]
    templates_dir: PathBuf,
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
        output,
        templates_dir,
        encoding,
    } = opts;

    let encoding = match encoding {
        Some(encoding) => match encoding.as_str() {
            "json" => Encoding::Json,
            _ => Encoding::Default,
        },
        None => Encoding::Default,
    };

    let template_path = |file_name: &str| templates_dir.join(file_name).to_string_lossy().into_owned();

    generate(
        verification_key,
        output.clone(),
        encoding,
        vec![
            template_path(VEERIFIER_TEMPLATE_FILE_PATH).as_ref(),
            template_path(PLONK_4_VERIFIER_FILE_PATH).as_ref(),
            template_path(TRANSCRIPT_LIB_FILE_PATH).as_ref(),
            template_path(PAIRING_BN_254_FILE_PATH).as_ref(),
            template_path(UNCHECKED_MATH_FILE_PATH).as_ref(),
        ],
    );

    eprintln!("Success!");
}
