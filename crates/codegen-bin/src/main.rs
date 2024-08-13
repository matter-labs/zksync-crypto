use std::path::PathBuf;
use structopt::StructOpt;
use zksync_solidity_vk_codegen::{generate, Encoding};

const DEFAULT_OUTPUT_FILE: &str = "./hardhat/contracts";

const PAIRING_BN_254_FILE_PATH: &str = "./codegen/template/PairingsBn254.sol";
const TRANSCRIPT_LIB_FILE_PATH: &str = "./codegen/template/TranscriptLib.sol";
const UNCHECKED_MATH_FILE_PATH: &str = "./codegen/template/UncheckedMath.sol";
const PLONK_4_VERIFIER_FILE_PATH: &str = "./codegen/template/Plonk4VerifierWithAccessToDNext.sol";
const VEERIFIER_TEMPLATE_FILE_PATH: &str = "./codegen/template/Verifier.sol";

#[derive(StructOpt, Debug)]
pub struct Opts {
    /// Path to verification key(required)
    #[structopt(long, parse(from_os_str))]
    verification_key: PathBuf,
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
        output.clone(),
        encoding,
        vec![
            VEERIFIER_TEMPLATE_FILE_PATH,
            PLONK_4_VERIFIER_FILE_PATH,
            TRANSCRIPT_LIB_FILE_PATH,
            PAIRING_BN_254_FILE_PATH,
            UNCHECKED_MATH_FILE_PATH,
        ],
    );

    eprintln!("Success!");
}
