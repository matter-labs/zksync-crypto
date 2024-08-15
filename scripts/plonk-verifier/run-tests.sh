#!/bin/bash
# set -e

# install foundry

if [[ ! -d $HOME/.foundry/bin ]]; then
    curl -L https://foundry.paradigm.xyz | bash
    . $HOME/.foundry/bin/foundryup 
fi

# run reference test first
PLONK_VERIFIER_DATA_DIR=$PWD/data/plonk-verifier
mkdir -p $PLONK_VERIFIER_DATA_DIR
cargo run --bin generate  --  --verification-key $PLONK_VERIFIER_DATA_DIR/reference_block_20_keccak.key --proof $PLONK_VERIFIER_DATA_DIR/reference_block_20_keccak.proof
cd -
cd crates/plonk-verifier-foundry
$HOME/.foundry/bin/forge test
cd -

# then check sample proofs and vks
mkdir -p $PLONK_VERIFIER_DATA_DIR/std
rm -rf $PLONK_VERIFIER_DATA_DIR/std/*
mkdir -p $PLONK_VERIFIER_DATA_DIR/optimized
rm -rf $PLONK_VERIFIER_DATA_DIR/optimized/*

cd crates/plonk-verifier-codegen
PLONK_VERIFIER_DATA_DIR=$PLONK_VERIFIER_DATA_DIR cargo test test_create_proof_for_all_circuits  --release -- --nocapture
cd -

for main_gate in "std" "optimized"
do    
    for vk_file in $(ls $PLONK_VERIFIER_DATA_DIR/$main_gate/*_vk.json)
    do
        proof_file=$(echo $vk_file | sed s/_vk/_proof/g)

        cargo run --bin generate  -- --encoding json --verification-key $vk_file --proof $proof_file
        cd crates/plonk-verifier-foundry
        $HOME/.foundry/bin/forge test
        cd -
    done
done
