// SPDX-License-Identifier: UNLICENSED
pragma solidity >=0.8.25 <0.9.0;

import {Test, console} from "forge-std/src/Test.sol";
import {Verifier} from "../src/Verifier.sol";
import {HardcodedValues} from "./HardcodedValues.sol";

contract VerifierTest is Verifier {
    Verifier public _verifier;

    function setUp() public {
        _verifier = new Verifier();
    }

    function test1FullProtocol() public view {
        uint256 _g0 = gasleft();
        (
            uint256[] memory serializedInputs,
            uint256[] memory serializedProof
        ) = HardcodedValues.hardcoded_serialized_proof();
        bool valid = verify_serialized_proof_helper(
            serializedInputs,
            serializedProof
        );
        uint256 _g1 = gasleft();
        require(valid, "proof is valid");
        console.log("test8FullProtocol gas cost: %d", _g0 - _g1);
    }

    function verify_serialized_proof_helper(
        uint256[] memory serializedInputs,
        uint256[] memory serializedProof
    ) public view  returns (bool) {
        console.log("verifying proof");
        return this.verify_serialized_proof(serializedInputs, serializedProof);
    }
}
