pragma solidity >=0.8.25 <0.9.0;
import {Proof} from "../src/Plonk4VerifierWithAccessToDNext.sol";

library HardcodedValues{
    function hardcoded_proof() internal pure returns(Proof memory proof){
    
        proof.input_values[0] = 0x0000000000000000000000000000000000000000000000000000000000000002;
    
    
        proof.state_polys_commitments[0].X = 0x125a96306030c34958b4ad672033277ea109a796adf28dc4976027a44c12ef22;
        proof.state_polys_commitments[0].Y = 0x00946c4bb159e3d8576b2b0be5bf7b3fc0b61f1c2e72da4374f9f0ad26f293b7;
    
        proof.state_polys_commitments[1].X = 0x25bdc8e25058040ad05096a4fd2cd0a81cbad973deda5f7415840e6593439de9;
        proof.state_polys_commitments[1].Y = 0x13d3fae4b771527d6283610904ee1c95f6d3e546ed861bd071f0c910af39d0e6;
    
        proof.state_polys_commitments[2].X = 0x0e0fb43d701a22c3fdd6670b69eca01b815ef2f6b3687191a4fb819d25c5e117;
        proof.state_polys_commitments[2].Y = 0x22c377922737d710453e0e6efb6d8b5b6b6fc346f4bd811802db23d4f7fc8019;
    
        proof.state_polys_commitments[3].X = 0x2a0edd1ead0847b61d9b9e6e026c47a758afbe347ae119486b9ad1719f6c367d;
        proof.state_polys_commitments[3].Y = 0x21e058312fc5a32072d821eadf55644320e7ca6e9ec236925edd632fc3ff463d;
    
        proof.copy_permutation_grand_product_commitment.X = 0x013bb8312203b14b634a689c6d4acef9d1f54e1068de1470e19279ede8bb991f;
        proof.copy_permutation_grand_product_commitment.Y = 0x0c7b763a6bd20f6fd5295c9e85846c4b132b26bb5d5acf71bcb29ee2a42ed5be;
    
        proof.quotient_poly_parts_commitments[0].X = 0x1837c3c200f0aebb25d7e1a2af2685892fc0c701433996b3b80953505f4efd9c;
        proof.quotient_poly_parts_commitments[0].Y = 0x1a9e03970fbb2ad6876f0f0bae91097eb19ac3face44eaae4396fc16ab9116ad;
    
        proof.quotient_poly_parts_commitments[1].X = 0x060df9de261271340128c4b3c77b2a41178e5b87db2355d3fe1930823abe56ac;
        proof.quotient_poly_parts_commitments[1].Y = 0x022f38e477d4a7a8ad6d91d56443d53628bbfc2b5d48e577f8a28bbd0f9a90a7;
    
        proof.quotient_poly_parts_commitments[2].X = 0x1af9ec6291b5a9da75e943391bcd283a51cc8a7c1a4825c01395e6b840bd16b4;
        proof.quotient_poly_parts_commitments[2].Y = 0x12b652f3ffe11917e28cd4731156df8a8143299d19c9a333a99e41e666c1baee;
    
        proof.quotient_poly_parts_commitments[3].X = 0x0daad57410d2430d709a9d5904c892b45189012336a110e6d89cf670d15a0b31;
        proof.quotient_poly_parts_commitments[3].Y = 0x2c064d6c18d55c0d370f56cbb0548e21f7866f02a52e63ec774805f32bf5c0f0;
    
        
        proof.opening_proof_at_z.X = 0x0d63c1990dcd7716af77977ce95bc0985cd416ca76bb65aebed69c79e85ae67f;
        proof.opening_proof_at_z.Y = 0x11e0cd7c46cf5343986e816559536280ffcc362957c9074aa6b0cbce69cfd1a5;
        proof.opening_proof_at_z_omega.X = 0x0c4d652ff36b77dc918b8df517d13a95c2150fd3163e03b19480bf733a6fdc39;
        proof.opening_proof_at_z_omega.Y = 0x20ded5c2838da911842c0ca7e65095623e3d8e2dd177fc58377119834084022b;

    
        proof.state_polys_openings_at_z[0].value = 0x260fb5a68c251e1d5fb776f1ae703a8c5e98ceb86cd57cfd03314b3a2dd78e44;
    
        proof.state_polys_openings_at_z[1].value = 0x22ca21b7257fb02fd9709e2c71c074e4fd8b2d3fbd969c2dc4624be529fd7cb2;
    
        proof.state_polys_openings_at_z[2].value = 0x1dc0e5fd84d3254929882dae104a0aa024caac4732894dedd1abb620c7d9f70a;
    
        proof.state_polys_openings_at_z[3].value = 0x28ae354f4b84775ea856462e682fda6e3489082b380423daadc97ed741a39eb9;
    
    
        proof.state_polys_openings_at_z_omega[0].value = 0x2f926ea1197a008ab9a586b00063fd9a1c61627211d62932780ef7b3ff9c58dd;
    
    
        proof.gate_selectors_openings_at_z[0].value = 0x0a6bf79bcd59025be43d0a022107876cf6adf94eb907448727c43ac75ec3c94e;
    
    
        proof.copy_permutation_polys_openings_at_z[0].value = 0x2067b179640d9fbca42c4814286b9211e4de99fa74dcbf892d9cc88ae9b13bf9;
    
        proof.copy_permutation_polys_openings_at_z[1].value = 0x2094e6d2b8ae0afa2e548e05d3fda117da3e595a78f2818d54406627e85e69f6;
    
        proof.copy_permutation_polys_openings_at_z[2].value = 0x2c4b1065f0312db8067a3a627e06a07a72dbbcaadea183f7ab12b3cabc6611d8;
    
        proof.copy_permutation_grand_product_opening_at_z_omega.value = 0x06856e549033f929f5e4433bb0d2f7c80b2c27ce5c84cbb6745d315ea2db10a7;
        proof.quotient_poly_opening_at_z.value = 0x201be4dd1d31bd3e0677b63dd095f21ecd194c240f9cbdb776ac3318fa1287cf;
        proof.linearization_poly_opening_at_z.value = 0x02935a5e887039918f79ed4e98be063a48c6af224f8cbfde3e52ea936e524a6f;
    
    }    

    function hardcoded_serialized_proof() public pure returns(uint256[] memory serializedInputs, uint256[] memory serializedProof){
        serializedInputs = new uint256[](1);
    
        serializedInputs[0] = 0x0000000000000000000000000000000000000000000000000000000000000002;
    
        serializedProof = new uint256[](34);
    
        serializedProof[0] = 0x125a96306030c34958b4ad672033277ea109a796adf28dc4976027a44c12ef22;
    
        serializedProof[1] = 0x00946c4bb159e3d8576b2b0be5bf7b3fc0b61f1c2e72da4374f9f0ad26f293b7;
    
        serializedProof[2] = 0x25bdc8e25058040ad05096a4fd2cd0a81cbad973deda5f7415840e6593439de9;
    
        serializedProof[3] = 0x13d3fae4b771527d6283610904ee1c95f6d3e546ed861bd071f0c910af39d0e6;
    
        serializedProof[4] = 0x0e0fb43d701a22c3fdd6670b69eca01b815ef2f6b3687191a4fb819d25c5e117;
    
        serializedProof[5] = 0x22c377922737d710453e0e6efb6d8b5b6b6fc346f4bd811802db23d4f7fc8019;
    
        serializedProof[6] = 0x2a0edd1ead0847b61d9b9e6e026c47a758afbe347ae119486b9ad1719f6c367d;
    
        serializedProof[7] = 0x21e058312fc5a32072d821eadf55644320e7ca6e9ec236925edd632fc3ff463d;
    
        serializedProof[8] = 0x013bb8312203b14b634a689c6d4acef9d1f54e1068de1470e19279ede8bb991f;
    
        serializedProof[9] = 0x0c7b763a6bd20f6fd5295c9e85846c4b132b26bb5d5acf71bcb29ee2a42ed5be;
    
        serializedProof[10] = 0x1837c3c200f0aebb25d7e1a2af2685892fc0c701433996b3b80953505f4efd9c;
    
        serializedProof[11] = 0x1a9e03970fbb2ad6876f0f0bae91097eb19ac3face44eaae4396fc16ab9116ad;
    
        serializedProof[12] = 0x060df9de261271340128c4b3c77b2a41178e5b87db2355d3fe1930823abe56ac;
    
        serializedProof[13] = 0x022f38e477d4a7a8ad6d91d56443d53628bbfc2b5d48e577f8a28bbd0f9a90a7;
    
        serializedProof[14] = 0x1af9ec6291b5a9da75e943391bcd283a51cc8a7c1a4825c01395e6b840bd16b4;
    
        serializedProof[15] = 0x12b652f3ffe11917e28cd4731156df8a8143299d19c9a333a99e41e666c1baee;
    
        serializedProof[16] = 0x0daad57410d2430d709a9d5904c892b45189012336a110e6d89cf670d15a0b31;
    
        serializedProof[17] = 0x2c064d6c18d55c0d370f56cbb0548e21f7866f02a52e63ec774805f32bf5c0f0;
    
        serializedProof[18] = 0x260fb5a68c251e1d5fb776f1ae703a8c5e98ceb86cd57cfd03314b3a2dd78e44;
    
        serializedProof[19] = 0x22ca21b7257fb02fd9709e2c71c074e4fd8b2d3fbd969c2dc4624be529fd7cb2;
    
        serializedProof[20] = 0x1dc0e5fd84d3254929882dae104a0aa024caac4732894dedd1abb620c7d9f70a;
    
        serializedProof[21] = 0x28ae354f4b84775ea856462e682fda6e3489082b380423daadc97ed741a39eb9;
    
        serializedProof[22] = 0x2f926ea1197a008ab9a586b00063fd9a1c61627211d62932780ef7b3ff9c58dd;
    
        serializedProof[23] = 0x0a6bf79bcd59025be43d0a022107876cf6adf94eb907448727c43ac75ec3c94e;
    
        serializedProof[24] = 0x2067b179640d9fbca42c4814286b9211e4de99fa74dcbf892d9cc88ae9b13bf9;
    
        serializedProof[25] = 0x2094e6d2b8ae0afa2e548e05d3fda117da3e595a78f2818d54406627e85e69f6;
    
        serializedProof[26] = 0x2c4b1065f0312db8067a3a627e06a07a72dbbcaadea183f7ab12b3cabc6611d8;
    
        serializedProof[27] = 0x06856e549033f929f5e4433bb0d2f7c80b2c27ce5c84cbb6745d315ea2db10a7;
    
        serializedProof[28] = 0x201be4dd1d31bd3e0677b63dd095f21ecd194c240f9cbdb776ac3318fa1287cf;
    
        serializedProof[29] = 0x02935a5e887039918f79ed4e98be063a48c6af224f8cbfde3e52ea936e524a6f;
    
        serializedProof[30] = 0x0d63c1990dcd7716af77977ce95bc0985cd416ca76bb65aebed69c79e85ae67f;
    
        serializedProof[31] = 0x11e0cd7c46cf5343986e816559536280ffcc362957c9074aa6b0cbce69cfd1a5;
    
        serializedProof[32] = 0x0c4d652ff36b77dc918b8df517d13a95c2150fd3163e03b19480bf733a6fdc39;
    
        serializedProof[33] = 0x20ded5c2838da911842c0ca7e65095623e3d8e2dd177fc58377119834084022b;
    
    }
}
