pragma solidity ^0.8.0;

import "./Plonk4VerifierWithAccessToDNext.sol";
import "./UncheckedMath.sol";

contract Verifier is Plonk4VerifierWithAccessToDNext {
    using UncheckedMath for uint256;

    function get_verification_key() public pure returns(VerificationKey memory vk) {
        vk.num_inputs = 1;
        vk.domain_size = 512;
        vk.omega = PairingsBn254.new_fr(0x0dd30b9ad8c173555d2a33029bc807ac165b61281e9054a173af7ff4e4fc88fc);
        // coefficients
        
        vk.gate_setup_commitments[0] = PairingsBn254.new_g1(0x026c0a977984e62a6f293349667d71480507b34051b957ad0ac785e73839028a,0x1ff4821b9192b1c3eb1f63a57c31231cfdf227c4fdf3e28ef5fec6a3d9105562);
        
        vk.gate_setup_commitments[1] = PairingsBn254.new_g1(0x1fc4f10763c1bdb57e8f26793eea070c1f2a9317735ad097973d5c291f3ebd0b,0x263fbbb9ead5cc24b1430a775a19a1d4c19c0dab7b9e5d44d25f2e49cd1b1cd0);
        
        vk.gate_setup_commitments[2] = PairingsBn254.new_g1(0x0ee3d81f2769d2c6656c420ddf09fcce84b4989f734e9d7ee7b2a9e5534a51dd,0x124cf4a1a38655564efa5997aee04ffba89dc6e6e1c90a0f449d24a920ca85c3);
        
        vk.gate_setup_commitments[3] = PairingsBn254.new_g1(0x2fc30bd1b2c14147665ab9d6d48200d608bb51e5906223e0280a0eb360c03d2a,0x08b01a5413d4cad971b20da06a41828f225dcbda78da225f91671101ab222700);
        
        vk.gate_setup_commitments[4] = PairingsBn254.new_g1(0x1008537e69268b91b6755639e94014b55d67b86b27be1aaba7c2cef49bb8f8cc,0x21f1c559b3ca173e6a7cac2f3c1c380538018da489baac83495302dd86ce9a7a);
        
        vk.gate_setup_commitments[5] = PairingsBn254.new_g1(0x00b3ca41775920fb0efe081cdf105c914de92eec51972ac49f827a98788c11d4,0x2b3f6aa0898924bb4e547c1038794907dd5dbe5a2f7cd6ed957a2722654bb1c2);
        
        vk.gate_setup_commitments[6] = PairingsBn254.new_g1(0x0e0f97ed437fce7d2c42b4f63f3ba78664d2518f377bbfc7399d9dddf33d47af,0x2353c06d2f7c3cf0113c53e64b6c6457360998a6dca6d00867c6cbbc2cf68228);
        
        vk.gate_setup_commitments[7] = PairingsBn254.new_g1(0x0000000000000000000000000000000000000000000000000000000000000000,0x0000000000000000000000000000000000000000000000000000000000000001);
        
        
        // gate selectors
        
        vk.gate_selectors_commitments[0] = PairingsBn254.new_g1(0x103243ba91656521ff89a14b1dd44cffa6aa98bb7aac1005c8d449feb47112b2,0x0eef3e086c6e923975f49579fbfa57c9033ce7136f8dc72c7df2490fef730760);
        
        vk.gate_selectors_commitments[1] = PairingsBn254.new_g1(0x0f741c0f15c95d7090c89895281d1da928d1b4500a1a89e843abaeffc7436efc,0x01147e8779348cfb5deca537d0801c5946d6d2129a7926a664357d43a95c53a6);
        
        
        // permutation
        
        vk.permutation_commitments[0] = PairingsBn254.new_g1(0x2e6f7e063c7973d29265ac3bd204298d3ffaa193491282f2bef6b94cf2007924,0x2dc76bd38ab72c5eae3c45a3a1100fa6facfef1c25e27e273bc5fde8d901f858);
        
        vk.permutation_commitments[1] = PairingsBn254.new_g1(0x2c9266ce160919359ba187b3b1f3ac9ff1b69d874aaa0a8cd38623bc87d3bded,0x29979cad34a229251a85df12da22d23f5b88625f07b59101c53b87be9b5ecfa2);
        
        vk.permutation_commitments[2] = PairingsBn254.new_g1(0x09e326950b3f812891d32e539b9d5f944960e63eb1c1c2f647a065574baa11f1,0x263a306bb6a7dbbb4a99c76b06073b1da64fe4b1cfc9de7c11a853c647a592d8);
        
        vk.permutation_commitments[3] = PairingsBn254.new_g1(0x2067ee34765855a14d66f93012a1301639dc27fe4984fa09478e74d83c0979aa,0x10ef67f3e4744a89933f3e10f1711e3e33b3ccbe452454e6edc00d816f13e0ae);
        
        
        // non residues
        
        vk.non_residues[0] = PairingsBn254.new_fr(0x0000000000000000000000000000000000000000000000000000000000000005);
        
        vk.non_residues[1] = PairingsBn254.new_fr(0x0000000000000000000000000000000000000000000000000000000000000007);
        
        vk.non_residues[2] = PairingsBn254.new_fr(0x000000000000000000000000000000000000000000000000000000000000000a);
        
        
        // g2 elements
        
        vk.g2_elements[0] = PairingsBn254.new_g2([0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2,0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed],[0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b,0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa]);
        
        vk.g2_elements[1] = PairingsBn254.new_g2([0x12740934ba9615b77b6a49b06fcce83ce90d67b1d0e2a530069e3a7306569a91,0x116da8c89a0d090f3d8644ada33a5f1c8013ba7204aeca62d66d931b99afe6e7],[0x25222d9816e5f86b4a7dedd00d04acc5c979c18bd22b834ea8c6d07c0ba441db,0x076441042e77b6309644b56251f059cf14befc72ac8a6157d30924e58dc4c172]);
        
    }

    function deserialize_proof(
        uint256[] calldata public_inputs, 
        uint256[] calldata serialized_proof
    ) internal pure returns(Proof memory proof) {
        require(serialized_proof.length == 34);
        proof.input_values = new uint256[](public_inputs.length);
        for (uint256 i = 0; i < public_inputs.length; i = i.uncheckedInc()) {
            proof.input_values[i] = public_inputs[i];
        }
 
        uint256 j;
        for (uint256 i = 0; i < STATE_WIDTH; i = i.uncheckedInc()) {
            proof.state_polys_commitments[i] = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j.uncheckedInc()]
            );

            j = j.uncheckedAdd(2);
        }
        proof.copy_permutation_grand_product_commitment = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j.uncheckedInc()]
        );
        j = j.uncheckedAdd(2);
        
        
        for (uint256 i = 0; i < proof.quotient_poly_parts_commitments.length; i = i.uncheckedInc()) {
            proof.quotient_poly_parts_commitments[i] = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j.uncheckedInc()]
            );
            j = j.uncheckedAdd(2);
        }

        for (uint256 i = 0; i < proof.state_polys_openings_at_z.length; i = i.uncheckedInc()) {
            proof.state_polys_openings_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j = j.uncheckedInc();
        }

        for (uint256 i = 0; i < proof.state_polys_openings_at_z_omega.length; i = i.uncheckedInc()) {
            proof.state_polys_openings_at_z_omega[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j = j.uncheckedInc();
        } 
        
        for (uint256 i = 0; i < proof.gate_selectors_openings_at_z.length; i = i.uncheckedInc()) {
            proof.gate_selectors_openings_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j = j.uncheckedInc();
        }        
        
        for (uint256 i = 0; i < proof.copy_permutation_polys_openings_at_z.length; i = i.uncheckedInc()) {
            proof.copy_permutation_polys_openings_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j = j.uncheckedInc();
        }
        proof.copy_permutation_grand_product_opening_at_z_omega = PairingsBn254.new_fr(
                serialized_proof[j]
            );

        j = j.uncheckedInc();
        
        proof.quotient_poly_opening_at_z = PairingsBn254.new_fr(
            serialized_proof[j]
        );
        j = j.uncheckedInc();
        proof.linearization_poly_opening_at_z = PairingsBn254.new_fr(
            serialized_proof[j]
        );
        j = j.uncheckedInc();
        proof.opening_proof_at_z = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j.uncheckedInc()]
        );
        j = j.uncheckedAdd(2);
        proof.opening_proof_at_z_omega = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j.uncheckedInc()]
        );
    }
    
    function verify_serialized_proof(
        uint256[] calldata public_inputs, 
        uint256[] calldata serialized_proof
    ) public view returns (bool) {
        VerificationKey memory vk = get_verification_key();
        require(vk.num_inputs == public_inputs.length);

        Proof memory proof = deserialize_proof(public_inputs, serialized_proof);

        return verify(proof, vk);
    }
}
