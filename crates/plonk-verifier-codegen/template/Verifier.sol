pragma solidity ^0.8.0;

import "./Plonk4VerifierWithAccessToDNext.sol";
import "./UncheckedMath.sol";

contract Verifier is Plonk4VerifierWithAccessToDNext {
    using UncheckedMath for uint256;

    function get_verification_key() public pure returns(VerificationKey memory vk) {
        vk.num_inputs = {{num_inputs}};
        vk.domain_size = {{domain_size}};
        vk.omega = {{domain_generator.el}};
        // coefficients
        {{#each gate_setup_commitments}}
        vk.gate_setup_commitments[{{@index}}] = {{this.g1}};
        {{/each}}
        {{#if has_rescue_custom_gate}}
        // gate selectors
        {{#each gate_selectors_commitments}}
        vk.gate_selectors_commitments[{{@index}}] = {{this.g1}};
        {{/each}}
        {{/if}}
        // permutation
        {{#each permutation_commitments}}
        vk.permutation_commitments[{{@index}}] = {{this.g1}};
        {{/each}}
        {{#if has_lookup}}
        // lookup table commitments
        vk.lookup_selector_commitment = {{lookup_selector_commitment.g1}};
        {{#each lookup_tables_commitments}}
        vk.lookup_tables_commitments[{{@index}}] = {{this.g1}};
        {{/each}}
        vk.lookup_table_type_commitment = {{lookup_table_type_commitment.g1}};
        {{/if}}
        // non residues
        {{#each non_residues}}
        vk.non_residues[{{@index}}] = {{this.el}};
        {{/each}}
        
        // g2 elements
        {{#each g2_elements}}
        vk.g2_elements[{{@index}}] = {{this.g2}};
        {{/each}}
    }

    function deserialize_proof(
        uint256[] calldata public_inputs, 
        uint256[] calldata serialized_proof
    ) internal pure returns(Proof memory proof) {
        require(serialized_proof.length == {{SERIALIZED_PROOF_LENGTH}});
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
        
        {{#if has_lookup}}
        proof.lookup_s_poly_commitment = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j.uncheckedInc()]
        );
        j = j.uncheckedAdd(2);

        proof.lookup_grand_product_commitment = PairingsBn254.new_g1_checked(
                serialized_proof[j],
                serialized_proof[j.uncheckedInc()]
        );
        j = j.uncheckedAdd(2);
        {{/if}}
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
        {{#if has_rescue_custom_gate}}
        for (uint256 i = 0; i < proof.gate_selectors_openings_at_z.length; i = i.uncheckedInc()) {
            proof.gate_selectors_openings_at_z[i] = PairingsBn254.new_fr(
                serialized_proof[j]
            );

            j = j.uncheckedInc();
        }        
        {{/if}}
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
        {{#if has_lookup}}
        proof.lookup_s_poly_opening_at_z_omega = PairingsBn254.new_fr(
                serialized_proof[j]
            );
        j = j.uncheckedInc();
        proof.lookup_grand_product_opening_at_z_omega = PairingsBn254.new_fr(
                serialized_proof[j]
            );

        j = j.uncheckedInc();
        proof.lookup_t_poly_opening_at_z = PairingsBn254.new_fr(
                serialized_proof[j]
            );

        j = j.uncheckedInc();
        proof.lookup_t_poly_opening_at_z_omega = PairingsBn254.new_fr(
                serialized_proof[j]
            );
        j = j.uncheckedInc();
        proof.lookup_selector_poly_opening_at_z = PairingsBn254.new_fr(
                serialized_proof[j]
            );
        j = j.uncheckedInc();
        proof.lookup_table_type_poly_opening_at_z = PairingsBn254.new_fr(
                serialized_proof[j]
            );
        j = j.uncheckedInc();
        {{/if}}
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
