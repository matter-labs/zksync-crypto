pragma solidity >=0.8.25 <0.9.0;
import {Proof} from "../src/Plonk4VerifierWithAccessToDNext.sol";

library HardcodedValues{
    function hardcoded_proof() internal pure returns(Proof memory proof){
    {{#each proof.inputs}}
        proof.input_values[{{@index}}] = {{this}};
    {{/each}}
    {{#each proof.state_polys_commitments}}
        proof.state_polys_commitments[{{@index}}].X = {{this.[0]}};
        proof.state_polys_commitments[{{@index}}].Y = {{this.[1]}};
    {{/each}}
        proof.copy_permutation_grand_product_commitment.X = {{proof.copy_permutation_grand_product_commitment.[0]}};
        proof.copy_permutation_grand_product_commitment.Y = {{proof.copy_permutation_grand_product_commitment.[1]}};
    {{#each proof.quotient_poly_parts_commitments}}
        proof.quotient_poly_parts_commitments[{{@index}}].X = {{this.[0]}};
        proof.quotient_poly_parts_commitments[{{@index}}].Y = {{this.[1]}};
    {{/each}}
        {{#if has_lookup}}
        proof.lookup_s_poly_commitment.X = {{proof.lookup_s_poly_commitment.[0]}};
        proof.lookup_s_poly_commitment.Y = {{proof.lookup_s_poly_commitment.[1]}};
        proof.lookup_grand_product_commitment.X = {{proof.lookup_grand_product_commitment.[0]}};
        proof.lookup_grand_product_commitment.Y = {{proof.lookup_grand_product_commitment.[1]}};
        {{/if}}
        proof.opening_proof_at_z.X = {{proof.opening_proof_at_z.[0]}};
        proof.opening_proof_at_z.Y = {{proof.opening_proof_at_z.[1]}};
        proof.opening_proof_at_z_omega.X = {{proof.opening_proof_at_z_omega.[0]}};
        proof.opening_proof_at_z_omega.Y = {{proof.opening_proof_at_z_omega.[1]}};

    {{#each proof.state_polys_openings_at_z}}
        proof.state_polys_openings_at_z[{{@index}}].value = {{this}};
    {{/each}}
    {{#each proof.state_polys_openings_at_z_omega}}
        proof.state_polys_openings_at_z_omega[{{@index}}].value = {{this}};
    {{/each}}
    {{#each proof.gate_selectors_openings_at_z}}
        proof.gate_selectors_openings_at_z[{{@index}}].value = {{this}};
    {{/each}}
    {{#each proof.copy_permutation_polys_openings_at_z}}
        proof.copy_permutation_polys_openings_at_z[{{@index}}].value = {{this}};
    {{/each}}
        proof.copy_permutation_grand_product_opening_at_z_omega.value = {{proof.copy_permutation_grand_product_opening_at_z_omega}};
        proof.quotient_poly_opening_at_z.value = {{proof.quotient_poly_opening_at_z}};
        proof.linearization_poly_opening_at_z.value = {{proof.linearization_poly_opening_at_z}};
    {{#if has_lookup}}
        proof.lookup_s_poly_opening_at_z_omega.value = {{proof.lookup_s_poly_opening_at_z_omega}};
        proof.lookup_grand_product_opening_at_z_omega.value = {{proof.lookup_grand_product_opening_at_z_omega}};
        proof.lookup_t_poly_opening_at_z.value = {{proof.lookup_t_poly_opening_at_z}};
        proof.lookup_selector_poly_opening_at_z.value = {{proof.lookup_selector_poly_opening_at_z}};
        proof.lookup_table_type_poly_opening_at_z.value = {{proof.lookup_table_type_poly_opening_at_z}};
    {{/if}}
    }    

    function hardcoded_serialized_proof() public pure returns(uint256[] memory serializedInputs, uint256[] memory serializedProof){
        serializedInputs = new uint256[]({{num_inputs}});
    {{#each serialized_inputs}}
        serializedInputs[{{@index}}] = {{this}};
    {{/each}}
        serializedProof = new uint256[]({{serialized_proof_length}});
    {{#each serialized_proof}}
        serializedProof[{{@index}}] = {{this}};
    {{/each}}
    }
}
