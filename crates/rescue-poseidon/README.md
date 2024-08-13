# Rescue and Poseidon
## Overview
This repo contains implementations of arithmetization oriented hash functions(Rescue, Poseidon, Rescue Prime) that constructed by a sponge construction over prime field for both out-of circuits and in-circuit usages. Each algebraic hash function uses same sponge construction with different round function or permutation function. Gadgets are optimal in the constraint systems while also supporting different scalar fields which supported by bellman. 

## Usage
Add dependency
```toml
rescue_poseidon = 0.1
```

```rust
    use franklin_crypto::bellman::bn256::Fr;
    use franklin_crypto::bellman::Field;
    use rescue_poseidon::rescue_hash;

    const L: usize = 2;
    let input = [Fr::one(); L]; // dummy input

    // fixed length rescue hash
    let result = rescue_hash::<Bn256, L>(&input);
    assert_eq!(result.len(), 2);
```
More examples can be found in `examples` folder.


## Testing
`cargo test -- --nocapture`

## Benchmarks & Constraint System Costs
`cargo bench -- --nocapture`


_CPU: 3,1 GHz Intel Core i5_

| hashes    | 1x permutation runtime (μs) | 1x permutation gates | number of rounds |
| --- | -------- | -------- | -------- |
| Poseidon   | 13     | 166     | 8f + 33p     |
| Rescue   | 680     | 266     | 44f     |
| Rescue Prime   | 300     | 104     | 9f     |



## References
- [1] [Cryprographic sponge functions](https://keccak.team/files/CSF-0.1.pdf)
- [2] [The sponge and duplex constructions](https://keccak.team/sponge_duplex.html)
- [3] [STARK Friendly Hash – Survey and Recommendation](https://eprint.iacr.org/2020/948.pdf)
- [4] [MARVELlous: a STARK-Friendly Family of Cryptographic Primitives](https://eprint.iacr.org/2018/1098.pdf)
- [5] [POSEIDON: A New Hash Function for Zero-Knowledge Proof Systems](https://eprint.iacr.org/2019/458.pdf)
- [6] [Rescue-Prime: a Standard Specification (SoK)](https://eprint.iacr.org/2020/1143.pdf)