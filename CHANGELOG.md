# Changelog

## [0.30.13](https://github.com/matter-labs/zksync-crypto/compare/v0.30.12...v0.30.13) (2025-01-11)


### Features

* **bellman:** CRS constructor with custom allocator ([#62](https://github.com/matter-labs/zksync-crypto/issues/62)) ([296045b](https://github.com/matter-labs/zksync-crypto/commit/296045bcd948051fdfbb22a13c8bc72c5f2378b8))
* **fflonk:** non-default allocator feature ([#61](https://github.com/matter-labs/zksync-crypto/issues/61)) ([debf9a3](https://github.com/matter-labs/zksync-crypto/commit/debf9a35c090eb58a862c3390827a097cdc43061))

## [0.30.12](https://github.com/matter-labs/zksync-crypto/compare/v0.30.11...v0.30.12) (2024-12-20)


### Features

* FFLONK ([#58](https://github.com/matter-labs/zksync-crypto/issues/58)) ([a57cf0f](https://github.com/matter-labs/zksync-crypto/commit/a57cf0fc5cee14fa6a361f0d5990ed2de8b094f5))

## [0.30.11](https://github.com/matter-labs/zksync-crypto/compare/v0.30.10...v0.30.11) (2024-12-18)


### Bug Fixes

* **snark-wrapper:** call range check with proper params ([#54](https://github.com/matter-labs/zksync-crypto/issues/54)) ([cba8e9c](https://github.com/matter-labs/zksync-crypto/commit/cba8e9c4334646ef85caedf1584fc7bbcd4656c5))
* **snark-wrapper:** range check bitlen in snark wrapper circuit ([#56](https://github.com/matter-labs/zksync-crypto/issues/56)) ([4549fb2](https://github.com/matter-labs/zksync-crypto/commit/4549fb2fb2648537e4ac967c6bf9d7001ac93b69))

## [0.30.10](https://github.com/matter-labs/zksync-crypto/compare/v0.30.9...v0.30.10) (2024-11-27)


### Features

* **franklin-crypto:** add `extern crate core;` to lib.rs ([#27](https://github.com/matter-labs/zksync-crypto/issues/27)) ([3084cb8](https://github.com/matter-labs/zksync-crypto/commit/3084cb821965382d63cc6f5bc074cf6dcfaff84d))


### Bug Fixes

* Disabled vectorisation for neon ([#52](https://github.com/matter-labs/zksync-crypto/issues/52)) ([3900a4b](https://github.com/matter-labs/zksync-crypto/commit/3900a4b4b225545e9cc05ad7ebca570aac3dd300))

## [0.30.9](https://github.com/matter-labs/zksync-crypto/compare/v0.30.8...v0.30.9) (2024-11-21)


### Bug Fixes

* **boojum,snark_wrapper:** revert PoW fix ([#47](https://github.com/matter-labs/zksync-crypto/issues/47)) ([477bf1e](https://github.com/matter-labs/zksync-crypto/commit/477bf1e72e63ea758fc7d8520178a832c2eea0e3))

## [0.30.8](https://github.com/matter-labs/zksync-crypto/compare/v0.30.7...v0.30.8) (2024-11-19)


### Bug Fixes

* **snark_wrapper:** compute the correct number of PoW seed challenges ([#45](https://github.com/matter-labs/zksync-crypto/issues/45)) ([c2999e1](https://github.com/matter-labs/zksync-crypto/commit/c2999e11a9643f0c5174d106849085d2f908a3e2))

## [0.30.7](https://github.com/matter-labs/zksync-crypto/compare/v0.30.6...v0.30.7) (2024-11-18)


### Features

* **bellman:** Remove allocator feature from default features ([#39](https://github.com/matter-labs/zksync-crypto/issues/39)) ([9ecb5cd](https://github.com/matter-labs/zksync-crypto/commit/9ecb5cdbfa4f4a1157c349b56cc807c4a842ad49))


### Bug Fixes

* **boojum:** compute the correct number of PoW seed challenges ([#43](https://github.com/matter-labs/zksync-crypto/issues/43)) ([8d2f5f7](https://github.com/matter-labs/zksync-crypto/commit/8d2f5f74c7a7a22db2ca3212a2fc6653b3ee0c76))

## [0.30.6](https://github.com/matter-labs/zksync-crypto/compare/v0.30.5...v0.30.6) (2024-10-31)


### Bug Fixes

* **franklin-crypto:** range check goldilocks with naive gate ([#37](https://github.com/matter-labs/zksync-crypto/issues/37)) ([450cdfe](https://github.com/matter-labs/zksync-crypto/commit/450cdfe4cb2d6f1ffd616744129d217d5cec0126))

## [0.30.5](https://github.com/matter-labs/zksync-crypto/compare/v0.30.4...v0.30.5) (2024-10-31)


### Bug Fixes

* **fflonk:** fix Cargo.toml ([#35](https://github.com/matter-labs/zksync-crypto/issues/35)) ([86cce2b](https://github.com/matter-labs/zksync-crypto/commit/86cce2b833f3a4da0ba2bb3fa1c994447b0389bc))

## [0.30.4](https://github.com/matter-labs/zksync-crypto/compare/v0.30.3...v0.30.4) (2024-10-30)


### Features

* fflonk protocol implementation  ([#11](https://github.com/matter-labs/zksync-crypto/issues/11)) ([a1485ce](https://github.com/matter-labs/zksync-crypto/commit/a1485ce53f1a92892c4845f02f0fc3416899bd92))

## [0.30.3](https://github.com/matter-labs/zksync-crypto/compare/v0.30.2...v0.30.3) (2024-10-29)


### Features

* **bellman:** declare naive main gate ([#28](https://github.com/matter-labs/zksync-crypto/issues/28)) ([5f563e0](https://github.com/matter-labs/zksync-crypto/commit/5f563e06a0c0c76c1c232ef041c359e7256d333c))
* **franklin-crypto:** naive main gate compatible gadgets ([#30](https://github.com/matter-labs/zksync-crypto/issues/30)) ([f44dd45](https://github.com/matter-labs/zksync-crypto/commit/f44dd45ce587326bb6f0a0b84ce6096e191ca298))
* **snark-wrapper:** wrapper circuit with naive main gate ([#29](https://github.com/matter-labs/zksync-crypto/issues/29)) ([235d0c8](https://github.com/matter-labs/zksync-crypto/commit/235d0c8481b7079a07ccb621745a230194bb00ce))

## [0.30.2](https://github.com/matter-labs/zksync-crypto/compare/v0.30.1...v0.30.2) (2024-10-29)


### Features

* **boojum:** add get_light_setup method to CSReferenceAssembly ([#31](https://github.com/matter-labs/zksync-crypto/issues/31)) ([6dde34c](https://github.com/matter-labs/zksync-crypto/commit/6dde34c119bf7f0ff91734d513adc8b265d17d16))

## [0.30.1](https://github.com/matter-labs/zksync-crypto/compare/v0.30.0...v0.30.1) (2024-09-05)


### Bug Fixes

* **snark-wrapper:** Revert "improvement of parameters" ([#24](https://github.com/matter-labs/zksync-crypto/issues/24)) ([cce1da3](https://github.com/matter-labs/zksync-crypto/commit/cce1da378761dd76271730ad154e6f5b8a7675bb))

## [0.30.0](https://github.com/matter-labs/zksync-crypto/compare/v0.29.0...v0.30.0) (2024-09-05)


### ⚠ BREAKING CHANGES

* Rename crypto crates and properly set metadata ([#21](https://github.com/matter-labs/zksync-crypto/issues/21))

### Features

* Rename crypto crates and properly set metadata ([#21](https://github.com/matter-labs/zksync-crypto/issues/21)) ([14f44d7](https://github.com/matter-labs/zksync-crypto/commit/14f44d7c3054e02fe8fbaa093a4548b4b5d2f5cf))

## [0.29.0](https://github.com/matter-labs/zksync-crypto/compare/v0.28.0...v0.29.0) (2024-09-04)


### Features

* **ci:** Introduce release-please and automatic publishing ([#19](https://github.com/matter-labs/zksync-crypto/issues/19)) ([a5444e3](https://github.com/matter-labs/zksync-crypto/commit/a5444e35f5074c0f0de6a9556c49682c228d92de))


### Bug Fixes

* Fix hard static analysis errors ([#18](https://github.com/matter-labs/zksync-crypto/issues/18)) ([29f0bdd](https://github.com/matter-labs/zksync-crypto/commit/29f0bddac058f0c460c36e914616252e9eee736e))
