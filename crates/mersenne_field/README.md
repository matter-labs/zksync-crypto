# Mersenne field

This implementation of the Mersenne field uses the modulus \(2^{31} - 1\).

On top of regular methods, it has optimized methods for AVX-512F and riscv32.

## RiscV32 Extensions

When compiling this crate for riscv32, you can enable extra features that reduce the number of instructions. Note that these options rely on custom opcodes and specific hardware or simulator support.

- **use_division**: Attempts to use division instead of subtraction.
- **modular_ops**: Works with use_division to replace the standard "remu" instruction with custom MOP.RR opcodes.
- **modular_ext4_ops**: Uses the control status register 0x7c5 to speed up certain operations (currently FMA) on the Mersenne31Quartic.



