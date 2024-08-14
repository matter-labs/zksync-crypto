mod test_large_field {
    use ff::*;
    #[PrimeFieldModulus = "21888242871839275222246405745257275088696311157297823662689037894645226208583"]
    #[PrimeFieldGenerator = "2"]
    #[UseADX = "true"]
    struct FrAsm(FrReprAsm);
}
