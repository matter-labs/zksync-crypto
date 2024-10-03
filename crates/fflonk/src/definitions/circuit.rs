use super::*;

#[derive(Clone, Debug)]
pub struct MockCircuitWidth3;
impl Circuit<Bn256> for MockCircuitWidth3 {
    type MainGate = NaiveMainGate;

    fn synthesize<CS: bellman::plonk::better_better_cs::cs::ConstraintSystem<Bn256> + 'static>(
        &self,
        cs: &mut CS,
    ) -> Result<(), bellman::SynthesisError> {
        todo!()
    }
    fn declare_used_gates() -> Result<
        Vec<Box<dyn bellman::plonk::better_better_cs::cs::GateInternal<Bn256>>>,
        bellman::SynthesisError,
    > {
        Ok(vec![Self::MainGate::default().into_internal()])
    }
}
#[derive(Clone, Debug)]
pub struct MockCircuitWidth4DNextCustomGate;
impl Circuit<Bn256> for MockCircuitWidth4DNextCustomGate {
    type MainGate = SelectorOptimizedWidth4MainGateWithDNext;

    fn synthesize<CS: bellman::plonk::better_better_cs::cs::ConstraintSystem<Bn256> + 'static>(
        &self,
        cs: &mut CS,
    ) -> Result<(), bellman::SynthesisError> {
        todo!()
    }
    fn declare_used_gates() -> Result<
        Vec<Box<dyn bellman::plonk::better_better_cs::cs::GateInternal<Bn256>>>,
        bellman::SynthesisError,
    > {
        Ok(vec![
            Self::MainGate::default().into_internal(),
            Rescue5CustomGate::default().into_internal(),
        ])
    }
}
