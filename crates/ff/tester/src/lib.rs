#![feature(llvm_asm)]
#![feature(asm)]

extern crate ff;
extern crate rand;

mod test_short_field;
mod test_large_field;
mod test_large_cios_field;
mod check_cios;
mod check_assembly_4;
pub mod mul_variant0;
pub mod assembly_4;

pub mod adx_4;