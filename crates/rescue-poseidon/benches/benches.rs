#[allow(dead_code)]
use criterion::Criterion;
mod hashes;

fn main() {
    let crit = &mut Criterion::default().configure_from_args();
    hashes::group(crit);
    crit.final_summary();
}
