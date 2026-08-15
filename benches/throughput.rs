//! Throughput benchmark comparing the Rust `Simulation` against the reference C
//! `treecode` binary. Run with `cargo bench`. Uses `harness = false` so it needs
//! no external benchmark crate; it simply times both implementations over a fixed
//! configuration and asserts the ratio stays within a sane window. The force
//! algorithm is byte-for-byte identical to the C reference, but Rust fans the
//! root's children out across threads, so it is expected (and intended) to be
//! substantially *faster* than C: the lower bound guards against a broken run
//! reporting an unrealistic speedup, the upper bound guards against a regression
//! that makes Rust much slower than C.

use std::{path::Path, process::Command, time::Instant};

use treecode::Simulation;

fn main() {
    let params = ["nbody=4096", "tstop=10.0", "dtout=0.025", "seed=123"];
    let iters = 3;

    // Warm up (also verifies the simulation runs end-to-end).
    Simulation::new(params).unwrap().run().unwrap();

    let start = Instant::now();
    for _ in 0..iters {
        Simulation::new(params).unwrap().run().unwrap();
    }
    let rust_ms = start.elapsed().as_secs_f64() * 1000.0 / iters as f64;
    println!("rust throughput : {rust_ms:8.2} ms/run");

    let cbin = concat!(env!("CARGO_MANIFEST_DIR"), "/external/treecode/treecode");
    if Path::new(cbin).exists() {
        let mut total = 0.0;
        for _ in 0..iters {
            let t = Instant::now();
            let out = Command::new(cbin)
                .args(params)
                .output()
                .expect("failed to run reference C treecode");
            assert!(out.status.success(), "C treecode failed");
            total += t.elapsed().as_secs_f64() * 1000.0;
        }
        let c_ms = total / iters as f64;
        println!("c    throughput : {c_ms:8.2} ms/run");

        let ratio = rust_ms / c_ms;
        println!("rust/c ratio    : {ratio:8.3}");
        assert!(
            (0.1..=2.0).contains(&ratio),
            "timing ratio {ratio:.3} outside the expected window (Rust should be \
             at most ~10x faster via parallel fan-out, and never more than 2x slower than C)"
        );
    } else {
        println!("reference C binary not found at {cbin}; skipping comparison");
    }

    println!(
        "interaction record: Interact={} B vs Cell={} B (≈{:.1}x smaller)",
        std::mem::size_of::<treecode::types::Interact>(),
        std::mem::size_of::<treecode::types::Cell>(),
        std::mem::size_of::<treecode::types::Cell>() as f64
            / std::mem::size_of::<treecode::types::Interact>() as f64,
    );
}
