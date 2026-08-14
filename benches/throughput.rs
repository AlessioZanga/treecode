//! Throughput benchmark comparing the Rust `Simulation` against the reference C
//! `treecode` binary. Run with `cargo bench`. Uses `harness = false` so it needs
//! no external benchmark crate; it simply times both implementations over a fixed
//! configuration and asserts the ratio stays within the 0.5x–2.0x window the
//! port is expected to hold (the algorithm is unchanged, so timings track C).

use std::path::Path;
use std::process::Command;
use std::time::Instant;

use treecode::Simulation;

fn main() {
    let params = ["nbody=1500", "tstop=0.05", "dtout=0.025", "seed=123"];
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
            (0.5..=2.0).contains(&ratio),
            "timing ratio {ratio:.3} outside the 0.5x–2.0x window"
        );
    } else {
        println!("reference C binary not found at {cbin}; skipping comparison");
    }
}
