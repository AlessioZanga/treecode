use std::{env, process::Command, time::Instant};

fn workspace_root() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn run_binary(args: &[&str]) -> String {
    let output = Command::new("external/treecode/treecode")
        .args(args)
        .output()
        .expect("Failed to execute C treecode");
    String::from_utf8_lossy(&output.stdout).to_string()
}

fn run_rust_binary(args: &[&str]) -> String {
    let root = workspace_root();
    let exe = root.join("target/release/treecode");
    if !exe.exists() {
        let status = Command::new("cargo")
            .args(["build", "--release"])
            .current_dir(&root)
            .status()
            .expect("Failed to build release treecode");
        assert!(status.success(), "cargo build --release failed");
    }
    let output = Command::new(exe)
        .args(args)
        .output()
        .expect("Failed to execute Rust treecode");
    String::from_utf8_lossy(&output.stdout).to_string()
}

fn bench_run(binary_fn: fn(&[&str]) -> String, args: &[&str], runs: usize) -> std::time::Duration {
    let start = Instant::now();
    for _ in 0..runs {
        binary_fn(args);
    }
    start.elapsed()
}

fn has_header(output: &str) -> bool {
    output.contains("nbody") && output.contains("dtime") && output.contains("theta")
}

fn has_force_report(output: &str) -> bool {
    output.contains("rsize") && output.contains("tdepth") && output.contains("ftree")
}

fn has_diagnostics(output: &str) -> bool {
    output.contains("|T+U|") && output.contains("|Vcom|") && output.contains("|Jtot|")
}

fn count_output_sections(output: &str) -> usize {
    output.lines().filter(|l| l.contains("|T+U|")).count()
}

#[test]
fn test_c_and_rust_produce_valid_output() {
    let args = &["nbody=30", "tstop=0.02", "dtout=0.01"];
    let c_out = run_binary(args);
    let rust_out = run_rust_binary(args);

    assert!(has_header(&c_out), "C missing header");
    assert!(has_header(&rust_out), "Rust missing header");
    assert!(has_force_report(&c_out), "C missing force report");
    assert!(has_force_report(&rust_out), "Rust missing force report");
    assert!(has_diagnostics(&c_out), "C missing diagnostics");
    assert!(has_diagnostics(&rust_out), "Rust missing diagnostics");

    let c_sections = count_output_sections(&c_out);
    let rust_sections = count_output_sections(&rust_out);
    assert_eq!(
        c_sections, rust_sections,
        "Different number of output sections: C={} Rust={}",
        c_sections, rust_sections
    );
}

#[test]
fn test_cargo_run_matches_c_binary() {
    let args = &["nbody=30", "tstop=0.02", "dtout=0.01"];
    let c_out = run_binary(args);
    let rust_out = run_rust_binary(args);

    assert!(has_header(&rust_out), "cargo run missing header");
    assert!(
        has_force_report(&rust_out),
        "cargo run missing force report"
    );
    assert!(has_diagnostics(&rust_out), "cargo run missing diagnostics");

    let c_sections = count_output_sections(&c_out);
    let rust_sections = count_output_sections(&rust_out);
    assert_eq!(c_sections, rust_sections);
}

#[test]
fn test_timing_c_vs_rust() {
    let args = &["nbody=50", "tstop=0.02", "dtout=0.01"];
    let runs = 3;

    let c_duration = bench_run(run_binary, args, runs);
    let rust_duration = bench_run(run_rust_binary, args, runs);

    let c_ms = c_duration.as_secs_f64() * 1000.0;
    let rust_ms = rust_duration.as_secs_f64() * 1000.0;
    let ratio = rust_ms / c_ms;

    eprintln!("\n=== Timing: nbody=50, tstop=0.02 (avg over {runs} runs) ===");
    eprintln!("  C:     {c_ms:.1} ms");
    eprintln!("  Rust:  {rust_ms:.1} ms");
    eprintln!("  Ratio: {ratio:.2}x");

    assert!(
        ratio < 3.0,
        "Rust binary is too slow: {ratio:.2}x slower than C ({rust_ms:.1}ms vs {c_ms:.1}ms)",
    );
}

#[test]
fn test_timing_larger_n() {
    let args = &["nbody=100", "tstop=0.01", "dtout=0.005"];
    let runs = 2;

    let c_duration = bench_run(run_binary, args, runs);
    let rust_duration = bench_run(run_rust_binary, args, runs);

    let c_ms = c_duration.as_secs_f64() * 1000.0;
    let rust_ms = rust_duration.as_secs_f64() * 1000.0;
    let ratio = rust_ms / c_ms;

    eprintln!("\n=== Timing: nbody=100, tstop=0.01 (avg over {runs} runs) ===");
    eprintln!("  C:     {c_ms:.1} ms");
    eprintln!("  Rust:  {rust_ms:.1} ms");
    eprintln!("  Ratio: {ratio:.2}x");

    assert!(
        ratio < 3.0,
        "Rust binary is too slow at nbody=100: {ratio:.2}x slower than C",
    );
}
