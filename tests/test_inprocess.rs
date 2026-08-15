use std::{path::Path, process::Command};

fn run_rust_in(dir: &Path, args: &[&str]) -> String {
    run_rust_with_prog(dir, "treecode", args)
}

fn run_rust_with_prog(dir: &Path, _prog: &str, args: &[&str]) -> String {
    // Spawn the built binary and capture its stdout. This is cross-platform
    // (unlike fd-1 redirection) and exercises the real `main` entry point.
    let output = Command::new(env!("CARGO_BIN_EXE_treecode"))
        .current_dir(dir)
        .args(args)
        .output()
        .expect("failed to run treecode binary");
    assert!(
        output.status.success(),
        "treecode failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8_lossy(&output.stdout).into_owned()
}

#[test]
fn inprocess_plummer_default() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(dir.path(), &["nbody=100", "tstop=0.05", "dtout=0.025"]);
    assert!(out.contains("Hierarchical N-body code"));
    assert!(out.contains("|T+U|"));
    assert!(out.contains("rsize"));
    assert!(out.contains("nbbtot"));
}

#[test]
fn inprocess_plummer_usequad() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(
        dir.path(),
        &["nbody=60", "usequad=true", "tstop=0.04", "dtout=0.02"],
    );
    assert!(out.contains("true"));
    assert!(out.contains("|T+U|"));
}

#[test]
fn inprocess_small_theta() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(
        dir.path(),
        &["nbody=40", "theta=0.5", "tstop=0.03", "dtout=0.015"],
    );
    assert!(out.contains("|T+U|"));
}

#[test]
fn inprocess_large_theta() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(
        dir.path(),
        &["nbody=40", "theta=1.5", "tstop=0.03", "dtout=0.015"],
    );
    assert!(out.contains("|T+U|"));
}

#[test]
fn inprocess_smallest_nbody() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(dir.path(), &["nbody=8", "tstop=0.01", "dtout=0.005"]);
    assert!(out.contains("|T+U|"));
}

#[test]
fn inprocess_different_seed() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(
        dir.path(),
        &["nbody=30", "seed=42", "tstop=0.02", "dtout=0.01"],
    );
    assert!(out.contains("|T+U|"));
}

#[test]
fn inprocess_fractional_dtime() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(
        dir.path(),
        &["nbody=20", "dtime=1/64", "tstop=0.02", "dtout=0.01"],
    );
    assert!(out.contains("|T+U|"));
}

#[test]
fn inprocess_usequad_with_options() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(
        dir.path(),
        &[
            "nbody=40",
            "usequad=true",
            "options=out=phi,out=acc",
            "tstop=0.02",
            "dtout=0.01",
        ],
    );
    assert!(out.contains("options:"));
    assert!(out.contains("|T+U|"));
}

#[test]
fn inprocess_input_from_file() {
    let dir = tempfile::TempDir::new().unwrap();
    let input = dir.path().join("snap_0.txt");
    write_snapshot_input(&input, 20);
    let out = run_rust_in(
        dir.path(),
        &[
            &format!("in={}", input.display()),
            "tstop=0.02",
            "dtout=0.01",
        ],
    );
    assert!(out.contains("|T+U|"));
}

#[test]
fn inprocess_save_restore_roundtrip() {
    let dir = tempfile::TempDir::new().unwrap();
    let save = dir.path().join("st.rst");
    let restore = dir.path().join("st.rst");
    let save_str = format!("save={}", save.display());
    let restore_str = format!("restore={}", restore.display());

    let out1 = run_rust_in(
        dir.path(),
        &["nbody=40", "tstop=0.02", "dtout=0.01", &save_str],
    );
    assert!(out1.contains("|T+U|"));
    assert!(save.exists(), "save file not created");

    // Restore with a different `VERSION` than the one recorded in the save file
    // so the "state file may be outdated" warning is emitted. This is portable:
    // unlike the Unix fd-redirect approach, both runs spawn the real binary.
    let out2 = run_rust_with_prog(
        dir.path(),
        "restored",
        &[
            "nbody=40",
            "tstop=0.04",
            "dtout=0.01",
            "VERSION=2.0",
            &restore_str,
        ],
    );
    assert!(out2.contains("|T+U|"));
    eprintln!("OUT2>>>\n{}<<<OUT2", out2);
    assert!(out2.contains("warning: state file may be outdated"));
}

#[test]
fn inprocess_energy_conservation() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(dir.path(), &["nbody=80", "tstop=0.05", "dtout=0.025"]);
    let energies = parse_energy_lines(&out);
    assert!(energies.len() >= 2);
    let drift = ((energies.last().unwrap().1 - energies[0].1) / energies[0].1).abs();
    assert!(drift < 0.01, "energy drift too large: {}", drift);
}

#[test]
fn inprocess_momentum_conservation() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(dir.path(), &["nbody=80", "tstop=0.05", "dtout=0.025"]);
    let diags = parse_diagnostics(&out);
    assert!(!diags.is_empty());
    let vcom_max = diags.iter().map(|d| d.0).fold(0.0f64, f64::max);
    assert!(vcom_max < 0.1, "COM velocity too large: {}", vcom_max);
}

fn parse_energy_lines(output: &str) -> Vec<(f64, f64, f64, f64)> {
    let mut energies = Vec::new();
    let lines: Vec<&str> = output.lines().collect();
    for i in 0..lines.len() {
        if lines[i].contains("time") && lines[i].contains("|T+U|") && i + 1 < lines.len() {
            let data_line = lines[i + 1];
            let parts: Vec<&str> = data_line.split_whitespace().collect();
            if parts.len() >= 6 {
                if let (Ok(time), Ok(t_plus_u), Ok(t), Ok(neg_u)) = (
                    parts[0].parse::<f64>(),
                    parts[1].parse::<f64>(),
                    parts[2].parse::<f64>(),
                    parts[3].parse::<f64>(),
                ) {
                    energies.push((time, t_plus_u, t, neg_u));
                }
            }
        }
    }
    energies
}

fn parse_diagnostics(output: &str) -> Vec<(f64, f64)> {
    let mut result = Vec::new();
    let lines: Vec<&str> = output.lines().collect();
    for i in 0..lines.len() {
        if lines[i].contains("time") && lines[i].contains("|T+U|") && i + 1 < lines.len() {
            let data_line = lines[i + 1];
            let parts: Vec<&str> = data_line.split_whitespace().collect();
            if parts.len() >= 8 {
                if let (Ok(_t), Ok(_e), Ok(_x), Ok(_y), Ok(_z), Ok(vcom), Ok(jtot)) = (
                    parts[0].parse::<f64>(),
                    parts[1].parse::<f64>(),
                    parts[2].parse::<f64>(),
                    parts[3].parse::<f64>(),
                    parts[4].parse::<f64>(),
                    parts[5].parse::<f64>(),
                    parts[6].parse::<f64>(),
                ) {
                    result.push((vcom, jtot));
                }
            }
        }
    }
    result
}

fn write_snapshot_input(path: &Path, nbody: usize) {
    use std::io::Write;
    let mut content = String::new();
    content.push_str(&format!("{}\n", nbody));
    content.push_str("3\n");
    content.push_str("0.0\n");
    for _ in 0..nbody {
        content.push_str("1.0\n");
    }
    for i in 0..nbody {
        let x = (i as f64 * 0.1371).fract();
        let y = (i as f64 * 0.1577).fract();
        let z = (i as f64 * 0.1733).fract();
        content.push_str(&format!("{:.6} {:.6} {:.6}\n", x, y, z));
    }
    for _ in 0..nbody {
        content.push_str("0.0 0.0 0.0\n");
    }
    let mut f = std::fs::File::create(path).unwrap();
    f.write_all(content.as_bytes()).unwrap();
}

#[test]
fn inprocess_real_snapshot() {
    let dir = tempfile::TempDir::new().unwrap();
    let input = dir.path().join("snap_0.txt");
    write_snapshot_input(&input, 20);
    let out = run_rust_in(
        dir.path(),
        &[
            &format!("in={}", input.display()),
            "tstop=0.01",
            "dtout=0.005",
        ],
    );
    assert!(out.contains("|T+U|"));
}

#[test]
fn inprocess_particle_output() {
    let dir = tempfile::TempDir::new().unwrap();
    let out_pat = dir.path().join("snap_%d.txt").display().to_string();
    let out = run_rust_in(
        dir.path(),
        &[
            "nbody=50",
            "tstop=0.04",
            "dtout=0.02",
            "usequad=true",
            "options=out-phi,out-acc",
            &format!("out={}", out_pat),
        ],
    );
    assert!(out.contains("data output to file"));
    let snap0 = dir.path().join("snap_0.txt");
    let snap1 = dir.path().join("snap_1.txt");
    assert!(snap0.exists(), "snap_0 not created");
    assert!(snap1.exists(), "snap_1 not created");
    let text = std::fs::read_to_string(&snap0).unwrap();
    assert!(text.trim_start().starts_with("50"));
    assert!(text.contains('E'), "expected scientific notation in output");
}

#[test]
fn inprocess_bh86_criterion() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(
        dir.path(),
        &["nbody=40", "options=bh86", "tstop=0.03", "dtout=0.015"],
    );
    assert!(out.contains("|T+U|"));
}

#[test]
fn inprocess_sw94_criterion() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(
        dir.path(),
        &["nbody=40", "options=sw94", "tstop=0.03", "dtout=0.015"],
    );
    assert!(out.contains("|T+U|"));
}

#[test]
fn inprocess_dtout_fraction() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(
        dir.path(),
        &["nbody=40", "dtime=1/64", "dtout=1/16", "tstop=0.02"],
    );
    assert!(out.contains("|T+U|"));
}

#[test]
fn inprocess_new_tout_option() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(
        dir.path(),
        &["nbody=40", "options=new-tout", "tstop=0.03", "dtout=0.015"],
    );
    assert!(out.contains("|T+U|"));
}

#[test]
fn inprocess_save_with_percent_pattern() {
    let dir = tempfile::TempDir::new().unwrap();
    let save_pat = dir.path().join("st_%d.rst").display().to_string();
    let out = run_rust_in(
        dir.path(),
        &[
            "nbody=40",
            "tstop=0.02",
            "dtout=0.01",
            &format!("save={}", save_pat),
        ],
    );
    assert!(out.contains("|T+U|"));
    assert!(dir.path().join("st_0.rst").exists() || dir.path().join("st_1.rst").exists());
}

#[test]
fn inprocess_restore_with_overrides() {
    let dir = tempfile::TempDir::new().unwrap();
    let save = dir.path().join("st.rst");
    let save_str = format!("save={}", save.display());
    run_rust_in(
        dir.path(),
        &["nbody=40", "tstop=0.02", "dtout=0.01", &save_str],
    );
    let restore_str = format!("restore={}", save.display());
    let out = run_rust_in(
        dir.path(),
        &[
            &restore_str,
            "eps=0.05",
            "theta=0.5",
            "usequad=true",
            "tstop=0.04",
            "dtout=1/16",
            "options=new-tout",
        ],
    );
    assert!(out.contains("|T+U|"));
}

#[test]
fn inprocess_particle_output_append() {
    let dir = tempfile::TempDir::new().unwrap();
    let outfile = dir.path().join("out.txt").display().to_string();
    let out = run_rust_in(
        dir.path(),
        &[
            "nbody=30",
            "tstop=0.03",
            "dtout=0.015",
            &format!("out={}", outfile),
        ],
    );
    assert!(out.contains("data output to file"));
    let text = std::fs::read_to_string(&outfile).unwrap();
    let times = text.matches("\n").count();
    assert!(times > 2, "expected multiple snapshot groups appended");
}

#[test]
fn inprocess_decimal_dtime() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(
        dir.path(),
        &["nbody=40", "dtime=0.001", "tstop=0.002", "dtout=0.001"],
    );
    assert!(out.contains("|T+U|"));
}

#[test]
fn inprocess_reset_time_option() {
    let dir = tempfile::TempDir::new().unwrap();
    let out = run_rust_in(
        dir.path(),
        &["nbody=40", "options=reset-time", "tstop=0.02", "dtout=0.01"],
    );
    assert!(out.contains("|T+U|"));
}
