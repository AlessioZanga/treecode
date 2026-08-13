use std::process::Command;

fn run_treecode(args: &[&str]) -> String {
    let output = Command::new("external/treecode/treecode")
        .args(args)
        .output()
        .expect("Failed to execute treecode");
    String::from_utf8_lossy(&output.stdout).to_string()
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

#[test]
fn test_nbody100_energy_conservation() {
    let output = run_treecode(&["nbody=100", "tstop=0.1", "dtout=0.05"]);
    let energies = parse_energy_lines(&output);

    assert!(
        energies.len() >= 2,
        "Not enough energy output lines: {}",
        energies.len()
    );

    let e0 = energies[0].1;
    let ef = energies.last().unwrap().1;
    let energy_drift = ((ef - e0) / e0).abs();
    assert!(
        energy_drift < 0.01,
        "Energy drift too large: {} (E0={}, Ef={})",
        energy_drift,
        e0,
        ef
    );
}

#[test]
fn test_nbody50_energy_conservation() {
    let output = run_treecode(&["nbody=50", "tstop=0.05", "dtout=0.025"]);
    let energies = parse_energy_lines(&output);

    assert!(energies.len() >= 2, "Not enough energy output lines");

    let e0 = energies[0].1;
    let ef = energies.last().unwrap().1;
    let energy_drift = ((ef - e0) / e0).abs();
    assert!(
        energy_drift < 0.01,
        "Energy drift too large: {}",
        energy_drift
    );
}

#[test]
fn test_nbody100_reference_output_matches() {
    let output = run_treecode(&["nbody=100", "tstop=0.1", "dtout=0.05"]);
    let reference = std::fs::read_to_string("tests/reference_outputs/nbody100_tstop0.1.txt")
        .expect("Failed to read reference output");

    let output_lines: Vec<&str> = output.lines().collect();
    let reference_lines: Vec<&str> = reference.lines().collect();

    assert_eq!(
        output_lines.len(),
        reference_lines.len(),
        "Output has different number of lines"
    );

    for (i, (out_line, ref_line)) in output_lines.iter().zip(reference_lines.iter()).enumerate() {
        if out_line.trim().is_empty() || ref_line.trim().is_empty() {
            continue;
        }
        if out_line.contains("time") && out_line.contains("|T+U|") {
            continue;
        }
        if out_line.contains("rsize") && out_line.contains("tdepth") {
            continue;
        }
        assert_eq!(
            out_line.trim(),
            ref_line.trim(),
            "Line {} differs:\n  got:      {:?}\n  expected: {:?}",
            i + 1,
            out_line.trim(),
            ref_line.trim()
        );
    }
}

#[test]
fn test_nbody50_reference_output_matches() {
    let output = run_treecode(&["nbody=50", "tstop=0.05", "dtout=0.025"]);
    let reference = std::fs::read_to_string("tests/reference_outputs/nbody50_tstop0.05.txt")
        .expect("Failed to read reference output");

    let output_lines: Vec<&str> = output.lines().collect();
    let reference_lines: Vec<&str> = reference.lines().collect();

    assert_eq!(
        output_lines.len(),
        reference_lines.len(),
        "Output has different number of lines"
    );

    for (i, (out_line, ref_line)) in output_lines.iter().zip(reference_lines.iter()).enumerate() {
        if out_line.trim().is_empty() || ref_line.trim().is_empty() {
            continue;
        }
        if out_line.contains("time") && out_line.contains("|T+U|") {
            continue;
        }
        if out_line.contains("rsize") && out_line.contains("tdepth") {
            continue;
        }
        assert_eq!(
            out_line.trim(),
            ref_line.trim(),
            "Line {} differs:\n  got:      {:?}\n  expected: {:?}",
            i + 1,
            out_line.trim(),
            ref_line.trim()
        );
    }
}

#[test]
fn test_treecode_header_output() {
    let output = run_treecode(&["nbody=10", "tstop=0.01", "dtout=0.005"]);
    assert!(
        output.contains("Hierarchical N-body code"),
        "Missing header in output"
    );
    assert!(output.contains("nbody"), "Missing nbody in output");
    assert!(output.contains("10"), "Missing nbody value in output");
}

#[test]
fn test_treecode_diagnostics_present() {
    let output = run_treecode(&["nbody=20", "tstop=0.01", "dtout=0.005"]);
    assert!(output.contains("|T+U|"), "Missing energy diagnostic");
    assert!(output.contains("rsize"), "Missing tree diagnostic");
}
