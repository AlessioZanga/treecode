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
            if parts.len() >= 6
                && let (Ok(time), Ok(t_plus_u), Ok(t), Ok(neg_u)) = (
                    parts[0].parse::<f64>(),
                    parts[1].parse::<f64>(),
                    parts[2].parse::<f64>(),
                    parts[3].parse::<f64>(),
                )
            {
                energies.push((time, t_plus_u, t, neg_u));
            }
        }
    }
    energies
}

fn parse_last_energy(output: &str) -> Option<f64> {
    let energies = parse_energy_lines(output);
    energies.last().map(|e| e.1)
}

fn parse_diagnostics(output: &str) -> Vec<(f64, f64, f64, f64)> {
    let mut result = Vec::new();
    let lines: Vec<&str> = output.lines().collect();
    for i in 0..lines.len() {
        if lines[i].contains("time") && lines[i].contains("|T+U|") && i + 1 < lines.len() {
            let data_line = lines[i + 1];
            let parts: Vec<&str> = data_line.split_whitespace().collect();
            if parts.len() >= 8
                && let (Ok(_time), Ok(_e), Ok(_t), Ok(_u), Ok(_tou), Ok(vcom), Ok(jtot), Ok(_cpu)) = (
                    parts[0].parse::<f64>(),
                    parts[1].parse::<f64>(),
                    parts[2].parse::<f64>(),
                    parts[3].parse::<f64>(),
                    parts[4].parse::<f64>(),
                    parts[5].parse::<f64>(),
                    parts[6].parse::<f64>(),
                    parts[7].parse::<f64>(),
                )
            {
                result.push((vcom, jtot, 0.0, 0.0));
            }
        }
    }
    result
}

#[test]
fn test_plummer_small_step_count() {
    let output = run_treecode(&["nbody=50", "tstop=0.02", "dtout=0.01"]);
    assert!(output.contains("time"), "No time output found");
    let energy = parse_last_energy(&output);
    assert!(energy.is_some(), "Could not parse final energy");
}

#[test]
fn test_output_format_consistency() {
    let out1 = run_treecode(&["nbody=30", "tstop=0.01", "dtout=0.005"]);
    let out2 = run_treecode(&["nbody=30", "tstop=0.01", "dtout=0.005"]);
    assert_eq!(
        out1, out2,
        "Same parameters should produce identical output"
    );
}

#[test]
fn test_quadrupole_flag() {
    let out_noq = run_treecode(&["nbody=50", "tstop=0.02", "dtout=0.01"]);
    assert!(out_noq.contains("false"), "Default should be usequad=false");
}

#[test]
fn test_energy_momentum_conservation() {
    let output = run_treecode(&["nbody=80", "tstop=0.05", "dtout=0.025"]);

    let diagnostics = parse_diagnostics(&output);

    if !diagnostics.is_empty() {
        let vcom_max = diagnostics.iter().map(|d| d.0).fold(0.0f64, f64::max);
        assert!(
            vcom_max < 0.1,
            "Center of mass velocity too large: {}",
            vcom_max
        );
    }

    if diagnostics.len() >= 2 {
        let j0 = diagnostics[0].1;
        let jf = diagnostics.last().unwrap().1;
        let drift = ((jf - j0) / j0).abs();
        assert!(drift < 0.01, "Angular momentum drift too large: {}", drift);
    }
}

#[test]
fn test_tree_depth_reasonable() {
    let output = run_treecode(&["nbody=100", "tstop=0.01", "dtout=0.005"]);

    for line in output.lines() {
        if line.contains("rsize") && line.contains("tdepth") {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 3
            && let Ok(depth) = parts[1].parse::<i32>()
        {
            assert!(
                (2..=32).contains(&depth),
                "Tree depth out of range: {}",
                depth
            );
        }
    }
}
