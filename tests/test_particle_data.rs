use std::process::Command;

use approx::assert_relative_eq;
use tempfile::TempDir;

fn workspace_root() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn run_treecode_with_output(args: &[&str], out_dir: &std::path::Path) -> String {
    let root = workspace_root();
    let out_pattern = out_dir.join("out_%d.txt").to_str().unwrap().to_string();

    let mut full_args: Vec<String> = args.iter().map(|s| s.to_string()).collect();
    full_args.push(format!("out={}", out_pattern));

    let str_refs: Vec<&str> = full_args.iter().map(|s| s.as_str()).collect();

    let output = Command::new(root.join("external/treecode/treecode"))
        .args(&str_refs)
        .output()
        .expect("Failed to execute treecode");
    String::from_utf8_lossy(&output.stdout).to_string()
}

struct ParticleData {
    nbody: usize,
    ndim: usize,
    time: f64,
    masses: Vec<f64>,
    positions: Vec<[f64; 3]>,
    velocities: Vec<[f64; 3]>,
}

fn parse_particle_file(path: &std::path::Path) -> ParticleData {
    let content = std::fs::read_to_string(path).expect("Failed to read particle file");
    let mut lines = content.lines().filter(|l| !l.trim().is_empty());

    let nbody: usize = lines
        .next()
        .expect("Missing nbody")
        .trim()
        .parse()
        .expect("Invalid nbody");
    let ndim: usize = lines
        .next()
        .expect("Missing ndim")
        .trim()
        .parse()
        .expect("Invalid ndim");
    let time: f64 = lines
        .next()
        .expect("Missing time")
        .trim()
        .parse()
        .expect("Invalid time");

    let mut masses = Vec::with_capacity(nbody);
    for _ in 0..nbody {
        let line = lines.next().expect("Missing mass");
        let m: f64 = line.trim().parse().expect("Invalid mass");
        masses.push(m);
    }

    let mut positions = Vec::with_capacity(nbody);
    for _ in 0..nbody {
        let line = lines.next().expect("Missing position");
        let parts: Vec<f64> = line
            .split_whitespace()
            .map(|s| s.parse().expect("Invalid position component"))
            .collect();
        assert_eq!(parts.len(), ndim, "Position vector has wrong dimensions");
        positions.push([parts[0], parts[1], parts[2]]);
    }

    let mut velocities = Vec::with_capacity(nbody);
    for _ in 0..nbody {
        let line = lines.next().expect("Missing velocity");
        let parts: Vec<f64> = line
            .split_whitespace()
            .map(|s| s.parse().expect("Invalid velocity component"))
            .collect();
        assert_eq!(parts.len(), ndim, "Velocity vector has wrong dimensions");
        velocities.push([parts[0], parts[1], parts[2]]);
    }

    ParticleData {
        nbody,
        ndim,
        time,
        masses,
        positions,
        velocities,
    }
}

#[test]
fn test_particle_count_matches() {
    let tmp = TempDir::new().unwrap();
    run_treecode_with_output(&["nbody=30", "tstop=0.02", "dtout=0.01"], tmp.path());
    let out_file = tmp.path().join("out_0.txt");
    assert!(out_file.exists(), "Output file not created");

    let data = parse_particle_file(&out_file);
    assert_eq!(data.nbody, 30);
    assert_eq!(data.ndim, 3);
    assert_eq!(data.masses.len(), 30);
    assert_eq!(data.positions.len(), 30);
    assert_eq!(data.velocities.len(), 30);
}

#[test]
fn test_masses_equal() {
    let tmp = TempDir::new().unwrap();
    let nbody = 30;
    run_treecode_with_output(&["nbody=30", "tstop=0.02", "dtout=0.01"], tmp.path());
    let data = parse_particle_file(&tmp.path().join("out_0.txt"));

    let expected_mass = 1.0 / nbody as f64;
    for m in &data.masses {
        assert_relative_eq!(*m, expected_mass, epsilon = 1e-6);
    }
}

#[test]
fn test_time_at_output() {
    let tmp = TempDir::new().unwrap();
    run_treecode_with_output(&["nbody=30", "tstop=0.04", "dtout=0.01"], tmp.path());
    let data0 = parse_particle_file(&tmp.path().join("out_0.txt"));
    let data1 = parse_particle_file(&tmp.path().join("out_1.txt"));

    assert_relative_eq!(data0.time, 0.0, epsilon = 1e-6);
    assert_relative_eq!(data1.time, 0.03125, epsilon = 1e-4);
}

#[test]
fn test_positions_within_bounds() {
    let tmp = TempDir::new().unwrap();
    run_treecode_with_output(&["nbody=40", "tstop=0.02", "dtout=0.01"], tmp.path());
    let data = parse_particle_file(&tmp.path().join("out_0.txt"));

    for (i, pos) in data.positions.iter().enumerate() {
        let r: f64 = pos.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(
            r < 10.0,
            "Body {} position magnitude {} out of bounds",
            i,
            r
        );
    }
}

#[test]
fn test_velocities_finite() {
    let tmp = TempDir::new().unwrap();
    run_treecode_with_output(&["nbody=30", "tstop=0.02", "dtout=0.01"], tmp.path());
    let data = parse_particle_file(&tmp.path().join("out_0.txt"));

    for (i, vel) in data.velocities.iter().enumerate() {
        for (j, v) in vel.iter().enumerate() {
            assert!(
                v.is_finite(),
                "Body {} velocity component {} is not finite: {}",
                i,
                j,
                v
            );
        }
    }
}

#[test]
fn test_center_of_mass_near_zero() {
    let tmp = TempDir::new().unwrap();
    run_treecode_with_output(&["nbody=50", "tstop=0.02", "dtout=0.01"], tmp.path());
    let data = parse_particle_file(&tmp.path().join("out_0.txt"));

    let mut cm = [0.0f64; 3];
    let mut total_mass = 0.0;
    for (pos, m) in data.positions.iter().zip(data.masses.iter()) {
        cm[0] += pos[0] * m;
        cm[1] += pos[1] * m;
        cm[2] += pos[2] * m;
        total_mass += m;
    }
    cm[0] /= total_mass;
    cm[1] /= total_mass;
    cm[2] /= total_mass;

    let cm_mag: f64 = cm.iter().map(|x| x * x).sum::<f64>().sqrt();
    assert!(cm_mag < 1e-5, "Center of mass not near zero: {}", cm_mag);
}

#[test]
fn test_particle_data_consistency() {
    let tmp = TempDir::new().unwrap();
    run_treecode_with_output(&["nbody=30", "tstop=0.02", "dtout=0.01"], tmp.path());
    let data = parse_particle_file(&tmp.path().join("out_0.txt"));

    assert_eq!(data.nbody, 30);
    assert_eq!(data.ndim, 3);

    let total_mass: f64 = data.masses.iter().sum();
    assert_relative_eq!(total_mass, 1.0, epsilon = 1e-5);
}

#[test]
fn test_deterministic_output() {
    let tmp1 = TempDir::new().unwrap();
    let tmp2 = TempDir::new().unwrap();

    run_treecode_with_output(
        &["nbody=30", "tstop=0.02", "dtout=0.01", "seed=42"],
        tmp1.path(),
    );
    run_treecode_with_output(
        &["nbody=30", "tstop=0.02", "dtout=0.01", "seed=42"],
        tmp2.path(),
    );

    let content1 = std::fs::read_to_string(tmp1.path().join("out_0.txt")).unwrap();
    let content2 = std::fs::read_to_string(tmp2.path().join("out_0.txt")).unwrap();
    assert_eq!(
        content1, content2,
        "Same seed should produce identical output"
    );
}

#[test]
fn test_multiple_output_files() {
    let tmp = TempDir::new().unwrap();
    run_treecode_with_output(&["nbody=40", "tstop=0.05", "dtout=0.01"], tmp.path());

    let mut files: Vec<_> = std::fs::read_dir(tmp.path())
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().is_some_and(|ext| ext == "txt"))
        .collect();
    files.sort_by_key(|e| e.path());

    assert!(
        files.len() >= 3,
        "Should have at least 3 output files, got {}",
        files.len()
    );

    let mut prev_time = -1.0;
    for file in &files {
        let data = parse_particle_file(&file.path());
        assert!(
            data.time > prev_time,
            "Output times should be increasing: {} <= {}",
            data.time,
            prev_time
        );
        prev_time = data.time;
    }
}
