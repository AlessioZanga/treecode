use std::{
    io::Cursor,
    os::unix::io::RawFd,
    path::{Path, PathBuf},
};

use treecode::vecmath::{Matrix, Vector, matrix_identity, matrix_zero, vector_length, vector_zero};

const STDOUT_FD: RawFd = 1;

struct StdoutCapture {
    saved: RawFd,
    _path: PathBuf,
}

impl StdoutCapture {
    fn new(dir: &Path) -> Self {
        let saved = unsafe { libc::dup(STDOUT_FD) };
        assert!(saved >= 0);
        let path = dir.join("stdout.txt");
        let file = std::fs::File::create(&path).unwrap();
        let fd = std::os::unix::io::IntoRawFd::into_raw_fd(file);
        unsafe {
            libc::dup2(fd, STDOUT_FD);
            libc::close(fd);
        }
        StdoutCapture { saved, _path: path }
    }
}

impl Drop for StdoutCapture {
    fn drop(&mut self) {
        unsafe {
            libc::dup2(self.saved, STDOUT_FD);
            libc::close(self.saved);
        }
    }
}

fn run_rust_in(dir: &Path, args: &[&str]) -> treecode::treecode::Tree {
    let cap = StdoutCapture::new(dir);
    let mut full: Vec<String> = vec!["treecode".to_string()];
    full.extend(args.iter().map(|s| s.to_string()));
    let refs: Vec<&str> = full.iter().map(|s| s.as_str()).collect();
    let tree = treecode::treecode::run(&refs).unwrap();
    drop(cap);
    tree
}

#[test]
fn tree_and_force_accessors() {
    let dir = tempfile::TempDir::new().unwrap();
    let tree = run_rust_in(
        dir.path(),
        &["nbody=60", "usequad=true", "tstop=0.02", "dtout=0.01"],
    );

    let depth = tree.tree_depth();
    assert!((2..=32).contains(&depth), "depth out of range: {}", depth);
    let ncell = tree.cell_count();
    assert!(ncell > 0);
    let build_time = tree.tree_build_time();
    assert!(build_time >= 0.0);

    let max_active = tree.force_max_active();
    assert!(max_active > 0);
    assert!(tree.force_bb_calc() >= 0);
    assert!(tree.force_bc_calc() >= 0);
    let cpu = tree.force_cpu_time();
    assert!(cpu >= 0.0);
}

#[test]
fn save_restore_roundtrip() {
    let dir = tempfile::TempDir::new().unwrap();
    let mut tree = run_rust_in(dir.path(), &["nbody=30", "tstop=0.02", "dtout=0.01"]);

    let state = dir.path().join("w.rst");
    let state_str = state.display().to_string();
    tree.savestate(&state_str).unwrap();
    assert!(state.exists());

    run_rust_in(dir.path(), &["nbody=30", "tstop=0.03", "dtout=0.01"]);
    tree.restorestate(&state_str).unwrap();
    assert_eq!(tree.nbody, 30);
}

#[test]
fn save_restore_roundtrip_via_writer() {
    let dir = tempfile::TempDir::new().unwrap();
    let tree = run_rust_in(dir.path(), &["nbody=30", "tstop=0.02", "dtout=0.01"]);

    let mass0 = tree.bodytab[0].bodynode.mass;

    // Exercise the injectable I/O core directly: serialize into a Vec<u8>
    // instead of opening a file on disk.
    let mut buf: Vec<u8> = Vec::new();
    tree.savestate_to(&mut buf).unwrap();

    let mut restored = treecode::treecode::Tree::new();
    restored.restorestate_from(&mut Cursor::new(buf)).unwrap();
    assert_eq!(restored.nbody, 30);
    assert_eq!(restored.bodytab[0].bodynode.mass, mass0);
}

#[test]
fn vector_matrix_helpers() {
    let mut v = Vector::from([1.0, 2.0, 3.0]);
    vector_zero(&mut v);
    assert_eq!(v, Vector::zero());
    assert!((vector_length(&Vector::from([3.0, 4.0, 0.0])) - 5.0).abs() < 1e-5);

    let mut m = Matrix::ones();
    matrix_zero(&mut m);
    assert_eq!(m, Matrix::zero());

    let mut i = Matrix::zero();
    matrix_identity(&mut i);
    assert_eq!(i[0][0], 1.0);
    assert_eq!(i[1][0], 0.0);
    assert_eq!(i[2][2], 1.0);
}

#[test]
fn simulation_api_new_and_run() {
    let mut sim = treecode::Simulation::new(["nbody=30", "tstop=0.01", "dtout=0.005"]).unwrap();
    sim.run().unwrap();
    assert!(sim.state().nstep > 0, "simulation should advance steps");
}

#[test]
fn simulation_api_step_and_output() {
    let mut sim = treecode::Simulation::new(["nbody=20", "tstop=0.02", "dtout=0.01"]).unwrap();
    sim.step().unwrap();
    sim.output().unwrap();
    assert!(sim.state().nstep > 0);
}
