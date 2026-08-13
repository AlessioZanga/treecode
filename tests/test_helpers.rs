use std::os::unix::io::RawFd;
use std::path::{Path, PathBuf};

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

fn run_rust_in(dir: &Path, args: &[&str]) {
    let cap = StdoutCapture::new(dir);
    let mut full: Vec<String> = vec!["treecode".to_string()];
    full.extend(args.iter().map(|s| s.to_string()));
    let refs: Vec<&str> = full.iter().map(|s| s.as_str()).collect();
    treecode::treecode::run(&refs);
    drop(cap);
}

#[test]
fn tree_and_force_accessors() {
    let dir = tempfile::TempDir::new().unwrap();
    run_rust_in(
        dir.path(),
        &["nbody=60", "usequad=true", "tstop=0.02", "dtout=0.01"],
    );

    let depth = treecode::treeload::tree_depth();
    assert!((2..=32).contains(&depth), "depth out of range: {}", depth);
    let ncell = treecode::treeload::cell_count();
    assert!(ncell > 0);
    let build_time = treecode::treeload::tree_build_time();
    assert!(build_time >= 0.0);

    let max_active = treecode::treegrav::force_max_active();
    assert!(max_active > 0);
    assert!(treecode::treegrav::force_bb_calc() >= 0);
    assert!(treecode::treegrav::force_bc_calc() >= 0);
    let cpu = treecode::treegrav::force_cpu_time();
    assert!(cpu >= 0.0);
}

#[test]
fn save_restore_roundtrip() {
    let dir = tempfile::TempDir::new().unwrap();
    run_rust_in(dir.path(), &["nbody=30", "tstop=0.02", "dtout=0.01"]);

    let state = dir.path().join("w.rst");
    let state_str = state.display().to_string();
    treecode::treeio::savestate(&state_str);
    assert!(state.exists());

    run_rust_in(dir.path(), &["nbody=30", "tstop=0.03", "dtout=0.01"]);
    treecode::treeio::restorestate(&state_str);
    let nbody = unsafe { treecode::types::nbody };
    assert_eq!(nbody, 30);
}

#[test]
fn vector_matrix_helpers() {
    let mut v = [1.0, 2.0, 3.0];
    treecode::types::vector_zero(&mut v);
    assert_eq!(v, [0.0; 3]);
    assert!((treecode::types::vector_length(&[3.0, 4.0, 0.0]) - 5.0).abs() < 1e-5);

    let mut m = [[1.0; 3]; 3];
    treecode::types::matrix_zero(&mut m);
    assert_eq!(m, [[0.0; 3]; 3]);

    let mut i = [[0.0; 3]; 3];
    treecode::types::matrix_identity(&mut i);
    assert_eq!(i[0][0], 1.0);
    assert_eq!(i[1][0], 0.0);
    assert_eq!(i[2][2], 1.0);
}
