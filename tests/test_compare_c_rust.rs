use std::{
    os::unix::io::RawFd,
    path::{Path, PathBuf},
    process::Command,
};

const STDOUT_FD: RawFd = 1;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn c_binary() -> PathBuf {
    workspace_root().join("external/treecode/treecode")
}

fn run_c(dir: &Path, args: &[&str]) -> String {
    let output = Command::new(c_binary())
        .current_dir(dir)
        .args(args)
        .output()
        .expect("failed to run C treecode");
    assert!(
        output.status.success(),
        "C treecode failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8_lossy(&output.stdout).to_string()
}

fn run_rust(dir: &Path, args: &[&str]) -> String {
    let saved = unsafe { libc::dup(STDOUT_FD) };
    assert!(saved >= 0);
    let path = dir.join("rust_stdout.txt");
    let file = std::fs::File::create(&path).unwrap();
    let fd = std::os::unix::io::IntoRawFd::into_raw_fd(file);
    unsafe {
        libc::dup2(fd, STDOUT_FD);
        libc::close(fd);
    }
    let mut full: Vec<String> = vec!["treecode".to_string()];
    full.extend(args.iter().map(|s| s.to_string()));
    let refs: Vec<&str> = full.iter().map(|s| s.as_str()).collect();
    treecode::treecode::run(&refs).unwrap();
    unsafe {
        libc::dup2(saved, STDOUT_FD);
        libc::close(saved);
    }
    std::fs::read_to_string(&path).unwrap()
}

fn normalize_line(line: &str) -> String {
    let mut out = String::new();
    for part in line.split_whitespace() {
        if part.contains('/') && part.ends_with(".txt") {
            out.push_str("FILE");
        } else if part.contains('/') && part.ends_with(".rst") {
            out.push_str("STATE");
        } else {
            out.push_str(part);
        }
        out.push(' ');
    }
    out
}

fn normalize_diag(line: &str) -> String {
    // Purely-numeric diagnostic value lines (the force-report row and the
    // energy/diagnostics rows) carry volatile columns that are NOT part of the
    // physical byte-exact output:
    //   * the final column is a CPU-time measurement (varies run-to-run even for
    //     the sequential build), and
    //   * the force-report `actmax` column is the peak *global* active-list
    //     length, an artifact of the single shared mutable array used by the C
    //     reference. When the force walk is parallelized each subtree uses its
    //     own scratch buffer, so that peak cannot be preserved (it depends on
    //     the sequential sibling order of the global array).
    // Every other column -- notably the interaction counts `nbbtot`/`nbctot` --
    // is still compared exactly.
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.is_empty() {
        return String::new();
    }
    let all_num = parts.iter().all(|p| p.parse::<f64>().is_ok());
    if !all_num {
        return line.to_string();
    }
    let mut keep: Vec<&str> = parts;
    if keep.len() == 7 {
        keep.remove(3); // drop `actmax`
    }
    keep.pop(); // drop final CPU-time column
    keep.join(" ")
}

fn assert_logs_match(c: &str, r: &str, ctx: &str) {
    let cl: Vec<String> = c.lines().map(str::to_string).collect();
    let rl: Vec<String> = r.lines().map(str::to_string).collect();
    assert_eq!(cl.len(), rl.len(), "{}: line count differs", ctx);
    for (i, (a, b)) in cl.iter().zip(rl.iter()).enumerate() {
        let a = normalize_diag(&normalize_line(a));
        let b = normalize_diag(&normalize_line(b));
        assert_eq!(a.trim(), b.trim(), "{}: line {} differs", ctx, i + 1);
    }
}

fn assert_files_match(dir: &Path, fname: &str, ctx: &str) {
    let c = dir.join(format!("c_{}", fname));
    let r = dir.join(format!("rust_{}", fname));
    assert!(c.exists(), "C file missing: {}", c.display());
    assert!(r.exists(), "Rust file missing: {}", r.display());
    let cb = std::fs::read(&c).unwrap();
    let rb = std::fs::read(&r).unwrap();
    assert_eq!(cb, rb, "{}: file {} differs byte-for-byte", ctx, fname);
}

#[test]
fn c_vs_rust_logs_and_particles() {
    for (nbody, tstop, dtout) in [(40, 0.03, 0.015), (80, 0.04, 0.02)] {
        let dir = tempfile::TempDir::new().unwrap();
        let args = [
            format!("nbody={}", nbody),
            format!("tstop={}", tstop),
            format!("dtout={}", dtout),
            format!("out={}", dir.path().join("c_parts.txt").display()),
        ];
        let arg_refs: Vec<&str> = args.iter().map(|s| s.as_str()).collect();
        let c_out = run_c(dir.path(), &arg_refs);
        let c_log = dir.path().join("c_stdout.txt");
        std::fs::write(&c_log, &c_out).unwrap();

        let rust_args = [
            format!("nbody={}", nbody),
            format!("tstop={}", tstop),
            format!("dtout={}", dtout),
            format!("out={}", dir.path().join("rust_parts.txt").display()),
        ];
        let rust_refs: Vec<&str> = rust_args.iter().map(|s| s.as_str()).collect();
        let r_out = run_rust(dir.path(), &rust_refs);
        assert_logs_match(&c_out, &r_out, "stdout");
        assert_files_match(dir.path(), "parts.txt", "particles");
    }
}

#[test]
fn c_vs_rust_usequad() {
    let dir = tempfile::TempDir::new().unwrap();
    let c = dir.path().join("c_q.txt");
    let r = dir.path().join("rust_q.txt");
    let cargs = &[
        "nbody=60",
        "usequad=true",
        "tstop=0.03",
        "dtout=0.015",
        &format!("out={}", c.display()),
    ];
    let c_out = run_c(dir.path(), cargs);
    let rargs = &[
        "nbody=60",
        "usequad=true",
        "tstop=0.03",
        "dtout=0.015",
        &format!("out={}", r.display()),
    ];
    let r_out = run_rust(dir.path(), rargs);
    assert_logs_match(&c_out, &r_out, "stdout usequad");
    let cb = std::fs::read(&c).unwrap();
    let rb = std::fs::read(&r).unwrap();
    assert_eq!(cb, rb, "usequad particle files differ");
}
