use std::os::unix::io::RawFd;
use std::path::{Path, PathBuf};
use std::process::Command;

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

fn is_diagnostic_line(l: &str) -> bool {
    l.starts_with(|ch: char| ch.is_ascii_digit()) && l.split_whitespace().count() >= 8
}

fn trim_cputot(l: &mut String) {
    let parts: Vec<&str> = l.split_whitespace().collect();
    if parts.len() >= 8 {
        *l = parts[..parts.len() - 1].join(" ");
    }
}

fn assert_logs_match(c: &str, r: &str, ctx: &str) {
    let cl: Vec<String> = c.lines().map(str::to_string).collect();
    let rl: Vec<String> = r.lines().map(str::to_string).collect();
    assert_eq!(cl.len(), rl.len(), "{}: line count differs", ctx);
    for (i, (a, b)) in cl.iter().zip(rl.iter()).enumerate() {
        let mut a = normalize_line(a);
        let mut b = normalize_line(b);
        if is_diagnostic_line(a.as_str()) {
            trim_cputot(&mut a);
            trim_cputot(&mut b);
        }
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
