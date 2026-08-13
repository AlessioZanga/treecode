use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    let exe = std::env::current_exe().unwrap();
    let debug_dir = exe.parent().unwrap().parent().unwrap();
    debug_dir.join("treecode")
}

#[test]
fn subprocess_oom_path() {
    let out = Command::new(bin())
        .args(["nbody=40", "theta=0.0", "tstop=0.01", "dtout=0.005"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("out of memory"), "stderr: {}", stderr);
}

#[test]
fn subprocess_unknown_param() {
    let out = Command::new(bin())
        .args(["nbody=40", "frobnicate=1", "tstop=0.01", "dtout=0.005"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("unknown"), "stderr: {}", stderr);
}

#[test]
fn subprocess_help() {
    let out = Command::new(bin()).args(["-help"]).output().unwrap();
    assert!(out.status.success());
}

#[test]
fn subprocess_bad_input_file() {
    let dir = tempfile::TempDir::new().unwrap();
    let bad = dir.path().join("bad.txt");
    std::fs::write(&bad, "3 3\n1.0\n").unwrap();
    let out = Command::new(bin())
        .args([
            "nbody=40",
            &format!("in={}", bad.display()),
            "tstop=0.01",
            "dtout=0.005",
        ])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(!stderr.is_empty());
}

fn run_args(args: &[&str]) -> (i32, String, String) {
    let out = Command::new(bin()).args(args).output().unwrap();
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

#[test]
fn subprocess_clue_help() {
    let (c, clue, _) = run_args(&["-clue"]);
    assert_eq!(c, 0);
    assert!(clue.contains("nbody"));
    let (h, help, _) = run_args(&["-help"]);
    assert_eq!(h, 0);
    assert!(help.contains("Hierarchical"));
}

#[test]
fn subprocess_duplicate_param() {
    let (c, _, err) = run_args(&["nbody=10", "nbody=20"]);
    assert_ne!(c, 0);
    assert!(err.contains("duplicated"), "{}", err);
}

#[test]
fn subprocess_too_many_args() {
    let positionals: Vec<String> = (0..20).map(|i| format!("arg{}", i)).collect();
    let refs: Vec<&str> = positionals.iter().map(|s| s.as_str()).collect();
    let (c, _, err) = run_args(&refs);
    assert_ne!(c, 0);
    assert!(err.contains("too many"), "{}", err);
}

#[test]
fn subprocess_nameless_arg_after_options() {
    let (c, _, err) = run_args(&["nbody=40", "bogus"]);
    assert_ne!(c, 0);
    assert!(err.contains("nameless"), "{}", err);
}

#[test]
fn subprocess_bad_bool() {
    let (c, _, err) = run_args(&["nbody=40", "usequad=xyz"]);
    assert_ne!(c, 0);
    assert!(err.contains("not bool"), "{}", err);
}

#[test]
fn subprocess_walktree_overflow() {
    let (c, _, err) = run_args(&["nbody=4"]);
    assert_ne!(c, 0);
    assert!(err.contains("overflow"), "{}", err);
}

#[test]
fn subprocess_theta_zero_oom() {
    let (c, _, err) = run_args(&["nbody=40", "theta=0.0"]);
    assert_ne!(c, 0);
    assert!(err.contains("out of memory"), "{}", err);
}
