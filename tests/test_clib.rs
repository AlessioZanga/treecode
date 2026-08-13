use treecode::error::TreeError;
use treecode::types::{allocate, cputime, eprintf, scanopt};

#[test]
fn test_allocate_non_null() {
    let ptr = allocate(100).unwrap();
    assert!(!ptr.is_null());
}

#[test]
fn test_allocate_zeroed() {
    let ptr = allocate(100).unwrap();
    let slice = unsafe { std::slice::from_raw_parts(ptr, 100) };
    assert!(slice.iter().all(|&b| b == 0));
}

#[test]
fn test_allocate_small() {
    let ptr = allocate(1).unwrap();
    assert!(!ptr.is_null());
}

#[test]
fn test_scanopt_found() {
    assert!(scanopt("a,b,c", "a"));
    assert!(scanopt("a,b,c", "b"));
    assert!(scanopt("a,b,c", "c"));
}

#[test]
fn test_scanopt_not_found() {
    assert!(!scanopt("a,b,c", "d"));
    assert!(!scanopt("a,b,c", "ab"));
    assert!(!scanopt("a,b,c", ""));
}

#[test]
fn test_scanopt_single() {
    assert!(scanopt("only", "only"));
    assert!(!scanopt("only", "two"));
}

#[test]
fn test_cputime_positive() {
    let t = cputime().unwrap();
    assert!(t >= 0.0);
}

#[test]
fn test_cputime_units() {
    let t1 = cputime().unwrap();
    std::thread::sleep(std::time::Duration::from_millis(10));
    let t2 = cputime().unwrap();
    assert!(t2 >= t1);
}

#[test]
fn test_eprintf() {
    eprintf("test message\n");
}

#[test]
fn test_tree_error_display() {
    let e = TreeError::AbsurdNbody(0);
    assert!(format!("{}", e).contains("nbody"));

    let e = TreeError::OutOfMemory(1024);
    assert!(format!("{}", e).contains("1024"));

    let e = TreeError::IncompatibleOptions;
    assert!(format!("{}", e).contains("incompatible"));
}
