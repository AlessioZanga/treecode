use treecode::types;

#[test]
fn test_allocate_non_null() {
    let ptr = types::allocate(100);
    assert!(!ptr.is_null());
}

#[test]
fn test_allocate_zeroed() {
    let ptr = types::allocate(100);
    let slice = unsafe { std::slice::from_raw_parts(ptr, 100) };
    assert!(slice.iter().all(|&b| b == 0));
}

#[test]
fn test_allocate_small() {
    let ptr = types::allocate(1);
    assert!(!ptr.is_null());
}

#[test]
fn test_scanopt_found() {
    assert!(types::scanopt("a,b,c", "a"));
    assert!(types::scanopt("a,b,c", "b"));
    assert!(types::scanopt("a,b,c", "c"));
}

#[test]
fn test_scanopt_not_found() {
    assert!(!types::scanopt("a,b,c", "d"));
    assert!(!types::scanopt("a,b,c", "ab"));
    assert!(!types::scanopt("a,b,c", ""));
}

#[test]
fn test_scanopt_single() {
    assert!(types::scanopt("only", "only"));
    assert!(!types::scanopt("only", "two"));
}

#[test]
fn test_cputime_positive() {
    let t = types::cputime();
    assert!(t >= 0.0);
}

#[test]
fn test_cputime_units() {
    let t1 = types::cputime();
    std::thread::sleep(std::time::Duration::from_millis(10));
    let t2 = types::cputime();
    assert!(t2 >= t1);
}

#[test]
fn test_eprintf() {
    types::eprintf("test message\n");
}

#[test]
fn test_error_function_exists() {
    // error() calls process::exit so we can't test it directly
    // Just verify the function compiles and is callable
    let _ = types::error;
}
