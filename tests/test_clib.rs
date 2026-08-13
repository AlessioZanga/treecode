use treecode::clib;

#[test]
fn test_allocate_non_null() {
    let ptr = clib::allocate(100);
    assert!(!ptr.is_null());
}

#[test]
fn test_allocate_zeroed() {
    let ptr = clib::allocate(100);
    let slice = unsafe { std::slice::from_raw_parts(ptr, 100) };
    assert!(slice.iter().all(|&b| b == 0));
}

#[test]
fn test_allocate_small() {
    let ptr = clib::allocate(1);
    assert!(!ptr.is_null());
}

#[test]
fn test_scanopt_found() {
    assert!(clib::scanopt("a,b,c", "a"));
    assert!(clib::scanopt("a,b,c", "b"));
    assert!(clib::scanopt("a,b,c", "c"));
}

#[test]
fn test_scanopt_not_found() {
    assert!(!clib::scanopt("a,b,c", "d"));
    assert!(!clib::scanopt("a,b,c", "ab"));
    assert!(!clib::scanopt("a,b,c", ""));
}

#[test]
fn test_scanopt_single() {
    assert!(clib::scanopt("only", "only"));
    assert!(!clib::scanopt("only", "two"));
}

#[test]
fn test_cputime_positive() {
    let t = clib::cputime();
    assert!(t >= 0.0);
}

#[test]
fn test_cputime_units() {
    let t1 = clib::cputime();
    std::thread::sleep(std::time::Duration::from_millis(10));
    let t2 = clib::cputime();
    assert!(t2 >= t1);
}

#[test]
fn test_eprintf() {
    clib::eprintf("test message\n");
}

#[test]
fn test_error_function_exists() {
    // error() calls process::exit so we can't test it directly
    // Just verify the function compiles and is callable
    let _ = clib::error;
}
