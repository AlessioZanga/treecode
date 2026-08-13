use approx::assert_relative_eq;
use treecode::wrapper;

#[test]
fn test_wrapper_rsqr() {
    assert_eq!(wrapper::rsqr(2.0), 4.0);
    assert_eq!(wrapper::rsqr(-3.0), 9.0);
}

#[test]
fn test_wrapper_rqbe() {
    assert_eq!(wrapper::rqbe(2.0), 8.0);
    assert_eq!(wrapper::rqbe(-3.0), -27.0);
}

#[test]
fn test_wrapper_rlog2() {
    assert_relative_eq!(wrapper::rlog2(8.0), 3.0);
}

#[test]
fn test_wrapper_rexp2() {
    assert_relative_eq!(wrapper::rexp2(3.0), 8.0);
}

#[test]
fn test_wrapper_rdex() {
    assert_relative_eq!(wrapper::rdex(2.0), 100.0, max_relative = 0.01);
}

#[test]
fn test_wrapper_fcbrt() {
    assert_relative_eq!(wrapper::fcbrt(8.0), 2.0);
}

#[test]
fn test_wrapper_xrandom() {
    let val = wrapper::xrandom(0.0, 1.0);
    assert!((0.0..=1.0).contains(&val));
}

#[test]
fn test_wrapper_grandom() {
    let val = wrapper::grandom(0.0, 1.0);
    assert!(val.is_finite());
}

#[test]
fn test_wrapper_pickshell() {
    let mut vec = [0.0f32; 3];
    wrapper::pickshell(&mut vec, 3, 1.0);
    let r: f32 = vec.iter().map(|x| x * x).sum();
    assert_relative_eq!(r.sqrt(), 1.0, max_relative = 1e-4);
}

#[test]
fn test_wrapper_pickball() {
    let mut vec = [0.0f32; 3];
    wrapper::pickball(&mut vec, 3, 1.0);
    let r: f32 = vec.iter().map(|x| x * x).sum();
    assert!(r <= 1.001);
}

#[test]
fn test_wrapper_pickbox() {
    let mut vec = [0.0f32; 3];
    wrapper::pickbox(&mut vec, 3, 1.0);
    for v in &vec {
        assert!(v.abs() <= 1.0);
    }
}

#[test]
fn test_wrapper_scanopt() {
    assert!(wrapper::scanopt("a,b,c", "b"));
    assert!(!wrapper::scanopt("a,b,c", "d"));
}

#[test]
fn test_wrapper_allocate() {
    let ptr = wrapper::allocate(100);
    assert!(!ptr.is_null());
}

#[test]
fn test_wrapper_cputime() {
    let t = wrapper::cputime();
    assert!(t >= 0.0);
}

#[test]
fn test_wrapper_initparam_and_getparam() {
    let mut argv = ["test", "nbody=50"];
    let mut defv = [";test", "nbody=10", "tstop=1.0"];
    wrapper::initparam(&mut argv, &mut defv);
    assert_eq!(wrapper::getparam("nbody"), "50");
    assert_eq!(wrapper::getparam("tstop"), "1.0");
}

#[test]
fn test_wrapper_getiparam() {
    let mut argv = ["test", "nbody=256"];
    let mut defv = [";test", "nbody=10"];
    wrapper::initparam(&mut argv, &mut defv);
    assert_eq!(wrapper::getiparam("nbody"), 256);
}

#[test]
fn test_wrapper_getdparam() {
    let mut argv = ["test", "eps=0.05"];
    let mut defv = [";test", "eps=0.025"];
    wrapper::initparam(&mut argv, &mut defv);
    let val = wrapper::getdparam("eps");
    assert_relative_eq!(val, 0.05);
}

#[test]
fn test_wrapper_getbparam() {
    let mut argv = ["test", "usequad=true"];
    let mut defv = [";test", "usequad=false"];
    wrapper::initparam(&mut argv, &mut defv);
    assert!(wrapper::getbparam("usequad"));
}

#[test]
fn test_wrapper_param_roundtrip() {
    let mut argv = ["test", "name=value", "count=42", "ratio=2.71", "flag=true"];
    let mut defv = [
        ";test",
        "name=default",
        "count=0",
        "ratio=0.0",
        "flag=false",
    ];
    wrapper::initparam(&mut argv, &mut defv);
    assert_eq!(wrapper::getparam("name"), "value");
    assert_eq!(wrapper::getiparam("count"), 42);
    assert_relative_eq!(wrapper::getdparam("ratio"), 2.71);
    assert!(wrapper::getbparam("flag"));
}

#[test]
fn test_wrapper_allocate_sizes() {
    for size in [1, 10, 100, 1000] {
        let ptr = wrapper::allocate(size);
        assert!(!ptr.is_null());
    }
}

#[test]
fn test_wrapper_math_chain() {
    let x = 2.0f32;
    let x2 = wrapper::rsqr(x);
    let x3 = wrapper::rqbe(x);
    assert_relative_eq!(x2, 4.0);
    assert_relative_eq!(x3, 8.0);
    let log2x = wrapper::rlog2(x);
    let exp2log = wrapper::rexp2(log2x);
    assert_relative_eq!(exp2log, x, max_relative = 1e-4);
}

#[test]
fn test_wrapper_cbrt_values() {
    let values = [(0.0, 0.0), (1.0, 1.0), (8.0, 2.0), (27.0, 3.0)];
    for (input, expected) in values {
        let result = wrapper::fcbrt(input);
        assert_relative_eq!(result, expected);
    }
}

#[test]
fn test_wrapper_random_distributions() {
    let mut sum = 0.0;
    for _ in 0..1000 {
        sum += wrapper::xrandom(0.0, 1.0);
    }
    let mean = sum / 1000.0;
    assert_relative_eq!(mean, 0.5, max_relative = 0.2);
}

#[test]
fn test_wrapper_normal_distribution() {
    let mut sum = 0.0;
    for _ in 0..1000 {
        sum += wrapper::grandom(5.0, 1.0);
    }
    let mean = sum / 1000.0;
    assert_relative_eq!(mean, 5.0, max_relative = 0.1);
}
