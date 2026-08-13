use approx::assert_relative_eq;
use treecode::wrapper;

#[test]
fn test_rsqr() {
    assert_eq!(wrapper::rsqr(0.0), 0.0);
    assert_eq!(wrapper::rsqr(1.0), 1.0);
    assert_eq!(wrapper::rsqr(2.0), 4.0);
    assert_eq!(wrapper::rsqr(-3.0), 9.0);
    assert_relative_eq!(wrapper::rsqr(0.5), 0.25);
}

#[test]
fn test_rqbe() {
    assert_eq!(wrapper::rqbe(0.0), 0.0);
    assert_eq!(wrapper::rqbe(1.0), 1.0);
    assert_eq!(wrapper::rqbe(2.0), 8.0);
    assert_eq!(wrapper::rqbe(-3.0), -27.0);
}

#[test]
fn test_rlog2() {
    assert_relative_eq!(wrapper::rlog2(1.0), 0.0);
    assert_relative_eq!(wrapper::rlog2(2.0), 1.0);
    assert_relative_eq!(wrapper::rlog2(4.0), 2.0);
    assert_relative_eq!(wrapper::rlog2(8.0), 3.0);
}

#[test]
fn test_rexp2() {
    assert_relative_eq!(wrapper::rexp2(0.0), 1.0);
    assert_relative_eq!(wrapper::rexp2(1.0), 2.0);
    assert_relative_eq!(wrapper::rexp2(2.0), 4.0);
    assert_relative_eq!(wrapper::rexp2(3.0), 8.0);
}

#[test]
fn test_rdex() {
    assert_relative_eq!(wrapper::rdex(0.0), 1.0);
    assert_relative_eq!(wrapper::rdex(1.0), 10.0);
    assert_relative_eq!(wrapper::rdex(2.0), 100.0);
}

#[test]
fn test_fcbrt() {
    assert_relative_eq!(wrapper::fcbrt(0.0), 0.0);
    assert_relative_eq!(wrapper::fcbrt(1.0), 1.0);
    assert_relative_eq!(wrapper::fcbrt(8.0), 2.0);
    assert_relative_eq!(wrapper::fcbrt(27.0), 3.0);
}

#[test]
fn test_xrandom_range() {
    for _ in 0..1000 {
        let val = wrapper::xrandom(-1.0, 1.0);
        assert!((-1.0..=1.0).contains(&val), "xrandom out of range: {}", val);
    }
}

#[test]
fn test_xrandom_distribution() {
    let n = 10000;
    let mut sum = 0.0;
    for _ in 0..n {
        sum += wrapper::xrandom(0.0, 1.0);
    }
    let mean = sum / n as f64;
    assert_relative_eq!(mean, 0.5, max_relative = 0.1);
}

#[test]
fn test_grandom_distribution() {
    let n = 10000;
    let mut sum = 0.0;
    for _ in 0..n {
        sum += wrapper::grandom(0.0, 1.0);
    }
    let mean = sum / n as f64;
    assert_relative_eq!(mean, 0.0, epsilon = 0.1);
}

#[test]
fn test_pickshell_in_range() {
    let mut vec = [0.0f32; 3];
    for _ in 0..100 {
        wrapper::pickshell(&mut vec, 3, 1.0);
        let r: f32 = vec.iter().map(|x| x * x).sum();
        assert_relative_eq!(r.sqrt(), 1.0, max_relative = 1e-4);
    }
}

#[test]
fn test_pickball_in_range() {
    let mut vec = [0.0f32; 3];
    for _ in 0..100 {
        wrapper::pickball(&mut vec, 3, 1.0);
        let r: f32 = vec.iter().map(|x| x * x).sum();
        assert!(r <= 1.001, "pickball outside unit ball: {:?}", vec);
    }
}

#[test]
fn test_pickbox_in_range() {
    let mut vec = [0.0f32; 3];
    for _ in 0..100 {
        wrapper::pickbox(&mut vec, 3, 1.0);
        for v in &vec {
            assert!(v.abs() <= 1.0, "pickbox outside cube: {:?}", vec);
        }
    }
}

#[test]
fn test_scanopt_found() {
    assert!(wrapper::scanopt("a,b,c", "a"));
    assert!(wrapper::scanopt("a,b,c", "b"));
    assert!(wrapper::scanopt("a,b,c", "c"));
}

#[test]
fn test_scanopt_not_found() {
    assert!(!wrapper::scanopt("a,b,c", "d"));
    assert!(!wrapper::scanopt("a,b,c", "ab"));
    assert!(!wrapper::scanopt("a,b,c", ""));
}

#[test]
fn test_scanopt_single() {
    assert!(wrapper::scanopt("only", "only"));
    assert!(!wrapper::scanopt("only", "two"));
}
