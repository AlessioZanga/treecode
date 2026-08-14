use approx::assert_relative_eq;
use treecode::mathfns;
use treecode::rng;
use treecode::types::{scanopt, Vector};

#[test]
fn test_rsqr() {
    assert_eq!(mathfns::rsqr(0.0), 0.0);
    assert_eq!(mathfns::rsqr(1.0), 1.0);
    assert_eq!(mathfns::rsqr(2.0), 4.0);
    assert_eq!(mathfns::rsqr(-3.0), 9.0);
    assert_relative_eq!(mathfns::rsqr(0.5), 0.25);
}

#[test]
fn test_rqbe() {
    assert_eq!(mathfns::rqbe(0.0), 0.0);
    assert_eq!(mathfns::rqbe(1.0), 1.0);
    assert_eq!(mathfns::rqbe(2.0), 8.0);
    assert_eq!(mathfns::rqbe(-3.0), -27.0);
}

#[test]
fn test_rlog2() {
    assert_relative_eq!(mathfns::rlog2(1.0), 0.0);
    assert_relative_eq!(mathfns::rlog2(2.0), 1.0);
    assert_relative_eq!(mathfns::rlog2(4.0), 2.0);
    assert_relative_eq!(mathfns::rlog2(8.0), 3.0);
}

#[test]
fn test_rexp2() {
    assert_relative_eq!(mathfns::rexp2(0.0), 1.0);
    assert_relative_eq!(mathfns::rexp2(1.0), 2.0);
    assert_relative_eq!(mathfns::rexp2(2.0), 4.0);
    assert_relative_eq!(mathfns::rexp2(3.0), 8.0);
}

#[test]
fn test_rdex() {
    assert_relative_eq!(mathfns::rdex(0.0), 1.0);
    assert_relative_eq!(mathfns::rdex(1.0), 10.0);
    assert_relative_eq!(mathfns::rdex(2.0), 100.0);
}

#[test]
fn test_fcbrt() {
    assert_relative_eq!(mathfns::fcbrt(0.0), 0.0);
    assert_relative_eq!(mathfns::fcbrt(1.0), 1.0);
    assert_relative_eq!(mathfns::fcbrt(8.0), 2.0);
    assert_relative_eq!(mathfns::fcbrt(27.0), 3.0);
}

#[test]
fn test_xrandom_range() {
    let mut rng = rng::RngState::new(42);
    for _ in 0..1000 {
        let val = rng::xrandom(&mut rng, -1.0, 1.0);
        assert!((-1.0..=1.0).contains(&val), "xrandom out of range: {}", val);
    }
}

#[test]
fn test_xrandom_distribution() {
    let mut rng = rng::RngState::new(42);
    let n = 10000;
    let mut sum = 0.0;
    for _ in 0..n {
        sum += rng::xrandom(&mut rng, 0.0, 1.0);
    }
    let mean = sum / n as f64;
    assert_relative_eq!(mean, 0.5, max_relative = 0.1);
}

#[test]
fn test_grandom_distribution() {
    let mut rng = rng::RngState::new(42);
    let n = 10000;
    let mut sum = 0.0;
    for _ in 0..n {
        sum += rng::grandom(&mut rng, 0.0, 1.0);
    }
    let mean = sum / n as f64;
    assert_relative_eq!(mean, 0.0, epsilon = 0.1);
}

#[test]
fn test_pickshell_in_range() {
    let mut rng = rng::RngState::new(42);
    let mut vec = Vector::zero();
    for _ in 0..100 {
        rng::pickshell(&mut rng, &mut vec, 3, 1.0);
        let r: f32 = vec.iter().map(|x| x * x).sum();
        assert_relative_eq!(r.sqrt(), 1.0, max_relative = 1e-4);
    }
}

#[test]
fn test_pickball_in_range() {
    let mut rng = rng::RngState::new(42);
    let mut vec = Vector::zero();
    for _ in 0..100 {
        rng::pickball(&mut rng, &mut vec, 3, 1.0);
        let r: f32 = vec.iter().map(|x| x * x).sum();
        assert!(r <= 1.001, "pickball outside unit ball: {:?}", vec);
    }
}

#[test]
fn test_pickbox_in_range() {
    let mut rng = rng::RngState::new(42);
    let mut vec = Vector::zero();
    for _ in 0..100 {
        rng::pickbox(&mut rng, &mut vec, 3, 1.0);
        for v in &vec {
            assert!(v.abs() <= 1.0, "pickbox outside cube: {:?}", vec);
        }
    }
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
