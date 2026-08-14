use treecode::error::TreeError;
use treecode::types::{cputime, eprintf, scanopt, Cell, Vector};

#[test]
fn test_arena_cells_default_zeroed() {
    // The tree's cell arena is a `Vec<Cell>`; `makecell` relies on
    // `Cell::default()` zeroing every field (what `allocate`'s calloc
    // used to guarantee).
    let arena: Vec<Cell> = vec![Cell::default(); 128];
    assert!(arena.iter().all(|c| {
        c.cellnode.mass == 0.0
            && c.cellnode.pos == Vector::zero()
            && c.rcrit2 == 0.0
            && c.more.is_none()
    }));
}

#[test]
fn test_arena_vec_growth() {
    // Cells live in a growable `Vec` indexed by `CellId`.
    let mut arena: Vec<Cell> = Vec::new();
    arena.push(Cell::default());
    let first = arena.len() - 1;
    arena.push(Cell::default());
    let second = arena.len() - 1;
    assert_ne!(first, second);
    assert_eq!(arena.len(), 2);
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
