use approx::assert_relative_eq;
use treecode::getparam;

#[test]
fn test_initparam_and_getparam() {
    let argv = ["test_program", "nbody=100", "tstop=0.5"];
    let defv = [
        ";Test program",
        "nbody=10",
        "tstop=1.0",
        "eps=0.025",
        "VERSION=1.0",
    ];
    getparam::initparam(&argv, &defv).unwrap();
    assert_eq!(getparam::getparam("argv0").unwrap(), "test_program");
    assert_eq!(getparam::getparam("nbody").unwrap(), "100");
    assert_eq!(getparam::getparam("tstop").unwrap(), "0.5");
    assert_eq!(getparam::getparam("eps").unwrap(), "0.025");
}

#[test]
fn test_getiparam() {
    let argv = ["test", "nbody=256"];
    let defv = [";test", "nbody=10"];
    getparam::initparam(&argv, &defv).unwrap();
    assert_eq!(getparam::getiparam("nbody").unwrap(), 256);
}

#[test]
fn test_getiparam_suffix_k() {
    let argv = ["test", "nbody=4k"];
    let defv = [";test", "nbody=10"];
    getparam::initparam(&argv, &defv).unwrap();
    assert_eq!(getparam::getiparam("nbody").unwrap(), 4096);
}

#[test]
fn test_getiparam_suffix_m() {
    let argv = ["test", "nbody=2m"];
    let defv = [";test", "nbody=10"];
    getparam::initparam(&argv, &defv).unwrap();
    assert_eq!(getparam::getiparam("nbody").unwrap(), 2097152);
}

#[test]
fn test_getdparam() {
    let argv = ["test", "eps=0.05"];
    let defv = [";test", "eps=0.025"];
    getparam::initparam(&argv, &defv).unwrap();
    let val = getparam::getdparam("eps").unwrap();
    assert_relative_eq!(val, 0.05);
}

#[test]
fn test_getdparam_fraction() {
    let argv = ["test", "dtime=1/32"];
    let defv = [";test", "dtime=0.03125"];
    getparam::initparam(&argv, &defv).unwrap();
    let val = getparam::getdparam("dtime").unwrap();
    assert_relative_eq!(val, 0.03125);
}

#[test]
fn test_getbparam_true() {
    let argv = ["test", "usequad=true"];
    let defv = [";test", "usequad=false"];
    getparam::initparam(&argv, &defv).unwrap();
    assert!(getparam::getbparam("usequad").unwrap());
}

#[test]
fn test_getbparam_false() {
    let argv = ["test", "usequad=false"];
    let defv = [";test", "usequad=true"];
    getparam::initparam(&argv, &defv).unwrap();
    assert!(!getparam::getbparam("usequad").unwrap());
}

#[test]
fn test_getbparam_variants() {
    for &val in &["t", "T", "y", "Y", "1"] {
        let argv = ["test", &format!("flag={}", val)];
        let defv = [";test", "flag=false"];
        getparam::initparam(&argv, &defv).unwrap();
        assert!(
            getparam::getbparam("flag").unwrap(),
            "Expected true for '{}'",
            val
        );
    }
    for &val in &["f", "F", "n", "N", "0"] {
        let argv = ["test", &format!("flag={}", val)];
        let defv = [";test", "flag=true"];
        getparam::initparam(&argv, &defv).unwrap();
        assert!(
            !getparam::getbparam("flag").unwrap(),
            "Expected false for '{}'",
            val
        );
    }
}

#[test]
fn test_getparamstat_default() {
    let argv = ["test", "nbody=100"];
    let defv = [";test", "nbody=10", "extra=5"];
    getparam::initparam(&argv, &defv).unwrap();
    let stat = getparam::getparamstat("extra");
    assert!(stat & 0o1 != 0);
}

#[test]
fn test_getparamstat_arg() {
    let argv = ["test", "nbody=100"];
    let defv = [";test", "nbody=10"];
    getparam::initparam(&argv, &defv).unwrap();
    let stat = getparam::getparamstat("nbody");
    assert!(stat & 0o4 != 0);
}

#[test]
fn test_getparamstat_unknown() {
    let argv = ["test"];
    let defv = [";test", "nbody=10"];
    getparam::initparam(&argv, &defv).unwrap();
    assert_eq!(getparam::getparamstat("nonexistent"), 0);
}

#[test]
fn test_positional_args() {
    let argv = ["test", "100", "0.5"];
    let defv = [";test", "nbody=10", "tstop=1.0"];
    getparam::initparam(&argv, &defv).unwrap();
    assert_eq!(getparam::getparam("nbody").unwrap(), "100");
    assert_eq!(getparam::getparam("tstop").unwrap(), "0.5");
}

#[test]
fn test_default_values() {
    let argv = ["test"];
    let defv = [";test", "nbody=42", "tstop=2.71"];
    getparam::initparam(&argv, &defv).unwrap();
    assert_eq!(getparam::getparam("nbody").unwrap(), "42");
    assert_eq!(getparam::getparam("tstop").unwrap(), "2.71");
}

#[test]
fn test_comments_in_defaults() {
    let argv = ["test"];
    let defv = [";main comment", "nbody=10", ";another comment", "tstop=1.0"];
    getparam::initparam(&argv, &defv).unwrap();
    assert_eq!(getparam::getparam("nbody").unwrap(), "10");
    assert_eq!(getparam::getparam("tstop").unwrap(), "1.0");
}

#[test]
fn test_getiparam_no_suffix() {
    let argv = ["test", "val=42"];
    let defv = [";test", "val=0"];
    getparam::initparam(&argv, &defv).unwrap();
    assert_eq!(getparam::getiparam("val").unwrap(), 42);
}

#[test]
fn test_getiparam_suffix_uppercase_k() {
    let argv = ["test", "val=2K"];
    let defv = [";test", "val=0"];
    getparam::initparam(&argv, &defv).unwrap();
    assert_eq!(getparam::getiparam("val").unwrap(), 2048);
}

#[test]
fn test_getiparam_suffix_uppercase_m() {
    let argv = ["test", "val=3M"];
    let defv = [";test", "val=0"];
    getparam::initparam(&argv, &defv).unwrap();
    assert_eq!(getparam::getiparam("val").unwrap(), 3145728);
}

#[test]
fn test_getiparam_unknown_suffix() {
    let argv = ["test", "val=99x"];
    let defv = [";test", "val=0"];
    getparam::initparam(&argv, &defv).unwrap();
    assert_eq!(getparam::getiparam("val").unwrap(), 99);
}

#[test]
fn test_getdparam_fraction_complex() {
    let argv = ["test", "val=3/4"];
    let defv = [";test", "val=0"];
    getparam::initparam(&argv, &defv).unwrap();
    let val = getparam::getdparam("val").unwrap();
    assert_relative_eq!(val, 0.75);
}

#[test]
fn test_getdparam_integer() {
    let argv = ["test", "val=42"];
    let defv = [";test", "val=0"];
    getparam::initparam(&argv, &defv).unwrap();
    let val = getparam::getdparam("val").unwrap();
    assert_relative_eq!(val, 42.0);
}

#[test]
fn test_getdparam_float() {
    let argv = ["test", "val=2.71"];
    let defv = [";test", "val=0"];
    getparam::initparam(&argv, &defv).unwrap();
    let val = getparam::getdparam("val").unwrap();
    assert_relative_eq!(val, 2.71);
}

#[test]
fn test_getparamstat_defparam() {
    let argv = ["test"];
    let defv = [";test", "nbody=10"];
    getparam::initparam(&argv, &defv).unwrap();
    let stat = getparam::getparamstat("nbody");
    assert!(stat & 0o1 != 0);
}

#[test]
fn test_multiple_params() {
    let argv = ["test", "a=1", "b=2", "c=3", "d=4", "e=5"];
    let defv = [";test", "a=0", "b=0", "c=0", "d=0", "e=0"];
    getparam::initparam(&argv, &defv).unwrap();
    assert_eq!(getparam::getiparam("a").unwrap(), 1);
    assert_eq!(getparam::getiparam("b").unwrap(), 2);
    assert_eq!(getparam::getiparam("c").unwrap(), 3);
    assert_eq!(getparam::getiparam("d").unwrap(), 4);
    assert_eq!(getparam::getiparam("e").unwrap(), 5);
}

#[test]
fn test_named_args_override_positional() {
    let argv = ["test", "100", "tstop=0.5"];
    let defv = [";test", "nbody=10", "tstop=1.0"];
    getparam::initparam(&argv, &defv).unwrap();
    assert_eq!(getparam::getparam("nbody").unwrap(), "100");
    assert_eq!(getparam::getparam("tstop").unwrap(), "0.5");
}

#[test]
fn test_all_default_values_preserved() {
    let argv = ["test"];
    let defv = [";test", "a=x", "b=y", "c=z"];
    getparam::initparam(&argv, &defv).unwrap();
    assert_eq!(getparam::getparam("a").unwrap(), "x");
    assert_eq!(getparam::getparam("b").unwrap(), "y");
    assert_eq!(getparam::getparam("c").unwrap(), "z");
}

#[test]
fn test_getparamstat_multiple() {
    let argv = ["test", "a=1"];
    let defv = [";test", "a=0", "b=0"];
    getparam::initparam(&argv, &defv).unwrap();
    let stat_a = getparam::getparamstat("a");
    let stat_b = getparam::getparamstat("b");
    assert!(stat_a & 0o4 != 0);
    assert!(stat_b & 0o1 != 0);
}
