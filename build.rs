fn main() {
    cc::Build::new()
        .file("external/treecode/clib.c")
        .file("external/treecode/getparam.c")
        .file("external/treecode/mathfns.c")
        .file("external/treecode/treeload.c")
        .file("external/treecode/treegrav.c")
        .file("external/treecode/treeio.c")
        .flag("-DLINUX")
        .flag("-DSINGLEPREC")
        .flag("-DTHREEDIM")
        .flag("-O3")
        .compile("treecode_c");

    println!("cargo:rerun-if-changed=external/treecode/");
    println!("cargo:rerun-if-changed=build.rs");
}
