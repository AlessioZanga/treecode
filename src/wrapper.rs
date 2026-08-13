use crate::types;

pub fn pickshell(vec: &mut [f32], ndim: i32, rad: f32) {
    crate::mathfns::pickshell(vec, ndim as usize, rad);
}

pub fn pickball(vec: &mut [f32], ndim: i32, rad: f32) {
    crate::mathfns::pickball(vec, ndim as usize, rad);
}

pub fn pickbox(vec: &mut [f32], ndim: i32, size: f32) {
    crate::mathfns::pickbox(vec, ndim as usize, size);
}

pub fn xrandom(xl: f64, xh: f64) -> f64 {
    crate::mathfns::xrandom(xl, xh)
}

pub fn grandom(mean: f64, sdev: f64) -> f64 {
    crate::mathfns::grandom(mean, sdev)
}

pub fn scanopt(opt: &str, key: &str) -> bool {
    crate::clib::scanopt(opt, key)
}

pub fn maketree(btab: &mut [types::Body], nbody: i32) {
    crate::treeload::maketree(btab, nbody)
}

pub fn gravcalc() {
    crate::treegrav::gravcalc()
}

pub fn inputdata() {
    crate::treeio::inputdata()
}

pub fn startoutput() {
    crate::treeio::startoutput()
}

pub fn forcereport() {
    crate::treeio::forcereport()
}

pub fn output() {
    crate::treeio::output()
}

pub fn savestate(pattern: &str) {
    crate::treeio::savestate(pattern)
}

pub fn restorestate(file: &str) {
    crate::treeio::restorestate(file)
}

pub fn initparam(argv: &mut [&str], defv: &mut [&str]) {
    crate::getparam::initparam(argv, defv);
}

pub fn getparam(name: &str) -> String {
    crate::getparam::getparam(name)
}

pub fn getiparam(name: &str) -> i32 {
    crate::getparam::getiparam(name)
}

pub fn getdparam(name: &str) -> f64 {
    crate::getparam::getdparam(name)
}

pub fn getbparam(name: &str) -> bool {
    crate::getparam::getbparam(name)
}

pub fn allocate(nb: i32) -> *mut std::ffi::c_void {
    crate::clib::allocate(nb as usize) as *mut std::ffi::c_void
}

pub fn cputime() -> f64 {
    crate::clib::cputime()
}

pub fn rsqr(x: f32) -> f32 {
    crate::mathfns::rsqr(x)
}

pub fn rqbe(x: f32) -> f32 {
    crate::mathfns::rqbe(x)
}

pub fn rlog2(x: f32) -> f32 {
    crate::mathfns::rlog2(x)
}

pub fn rexp2(x: f32) -> f32 {
    crate::mathfns::rexp2(x)
}

pub fn rdex(x: f32) -> f32 {
    crate::mathfns::rdex(x)
}

pub fn fcbrt(x: f32) -> f32 {
    crate::mathfns::fcbrt(x)
}
