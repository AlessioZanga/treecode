use std::io::Write;

pub fn allocate(nb: usize) -> *mut u8 {
    unsafe {
        let ptr = libc::calloc(nb as libc::size_t, 1) as *mut u8;
        if ptr.is_null() {
            eprintln!("allocate: not enough memory ({} bytes)", nb);
            std::process::exit(1);
        }
        ptr
    }
}

pub fn cputime() -> f64 {
    unsafe {
        let mut buffer: libc::tms = std::mem::zeroed();
        if libc::times(&mut buffer) == -1 {
            eprintln!("cputime: times() call failed");
            std::process::exit(1);
        }
        let hz = libc::sysconf(libc::_SC_CLK_TCK) as f64;
        (buffer.tms_utime + buffer.tms_stime) as f64 / (60.0 * hz)
    }
}

pub fn eprintf(fmt: &str) {
    eprint!("{}", fmt);
    let _ = std::io::stderr().flush();
}

pub fn error(fmt: &str) {
    eprint!("{}", fmt);
    let _ = std::io::stderr().flush();
    std::process::exit(1);
}

pub fn scanopt(opt: &str, key: &str) -> bool {
    for word in opt.split(',') {
        if word == key {
            return true;
        }
    }
    false
}
