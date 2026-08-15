// Cross-platform stdout handling for the in-process integration tests.
//
// On Unix the library prints diagnostics via `print!`; the tests redirect
// file-descriptor 1 to a throwaway sink so the test output stays clean. On
// Windows there is no `dup`/`dup2`, so we simply let the output reach the
// console (the capture is only used to silence output, never to read it).
//
// The byte-exact comparison test still requires the C reference binary and
// stays Unix-only; `test_inprocess.rs` instead spawns the built binary via
// `Command`, which is cross-platform.

#[cfg(unix)]
mod imp {
    use std::{
        os::unix::io::{IntoRawFd, RawFd},
        path::Path,
    };

    const STDOUT_FD: RawFd = 1;

    pub struct StdoutCapture {
        saved: RawFd,
    }

    impl StdoutCapture {
        pub fn new(dir: &Path) -> Self {
            let saved = unsafe { libc::dup(STDOUT_FD) };
            assert!(saved >= 0, "dup failed");
            let path = dir.join("stdout.txt");
            let file = std::fs::File::create(&path).unwrap();
            let fd = file.into_raw_fd();
            unsafe {
                libc::dup2(fd, STDOUT_FD);
                libc::close(fd);
            }
            StdoutCapture { saved }
        }
    }

    impl Drop for StdoutCapture {
        fn drop(&mut self) {
            unsafe {
                libc::dup2(self.saved, STDOUT_FD);
                libc::close(self.saved);
            }
        }
    }
}

#[cfg(windows)]
mod imp {
    use std::path::Path;

    pub struct StdoutCapture;

    impl StdoutCapture {
        pub fn new(_dir: &Path) -> Self {
            StdoutCapture
        }
    }
}

pub use imp::StdoutCapture;
