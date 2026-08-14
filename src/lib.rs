#![forbid(unsafe_code)]
#![forbid(clippy::unwrap_used)]
#![forbid(clippy::expect_used)]

//! Faithful, safe Rust port of Joshua E. Barnes' hierarchical N-body treecode
//! (`treecode`), refactored from the original C into idiomatic Rust while
//! preserving a 1:1 mapping of function names so the implementation stays
//! auditable against the reference.
//!
//! Two API layers are exposed:
//!
//! * The **high-level** [`Simulation`] handle: parse parameters and run the
//!   whole integration, or drive it step-by-step.
//! * The **low-level** 1:1 functions (`maketree`, `gravcalc`, `stepsystem`,
//!   `outputdata`, `savestate`, …), kept public on their respective modules so
//!   the C mapping remains exact and inspectable.
//!
//! Parameter parsing lives in [`getparam`], the math helpers in [`mathfns`] and
//! [`vecmath`], the PRNG in [`rng`], and all simulation state in [`Tree`].

pub mod error;
pub mod getparam;
pub mod mathfns;
pub mod rng;
pub mod treecode;
pub mod treegrav;
pub mod treeio;
pub mod treeload;
pub mod types;
pub mod vecmath;

pub use treecode::Simulation;
