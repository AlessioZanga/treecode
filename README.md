# treecode

A faithful, **100% safe** Rust port of Joshua E. Barnes' hierarchical N‑body
**treecode** (Version 1.4). The force algorithm is a line‑for‑line mapping of
the original C `treecode` (`external/treecode/`), refactored into idiomatic,
auditable Rust while preserving a **byte‑exact** match with the reference
binary.

- **Hierarchical gravity.** O(N log N) tree‑based approximation of the
  gravitational force (and potential) on every body, using the Barnes & Hut
  opening criterion and optional quadrupole moments for higher accuracy.
- **Byte‑exact output.** Particle files and interaction counts are
  bit‑for‑bit identical to the C reference. The only exempt columns are
  run‑to‑run CPU timings and the parallel force walk's peak `actmax`.
- **Safe & clean.** `#![forbid(unsafe_code)]`, no `unwrap`/`expect`
  (`clippy::unwrap_used`/`clippy::expect_used` are forbidden), and a
  `clippy --all-targets -D warnings` / `rustfmt` clean codebase.
- **Parallel fan‑out.** The force walk fans the tree's root children across
  OS threads (`std::thread::scope`); each subtree uses private scratch so the
  per‑body float reduction order is unchanged → physics stays byte‑exact while
  Rust runs ~3× faster than C at scale.

The original C is vendored under [`external/treecode/`](external/treecode) so
the port can be diffed against it and validated by the byte‑exact test suite.

## Building

```sh
cargo build --release
```

No external runtime dependencies beyond `thiserror`. The build is tuned via
`.cargo/config.toml` with `-C target-cpu=native -C target-feature=-fma`; FMA is
disabled because enabling it changes `f32` rounding and would break byte‑exact
parity with the SSE2 C reference.

## Running

The binary mirrors the C command‑line interface (name=value parameters):

```sh
cargo run --release -- nbody=4096 tstop=10.0 dtout=0.025 seed=123
# or, to write a particle snapshot:
cargo run --release -- nbody=1000 tstop=1.0 out=parts.txt
```

Common parameters: `in=`, `out=`, `save=`, `restore=`, `dtime=`, `eps=`,
`theta=`, `usequad=`, `options=`, `tstop=`, `dtout=`, `nbody=`, `seed=`.
Run with `-help` / `-clue` to print the parameter table. When no `in=` file is
given, an initial Plummer‑model test distribution is generated internally.

### Library usage

```rust
use treecode::Simulation;

let mut sim = Simulation::new(["nbody=16", "tstop=0.01", "dtout=0.005"])?;
sim.run()?;
assert!(sim.state().nstep > 0);
# Ok::<(), treecode::error::TreeError>(())
```

The high‑level [`Simulation`](src/treecode.rs) handle parses parameters and
runs the whole integration, or lets you drive it step‑by‑step (`step` /
`output`). The low‑level 1:1 functions (`maketree`, `gravcalc`, `stepsystem`,
`outputdata`, `savestate`, …) are also public on their modules for inspection.

## Testing & verification

The central guarantee is **byte‑exactness vs the C reference**:

```sh
cargo test                                  # all unit/integration tests
cargo test --release --test test_compare_c_rust   # C vs Rust byte-exact gate
cargo bench                                  # Rust/C throughput ratio
cargo clippy --all-targets -- -D warnings    # lint gate
cargo fmt --check                            # format gate
```

`tests/test_compare_c_rust.rs` runs the vendored C binary and the Rust port
over identical parameters and asserts their particle files match
byte‑for‑byte and their stdout diagnostics match (with the volatile CPU‑timing
and `actmax` columns normalized out).

Optional quality gates used historically: `cargo llvm-cov` (line coverage
≥ 90%, measured 94.25%) and `pre-commit run --all-files` (fmt + clippy + test).

## Project layout

| Path | Contents |
|------|----------|
| `src/vecmath.rs` | `#[repr(C)]` `Vector`/`Matrix` (`f32`, 3‑D), mirroring `vectmath.h`. |
| `src/types.rs` | Core types: `Node`, `Body`, `Cell`, `Sorq`, `NodeRef`, CPU‑time reading. |
| `src/mathfns.rs` | Scalar math helpers (`rsqr`, `rlog2`, …). |
| `src/rng.rs` | Pure‑Rust reimplementation of glibc `random()` (byte‑identical stream). |
| `src/getparam.rs` | Parameter parsing (`initparam` + typed getters), replacing C globals. |
| `src/treeload.rs` | Tree construction: `maketree`, `loadbody`, `hackcofm`, `hackquad`, `threadtree`. |
| `src/treegrav.rs` | Force calculation: `gravcalc`, `walktree`/`walksub`, `gravsum`, parallel fan‑out. |
| `src/treeio.rs` | Input/output, diagnostics, save/restore state (byte‑exact serialization). |
| `src/treecode.rs` | `Tree` state, `stepsystem`, `Simulation` handle, `run`. |
| `src/error.rs` | `thiserror`‑based `TreeError`. |
| `tests/` | `test_compare_c_rust` (byte‑exact gate) and unit/integration suites. |
| `benches/throughput.rs` | Rust vs C throughput benchmark (`harness = false`). |
| `external/treecode/` | Vendored C reference (v1.4) used for the byte‑exact comparison. |

## Performance notes

- Interaction buffer is a lightweight `Interact { mass, pos, quad }` (~56 B)
  instead of copying whole `Cell`s (~192 B) into the walk — ~3.4× smaller.
- Scratch buffers are retained across timesteps; the traversal copies `Node`s
  by value to avoid repeated enum dispatch; hot leaves are `#[inline]`.
- `gravcalc` fans the root's `NSUB` children across threads via
  `std::thread::scope` with private scratch → byte‑exact physics and a
  Rust/C throughput ratio of ≈0.27 (Rust ~3.7× faster) at nbody=1500.
