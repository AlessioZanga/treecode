# AGENTS.md

Guidance for AI agents (and humans) working in this repository.

## What this repo is

A **faithful, 100% safe Rust port** of Joshua E. Barnes' hierarchical N‑body
**treecode** v1.4. The C reference is vendored under `external/treecode/`. The
port must stay **byte‑exact** with that C binary: particle files and
interaction counts must match bit‑for‑bit. See `README.md`.

## The single most important invariant

> **Output must remain byte‑exact vs the C reference.**

Any change that **reorders floating‑point operations** on `f32` will change the
results and fail `tests/test_compare_c_rust.rs`. The float‑order‑sensitive
functions are: `sumnode`, `sumcell`, `gravsum`, `hackcofm`,
`accumulate_moment`, `quadrupole_tensor`, `stepsystem`, `testdata`, and all the
`vecmath`/`rng` math.

**Safe changes** (structure/layout/control‑flow only, float order untouched):
allocator reuse, copying values by value instead of accessor calls, `#[inline]`,
branchless rewrites that are bit‑identical (`abs`, sign tricks), extracting
injectable I/O cores (`*_to` writers).

**Gated changes** (must be verified against the byte‑exact test, or kept behind
a feature flag): any parallelism that reorders per‑body reductions, FMA. (SIMD is
now on by default via the `wide` crate and is itself gated behind the `simd`
feature / `--no-default-features`, so it can still be disabled.)

Never enable FMA: it changes `f32` rounding vs the SSE2 C reference. The build
already sets `-C target-feature=-fma` in `.cargo/config.toml`.

## Hard constraints

- `#![forbid(unsafe_code)]` (in `src/lib.rs`, `src/main.rs`).
- `#![forbid(clippy::unwrap_used)]` and `#![forbid(clippy::expect_used)]`.
  Handle errors via the `Result`/`TreeError` types in `src/error.rs`.
- No `#[allow(...)]` lint suppressions anywhere in the crate. Fix the lint
  instead of allowing it.
- Keep the 1:1 function mapping to the C source so the port stays auditable
  (e.g. `maketree`, `gravcalc`, `hackcofm`, `threadtree`).

## Layout & key modules

- `src/vecmath.rs` — `Vector`/`Matrix` (`f32`, `#[repr(C)]`, 3‑D). Layout‑exact
  with C; never change field ordering or precision (`Real = f32`).
- `src/types.rs` — `Node`, `Body`, `Cell`, `Sorq` (subcells **or** quadrupole),
  `NodeRef` (body/cell index discriminant), `Interact`, `cputime`.
- `src/rng.rs` — glibc `random()` reimplementation; the `unit()` stream MUST
  stay byte‑identical to C (do not swap PRNGs).
- `src/getparam.rs` — parameter table + typed getters (`getparam`, `getiparam`,
  `getdparam`, `getbparam`, `getparamstat`).
- `src/treeload.rs` — tree build: `maketree`, `loadbody`, `hackcofm`,
  `hackquad`, `threadtree`.
- `src/treegrav.rs` — force walk: `gravcalc`, `walktree`/`walksub`/`gravsum`,
  parallel fan‑out (`walk_parallel`/`run_child_walk`). Walk functions take a
  single `WalkContext<'a>` to stay under the clippy `too_many_arguments` limit.
- `src/treeio.rs` — `inputdata`, `output`/`outputdata`, `diagnostics`,
  `savestate`/`restorestate`. Prefer extending the injectable `*_to` writers
  over hardcoding `Write`/`Read` sinks.
- `src/treecode.rs` — `Tree` (all formerly‑global state), `stepsystem`,
  `Simulation` handle, `run`.
- `tests/test_compare_c_rust.rs` — **the byte‑exact gate** (runs the C binary +
  Rust port, compares particle files and normalized stdout).
- `benches/throughput.rs` — Rust vs C ratio (`harness = false`).

## Commands (run before committing)

```sh
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo test
cargo test --release --test test_compare_c_rust   # must stay byte-exact
cargo bench                                        # rust/c ratio in 0.1..=2.0
```

Note: `.cargo/config.toml` forces `RUST_TEST_THREADS=1` and
`RUST_TEST_NOCAPTURE=1` (the byte‑exact test captures stdout via `dup2`, so
tests must run single‑threaded). `pre-commit` (see `.pre-commit-config.yaml`)
runs fmt + clippy + `cargo test -- --test-threads=1`.

## Editing patterns

- **Adding a parameter:** extend `default_defv()` in `src/treecode.rs` and add
  parsing in `startrun`/`restorestate`; mirror C semantics exactly.
- **New I/O:** add an injectable `*_to(&mut impl Write/Read)` core, keep the
  file‑opening wrapper thin, and preserve the `f32` `to_ne_bytes` /
  `%14.7E` formatting in `treeio.rs` (`fmt_e14`) — that formatting is required
  for byte‑exact output.
- **Parallelism:** only the tree‑walk **fan‑out** form (private scratch per
  subtree, identical per‑body reduction order) is allowed. Never do per‑body
  `rayon` parallelism — it reorders interaction lists and breaks byte‑exactness.
- **SIMD:** **on by default** (the `simd` cargo feature is in `default`); built
  on the `wide` crate so it runs on **stable Rust** and cross‑compiles for
  Windows (verified for `x86_64-pc-windows-gnu` and `-msvc`). `wide::f32x8` uses
  AVX2 when present and transparently falls back to two `f32x4` on SSE2, with
  identical per‑lane rounding either way. The kernels in `src/treegrav.rs`
  (`sumnode_simd` / `sumcell_simd`) vectorize the *independent per‑interaction*
  map (`sqrt`, divide, `dr*mr3i`, the quadrupole matvec) and fold back into the
  running `f32` accumulators **one lane at a time in list order** — no FP
  reassociation, so they stay bit‑identical to the scalar path and the C
  reference. FMA is still avoided (the kernels use separate mul/add to match the
  `-fma`-disabled SSE2 reference). The same pattern is applied in `stepsystem`
  (a pure per‑body map, so lane order is irrelevant) and in the tree‑build
  moment reductions `accumulate_subnodes`/`hackcofm` and `accumulate_moment`/
  `hackquad` (NSUB=8 descendants map perfectly onto `f32x8`; `None` subnodes
  contribute exactly 0). Disable with `--no-default-features` to fall back to the
  scalar path. Verify any SIMD change with
  `cargo test --release --test test_compare_c_rust` (the default build already
  enables SIMD).
- **Tests:** add unit tests next to the code (`#[cfg(test)]` modules) and
  integration tests under `tests/`. The `test_compare_c_rust` gate must remain
  green; do not weaken its byte‑for‑byte `assert_eq!` on particle files.

## Reference docs

- `external/treecode/README.md` and `external/treecode/treeguide.html` — the
  original C documentation and algorithm guide.
