# Contributing to archi-kernel

Issues and pull requests are welcome — in Japanese or English.

## Before you start

Read [DESIGN.md](DESIGN.md) (Japanese). It is the project's source of truth: every
non-obvious decision in the code traces back to a section there, and changes that
contradict it need a design argument, not just working code. The short version of
the philosophy:

- **Closed form or explicit error.** No numerical marching, no silent
  approximation. If an operation cannot be answered exactly within the
  plane + cylinder vocabulary, it must return a machine-readable error.
- **Failures are member-local and detectable.** A boolean that cannot produce a
  watertight, Euler-consistent result returns `EvalError` — it never hands back
  corrupted geometry.
- **Topology never touches coordinates.** `src/topo/` references geometry only
  through ID handles; CI enforces this with an import check.
- **One absolute tolerance.** All predicates go through `Tol` / the
  `predicate/` facade. No bare float comparisons, no ad-hoc epsilons.

## Gates (all must pass)

```bash
cargo test --all-features
cargo clippy --all-targets --all-features -- -D warnings
cargo fmt --all -- --check
cargo check --no-default-features
```

## Conventions

- Public types carry rustdoc; public enums are `#[non_exhaustive]`; public
  constructors return `Result` (the library never panics on user input).
- Numeric literals in tests are `f64`-annotated with explicit tolerances, and
  expected values are hand-computed analytic references where possible.
- Degenerate configurations (coincident faces, tangencies, edge-flush openings)
  are first-class test cases, not afterthoughts — if you fix a bug, land its
  minimal reproduction as a regression test.
- SI units only (metres, radians). Unit conversion belongs to calling adapters.

## Viewer

The Three.js viewer (`viewer/`, wasm bindings in `wasm/`) doubles as an
integration test bed — its demo building has caught real kernel bugs. If your
change affects booleans, sections or tessellation, load it and look:

```bash
wasm-pack build wasm --target web --out-dir ../viewer/pkg --release
cd viewer && python3 -m http.server 8741
```
