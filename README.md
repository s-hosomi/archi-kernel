# archi-kernel

A domain-specific B-rep geometry kernel for building simulation, written in Rust.

Instead of competing with general-purpose kernels (Parasolid, ACIS, Open CASCADE), this kernel restricts its surface vocabulary to what actually occurs in building structures — **planes and circular cylinders** — so that every surface-surface intersection has a closed-form solution. No NURBS, no numerical curve marching, and a single absolute tolerance for the whole model. This removes the two main sources of fragility in from-scratch kernels.

設計方針・調査記録・ロードマップは [DESIGN.md](DESIGN.md) (日本語) を参照。

## Status

v0.2.0 — Phase 0.5 complete (API hygiene):

- Analytic primitives: `Plane`, `Cylinder`, `Line3`, `Circle3`, `Ellipse3`
- Self-contained 3-D math module (`Point3`, `Vec3`, `Unit3`) — no third-party linear algebra dependency
- All constructors return `Result` — panic-free public API
- Private fields with value-passing accessors; `#[non_exhaustive]` on all public enums
- `KernelError` with hand-written `Display` / `std::error::Error`
- Optional `serde` feature for serialization of math types
- Tolerant predicates (`Tol`, 3-value classification `Sign3`)
- Closed-form intersections: plane × plane, plane × cylinder (all five cases), line × plane
- Verification tests against hand-computed analytic solutions

Planned next (see [DESIGN.md](DESIGN.md) §10): half-edge topology → extruded solids → half-space cuts → boolean difference (2.5D-first) → section drawings → mass properties / quantity takeoff → tessellation.

## Usage

```rust
use archi_kernel::intersect::{plane_cylinder, PlaneCylinder};
use archi_kernel::{Cylinder, Line3, Plane, Point3, Tol, Vec3};

let tol = Tol::default();
let column = Cylinder::new(
    Line3::new(Point3::origin(), Vec3::Z).expect("valid axis"),
    0.3,
)
.expect("valid cylinder");
let slab = Plane::new(Point3::new(0.0, 0.0, 2.8), Vec3::Z).expect("valid plane");

match plane_cylinder(&slab, &column, &tol) {
    PlaneCylinder::Circle(c) => println!("section circle r = {} m", c.radius()),
    other => println!("{other:?}"),
}
```

## Development

```bash
cargo test --all-features                                   # analytic verification + serde smoke tests
cargo clippy --all-targets --all-features -- -D warnings    # zero-warning policy
cargo fmt --all -- --check
```

All lengths are SI metres, angles in radians. Unit conversion for display is the caller's responsibility.

## License

MIT — see `LICENSE`.
