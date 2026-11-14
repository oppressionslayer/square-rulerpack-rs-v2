# Square Rulerpack

An optimized self-initializing quadratic sieve in Rust by Lars Rocha. Uses
Rulerpack for modular arithmetic, primality tests and matrix dependencies.

```sh
SIQS_Q2=1 SIQS_SMALL=512 RAYON_NUM_THREADS=8 cargo run --release --example factor -- 1812159225242518635471100871179858454432957403510435424598139804492783352129850030361426871773756197
```

The sibling `../rulerpack-rs` source directory is included in the source bundle.
Both directories are needed when building this checkout. Registry packages use
the declared `rulerpack` version instead of the local path.

## Command-line installation

From this directory:

```sh
cargo install --path . --locked
SIQS_Q2=1 SIQS_SMALL=512 RAYON_NUM_THREADS=8 square-rulerpack 1009732533765211
```

The example and installed command share one entry point. Both print factors,
check their product, test primality, and report elapsed seconds. The factors of
the small example are `11344301` and `89007911`.

## Rust API

```rust
use num_bigint::BigUint;
use square_rulerpack::siqs;

let n = BigUint::from(1009732533765211u64);
let factors = siqs::factor_completely(&n);
assert_eq!(factors, vec![BigUint::from(11344301u32), BigUint::from(89007911u32)]);
```

The result is sorted and retains repeated factors. If the sieve cannot split a
composite, the library retains it as an unresolved cofactor. Verify completeness
as well as the product. The CLI performs both checks and exits nonzero on failure.
Large factors pass probable-prime tests; these are not primality certificates.

## Settings

| Setting | Meaning |
| --- | --- |
| `RAYON_NUM_THREADS=8` | Use eight factoring workers |
| `SIQS_Q2=1` | Enable Q2 when the multiplied input is 1 modulo 8 |
| `SIQS_SMALL=512` | Small-prime bound for the exact-threshold path |
| `SIQS_PROGRESS=0` | Suppress periodic progress output |

Q1 is used when Q2 is ineligible. The progress header distinguishes the input
from the multiplied sieve target. SIMD is selected at runtime; portable paths
remain available. The `profile` feature enables additional stage counters.

Always use a release build when measuring speed. The supplied 100-digit example
was reported at 794.20 seconds with eight workers on the development machine;
performance varies with hardware, load and input.

## License

MIT. See LICENSE. Related techniques and implementations are credited in the
source comments, including YAFU, msieve and Montgomery's block Lanczos method.
