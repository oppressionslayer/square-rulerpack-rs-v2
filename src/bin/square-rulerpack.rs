// `cargo run --release --example factor <n>`: the complete factorisation of n,
// small primes by trial division, perfect powers, then the self-initialising quadratic sieve
// on every composite piece until all are prime.
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args == ["--help"] || args == ["-h"] {
        println!("Usage: square-rulerpack <positive-decimal-integer>\n\nUse RAYON_NUM_THREADS to select the worker count. SIQS_Q2 and SIQS_SMALL\nretain the factoring engine's environment settings.");
        return;
    }
    if args == ["--version"] {
        println!("square-rulerpack {}", env!("CARGO_PKG_VERSION"));
        return;
    }
    if args.len() != 1 || args[0].is_empty() || !args[0].bytes().all(|c| c.is_ascii_digit()) {
        eprintln!("Usage: square-rulerpack <positive-decimal-integer>");
        std::process::exit(2);
    }
    let arg = &args[0];
    // progress lines with an estimate on stderr, unless the caller turned them off
    if std::env::var("SIQS_PROGRESS").is_err() { std::env::set_var("SIQS_PROGRESS", "1"); }
    let n = BigUint::parse_bytes(arg.as_bytes(), 10).expect("a decimal number");
    if n.is_zero() {
        eprintln!("The input must be positive.");
        std::process::exit(2);
    }
    let t = Instant::now();
    let factors = square_rulerpack::siqs::factor_completely(&n);
    let shown: Vec<String> = factors.iter().map(|f| f.to_string()).collect();
    println!("{} = {}", arg, if shown.is_empty() { "1".to_string() } else { shown.join(" × ") });
    let back = factors.iter().fold(BigUint::one(), |acc, f| acc * f);
    let product_ok = back == n;
    let all_probable_prime = factors.iter().all(|p| p.to_u64().map(square_rulerpack::arith64::is_prime_u64).unwrap_or_else(|| rulerpack::big::is_prime_big(p)));
    println!("{} factors, product checks: {}, all probable prime: {}, {:.2} s",
        factors.len(), product_ok, all_probable_prime, t.elapsed().as_secs_f64());
    // A deliberately short collection can leave an unresolved composite. Multiplication
    // alone does not certify a complete factorization, so do not label that result prime.
    if !product_ok || !all_probable_prime { std::process::exit(1); }
}
