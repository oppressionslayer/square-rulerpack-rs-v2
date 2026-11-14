//! Self-initializing quadratic sieve with double large primes.
//!
//! Q1 uses g(x)=a*x²+2*b*x+(b²-K)/a and witness r=a*x+b.
//! Q2 uses g(x)=a*x²+b*x+(b²-K)/(4*a) and witness r=2*a*x+b.
//! The full square offset is respectively a*g(x) or 4*a*g(x).
//! Families share a square-free a and visit b values in Gray-code order.
//! Modular root updates, logarithmic sieving, exact candidate division,
//! large-prime cycles and GF(2) dependencies lead to verified GCD factors.

use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

use rulerpack::qs::{dependencies_sparse, legendre, sqrt_mod};
use rulerpack::smooth::primes_upto;
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

/// Nanoseconds spent in each stage, summed over threads, and the candidates seen: the sieve
/// fill, the scan for hits, the exact division of hits, the polynomial switches. Read by
/// `examples/siqs_sweep.rs`; reset with `stats_reset`.
pub static STATS: [AtomicU64; 24] = [const { AtomicU64::new(0) }; 24];
pub const STAT_NAMES: [&str; 24] = ["sieve fill", "scan", "exact division", "switch", "candidates", "direct sieve ns", "cofactor split ns", "dlp range", "dlp prp dropped", "dlp split", "fill p<8192 ns", "fill p>=8192 ns", "cand: Q eval ns", "cand: small loop ns", "cand: divide out ns", "cand: cofactor+y ns", "reports (exact mode)", "fill setup ns", "tag ns (all)", "tag wide ns", "shift ns", "mid add ns", "wide add ns", ""];
pub fn stats_reset() { for s in &STATS { s.store(0, Relaxed); } }
/// The clock: `std::time::Instant` natively; in a browser (`wasm32-unknown-unknown`) there is
/// no monotonic clock in `std`, so the stage counters and the progress estimate read zero and
/// the sieve itself runs unchanged.
pub mod clock {
    #[cfg(not(target_arch = "wasm32"))]
    pub use std::time::Instant;
    #[cfg(target_arch = "wasm32")]
    #[derive(Clone, Copy, Debug)]
    pub struct Instant;
    #[cfg(target_arch = "wasm32")]
    impl Instant {
        pub fn now() -> Instant { Instant }
        pub fn elapsed(&self) -> std::time::Duration { std::time::Duration::ZERO }
        pub fn duration_since(&self, _earlier: Instant) -> std::time::Duration { std::time::Duration::ZERO }
    }
}

/// Whether the stage counters are kept. Off by default: the `profile` feature turns them on for
/// `siqs_sweep`. An earlier parameter audit asked for the cost of measuring to be measured.
pub const PROFILE: bool = cfg!(feature = "profile");
/// A time stamp that costs nothing unless `PROFILE` is on.
#[derive(Clone, Copy)]
struct Stamp(Option<clock::Instant>);
impl Stamp {
    #[inline(always)]
    fn now() -> Stamp { Stamp(if PROFILE { Some(clock::Instant::now()) } else { None }) }
    #[inline(always)]
    fn ns(self) -> u64 { self.0.map_or(0, |t| t.elapsed().as_nanos() as u64) }
    #[inline(always)]
    fn since(self, earlier: Stamp) -> u64 { match (self.0, earlier.0) { (Some(a), Some(b)) => a.duration_since(b).as_nanos() as u64, _ => 0 } }
}
/// Optional observer for a page's worker (the same shape as rulerpack 0.3.0's browser build, so
/// Spindexin's Rulerpack Lab can drive either engine): no clock, no allocation in the hook.
/// Stages: 0 preparing, 1 relation collection, 2 matrix, 3 square reconstruction.
pub type ProgressHook = fn(u32, usize, usize, usize, usize);
std::thread_local! {
    static PROGRESS_HOOK: std::cell::Cell<Option<ProgressHook>> = const { std::cell::Cell::new(None) };
}
pub fn set_progress_hook(hook: Option<ProgressHook>) { PROGRESS_HOOK.with(|h| h.set(hook)); }
fn report_progress(stage: u32, full: usize, need: usize, partials: usize, polys: usize) {
    if let Some(hook) = PROGRESS_HOOK.with(|h| h.get()) { hook(stage, full, need, partials, polys); }
}
#[inline(always)]
fn tick(i: usize, t: Stamp) { if PROFILE { STATS[i].fetch_add(t.ns(), Relaxed); } }
#[inline(always)]
fn count(i: usize, v: u64) { if PROFILE { STATS[i].fetch_add(v, Relaxed); } }

/// A factor-base prime with the root of n, its log in the byte scale, and the pair for the
/// Granlund–Montgomery exact-divisibility test: p | d iff d·inv (mod 2³²) ≤ limit, no division.
#[derive(Clone, Debug)]
pub struct Prime { pub p: u64, pub root: u64, pub log: u8, pub inv: u32, pub limit: u32 }

impl Prime {
    fn new(p: u64, root: u64) -> Prime {
        let inv = if p % 2 == 1 { inverse_mod_2_32(p as u32) } else { 0 };
        Prime { p, root, log: ((p as f64).log2() + 0.5) as u8, inv, limit: (u32::MAX / p as u32) }
    }
    /// Does p divide d? For odd p by the multiply-and-compare test, for 2 by the bit.
    #[inline]
    pub fn divides(&self, d: u32) -> bool { if self.p == 2 { d & 1 == 0 } else { d.wrapping_mul(self.inv) <= self.limit } }
}

/// The inverse of an odd number modulo 2³², by Newton's iteration.
fn inverse_mod_2_32(p: u32) -> u32 {
    let mut x = p; // correct to 3 bits for odd p
    for _ in 0..5 { x = x.wrapping_mul(2u32.wrapping_sub(p.wrapping_mul(x))); }
    x
}

/// Primes below this are left out of the sieve fill (their strides touch every few cells and
/// cost more than they tell); the threshold is lowered by their expected contribution instead.
/// `SIQS_SMALL` in the environment overrides it, for the sweeps.
/// The sieve block, in cells. Measured on the 60-digit input: the whole 2^17 interval at once
/// is fastest (1.28 s), 32 KiB blocks slower (1.39 s), 8 KiB much slower (2.2 s), so the fill is
/// not cache-latency bound and the default is one block; the blocks stay for a factor base that
/// outgrows one. `SIQS_SEGMENT` in the environment overrides it, for the sweeps.
pub fn segment_len(bucket: bool, primes: usize) -> usize {
    // a bucket entry packs the prime index above the position in the block, so a base past
    // 131,072 primes takes the 16 KB block
    let seg: usize = std::env::var("SIQS_SEGMENT").ok().and_then(|v| v.parse().ok()).unwrap_or(if !bucket { 1 << 20 } else if primes < 1 << 17 { 1 << 15 } else { 1 << 14 });
    seg.next_power_of_two()
}

/// The direct sieve for the primes from `SIQS_DIRECT_FROM` up (the bucket sieve's place in
/// 0.3.0, hence the name): `SIQS_BUCKET=1` forces it on, `0` off, and unset leaves it to the
/// rule in the sieve (on once the largest prime is wider than the interval). With it the block
/// is 32 KB, the L1 size, instead of the whole interval.
pub fn bucket_sieve() -> Option<bool> { std::env::var("SIQS_BUCKET").ok().map(|v| v != "0") }

/// The DLP exponent in thousandths, set from the table row by `Params::for_bits` (YAFU's
/// `dlp_exp` column: 1.72 at 268 bits, 1.8 at 330) and read by `split_cofactor`;
/// `SIQS_DLP_EXP` overrides it.
pub static DLP_EXP_MILLI: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(1800);
pub fn dlp_exp() -> f64 {
    static OVERRIDE: std::sync::OnceLock<Option<f64>> = std::sync::OnceLock::new();
    OVERRIDE.get_or_init(|| std::env::var("SIQS_DLP_EXP").ok().and_then(|v| v.parse().ok())).unwrap_or(DLP_EXP_MILLI.load(Relaxed) as f64 / 1000.0)
}

/// Two large primes on (the double-large-prime variation), unless `SIQS_DLP=0`.
pub fn double_large() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("SIQS_DLP").map(|v| v != "0").unwrap_or(true))
}

/// Miller–Rabin for u64, deterministic with the first twelve prime bases.
pub fn is_prime_u64(n: u64) -> bool {
    if n < 2 { return false; }
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] { if n % p == 0 { return n == p; } }
    let (mut d, mut s) = (n - 1, 0);
    while d % 2 == 0 { d /= 2; s += 1; }
    'outer: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = pow_mod_u64(a, d, n);
        if x == 1 || x == n - 1 { continue; }
        for _ in 1..s { x = mul_mod_u64(x, x, n); if x == n - 1 { continue 'outer; } }
        return false;
    }
    true
}
fn mul_mod_u64(a: u64, b: u64, m: u64) -> u64 { ((a as u128 * b as u128) % m as u128) as u64 }
fn pow_mod_u64(mut b: u64, mut e: u64, m: u64) -> u64 { let mut r = 1u64; b %= m; while e > 0 { if e & 1 == 1 { r = mul_mod_u64(r, b, m); } b = mul_mod_u64(b, b, m); e >>= 1; } r }

/// A nontrivial factor of a composite u64 by Pollard's rho with Brent's cycle finding.
pub fn rho_u64(n: u64) -> Option<u64> {
    if n % 2 == 0 { return Some(2); }
    for c in 1..20u64 {
        let (mut y, mut r, mut q, mut g) = (2u64, 1u64, 1u64, 1u64);
        let (mut x, mut ys) = (y, y);
        let f = |v: u64| (mul_mod_u64(v, v, n) + c) % n;
        while g == 1 {
            x = y;
            for _ in 0..r { y = f(y); }
            let mut k = 0;
            while k < r && g == 1 {
                ys = y;
                for _ in 0..128.min(r - k) { y = f(y); q = mul_mod_u64(q, x.abs_diff(y), n); }
                g = gcd_u64(q, n);
                k += 128;
            }
            r *= 2;
            if r > 1 << 20 { break; }
        }
        if g == n { loop { ys = f(ys); g = gcd_u64(x.abs_diff(ys), n); if g > 1 { break; } } }
        if g > 1 && g < n { return Some(g); }
    }
    None
}
fn gcd_u64(mut a: u64, mut b: u64) -> u64 { while b != 0 { let t = a % b; a = b; b = t; } a }

/// Primes from this size up to the block are resieved: after the scan each walks the interval
/// once and tags the candidates it hits, so a candidate divides only by primes known to divide
/// it instead of testing every block-sized prime (YAFU resieves 13-16-bit primes the same way).
/// `SIQS_RESIEVE=0` turns it off; default 8192.
pub fn resieve_min() -> u64 { std::env::var("SIQS_RESIEVE").ok().and_then(|v| v.parse().ok()).unwrap_or(8192) }

pub fn small_bound() -> u64 { small_bound_with(exact_mode()) }
/// The unsieved bound: 512 with the exact threshold, 30 with the plain threshold.
fn small_bound_with(exact: bool) -> u64 { std::env::var("SIQS_SMALL").ok().and_then(|v| v.parse().ok()).unwrap_or(if exact { 512 } else { 30 }) }
/// Inputs below this many bits (about 65 digits) keep the plain threshold: a polynomial is so
/// cheap there that report validation can outweigh the benefit of additional relations.
pub const EXACT_MIN_BITS: u64 = 215;

/// **The exact threshold** (on by default from `EXACT_MIN_BITS`; `SIQS_EXACT=0` restores the
/// plain rule), after YAFU's small-prime
/// variation and trial-division cutoff (`tdiv_small.c`, Ben Buhrow, public domain).
///
/// The plain rule keeps a sieve report when its sum of logs reaches one constant, having
/// assumed that the unsieved small primes contribute their *average* and that Q(x) has its
/// *largest* size. Both assumptions lose relations: a value rich in small primes or in a power
/// of 2 has a sieve sum several bits short of what its factor-base part really is, and a
/// position near a root of the polynomial has a value several bits smaller than the largest.
///
/// With the exact threshold a report is kept when
/// `sieve sum + exact log of its small-prime part ≥ log2|g(x)| − allowed remainder`,
/// where the small-prime part is counted with multiplicity (`SmallPart`), log2|g(x)| is that
/// position's own, and the allowed remainder is the double-large-prime ceiling plus
/// `SIQS_MARGIN` bits (default 1) for the rounding of the byte logs. The sieve itself starts
/// each byte from a per-position template instead of one constant, lowered by `SIQS_ECAP` bits
/// (default 32) so that a report rich in small primes reaches the test at all.
/// Because the small primes are then accounted for exactly, they need not be sieved: the
/// unsieved bound defaults to 512 in exact mode (`SIQS_SMALL`).
/// Relations are found by exact division as before, so the mode changes only which positions are tried.
pub fn exact_mode() -> bool {
    if let Some(v) = EXACT_OVERRIDE.with(|o| o.get()) { return v; }
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("SIQS_EXACT").map(|v| v != "0").unwrap_or(true))
}
std::thread_local! {
    /// Per-thread override of `SIQS_EXACT`, for tests (the environment is process-wide).
    static EXACT_OVERRIDE: std::cell::Cell<Option<bool>> = const { std::cell::Cell::new(None) };
}
/// Bits the sieve threshold is lowered by so that reports rich in small primes are seen (`SIQS_ECAP`).
pub fn exact_ecap() -> f32 {
    static V: std::sync::OnceLock<f32> = std::sync::OnceLock::new();
    *V.get_or_init(|| std::env::var("SIQS_ECAP").ok().and_then(|v| v.parse().ok()).unwrap_or(32.0))
}
/// Bits of remainder allowed beyond the large-prime ceiling (`SIQS_MARGIN`).
pub fn exact_margin() -> f32 {
    static V: std::sync::OnceLock<f32> = std::sync::OnceLock::new();
    *V.get_or_init(|| std::env::var("SIQS_MARGIN").ok().and_then(|v| v.parse().ok()).unwrap_or(1.0))
}

/// log2 of a positive finite double to within 0.01: the exponent, plus the mantissa fraction f
/// with the correction f + 0.344·f·(1 − f) for log2(1 + f). The report test compares against a
/// margin of bits, so a library logarithm per report (about 20 ns) is not needed.
#[inline(always)]
fn log2_fast(g: f64) -> f32 {
    let b = g.to_bits();
    let e = ((b >> 52) & 0x7ff) as i32 - 1023;
    let f = ((b >> 29) & 0x7f_ffff) as f32 * (1.0 / (1u32 << 23) as f32);
    e as f32 + f + 0.344 * f * (1.0 - f)
}

/// v mod m for a signed v, in [0, m), without allocating.
#[inline]
fn rem_floor_u32(v: &BigInt, m: u32) -> u32 {
    let r = rem_u32(v.magnitude(), m);
    if v.is_negative() && r != 0 { m - r } else { r }
}
/// The low 64 bits of a signed v in two's complement.
#[inline]
fn low64(v: &BigInt) -> u64 {
    let lo = v.magnitude().iter_u64_digits().next().unwrap_or(0);
    if v.is_negative() { lo.wrapping_neg() } else { lo }
}

/// The unsieved small primes of one polynomial family, for the exact threshold: the polynomial
/// g(x) = (a·x + lin)·x + c reduced modulo a power of each of them and modulo 2^64, which is
/// enough to read off, in machine words, which of them divide g(x) and how many times.
pub struct SmallPart {
    /// per odd unsieved prime: its index in the base, the prime, p^e ≤ 2^31, log2 p, and the
    /// multiply-and-compare divisibility pair of `Prime`
    idx: Vec<u32>, p: Vec<u32>, m: Vec<u32>, lg: Vec<f32>, inv: Vec<u32>, lim: Vec<u32>,
    /// e of p^e, and the 64-bit divisibility pair: p⁻¹ mod 2^64 and ⌊(2^64 − 1)/p⌋, with which
    /// "p divides g" is one multiply and a compare and the quotient is the product itself
    e: Vec<u32>, pinv: Vec<u64>, plim: Vec<u64>,
    /// ⌊(2^64 − 1)/p^e⌋ for Barrett reduction modulo p^e
    mu: Vec<u64>,
    /// fixed for the family: 2a, the discriminant lin² − 4ac and (4a)⁻¹, all mod p^e (the
    /// inverse is 0 for a prime that divides a, which is then evaluated from a, lin and c)
    a2: Vec<u32>, disc: Vec<u32>, i4a: Vec<u32>,
    /// a mod p^e (fixed for the family), then lin mod p^e (set per polynomial), and c mod p^e
    /// for the primes that divide a only
    a: Vec<u32>, lin: Vec<u32>, c: Vec<u32>,
    /// the same three modulo 2^64
    a_lo: u64, lin_lo: u64, c_lo: u64,
    /// whether 2 is among the unsieved primes
    two: bool,
}

impl SmallPart {
    /// The base primes below `small`, with `a` and the family's discriminant `disc` = lin² − 4ac
    /// (the same for every polynomial of the family) reduced for each.
    pub fn new(base: &[Prime], small: u64, a: &BigInt, disc: &BigInt) -> SmallPart {
        let mut s = SmallPart { idx: vec![], p: vec![], m: vec![], lg: vec![], inv: vec![], lim: vec![], e: vec![], pinv: vec![], plim: vec![], mu: vec![], a2: vec![], disc: vec![], i4a: vec![], a: vec![], lin: vec![], c: vec![], a_lo: low64(a), lin_lo: 0, c_lo: 0, two: false };
        for (j, q) in base.iter().enumerate() {
            if q.p >= small { break; }
            if q.p == 2 { s.two = true; continue; }
            let p = q.p as u32;
            let (mut m, mut e) = (p, 1u32);
            while (m as u64) * (p as u64) <= 1 << 31 { m *= p; e += 1; }
            let mut pinv = q.p;   // Newton: correct to 3 bits for odd p, doubling each step
            for _ in 0..6 { pinv = pinv.wrapping_mul(2u64.wrapping_sub(q.p.wrapping_mul(pinv))); }
            s.idx.push(j as u32); s.p.push(p); s.m.push(m); s.lg.push((q.p as f32).log2()); s.inv.push(q.inv); s.lim.push(q.limit);
            s.e.push(e); s.pinv.push(pinv); s.plim.push(u64::MAX / q.p);
            s.mu.push(u64::MAX / m as u64);
            let am = rem_floor_u32(a, m) as u64;
            s.a.push(am as u32);
            s.a2.push((2 * am % m as u64) as u32);
            s.disc.push(rem_floor_u32(disc, m));
            // (4a)⁻¹ mod p^e when p does not divide a: the i64 extended Euclid of `modinv`
            s.i4a.push(if am % q.p == 0 { 0 } else { modinv(4 * am % m as u64, m as u64) as u32 });
        }
        s.lin = vec![0; s.p.len()];
        s.c = vec![0; s.p.len()];
        s
    }
    /// The coefficients of this polynomial: g(x) = (a·x + lin)·x + c. Only lin is reduced per
    /// prime: 4a·g(x) = (2a·x + lin)² − disc gives g modulo p^e without c.
    pub fn set_poly(&mut self, lin: &BigInt, c: &BigInt) {
        self.lin_lo = low64(lin);
        self.c_lo = low64(c);
        for t in 0..self.p.len() {
            self.lin[t] = rem_floor_u32(lin, self.m[t]);
            if self.i4a[t] == 0 { self.c[t] = rem_floor_u32(c, self.m[t]); }
        }
    }
    /// v mod p^e of the t-th prime for v < 2^63, by Barrett: one multiply-high instead of a division.
    #[inline(always)]
    fn reduce(&self, t: usize, v: u64) -> u64 {
        let m = self.m[t] as u64;
        let q = ((v as u128 * self.mu[t] as u128) >> 64) as u64;
        let mut r = v - q * m;   // the quotient estimate is short by at most 2
        if r >= m { r -= m; }
        if r >= m { r -= m; }
        r
    }
    /// log2 of the part of g(x) made of the unsieved primes, each with its multiplicity.
    /// `pos` is the sieve position of x and `starts` the position roots of every prime, used
    /// to skip the primes that do not divide g(x); a prime with no roots here (a factor of a,
    /// marked u32::MAX) is evaluated directly.
    #[inline]
    pub fn bits(&self, xi: i64, pos: u32, starts: &[[u32; 2]]) -> f32 {
        let mut e = 0f32;
        if self.two {
            let x = xi as u64;
            let g = self.a_lo.wrapping_mul(x).wrapping_add(self.lin_lo).wrapping_mul(x).wrapping_add(self.c_lo);
            e += g.trailing_zeros().min(48) as f32;
        }
        // which primes divide g(x): one bit each, without a branch per prime (a mispredict
        // costs more than the test), then the multiplicity of the few that do
        let n = self.p.len();
        let mut t0 = 0;
        while t0 < n {
            let t1 = (t0 + 64).min(n);
            let mut mask = 0u64;
            for t in t0..t1 {
                let [s0, s1] = starts[self.idx[t] as usize];
                let (inv, lim) = (self.inv[t], self.lim[t]);
                let hit = (s0 == u32::MAX) | (pos.abs_diff(s0).wrapping_mul(inv) <= lim) | (pos.abs_diff(s1).wrapping_mul(inv) <= lim);
                mask |= (hit as u64) << (t - t0);
            }
            while mask != 0 {
                let t = t0 + mask.trailing_zeros() as usize;
                e += self.valuation(t, xi) as f32 * self.lg[t];
                mask &= mask - 1;
            }
            t0 = t1;
        }
        debug_assert!((e - self.bits_direct(xi)).abs() < 1e-3, "root test and direct evaluation disagree at x = {xi}");
        e
    }
    /// `bits` without the root test: every unsieved prime evaluated directly. The reference
    /// for the tests and for the debug check in `bits`.
    pub fn bits_direct(&self, xi: i64) -> f32 {
        let mut e = 0f32;
        if self.two {
            let x = xi as u64;
            let g = self.a_lo.wrapping_mul(x).wrapping_add(self.lin_lo).wrapping_mul(x).wrapping_add(self.c_lo);
            e += g.trailing_zeros().min(48) as f32;
        }
        for t in 0..self.p.len() { e += self.valuation(t, xi) as f32 * self.lg[t]; }
        e
    }
    /// How many times the t-th unsieved odd prime divides g(x) (capped at its stored power).
    #[inline]
    fn valuation(&self, t: usize, xi: i64) -> u32 {
        let m = self.m[t] as u64;
        // x mod p^e: p^e exceeds 2^23 for p < 256, far above any |x| of the interval
        let x = if xi.unsigned_abs() < m { if xi < 0 { (xi + m as i64) as u64 } else { xi as u64 } } else { xi.rem_euclid(m as i64) as u64 };
        let mut g = if self.i4a[t] != 0 {
            // g ≡ ((2a·x + lin)² − disc)·(4a)⁻¹ (mod p^e)
            let y = self.reduce(t, self.a2[t] as u64 * x + self.lin[t] as u64);
            let w = self.reduce(t, y * y + m - self.disc[t] as u64);
            self.reduce(t, w * self.i4a[t] as u64)
        } else {
            ((self.a[t] as u64 * x + self.lin[t] as u64) % m * x + self.c[t] as u64) % m
        };
        if g == 0 { return self.e[t]; }
        // exact division by an odd p is multiplication by its inverse mod 2^64, and the product
        // is at most ⌊(2^64 − 1)/p⌋ exactly when p divided: no division in the loop
        let (pinv, lim) = (self.pinv[t], self.plim[t]);
        let mut k = 0;
        loop { let q = g.wrapping_mul(pinv); if q > lim { return k; } g = q; k += 1; }
    }
}

/// The knobs: factor-base size, half the sieve interval, large-prime multiple, threshold slack
/// in largest-prime logs (0.7 everywhere: stricter was faster on every input tried, as also
/// found on Alpern's sieve; 1.3 cost a third more and 2.0 four times), relations past the base.
#[derive(Clone, Copy, Debug)]
pub struct Params { pub base: usize, pub half: usize, pub large: u64, pub slack: f32, pub extra: usize }

impl Params {
    /// By the size of n in digits, from the usual tables, measured here where it says so.
    pub fn for_digits(d: usize) -> Params {
        Params::for_bits((d as f64 / std::f64::consts::LOG10_2).round() as u32, d)
    }
    /// By the exact size of n: `for_bits(n.bits(), digits)`. An earlier parameter audit found the
    /// digit estimate landing on 332 bits for a 330-bit input, which picked a 25 % wider window
    /// from the YAFU tables; `factor` now passes the real bit count, as YAFU does.
    pub fn for_n(n: &BigUint) -> Params { Params::for_bits(n.bits() as u32, n.to_string().len()) }
    pub fn for_bits(bits: u32, d: usize) -> Params {
        // `SIQS_TABLE=yafu` takes the row from YAFU's direct-sieve table (crate::params), which
        // assumes double large primes: blocks per side of 32 KB become the half interval, the
        // large-prime multiplier carries over, and the threshold slack is `SIQS_SLACK` (1.3).
        // `SIQS_TABLE=yafu512` uses the AVX-512 table's rows.
        // The YAFU-derived AVX2 rows are the default from 2026-11-12 (908.9 s against the hand
        // table's slower rows at 100 digits on the desktop, and the browser has no environment
        // to set them); `SIQS_TABLE=own` keeps the hand table below, `yafu512` the AVX-512 rows.
        let t = std::env::var("SIQS_TABLE").unwrap_or_else(|_| "yafu".to_string());
        {
            let table = match t.as_str() { "yafu" => Some(crate::params::YAFU_DIRECT), "yafu512" => Some(crate::params::YAFU_DIRECT_AVX512), _ => None };
            // the tables start where YAFU's do; smaller inputs use the hand rows
            let table = table.filter(|tb| tb.first().is_some_and(|r| bits >= r.bits));
            if let Some(table) = table {
                let r = crate::params::for_bits(table, bits);
                if r.dlp_exp > 0.0 { DLP_EXP_MILLI.store((r.dlp_exp * 1000.0).round() as u32, Relaxed); }
                let slack = std::env::var("SIQS_SLACK").ok().and_then(|v| v.parse().ok()).unwrap_or(1.3);
                return Params { base: r.base as usize, half: r.blocks as usize * 32768, large: r.lp_mult as u64, slack, extra: 64 };
            }
        }
        match d {
            0..=30 => Params { base: 150, half: 1 << 13, large: 30, slack: 0.7, extra: 24 },
            31..=40 => Params { base: 400, half: 1 << 14, large: 40, slack: 0.7, extra: 32 },
            41..=50 => Params { base: 1100, half: 1 << 15, large: 50, slack: 0.7, extra: 48 },
            51..=56 => Params { base: 1800, half: 1 << 15, large: 60, slack: 0.7, extra: 64 },   // measured on the 52-digit case
            57..=62 => Params { base: 3000, half: 1 << 16, large: 60, slack: 0.7, extra: 64 },
            63..=67 => Params { base: 4500, half: 1 << 16, large: 80, slack: 0.7, extra: 80 },
            68..=72 => Params { base: 15000, half: 1 << 16, large: 100, slack: 0.7, extra: 96 },   // 69-digit sweep with the bucket: 7,000/2^17 3.8 s to the count, 15,000/2^16 3.2 s
            73..=77 => Params { base: 30000, half: 1 << 17, large: 120, slack: 0.7, extra: 120 },   // 75-digit sweep: 10,000 26 s to the count, 25,000 13 s, 40,000 12 s
            // with the bucket sieve the exact division no longer grows with the base, and the
            // budgeted sweeps of 2026-11-09 moved to much larger bases: 80 digits 20,000/2^18
            // 217 relations/s → 50,000/2^17 1,200/s (41 s to the count against 92); 90 digits
            // 50,000 at 35/s → 80,000 at 141/s; 100 digits 120,000 at 14.6/s → 200,000 at 27/s
            78..=84 => Params { base: 50000, half: 1 << 17, large: 100, slack: 0.7, extra: 150 },
            85..=92 => Params { base: 80000, half: 1 << 17, large: 200, slack: 0.7, extra: 200 },
            93..=97 => Params { base: 120000, half: 1 << 17, large: 250, slack: 0.7, extra: 250 },
            _ => Params { base: 200000, half: 1 << 17, large: 300, slack: 0.7, extra: 300 },
        }
    }
}

/// **The Knuth–Schroeppel multiplier**: the small square-free k for which k·n has the most
/// small primes in its factor base, scored as Σ g(p)·ln p over the primes to 1000 less ½ ln k,
/// with g = 2/(p−1) for a prime where k·n is a square, 1/(p−1) for a prime dividing k, and the
/// usual cases at 2. The sieve then factors k·n; the gcd at the end is taken with n.
pub fn multiplier(n: &BigUint) -> u64 {
    if let Some(k) = std::env::var("SIQS_MULT").ok().and_then(|v| v.parse().ok()) { return k; }
    let candidates = [1u64, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38, 39, 41, 42, 43, 46, 47, 51, 53, 55, 57, 58, 59, 61, 62, 65, 66, 67, 69, 70, 71, 73];
    let primes = primes_upto(1000);
    let mut best = (f64::MIN, 1u64);
    for &k in &candidates {
        let kn = n * k;
        let mut score = -0.5 * (k as f64).ln();
        let m8 = (&kn % 8u32).to_u64().unwrap();
        score += match m8 { 1 => 2.0 * 2f64.ln(), 5 => 2f64.ln(), 3 | 7 => 0.5 * 2f64.ln(), _ => 0.0 };
        for &p in primes.iter().skip(1) {
            let r = (&kn % p).to_u64().unwrap();
            let lp = (p as f64).ln();
            if k % p == 0 { score += lp / (p - 1) as f64; } else if legendre(r, p) == 1 { score += 2.0 * lp / (p - 1) as f64; }
        }
        if score > best.0 { best = (score, k); }
    }
    best.1
}

/// The first `count` primes modulo which n is a square (2 handled as a stride of its own).
pub fn factor_base(n: &BigUint, count: usize) -> Vec<Prime> {
    let mut out = vec![Prime::new(2, (n % 2u32).to_u64().unwrap())];
    let mut bound = 1000u64;
    let mut checked_through = 2u64;
    loop {
        for p in primes_upto(bound).into_iter().skip(1) {
            if out.len() >= count { return out; }
            if p <= checked_through { continue; }
            let np = (n % p).to_u64().unwrap();
            if np == 0 { out.push(Prime::new(p, 0)); } else if legendre(np, p) == 1 { out.push(Prime::new(p, sqrt_mod(np, p))); }
        }
        checked_through = bound;
        bound *= 2;
    }
}

/// One relation: (a·x + b) mod n, the exponents over the base as `(column, exponent)` pairs
/// (column 0 is the sign, column j + 1 the j-th base prime), the large primes outside the base
/// with their exponents (empty for a full relation).
#[derive(Clone, Debug)]
pub struct Relation { pub y: BigUint, pub exps: Vec<(u32, u32)>, pub large: Vec<(u64, u32)> }

impl Relation {
    /// The product of two relations: exponents summed, large primes summed.
    pub fn times(&self, other: &Relation, n: &BigUint) -> Relation {
        let mut exps: Vec<(u32, u32)> = self.exps.clone();
        for &(c, e) in &other.exps { match exps.iter_mut().find(|(d, _)| *d == c) { Some(t) => t.1 += e, None => exps.push((c, e)) } }
        let mut large = self.large.clone();
        for &(p, e) in &other.large { match large.iter_mut().find(|(q, _)| *q == p) { Some(t) => t.1 += e, None => large.push((p, e)) } }
        Relation { y: &self.y * &other.y % n, exps, large }
    }
    /// Is every exponent even, base and large primes alike?
    pub fn is_square(&self) -> bool { self.exps.iter().all(|(_, e)| e % 2 == 0) && self.large.iter().all(|(_, e)| e % 2 == 0) }
}

/// A polynomial family: a = ∏ q_i, the B_l such that b = Σ ±B_l, and the inverse
/// of a (Q1) or 2a (Q2) mod p per prime. Q2 additionally makes b odd.
struct Family { a: BigUint, qs: Vec<usize>, bs: Vec<BigUint>, ainv: Vec<u64>, q2: bool }

/// YAFU's Q2 variation: for kN = 1 (mod 8), sieve ((2ax+b)^2-kN)/(4a)
/// with odd b and half-sized a. Enabled by default for eligible K; SIQS_Q2=0 selects Q1.
fn q2_polynomials(n: &BigUint) -> bool {
    // on by default since 2026-11-13 (11 % fewer polynomials at 80 digits, 0.199 -> 0.222 partials
    // a polynomial at 100); `SIQS_Q2=0` restores the Q1 form. In the browser there is no
    // environment, so the default is what runs there.
    !std::env::var("SIQS_Q2").is_ok_and(|v| v == "0") && (n % 8u32) == BigUint::one()
}

/// The family setup in machine words (the default) or, with `SIQS_FAST_SETUP=0`, by the
/// big-integer remainders and 128-bit arithmetic used before 2026-11-13, for A/B runs.
///
/// Per family the setup computes, for every base prime, a⁻¹ mod p, one root shift per B_l and
/// the two starting roots. With `&a % p` (which builds a BigUint), a 128-bit extended Euclid
/// and 128-bit modular products that came to about 240 ms a family at 100 digits (90,000
/// primes, 14 B values): the 8.5 % of worker time the stage timers did not account for in
/// a stage profile. Every operand fits in 32 bits and every product in 64, so
/// the same values come from `rem_u32` and 64-bit arithmetic. Results are identical.
fn fast_setup() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var("SIQS_FAST_SETUP").map(|v| v != "0").unwrap_or(true))
}

fn modinv(a: u64, p: u64) -> u64 {
    if p < 1 << 62 && fast_setup() {
        // every intermediate of the extended Euclid is below 2p in absolute value: i64 is enough
        let (mut t, mut nt, mut r, mut nr) = (0i64, 1i64, p as i64, (a % p) as i64);
        while nr != 0 { let q = r / nr; (t, nt) = (nt, t - q * nt); (r, nr) = (nr, r - q * nr); }
        return t.rem_euclid(p as i64) as u64;
    }
    let (mut t, mut nt, mut r, mut nr) = (0i128, 1i128, p as i128, (a % p) as i128);
    while nr != 0 { let q = r / nr; (t, nt) = (nt, t - q * nt); (r, nr) = (nr, r - q * nr); }
    t.rem_euclid(p as i128) as u64
}

/// Pick a's primes: s of them from the middle of the base so that a ≈ √(2n)/M, chosen by a
/// seeded walk so every thread and every round takes a different set.
fn choose_a(n: &BigUint, base: &[Prime], half: usize, seed: u64) -> Family {
    // The target-A repair is on unless `SIQS_CLOSE_A=0`: measured here, 60 digits
    // 1.27 s → 0.74 s over three alternating pairs, 80 digits 121 s → 91 s
    choose_a_mode(n, base, half, seed, std::env::var("SIQS_CLOSE_A").map(|v| v != "0").unwrap_or(true))
}

fn choose_a_mode(n: &BigUint, base: &[Prime], half: usize, seed: u64, close: bool) -> Family {
    choose_a_form(n, base, half, seed, close, q2_polynomials(n))
}

fn choose_a_form(n: &BigUint, base: &[Prime], half: usize, seed: u64, close: bool, q2: bool) -> Family {
    assert!(!q2 || (n % 8u32) == BigUint::one());
    let step = if q2 { 2 } else { 1 };
    let target = (n.to_f64().unwrap() * 2.0).sqrt() / (half as f64 * step as f64);
    let target_log = target.ln();
    // candidates: primes around 2000, the usual choice; s from the target
    let lo = base.iter().position(|q| q.p > 1000).unwrap_or(base.len() / 2).min(base.len().saturating_sub(10));
    let hi = base.len().min(lo + 300).max(lo + 10).min(base.len());
    let pool: Vec<usize> = (lo..hi).collect();
    let mut rng = seed.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(0xD1B54A32D192ED03);
    let mut next = || { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; rng };
    let mut qs: Vec<usize> = Vec::new();
    let mut log_a = 0.0f64;
    while log_a < target_log - (base[pool[pool.len() / 2]].p as f64).ln() / 2.0 && qs.len() < 20 {
        let i = pool[(next() % pool.len() as u64) as usize];
        if qs.contains(&i) { continue; }
        qs.push(i); log_a += (base[i].p as f64).ln();
    }
    if qs.len() < 2 { qs.push(pool[0]); qs.push(pool[1]); }
    if close { repair_a(&mut qs, base, target_log); }
    qs.sort_unstable();
    let a: BigUint = qs.iter().fold(BigUint::one(), |acc, &i| acc * base[i].p);
    // B_l = (a/q_l) · (t_l · (a/q_l)⁻¹ mod q_l), so B_l ≡ t_l (mod q_l), B_l ≡ 0 mod the others
    let bs: Vec<BigUint> = qs.iter().map(|&i| {
        let q = base[i].p;
        let aq = &a / q;
        let aq_mod = (&aq % q).to_u64().unwrap();
        let gamma = (base[i].root as u128 * modinv(aq_mod, q) as u128 % q as u128) as u64;
        let gamma = if gamma > q / 2 { q - gamma } else { gamma };
        aq * gamma
    }).collect();
    let word = fast_setup() && base.last().is_some_and(|q| q.p <= u32::MAX as u64);
    let ainv: Vec<u64> = if word {
        base.iter().map(|q| { let am = rem_u32(&a, q.p as u32) as u64 * step % q.p; if am == 0 { 0 } else { modinv(am, q.p) } }).collect()
    } else {
        base.iter().map(|q| { let am = ((&a % q.p).to_u64().unwrap() * step) % q.p; if am == 0 { 0 } else { modinv(am, q.p) } }).collect()
    };
    Family { a, qs, bs, ainv, q2 }
}

/// Keep the number of factors and replace at most one by an unused odd base
/// prime to minimize the logarithmic distance from the usual SIQS target.
/// Selection is approximate; a, b and all relations are subsequently exact.
fn repair_a(qs: &mut [usize], base: &[Prime], target_log: f64) {
    let log_a: f64 = qs.iter().map(|&i| (base[i].p as f64).ln()).sum();
    let mut best = ((log_a - target_log).abs(), None);
    for (slot, &old) in qs.iter().enumerate() {
        let rest = log_a - (base[old].p as f64).ln();
        let wanted = (target_log - rest).exp();
        let at = base.partition_point(|q| (q.p as f64) < wanted);
        // Enough neighbors to step past the already-used factors.
        let radius = qs.len() + 2;
        for (j, q) in base.iter().enumerate().take((at + radius).min(base.len())).skip(at.saturating_sub(radius)) {
            if q.p == 2 || q.root == 0 || (j != old && qs.contains(&j)) { continue; }
            let error = (rest + (q.p as f64).ln() - target_log).abs();
            if error < best.0 { best = (error, Some((slot, j))); }
        }
    }
    if let Some((slot, j)) = best.1 { qs[slot] = j; }
}

/// Move a known reduced residue by a known reduced delta, with one correction.
fn shift_residue(r: u32, delta: u32, p: u32, add: bool) -> u32 {
    debug_assert!(p > 0 && r < p && delta < p);
    if add {
        if r >= p - delta { r - (p - delta) } else { r + delta }
    } else if r >= delta { r - delta } else { p - (delta - r) }
}


/// Shift every position by its delta modulo its prime, branch-free, over whole arrays.
///
/// The crate is built for baseline x86-64, so left to itself the compiler vectorises this
/// loop with SSE2: four lanes, and an unsigned minimum pieced together from compares because
/// SSE2 has none (checked in the shipped binary, 2026-11-13: no ymm register outside the
/// hand-written kernels). `shift_all_avx2` is the same loop compiled with AVX2 enabled, chosen
/// at run time: eight lanes and `vpminud`. `SIQS_NO_AVX2=1` keeps the baseline loop.
#[inline]
fn shift_all(r: &mut [u32], d: &[u32], p: &[u32], add: bool, avx2: bool) {
    #[cfg(target_arch = "x86_64")]
    if avx2 && avx512() {
        // SAFETY: `avx512()` checked avx512f at run time
        unsafe { return shift_all_avx512(r, d, p, add); }
    }
    #[cfg(target_arch = "x86_64")]
    if avx2 {
        // SAFETY: `avx2` was set from `is_x86_feature_detected!("avx2")`
        unsafe { return shift_all_avx2(r, d, p, add); }
    }
    let _ = avx2;
    shift_all_body(r, d, p, add)
}
#[cfg_attr(not(target_arch = "wasm32"), inline(always))]
#[cfg_attr(target_arch = "wasm32", inline(never))]   // a profile line of its own in the browser build
fn shift_all_body(r: &mut [u32], d: &[u32], p: &[u32], add: bool) {
    debug_assert!(r.len() == d.len() && r.len() == p.len());
    if add { for ((r, &d), &p) in r.iter_mut().zip(d).zip(p) { let t = r.wrapping_add(d); *r = t.min(t.wrapping_sub(p)); } }
    else { for ((r, &d), &p) in r.iter_mut().zip(d).zip(p) { let t = r.wrapping_sub(d); *r = t.min(t.wrapping_add(p)); } }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn shift_all_avx2(r: &mut [u32], d: &[u32], p: &[u32], add: bool) { shift_all_body(r, d, p, add) }

fn have_avx2() -> bool {
    #[cfg(target_arch = "x86_64")]
    { std::env::var("SIQS_NO_AVX2").is_err() && is_x86_feature_detected!("avx2") }
    #[cfg(not(target_arch = "x86_64"))]
    { false }
}

/// AVX-512 (F + BW + VL) for the hot vector kernels: the root shifts, the wide-prime compaction
/// and the mid-prime block walk run sixteen lanes at a time, and the compaction becomes one
/// `vpcompressd` into a register plus a full store instead of a table permute. Detected at run
/// time; `SIQS_AVX512=0` keeps the AVX2 kernels (the hit order and so every relation are the
/// same either way). Cached, since the sieve asks once per chunk.
#[cfg_attr(not(target_arch = "x86_64"), allow(dead_code))]
fn avx512() -> bool {
    static V: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *V.get_or_init(|| {
        #[cfg(target_arch = "x86_64")]
        { std::env::var("SIQS_AVX512").map(|v| v != "0").unwrap_or(true) && std::env::var("SIQS_NO_AVX2").is_err()
            && is_x86_feature_detected!("avx512f") && is_x86_feature_detected!("avx512bw") && is_x86_feature_detected!("avx512vl") }
        #[cfg(not(target_arch = "x86_64"))]
        { false }
    })
}

/// `shift_all` with sixteen lanes: the same sum (or difference) and wrap, the smaller kept.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn shift_all_avx512(r: &mut [u32], d: &[u32], p: &[u32], add: bool) {
    use std::arch::x86_64::*;
    let n = r.len();
    let mut i = 0;
    while i + 16 <= n {
        let rv = _mm512_loadu_si512(r.as_ptr().add(i) as *const _);
        let dv = _mm512_loadu_si512(d.as_ptr().add(i) as *const _);
        let pv = _mm512_loadu_si512(p.as_ptr().add(i) as *const _);
        let t = if add { _mm512_add_epi32(rv, dv) } else { _mm512_sub_epi32(rv, dv) };
        let w = if add { _mm512_sub_epi32(t, pv) } else { _mm512_add_epi32(t, pv) };
        _mm512_storeu_si512(r.as_mut_ptr().add(i) as *mut _, _mm512_min_epu32(t, w));
        i += 16;
    }
    shift_all_body(&mut r[i..], &d[i..], &p[i..], add);
}

/// v mod m without allocating (num-bigint's `%` builds a BigUint for the result).
#[inline]
fn rem_u32(v: &BigUint, m: u32) -> u32 {
    v.iter_u32_digits().rev().fold(0u64, |r, d| ((r << 32) | d as u64) % m as u64) as u32
}
/// Divide the factor-base primes in `js` out of `v`: the multiplicity of each from remainders
/// (no allocation), then one division per 64-bit chunk of their product.
#[inline]
fn divide_out(v: &mut BigUint, js: &[u32], base: &[Prime], exps: &mut Vec<(u32, u32)>) {
    let mut acc: u64 = 1;
    for &j in js {
        let p = base[j as usize].p;
        if p == 2 { continue; }
        let p32 = p as u32;
        let mut e = 0u32;
        let mut pk = p32 as u64;
        while pk <= u32::MAX as u64 && rem_u32(v, pk as u32) == 0 { e += 1; pk *= p; }
        if pk > u32::MAX as u64 {
            // a high power: finish with the slow path
            let mut q = p32 as u64; for _ in 1..e { q *= p; }
            if acc > 1 { *v /= acc; acc = 1; }   // (the first version dropped `acc` here and lost relations)
            *v /= q;
            while (&*v % p).is_zero() { *v /= p; e += 1; }
            for _ in 0..e { bump_exp(exps, j + 1); }
            continue;
        }
        if e == 0 { continue; }
        let pe = pk / p;
        for _ in 0..e { bump_exp(exps, j + 1); }
        match acc.checked_mul(pe) { Some(a) => acc = a, None => { *v /= acc; acc = pe; } }
    }
    if acc > 1 { *v /= acc; }
}
fn bump_exp(exps: &mut Vec<(u32, u32)>, c: u32) { match exps.iter_mut().find(|(d, _)| *d == c) { Some(t) => t.1 += 1, None => exps.push((c, 1)) } }

/// The sink bytes past the interval that a stride falling outside it adds to, so the add loops
/// have no branch. One sink byte made every miss wait for the previous miss's store (a serial
/// read-modify-write chain through memory, about 6 cycles a miss, and 30 % of the mid-band
/// strides miss); `SINK` bytes chosen by the low bits of the position break the chain.
/// `SIQS_SINK_MASK=0` restores the single byte, for A/B runs.
pub const SINK: usize = 64;
pub fn sink_mask() -> u32 {
    static V: std::sync::OnceLock<u32> = std::sync::OnceLock::new();
    *V.get_or_init(|| std::env::var("SIQS_SINK_MASK").ok().and_then(|v| v.parse().ok()).unwrap_or(SINK as u32 - 1) & (SINK as u32 - 1))
}

/// acc[pos] += log for every root below `width`; a root at or past it adds to one of the sink
/// bytes acc[width + (pos & mask)], so the loop has no branch. `acc.len() == width + SINK`.
#[inline(never)]
fn direct_add(acc: &mut [u8], roots: &[u32], logs: &[u8], width: u32, mask: u32) {
    assert!(acc.len() == width as usize + SINK && roots.len() == logs.len() && (mask as usize) < SINK);
    let a = acc.as_mut_ptr();
    for (&pos, &lg) in roots.iter().zip(logs) {
        let idx = if pos < width { pos } else { width + (pos & mask) } as usize;
        // SAFETY: idx < width + SINK = acc.len()
        unsafe { *a.add(idx) = (*a.add(idx)).wrapping_add(lg); }
    }
}

/// `direct_add` for primes that land at most `k` times per root: k fixed strides each, the
/// ones past the interval into the sink byte.
#[inline(never)]
fn direct_add_k(acc: &mut [u8], roots: &[u32], primes: &[u32], logs: &[u8], width: u32, k: u32, mask: u32) {
    assert!(acc.len() == width as usize + SINK && roots.len() == logs.len() && roots.len() == primes.len() && (mask as usize) < SINK);
    let a = acc.as_mut_ptr();
    for ((&r, &p), &lg) in roots.iter().zip(primes).zip(logs) {
        let mut pos = r;
        for _ in 0..k {
            let idx = if pos < width { pos } else { width + (pos & mask) } as usize;
            // SAFETY: idx < width + SINK = acc.len()
            unsafe { *a.add(idx) = (*a.add(idx)).wrapping_add(lg); }
            pos = pos.wrapping_add(p);
        }
    }
}

/// The wide primes' hits of one chunk, compacted: each in-range root as one word
/// `pos | i << POS_BITS` (`i` the index within the chunk of at most 2^13 primes). Under AVX2
/// eight roots are compared at once and the in-range lanes permuted to the front with a
/// 256-entry lane table (YAFU's `creg` kernel does the same), so a miss costs no store and no
/// read-modify-write of a sink byte: the adds then run over the hits alone. At 100 digits
/// about 117,000 of the 130,000 wide roots per polynomial miss, and each miss used to be a
/// dependent read-modify-write on the sink byte (about 1 ns).
const POS_BITS: u32 = 19;
#[cfg_attr(not(target_arch = "x86_64"), allow(dead_code))]
static LANE_TABLE: [[u32; 8]; 256] = {
    let mut t = [[0u32; 8]; 256];
    let mut m = 0;
    while m < 256 {
        let (mut j, mut b) = (0, 0);
        while b < 8 { if m >> b & 1 == 1 { t[m][j] = b as u32; j += 1; } b += 1; }
        m += 1;
    }
    t
};
/// WebAssembly SIMD uses four lanes for the kernels that use eight or sixteen on x86.
/// Enabled with `-C target-feature=+simd128`; the runtime must support SIMD128.
#[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
mod wasm_simd {
    use core::arch::wasm32::*;
    use super::POS_BITS;
    /// Byte shuffles that move the set lanes of a 4-bit mask to the front (the 4-lane LANE_TABLE).
    static SWZ: [[u8; 16]; 16] = {
        let mut t = [[0u8; 16]; 16];
        let mut m = 0;
        while m < 16 {
            let (mut j, mut b) = (0, 0);
            while b < 4 {
                if m >> b & 1 == 1 { let mut i = 0; while i < 4 { t[m][4 * j + i] = (4 * b + i) as u8; i += 1; } j += 1; }
                b += 1;
            }
            m += 1;
        }
        t
    };
    #[inline(always)]
    unsafe fn lane_tag() -> v128 { u32x4(0, 1 << POS_BITS, 2 << POS_BITS, 3 << POS_BITS) }
    /// `compact_hits` with four lanes; the same entries in the same order as the scalar loop.
    #[inline(always)]
    pub unsafe fn compact(roots: &[u32], width: u32, out: *mut u32) -> usize {
        let w = u32x4_splat(width);
        let step = u32x4_splat(4 << POS_BITS);
        let mut tag = lane_tag();
        let (mut cnt, mut k) = (0usize, 0usize);
        while k + 4 <= roots.len() {
            let v = v128_load(roots.as_ptr().add(k) as *const v128);
            let m = i32x4_bitmask(u32x4_lt(v, w)) as usize;
            let sh = v128_load(SWZ[m].as_ptr() as *const v128);
            v128_store(out.add(cnt) as *mut v128, i8x16_swizzle(v128_or(v, tag), sh));
            cnt += m.count_ones() as usize;
            tag = u32x4_add(tag, step);
            k += 4;
        }
        while k < roots.len() { let pos = *roots.get_unchecked(k); if pos < width { *out.add(cnt) = pos | (k as u32) << POS_BITS; cnt += 1; } k += 1; }
        cnt
    }
    /// `compact_hits_k` with four lanes (a group's strides stride-major, as the x86 kernels).
    #[inline(always)]
    pub unsafe fn compact_k(roots: &[u32], primes: &[u32], kk: u32, width: u32, out: *mut u32) -> usize {
        let w = u32x4_splat(width);
        let step = u32x4_splat(4 << POS_BITS);
        let mut tag = lane_tag();
        let (mut cnt, mut i) = (0usize, 0usize);
        while i + 4 <= roots.len() {
            let mut v = v128_load(roots.as_ptr().add(i) as *const v128);
            let p = v128_load(primes.as_ptr().add(i) as *const v128);
            for _ in 0..kk {
                let m = i32x4_bitmask(u32x4_lt(v, w)) as usize;
                let sh = v128_load(SWZ[m].as_ptr() as *const v128);
                v128_store(out.add(cnt) as *mut v128, i8x16_swizzle(v128_or(v, tag), sh));
                cnt += m.count_ones() as usize;
                v = u32x4_add(v, p);
            }
            tag = u32x4_add(tag, step);
            i += 4;
        }
        while i < roots.len() {
            let mut pos = *roots.get_unchecked(i);
            let p = *primes.get_unchecked(i);
            for _ in 0..kk { if pos < width { *out.add(cnt) = pos | (i as u32) << POS_BITS; cnt += 1; } pos = pos.wrapping_add(p); }
            i += 1;
        }
        cnt
    }
    /// `divides_flags` with four lanes: |x − r|·inv ≤ lim for either root, eight primes a byte.
    #[inline(always)]
    pub unsafe fn divides(x: u32, r0: &[u32], r1: &[u32], inv: &[u32], lim: &[u32], flags: &mut [u8]) {
        let n = r0.len();
        let xv = u32x4_splat(x);
        let mut j = 0usize;
        let mut w = 0usize;
        while j + 8 <= n {
            let mut bits = 0u8;
            for half in 0..2 {
                let o = j + 4 * half;
                let iv = v128_load(inv.as_ptr().add(o) as *const v128);
                let lv = v128_load(lim.as_ptr().add(o) as *const v128);
                let a = v128_load(r0.as_ptr().add(o) as *const v128);
                let b = v128_load(r1.as_ptr().add(o) as *const v128);
                let da = u32x4_sub(u32x4_max(xv, a), u32x4_min(xv, a));
                let db = u32x4_sub(u32x4_max(xv, b), u32x4_min(xv, b));
                let ha = u32x4_le(i32x4_mul(da, iv), lv);
                let hb = u32x4_le(i32x4_mul(db, iv), lv);
                bits |= (i32x4_bitmask(v128_or(ha, hb)) as u8) << (4 * half);
            }
            *flags.get_unchecked_mut(w) = bits;
            w += 1;
            j += 8;
        }
        // the tail byte, scalar
        if j < n {
            let mut bits = 0u8;
            for t in 0..8 {
                let q = j + t;
                if q >= n { break; }
                let hit = x.abs_diff(r0[q]).wrapping_mul(inv[q]) <= lim[q] || x.abs_diff(r1[q]).wrapping_mul(inv[q]) <= lim[q];
                bits |= (hit as u8) << t;
            }
            flags[w] = bits;
        }
    }
    /// The high-bit scan sixteen bytes at a time (`i8x16.bitmask`), every position in order,
    /// tails scalar; the counterpart of `scan_block_avx2`.
    #[inline(always)]
    pub unsafe fn scan(block: &[u8], offset: usize, hits: &mut Vec<usize>) {
        let mut i = 0;
        while block.len() - i >= 16 {
            let v = v128_load(block.as_ptr().add(i) as *const v128);
            let mut m = i8x16_bitmask(v) as u32;
            while m != 0 { let b = m.trailing_zeros() as usize; hits.push(offset + i + b); m &= m - 1; }
            i += 16;
        }
        while i < block.len() { if *block.get_unchecked(i) & 0x80 != 0 { hits.push(offset + i); } i += 1; }
    }
}

#[inline(always)]
fn compact_hits(roots: &[u32], width: u32, hits: &mut Vec<u32>, avx2: bool) {
    debug_assert!(roots.len() <= 1 << (32 - POS_BITS) && width < 1 << POS_BITS);
    let n = roots.len();
    hits.reserve(n + 16);
    let base_len = hits.len();
    #[cfg(target_arch = "x86_64")]
    if avx2 && avx512() {
        // SAFETY: avx512f detected; `n + 16` slots of room for the last full-width store
        unsafe {
            let cnt = compact_hits_avx512(roots, width, hits.as_mut_ptr().add(base_len));
            hits.set_len(base_len + cnt);
        }
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if avx2 {
        // SAFETY: `avx2` was set from `is_x86_feature_detected!("avx2")`; `hits` has `n + 8`
        // slots of room so the full-width store of the last group stays in bounds
        unsafe {
            let cnt = compact_hits_avx2(roots, width, hits.as_mut_ptr().add(base_len));
            hits.set_len(base_len + cnt);
        }
        return;
    }
    let _ = (avx2, base_len);
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    {
        // SAFETY: `n + 16` slots of room for the last full-width store
        unsafe { let cnt = wasm_simd::compact(roots, width, hits.as_mut_ptr().add(base_len)); hits.set_len(base_len + cnt); }
        return;
    }
    #[allow(unreachable_code)]
    for (i, &pos) in roots.iter().enumerate() { if pos < width { hits.push(pos | (i as u32) << POS_BITS); } }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn compact_hits_avx2(roots: &[u32], width: u32, out: *mut u32) -> usize {
    use std::arch::x86_64::*;
    let wm1 = _mm256_set1_epi32(width.wrapping_sub(1) as i32);
    let lane_step = _mm256_set1_epi32(8 << POS_BITS);
    let mut tag = _mm256_slli_epi32(_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7), POS_BITS as i32);
    let mut cnt = 0usize;
    let mut k = 0;
    while k + 8 <= roots.len() {
        let v = _mm256_loadu_si256(roots.as_ptr().add(k) as *const __m256i);
        // pos < width, unsigned (the skipped primes' sentinel root u32::MAX − 1 is negative
        // as i32 and a signed compare took it for a hit): min(pos, width − 1) == pos
        let m = _mm256_cmpeq_epi32(_mm256_min_epu32(v, wm1), v);
        let bits = _mm256_movemask_ps(_mm256_castsi256_ps(m)) as usize;
        let perm = _mm256_loadu_si256(LANE_TABLE[bits].as_ptr() as *const __m256i);
        let packed = _mm256_or_si256(v, tag);
        _mm256_storeu_si256(out.add(cnt) as *mut __m256i, _mm256_permutevar8x32_epi32(packed, perm));
        cnt += bits.count_ones() as usize;
        tag = _mm256_add_epi32(tag, lane_step);
        k += 8;
    }
    while k < roots.len() { let pos = *roots.get_unchecked(k); if pos < width { *out.add(cnt) = pos | (k as u32) << POS_BITS; cnt += 1; } k += 1; }
    cnt
}
/// `compact_hits_avx2` with sixteen lanes: the in-range lanes are compressed to the front of a
/// register (`vpcompressd`) and the whole register stored, so a group costs one compare, one
/// compress and one store, no table lookup. Same entries in the same order as the 8-lane loop.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn compact_hits_avx512(roots: &[u32], width: u32, out: *mut u32) -> usize {
    use std::arch::x86_64::*;
    let wm1 = _mm512_set1_epi32(width.wrapping_sub(1) as i32);
    let lane_step = _mm512_set1_epi32(16 << POS_BITS);
    let mut tag = _mm512_slli_epi32::<POS_BITS>(_mm512_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));
    let mut cnt = 0usize;
    let mut k = 0;
    while k + 16 <= roots.len() {
        let v = _mm512_loadu_si512(roots.as_ptr().add(k) as *const _);
        let m = _mm512_cmple_epu32_mask(v, wm1);   // pos < width, unsigned
        _mm512_storeu_si512(out.add(cnt) as *mut _, _mm512_maskz_compress_epi32(m, _mm512_or_si512(v, tag)));
        cnt += m.count_ones() as usize;
        tag = _mm512_add_epi32(tag, lane_step);
        k += 16;
    }
    while k < roots.len() { let pos = *roots.get_unchecked(k); if pos < width { *out.add(cnt) = pos | (k as u32) << POS_BITS; cnt += 1; } k += 1; }
    cnt
}
/// `compact_hits_k_avx2` with sixteen lanes, the strides past the interval masked out of the
/// compress and the position advance.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn compact_hits_k_avx512(roots: &[u32], primes: &[u32], k: u32, width: u32, out: *mut u32) -> usize {
    use std::arch::x86_64::*;
    let wm1 = _mm512_set1_epi32(width.wrapping_sub(1) as i32);
    let lane_step = _mm512_set1_epi32(16 << POS_BITS);
    let mut tag = _mm512_slli_epi32::<POS_BITS>(_mm512_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));
    let mut cnt = 0usize;
    let mut i = 0;
    while i + 16 <= roots.len() {
        let mut v = _mm512_loadu_si512(roots.as_ptr().add(i) as *const _);
        let p = _mm512_loadu_si512(primes.as_ptr().add(i) as *const _);
        for _ in 0..k {
            let m = _mm512_cmple_epu32_mask(v, wm1);
            _mm512_storeu_si512(out.add(cnt) as *mut _, _mm512_maskz_compress_epi32(m, _mm512_or_si512(v, tag)));
            cnt += m.count_ones() as usize;
            v = _mm512_add_epi32(v, p);
        }
        tag = _mm512_add_epi32(tag, lane_step);
        i += 16;
    }
    while i < roots.len() {
        let mut pos = *roots.get_unchecked(i);
        let p = *primes.get_unchecked(i);
        for _ in 0..k { if pos < width { *out.add(cnt) = pos | (i as u32) << POS_BITS; cnt += 1; } pos = pos.wrapping_add(p); }
        i += 1;
    }
    cnt
}
/// `mid_block_avx2` with sixteen lanes: roots below `hi` compressed out and advanced by their
/// prime, the rest left where they are, until no lane is inside the block.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
unsafe fn mid_block_avx512(roots: &mut [u32], primes: &[u32], hi: u32, kb: u32, out: *mut u32) -> usize {
    use std::arch::x86_64::*;
    let hm1 = _mm512_set1_epi32(hi.wrapping_sub(1) as i32);
    let lane_step = _mm512_set1_epi32(16 << POS_BITS);
    let mut tag = _mm512_slli_epi32::<POS_BITS>(_mm512_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));
    let mut cnt = 0usize;
    let mut i = 0;
    while i + 16 <= roots.len() {
        let rp = roots.as_mut_ptr().add(i);
        let mut v = _mm512_loadu_si512(rp as *const _);
        let p = _mm512_loadu_si512(primes.as_ptr().add(i) as *const _);
        for _ in 0..kb {
            let m = _mm512_cmple_epu32_mask(v, hm1);   // v < hi, unsigned
            if m == 0 { break; }
            _mm512_storeu_si512(out.add(cnt) as *mut _, _mm512_maskz_compress_epi32(m, _mm512_or_si512(v, tag)));
            cnt += m.count_ones() as usize;
            v = _mm512_mask_add_epi32(v, m, v, p);
        }
        _mm512_storeu_si512(rp as *mut _, v);
        tag = _mm512_add_epi32(tag, lane_step);
        i += 16;
    }
    while i < roots.len() {
        let r = roots.get_unchecked_mut(i);
        let p = *primes.get_unchecked(i);
        while *r < hi { *out.add(cnt) = *r | (i as u32) << POS_BITS; cnt += 1; *r += p; }
        i += 1;
    }
    cnt
}

/// `compact_hits` for a band whose primes land up to `k` times per root: stride `j` of every
/// root is `r + j·p`, compared and compacted like the single stride. Scalar fallback included.
#[inline(always)]
fn compact_hits_k(roots: &[u32], primes: &[u32], k: u32, width: u32, hits: &mut Vec<u32>, avx2: bool) {
    debug_assert!(roots.len() == primes.len() && roots.len() <= 1 << (32 - POS_BITS) && width < 1 << POS_BITS);
    let n = roots.len();
    hits.reserve(n * k as usize + 16);
    let base_len = hits.len();
    #[cfg(target_arch = "x86_64")]
    if avx2 && avx512() {
        // SAFETY: avx512f detected; room for `n·k + 16` entries
        unsafe {
            let cnt = compact_hits_k_avx512(roots, primes, k, width, hits.as_mut_ptr().add(base_len));
            hits.set_len(base_len + cnt);
        }
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if avx2 {
        // SAFETY: as for `compact_hits`; room for `n·k + 8` entries
        unsafe {
            let cnt = compact_hits_k_avx2(roots, primes, k, width, hits.as_mut_ptr().add(base_len));
            hits.set_len(base_len + cnt);
        }
        return;
    }
    let _ = (avx2, base_len);
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    {
        // SAFETY: `n·k + 16` slots of room
        unsafe { let cnt = wasm_simd::compact_k(roots, primes, k, width, hits.as_mut_ptr().add(base_len)); hits.set_len(base_len + cnt); }
        return;
    }
    #[allow(unreachable_code)]
    for (i, (&r, &p)) in roots.iter().zip(primes).enumerate() {
        let mut pos = r;
        for _ in 0..k { if pos < width { hits.push(pos | (i as u32) << POS_BITS); } pos = pos.wrapping_add(p); }
    }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn compact_hits_k_avx2(roots: &[u32], primes: &[u32], k: u32, width: u32, out: *mut u32) -> usize {
    use std::arch::x86_64::*;
    let wm1 = _mm256_set1_epi32(width.wrapping_sub(1) as i32);
    let lane_step = _mm256_set1_epi32(8 << POS_BITS);
    let mut tag = _mm256_slli_epi32(_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7), POS_BITS as i32);
    let mut cnt = 0usize;
    let mut i = 0;
    while i + 8 <= roots.len() {
        let mut v = _mm256_loadu_si256(roots.as_ptr().add(i) as *const __m256i);
        let p = _mm256_loadu_si256(primes.as_ptr().add(i) as *const __m256i);
        for _ in 0..k {
            let m = _mm256_cmpeq_epi32(_mm256_min_epu32(v, wm1), v);
            let bits = _mm256_movemask_ps(_mm256_castsi256_ps(m)) as usize;
            let perm = _mm256_loadu_si256(LANE_TABLE[bits].as_ptr() as *const __m256i);
            _mm256_storeu_si256(out.add(cnt) as *mut __m256i, _mm256_permutevar8x32_epi32(_mm256_or_si256(v, tag), perm));
            cnt += bits.count_ones() as usize;
            v = _mm256_add_epi32(v, p);
        }
        tag = _mm256_add_epi32(tag, lane_step);
        i += 8;
    }
    while i < roots.len() {
        let mut pos = *roots.get_unchecked(i);
        let p = *primes.get_unchecked(i);
        for _ in 0..k { if pos < width { *out.add(cnt) = pos | (i as u32) << POS_BITS; cnt += 1; } pos = pos.wrapping_add(p); }
        i += 1;
    }
    cnt
}
/// One 32 KB block of the mid primes (YAFU's "sparse" strategy): each root is advanced
/// through the block, every position below `hi` compacted as a hit word, and the root left at
/// its first position at or past `hi` for the next block. `kb` bounds the strides per root in a
/// block (ceil(block / smallest prime of the slice)). The adds then land in one block, in L1,
/// instead of anywhere in a 256 KB interval.
#[inline(always)]
fn mid_block(roots: &mut [u32], primes: &[u32], hi: u32, kb: u32, hits: &mut Vec<u32>, avx2: bool) {
    debug_assert!(roots.len() == primes.len() && roots.len() <= 1 << (32 - POS_BITS) && hi <= 1 << POS_BITS);
    let n = roots.len();
    hits.reserve(n * kb as usize + 16);
    let base_len = hits.len();
    #[cfg(target_arch = "x86_64")]
    if avx2 && avx512() {
        // SAFETY: avx512f detected; room for `n·kb + 16` entries
        unsafe {
            let cnt = mid_block_avx512(roots, primes, hi, kb, hits.as_mut_ptr().add(base_len));
            hits.set_len(base_len + cnt);
        }
        return;
    }
    #[cfg(target_arch = "x86_64")]
    if avx2 {
        // SAFETY: as for `compact_hits`; room for `n·kb + 8` entries
        unsafe {
            let cnt = mid_block_avx2(roots, primes, hi, kb, hits.as_mut_ptr().add(base_len));
            hits.set_len(base_len + cnt);
        }
        return;
    }
    let _ = (avx2, base_len);
    for (i, (r, &p)) in roots.iter_mut().zip(primes).enumerate() {
        while *r < hi { hits.push(*r | (i as u32) << POS_BITS); *r += p; }
    }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn mid_block_avx2(roots: &mut [u32], primes: &[u32], hi: u32, kb: u32, out: *mut u32) -> usize {
    use std::arch::x86_64::*;
    let hm1 = _mm256_set1_epi32(hi.wrapping_sub(1) as i32);
    let lane_step = _mm256_set1_epi32(8 << POS_BITS);
    let mut tag = _mm256_slli_epi32(_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7), POS_BITS as i32);
    let mut cnt = 0usize;
    let mut i = 0;
    while i + 8 <= roots.len() {
        let rp = roots.as_mut_ptr().add(i) as *mut __m256i;
        let mut v = _mm256_loadu_si256(rp);
        let p = _mm256_loadu_si256(primes.as_ptr().add(i) as *const __m256i);
        for _ in 0..kb {
            let m = _mm256_cmpeq_epi32(_mm256_min_epu32(v, hm1), v);   // v < hi, unsigned
            let bits = _mm256_movemask_ps(_mm256_castsi256_ps(m)) as usize;
            if bits == 0 { break; }
            let perm = _mm256_loadu_si256(LANE_TABLE[bits].as_ptr() as *const __m256i);
            _mm256_storeu_si256(out.add(cnt) as *mut __m256i, _mm256_permutevar8x32_epi32(_mm256_or_si256(v, tag), perm));
            cnt += bits.count_ones() as usize;
            v = _mm256_add_epi32(v, _mm256_and_si256(p, m));
        }
        _mm256_storeu_si256(rp, v);
        tag = _mm256_add_epi32(tag, lane_step);
        i += 8;
    }
    while i < roots.len() {
        let r = roots.get_unchecked_mut(i);
        let p = *primes.get_unchecked(i);
        while *r < hi { *out.add(cnt) = *r | (i as u32) << POS_BITS; cnt += 1; *r += p; }
        i += 1;
    }
    cnt
}
/// The recorded hits of one chunk whose position is a candidate (one bit each in `hitbits`),
/// reported as (position, prime index `idx0 + lane`).
#[cfg_attr(target_arch = "wasm32", inline(never))]
fn tag_bitmap(run: &[u32], idx0: u32, hitbits: &[u64], f: &mut impl FnMut(usize, u32)) {
    let pmask = (1u32 << POS_BITS) - 1;
    for &e in run {
        let pos = (e & pmask) as usize;
        if hitbits[pos >> 6] >> (pos & 63) & 1 == 1 { f(pos, idx0 + (e >> POS_BITS)); }
    }
}

/// Add each recorded hit's log at its position.
#[inline(never)]
fn add_hits(acc: &mut [u8], hits: &[u32], logs: &[u8]) {
    let mask = (1u32 << POS_BITS) - 1;
    for &e in hits {
        let pos = (e & mask) as usize;
        let i = (e >> POS_BITS) as usize;
        acc[pos] = acc[pos].wrapping_add(logs[i]);
    }
}

#[inline(always)]
#[cfg(target_arch = "x86_64")]
fn append_scan_mask(mut mask: u32, offset: usize, hits: &mut Vec<usize>) {
    while mask != 0 {
        hits.push(offset + mask.trailing_zeros() as usize);
        mask &= mask - 1;
    }
}

/// Adapt YAFU's tdiv_scan.c bulk high-bit scan (64/128-byte variants), decoding
/// masks only when a high bit is set. Preserve every byte position and order,
/// including the endpoints and arbitrary tails of our [-M, M] interval.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn scan_block_avx2(block: &[u8], offset: usize, hits: &mut Vec<usize>) {
    use std::arch::x86_64::*;
    let mut i = 0;
    while block.len() - i >= 128 {
        let a = _mm256_loadu_si256(block.as_ptr().add(i).cast());
        let b = _mm256_loadu_si256(block.as_ptr().add(i + 32).cast());
        let c = _mm256_loadu_si256(block.as_ptr().add(i + 64).cast());
        let d = _mm256_loadu_si256(block.as_ptr().add(i + 96).cast());
        let any = _mm256_or_si256(_mm256_or_si256(a, b), _mm256_or_si256(c, d));
        if _mm256_movemask_epi8(any) != 0 {
            append_scan_mask(_mm256_movemask_epi8(a) as u32, offset + i, hits);
            append_scan_mask(_mm256_movemask_epi8(b) as u32, offset + i + 32, hits);
            append_scan_mask(_mm256_movemask_epi8(c) as u32, offset + i + 64, hits);
            append_scan_mask(_mm256_movemask_epi8(d) as u32, offset + i + 96, hits);
        }
        i += 128;
    }
    while block.len() - i >= 32 {
        let v = _mm256_loadu_si256(block.as_ptr().add(i).cast());
        append_scan_mask(_mm256_movemask_epi8(v) as u32, offset + i, hits);
        i += 32;
    }
    for (j, &v) in block[i..].iter().enumerate() {
        if v & 0x80 != 0 { hits.push(offset + i + j); }
    }
}

/// One bit per prime: does p divide (x − r0) or (x − r1)? By the multiply-and-compare test
/// (`Prime::divides`), eight primes a lane group under AVX2. `flags.len() == ceil(n / 8)`.
#[inline(always)]
fn divides_flags(x: u32, r0: &[u32], r1: &[u32], inv: &[u32], lim: &[u32], flags: &mut [u8], avx2: bool) {
    let n = r0.len();
    debug_assert!(r1.len() == n && inv.len() == n && lim.len() == n && flags.len() == n.div_ceil(8));
    #[cfg(target_arch = "x86_64")]
    if avx2 {
        // SAFETY: `avx2` was set from `is_x86_feature_detected!("avx2")`; loads are within the slices
        unsafe { return divides_flags_avx2(x, r0, r1, inv, lim, flags); }
    }
    let _ = avx2;
    #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
    {
        // SAFETY: all loads are within the slices (whole groups of eight, the tail scalar)
        unsafe { wasm_simd::divides(x, r0, r1, inv, lim, flags); }
        return;
    }
    #[allow(unreachable_code)]
    for (w, f) in flags.iter_mut().enumerate() {
        let mut bits = 0u8;
        for t in 0..8 {
            let j = w * 8 + t;
            if j >= n { break; }
            let hit = x.abs_diff(r0[j]).wrapping_mul(inv[j]) <= lim[j] || x.abs_diff(r1[j]).wrapping_mul(inv[j]) <= lim[j];
            bits |= (hit as u8) << t;
        }
        *f = bits;
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn divides_flags_avx2(x: u32, r0: &[u32], r1: &[u32], inv: &[u32], lim: &[u32], flags: &mut [u8]) {
    use std::arch::x86_64::*;
    let n = r0.len();
    let xv = _mm256_set1_epi32(x as i32);
    let mut j = 0;
    let mut w = 0;
    while j + 8 <= n {
        let a = _mm256_loadu_si256(r0.as_ptr().add(j) as *const __m256i);
        let b = _mm256_loadu_si256(r1.as_ptr().add(j) as *const __m256i);
        let iv = _mm256_loadu_si256(inv.as_ptr().add(j) as *const __m256i);
        let lm = _mm256_loadu_si256(lim.as_ptr().add(j) as *const __m256i);
        // |x − r| = max(x, r) − min(x, r), unsigned
        let da = _mm256_sub_epi32(_mm256_max_epu32(xv, a), _mm256_min_epu32(xv, a));
        let db = _mm256_sub_epi32(_mm256_max_epu32(xv, b), _mm256_min_epu32(xv, b));
        let pa = _mm256_mullo_epi32(da, iv);
        let pb = _mm256_mullo_epi32(db, iv);
        // p ≤ lim  ⇔  max(p, lim) == lim, unsigned
        let ma = _mm256_cmpeq_epi32(_mm256_max_epu32(pa, lm), lm);
        let mb = _mm256_cmpeq_epi32(_mm256_max_epu32(pb, lm), lm);
        let m = _mm256_or_si256(ma, mb);
        flags[w] = _mm256_movemask_ps(_mm256_castsi256_ps(m)) as u8;
        j += 8;
        w += 1;
    }
    if j < n {
        let mut bits = 0u8;
        for t in 0..(n - j) {
            let k = j + t;
            let hit = x.abs_diff(r0[k]).wrapping_mul(inv[k]) <= lim[k] || x.abs_diff(r1[k]).wrapping_mul(inv[k]) <= lim[k];
            bits |= (hit as u8) << t;
        }
        flags[w] = bits;
    }
}

/// Call `f` on every bucket entry whose position (entry & seg_mask) is one of the up-to-eight
/// survivor positions in `surv` (unused slots are u32::MAX, which no position equals).
#[inline(always)]
fn run_hits(run: &[u32], seg_mask: u32, surv: &[u32; 8], avx2: bool, mut f: impl FnMut(usize, u32)) {
    #[cfg(target_arch = "x86_64")]
    if avx2 {
        // SAFETY: `avx2` was set from `is_x86_feature_detected!("avx2")`; loads are unaligned
        // and within `run`.
        unsafe { return run_hits_avx2(run, seg_mask, surv, &mut f); }
    }
    let _ = avx2;
    for (k, &e) in run.iter().enumerate() { let c = e & seg_mask; if surv.contains(&c) { f(k, e); } }
}
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn run_hits_avx2(run: &[u32], seg_mask: u32, surv: &[u32; 8], f: &mut impl FnMut(usize, u32)) {
    use std::arch::x86_64::*;
    let mask = _mm256_set1_epi32(seg_mask as i32);
    let s0 = _mm256_set1_epi32(surv[0] as i32);
    let s1 = _mm256_set1_epi32(surv[1] as i32);
    let s2 = _mm256_set1_epi32(surv[2] as i32);
    let s3 = _mm256_set1_epi32(surv[3] as i32);
    let s4 = _mm256_set1_epi32(surv[4] as i32);
    let s5 = _mm256_set1_epi32(surv[5] as i32);
    let s6 = _mm256_set1_epi32(surv[6] as i32);
    let s7 = _mm256_set1_epi32(surv[7] as i32);
    let four = surv[4] == u32::MAX;   // unused slots are u32::MAX, which no position equals
    let mut k = 0;
    while k + 8 <= run.len() {
        let v = _mm256_and_si256(_mm256_loadu_si256(run.as_ptr().add(k) as *const __m256i), mask);
        let mut m = _mm256_or_si256(_mm256_or_si256(_mm256_cmpeq_epi32(v, s0), _mm256_cmpeq_epi32(v, s1)), _mm256_or_si256(_mm256_cmpeq_epi32(v, s2), _mm256_cmpeq_epi32(v, s3)));
        if !four { m = _mm256_or_si256(m, _mm256_or_si256(_mm256_or_si256(_mm256_cmpeq_epi32(v, s4), _mm256_cmpeq_epi32(v, s5)), _mm256_or_si256(_mm256_cmpeq_epi32(v, s6), _mm256_cmpeq_epi32(v, s7)))); }
        let mut bits = _mm256_movemask_ps(_mm256_castsi256_ps(m)) as u32;
        while bits != 0 { let j = bits.trailing_zeros() as usize; f(k + j, run[k + j]); bits &= bits - 1; }
        k += 8;
    }
    while k < run.len() { let c = run[k] & seg_mask; if surv.contains(&c) { f(k, run[k]); } k += 1; }
}


/// Add `lg` at every position of two strides of `p` inside the block `[lo, hi)`, each stride
/// starting at its absolute `next` position, and advance both to their first position at or
/// past `hi`. A root already past the block is left alone (an earlier implementation
/// version clamped it to the block end and lost it). Equal roots, a prime dividing k·n, are
/// one stride.
#[inline]
fn fill_block(block: &mut [u8], lo: usize, hi: usize, p: usize, lg: u8, next: &mut [u32; 2]) {
    let len = hi - lo;
    let block = &mut block[..len];
    let (n0, n1) = (next[0] as usize, next[1] as usize);
    if n0 == n1 {
        if n0 < hi { let mut c = n0 - lo; while c < len { block[c] = block[c].wrapping_add(lg); c += p; } next[0] = (lo + c) as u32; next[1] = next[0]; }
        return;
    }
    let (mut c0, mut c1) = (if n0 < hi { n0 - lo } else { len }, if n1 < hi { n1 - lo } else { len });
    // the two roots ordered, `c1` the later one: while it is in range both are (YAFU's medium
    // sieve does the same, with the bulk of the strides unrolled eight pairs at a time)
    let swapped = c0 > c1;
    if swapped { std::mem::swap(&mut c0, &mut c1); }
    let b = block.as_mut_ptr();
    // SAFETY: every index below is checked against `len` before use: in the unrolled loop
    // c1 + 7p < len bounds all sixteen, in the pair loop c0 ≤ c1 < len
    unsafe {
        let step8 = 7 * p;
        while c1 + step8 < len {
            let (mut a0, mut a1) = (c0, c1);
            for _ in 0..8 {
                *b.add(a0) = (*b.add(a0)).wrapping_add(lg);
                *b.add(a1) = (*b.add(a1)).wrapping_add(lg);
                a0 += p; a1 += p;
            }
            c0 += 8 * p; c1 += 8 * p;
        }
        while c1 < len {
            *b.add(c0) = (*b.add(c0)).wrapping_add(lg);
            *b.add(c1) = (*b.add(c1)).wrapping_add(lg);
            c0 += p; c1 += p;
        }
        if c0 < len { *b.add(c0) = (*b.add(c0)).wrapping_add(lg); c0 += p; }
    }
    if swapped { std::mem::swap(&mut c0, &mut c1); }
    if n0 < hi { next[0] = (lo + c0) as u32; }
    if n1 < hi { next[1] = (lo + c1) as u32; }
}

/// A cofactor between the large-prime bound and its square: two primes each below the bound,
/// or nothing. Below pmax² it is prime and too big; above, it must be composite and split.
fn split_cofactor(v: u128, bound: u128, pmax: u128) -> Option<Vec<(u64, u32)>> {
    if !double_large() || v < pmax * pmax || v > u64::MAX as u128 { return None; }
    let top = (bound as f64).powf(dlp_exp());
    if (v as f64) >= top || v >= 1u128 << 63 { return None; }
    let v = v as u64;
    count(7, 1);
    // a prime squared defeats every curve (both factors appear at once): take the root first
    let r = (v as f64).sqrt() as u64;
    for c in [r.saturating_sub(1), r, r + 1] { if c * c == v { return if (c as u128) < bound { count(9, 1); Some(vec![(c, 2)]) } else { None }; } }
    // a probable prime is too big for one large prime and cannot be two: drop it after one
    // Montgomery Fermat test (microseconds, not the twelve-base Miller-Rabin)
    if crate::arith64::fermat2(v) { count(8, 1); return None; }
    let t = Stamp::now();
    let mut rng = crate::arith64::Lcg(v ^ 0x9E37_79B9_7F4A_7C15);
    let a = crate::uecm::split(v, &mut rng);
    tick(6, t);
    let a = a?;
    count(9, 1);
    let b = v / a;
    if a as u128 >= bound || b as u128 >= bound { return None; }
    if a == b { Some(vec![(a, 2)]) } else { Some(vec![(a.min(b), 1), (a.max(b), 1)]) }
}

/// Sieve every b of one a. Returns the relations and the number of polynomials used.
fn sieve_family(n: &BigUint, base: &[Prime], fam: &Family, params: &Params, threshold_log: f32) -> (Vec<Relation>, usize) {
    let half = params.half as i64;
    let width = 2 * params.half + 1;
    let s = fam.qs.len();
    let pmax = base.last().unwrap().p;
    let large_bound = BigUint::from(pmax * params.large);
    let nb = BigInt::from(n.clone());
    let ab = BigInt::from(fam.a.clone());
    let y_a = if fam.q2 { &ab * 2 } else { ab.clone() };
    // start with b = Σ B_l
    let mut b: BigInt = fam.bs.iter().fold(BigInt::zero(), |acc, x| acc + BigInt::from(x.clone()));
    if fam.q2 && b.is_even() { b += &ab; }
    // Roots x ≡ (±t − b)/(step*a) and Gray-code deltas 2*B_l/(step*a),
    // where step is 1 for Q1 or 2 for Q2.
    let mut roots: Vec<[u32; 2]> = vec![[0; 2]; base.len()];
    // deltas[v][i]: for Gray bit v, the root shift 2·B_v·a⁻¹ mod p of every prime, contiguous
    // over the primes so a polynomial switch streams through one array (a large-input
    // review: the per-prime layout cost a cache miss a prime and a third of the time at a
    // 60,000-prime base)
    let word = fast_setup() && pmax <= u32::MAX as u64;
    let deltas: Vec<Vec<u32>> = if word {
        // 2·(B_l mod p)·a⁻¹ mod p with p < 2^32: the reduced double is below p, the product below 2^64
        fam.bs.iter().map(|bl| base.iter().enumerate().map(|(i, q)| ((2 * rem_u32(bl, q.p as u32) as u64 % q.p) * fam.ainv[i] % q.p) as u32).collect()).collect()
    } else {
        fam.bs.iter().map(|bl| base.iter().enumerate().map(|(i, q)| ((2u128 * (bl % q.p).to_u64().unwrap() as u128 * fam.ainv[i] as u128) % q.p as u128) as u32).collect()).collect()
    };
    let set_roots = |b: &BigInt, roots: &mut Vec<[u32; 2]>| {
        for (i, q) in base.iter().enumerate() {
            if fam.ainv[i] == 0 { continue; }
            if word {
                let bm = rem_floor_u32(b, q.p as u32) as u64;
                let r1 = (q.root + q.p - bm) % q.p * fam.ainv[i] % q.p;
                let r2 = (2 * q.p - q.root - bm) % q.p * fam.ainv[i] % q.p;
                roots[i] = [r1 as u32, r2 as u32];
                continue;
            }
            let bm = b.mod_floor(&BigInt::from(q.p)).to_u64().unwrap();
            let r1 = ((q.root + q.p - bm) % q.p) as u128 * fam.ainv[i] as u128 % q.p as u128;
            let r2 = ((2 * q.p - q.root - bm) % q.p) as u128 * fam.ainv[i] as u128 % q.p as u128;
            roots[i] = [r1 as u32, r2 as u32];
        }
    };
    set_roots(&b, &mut roots);
    let mut out = Vec::new();
    let mut acc = vec![0u8; width + SINK + 8];
    let sink = sink_mask();
    let half_mod: Vec<u32> = base.iter().map(|q| (params.half as u64 % q.p) as u32).collect();
    let mut starts: Vec<[u32; 2]> = vec![[0; 2]; base.len()];
    // Enable the direct large-prime stage when primes or the interval outgrow
    // the small-prime fill. The historical name bucket_on remains in this API;
    // this implementation uses direct root arrays, not stored sieve buckets.
    let bucket_on = bucket_sieve().unwrap_or(pmax as usize > width || width > 1 << 17);
    let seg = segment_len(bucket_on, base.len());
    // Below SIQS_DIRECT_FROM, fill cache-sized blocks by prime strides.
    // At and above that cutoff, add directly from the polynomial's root arrays.
    let direct_from: usize = std::env::var("SIQS_DIRECT_FROM").ok().and_then(|v| v.parse().ok()).unwrap_or(4096).min(seg);
    let first_large = if bucket_on { base.partition_point(|q| (q.p as usize) <= direct_from) } else { base.len() };
    // the Gray-code step from polynomial t−1 to t: which B_v flips, and whether it is added
    let gray_step = |t: usize| { let v = t.trailing_zeros() as usize; (v, (t >> (v + 1)) & 1 == 0) };
    // **the direct primes as arrays**: positions in the interval (start = root + half mod p),
    // primes, logs, the exact-division pair and per-switch deltas, as separate contiguous arrays
    // so a polynomial switch is one straight pass the compiler vectorises and the tagging tests
    // run eight lanes at a time. Primes that cannot hit (factors of a, or the rare prime dividing
    // kN with one root) get a sentinel position, zero deltas and a test that never passes.
    let nd = base.len() - first_large;
    let mut dr0: Vec<u32> = Vec::with_capacity(nd);
    let mut dr1: Vec<u32> = Vec::with_capacity(nd);
    let mut dp: Vec<u32> = Vec::with_capacity(nd);
    let mut dlog: Vec<u8> = Vec::with_capacity(nd);
    let mut dinv: Vec<u32> = Vec::with_capacity(nd);
    let mut dlim: Vec<u32> = Vec::with_capacity(nd);
    let mut dskip: Vec<bool> = Vec::with_capacity(nd);
    for i in first_large..base.len() {
        let q = &base[i];
        let skip = fam.ainv[i] == 0 || roots[i][0] == roots[i][1];
        dskip.push(skip);
        if skip { dr0.push(u32::MAX - 1); dr1.push(u32::MAX - 1); dp.push(u32::MAX); dlog.push(0); dinv.push(1); dlim.push(0); }
        else {
            dr0.push(shift_residue(roots[i][0], half_mod[i], q.p as u32, true));
            dr1.push(shift_residue(roots[i][1], half_mod[i], q.p as u32, true));
            dp.push(q.p as u32); dlog.push(q.log); dinv.push(q.inv); dlim.push(q.limit);
        }
    }
    let dd: Vec<Vec<u32>> = (0..s).map(|v| (first_large..base.len()).map(|i| if dskip[i - first_large] { 0 } else { deltas[v][i] }).collect()).collect();
    let avx2 = have_avx2();
    // The bulk mask scan won complete 70/80-digit comparisons at bases >=8192.
    // Smaller-base timings did not establish a gain; keep their existing scan.
    // SIQS_VECTOR_SCAN=0/1 overrides that size policy for controlled comparisons.
    let vector_scan = avx2 && std::env::var("SIQS_VECTOR_SCAN")
        .map(|v| v != "0").unwrap_or(base.len() >= 8192);
    // `SIQS_SHIFT_SIMD=0`: the root shifts by the baseline (SSE2) loop, for A/B runs
    let shift_simd = std::env::var("SIQS_SHIFT_SIMD").map(|v| v != "0").unwrap_or(true);
    // `SIQS_REC=0`: tag the wide primes by re-reading the root arrays instead of recording hits
    let rec = std::env::var("SIQS_REC").map(|v| v != "0").unwrap_or(true);
    // `SIQS_WIDE_TAG=bitmap`: tag the wide primes by one bitmap test per recorded hit instead of
    // `run_hits` per eight survivors. Measured 2026-11-13 at 100 digits, 8 threads: the bitmap
    // loop cost 87 thread-us a polynomial against 40 (a random 40 KB bitmap load per entry, where
    // `run_hits` streams the records against eight survivors in registers), so `runs` stays.
    // Off x86 (the browser) there is no vector `run_hits`, and the scalar one compares every
    // record with eight survivors; the bitmap test is the default there.
    let wide_tag_bitmap = std::env::var("SIQS_WIDE_TAG").map(|v| v == "bitmap").unwrap_or(cfg!(not(target_arch = "x86_64")));
    // `SIQS_REC_MID=1`: the same compaction for the mid bands (k strides per root)
    let rec_mid = std::env::var("SIQS_REC_MID").map(|v| v == "1").unwrap_or(false);
    // `SIQS_MID=block`: the mid primes sieved block by block (adds in L1), `SIQS_MID_BLOCK` the block size
    let mid_blockwise = std::env::var("SIQS_MID").map(|v| v == "block").unwrap_or(false);
    let mid_bs: usize = std::env::var("SIQS_MID_BLOCK").ok().and_then(|v| v.parse().ok()).unwrap_or(32768);
    // bands of the direct primes by strides per root: a prime in [width/2^k, width/2^(k-1))
    // lands at most 2^k times, so each band is a fixed-count loop with the misses sunk; the
    // primes wider than the interval (one stride) are the last band, index `first_wide` on
    let first_wide = base.partition_point(|q| (q.p as usize) < width).max(first_large);
    let fine_bands = std::env::var("SIQS_BANDS").map(|v| v != "pow2").unwrap_or(true);
    let mut bands: Vec<(usize, usize, u32)> = vec![(first_wide - first_large, nd, 1)];
    {
        let mut hi = first_wide;
        let mut k = 2u32;
        while hi > first_large {
            let pmin = width.div_ceil(k as usize);
            let lo_i = base.partition_point(|q| (q.p as usize) < pmin).max(first_large).min(hi);
            if lo_i < hi { bands.push((lo_i - first_large, hi - first_large, k)); }
            hi = lo_i;
            // the next band's stride count. Doubling (`SIQS_BANDS=pow2`, the rule until
            // 2026-11-13) gives every root of a band k strides although its prime may need only
            // k/2, so about 28 % of the strides fell outside the interval; steps of one up to 8
            // and of a quarter beyond leave about 8 %
            k = if fine_bands { if k < 8 { k + 1 } else { k + k / 4 } } else { k * 2 };
        }
    }
    let nmid = first_wide - first_large;
    let mut md_flags: Vec<u8> = vec![0; nmid.div_ceil(8)];   // one bit per mid prime
    // the primes of a that lie beyond the interval, checked by division on every candidate
    let large_a: Vec<usize> = fam.qs.iter().copied().filter(|&qi| qi >= first_large).collect();
    let polys = 1usize << (s - 1);
    let rmin = resieve_min();
    let band = base.partition_point(|q| q.p < 8192).min(first_large);
    let resieve_from = if rmin == 0 { first_large } else { base.partition_point(|q| q.p < rmin).min(first_large) };
    // which positions are candidates, one bit each (2-4 KB, stays in L1 across 16 threads where a
    // u32 index per position was 512 KB a thread and thrashed the shared cache); the candidate's
    // number is then a binary search in the short sorted `hits`
    let mut hitbits: Vec<u64> = vec![0; width.div_ceil(64)];
    let mut hpos: Vec<u32> = Vec::new();   // the wide primes' hits of this polynomial, packed
    let mut hchunk: Vec<(u32, u32)> = Vec::new();   // (start in hpos, first prime index) per chunk
    let mut mpos: Vec<u32> = Vec::new();   // scratch for one chunk of a mid band's hits
    let mut mr0: Vec<u32> = vec![0; nmid];   // the mid roots advanced block by block (`SIQS_MID=block`)
    let mut mr1: Vec<u32> = vec![0; nmid];
    // a's primes inside the resieve band have no sieve roots, so the resieve cannot tag them;
    // they are divided out of every candidate explicitly, like `large_a`
    let mid_a: Vec<usize> = fam.qs.iter().copied().filter(|&qi| qi >= resieve_from && qi < first_large).collect();
    let mid_a32: Vec<u32> = mid_a.iter().map(|&j| j as u32).collect();
    let large_a32: Vec<u32> = large_a.iter().map(|&j| j as u32).collect();
    let mut sp_s0: Vec<u32> = vec![0; resieve_from];
    let mut sp_s1: Vec<u32> = vec![0; resieve_from];
    let sp_special: Vec<usize> = (0..resieve_from).filter(|&j| base[j].p < 5 || fam.ainv[j] == 0).collect();
    // a special prime gets limit 0 and inverse 0: 0·d = 0 ≤ 0 would flag it, so use u32::MAX as the
    // "never" limit trick in reverse: inverse 1, limit 0 flags only d = 0, and the scalar pass overwrites
    let sp_inv: Vec<u32> = (0..resieve_from).map(|j| if sp_special.contains(&j) { 1 } else { base[j].inv }).collect();
    let sp_lim: Vec<u32> = (0..resieve_from).map(|j| if sp_special.contains(&j) { 0 } else { base[j].limit }).collect();
    // the small primes are not sieved: take their expected log off the threshold, then start
    // every cell at 128 − threshold so a hit is a byte with its high bit set, found a word at a time
    // whether this sieve runs the exact threshold (see `exact_mode`); the unsieved bound follows it
    let exact = match EXACT_OVERRIDE.with(|o| o.get()) { Some(v) => v, None => exact_mode() && n.bits() >= EXACT_MIN_BITS };
    let small = small_bound_with(exact);
    let small_log: f32 = base.iter().filter(|q| q.p < small).map(|q| 2.0 * (q.p as f32).log2() / q.p as f32).sum();
    // `SIQS_SPV_BITS=b` (YAFU's small-prime-variation correction, 20 bits there): scan with the
    // threshold lowered by b, then keep a hit only if its sieve sum plus the exact logs of the
    // unsieved small primes (found by the divisibility test on the roots) reaches the full
    // threshold. Unset: no filter, today's threshold, identical relations.
    let spv: Option<f32> = std::env::var("SIQS_SPV_BITS").ok().and_then(|v| v.parse().ok());
    let thr = (threshold_log - small_log - spv.unwrap_or(0.0)).clamp(1.0, 120.0) as u8;
    let init = 128u8 - thr;
    let small_idx: Vec<usize> = (0..first_large).filter(|&j| base[j].p < small && fam.ainv[j] != 0).collect();
    let small_logs: Vec<f32> = small_idx.iter().map(|&j| (base[j].p as f32).log2()).collect();
    // the polynomial's coefficients from b: g(x) = (a·x + lin)·x + c.
    // Q1: lin = 2b, c = (b² − n)/a.  Q2 (following YAFU): lin = b, c = (b² − n)/4a.
    let coeffs = |b: &BigInt| -> (BigInt, BigInt) {
        let c = (b * b - &nb) / &ab;
        if fam.q2 { (b.clone(), c / 4) } else { (b * 2, c) }
    };
    // **the exact threshold** (`SIQS_EXACT=1`, see `exact_mode`): the per-position template the
    // interval starts from, the remainder a report may keep, and the unsieved primes' reductions
    let af = fam.a.to_f64().unwrap();
    let (mut sp, tmpl, runs, allow, knoa) = if exact {
        let (lin0, c0) = coeffs(&b);
        // g(x) = a·(x + lin/2a)² − disc/4a with disc = lin² − 4ac (4kN in this form), the same
        // for every b of the family
        let disc = &lin0 * &lin0 - BigInt::from(4) * &ab * &c0;
        let knoa = disc.to_f64().unwrap() / (4.0 * af);
        let lp_bits = ((pmax * params.large) as f64).log2();
        let allow = (if double_large() { lp_bits * dlp_exp() } else { lp_bits }) as f32 + exact_margin();
        let ecap = exact_ecap();
        let lg_max = knoa.log2() as f32;
        // near a root of g the value is tiny and the rule would report every position: hold the
        // threshold at 12 bits under the interval's highest, which costs a few dozen positions
        let floor = (lg_max - allow - ecap - 12.0).clamp(1.0, 120.0);
        // one starting byte per 64 positions for the report test, and the same as runs of equal
        // bytes for filling the interval (a few dozen runs: the threshold moves a bit at a time),
        // so no interval-sized template sits in each thread's cache
        let mut tmpl: Vec<u8> = Vec::with_capacity(width.div_ceil(64));
        let mut runs: Vec<(usize, u8)> = Vec::new();   // (end position, byte)
        const STEP: usize = 64;
        let root = (knoa / af).sqrt();
        let g_abs = |x: f64| (af * x * x - knoa).abs();
        let mut lo = 0usize;
        while lo < width {
            let hi = (lo + STEP).min(width);
            // the smallest |g| over the chunk, taken a few positions wider for the shift lin/2a
            let (x0, x1) = (lo as f64 - half as f64 - 8.0, hi as f64 - half as f64 + 8.0);
            let crosses = (x0 - root) * (x1 - root) <= 0.0 || (x0 + root) * (x1 + root) <= 0.0;
            let lg_min = if crosses { 0.0 } else { g_abs(x0).min(g_abs(x1)).max(1.0).log2() as f32 };
            let thr = (lg_min - allow - ecap).floor().clamp(floor, 120.0) as u8;
            tmpl.push(128 - thr);
            match runs.last_mut() { Some(r) if r.1 == 128 - thr => r.0 = hi, _ => runs.push((hi, 128 - thr)) }
            lo = hi;
        }
        (Some(SmallPart::new(base, small, &ab, &disc)), tmpl, runs, allow, knoa)
    } else { (None, Vec::new(), Vec::new(), 0.0, 0.0) };
    for poly in 0..polys {
        // sieve: position i is x = i − half; p | Q(x) iff x ≡ root (mod p)
        let t0 = Stamp::now();
        // in blocks that fit L1: every prime's next position is carried from block to block,
        // so the random byte adds land in cache; the hits of a block are read while it is there
        // the first polynomial sets every start from its roots; after that the Gray-code switch
        // below carries the starts along with the roots, so this is one copy, not 7,000 mod-adds
        if poly == 0 {
            for (i, q) in base.iter().enumerate().take(first_large) {
                if fam.ainv[i] == 0 { starts[i] = [u32::MAX; 2]; continue; }
                starts[i] = [shift_residue(roots[i][0], half_mod[i], q.p as u32, true), shift_residue(roots[i][1], half_mod[i], q.p as u32, true)];
            }
        }
        let mut next: Vec<[u32; 2]> = starts[..first_large].to_vec();
        if exact { let mut lo = 0; for &(hi, v) in &runs { acc[lo..hi].fill(v); lo = hi; } } else { acc[..width].iter_mut().for_each(|v| *v = init); }
        // this polynomial's coefficients: the candidate stage evaluates g with them, and the
        // exact threshold reduces them modulo the unsieved primes
        let (lin, c) = coeffs(&b);
        let boa = if exact { lin.to_f64().unwrap() / (2.0 * af) } else { 0.0 };
        if let Some(sp) = sp.as_mut() { sp.set_poly(&lin, &c); }
        count(17, t0.ns());
        // the direct primes: shift every position by this switch's delta in one pass, then add
        // each prime's log at its strides straight into the interval (no bucket entry written,
        // read back and applied). A stride past the end adds to the sink byte at acc[width].
        // The shift and the adds run chunk by chunk so a chunk's roots are still in L1 when
        // they are added: one trip through memory per polynomial instead of two. At 16
        // threads on the desktop the separate passes were memory-bandwidth bound (polynomials
        // per second fell from 10,721 at 8 threads to 8,320 at 16, 2026-11-12).
        let tb = Stamp::now();
        let step = if poly > 0 { Some(gray_step(poly)) } else { None };
        hpos.clear();
        hchunk.clear();
        const CHUNK: usize = 1024;
        for &(a, b, k) in &bands {
            let mut lo = a;
            while lo < b {
                let hi = (lo + CHUNK).min(b);
                let tsh = Stamp::now();
                if let Some((v, add)) = step { shift_all(&mut dr0[lo..hi], &dd[v][lo..hi], &dp[lo..hi], add, avx2 && shift_simd); shift_all(&mut dr1[lo..hi], &dd[v][lo..hi], &dp[lo..hi], add, avx2 && shift_simd); }
                let tad = Stamp::now();
                count(20, tad.since(tsh));
                if k == 1 && rec {
                    for r in [&dr0[lo..hi], &dr1[lo..hi]] {
                        let start = hpos.len();
                        hchunk.push((start as u32, (first_large + lo) as u32));
                        compact_hits(r, width as u32, &mut hpos, avx2);
                        add_hits(&mut acc[..width], &hpos[start..], &dlog[lo..hi]);
                    }
                } else if k == 1 {
                    direct_add(&mut acc[..width + SINK], &dr0[lo..hi], &dlog[lo..hi], width as u32, sink);
                    direct_add(&mut acc[..width + SINK], &dr1[lo..hi], &dlog[lo..hi], width as u32, sink);
                } else if mid_blockwise {
                    mr0[lo..hi].copy_from_slice(&dr0[lo..hi]);
                    mr1[lo..hi].copy_from_slice(&dr1[lo..hi]);
                } else if rec_mid {
                    for r in [&dr0[lo..hi], &dr1[lo..hi]] {
                        mpos.clear();
                        compact_hits_k(r, &dp[lo..hi], k, width as u32, &mut mpos, avx2);
                        add_hits(&mut acc[..width], &mpos, &dlog[lo..hi]);
                    }
                } else {
                    direct_add_k(&mut acc[..width + SINK], &dr0[lo..hi], &dp[lo..hi], &dlog[lo..hi], width as u32, k, sink);
                    direct_add_k(&mut acc[..width + SINK], &dr1[lo..hi], &dp[lo..hi], &dlog[lo..hi], width as u32, k, sink);
                }
                // the adds of this chunk: wide primes (one stride) into slot 22, mid primes into 21
                count(if k == 1 { 22 } else { 21 }, tad.ns());
                lo = hi;
            }
        }
        let tmb = Stamp::now();
        if mid_blockwise {
            let mut blo = 0usize;
            while blo < width {
                let bhi = (blo + mid_bs).min(width);
                for &(a, b, k) in &bands {
                    if k == 1 { continue; }
                    let mut lo = a;
                    while lo < b {
                        let hi = (lo + (1 << (32 - POS_BITS))).min(b);
                        // strides per block from the band's smallest possible prime (width / k), not
                        // from dp[lo], which may be a skipped prime's u32::MAX sentinel
                        let kb = (mid_bs as u32).div_ceil(width.div_ceil(k as usize) as u32).min(k);
                        for r in [&mut mr0[lo..hi], &mut mr1[lo..hi]] {
                            mpos.clear();
                            mid_block(r, &dp[lo..hi], bhi as u32, kb, &mut mpos, avx2);
                            add_hits(&mut acc[..width], &mpos, &dlog[lo..hi]);
                        }
                        lo = hi;
                    }
                }
                blo = bhi;
            }
        }
        count(21, tmb.ns());
        let bucket_ns = tb.ns();
        count(5, bucket_ns);
        let mut hits: Vec<usize> = Vec::new();
        // reports of this polynomial, added to the shared counter once at the end: one atomic
        // add per report on a line every thread writes made the test look three times dearer
        let reports = std::cell::Cell::new(0u64);
        let mut scan_ns = 0u64;
        let mut lo = 0usize;
        while lo < width {
            let hi = (lo + seg).min(width);
            let block = &mut acc[lo..hi];
            let tf = Stamp::now();
            // (a fixed-count variant with a sink byte for the primes from 2048 up, to avoid the
            // loop-exit mispredicts, was 5 % slower than this loop, 2026-11-12: not kept)
            for (i, q) in base.iter().enumerate().take(band) {
                if fam.ainv[i] == 0 || q.p < small { continue; }
                let p = q.p as usize;
                let lg = q.log;
                fill_block(block, lo, hi, p, lg, &mut next[i]);
            }
            let tf2 = Stamp::now();
            count(10, tf2.since(tf));
            for (i, q) in base.iter().enumerate().take(first_large).skip(band) {
                if fam.ainv[i] == 0 { continue; }
                fill_block(block, lo, hi, q.p as usize, q.log, &mut next[i]);
            }
            count(11, tf2.ns());
            // one byte per prime, 20-200 KB, so the log lookup stays in cache whatever order the
            // entries arrive in (reading it from the 32-byte Prime cost 1 ns an entry)
            let ts = Stamp::now();
            let full = block.len() / 8 * 8;
            let accept = |pos: usize, v: u8| -> bool {
                if let Some(sp) = &sp {
                    // the exact threshold: this position's own size and sieve sum, and the exact
                    // small-prime part only when the sum alone falls short
                    reports.set(reports.get() + 1);
                    let xi = pos as i64 - half;
                    let t = xi as f64 + boa;
                    let g = (af * t * t - knoa).abs();
                    if g <= 1.0 { return true; }
                    let need = log2_fast(g) - allow;
                    let sum = v.wrapping_sub(tmpl[pos >> 6]) as f32;
                    return sum >= need || sum + sp.bits(xi, pos as u32, &starts) >= need;
                }
                let Some(_) = spv else { return true };
                let iu = pos as u32;
                let mut sum = (v - init) as f32;
                for (&j, &lg) in small_idx.iter().zip(&small_logs) {
                    let [s0, s1] = starts[j];
                    if base[j].divides(iu.abs_diff(s0)) || base[j].divides(iu.abs_diff(s1)) { sum += lg; }
                }
                sum >= threshold_log
            };
            if vector_scan || cfg!(all(target_arch = "wasm32", target_feature = "simd128")) {
                let start = hits.len();
                #[cfg(target_arch = "x86_64")]
                // SAFETY: AVX2 was detected; the scanner bounds every load
                // and handles the tail without reading outside this block.
                unsafe { scan_block_avx2(block, lo, &mut hits); }
                #[cfg(all(target_arch = "wasm32", target_feature = "simd128"))]
                // SAFETY: whole sixteen-byte loads inside the block, the tail scalar
                unsafe { wasm_simd::scan(block, lo, &mut hits); }
                if spv.is_some() || exact {
                    let mut write = start;
                    for read in start..hits.len() {
                        let pos = hits[read];
                        if accept(pos, block[pos - lo]) { hits[write] = pos; write += 1; }
                    }
                    hits.truncate(write);
                }
            } else {
                for (w, chunk) in block[..full].chunks_exact(8).enumerate() {
                    let word = u64::from_le_bytes(chunk.try_into().unwrap());
                    if word & 0x8080_8080_8080_8080 != 0 { for (k, &v) in chunk.iter().enumerate() { if v & 0x80 != 0 && accept(lo + w * 8 + k, v) { hits.push(lo + w * 8 + k); } } }
                }
                for (k, &v) in block[full..].iter().enumerate() { if v & 0x80 != 0 && accept(lo + full + k, v) { hits.push(lo + full + k); } }
            }
            scan_ns += ts.ns();
            lo = hi;
        }
        count(0, t0.ns() - scan_ns - bucket_ns);
        count(1, scan_ns);
        count(16, reports.get());

        count(4, hits.len() as u64);
        let t2 = Stamp::now();
        let tr = Stamp::now();
        // resieve: the medium primes walk the interval once more and tag the candidates they hit
        let mut cand_primes: Vec<Vec<u32>> = vec![Vec::new(); hits.len()];
        if !hits.is_empty() {
            for &i in &hits { hitbits[i >> 6] |= 1u64 << (i & 63); }
            let is_hit = |pos: usize| hitbits[pos >> 6] >> (pos & 63) & 1 == 1;
            let which = |pos: usize| hits.binary_search(&pos).unwrap();
            if resieve_from < first_large {
                for j in resieve_from..first_large {
                    if fam.ainv[j] == 0 { continue; }
                    let p = base[j].p as usize;
                    let [s0, s1] = starts[j];
                    let mut pos = s0 as usize;
                    while pos < width { if is_hit(pos) { cand_primes[which(pos)].push(j as u32); } pos += p; }
                    if s1 != s0 { let mut pos = s1 as usize; while pos < width { if is_hit(pos) { cand_primes[which(pos)].push(j as u32); } pos += p; } }
                }
            }
            // the direct primes that landed on a candidate. Mid primes (several strides): the
            // divisibility test of the candidate against both roots, eight lanes at a time, per
            // candidate. Wide primes (one stride): one compare pass over the root arrays per eight
            // survivors.
            for (k, &x) in hits.iter().enumerate() {
                divides_flags(x as u32, &dr0[..nmid], &dr1[..nmid], &dinv[..nmid], &dlim[..nmid], &mut md_flags, avx2);
                for (w, &bits) in md_flags.iter().enumerate() {
                    let mut bits = bits;
                    while bits != 0 { let t = bits.trailing_zeros() as usize; cand_primes[k].push((first_large + w * 8 + t) as u32); bits &= bits - 1; }
                }
            }
            let twt = Stamp::now();
            if rec && wide_tag_bitmap {
                // one pass over the recorded wide hits, each tested against the candidate bitmap
                // (`SIQS_WIDE_TAG=bitmap`, opt-in: slower here, see the switch)
                for (c, &(start, idx0)) in hchunk.iter().enumerate() {
                    let end = hchunk.get(c + 1).map_or(hpos.len(), |&(e, _)| e as usize);
                    tag_bitmap(&hpos[start as usize..end], idx0, &hitbits, &mut |pos, j| cand_primes[which(pos)].push(j));
                }
            }
            for grp in hits.chunks(8) {
                if rec && wide_tag_bitmap { break; }
                let mut surv = [u32::MAX; 8];
                for (j, &h) in grp.iter().enumerate() { surv[j] = h as u32; }
                if rec {
                    for (c, &(start, idx0)) in hchunk.iter().enumerate() {
                        let end = hchunk.get(c + 1).map_or(hpos.len(), |&(e, _)| e as usize);
                        run_hits(&hpos[start as usize..end], (1 << POS_BITS) - 1, &surv, avx2, |_, e| { cand_primes[which((e & ((1 << POS_BITS) - 1)) as usize)].push(idx0 + (e >> POS_BITS)); });
                    }
                } else {
                    for r in [&dr0[nmid..], &dr1[nmid..]] {
                        run_hits(r, u32::MAX, &surv, avx2, |k, pos| { cand_primes[which(pos as usize)].push((first_wide + k) as u32); });
                    }
                }
            }
            count(19, twt.ns());
            for &i in &hits { hitbits[i >> 6] = 0; }
        }
        count(18, tr.ns());
        // c = (b² − n)/a; Q(x)/a = a x² + 2 b x + c. Below about 57 digits this fits i128 and
        // the whole candidate step runs in machine integers; above, the same in BigInt.
        let small = ab.bits() + 2 * (params.half as u64).ilog2() as u64 + 4 < 126 && c.bits() < 125 && b.bits() + (params.half as u64).ilog2() as u64 + 3 < 126;
        let (a128, b128, c128) = if small { (ab.to_i128().unwrap(), lin.to_i128().unwrap(), c.to_i128().unwrap()) } else { (0, 0, 0) };
        // the small primes' roots for this polynomial as flat arrays, so the divisibility test
        // per candidate runs eight lanes at a time (the loop over the 40-byte `Prime` structs
        // was 1.7 ns a prime, 2026-11-12); the special cases (2, primes dividing a or k·n,
        // 3) stay scalar in `sp_special`
        if !hits.is_empty() && !small {
            for (j, &[s0, s1]) in starts.iter().enumerate().take(resieve_from) { sp_s0[j] = s0; sp_s1[j] = s1; }
        }
        let mut sp_flags: Vec<u8> = vec![0; resieve_from];
        for (k, &i) in hits.iter().enumerate() {
            let xi = i as i64 - half;
            let mut exps: Vec<(u32, u32)> = fam.qs.iter().map(|&qi| (qi as u32 + 1, 1)).collect();
            // Q2 divides out 4a rather than a. The two powers of 2 must be restored
            // in every relation, including those that enter the large-prime graph.
            if fam.q2 { exps.push((1, 2)); }
            let bump = |exps: &mut Vec<(u32, u32)>, c: u32| { match exps.iter_mut().find(|(d, _)| *d == c) { Some(t) => t.1 += 1, None => exps.push((c, 1)) } };
            let large;
            if small {
                let x = xi as i128;
                let qa = (a128 * x + b128) * x + c128;
                if qa == 0 { continue; }
                if qa < 0 { exps.push((0, 1)); }
                let mut v = qa.unsigned_abs();
                let iu = i as u32;
                for (j, q) in base.iter().enumerate().take(resieve_from) {
                    let p = q.p as u128;
                    let divides = if fam.ainv[j] == 0 || q.p < 5 { v % p == 0 } else {
                        let [s0, s1] = starts[j];
                        q.divides(iu.abs_diff(s0)) || q.divides(iu.abs_diff(s1))
                    };
                    if divides { while v % p == 0 { v /= p; bump(&mut exps, j as u32 + 1); } }
                }
                for &j in &cand_primes[k] { let p = base[j as usize].p as u128; while v % p == 0 { v /= p; bump(&mut exps, j + 1); } }
                for &j in &mid_a { let p = base[j].p as u128; while v % p == 0 { v /= p; bump(&mut exps, j as u32 + 1); } }
                for &j in &large_a { let p = base[j].p as u128; while v % p == 0 { v /= p; bump(&mut exps, j as u32 + 1); } }
                let lb = (pmax * params.large) as u128;
                if v == 1 { large = vec![]; }
                else if v < lb { large = vec![(v as u64, 1)]; }
                else if let Some(pair) = split_cofactor(v, lb, pmax as u128) { large = pair; }
                else { continue; }
            } else {
                let tq = Stamp::now();
                let x = BigInt::from(xi);
                let qa: BigInt = (&ab * &x + &lin) * &x + &c;
                if qa.is_zero() { continue; }
                if qa.is_negative() { exps.push((0, 1)); }
                let mut v = qa.magnitude().clone();
                let iu = i as u32;
                let ts = Stamp::now();
                count(12, ts.since(tq));
                let mut hitlist: Vec<u32> = Vec::with_capacity(16);
                for ((((f, &s0), &s1), &inv), &lim) in sp_flags.iter_mut().zip(&sp_s0).zip(&sp_s1).zip(&sp_inv).zip(&sp_lim) {
                    *f = ((iu.abs_diff(s0).wrapping_mul(inv) <= lim) | (iu.abs_diff(s1).wrapping_mul(inv) <= lim)) as u8;
                }
                for &j in &sp_special {
                    let q = &base[j];
                    sp_flags[j] = if q.p == 2 { v.iter_u32_digits().next().is_none_or(|d| d & 1 == 0) } else { rem_u32(&v, q.p as u32) == 0 } as u8;
                }
                for (j, &f) in sp_flags.iter().enumerate() { if f != 0 { hitlist.push(j as u32); } }
                let td = Stamp::now();
                count(13, td.since(ts));
                // 2 first (the remainder trick is for odd moduli), then everything else in batches
                if base[0].p == 2 && hitlist.first() == Some(&0) { while (&v % 2u32).is_zero() { v /= 2u32; bump(&mut exps, 1); } }
                divide_out(&mut v, &hitlist, base, &mut exps);
                divide_out(&mut v, &cand_primes[k], base, &mut exps);
                divide_out(&mut v, &mid_a32, base, &mut exps);
                divide_out(&mut v, &large_a32, base, &mut exps);
                let tc = Stamp::now();
                count(14, tc.since(td));
                let lb = (pmax * params.large) as u128;
                if v.is_one() { large = vec![]; }
                else if v < large_bound { large = vec![(v.to_u64().unwrap(), 1)]; }
                else if v.bits() <= 63 { match split_cofactor(v.to_u128().unwrap(), lb, pmax as u128) { Some(pair) => large = pair, None => { tick(15, tc); continue } } }
                else { tick(15, tc); continue; }
                let y = (&y_a * BigInt::from(xi) + &b).mod_floor(&nb).to_biguint().unwrap();
                out.push(Relation { y, exps, large });
                tick(15, tc);
                continue;
            }
            let y = (&y_a * BigInt::from(xi) + &b).mod_floor(&nb).to_biguint().unwrap();
            out.push(Relation { y, exps, large });
        }
        tick(2, t2);
        // next b by the Gray code: flip B_l where l is the lowest set bit of poly+1
        let t3 = Stamp::now();
        if poly + 1 < polys {
            let v = (poly + 1).trailing_zeros() as usize;
            let sign_plus = ((poly + 1) >> (v + 1)) & 1 == 0;   // whether B_v is added or subtracted next
            let bl = BigInt::from(fam.bs[v].clone());
            if sign_plus { b -= 2 * &bl; } else { b += 2 * &bl; }
            let dv = &deltas[v];
            for (i, q) in base.iter().enumerate().take(first_large) {
                if fam.ainv[i] == 0 { continue; }
                let d = dv[i];
                let p = q.p as u32;
                for r in roots[i].iter_mut() { *r = shift_residue(*r, d, p, sign_plus); }
                for r in starts[i].iter_mut() { *r = shift_residue(*r, d, p, sign_plus); }
            }
        }
        tick(3, t3);
    }
    (out, polys)
}

/// **The cycle finder** for partial relations. A partial with large primes p (and q) is an
/// edge p–q of a graph whose vertex 1 stands for "no second prime"; a cycle multiplies to a
/// relation in which every large prime appears an even number of times. A spanning forest is
/// kept as union-find plus the tree edges; an edge that lands inside one tree closes a cycle,
/// and the path between its ends through the tree edges is the rest of it.
#[derive(Default)]
pub struct Cycles { pub edges: Vec<Relation>, up: std::collections::HashMap<u64, (u64, usize)>, size: std::collections::HashMap<u64, usize> }

impl Cycles {
    fn root(&self, mut v: u64) -> (u64, usize) {
        let mut depth = 0;
        while let Some(&(p, _)) = self.up.get(&v) { v = p; depth += 1; }
        (v, depth)
    }
    /// Make `x` the root of its tree by reversing the parent pointers on its path to the root.
    fn reroot(&mut self, x: u64) {
        let mut path = Vec::new();
        let mut v = x;
        while let Some(&(p, e)) = self.up.get(&v) { path.push((v, p, e)); v = p; }
        for &(child, parent, e) in &path { self.up.remove(&child); self.up.insert(parent, (child, e)); }
        self.up.remove(&x);
    }
    /// Add a partial; when it closes a cycle, the relations of that cycle (this one included).
    pub fn add(&mut self, r: Relation) -> Option<Vec<Relation>> {
        // Only odd exponents need partners. Preserve even powers in the relation for
        // square-root reconstruction, but never put them into the partial forest.
        let mut odd = r.large.iter().filter(|(_, e)| e % 2 == 1).map(|(p, _)| *p);
        let (u, v) = match (odd.next(), odd.next(), odd.next()) {
            (None, None, None) => return Some(vec![r]),
            (Some(p), None, None) => (1u64, p),
            (Some(p), Some(q), None) => (p, q),
            _ => return None,
        };
        let idx = self.edges.len();
        self.edges.push(r);
        let (ru, _) = self.root(u);
        let (rv, _) = self.root(v);
        if ru != rv {
            // join: re-root the smaller tree at its end of the edge and hang it under the other
            let (su, sv) = (*self.size.get(&ru).unwrap_or(&1), *self.size.get(&rv).unwrap_or(&1));
            let (small, big, small_root, big_root) = if su <= sv { (u, v, ru, rv) } else { (v, u, rv, ru) };
            self.reroot(small);
            self.up.insert(small, (big, idx));
            self.size.remove(&small_root);
            self.size.insert(big_root, su + sv);
            return None;
        }
        // the cycle: both paths to the root, cut at the common ancestor
        let mut seen: std::collections::HashMap<u64, usize> = Default::default();
        let mut x = u;
        let mut i = 0;
        loop { seen.insert(x, i); i += 1; match self.up.get(&x) { Some(&(p, _)) => x = p, None => break } }
        let mut path = vec![self.edges[idx].clone()];
        let mut y = v;
        let mut up_from_v: Vec<usize> = Vec::new();
        while !seen.contains_key(&y) { let &(p, e) = self.up.get(&y).unwrap(); up_from_v.push(e); y = p; }
        let meet = y;
        let mut x = u;
        while x != meet { let &(p, e) = self.up.get(&x).unwrap(); path.push(self.edges[e].clone()); x = p; }
        for e in up_from_v { path.push(self.edges[e].clone()); }
        Some(path)
    }
}

/// Collect relations across threads until the base is covered, merging partials.
pub fn relations(n: &BigUint, base: &[Prime], params: &Params, need: usize) -> (Vec<Relation>, usize, usize) {
    // threshold: log2 of the typical value M·√(n/2) less the slack and the large-prime room
    let size = (params.half as f64 * (n.to_f64().unwrap() / 2.0).sqrt()).log2()
        - if q2_polynomials(n) { 1.0 } else { 0.0 };
    let pmax = base.last().unwrap().p as f64;
    let threshold = (size - params.slack as f64 * pmax.log2() - (pmax * params.large as f64).log2()) as f32;
    let mut full: Vec<Relation> = Vec::new();
    let mut graph = Cycles::default();
    let mut round = 0u64;
    // progress with an estimate, like Alpern's countdown: remaining relations over the rate of
    // the last stretch, printed to stderr every few seconds when `SIQS_PROGRESS` is set
    let progress = std::env::var("SIQS_PROGRESS").map(|v| v != "0").unwrap_or(false);
    // `SIQS_BUDGET_SECS` stops the collection after that many seconds, for sweeps that compare
    // relations a second across parameter rows without running each to the end
    let budget: Option<f64> = std::env::var("SIQS_BUDGET_SECS").ok().and_then(|v| v.parse().ok());
    // Fixed work for reproducible kernel comparisons. Workers never start a family past
    // this bound, so counters include exactly the same work at any thread count.
    let family_limit: Option<u64> = std::env::var("SIQS_FAMILIES").ok()
        .and_then(|v| v.parse().ok()).filter(|&v| v > 0);
    let started = clock::Instant::now();
    let mut last_print = started;
    let mut mark = (started, 0usize);
    let mut polys = 0usize;
    // Families stream through a channel from `threads` workers and are merged in family order
    // on this thread, so no worker waits at a round barrier for the slowest core (the
    // round-of-`threads` version idled 35 % of 8 threads on this desktop, 2026-11-12). Merging in
    // index order keeps the relation set a deterministic function of the seeds: the stop
    // decision is taken after family k from the first k families alone, whichever finished first.
    let threads = rayon::current_num_threads().max(1);
    // one family's relations into the graph and the full list; true when collection is done
    report_progress(1, 0, need, 0, 0);
    let mut absorb = |rels: Vec<Relation>, np: usize| -> bool {
        polys += np;
        for r in rels {
            if r.large.is_empty() { full.push(r); continue; }
            if let Some(cycle) = graph.add(r) {
                let mut acc = cycle[0].clone();
                for r in &cycle[1..] { acc = acc.times(r, n); }
                full.push(acc);
            }
        }
        round += 1;
        report_progress(1, full.len(), need, graph.edges.len(), polys);
        if full.len() >= need || family_limit.is_some_and(|limit| round >= limit)
            || budget.is_some_and(|b| started.elapsed().as_secs_f64() >= b) { return true; }
        if progress && last_print.elapsed().as_secs_f64() >= 5.0 {
            let now = clock::Instant::now();
            let elapsed = now.duration_since(started).as_secs_f64();
            let (t0, c0) = mark;
            let recent = (full.len().saturating_sub(c0)) as f64 / now.duration_since(t0).as_secs_f64().max(1e-9);
            let overall = full.len() as f64 / elapsed.max(1e-9);
            let rate = if full.len() > c0 + 20 { recent } else { overall };
            let left = if rate > 0.0 { (need.saturating_sub(full.len())) as f64 / rate } else { f64::INFINITY };
            eprintln!("  {}/{} relations, {} partial edges, {} polynomials, {:.0} s elapsed, about {:.0} s to go", full.len(), need, graph.edges.len(), polys, elapsed, left);
            if now.duration_since(t0).as_secs_f64() >= 20.0 { mark = (now, full.len()); }
            last_print = now;
        }
        false
    };
    if threads == 1 {
        // one thread (and the browser, where std threads do not exist): the same families in
        // the same order on the calling thread, so the relation set is identical
        let mut idx = 1u64;
        loop {
            let fam = choose_a(n, base, params.half, idx);
            let (rels, np) = sieve_family(n, base, &fam, params, threshold);
            if absorb(rels, np) { break; }
            idx += 1;
        }
    } else if cfg!(target_arch = "wasm32") {
        // the browser: the only threads are rayon's pool (wasm-bindgen-rayon's Web Workers), and
        // the caller is outside it, so `in_place_scope` keeps the merge on the caller and puts
        // one long sieve task on each worker (a plain `rayon::scope` would park the merge on a
        // worker, and with one worker deadlock)
        let next = std::sync::atomic::AtomicU64::new(1);
        let stop = std::sync::atomic::AtomicBool::new(false);
        let (tx, rx) = std::sync::mpsc::channel::<(u64, Vec<Relation>, usize)>();
        rayon::in_place_scope(|sc| {
            for _ in 0..threads {
                let tx = tx.clone();
                let (next, stop) = (&next, &stop);
                sc.spawn(move |_| {
                    while !stop.load(Relaxed) {
                        let idx = next.fetch_add(1, Relaxed);
                        if family_limit.is_some_and(|limit| idx > limit) { break; }
                        let fam = choose_a(n, base, params.half, idx);
                        let (rels, np) = sieve_family(n, base, &fam, params, threshold);
                        if tx.send((idx, rels, np)).is_err() { break; }
                    }
                });
            }
            drop(tx);
            let mut pending: std::collections::BTreeMap<u64, (Vec<Relation>, usize)> = std::collections::BTreeMap::new();
            let mut want = 1u64;
            'merge: while let Ok((idx, rels, np)) = rx.recv() {
                pending.insert(idx, (rels, np));
                while let Some((rels, np)) = pending.remove(&want) {
                    want += 1;
                    if absorb(rels, np) { stop.store(true, Relaxed); break 'merge; }
                }
            }
            stop.store(true, Relaxed);
            for _ in rx.try_iter() {}
        });
    } else {
        let next = std::sync::atomic::AtomicU64::new(1);
        let stop = std::sync::atomic::AtomicBool::new(false);
        let (tx, rx) = std::sync::mpsc::channel::<(u64, Vec<Relation>, usize)>();
        // plain threads, not rayon::scope: with a one-thread pool the scope body would hold the
        // only worker while waiting on the channel (deadlock, 2026-11-12); the count still
        // follows RAYON_NUM_THREADS through rayon's pool size
        std::thread::scope(|sc| {
            for _ in 0..threads {
                let tx = tx.clone();
                let (next, stop) = (&next, &stop);
                sc.spawn(move || {
                    while !stop.load(Relaxed) {
                        let idx = next.fetch_add(1, Relaxed);
                        if family_limit.is_some_and(|limit| idx > limit) { break; }
                        let fam = choose_a(n, base, params.half, idx);
                        let (rels, np) = sieve_family(n, base, &fam, params, threshold);
                        if tx.send((idx, rels, np)).is_err() { break; }
                    }
                });
            }
            drop(tx);
            let mut pending: std::collections::BTreeMap<u64, (Vec<Relation>, usize)> = std::collections::BTreeMap::new();
            let mut want = 1u64;
            'merge: while let Ok((idx, rels, np)) = rx.recv() {
                pending.insert(idx, (rels, np));
                while let Some((rels, np)) = pending.remove(&want) {
                    want += 1;
                    if absorb(rels, np) { stop.store(true, Relaxed); break 'merge; }
                }
            }
            stop.store(true, Relaxed);
            // drain what the workers still send so the scope can end
            for _ in rx.try_iter() {}
        });
    }
    let _ = round;
    (full, polys, graph.edges.len())
}

/// A nontrivial factor of the odd composite n, or None.
pub fn factor(n: &BigUint) -> Option<BigUint> { factor_with(n, Params::for_n(n)) }

pub fn factor_with(n: &BigUint, params: Params) -> Option<BigUint> {
    report_progress(0, 0, 0, 0, 0);
    if *n < BigUint::from(4u32) { return None; }
    for p in primes_upto(10_000) { if (n % p).is_zero() { return if *n == BigUint::from(p) { None } else { Some(BigUint::from(p)) }; } }
    let r = n.sqrt();
    if &r * &r == *n { return Some(r); }
    if rulerpack::big::is_prime_big(n) { return None; }
    let k = multiplier(n);
    let kn = n * k;
    if std::env::var("SIQS_QUIET").is_err() {
        eprintln!("  siqs: input {} digits / {} bits; kN {} digits / {} bits; multiplier {k}",
            n.to_string().len(), n.bits(), kn.to_string().len(), kn.bits());
        eprintln!("  base {} half {} ({} blocks/side) large {} slack {} extra {}; {} threads; {}; {}",
            params.base, params.half, params.half / 32768, params.large, params.slack, params.extra,
            rayon::current_num_threads(), if have_avx2() { "AVX2" } else { "portable" },
            if q2_polynomials(&kn) { "Q2" } else { "Q1" });
    }
    let base = factor_base(&kn, params.base);
    let need = base.len() + 1 + params.extra;
    let (rels, polys, partials) = relations(&kn, &base, &params, need);
    report_progress(2, rels.len(), need, partials, polys);
    let n_only = n;
    let n = &kn;
    let cols = base.len() + 1;
    let odd: Vec<Vec<u32>> = rels.iter().map(|r| { let mut v: Vec<u32> = r.exps.iter().filter(|(_, e)| e % 2 == 1).map(|(c, _)| *c).collect(); v.sort_unstable(); v }).collect();
    let keep = prune_singletons(&odd, cols);
    let kept_rows: Vec<Vec<u32>> = keep.iter().map(|&i| odd[i].clone()).collect();
    let deps = solve(&kept_rows, cols);
    report_progress(3, rels.len(), need, partials, polys);
    for dep in deps {
        if dep.is_empty() { continue; }
        let mut x = BigUint::one();
        let mut sum = vec![0u32; cols];
        let mut larges: std::collections::HashMap<u64, u32> = Default::default();
        for &k in &dep {
            let r = &rels[keep[k]];
            x = &x * &r.y % n;
            for &(c, e) in &r.exps { sum[c as usize] += e; }
            for &(p, e) in &r.large { *larges.entry(p).or_default() += e; }
        }
        if sum.iter().any(|e| e % 2 == 1) || larges.values().any(|e| e % 2 == 1) { continue; }
        let mut y = BigUint::one();
        for (j, q) in base.iter().enumerate() { if sum[j + 1] > 0 { y = y * BigUint::from(q.p).modpow(&BigUint::from(sum[j + 1] / 2), n) % n; } }
        for (&p, &e) in &larges { y = y * BigUint::from(p).modpow(&BigUint::from(e / 2), n) % n; }
        let diff = if x >= y { &x - &y } else { &y - &x };
        for g in [diff.gcd(n_only), (&x + &y).gcd(n_only)] { if g > BigUint::one() && g < *n_only { return Some(g); } }
    }
    None
}

/// Rows to keep after removing, repeatedly, every row that is the only one with an odd
/// exponent in some column: such a row can be in no dependency, and dropping it shrinks the
/// matrix the elimination sees. The count-and-XOR form: per column the
/// number of live rows and the XOR of their ids, so when a count reaches one the XOR *is* the
/// remaining row, and the cascade runs in one pass over the nonzeros with a queue.
pub fn prune_singletons(odd: &[Vec<u32>], cols: usize) -> Vec<usize> {
    let mut count = vec![0u32; cols];
    let mut owner = vec![0usize; cols];
    for (i, r) in odd.iter().enumerate() { for &c in r { count[c as usize] += 1; owner[c as usize] ^= i; } }
    let mut alive = vec![true; odd.len()];
    let mut queue: std::collections::VecDeque<usize> = (0..cols).filter(|&c| count[c] == 1).collect();
    while let Some(c) = queue.pop_front() {
        if count[c] != 1 { continue; }
        let i = owner[c];
        if !alive[i] { continue; }
        alive[i] = false;
        for &d in &odd[i] { let d = d as usize; count[d] -= 1; owner[d] ^= i; if count[d] == 1 { queue.push_back(d); } }
    }
    (0..odd.len()).filter(|&i| alive[i]).collect()
}

/// The dependencies of the pruned matrix: dense elimination up to `LANCZOS_FROM` rows, block
/// Lanczos above it (or always with `SIQS_LANCZOS=1`, never with `SIQS_LANCZOS=0`).
pub const LANCZOS_FROM: usize = 6000;
pub fn solve(rows: &[Vec<u32>], cols: usize) -> Vec<Vec<usize>> {
    let force = std::env::var("SIQS_LANCZOS").ok();
    let use_lanczos = match force.as_deref() { Some("1") => true, Some("0") => false, _ => rows.len() >= LANCZOS_FROM };
    if !use_lanczos { return dependencies_sparse(rows, cols); }
    let mut out = Vec::new();
    for seed in 1..=4u64 {
        out.extend(rulerpack::lanczos::dependencies(rows, cols, seed));
        if out.len() >= 32 { break; }
    }
    if out.is_empty() { dependencies_sparse(rows, cols) } else { out }
}

/// The complete factorisation of `n`: trial division to 100,000, perfect squares, primality
/// tests, and the sieve on every composite piece until all are prime. Sorted, with
/// multiplicity. A piece the sieve cannot split is kept whole.
pub fn factor_completely(n: &BigUint) -> Vec<BigUint> {
    report_progress(0, 0, 0, 0, 0);
    fn split(n: &BigUint, out: &mut Vec<BigUint>) {
        if n.is_one() { return; }
        if rulerpack::big::is_prime_big(n) { out.push(n.clone()); return; }
        let r = n.sqrt();
        if &r * &r == *n { split(&r, out); split(&r, out); return; }
        match factor(n) {
            Some(f) => { split(&f, out); split(&(n / &f), out); }
            None => out.push(n.clone()),
        }
    }
    let mut n = n.clone();
    let mut factors = Vec::new();
    if n.is_zero() { return factors; }
    for p in rulerpack::smooth::primes_upto(100_000) {
        while (&n % p).is_zero() { n /= p; factors.push(BigUint::from(p)); }
    }
    if !n.is_one() { split(&n, &mut factors); }
    factors.sort();
    factors
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn incremental_factor_base_matches_the_original_scan() {
fn original(n: &BigUint, count: usize) -> Vec<Prime> {
    let mut out = vec![Prime::new(2, (n % 2u32).to_u64().unwrap())];
    let mut bound = 1000u64;
    loop {
        for p in primes_upto(bound).into_iter().skip(1) {
            if out.len() >= count { return out; }
            if out.iter().any(|q| q.p == p) { continue; }
            let np = (n % p).to_u64().unwrap();
            if np == 0 { out.push(Prime::new(p, 0)); } else if legendre(np, p) == 1 { out.push(Prime::new(p, sqrt_mod(np, p))); }
        }
        bound *= 2;
    }
}

        for n in [BigUint::from(11u32), BigUint::from(1234567891011121314u64), (BigUint::one()<<180usize)+BigUint::from(37u32)] {
            for count in [1,150,400,1800] {
                let a=factor_base(&n,count);let b=original(&n,count);
                assert_eq!(a.iter().map(|q|(q.p,q.root,q.log,q.inv,q.limit)).collect::<Vec<_>>(),b.iter().map(|q|(q.p,q.root,q.log,q.inv,q.limit)).collect::<Vec<_>>());
            }
        }
    }
    #[test]
    fn repaired_families_keep_exact_crt_contracts() {
        let n = BigUint::parse_bytes(b"13729273257095901138749901421978803731471830970253115531", 10).unwrap();
        let base = factor_base(&n, 1800);
        let half = 1 << 15;
        let target = (n.to_f64().unwrap() * 2.0).sqrt() / half as f64;
        for seed in 1..=64 {
            let old = choose_a_mode(&n, &base, half, seed, false);
            let new = choose_a_mode(&n, &base, half, seed, true);
            assert_eq!(old.qs.len(), new.qs.len());
            assert!(new.qs.windows(2).all(|w| w[0] < w[1]));
            assert!((new.a.to_f64().unwrap() / target).ln().abs() <= (old.a.to_f64().unwrap() / target).ln().abs() + 1e-12);
            assert_eq!(new.a, new.qs.iter().fold(BigUint::one(), |a, &i| a * base[i].p));
            for (slot, bl) in new.bs.iter().enumerate() {
                for (j, &i) in new.qs.iter().enumerate() {
                    if slot == j { assert_eq!((bl * bl) % base[i].p, &n % base[i].p); }
                    else { assert!((bl % base[i].p).is_zero()); }
                }
            }
            let mut b: BigInt = new.bs.iter().fold(BigInt::zero(), |acc, v| acc + BigInt::from(v.clone()));
            let a = BigInt::from(new.a.clone()); let ni = BigInt::from(n.clone());
            for poly in 0..(1usize << (new.qs.len() - 1)) {
                assert!(((&b * &b - &ni) % &a).is_zero());
                if poly + 1 < (1usize << (new.qs.len() - 1)) {
                    let v = (poly + 1).trailing_zeros() as usize;
                    if ((poly + 1) >> (v + 1)) & 1 == 0 { b -= BigInt::from(new.bs[v].clone()) * 2; }
                    else { b += BigInt::from(new.bs[v].clone()) * 2; }
                }
            }
        }
    }
    #[test]
    fn bounded_stride_updates_match_modulo_without_overflow() {
        for p in [1u32, 2, 3, 5, 7, 31, 65537, 1000003, u32::MAX - 4, u32::MAX] {
            for r in [0, p / 2, p - 1] {
                for delta in [0, p / 2, p - 1] {
                    assert_eq!(shift_residue(r, delta, p, true) as u64, (r as u64 + delta as u64) % p as u64);
                    assert_eq!(shift_residue(r, delta, p, false) as u64, (r as u64 + p as u64 - delta as u64) % p as u64);
                }
            }
            let mut state = 20261109u32;
            for _ in 0..10000 {
                state = state.wrapping_mul(1664525).wrapping_add(1013904223);
                let r = state % p;
                state = state.wrapping_mul(1664525).wrapping_add(1013904223);
                let delta = state % p;
                assert_eq!(shift_residue(r, delta, p, true) as u64, (r as u64 + delta as u64) % p as u64);
                assert_eq!(shift_residue(r, delta, p, false) as u64, (r as u64 + p as u64 - delta as u64) % p as u64);
            }
        }
    }

    #[test]
    fn the_divisibility_test_agrees_with_the_remainder() {
        for p in [3u64, 5, 7, 11, 65537, 1_000_003, 4_294_967_291] {
            let q = Prime::new(p, 1);
            assert_eq!(q.inv.wrapping_mul(p as u32), 1, "inverse mod 2^32 of {p}");
            let mut state = 20261109u32;
            for d in (0..5000u32).chain((0..5000).map(|_| { state = state.wrapping_mul(1664525).wrapping_add(1013904223); state })) {
                assert_eq!(q.divides(d), d as u64 % p == 0, "{d} by {p}");
            }
        }
    }
    #[test]
    fn small_prime_tests_and_rho() {
        assert!(is_prime_u64(1_000_000_007) && !is_prime_u64(1_000_000_007 * 3) && is_prime_u64(18446744073709551557));
        let n = 1_000_003u64 * 1_000_033;
        let f = rho_u64(n).unwrap(); assert!(n % f == 0 && f > 1 && f < n);
        // bound 5,000,000: the DLP ceiling is bound^1.8 ≈ 1.08e12
        assert_eq!(split_cofactor((1_000_003u128) * 1_000_033, 5_000_000, 100_000), Some(vec![(1_000_003, 1), (1_000_033, 1)]));
        assert_eq!(split_cofactor(1_000_003u128 * 1_000_003, 5_000_000, 100_000), Some(vec![(1_000_003, 2)]));
        assert_eq!(split_cofactor(100_003u128 * 6_000_011, 5_000_000, 100_000), None, "one prime over the bound");
        assert_eq!(split_cofactor(1_000_003u128 * 1_100_033, 1_000_000, 100_000), None, "over the ceiling");
    }
    #[test]
    fn square_large_cofactors_are_immediately_usable_and_do_not_join_odd_edges() {
        let rel = |large| Relation { y: BigUint::from(101u32), exps: vec![], large };
        let mut graph = Cycles::default();
        assert!(graph.add(rel(vec![(101, 1)])).is_none());
        let square = graph.add(rel(vec![(101, 2)])).expect("an even large exponent needs no partner");
        assert_eq!(square.len(), 1);
        assert_eq!(square[0].large, vec![(101, 2)]);
        let n = BigUint::from(1009u32);
        assert_eq!(&square[0].y * &square[0].y % &n, BigUint::from(101u32).pow(2) % &n);
        assert_eq!(graph.edges.len(), 1, "the square must not enter the forest");
        let pair = graph.add(rel(vec![(101, 1)])).unwrap();
        assert_eq!(pair.len(), 2);
        assert!(pair[0].times(&pair[1], &n).large.iter().all(|(_, e)| e % 2 == 0));
        assert!(graph.add(rel(vec![(103, 2), (107, 1)])).is_none());
        let pair = graph.add(rel(vec![(107, 1)])).unwrap();
        assert!(pair[0].times(&pair[1], &n).large.iter().all(|(_, e)| e % 2 == 0));
    }
    #[test]
    fn cycles_close_and_multiply_to_squares() {
        let n = BigUint::from(1_000_000_007u64);
        let rel = |l: Vec<(u64, u32)>| Relation { y: BigUint::from(3u32), exps: vec![(1, 1)], large: l };
        let mut g = Cycles::default();
        assert!(g.add(rel(vec![(101, 1)])).is_none());
        assert!(g.add(rel(vec![(101, 1), (103, 1)])).is_none());
        let c = g.add(rel(vec![(103, 1)])).expect("1–101–103–1 closes");
        assert_eq!(c.len(), 3);
        let mut acc = c[0].clone(); for r in &c[1..] { acc = acc.times(r, &n); }
        assert!(acc.large.iter().all(|(_, e)| e % 2 == 0), "{:?}", acc.large);
        assert!(g.add(rel(vec![(107, 1), (109, 1)])).is_none());
        assert!(g.add(rel(vec![(107, 1), (109, 1)])).is_some(), "a repeated pair is a two-cycle");
    }
    #[test]
    fn blocked_fill_equals_whole_fill() {
        let mut state = 12345u64;
        let mut rnd = |m: u64| { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state % m };
        for _ in 0..300 {
            let width = 50 + rnd(500) as usize;
            let p = 2 + rnd(60) as usize;
            let r0 = rnd(p as u64) as u32;
            let r1 = if rnd(5) == 0 { r0 } else { rnd(p as u64) as u32 };
            let seg = 1 + rnd(64) as usize;
            let mut whole = vec![0u8; width];
            let mut nx = [r0, r1];
            fill_block(&mut whole, 0, width, p, 3, &mut nx);
            let mut blocked = vec![0u8; width];
            let mut nb = [r0, r1];
            let mut lo = 0;
            while lo < width { let hi = (lo + seg).min(width); fill_block(&mut blocked[lo..hi], lo, hi, p, 3, &mut nb); lo = hi; }
            assert_eq!(whole, blocked, "p={p} roots={r0},{r1} width={width} seg={seg}");
            assert_eq!(nx, nb, "the carried positions agree too");
        }
    }
    #[test]
    fn pruning_keeps_every_dependency() {
        // every 4-by-3 binary matrix: the rows pruned are never part of any dependency
        for bits in 0..4096u32 {
            let odd: Vec<Vec<u32>> = (0..4).map(|r| (0..3).filter(|c| bits >> (3 * r + c) & 1 == 1).collect()).collect();
            let keep = prune_singletons(&odd, 3);
            let kept: std::collections::HashSet<usize> = keep.iter().copied().collect();
            for subset in 1..16u32 {
                let rows: Vec<usize> = (0..4).filter(|r| subset >> r & 1 == 1).collect();
                let mut sum = [0u8; 3];
                for &r in &rows { for &c in &odd[r] { sum[c as usize] ^= 1; } }
                if sum == [0, 0, 0] { assert!(rows.iter().all(|r| kept.contains(r)), "{bits:012b} {rows:?}"); }
            }
        }
        // a 15,000-row chain cascades to nothing in one pass
        let chain: Vec<Vec<u32>> = (0..15000u32).map(|i| vec![i, i + 1]).collect();
        assert!(prune_singletons(&chain, 15001).is_empty());
    }
    #[test]
    fn the_polynomials_are_consistent() {
        let n = BigUint::parse_bytes(b"4534149033143998129927133103", 10).unwrap(); // 25-digit-ish composite
        let base = factor_base(&n, 150);
        assert!(base.iter().skip(1).all(|q| (q.root * q.root) % q.p == (&n % q.p).to_u64().unwrap()));
        let fam = choose_a(&n, &base, 1 << 13, 7);
        // b² ≡ n (mod a) for the first b
        let b: BigUint = fam.bs.iter().fold(BigUint::zero(), |acc, x| acc + x);
        assert_eq!((&b * &b) % &fam.a, &n % &fam.a);
    }
    #[test]
    fn factors_a_few() {
        let (p, q) = (BigUint::from(206158430209u64), BigUint::from(21990232567889u64));
        let f = factor(&(&p * &q)).unwrap(); assert!(f == p || f == q, "25 digits");
        let (p, q) = (BigUint::from(26388279066671u64), BigUint::from(22517998136864981u64));
        let f = factor(&(&p * &q)).unwrap(); assert!(f == p || f == q, "30 digits");
    }
}

#[cfg(test)]
#[cfg(target_arch = "x86_64")]
mod scan_tests {
    use super::*;

    #[test]
    fn vector_scan_preserves_all_positions_order_and_tails() {
        if !is_x86_feature_detected!("avx2") { return; }
        let mut rng = crate::arith64::Lcg(20261113);
        for len in [0, 1, 7, 8, 31, 32, 33, 63, 127, 128, 129, 255, 256, 32768, 32769] {
            for mode in 0..4 {
                let mut buf: Vec<u8> = (0..len + 2).map(|_| match mode {
                    0 => 0x7f, 1 => 0x80, 2 => rng.next_u64() as u8,
                    _ => if rng.next_u64() % 1024 == 0 { 0x80 } else { 0x7f },
                }).collect();
                buf[0] = 0xff;
                buf[len + 1] = 0xff;
                let block = &buf[1..len + 1];  // intentionally unaligned
                for offset in [0, 32768] {
                    let mut expected = vec![999];
                    expected.extend(block.iter().enumerate().filter(|(_, v)| **v & 0x80 != 0)
                        .map(|(i, _)| offset + i));
                    let mut got = vec![999];
                    unsafe { scan_block_avx2(block, offset, &mut got); }
                    assert_eq!(got, expected, "len={len}, mode={mode}, offset={offset}");
                }
            }
        }
    }

    /// The sixteen-lane kernels produce the same entries in the same order as the scalar loops
    /// (and so as the AVX2 ones), on lengths that exercise the tails, with sentinel roots mixed in.
    #[test]
    fn avx512_kernels_match_scalar() {
        if !(is_x86_feature_detected!("avx512f") && is_x86_feature_detected!("avx512bw") && is_x86_feature_detected!("avx512vl")) { return; }
        let mut rng = crate::arith64::Lcg(2026_11_13);
        let width: u32 = 2 * 163840 + 1;
        for n in [0usize, 1, 15, 16, 17, 100, 1024, 1031] {
            let primes: Vec<u32> = (0..n).map(|_| 4096 + (rng.next_u64() % 3_000_000) as u32 | 1).collect();
            let roots: Vec<u32> = primes.iter().map(|&p| if rng.next_u64() % 23 == 0 { u32::MAX - 1 } else { (rng.next_u64() % (p as u64 * 2)) as u32 }).collect();
            let d: Vec<u32> = primes.iter().map(|&p| (rng.next_u64() % p as u64) as u32).collect();
            // shifts
            for add in [false, true] {
                let (mut a, mut b) = (roots.clone(), roots.clone());
                shift_all_body(&mut a, &d, &primes, add);
                unsafe { shift_all_avx512(&mut b, &d, &primes, add); }
                assert_eq!(a, b, "shift n={n} add={add}");
            }
            // single-stride compaction
            let mut want = Vec::new();
            for (i, &pos) in roots.iter().enumerate() { if pos < width { want.push(pos | (i as u32) << POS_BITS); } }
            let mut got = vec![0u32; n + 16];
            let cnt = unsafe { compact_hits_avx512(&roots, width, got.as_mut_ptr()) };
            got.truncate(cnt);
            assert_eq!(got, want, "compact n={n}");
            // k strides
            for k in [1u32, 2, 3, 7] {
                let mut want = Vec::new();
                for (i, (&r, &p)) in roots.iter().zip(&primes).enumerate() {
                    let mut pos = r;
                    for _ in 0..k { if pos < width { want.push(pos | (i as u32) << POS_BITS); } pos = pos.wrapping_add(p); }
                }
                let mut got = vec![0u32; n * k as usize + 16];
                let cnt = unsafe { compact_hits_k_avx512(&roots, &primes, k, width, got.as_mut_ptr()) };
                got.truncate(cnt);
                // the vector kernels emit a group's strides stride-major, the scalar loop root-major:
                // the same set, and only the set matters to the adds and the tagging
                got.sort_unstable(); want.sort_unstable();
                assert_eq!(got, want, "compact_k n={n} k={k}");
            }
            // one block of the mid walk, roots advanced
            let hi = 32768u32;
            let kb = 9;
            let mut ra = roots.iter().map(|&r| r % 40000).collect::<Vec<u32>>();
            let mut rb = ra.clone();
            let mut want = Vec::new();
            for (i, (r, &p)) in ra.iter_mut().zip(&primes).enumerate() { while *r < hi { want.push(*r | (i as u32) << POS_BITS); *r += p; } }
            let mut got = vec![0u32; n * kb as usize + 16];
            let cnt = unsafe { mid_block_avx512(&mut rb, &primes, hi, kb, got.as_mut_ptr()) };
            got.truncate(cnt);
            got.sort_unstable(); want.sort_unstable();
            assert_eq!(got, want, "mid_block n={n}");
            assert_eq!(ra, rb, "mid_block roots n={n}");
        }
    }
}

#[cfg(test)]
mod q2_tests {
    use super::*;

    #[test]
    fn q2_relations_preserve_the_exact_square_identity_in_both_integer_paths() {
        for (number, k, base_size, half) in [
            ("4533471823462034378103601", 1u32, 150, 8192),
            ("277033964293258511963625724654703158138482725303722794254981", 29, 5472, 32768),
            ("4993025603380020620960810892041035826574148743425607450811322494363051", 3, 8208, 32768),
        ] {
            let n = BigUint::parse_bytes(number.as_bytes(), 10).unwrap() * k;
            assert_eq!(&n % 8u32, BigUint::one());
            let base = factor_base(&n, base_size);
            let params = Params { base: base_size, half, large: 225, slack: 1.3, extra: 64 };
            let size = (half as f64 * (n.to_f64().unwrap() / 2.0).sqrt()).log2() - 1.0;
            let pmax = base.last().unwrap().p as f64;
            let threshold = (size - params.slack as f64 * pmax.log2() - (pmax * params.large as f64).log2()) as f32;
            for seed in [1, 7] {
                let fam = choose_a_form(&n, &base, half, seed, true, true);
                let a = BigInt::from(fam.a.clone());
                let mut b = fam.bs.iter().fold(BigInt::zero(), |acc, bl| acc + BigInt::from(bl.clone()));
                if b.is_even() { b += &a; }
                let four_a: BigInt = &a * 4;
                let ni = BigInt::from(n.clone());
                assert!(((&b * &b - &ni) % &four_a).is_zero());
                let c = (&b * &b - &ni) / &four_a;
                for x in [-(half as i64), -1, 0, 1, half as i64] {
                    let x = BigInt::from(x);
                    let y = &a * &x * 2 + &b;
                    assert_eq!(&y * &y - &ni, &four_a * ((&a * &x + &b) * &x + &c));
                }
                let (rels, _) = sieve_family(&n, &base, &fam, &params, threshold);
                assert!(!rels.is_empty(), "the fixture must exercise returned relations");
                for r in rels {
                    let mut rhs = BigUint::one();
                    let mut negative = false;
                    for (column, e) in r.exps {
                        if column == 0 { negative = e % 2 == 1; }
                        else { rhs = rhs * BigUint::from(base[column as usize - 1].p)
                            .modpow(&BigUint::from(e), &n) % &n; }
                    }
                    for (p, e) in r.large { rhs = rhs * BigUint::from(p).modpow(&BigUint::from(e), &n) % &n; }
                    if negative && !rhs.is_zero() { rhs = &n - rhs; }
                    assert_eq!(&r.y * &r.y % &n, rhs, "Q2 must restore the removed factor 4a");
                }
            }
        }
    }
}

#[cfg(test)]
mod param_tests {
    use super::*;
    #[test]
    fn exact_bits_pick_the_330_row() {
        // Parameter regression: the digit estimate gave 332 bits and a 5-block window here.
        let n = BigUint::parse_bytes(b"1662411017553075237233839631781649972371163397777594492538570355101974205404281730101312377184050959", 10).unwrap();
        assert_eq!(n.bits(), 330);
        let r = crate::params::for_bits(crate::params::YAFU_DIRECT_AVX512, n.bits() as u32);
        assert_eq!((r.base, r.blocks, r.lp_mult), (87500, 4, 350));
        let r2 = crate::params::for_bits(crate::params::YAFU_DIRECT_AVX512, 332);
        assert_eq!(r2.blocks, 5);
    }
}

#[cfg(test)]
mod divide_tests {
    use super::*;
    #[test]
    fn divide_out_matches_the_plain_loop() {
        let base: Vec<Prime> = [2u64, 3, 5, 7, 11, 13, 65537, 100003, 600011].iter().map(|&p| Prime::new(p, 0)).collect();
        let mut rng = crate::arith64::Lcg(7);
        for _ in 0..200 {
            let mut v = BigUint::from(rng.next_u64() | 1) * BigUint::from(rng.next_u64() | 1) * BigUint::from(rng.next_u64() | 1);
            let mut want: Vec<(u32, u32)> = vec![];
            let mut js = vec![];
            for (j, q) in base.iter().enumerate().skip(1) {
                let e = (rng.next_u64() % 4) as u32;
                if e > 0 || rng.next_u64() % 3 == 0 { js.push(j as u32); }
                for _ in 0..e { v *= q.p; }
            }
            let mut w = v.clone();
            for &j in &js { let p = base[j as usize].p; while (&w % p).is_zero() { w /= p; bump_exp(&mut want, j + 1); } }
            let mut got: Vec<(u32, u32)> = vec![];
            divide_out(&mut v, &js, &base, &mut got);
            assert_eq!(v, w);
            got.sort(); want.sort();
            assert_eq!(got, want);
        }
    }
}

#[cfg(test)]
mod exact_tests {
    use super::*;

    /// y² ≡ ± ∏ p^e · ∏ large^e (mod n): what every relation must satisfy.
    fn holds(r: &Relation, base: &[Prime], n: &BigUint) -> bool {
        let mut rhs = BigUint::one();
        let mut negative = false;
        for &(column, e) in &r.exps {
            if column == 0 { negative = e % 2 == 1; }
            else { rhs = rhs * BigUint::from(base[column as usize - 1].p).modpow(&BigUint::from(e), n) % n; }
        }
        for &(p, e) in &r.large { rhs = rhs * BigUint::from(p).modpow(&BigUint::from(e), n) % n; }
        if negative && !rhs.is_zero() { rhs = n - rhs; }
        &r.y * &r.y % n == rhs
    }

    #[test]
    fn small_part_counts_every_small_prime_with_its_multiplicity() {
        // random polynomials with signed coefficients of about 150 bits; the truth by big integers
        let primes: Vec<Prime> = [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263].iter().map(|&p| Prime::new(p, 1)).collect();
        let mut rng = crate::arith64::Lcg(20261113);
        let big = |rng: &mut crate::arith64::Lcg| -> BigInt {
            let m = BigUint::from(rng.next_u64()) << 96usize | BigUint::from(rng.next_u64()) << 32usize | BigUint::from(rng.next_u64() >> 32);
            if rng.next_u64() & 1 == 0 { BigInt::from(m) } else { -BigInt::from(m) }
        };
        let mut checked_powers = 0;
        for round in 0..40 {
            let a = big(&mut rng);
            let sp_starts = vec![[u32::MAX; 2]; primes.len()];   // no roots: every prime evaluated directly
            for _ in 0..4 {
                let (lin, c) = (big(&mut rng), big(&mut rng));
                // the discriminant belongs to the family; in this test every polynomial is its own family
                let mut sp = SmallPart::new(&primes, 256, &a, &(&lin * &lin - BigInt::from(4) * &a * &c));
                assert!(sp.two && sp.p.len() == 53, "the 53 odd primes below 256");
                sp.set_poly(&lin, &c);
                // positions spread over a 2^18 interval, plus ones forced onto high powers
                let mut xs: Vec<i64> = (0..300).map(|_| (rng.next_u64() % (1 << 18)) as i64 - (1 << 17)).collect();
                if round % 2 == 0 {
                    // make g(x0) ≡ 0 mod 3^5·5^3·2^7 by choosing c for a given x0: c = −(a·x0 + lin)·x0 + t·M
                    let x0 = xs[0];
                    let m = BigInt::from(243 * 125 * 128);
                    let c2 = -(&a * x0 + &lin) * x0 + &m * BigInt::from(rng.next_u64() | 1);
                    let mut sp = SmallPart::new(&primes, 256, &a, &(&lin * &lin - BigInt::from(4) * &a * &c2));
                    sp.set_poly(&lin, &c2);
                    let g = (&a * x0 + &lin) * x0 + &c2;
                    assert!((&g % &m).is_zero());
                    assert!(sp.bits(x0, 0, &sp_starts) >= 7.0 + 5.0 * 3f32.log2() + 3.0 * 5f32.log2() - 1e-3);
                    checked_powers += 1;
                    xs.truncate(40);
                    for &x in &xs { check(&sp, &a, &lin, &c2, x, &primes, &sp_starts); }
                    continue;
                }
                for &x in &xs { check(&sp, &a, &lin, &c, x, &primes, &sp_starts); }
            }
        }
        assert!(checked_powers >= 20);

        fn check(sp: &SmallPart, a: &BigInt, lin: &BigInt, c: &BigInt, x: i64, primes: &[Prime], starts: &[[u32; 2]]) {
            let g = (a * x + lin) * x + c;
            let mut want = 0f32;
            if !g.is_zero() {
                let mut v = g.magnitude().clone();
                let mut twos = 0u32;
                while (&v % 2u32).is_zero() { v /= 2u32; twos += 1; }
                want += twos.min(48) as f32;
                for q in primes.iter().filter(|q| q.p > 2 && q.p < 256) {
                    let mut cap = 0u32; { let mut m = q.p; while m * q.p <= 1 << 31 { m *= q.p; cap += 1; } cap += 1; }
                    let mut e = 0u32;
                    while (&v % q.p).is_zero() { v /= q.p; e += 1; }
                    want += e.min(cap) as f32 * (q.p as f32).log2();
                }
            } else { return; }
            let got = sp.bits(x, 0, starts);
            assert!((got - want).abs() < 1e-2, "x = {x}: got {got}, want {want}");
        }
    }

    #[test]
    fn exact_threshold_relations_hold_and_are_not_fewer() {
        // a 60-digit input through `sieve_family` both ways, same families. In this (debug) build
        // `SmallPart::bits` also checks its root test against direct evaluation at every report.
        let n0 = BigUint::parse_bytes(b"277033964293258511963625724654703158138482725303722794254981", 10).unwrap();
        let k = multiplier(&n0);
        let n = &n0 * k;
        assert!(n.bits() >= 180);
        let params = Params { base: 5472, half: 32768, large: 180, slack: 1.3, extra: 64 };
        let base = factor_base(&n, params.base);
        let size = (params.half as f64 * (n.to_f64().unwrap() / 2.0).sqrt()).log2();
        let pmax = base.last().unwrap().p as f64;
        let threshold = (size - params.slack as f64 * pmax.log2() - (pmax * params.large as f64).log2()) as f32;
        let (mut plain, mut exact) = (0usize, 0usize);
        for seed in [1u64, 2, 3] {
            let fam = choose_a(&n, &base, params.half, seed);
            EXACT_OVERRIDE.with(|o| o.set(Some(false)));
            let (a, _) = sieve_family(&n, &base, &fam, &params, threshold);
            EXACT_OVERRIDE.with(|o| o.set(Some(true)));
            let (b, _) = sieve_family(&n, &base, &fam, &params, threshold);
            EXACT_OVERRIDE.with(|o| o.set(None));
            assert!(!a.is_empty() && !b.is_empty());
            for r in a.iter().chain(&b) { assert!(holds(r, &base, &n), "a relation that does not hold"); }
            plain += a.len();
            exact += b.len();
        }
        eprintln!("60 digits, 3 families: {plain} relations with the plain threshold, {exact} with the exact one");
        assert!(exact * 100 >= plain * 95, "the exact threshold lost relations: {exact} against {plain}");
    }
}

#[cfg(test)]
mod vector_kernel_tests {
    use super::*;
    #[test]
    fn log2_fast_is_within_a_hundredth_of_a_bit() {
        let mut rng = crate::arith64::Lcg(99);
        for _ in 0..200_000 {
            let g = (rng.next_u64() >> 11) as f64 * 2f64.powi((rng.next_u64() % 400) as i32 - 100);
            if g < 1.0 { continue; }
            assert!((log2_fast(g) as f64 - g.log2()).abs() < 0.011, "{g}");
        }
        for g in [1.0f64, 2.0, 3.0, 1e30, 2f64.powi(185), 1.5 * 2f64.powi(60)] { assert!((log2_fast(g) as f64 - g.log2()).abs() < 0.011); }
    }
    #[test]
    fn word_sized_setup_arithmetic_equals_the_wide_arithmetic() {
        // the i64 extended Euclid and the 64-bit modular products against the 128-bit originals
        let mut rng = crate::arith64::Lcg(7);
        for p in [3u64, 5, 7, 65537, 1_000_003, 2_461_973, 4_294_967_291] {
            for _ in 0..2000 {
                let a = 1 + rng.next_u64() % (p - 1);
                let inv = modinv(a, p);
                assert_eq!((a as u128 * inv as u128 % p as u128) as u64, 1, "{a}^-1 mod {p}");
                let wide = { let (mut t, mut nt, mut r, mut nr) = (0i128, 1i128, p as i128, a as i128); while nr != 0 { let q = r / nr; (t, nt) = (nt, t - q * nt); (r, nr) = (nr, r - q * nr); } t.rem_euclid(p as i128) as u64 };
                assert_eq!(inv, wide);
                let x = rng.next_u64() % p;
                assert_eq!((2 * x % p) * inv % p, (2u128 * x as u128 * inv as u128 % p as u128) as u64);
            }
        }
        let big = BigUint::parse_bytes(b"1266849965389271886888480969488083293990717800033718809356376758077582106924382257211996775138403103", 10).unwrap();
        for p in [3u32, 65537, 2_461_973, 4_294_967_291] {
            assert_eq!(rem_u32(&big, p) as u64, (&big % p).to_u64().unwrap());
            let neg = -BigInt::from(big.clone());
            assert_eq!(rem_floor_u32(&neg, p) as u64, neg.mod_floor(&BigInt::from(p)).to_u64().unwrap());
            assert_eq!(rem_floor_u32(&BigInt::from(big.clone()), p) as u64, (&big % p).to_u64().unwrap());
        }
    }
    #[test]
    fn spread_sink_adds_the_same_bytes_inside_the_interval() {
        let mut rng = crate::arith64::Lcg(5);
        let width = 5000u32;
        let roots: Vec<u32> = (0..3000).map(|_| (rng.next_u64() % 9000) as u32).collect();
        let primes: Vec<u32> = (0..3000).map(|_| 600 + (rng.next_u64() % 900) as u32).collect();
        let logs: Vec<u8> = (0..3000).map(|_| 1 + (rng.next_u64() % 20) as u8).collect();
        for k in [1u32, 2, 4, 8] {
            let mut one = vec![0u8; width as usize + SINK];
            let mut many = vec![0u8; width as usize + SINK];
            if k == 1 { direct_add(&mut one, &roots, &logs, width, 0); direct_add(&mut many, &roots, &logs, width, SINK as u32 - 1); }
            else { direct_add_k(&mut one, &roots, &primes, &logs, width, k, 0); direct_add_k(&mut many, &roots, &primes, &logs, width, k, SINK as u32 - 1); }
            assert_eq!(one[..width as usize], many[..width as usize], "k = {k}");
            assert!(many[width as usize..].iter().filter(|&&b| b != 0).count() > 8, "the misses are spread");
        }
    }
    #[test]
    fn shift_all_avx2_equals_the_portable_loop() {
        let mut rng = crate::arith64::Lcg(11);
        for n in [0usize, 1, 7, 8, 9, 33, 1024] {
            let p: Vec<u32> = (0..n).map(|j| [8209u32, 65537, 2_461_973, u32::MAX][j % 4]).collect();
            let r: Vec<u32> = p.iter().map(|&p| if p == u32::MAX { u32::MAX - 1 } else { rng.next_u64() as u32 % p }).collect();
            let d: Vec<u32> = p.iter().map(|&p| if p == u32::MAX { 0 } else { rng.next_u64() as u32 % p }).collect();
            for add in [true, false] {
                let (mut a, mut b) = (r.clone(), r.clone());
                shift_all(&mut a, &d, &p, add, false);
                shift_all(&mut b, &d, &p, add, have_avx2());
                assert_eq!(a, b);
                // a skipped prime's sentinel (root u32::MAX − 1, delta 0, prime u32::MAX) may drift down by one
                // on a subtracting step (existing behaviour); what matters is that it stays out of any interval
                for j in 0..n { if p[j] != u32::MAX { assert_eq!(a[j], shift_residue(r[j], d[j], p[j], add)); } else { assert!(a[j] >= u32::MAX - 2); } }
            }
        }
    }
}

/// Timings of the direct-sieve kernels on arrays shaped like a 100-digit polynomial (5 blocks a
/// side, 90,432 primes): `cargo test --release bench_direct_kernels -- --ignored --nocapture`.
/// One thread; prints nanoseconds per element and microseconds per polynomial for each kernel.
#[cfg(test)]
mod vector_kernel_bench {
    use super::*;
    use std::time::Instant;

    #[test]
    #[ignore]
    fn bench_direct_kernels() {
        let width: u32 = 2 * 163840 + 1;
        let mut rng = crate::arith64::Lcg(2026);
        // base primes are about half of all primes: take every other prime from 8192 to 2.46 M
        let all = primes_upto(2_462_000);
        let direct: Vec<u32> = all.iter().filter(|&&p| p >= 8192).step_by(2).map(|&p| p as u32).collect();
        let nd = direct.len();
        let first_wide = direct.partition_point(|&p| p < width);
        let r0: Vec<u32> = direct.iter().map(|&p| rng.next_u64() as u32 % p).collect();
        let r1: Vec<u32> = direct.iter().map(|&p| rng.next_u64() as u32 % p).collect();
        let dl: Vec<u32> = direct.iter().map(|&p| rng.next_u64() as u32 % p).collect();
        let logs: Vec<u8> = direct.iter().map(|&p| ((p as f64).log2() + 0.5) as u8).collect();
        let avx2 = have_avx2();
        let reps = 300;
        println!("direct primes {nd} ({} mid below the interval, {} wide), interval {width}, AVX2 {avx2}", first_wide, nd - first_wide);
        let time = |name: &str, elems: usize, f: &mut dyn FnMut()| {
            f();
            let t = Instant::now();
            for _ in 0..reps { f(); }
            let ns = t.elapsed().as_nanos() as f64 / reps as f64;
            println!("  {name:58} {:8.1} us per polynomial, {:6.2} ns per element ({elems} elements)", ns / 1e3, ns / elems.max(1) as f64);
        };
        // root shifts, both arrays, by chunks of 1024 as the sieve does
        for (label, simd) in [("shift, baseline loop (SSE2)", false), ("shift, AVX2 loop", avx2)] {
            let (mut a, mut b) = (r0.clone(), r1.clone());
            let mut add = true;
            time(label, 2 * nd, &mut || {
                let mut lo = 0;
                while lo < nd { let hi = (lo + 1024).min(nd); shift_all(&mut a[lo..hi], &dl[lo..hi], &direct[lo..hi], add, simd); shift_all(&mut b[lo..hi], &dl[lo..hi], &direct[lo..hi], add, simd); lo = hi; }
                add = !add;
            });
        }
        // wide band: compaction, then the adds over the hits
        let mut acc = vec![0u8; width as usize + SINK + 8];
        let mut hpos: Vec<u32> = Vec::new();
        let nw = nd - first_wide;
        time("wide: compact hits (both roots)", 2 * nw, &mut || {
            hpos.clear();
            let mut lo = first_wide;
            while lo < nd { let hi = (lo + 1024).min(nd); compact_hits(&r0[lo..hi], width, &mut hpos, avx2); compact_hits(&r1[lo..hi], width, &mut hpos, avx2); lo = hi; }
        });
        let nh = hpos.len();
        time("wide: compact and add the hits, chunk by chunk", 2 * nw, &mut || {
            hpos.clear();
            let mut lo = first_wide;
            while lo < nd {
                let hi = (lo + 1024).min(nd);
                for r in [&r0[lo..hi], &r1[lo..hi]] { let start = hpos.len(); compact_hits(r, width, &mut hpos, avx2); add_hits(&mut acc[..width as usize], &hpos[start..], &logs[lo..hi]); }
                lo = hi;
            }
        });
        println!("    ({nh} wide hits per polynomial)");
        // mid bands: fixed strides, the old doubling bands and the fine ones, one sink byte and 64
        for (label, fine, mask) in [("mid: doubling bands, one sink byte", false, 0u32), ("mid: doubling bands, 64 sink bytes", false, 63), ("mid: fine bands, one sink byte", true, 0), ("mid: fine bands, 64 sink bytes", true, 63)] {
            let mut bands: Vec<(usize, usize, u32)> = Vec::new();
            let (mut hi, mut k) = (first_wide, 2u32);
            while hi > 0 {
                let pmin = (width as usize).div_ceil(k as usize) as u32;
                let lo_i = direct.partition_point(|&p| p < pmin).min(hi);
                if lo_i < hi { bands.push((lo_i, hi, k)); }
                hi = lo_i;
                k = if fine { if k < 8 { k + 1 } else { k + k / 4 } } else { k * 2 };
            }
            let strides: usize = bands.iter().map(|&(a, b, k)| 2 * (b - a) * k as usize).sum();
            let real: usize = (0..first_wide).map(|i| { let mut c = 0; for r in [r0[i], r1[i]] { let mut pos = r; while pos < width { c += 1; pos += direct[i]; } } c }).sum();
            time(label, strides, &mut || {
                for &(a, b, k) in &bands {
                    let mut lo = a;
                    while lo < b { let hi = (lo + 1024).min(b); direct_add_k(&mut acc[..width as usize + SINK], &r0[lo..hi], &direct[lo..hi], &logs[lo..hi], width, k, mask); direct_add_k(&mut acc[..width as usize + SINK], &r1[lo..hi], &direct[lo..hi], &logs[lo..hi], width, k, mask); lo = hi; }
                }
            });
            println!("    ({} bands, {strides} strides of which {real} land in the interval: {:.0} % wasted)", bands.len(), 100.0 * (strides - real) as f64 / strides as f64);
        }
        // the small-prime block fill for comparison: primes 256..8192 and 30..256, 32 KB blocks
        let small: Vec<u32> = all.iter().filter(|&&p| (30..8192).contains(&p)).step_by(2).map(|&p| p as u32).collect();
        for (label, from) in [("fill: primes 30..8192 in 32 KB blocks", 30u32), ("fill: primes 256..8192 in 32 KB blocks", 256), ("fill: primes 1024..8192 in 32 KB blocks", 1024)] {
            let ps: Vec<u32> = small.iter().copied().filter(|&p| p >= from).collect();
            let starts: Vec<[u32; 2]> = ps.iter().map(|&p| [rng.next_u64() as u32 % p, rng.next_u64() as u32 % p]).collect();
            let adds: usize = ps.iter().map(|&p| 2 * (width / p) as usize).sum();
            time(label, adds, &mut || {
                let mut next = starts.clone();
                let mut lo = 0usize;
                while lo < width as usize {
                    let hi = (lo + 32768).min(width as usize);
                    let block = &mut acc[lo..hi];
                    for (i, &p) in ps.iter().enumerate() { fill_block(block, lo, hi, p as usize, 11, &mut next[i]); }
                    lo = hi;
                }
            });
        }
        std::hint::black_box(&acc);
    }
}
