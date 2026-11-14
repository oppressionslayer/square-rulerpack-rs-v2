//! Splitting a 40-64-bit composite cofactor into two large primes: the double-large-prime step.
//!
//! Pollard P-1 first (cheap, catches a good share), then ECM on Montgomery curves with Suyama's
//! parametrisation, stage 1 by a prime-power ladder and a simple stage 2 over primes to B2. The
//! curve counts and bounds per bit size follow the measured choices in YAFU's microecm (Buhrow),
//! which were tuned for exactly this job; the arithmetic is `arith64::Mont`.
use crate::arith64::{gcd, Lcg, Mont};

/// Bounds for a cofactor of the given bit size: (B1, curves). B2 = 25·B1.
fn plan(bits: u32) -> (u64, u32) {
    match bits {
        0..=36 => (20, 32),
        37..=40 => (27, 32),
        41..=44 => (47, 32),
        45..=48 => (70, 32),
        49..=52 => (85, 32),
        53..=56 => (125, 48),
        57..=60 => (165, 80),
        _ => (205, 80),
    }
}

/// Small primes for the stage-1 ladder and stage 2.
fn small_primes(limit: u64) -> Vec<u64> {
    let mut sieve = vec![true; limit as usize + 1];
    let mut out = Vec::new();
    for p in 2..=limit as usize {
        if sieve[p] {
            out.push(p as u64);
            let mut q = p * p;
            while q <= limit as usize { sieve[q] = false; q += p; }
        }
    }
    out
}

thread_local! {
    static PRIMES: Vec<u64> = small_primes(205 * 25 + 1);
}

/// Pollard P-1 with B1 and B2 = 25·B1. Returns a proper factor or None.
pub fn pm1(n: u64, b1: u64) -> Option<u64> {
    let m = Mont::new(n);
    let mut a = m.to_mont(3);
    PRIMES.with(|primes| {
        for &p in primes.iter().take_while(|&&p| p <= b1) {
            let mut q = p;
            while q * p <= b1 { q *= p; }
            a = m.pow(a, q);
        }
        let g = gcd(m.sub(a, m.one), n);
        if g > 1 && g < n { return Some(g); }
        if g == n { return None; }
        // stage 2: multiply (a^p - 1) over primes in (B1, B2]; accumulate, gcd at the end
        let b2 = b1 * 25;
        let mut acc = m.one;
        for &p in primes.iter().skip_while(|&&p| p <= b1).take_while(|&&p| p <= b2) {
            let ap = m.pow(a, p);
            acc = m.mul(acc, m.sub(ap, m.one));
        }
        let g = gcd(m.from_mont(acc), n);
        if g > 1 && g < n { Some(g) } else { None }
    })
}

/// A point on a Montgomery curve By^2 = x^3 + Ax^2 + x in projective (X : Z), with the curve
/// constant stored as a24 = (A + 2)/4 in Montgomery form.
#[derive(Clone, Copy)]
struct Pt { x: u64, z: u64 }

#[inline(always)]
fn dbl(m: &Mont, a24: u64, p: Pt) -> Pt {
    let s = m.add(p.x, p.z);
    let d = m.sub(p.x, p.z);
    let s2 = m.sqr(s);
    let d2 = m.sqr(d);
    let t = m.sub(s2, d2);
    Pt { x: m.mul(s2, d2), z: m.mul(t, m.add(d2, m.mul(a24, t))) }
}

/// Differential addition: P + Q given P - Q.
#[inline(always)]
fn dadd(m: &Mont, p: Pt, q: Pt, diff: Pt) -> Pt {
    let u = m.mul(m.sub(p.x, p.z), m.add(q.x, q.z));
    let v = m.mul(m.add(p.x, p.z), m.sub(q.x, q.z));
    let s = m.add(u, v);
    let d = m.sub(u, v);
    Pt { x: m.mul(diff.z, m.sqr(s)), z: m.mul(diff.x, m.sqr(d)) }
}

/// Montgomery ladder: k·P.
fn ladder(m: &Mont, a24: u64, p: Pt, k: u64) -> Pt {
    if k == 0 { return Pt { x: m.one, z: 0 }; }
    let mut r0 = p;
    let mut r1 = dbl(m, a24, p);
    let bits = 63 - k.leading_zeros();
    for i in (0..bits).rev() {
        if (k >> i) & 1 == 1 {
            r0 = dadd(m, r1, r0, p);
            r1 = dbl(m, a24, r1);
        } else {
            r1 = dadd(m, r0, r1, p);
            r0 = dbl(m, a24, r0);
        }
    }
    r0
}

/// Suyama's parametrisation from sigma: gives a24 and a starting point, both in Montgomery
/// form, or None if an inversion fails (which itself reveals a factor, returned via Err).
fn suyama(m: &Mont, sigma: u64) -> Result<(u64, Pt), u64> {
    let s = m.to_mont(sigma);
    let u = m.sub(m.sqr(s), m.to_mont(5));
    let v = m.mul(m.to_mont(4), s);
    let u3 = m.mul(m.sqr(u), u);
    let v3 = m.mul(m.sqr(v), v);
    // a24 = (v - u)^3 (3u + v) / (16 u^3 v)
    let vmu = m.sub(v, u);
    let num = m.mul(m.mul(m.sqr(vmu), vmu), m.add(m.mul(m.to_mont(3), u), v));
    let den = m.mul(m.to_mont(16), m.mul(u3, v));
    let den_n = m.from_mont(den);
    let g = gcd(den_n, m.n);
    if g != 1 { return Err(g); }
    let inv = m.pow(den, m.n - 2); // n is treated as prime for the inverse; a wrong inverse only
                                   // happens when gcd != 1, caught above for the common case
    // but n is composite here, so use the extended Euclid instead
    let inv = mont_inverse(m, den).unwrap_or(inv);
    let a24 = m.mul(num, inv);
    Ok((a24, Pt { x: u3, z: v3 }))
}

/// Modular inverse of a Montgomery-form value, via extended Euclid on the plain value.
fn mont_inverse(m: &Mont, a: u64) -> Option<u64> {
    let a_plain = m.from_mont(a);
    let (mut old_r, mut r) = (a_plain as i128, m.n as i128);
    let (mut old_s, mut s) = (1i128, 0i128);
    while r != 0 {
        let q = old_r / r;
        (old_r, r) = (r, old_r - q * r);
        (old_s, s) = (s, old_s - q * s);
    }
    if old_r != 1 { return None; }
    let inv = old_s.rem_euclid(m.n as i128) as u64;
    Some(m.to_mont(inv))
}

/// One ECM curve: stage 1 to B1 by the prime-power ladder, stage 2 to 25·B1 by baby steps and
/// giant steps on the x-only curve. Returns a proper factor or None.
fn ecm_curve(m: &Mont, sigma: u64, b1: u64) -> Option<u64> {
    let (a24, mut p) = match suyama(m, sigma) {
        Ok(v) => v,
        Err(g) => return if g < m.n { Some(g) } else { None },
    };
    PRIMES.with(|primes| {
        for &q in primes.iter().take_while(|&&q| q <= b1) {
            let mut e = q;
            while e * q <= b1 { e *= q; }
            p = ladder(m, a24, p, e);
        }
        let g = gcd(m.from_mont(p.z), m.n);
        if g > 1 && g < m.n { return Some(g); }
        if g == m.n { return None; }
        // stage 2: every prime q in (B1, B2] is m·D ± j with D even and j odd below D/2, so
        // q·P = O mod p exactly when (m·D)·P and j·P share an x-coordinate mod p. Baby steps:
        // T[j] = j·P for odd j; giant steps: G_m = (m·D)·P. One cross product per prime.
        let b2 = b1 * 25;
        let d = ((b2 as f64).sqrt() as u64 * 2 + 1) & !1;
        let half = d / 2;
        let p2 = dbl(m, a24, p);
        let mut table: Vec<Pt> = Vec::with_capacity((half / 2 + 1) as usize);
        table.push(p);                       // j = 1
        if half > 1 { table.push(dadd(m, p, p2, p)); } // j = 3: P + 2P with difference P
        let mut j = 5;
        while j < half { let k = table.len(); table.push(dadd(m, table[k - 1], p2, table[k - 2])); j += 2; }
        let pd = ladder(m, a24, p, d);
        let m_last = b2.div_ceil(d) + 1;
        let mut giants: Vec<Pt> = Vec::with_capacity(m_last as usize + 1);
        giants.push(Pt { x: m.one, z: 0 }); // m = 0, unused
        giants.push(pd);
        if m_last >= 2 { giants.push(dbl(m, a24, pd)); }
        for _ in 3..=m_last { let k = giants.len(); giants.push(dadd(m, giants[k - 1], pd, giants[k - 2])); }
        let mut acc = m.one;
        for &q in primes.iter().skip_while(|&&q| q <= b1).take_while(|&&q| q <= b2) {
            let mi = (q + half) / d;
            let jj = q.abs_diff(mi * d);
            if jj == 0 || jj % 2 == 0 { continue; } // cannot happen for odd q and even d, kept for safety
            let g = giants[mi as usize];
            let t = table[(jj / 2) as usize];
            let cross = m.sub(m.mul(g.x, t.z), m.mul(t.x, g.z));
            acc = m.mul(acc, cross);
        }
        let g = gcd(m.from_mont(acc), m.n);
        if g > 1 && g < m.n { Some(g) } else { None }
    })
}

/// Split an odd composite n < 2^63 that has no factor below the sieve's factor-base bound.
/// Returns a proper factor or None after the planned effort. `rng` seeds the curves so retries
/// on the same n take different curves.
pub fn split(n: u64, rng: &mut Lcg) -> Option<u64> {
    if n < 4 { return None; }
    if n & 1 == 0 { return Some(2); }
    let bits = 64 - n.leading_zeros();
    let (b1, curves) = plan(bits);
    if bits > 54 {
        if let Some(f) = pm1(n, 333) { return Some(f); }
    }
    let m = Mont::new(n);
    for _ in 0..curves {
        let sigma = 6 + rng.next_u64() % (n - 7).max(1);
        if let Some(f) = ecm_curve(&m, sigma, b1) { return Some(f); }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arith64::is_prime_u64;

    fn semiprime(rng: &mut Lcg, bits_each: u32) -> (u64, u64, u64) {
        let next_prime = |mut x: u64| { x |= 1; while !is_prime_u64(x) { x += 2; } x };
        let p = next_prime((rng.next_u64() >> (64 - bits_each)) | (1 << (bits_each - 1)));
        let q = next_prime((rng.next_u64() >> (64 - bits_each)) | (1 << (bits_each - 1)));
        (p, q, p * q)
    }

    #[test]
    fn pm1_finds_smooth_order() {
        // p = 51481: p - 1 = 51480 = 2^3·3^2·5·11·13 is 19-smooth; q = 1_000_003
        assert!(is_prime_u64(51481));
        let n = 51481u64 * 1_000_003;
        let f = pm1(n, 20).unwrap();
        assert!(n % f == 0 && f > 1 && f < n);
    }

    #[test]
    fn splits_random_semiprimes_up_to_62_bits() {
        let mut rng = Lcg(0x9E37_79B9_7F4A_7C15);
        for &bits in &[20u32, 24, 26, 28, 30, 31] {
            let mut ok = 0;
            let trials = 60;
            for _ in 0..trials {
                let (p, q, n) = semiprime(&mut rng, bits);
                if p == q { ok += 1; continue; }
                if let Some(f) = split(n, &mut rng) {
                    assert!(f == p || f == q, "n={n} f={f}");
                    ok += 1;
                }
            }
            // the plan is sized for a high hit rate; demand at least 90 % here
            assert!(ok * 10 >= trials * 9, "{bits}-bit halves: {ok}/{trials}");
        }
    }

    #[test]
    fn never_returns_a_non_factor() {
        let mut rng = Lcg(12345);
        for _ in 0..200 {
            let (_, _, n) = semiprime(&mut rng, 22);
            if let Some(f) = split(n, &mut rng) { assert_eq!(n % f, 0); assert!(f > 1 && f < n); }
        }
    }
}
