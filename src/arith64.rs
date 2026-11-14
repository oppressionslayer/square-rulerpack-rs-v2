//! 64-bit modular arithmetic in Montgomery form, for the cofactor work of the sieve: Fermat and
//! Miller-Rabin tests and the tiny ECM in `uecm`. Everything here is branch-light and allocation
//! free; a sieve calls it millions of times.

/// Montgomery context for an odd modulus `n < 2^63`.
#[derive(Clone, Copy, Debug)]
pub struct Mont {
    pub n: u64,
    /// -n^{-1} mod 2^64
    pub ninv: u64,
    /// R mod n with R = 2^64
    pub one: u64,
    /// R^2 mod n
    pub r2: u64,
}

impl Mont {
    pub fn new(n: u64) -> Mont {
        debug_assert!(n & 1 == 1 && n < (1 << 63));
        // Newton iteration for n^{-1} mod 2^64, then negate.
        let mut inv: u64 = 1;
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(n.wrapping_mul(inv)));
        }
        let ninv = inv.wrapping_neg();
        let one = ((1u128 << 64) % n as u128) as u64;
        let r2 = ((one as u128 * one as u128) % n as u128) as u64;
        Mont { n, ninv, one, r2 }
    }

    /// Montgomery product a·b·R^{-1} mod n, inputs below n.
    #[inline(always)]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        let t = a as u128 * b as u128;
        let m = (t as u64).wrapping_mul(self.ninv);
        let u = ((t + m as u128 * self.n as u128) >> 64) as u64;
        // u < 2n since a, b < n < 2^63; one conditional subtraction.
        if u >= self.n { u - self.n } else { u }
    }

    #[inline(always)]
    pub fn sqr(&self, a: u64) -> u64 { self.mul(a, a) }

    #[inline(always)]
    pub fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.n { s - self.n } else { s }
    }

    #[inline(always)]
    pub fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b { a - b } else { a + self.n - b }
    }

    /// Bring x (< n) into Montgomery form.
    #[inline(always)]
    pub fn to_mont(&self, x: u64) -> u64 { self.mul(x % self.n, self.r2) }

    /// Take x out of Montgomery form.
    #[inline(always)]
    pub fn from_mont(&self, x: u64) -> u64 { self.mul(x, 1) }

    /// a^e in Montgomery form, a already in Montgomery form.
    pub fn pow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut r = self.one;
        while e > 0 {
            if e & 1 == 1 { r = self.mul(r, a); }
            a = self.sqr(a);
            e >>= 1;
        }
        r
    }
}

/// Base-2 Fermat test for odd n < 2^63. A composite that passes is rare and harmless here: the
/// sieve only uses it to skip cofactors that are (probably) prime, as YAFU does.
pub fn fermat2(n: u64) -> bool {
    if n < 4 { return n >= 2; }
    if n & 1 == 0 { return false; }
    let m = Mont::new(n);
    let two = m.to_mont(2);
    m.pow(two, n - 1) == m.one
}

/// Deterministic Miller-Rabin for all n < 2^63 (bases 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
/// are enough below 3.3·10^24).
pub fn is_prime_u64(n: u64) -> bool {
    if n < 2 { return false; }
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n % p == 0 { return n == p; }
    }
    if n >= (1 << 63) { return is_prime_u64_slow(n); }
    let m = Mont::new(n);
    let nm1 = n - 1;
    let s = nm1.trailing_zeros();
    let d = nm1 >> s;
    let minus_one = m.sub(0, m.one);
    'base: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = m.pow(m.to_mont(a), d);
        if x == m.one || x == minus_one { continue; }
        for _ in 1..s {
            x = m.sqr(x);
            if x == minus_one { continue 'base; }
        }
        return false;
    }
    true
}

fn is_prime_u64_slow(n: u64) -> bool {
    fn mulmod(a: u64, b: u64, n: u64) -> u64 { ((a as u128 * b as u128) % n as u128) as u64 }
    fn powmod(mut a: u64, mut e: u64, n: u64) -> u64 {
        let mut r = 1u64;
        a %= n;
        while e > 0 { if e & 1 == 1 { r = mulmod(r, a, n); } a = mulmod(a, a, n); e >>= 1; }
        r
    }
    let nm1 = n - 1;
    let s = nm1.trailing_zeros();
    let d = nm1 >> s;
    'base: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(a, d, n);
        if x == 1 || x == nm1 { continue; }
        for _ in 1..s { x = mulmod(x, x, n); if x == nm1 { continue 'base; } }
        return false;
    }
    true
}

pub fn gcd(mut a: u64, mut b: u64) -> u64 {
    if a == 0 { return b; }
    if b == 0 { return a; }
    let shift = (a | b).trailing_zeros();
    a >>= a.trailing_zeros();
    loop {
        b >>= b.trailing_zeros();
        if a > b { std::mem::swap(&mut a, &mut b); }
        b -= a;
        if b == 0 { return a << shift; }
    }
}

/// A small xorshift generator; the sieve threads each carry one.
#[derive(Clone, Debug)]
pub struct Lcg(pub u64);
impl Lcg {
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        self.0 = x;
        x
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn montgomery_roundtrip_and_mul() {
        for &n in &[3u64, 101, 1_000_003, 4_294_967_291, (1u64 << 62) + 135] {
            let m = Mont::new(n);
            for &x in &[0u64, 1, 2, 12345, n - 1] {
                assert_eq!(m.from_mont(m.to_mont(x)), x % n);
            }
            let a = 987_654_321 % n; let b = 123_456_789 % n;
            let want = ((a as u128 * b as u128) % n as u128) as u64;
            assert_eq!(m.from_mont(m.mul(m.to_mont(a), m.to_mont(b))), want);
            assert_eq!(m.from_mont(m.pow(m.to_mont(3), 100)), {
                let mut r = 1u128; for _ in 0..100 { r = r * 3 % n as u128; } r as u64 });
        }
    }
    #[test]
    fn primality() {
        let primes = [2u64, 3, 5, 97, 7919, 1_000_000_007, 4_294_967_291, 2_305_843_009_213_693_951, 9_223_372_036_854_775_783];
        for &p in &primes { assert!(is_prime_u64(p), "{p}"); assert!(p == 2 || fermat2(p), "{p}"); }
        let composites = [1u64, 4, 91, 561, 1_000_000_007 * 3, 4_294_967_291 * 3, 3_215_031_751, 1_373_653 * 1_373_653];
        for &c in &composites { assert!(!is_prime_u64(c), "{c}"); }
        // 2047 = 23·89 is a base-2 Fermat liar: the test must say "probably prime" and MR must not.
        assert!(fermat2(2047)); assert!(!is_prime_u64(2047));
        let mut count = 0; for n in 1..10_000u64 { if is_prime_u64(n) { count += 1; } }
        assert_eq!(count, 1229);
    }
    #[test]
    fn gcd_works() {
        assert_eq!(gcd(0, 5), 5); assert_eq!(gcd(12, 18), 6); assert_eq!(gcd(1 << 40, 1 << 20), 1 << 20);
        assert_eq!(gcd(7919 * 104_729, 104_729 * 3), 104_729);
    }
}
