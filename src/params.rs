//! Sieve parameters by input size. Three tables: Lars's rows measured in the rulerpack-rs sweeps
//! (single large prime), and the two YAFU 3.x tables for its direct sieve (plain/AVX2 and
//! AVX-512), read from `siqs_params.c` in Ben Buhrow's public-domain source. The YAFU rows are a
//! starting point for the double-large-prime engine; every row that ends up in the default table
//! has to win a race here first.

/// One row of a parameter table.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Row {
    /// input size in bits
    pub bits: u32,
    /// factor-base primes
    pub base: u32,
    /// single large prime bound as a multiple of the largest base prime
    pub lp_mult: u32,
    /// 32 KB sieve blocks per side
    pub blocks: u32,
    /// trial-division cutoff in bits (0 = formula)
    pub tf_bits: u32,
    /// double-large-prime exponent: cofactor bound = (lp_mult * pmax)^dlp_exp (0 = single)
    pub dlp_exp: f32,
}

/// YAFU `siqs_table_direct` (plain and AVX2 builds), 60-110 digits.
pub const YAFU_DIRECT: &[Row] = &[
    Row { bits: 198, base: 5472, lp_mult: 180, blocks: 1, tf_bits: 66, dlp_exp: 0.0 },
    Row { bits: 215, base: 10096, lp_mult: 225, blocks: 1, tf_bits: 73, dlp_exp: 1.62 },
    Row { bits: 232, base: 8208, lp_mult: 225, blocks: 1, tf_bits: 68, dlp_exp: 1.62 },
    Row { bits: 248, base: 8416, lp_mult: 300, blocks: 1, tf_bits: 72, dlp_exp: 1.72 },
    Row { bits: 265, base: 11856, lp_mult: 300, blocks: 1, tf_bits: 77, dlp_exp: 1.76 },
    Row { bits: 281, base: 16128, lp_mult: 300, blocks: 2, tf_bits: 85, dlp_exp: 1.74 },
    Row { bits: 298, base: 35232, lp_mult: 350, blocks: 2, tf_bits: 89, dlp_exp: 1.74 },
    Row { bits: 310, base: 56071, lp_mult: 350, blocks: 3, tf_bits: 94, dlp_exp: 1.74 },
    Row { bits: 320, base: 71745, lp_mult: 350, blocks: 3, tf_bits: 97, dlp_exp: 1.78 },
    Row { bits: 330, base: 87500, lp_mult: 350, blocks: 4, tf_bits: 102, dlp_exp: 1.8 },
    Row { bits: 340, base: 102120, lp_mult: 350, blocks: 6, tf_bits: 104, dlp_exp: 1.8 },
    Row { bits: 350, base: 119102, lp_mult: 350, blocks: 8, tf_bits: 106, dlp_exp: 1.8 },
    Row { bits: 360, base: 132742, lp_mult: 350, blocks: 8, tf_bits: 111, dlp_exp: 1.8 },
    Row { bits: 370, base: 239500, lp_mult: 150, blocks: 16, tf_bits: 0, dlp_exp: 1.9 },
];

/// YAFU `siqs_table_direct_avx512`, 60-110 digits.
pub const YAFU_DIRECT_AVX512: &[Row] = &[
    Row { bits: 198, base: 4912, lp_mult: 120, blocks: 1, tf_bits: 64, dlp_exp: 0.0 },
    Row { bits: 215, base: 6512, lp_mult: 300, blocks: 1, tf_bits: 69, dlp_exp: 1.70 },
    Row { bits: 232, base: 9920, lp_mult: 200, blocks: 1, tf_bits: 69, dlp_exp: 1.70 },
    Row { bits: 248, base: 11208, lp_mult: 225, blocks: 2, tf_bits: 74, dlp_exp: 1.70 },
    Row { bits: 265, base: 24624, lp_mult: 250, blocks: 2, tf_bits: 79, dlp_exp: 1.72 },
    Row { bits: 281, base: 32256, lp_mult: 250, blocks: 3, tf_bits: 85, dlp_exp: 1.74 },
    Row { bits: 298, base: 43000, lp_mult: 250, blocks: 4, tf_bits: 89, dlp_exp: 1.74 },
    Row { bits: 310, base: 56071, lp_mult: 250, blocks: 4, tf_bits: 94, dlp_exp: 1.74 },
    Row { bits: 320, base: 71745, lp_mult: 300, blocks: 4, tf_bits: 97, dlp_exp: 1.78 },
    Row { bits: 330, base: 87500, lp_mult: 350, blocks: 4, tf_bits: 102, dlp_exp: 1.8 },
    Row { bits: 340, base: 102120, lp_mult: 350, blocks: 6, tf_bits: 104, dlp_exp: 1.8 },
    Row { bits: 350, base: 119102, lp_mult: 350, blocks: 8, tf_bits: 106, dlp_exp: 1.8 },
    Row { bits: 360, base: 132742, lp_mult: 350, blocks: 8, tf_bits: 111, dlp_exp: 1.8 },
    Row { bits: 370, base: 239500, lp_mult: 150, blocks: 16, tf_bits: 0, dlp_exp: 1.9 },
];

/// Lars's rows from the rulerpack-rs 0.3.0 sweeps (single large prime, bucket sieve). The
/// interval there is a half-width 2^16 or 2^17, written here as 32 KB blocks per side: 2^17 half
/// = 4 blocks. lp_mult is the "large" parameter; tf_bits 0 means the slack-based threshold.
pub const RULERPACK_030: &[Row] = &[
    Row { bits: 226, base: 15000, lp_mult: 100, blocks: 2, tf_bits: 0, dlp_exp: 0.0 },  // 68-72 d
    Row { bits: 243, base: 30000, lp_mult: 120, blocks: 4, tf_bits: 0, dlp_exp: 0.0 },  // 73-77
    Row { bits: 259, base: 50000, lp_mult: 100, blocks: 4, tf_bits: 0, dlp_exp: 0.0 },  // 78-84
    Row { bits: 283, base: 80000, lp_mult: 200, blocks: 4, tf_bits: 0, dlp_exp: 0.0 },  // 85-92
    Row { bits: 309, base: 120000, lp_mult: 250, blocks: 4, tf_bits: 0, dlp_exp: 0.0 }, // 93-97
    Row { bits: 326, base: 200000, lp_mult: 300, blocks: 4, tf_bits: 0, dlp_exp: 0.0 }, // 98+
];

/// Parameters for an input of `bits` bits from a table, interpolating between rows the way
/// YAFU's `get_params` does: base linear in bits; lp_mult and blocks the rounded average of the
/// two bracketing rows; tf_bits and dlp_exp linear when both rows have them. Below the first row
/// the first row is used; above the last, the last row's slope continues for the base.
pub fn for_bits(table: &[Row], bits: u32) -> Row {
    let first = table[0];
    let last = table[table.len() - 1];
    if bits <= first.bits { return Row { bits, ..first }; }
    if bits >= last.bits {
        let prev = table[table.len() - 2];
        let slope = (last.base as f64 - prev.base as f64) / (last.bits as f64 - prev.bits as f64);
        let base = (last.base as f64 + slope * (bits - last.bits) as f64).round() as u32;
        return Row { bits, base: round16(base), ..last };
    }
    if let Some(r) = table.iter().find(|r| r.bits == bits) { return *r; }
    let i = table.iter().position(|r| r.bits > bits).unwrap();
    let (lo, hi) = (table[i - 1], table[i]);
    let t = (bits - lo.bits) as f64 / (hi.bits - lo.bits) as f64;
    let lerp = |a: f64, b: f64| a + (b - a) * t;
    Row {
        bits,
        base: round16(lerp(lo.base as f64, hi.base as f64).round() as u32),
        lp_mult: ((lo.lp_mult + hi.lp_mult) as f64 / 2.0).round() as u32,
        blocks: ((lo.blocks + hi.blocks) as f64 / 2.0).round() as u32,
        tf_bits: if lo.tf_bits > 0 && hi.tf_bits > 0 { lerp(lo.tf_bits as f64, hi.tf_bits as f64).round() as u32 } else { 0 },
        dlp_exp: if lo.dlp_exp > 0.0 && hi.dlp_exp > 0.0 { lerp(lo.dlp_exp as f64, hi.dlp_exp as f64) as f32 } else { 0.0 },
    }
}

fn round16(b: u32) -> u32 { b.div_ceil(16) * 16 }

/// The decimal digit count of an n-bit number, as the tables' comments use it.
pub fn digits_of_bits(bits: u32) -> u32 { (bits as f64 * std::f64::consts::LOG10_2).ceil() as u32 }

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_rows_come_back() {
        let r = for_bits(YAFU_DIRECT, 330);
        assert_eq!((r.base, r.lp_mult, r.blocks, r.tf_bits), (87500, 350, 4, 102));
        assert!((r.dlp_exp - 1.8).abs() < 1e-6);
    }
    #[test]
    fn interpolation_matches_lars_100_digit_log() {
        // The 4886...7869 input is 332 bits. YAFU selects its table row before multiplying
        // by 29; its later "102 digits / 337 bits" diagnostic describes kN, not this lookup.
        // YAFU's final factor-base padding adds another 16 entries in the saved log.
        let r = for_bits(YAFU_DIRECT_AVX512, 332);
        assert_eq!((r.base, r.lp_mult, r.blocks, r.tf_bits), (90432, 350, 5, 102));
    }
    #[test]
    fn edges() {
        assert_eq!(for_bits(YAFU_DIRECT, 100).base, 5472);
        let r = for_bits(YAFU_DIRECT, 380);
        assert!(r.base > 239500);
        assert_eq!(digits_of_bits(337), 102);
        assert_eq!(digits_of_bits(265), 80);
    }
}
