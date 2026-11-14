# Square Rulerpack

An optimized self-initializing quadratic sieve in Rust by Lars Rocha. Uses
Rulerpack for modular arithmetic, primality tests and matrix dependencies.

This factors a 100 digit number in under 700 seconds on an Vendor ID:                  
GenuineIntel
  Model name:                11th Gen Intel(R) Core(TM) i9-11900H @ 2.50GHz


```sh
curl -L -o sq.crate https://crates.io/api/v1/crates/square-rulerpack/0.2.0/download
tar xzf sq.crate && cd square-rulerpack-0.2.0
SIQS_Q2=1 SIQS_SMALL=512 RAYON_NUM_THREADS=8 cargo run --release --example factor -- 4886189752164427378644372029277135326505927538331605193174055598646955377478099751972989945355617869


     Running `target/release/examples/factor 4886189752164427378644372029277135326505927538331605193174055598646955377478099751972989945355617869`
  siqs: input 100 digits / 332 bits; kN 102 digits / 337 bits; multiplier 29
  base 90432 half 163840 (5 blocks/side) large 350 slack 1.3 extra 64; 8 threads; AVX2; Q2
  88/90497 relations, 7744 partial edges, 36864 polynomials, 8 s elapsed, about 7801 s to go
  226/90497 relations, 21612 partial edges, 102400 polynomials, 15 s elapsed, about 5997 s to go
  374/90497 relations, 35397 partial edges, 167936 polynomials, 22 s elapsed, about 5419 s to go
  491/90497 relations, 48486 partial edges, 229376 polynomials, 28 s elapsed, about 3942 s to go
  585/90497 relations, 59910 partial edges, 282624 polynomials, 33 s elapsed, about 4375 s to go
  768/90497 relations, 77280 partial edges, 364544 polynomials, 38 s elapsed, about 3533 s to go
  1023/90497 relations, 98668 partial edges, 464896 polynomials, 44 s elapsed, about 2918 s to go
  1264/90497 relations, 119396 partial edges, 563200 polynomials, 49 s elapsed, about 2053 s to go
  1531/90497 relations, 140461 partial edges, 661504 polynomials, 55 s elapsed, about 1953 s to go
  1750/90497 relations, 161133 partial edges, 759808 polynomials, 60 s elapsed, about 2040 s to go
  1998/90497 relations, 182025 partial edges, 858112 polynomials, 66 s elapsed, about 2021 s to go
  2254/90497 relations, 203009 partial edges, 956416 polynomials, 72 s elapsed, about 1938 s to go
  2514/90497 relations, 223954 partial edges, 1054720 polynomials, 77 s elapsed, about 1896 s to go
  2778/90497 relations, 244666 partial edges, 1153024 polynomials, 83 s elapsed, about 1898 s to go
  2955/90497 relations, 259334 partial edges, 1222656 polynomials, 88 s elapsed, about 2048 s to go
  3277/90497 relations, 285520 partial edges, 1345536 polynomials, 94 s elapsed, about 1436 s to go
  3493/90497 relations, 303730 partial edges, 1431552 polynomials, 99 s elapsed, about 1667 s to go
  3705/90497 relations, 320985 partial edges, 1513472 polynomials, 104 s elapsed, about 1783 s to go
  3967/90497 relations, 340851 partial edges, 1607680 polynomials, 109 s elapsed, about 1760 s to go
  4262/90497 relations, 361806 partial edges, 1705984 polynomials, 115 s elapsed, about 1654 s to go
  4470/90497 relations, 378132 partial edges, 1783808 polynomials, 121 s elapsed, about 2024 s to go
  4813/90497 relations, 402273 partial edges, 1898496 polynomials, 126 s elapsed, about 1712 s to go
  5102/90497 relations, 422925 partial edges, 1996800 polynomials, 131 s elapsed, about 1694 s to go
  5394/90497 relations, 442942 partial edges, 2091008 polynomials, 136 s elapsed, about 1463 s to go
  5635/90497 relations, 460403 partial edges, 2172928 polynomials, 142 s elapsed, about 1617 s to go
  5964/90497 relations, 480353 partial edges, 2267136 polynomials, 147 s elapsed, about 1515 s to go
  6243/90497 relations, 498333 partial edges, 2353152 polynomials, 152 s elapsed, about 1528 s to go
  6538/90497 relations, 518082 partial edges, 2447360 polynomials, 157 s elapsed, about 1478 s to go
  6735/90497 relations, 530965 partial edges, 2508800 polynomials, 163 s elapsed, about 1775 s to go
  7063/90497 relations, 550191 partial edges, 2598912 polynomials, 168 s elapsed, about 1591 s to go
  7470/90497 relations, 574573 partial edges, 2713600 polynomials, 173 s elapsed, about 1409 s to go
  7825/90497 relations, 594729 partial edges, 2807808 polynomials, 178 s elapsed, about 1233 s to go
  8167/90497 relations, 613905 partial edges, 2897920 polynomials, 185 s elapsed, about 1392 s to go
  8554/90497 relations, 637371 partial edges, 3008512 polynomials, 190 s elapsed, about 1280 s to go
  8923/90497 relations, 657202 partial edges, 3102720 polynomials, 195 s elapsed, about 1242 s to go
  9271/90497 relations, 676564 partial edges, 3192832 polynomials, 200 s elapsed, about 1202 s to go
  9680/90497 relations, 697216 partial edges, 3291136 polynomials, 207 s elapsed, about 1289 s to go
  10144/90497 relations, 720692 partial edges, 3401728 polynomials, 212 s elapsed, about 1127 s to go
  10532/90497 relations, 739798 partial edges, 3491840 polynomials, 217 s elapsed, about 1111 s to go
  10920/90497 relations, 759763 partial edges, 3586048 polynomials, 223 s elapsed, about 1088 s to go
  11263/90497 relations, 777078 partial edges, 3667968 polynomials, 228 s elapsed, about 1119 s to go
  11641/90497 relations, 796283 partial edges, 3758080 polynomials, 233 s elapsed, about 1097 s to go
  12080/90497 relations, 815588 partial edges, 3848192 polynomials, 238 s elapsed, about 1043 s to go
  12534/90497 relations, 835625 partial edges, 3942400 polynomials, 245 s elapsed, about 1261 s to go
  13104/90497 relations, 861536 partial edges, 4065280 polynomials, 251 s elapsed, about 952 s to go
  13533/90497 relations, 879732 partial edges, 4151296 polynomials, 256 s elapsed, about 936 s to go
  13915/90497 relations, 896383 partial edges, 4229120 polynomials, 261 s elapsed, about 978 s to go
  14467/90497 relations, 918825 partial edges, 4335616 polynomials, 267 s elapsed, about 720 s to go
  14978/90497 relations, 939513 partial edges, 4433920 polynomials, 274 s elapsed, about 873 s to go
  15605/90497 relations, 964747 partial edges, 4552704 polynomials, 280 s elapsed, about 833 s to go
  16248/90497 relations, 987334 partial edges, 4659200 polynomials, 285 s elapsed, about 761 s to go
  16795/90497 relations, 1007227 partial edges, 4753408 polynomials, 291 s elapsed, about 715 s to go
  17393/90497 relations, 1027369 partial edges, 4847616 polynomials, 296 s elapsed, about 664 s to go
  17863/90497 relations, 1044966 partial edges, 4929536 polynomials, 301 s elapsed, about 702 s to go
  18411/90497 relations, 1064177 partial edges, 5019648 polynomials, 306 s elapsed, about 697 s to go
  18990/90497 relations, 1083332 partial edges, 5109760 polynomials, 311 s elapsed, about 639 s to go
  19551/90497 relations, 1103172 partial edges, 5203968 polynomials, 316 s elapsed, about 637 s to go
  20141/90497 relations, 1120431 partial edges, 5285888 polynomials, 322 s elapsed, about 628 s to go
  20634/90497 relations, 1136145 partial edges, 5359616 polynomials, 327 s elapsed, about 664 s to go
  21443/90497 relations, 1160445 partial edges, 5474304 polynomials, 332 s elapsed, about 431 s to go
  22079/90497 relations, 1179493 partial edges, 5564416 polynomials, 338 s elapsed, about 483 s to go
  22685/90497 relations, 1198596 partial edges, 5654528 polynomials, 343 s elapsed, about 509 s to go
  23129/90497 relations, 1211550 partial edges, 5715968 polynomials, 348 s elapsed, about 559 s to go
  24012/90497 relations, 1235744 partial edges, 5830656 polynomials, 353 s elapsed, about 381 s to go
  24739/90497 relations, 1255687 partial edges, 5924864 polynomials, 358 s elapsed, about 425 s to go
  25447/90497 relations, 1274810 partial edges, 6014976 polynomials, 364 s elapsed, about 441 s to go
  26092/90497 relations, 1292163 partial edges, 6096896 polynomials, 369 s elapsed, about 452 s to go
  26737/90497 relations, 1308009 partial edges, 6170624 polynomials, 375 s elapsed, about 570 s to go
  27678/90497 relations, 1331936 partial edges, 6285312 polynomials, 380 s elapsed, about 435 s to go
  28415/90497 relations, 1350987 partial edges, 6375424 polynomials, 385 s elapsed, about 430 s to go
  29127/90497 relations, 1369161 partial edges, 6461440 polynomials, 390 s elapsed, about 426 s to go
  29908/90497 relations, 1387339 partial edges, 6547456 polynomials, 395 s elapsed, about 405 s to go
  30785/90497 relations, 1407246 partial edges, 6641664 polynomials, 400 s elapsed, about 380 s to go
  31667/90497 relations, 1426439 partial edges, 6731776 polynomials, 406 s elapsed, about 362 s to go
  32422/90497 relations, 1442922 partial edges, 6809600 polynomials, 411 s elapsed, about 365 s to go
  33333/90497 relations, 1462099 partial edges, 6899712 polynomials, 416 s elapsed, about 320 s to go
  34172/90497 relations, 1480378 partial edges, 6985728 polynomials, 421 s elapsed, about 330 s to go
  35074/90497 relations, 1498689 partial edges, 7071744 polynomials, 426 s elapsed, about 331 s to go
  36096/90497 relations, 1518746 partial edges, 7165952 polynomials, 432 s elapsed, about 311 s to go
  36939/90497 relations, 1536030 partial edges, 7247872 polynomials, 437 s elapsed, about 345 s to go
  37913/90497 relations, 1555097 partial edges, 7337984 polynomials, 442 s elapsed, about 305 s to go
  38836/90497 relations, 1573049 partial edges, 7424000 polynomials, 447 s elapsed, about 294 s to go
  39778/90497 relations, 1590252 partial edges, 7505920 polynomials, 452 s elapsed, about 287 s to go
  40806/90497 relations, 1609262 partial edges, 7596032 polynomials, 458 s elapsed, about 269 s to go
  41923/90497 relations, 1629321 partial edges, 7690240 polynomials, 465 s elapsed, about 283 s to go
  43313/90497 relations, 1652862 partial edges, 7800832 polynomials, 470 s elapsed, about 236 s to go
  44323/90497 relations, 1670083 partial edges, 7882752 polynomials, 475 s elapsed, about 231 s to go
  45263/90497 relations, 1686184 partial edges, 7960576 polynomials, 480 s elapsed, about 243 s to go
  46365/90497 relations, 1704392 partial edges, 8046592 polynomials, 485 s elapsed, about 220 s to go
  47541/90497 relations, 1723360 partial edges, 8136704 polynomials, 490 s elapsed, about 204 s to go
  48397/90497 relations, 1737282 partial edges, 8202240 polynomials, 496 s elapsed, about 215 s to go
  49863/90497 relations, 1760589 partial edges, 8312832 polynomials, 501 s elapsed, about 141 s to go
  51001/90497 relations, 1777967 partial edges, 8394752 polynomials, 506 s elapsed, about 153 s to go
  52228/90497 relations, 1796437 partial edges, 8480768 polynomials, 512 s elapsed, about 155 s to go
  53543/90497 relations, 1814838 partial edges, 8566784 polynomials, 517 s elapsed, about 149 s to go
  54881/90497 relations, 1834149 partial edges, 8656896 polynomials, 522 s elapsed, about 139 s to go
  56215/90497 relations, 1852548 partial edges, 8744960 polynomials, 528 s elapsed, about 149 s to go
  57573/90497 relations, 1870556 partial edges, 8830976 polynomials, 534 s elapsed, about 141 s to go
  59511/90497 relations, 1896778 partial edges, 8953856 polynomials, 539 s elapsed, about 117 s to go
  61076/90497 relations, 1916811 partial edges, 9048064 polynomials, 544 s elapsed, about 99 s to go
  62557/90497 relations, 1935795 partial edges, 9138176 polynomials, 550 s elapsed, about 97 s to go
  64027/90497 relations, 1954003 partial edges, 9224192 polynomials, 555 s elapsed, about 92 s to go
  65559/90497 relations, 1973309 partial edges, 9314304 polynomials, 560 s elapsed, about 86 s to go
  67076/90497 relations, 1991619 partial edges, 9400320 polynomials, 565 s elapsed, about 78 s to go
  68353/90497 relations, 2007208 partial edges, 9474048 polynomials, 570 s elapsed, about 80 s to go
  70137/90497 relations, 2027630 partial edges, 9570304 polynomials, 575 s elapsed, about 67 s to go
  71503/90497 relations, 2043775 partial edges, 9648128 polynomials, 581 s elapsed, about 68 s to go
  73279/90497 relations, 2063829 partial edges, 9742336 polynomials, 586 s elapsed, about 49 s to go
  75441/90497 relations, 2088006 partial edges, 9857024 polynomials, 592 s elapsed, about 40 s to go
  77177/90497 relations, 2106832 partial edges, 9947136 polynomials, 597 s elapsed, about 36 s to go
  78854/90497 relations, 2124907 partial edges, 10033152 polynomials, 603 s elapsed, about 34 s to go
  81231/90497 relations, 2149170 partial edges, 10147840 polynomials, 608 s elapsed, about 20 s to go
  82574/90497 relations, 2163140 partial edges, 10213376 polynomials, 614 s elapsed, about 22 s to go
  84955/90497 relations, 2187547 partial edges, 10328064 polynomials, 619 s elapsed, about 14 s to go
  86923/90497 relations, 2207359 partial edges, 10422272 polynomials, 624 s elapsed, about 9 s to go
  88712/90497 relations, 2224717 partial edges, 10504192 polynomials, 629 s elapsed, about 5 s to go
4886189752164427378644372029277135326505927538331605193174055598646955377478099751972989945355617869 = 59938415729340155143870756798596946533779694910219 × 81520168538131931077736000604314114233924090649351
2 factors, product checks: true, all probable prime: true, 661.57 s

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
