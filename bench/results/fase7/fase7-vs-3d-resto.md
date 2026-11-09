# Comparativa de benchmarks alcp

Base: **f7i534**. Speedup = mediana(base) / mediana(variante); > 1 significa más rápido.
IQR = rango intercuartílico de las muestras. Las medias geométricas excluyen escenarios con error.

- **f7i534**: modo interleave, 2026-11-09T15:43:28+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -fno-ident -Wl,--no-insert-timestamp`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-93-g0877c55, sha256 1d935c5161da  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: no (batería)
- **f7resto3**: modo interleave, 2026-11-09T15:43:54+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -fno-ident -Wl,--no-insert-timestamp -Wa,-mbig-obj -Wl,--build-id`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-99-gd3c455e (con cambios sin commit), sha256 ce0110a2d72b  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: no (batería)

Checksums: todos los escenarios coinciden con la base.

## Media geométrica del speedup

| ámbito | f7resto3 |
|---|---|
| **global** | 1.025x (n=128) |
| grupo `fast` | 0.994x (n=10) |
| grupo `fp.m127` | 1.006x (n=14) |
| grupo `fp.m61` | 0.993x (n=14) |
| grupo `fp.ntt` | 1.112x (n=16) |
| grupo `fp.p101` | 1.068x (n=14) |
| grupo `fp.p256` | 1.004x (n=16) |
| grupo `fp.p7` | 1.048x (n=14) |
| grupo `fq.build` | 1.090x (n=2) |
| grupo `fq.m61e2` | 0.995x (n=10) |
| grupo `fq.q2e8` | 1.008x (n=6) |
| grupo `fq.q3e30` | 0.970x (n=6) |
| grupo `fq.q7e4` | 0.992x (n=6) |
| operación `build` | 1.090x (n=2) |
| operación `dft` | 0.982x (n=3) |
| operación `divmod` | 1.122x (n=8) |
| operación `elem_inv` | 1.007x (n=10) |
| operación `elem_muladd` | 1.006x (n=10) |
| operación `elem_pow` | 1.006x (n=6) |
| operación `factor` | 0.990x (n=10) |
| operación `gcd` | 1.177x (n=8) |
| operación `gcd_ext` | 1.059x (n=6) |
| operación `idft` | 1.002x (n=1) |
| operación `is_irreducible` | 1.082x (n=7) |
| operación `mul_auto` | 1.000x (n=7) |
| operación `mul_fft` | 1.026x (n=6) |
| operación `mul_karatsuba` | 0.984x (n=16) |
| operación `mul_school` | 1.004x (n=16) |
| operación `pow_mod` | 0.993x (n=6) |
| operación `toeplitz_divmod` | 1.012x (n=2) |
| operación `toeplitz_lower_inverse` | 0.994x (n=3) |
| operación `toeplitz_lower_vec` | 0.981x (n=1) |

## Detalle por escenario

| escenario | base mediana | base IQR | f7resto3 mediana | f7resto3 mín | f7resto3 IQR | speedup | checksum |
|---|---:|---:|---:|---:|---:|---:|---|
| `fp.p7.elem_muladd.n400000` | 4.080 ms | 124.6 us | 4.076 ms | 4.062 ms | 176.9 us | 1.001x | ok `cf4797a827e7f9ec` |
| `fp.p7.elem_inv.n100000` | 815.1 us | 9.7 us | 815.3 us | 808.8 us | 13.4 us | 1.000x | ok `7beb5a3238b919cb` |
| `fp.p7.elem_pow.n100000` | 2.253 ms | 206.0 us | 2.258 ms | 2.245 ms | 193.7 us | 0.998x | ok `a7f8c429d0d81b73` |
| `fp.p7.mul_school.d64` | 52.4 us | 0.4 us | 52.5 us | 51.9 us | 0.6 us | 0.998x | ok `559357b5bc3c57c2` |
| `fp.p7.mul_school.d512` | 127.4 us | 2.7 us | 127.4 us | 126.6 us | 0.9 us | 1.000x | ok `5936b7480b151eb3` |
| `fp.p7.mul_karatsuba.d512` | 195.2 us | 44.2 us | 211.7 us | 178.7 us | 60.9 us | 0.922x | ok `d72c59b658ac0379` |
| `fp.p7.mul_karatsuba.d4096` | 3.171 ms | 476.0 us | 2.941 ms | 2.250 ms | 696.9 us | 1.078x | ok `00f3dde93cbd25a4` |
| `fp.p7.mul_auto.d2048` | 121.9 us | 1.8 us | 121.8 us | 121.2 us | 1.0 us | 1.001x | ok `2b256c747e0070c7` |
| `fp.p7.divmod.d2000_1000` | 212.6 us | 14.3 us | 212.9 us | 206.5 us | 25.2 us | 0.999x | ok `958d383cda71b699` |
| `fp.p7.gcd.d500` | 620.5 us | 13.7 us | 456.5 us | 455.1 us | 8.4 us | 1.359x | ok `cdf18334cf9b13d7` |
| `fp.p7.gcd_ext.d400` | 798.6 us | 147.9 us | 713.5 us | 668.3 us | 57.8 us | 1.119x | ok `984ac156886a3f82` |
| `fp.p7.pow_mod.d128_k48` | 1.972 ms | 207.4 us | 1.990 ms | 1.912 ms | 383.1 us | 0.991x | ok `e3cedb8421cd962f` |
| `fp.p7.is_irreducible.d64` | 915.3 us | 58.1 us | 743.6 us | 721.9 us | 92.6 us | 1.231x | ok `9cb1bec87bc433d0` |
| `fp.p7.factor.d62` | 356.5 us | 28.4 us | 338.0 us | 319.8 us | 29.6 us | 1.055x | ok `3425aecb393956ba` |
| `fp.p101.elem_muladd.n400000` | 4.097 ms | 126.4 us | 4.197 ms | 4.063 ms | 405.4 us | 0.976x | ok `e51aa870ce6269eb` |
| `fp.p101.elem_inv.n100000` | 1.616 ms | 36.0 us | 1.615 ms | 1.606 ms | 47.4 us | 1.001x | ok `69d137d626075364` |
| `fp.p101.elem_pow.n100000` | 4.711 ms | 526.5 us | 4.794 ms | 4.679 ms | 364.5 us | 0.983x | ok `09fec44561c97afa` |
| `fp.p101.mul_school.d64` | 56.9 us | 0.7 us | 56.9 us | 56.4 us | 0.8 us | 1.000x | ok `5a0fa852855ea8b0` |
| `fp.p101.mul_school.d512` | 147.7 us | 1.7 us | 147.9 us | 147.4 us | 2.7 us | 0.999x | ok `8289633c899317a1` |
| `fp.p101.mul_karatsuba.d512` | 172.5 us | 22.6 us | 176.9 us | 162.4 us | 27.3 us | 0.975x | ok `fea453dad54c0188` |
| `fp.p101.mul_karatsuba.d4096` | 3.331 ms | 1.089 ms | 3.178 ms | 2.619 ms | 661.0 us | 1.048x | ok `82f83f0b7fe7713b` |
| `fp.p101.mul_auto.d2048` | 180.7 us | 12.1 us | 180.3 us | 178.7 us | 9.8 us | 1.002x | ok `3631fff5520b84d9` |
| `fp.p101.divmod.d2000_1000` | 314.5 us | 37.1 us | 310.3 us | 300.7 us | 22.2 us | 1.014x | ok `ce32beed1dda91cd` |
| `fp.p101.gcd.d500` | 746.1 us | 96.3 us | 468.6 us | 461.4 us | 18.1 us | 1.592x | ok `5cad0c7c5f168edb` |
| `fp.p101.gcd_ext.d400` | 907.7 us | 94.1 us | 778.4 us | 729.0 us | 219.2 us | 1.166x | ok `4fe0534b04a103a1` |
| `fp.p101.pow_mod.d128_k20` | 2.577 ms | 356.2 us | 2.500 ms | 2.461 ms | 199.1 us | 1.031x | ok `bc5334d07b025488` |
| `fp.p101.is_irreducible.d48` | 675.6 us | 45.7 us | 540.9 us | 526.4 us | 26.1 us | 1.249x | ok `5146dc1467eee05b` |
| `fp.p101.factor.d100` | 1.600 ms | 168.5 us | 1.521 ms | 1.468 ms | 76.4 us | 1.052x | ok `48a7604aa4ae63ee` |
| `fp.ntt.elem_muladd.n400000` | 6.155 ms | 463.2 us | 6.148 ms | 5.960 ms | 415.9 us | 1.001x | ok `682b833fe4edd17b` |
| `fp.ntt.elem_inv.n100000` | 7.622 ms | 691.7 us | 7.623 ms | 7.150 ms | 587.9 us | 1.000x | ok `d275b741e6b3050e` |
| `fp.ntt.elem_pow.n50000` | 17.270 ms | 1.162 ms | 17.125 ms | 16.385 ms | 899.7 us | 1.008x | ok `c27b88ab2935c700` |
| `fp.ntt.mul_school.d64` | 127.9 us | 5.8 us | 126.9 us | 125.2 us | 6.1 us | 1.008x | ok `67c7e1645a9cc6ae` |
| `fp.ntt.mul_school.d512` | 347.6 us | 37.7 us | 344.7 us | 343.6 us | 14.4 us | 1.008x | ok `e296393a6a637046` |
| `fp.ntt.mul_karatsuba.d512` | 256.8 us | 36.7 us | 236.4 us | 227.1 us | 29.6 us | 1.086x | ok `964e15063fb16df4` |
| `fp.ntt.mul_karatsuba.d4096` | 3.881 ms | 1.049 ms | 4.325 ms | 3.472 ms | 1.151 ms | 0.897x | ok `d656b678f0c2ca16` |
| `fp.ntt.mul_fft.d4096` | 513.8 us | 107.3 us | 500.1 us | 394.0 us | 86.3 us | 1.027x | ok `f4c81b85220fb683` |
| `fp.ntt.mul_fft.d32768` | 2.794 ms | 285.9 us | 2.582 ms | 2.250 ms | 361.9 us | 1.082x | ok `8228e76950faded3` |
| `fp.ntt.mul_auto.d2048` | 341.5 us | 106.6 us | 343.8 us | 309.8 us | 74.8 us | 0.993x | ok `02e67c888996018d` |
| `fp.ntt.divmod.d2000_1000` | 1.051 ms | 81.2 us | 441.5 us | 433.6 us | 20.6 us | 2.381x | ok `7adfbb2614325099` |
| `fp.ntt.gcd.d500` | 794.6 us | 46.8 us | 498.0 us | 495.8 us | 22.1 us | 1.596x | ok `678962135e7a2447` |
| `fp.ntt.gcd_ext.d400` | 1.070 ms | 251.9 us | 910.0 us | 838.8 us | 99.8 us | 1.176x | ok `f9dc4e400a7954ee` |
| `fp.ntt.pow_mod.d128_k5` | 7.363 ms | 1.337 ms | 7.540 ms | 6.923 ms | 852.5 us | 0.977x | ok `0eb2514aceb37261` |
| `fp.ntt.is_irreducible.d32` | 506.7 us | 46.8 us | 438.1 us | 425.2 us | 25.8 us | 1.157x | ok `1f2827a0d8aeb65e` |
| `fp.ntt.factor.d40` | 1.630 ms | 94.7 us | 1.655 ms | 1.558 ms | 429.4 us | 0.985x | ok `297c76c423fe55aa` |
| `fp.m61.elem_muladd.n400000` | 5.326 ms | 437.6 us | 5.317 ms | 5.241 ms | 275.6 us | 1.002x | ok `b8e325f73ad90582` |
| `fp.m61.elem_inv.n100000` | 20.691 ms | 839.6 us | 20.478 ms | 19.934 ms | 699.0 us | 1.010x | ok `41e7fe38de6fdcc1` |
| `fp.m61.elem_pow.n20000` | 15.155 ms | 569.6 us | 15.082 ms | 14.589 ms | 535.7 us | 1.005x | ok `a3680ece7fc07932` |
| `fp.m61.mul_school.d64` | 186.2 us | 7.3 us | 186.5 us | 185.4 us | 4.2 us | 0.998x | ok `bed6e9e86907de79` |
| `fp.m61.mul_school.d512` | 487.7 us | 41.9 us | 485.4 us | 481.9 us | 27.8 us | 1.005x | ok `43c3884710f7d468` |
| `fp.m61.mul_karatsuba.d512` | 325.3 us | 48.2 us | 316.3 us | 296.3 us | 25.3 us | 1.028x | ok `a2dcbeaa6574fcb8` |
| `fp.m61.mul_karatsuba.d4096` | 4.281 ms | 955.2 us | 4.589 ms | 3.860 ms | 1.550 ms | 0.933x | ok `37d5deaa1cfeb600` |
| `fp.m61.mul_auto.d2048` | 1.639 ms | 308.2 us | 1.654 ms | 1.300 ms | 437.9 us | 0.991x | ok `73174ffa4e90d6fa` |
| `fp.m61.divmod.d2000_1000` | 2.014 ms | 115.3 us | 2.013 ms | 1.984 ms | 122.3 us | 1.000x | ok `a8710fc5aaa24d09` |
| `fp.m61.gcd.d500` | 852.5 us | 27.8 us | 825.6 us | 811.6 us | 113.5 us | 1.033x | ok `fda394dde025c1a2` |
| `fp.m61.gcd_ext.d400` | 1.596 ms | 122.6 us | 1.689 ms | 1.625 ms | 192.4 us | 0.945x | ok `be3d81e8acc44b48` |
| `fp.m61.pow_mod.d64_k2` | 3.637 ms | 435.4 us | 3.777 ms | 3.539 ms | 838.9 us | 0.963x | ok `394fb320787700fe` |
| `fp.m61.is_irreducible.d24` | 646.0 us | 27.6 us | 654.4 us | 638.2 us | 41.8 us | 0.987x | ok `9047cf1d3fc29abf` |
| `fp.m61.factor.d26` | 2.130 ms | 218.3 us | 2.100 ms | 2.017 ms | 242.7 us | 1.014x | ok `931b8d0b8b91b0b4` |
| `fp.m127.elem_muladd.n200000` | 7.616 ms | 598.3 us | 7.409 ms | 7.244 ms | 591.2 us | 1.028x | ok `16fab7b139aa9aab` |
| `fp.m127.elem_inv.n50000` | 73.900 ms | 3.511 ms | 76.888 ms | 70.533 ms | 4.723 ms | 0.961x | ok `c7d67d6018036bf1` |
| `fp.m127.elem_pow.n10000` | 47.062 ms | 1.708 ms | 47.048 ms | 44.261 ms | 1.662 ms | 1.000x | ok `308896e55b40c815` |
| `fp.m127.mul_school.d64` | 675.7 us | 31.2 us | 675.5 us | 667.3 us | 26.4 us | 1.000x | ok `a9664408b0fb88e7` |
| `fp.m127.mul_school.d512` | 1.886 ms | 111.7 us | 1.884 ms | 1.861 ms | 112.1 us | 1.001x | ok `98eca357802203c0` |
| `fp.m127.mul_karatsuba.d512` | 1.003 ms | 81.6 us | 991.3 us | 942.6 us | 74.8 us | 1.012x | ok `f773c153fde34172` |
| `fp.m127.mul_karatsuba.d2048` | 4.643 ms | 591.1 us | 4.297 ms | 3.490 ms | 572.1 us | 1.081x | ok `695e0d0f13c14e87` |
| `fp.m127.mul_auto.d2048` | 3.594 ms | 773.6 us | 3.431 ms | 3.199 ms | 503.2 us | 1.047x | ok `1bc667b834b9e559` |
| `fp.m127.divmod.d2000_1000` | 7.518 ms | 707.1 us | 7.268 ms | 7.096 ms | 660.5 us | 1.034x | ok `a8babdcd119eb404` |
| `fp.m127.gcd.d400` | 4.119 ms | 432.3 us | 4.034 ms | 3.944 ms | 474.7 us | 1.021x | ok `1b15fac38c54accc` |
| `fp.m127.gcd_ext.d300` | 6.713 ms | 555.0 us | 7.052 ms | 6.351 ms | 881.7 us | 0.952x | ok `45a1b9001453c6ef` |
| `fp.m127.pow_mod.d64_k1` | 17.732 ms | 830.4 us | 17.537 ms | 16.765 ms | 883.8 us | 1.011x | ok `9e9e535b3c4e4e7d` |
| `fp.m127.is_irreducible.d16` | 2.057 ms | 214.5 us | 2.035 ms | 1.961 ms | 321.3 us | 1.011x | ok `2cfb69ee2746aa9d` |
| `fp.m127.factor.d19` | 12.924 ms | 1.819 ms | 13.891 ms | 12.305 ms | 1.146 ms | 0.930x | ok `7b0defab2417737c` |
| `fp.p256.elem_muladd.n200000` | 83.995 ms | 6.625 ms | 81.650 ms | 78.821 ms | 5.938 ms | 1.029x | ok `cc3cab838e9b19f1` |
| `fp.p256.elem_inv.n50000` | 171.982 ms | 7.157 ms | 167.874 ms | 162.853 ms | 8.891 ms | 1.024x | ok `89f920a8eba0c743` |
| `fp.p256.elem_pow.n2000` | 49.963 ms | 4.192 ms | 48.002 ms | 46.730 ms | 1.532 ms | 1.041x | ok `40472b256a3c877d` |
| `fp.p256.mul_school.d64` | 6.236 ms | 1.216 ms | 6.029 ms | 5.759 ms | 807.7 us | 1.034x | ok `f5d3d87e3ed0977e` |
| `fp.p256.mul_school.d512` | 16.117 ms | 1.424 ms | 15.649 ms | 14.501 ms | 1.409 ms | 1.030x | ok `6fd19bd7b3ceab69` |
| `fp.p256.mul_karatsuba.d512` | 5.510 ms | 1.314 ms | 5.829 ms | 4.910 ms | 442.4 us | 0.945x | ok `4b2f43bf0f878408` |
| `fp.p256.mul_karatsuba.d2048` | 25.760 ms | 2.706 ms | 27.881 ms | 24.008 ms | 4.830 ms | 0.924x | ok `b4b89a80c2a3aaa9` |
| `fp.p256.mul_fft.d4096` | 42.051 ms | 8.065 ms | 42.065 ms | 35.273 ms | 5.790 ms | 1.000x | ok `60ed510564513dfb` |
| `fp.p256.mul_fft.d8192` | 88.996 ms | 7.895 ms | 83.038 ms | 78.627 ms | 5.936 ms | 1.072x | ok `47812935df4de5fe` |
| `fp.p256.mul_auto.d2048` | 8.784 ms | 881.0 us | 8.992 ms | 8.279 ms | 867.4 us | 0.977x | ok `3d89402081500506` |
| `fp.p256.divmod.d1000_500` | 9.209 ms | 864.4 us | 9.144 ms | 7.770 ms | 745.2 us | 1.007x | ok `0bebce51d98eb3f4` |
| `fp.p256.gcd.d300` | 27.273 ms | 1.183 ms | 27.397 ms | 25.950 ms | 943.8 us | 0.995x | ok `f69e4050cb1a5861` |
| `fp.p256.gcd_ext.d200` | 31.174 ms | 1.238 ms | 30.538 ms | 29.367 ms | 1.282 ms | 1.021x | ok `1a9df499c28b9f49` |
| `fp.p256.pow_mod.d64_k1` | 166.922 ms | 11.425 ms | 169.143 ms | 162.073 ms | 9.653 ms | 0.987x | ok `dbe2b595567fbb9c` |
| `fp.p256.is_irreducible.d12` | 17.256 ms | 1.406 ms | 17.522 ms | 16.884 ms | 1.130 ms | 0.985x | ok `7b839a40c98b0ef0` |
| `fp.p256.factor.d12` | 37.008 ms | 1.522 ms | 36.854 ms | 34.714 ms | 1.685 ms | 1.004x | ok `c33e3512ed388a7e` |
| `fq.q3e30.elem_muladd.n4000` | 5.188 ms | 619.6 us | 5.021 ms | 4.892 ms | 746.9 us | 1.033x | ok `16edfcc06a08028a` |
| `fq.q3e30.elem_inv.n100` | 360.2 us | 31.2 us | 358.5 us | 348.4 us | 26.7 us | 1.005x | ok `a8167c01d1cfc5dd` |
| `fq.q3e30.mul_school.d32` | 368.8 us | 30.4 us | 366.8 us | 360.9 us | 18.7 us | 1.005x | ok `ed0ca8eb7d325de3` |
| `fq.q3e30.mul_karatsuba.d128` | 2.310 ms | 301.1 us | 2.812 ms | 2.193 ms | 444.0 us | 0.822x | ok `2464b901bd242c7d` |
| `fq.q3e30.divmod.d64_32` | 399.2 us | 13.9 us | 403.5 us | 394.8 us | 17.8 us | 0.989x | ok `a2b35f414a139741` |
| `fq.q3e30.factor.d5` | 2.802 ms | 173.9 us | 2.856 ms | 2.751 ms | 229.8 us | 0.981x | ok `417a5b753c7a2a07` |
| `fq.q2e8.elem_muladd.n50000` | 347.9 us | 8.1 us | 351.9 us | 346.6 us | 7.4 us | 0.989x | ok `c0ffe16ee54f128c` |
| `fq.q2e8.elem_inv.n5000` | 13.7 us | 0.3 us | 12.8 us | 12.5 us | 0.6 us | 1.070x | ok `81d71b5687e3e735` |
| `fq.q2e8.mul_school.d128` | 30.1 us | 0.8 us | 30.0 us | 29.3 us | 0.6 us | 1.003x | ok `c056d9d992271db9` |
| `fq.q2e8.mul_karatsuba.d512` | 225.1 us | 11.3 us | 225.6 us | 215.8 us | 19.8 us | 0.998x | ok `de1089fdce26d02c` |
| `fq.q2e8.is_irreducible.d16` | 92.6 us | 4.0 us | 93.3 us | 90.9 us | 3.2 us | 0.992x | ok `9478d590a29a9cbf` |
| `fq.q2e8.factor.d28` | 309.4 us | 21.4 us | 310.4 us | 300.6 us | 22.5 us | 0.997x | ok `d2fc95c10ececedd` |
| `fq.q7e4.elem_muladd.n50000` | 755.4 us | 25.1 us | 759.5 us | 746.9 us | 22.8 us | 0.995x | ok `a6cf18a43d40c5c0` |
| `fq.q7e4.elem_inv.n5000` | 52.1 us | 0.7 us | 52.1 us | 51.6 us | 0.6 us | 1.000x | ok `8664ac6d09409f94` |
| `fq.q7e4.mul_school.d128` | 187.3 us | 22.5 us | 189.5 us | 158.4 us | 22.0 us | 0.989x | ok `71d186560b1d8153` |
| `fq.q7e4.mul_karatsuba.d512` | 1.287 ms | 72.4 us | 1.292 ms | 1.253 ms | 52.5 us | 0.996x | ok `6a2e9b751ab48e99` |
| `fq.q7e4.gcd.d200` | 564.6 us | 33.8 us | 565.8 us | 550.5 us | 33.0 us | 0.998x | ok `8cf19ac90e512feb` |
| `fq.q7e4.factor.d20` | 330.3 us | 42.5 us | 338.2 us | 303.7 us | 41.9 us | 0.977x | ok `ce23318b1b5fd48e` |
| `fq.m61e2.elem_muladd.n100000` | 23.379 ms | 965.2 us | 23.205 ms | 22.719 ms | 991.7 us | 1.008x | ok `1814347c960518d8` |
| `fq.m61e2.elem_inv.n1000` | 962.3 us | 25.1 us | 956.6 us | 931.4 us | 25.4 us | 1.006x | ok `45326d8095032432` |
| `fq.m61e2.mul_school.d128` | 286.1 us | 10.8 us | 288.2 us | 281.8 us | 15.2 us | 0.993x | ok `2cb38bd639886e3c` |
| `fq.m61e2.mul_karatsuba.d512` | 3.706 ms | 699.0 us | 3.568 ms | 3.093 ms | 719.0 us | 1.039x | ok `f138ad86c0a7a838` |
| `fq.m61e2.mul_fft.d1024` | 35.648 ms | 1.586 ms | 35.471 ms | 34.032 ms | 1.256 ms | 1.005x | ok `a45642d67f5a9912` |
| `fq.m61e2.mul_fft.d4096` | 89.030 ms | 1.640 ms | 91.378 ms | 87.980 ms | 4.562 ms | 0.974x | ok `51844e8f95a9c0ab` |
| `fq.m61e2.mul_auto.d2048` | 6.330 ms | 322.1 us | 6.384 ms | 6.063 ms | 479.6 us | 0.992x | ok `76e66ad22fb3fa2e` |
| `fq.m61e2.divmod.d500_250` | 1.093 ms | 31.4 us | 1.083 ms | 1.059 ms | 39.3 us | 1.009x | ok `53509a111940ee2a` |
| `fq.m61e2.gcd.d150` | 3.353 ms | 307.0 us | 3.295 ms | 3.221 ms | 206.1 us | 1.017x | ok `bd8c3c98cc56df84` |
| `fq.m61e2.factor.d7` | 6.146 ms | 506.9 us | 6.752 ms | 6.157 ms | 590.1 us | 0.910x | ok `d61430d11ee187c3` |
| `fq.build.q3e30` | 266.7 us | 16.8 us | 246.0 us | 234.4 us | 19.8 us | 1.084x | ok `4973a3e567e2fc74` |
| `fq.build.q7e20` | 258.2 us | 17.5 us | 235.7 us | 224.8 us | 76.7 us | 1.095x | ok `9c377cd554aa108e` |
| `fast.ntt.dft.k10` | 2.925 ms | 460.5 us | 3.183 ms | 2.855 ms | 999.5 us | 0.919x | ok `b61f454e63ae6d8d` |
| `fast.ntt.dft.k11` | 2.975 ms | 491.7 us | 2.941 ms | 2.880 ms | 340.8 us | 1.012x | ok `311d4fed0f42034c` |
| `fast.ntt.dft.k12` | 3.076 ms | 593.5 us | 3.017 ms | 2.887 ms | 787.5 us | 1.020x | ok `70691e178cfc0c9d` |
| `fast.ntt.idft.k12` | 3.023 ms | 412.5 us | 3.017 ms | 2.924 ms | 509.4 us | 1.002x | ok `e0677b59c4e33615` |
| `fast.ntt.toeplitz_lower_inverse.n1024` | 103.5 us | 11.6 us | 102.7 us | 97.6 us | 10.8 us | 1.008x | ok `7749eb8ff8cac3d0` |
| `fast.ntt.toeplitz_lower_inverse.n8192` | 2.133 ms | 463.1 us | 2.249 ms | 1.962 ms | 473.0 us | 0.949x | ok `0d31554b19d8b996` |
| `fast.ntt.toeplitz_lower_vec.n4096` | 362.1 us | 77.0 us | 369.1 us | 329.7 us | 69.8 us | 0.981x | ok `86e36ad2a1eacb05` |
| `fast.ntt.toeplitz_divmod.d2000_1000` | 6.840 ms | 448.3 us | 6.354 ms | 6.157 ms | 292.7 us | 1.077x | ok `b138cdcceebcc32e` |
| `fast.p7.toeplitz_lower_inverse.n1024` | 142.8 us | 27.5 us | 138.8 us | 131.3 us | 20.4 us | 1.029x | ok `5df08c92a814155a` |
| `fast.p7.toeplitz_divmod.d2000_1000` | 2.838 ms | 357.8 us | 2.984 ms | 2.945 ms | 302.7 us | 0.951x | ok `433d13716b2f95e8` |
