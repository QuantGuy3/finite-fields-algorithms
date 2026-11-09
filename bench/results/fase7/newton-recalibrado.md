# Comparativa de benchmarks alcp

Base: **f7resto3**. Speedup = mediana(base) / mediana(variante); > 1 significa más rápido.
IQR = rango intercuartílico de las muestras. Las medias geométricas excluyen escenarios con error.

- **f7resto3**: modo interleave, 2026-11-09T16:03:07+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -fno-ident -Wl,--no-insert-timestamp -Wa,-mbig-obj -Wl,--build-id`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-99-gd3c455e (con cambios sin commit), sha256 ce0110a2d72b  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: no (batería)
- **f7newton2**: modo interleave, 2026-11-09T16:03:30+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -fno-ident -Wl,--no-insert-timestamp -Wl,--build-id=sha1`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-100-g5cd0529 (con cambios sin commit), sha256 d7c52180bb2c  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: no (batería)

Checksums: todos los escenarios coinciden con la base.

## Media geométrica del speedup

| ámbito | f7newton2 |
|---|---|
| **global** | 1.011x (n=128) |
| grupo `fast` | 1.023x (n=10) |
| grupo `fp.m127` | 1.037x (n=14) |
| grupo `fp.m61` | 1.002x (n=14) |
| grupo `fp.ntt` | 1.061x (n=16) |
| grupo `fp.p101` | 1.004x (n=14) |
| grupo `fp.p256` | 0.971x (n=16) |
| grupo `fp.p7` | 1.006x (n=14) |
| grupo `fq.build` | 1.023x (n=2) |
| grupo `fq.m61e2` | 1.008x (n=10) |
| grupo `fq.q2e8` | 0.997x (n=6) |
| grupo `fq.q3e30` | 0.987x (n=6) |
| grupo `fq.q7e4` | 0.998x (n=6) |
| operación `build` | 1.023x (n=2) |
| operación `dft` | 1.015x (n=3) |
| operación `divmod` | 0.989x (n=8) |
| operación `elem_inv` | 0.991x (n=10) |
| operación `elem_muladd` | 1.005x (n=10) |
| operación `elem_pow` | 1.005x (n=6) |
| operación `factor` | 1.011x (n=10) |
| operación `gcd` | 1.051x (n=8) |
| operación `gcd_ext` | 0.978x (n=6) |
| operación `idft` | 1.040x (n=1) |
| operación `is_irreducible` | 1.016x (n=7) |
| operación `mul_auto` | 1.001x (n=7) |
| operación `mul_fft` | 0.995x (n=6) |
| operación `mul_karatsuba` | 0.993x (n=16) |
| operación `mul_school` | 1.002x (n=16) |
| operación `pow_mod` | 1.150x (n=6) |
| operación `toeplitz_divmod` | 1.028x (n=2) |
| operación `toeplitz_lower_inverse` | 1.012x (n=3) |
| operación `toeplitz_lower_vec` | 1.051x (n=1) |

## Detalle por escenario

| escenario | base mediana | base IQR | f7newton2 mediana | f7newton2 mín | f7newton2 IQR | speedup | checksum |
|---|---:|---:|---:|---:|---:|---:|---|
| `fp.p7.elem_muladd.n400000` | 4.087 ms | 203.5 us | 4.142 ms | 4.065 ms | 178.7 us | 0.987x | ok `cf4797a827e7f9ec` |
| `fp.p7.elem_inv.n100000` | 815.2 us | 21.2 us | 825.2 us | 808.8 us | 121.0 us | 0.988x | ok `7beb5a3238b919cb` |
| `fp.p7.elem_pow.n100000` | 2.280 ms | 209.4 us | 2.294 ms | 2.245 ms | 348.6 us | 0.994x | ok `a7f8c429d0d81b73` |
| `fp.p7.mul_school.d64` | 52.5 us | 1.0 us | 52.8 us | 51.9 us | 1.6 us | 0.994x | ok `559357b5bc3c57c2` |
| `fp.p7.mul_school.d512` | 127.5 us | 6.8 us | 127.7 us | 127.0 us | 24.3 us | 0.998x | ok `5936b7480b151eb3` |
| `fp.p7.mul_karatsuba.d512` | 214.3 us | 69.2 us | 236.2 us | 180.0 us | 119.5 us | 0.907x | ok `d72c59b658ac0379` |
| `fp.p7.mul_karatsuba.d4096` | 2.889 ms | 949.1 us | 2.705 ms | 2.240 ms | 513.0 us | 1.068x | ok `00f3dde93cbd25a4` |
| `fp.p7.mul_auto.d2048` | 121.8 us | 6.8 us | 122.6 us | 121.2 us | 9.7 us | 0.993x | ok `2b256c747e0070c7` |
| `fp.p7.divmod.d2000_1000` | 215.1 us | 41.8 us | 221.3 us | 208.2 us | 35.6 us | 0.972x | ok `958d383cda71b699` |
| `fp.p7.gcd.d500` | 461.4 us | 36.7 us | 411.9 us | 407.1 us | 31.9 us | 1.120x | ok `cdf18334cf9b13d7` |
| `fp.p7.gcd_ext.d400` | 721.8 us | 107.6 us | 728.7 us | 658.1 us | 172.0 us | 0.991x | ok `984ac156886a3f82` |
| `fp.p7.pow_mod.d128_k48` | 2.004 ms | 497.0 us | 2.002 ms | 1.908 ms | 250.5 us | 1.001x | ok `e3cedb8421cd962f` |
| `fp.p7.is_irreducible.d64` | 747.7 us | 103.6 us | 735.7 us | 699.1 us | 109.2 us | 1.016x | ok `9cb1bec87bc433d0` |
| `fp.p7.factor.d62` | 367.4 us | 186.3 us | 345.5 us | 319.8 us | 50.2 us | 1.063x | ok `3425aecb393956ba` |
| `fp.p101.elem_muladd.n400000` | 4.078 ms | 42.0 us | 4.183 ms | 4.062 ms | 325.0 us | 0.975x | ok `e51aa870ce6269eb` |
| `fp.p101.elem_inv.n100000` | 1.615 ms | 41.0 us | 1.661 ms | 1.606 ms | 160.9 us | 0.972x | ok `69d137d626075364` |
| `fp.p101.elem_pow.n100000` | 4.836 ms | 249.5 us | 4.973 ms | 4.681 ms | 1.015 ms | 0.972x | ok `09fec44561c97afa` |
| `fp.p101.mul_school.d64` | 57.1 us | 2.5 us | 56.9 us | 56.4 us | 1.3 us | 1.004x | ok `5a0fa852855ea8b0` |
| `fp.p101.mul_school.d512` | 147.9 us | 10.5 us | 148.0 us | 147.2 us | 8.9 us | 0.999x | ok `8289633c899317a1` |
| `fp.p101.mul_karatsuba.d512` | 170.8 us | 36.3 us | 168.5 us | 160.2 us | 25.5 us | 1.013x | ok `fea453dad54c0188` |
| `fp.p101.mul_karatsuba.d4096` | 2.923 ms | 582.6 us | 3.207 ms | 2.528 ms | 727.0 us | 0.911x | ok `82f83f0b7fe7713b` |
| `fp.p101.mul_auto.d2048` | 180.7 us | 14.8 us | 181.8 us | 178.9 us | 21.7 us | 0.994x | ok `3631fff5520b84d9` |
| `fp.p101.divmod.d2000_1000` | 322.8 us | 52.2 us | 326.6 us | 302.7 us | 76.3 us | 0.988x | ok `ce32beed1dda91cd` |
| `fp.p101.gcd.d500` | 466.1 us | 38.5 us | 398.8 us | 396.1 us | 16.2 us | 1.169x | ok `5cad0c7c5f168edb` |
| `fp.p101.gcd_ext.d400` | 769.1 us | 80.4 us | 757.4 us | 704.4 us | 110.6 us | 1.016x | ok `4fe0534b04a103a1` |
| `fp.p101.pow_mod.d128_k20` | 2.626 ms | 578.6 us | 2.568 ms | 2.479 ms | 214.0 us | 1.022x | ok `bc5334d07b025488` |
| `fp.p101.is_irreducible.d48` | 556.2 us | 183.7 us | 535.4 us | 512.5 us | 53.9 us | 1.039x | ok `5146dc1467eee05b` |
| `fp.p101.factor.d100` | 1.558 ms | 328.4 us | 1.548 ms | 1.450 ms | 289.1 us | 1.006x | ok `48a7604aa4ae63ee` |
| `fp.ntt.elem_muladd.n400000` | 6.296 ms | 569.1 us | 6.134 ms | 5.957 ms | 520.7 us | 1.026x | ok `682b833fe4edd17b` |
| `fp.ntt.elem_inv.n100000` | 7.705 ms | 530.4 us | 7.661 ms | 7.213 ms | 1.001 ms | 1.006x | ok `d275b741e6b3050e` |
| `fp.ntt.elem_pow.n50000` | 17.178 ms | 643.8 us | 17.102 ms | 16.568 ms | 595.1 us | 1.004x | ok `c27b88ab2935c700` |
| `fp.ntt.mul_school.d64` | 131.0 us | 31.8 us | 130.0 us | 125.6 us | 15.4 us | 1.007x | ok `67c7e1645a9cc6ae` |
| `fp.ntt.mul_school.d512` | 346.5 us | 19.9 us | 347.2 us | 343.7 us | 30.6 us | 0.998x | ok `e296393a6a637046` |
| `fp.ntt.mul_karatsuba.d512` | 253.4 us | 40.4 us | 249.1 us | 227.5 us | 46.0 us | 1.017x | ok `964e15063fb16df4` |
| `fp.ntt.mul_karatsuba.d4096` | 3.781 ms | 721.2 us | 3.936 ms | 2.687 ms | 1.122 ms | 0.961x | ok `d656b678f0c2ca16` |
| `fp.ntt.mul_fft.d4096` | 523.8 us | 92.8 us | 524.0 us | 386.8 us | 93.2 us | 1.000x | ok `f4c81b85220fb683` |
| `fp.ntt.mul_fft.d32768` | 2.583 ms | 347.8 us | 2.589 ms | 2.205 ms | 317.8 us | 0.998x | ok `8228e76950faded3` |
| `fp.ntt.mul_auto.d2048` | 402.3 us | 139.7 us | 406.4 us | 316.5 us | 98.1 us | 0.990x | ok `02e67c888996018d` |
| `fp.ntt.divmod.d2000_1000` | 451.3 us | 50.0 us | 449.1 us | 434.7 us | 45.0 us | 1.005x | ok `7adfbb2614325099` |
| `fp.ntt.gcd.d500` | 498.7 us | 34.4 us | 435.0 us | 428.9 us | 28.8 us | 1.147x | ok `678962135e7a2447` |
| `fp.ntt.gcd_ext.d400` | 878.8 us | 87.2 us | 910.4 us | 805.1 us | 166.4 us | 0.965x | ok `f9dc4e400a7954ee` |
| `fp.ntt.pow_mod.d128_k5` | 7.549 ms | 932.0 us | 3.367 ms | 3.237 ms | 395.8 us | 2.242x | ok `0eb2514aceb37261` |
| `fp.ntt.is_irreducible.d32` | 448.6 us | 51.6 us | 435.2 us | 420.5 us | 59.4 us | 1.031x | ok `1f2827a0d8aeb65e` |
| `fp.ntt.factor.d40` | 1.628 ms | 228.0 us | 1.642 ms | 1.529 ms | 350.4 us | 0.991x | ok `297c76c423fe55aa` |
| `fp.m61.elem_muladd.n400000` | 5.375 ms | 302.3 us | 5.362 ms | 5.240 ms | 215.1 us | 1.002x | ok `b8e325f73ad90582` |
| `fp.m61.elem_inv.n100000` | 20.680 ms | 245.7 us | 20.754 ms | 20.116 ms | 1.394 ms | 0.996x | ok `41e7fe38de6fdcc1` |
| `fp.m61.elem_pow.n20000` | 15.122 ms | 222.6 us | 14.796 ms | 14.578 ms | 691.4 us | 1.022x | ok `a3680ece7fc07932` |
| `fp.m61.mul_school.d64` | 188.0 us | 18.1 us | 188.6 us | 185.8 us | 9.5 us | 0.997x | ok `bed6e9e86907de79` |
| `fp.m61.mul_school.d512` | 486.9 us | 26.8 us | 486.3 us | 481.8 us | 32.5 us | 1.001x | ok `43c3884710f7d468` |
| `fp.m61.mul_karatsuba.d512` | 309.9 us | 56.2 us | 309.6 us | 295.9 us | 42.2 us | 1.001x | ok `a2dcbeaa6574fcb8` |
| `fp.m61.mul_karatsuba.d4096` | 4.264 ms | 935.7 us | 4.185 ms | 3.188 ms | 523.7 us | 1.019x | ok `37d5deaa1cfeb600` |
| `fp.m61.mul_auto.d2048` | 1.744 ms | 297.1 us | 1.733 ms | 1.534 ms | 252.3 us | 1.006x | ok `73174ffa4e90d6fa` |
| `fp.m61.divmod.d2000_1000` | 2.036 ms | 238.7 us | 2.099 ms | 1.991 ms | 316.5 us | 0.970x | ok `a8710fc5aaa24d09` |
| `fp.m61.gcd.d500` | 826.6 us | 46.0 us | 831.2 us | 811.9 us | 71.4 us | 0.995x | ok `fda394dde025c1a2` |
| `fp.m61.gcd_ext.d400` | 1.692 ms | 218.9 us | 1.737 ms | 1.622 ms | 291.3 us | 0.974x | ok `be3d81e8acc44b48` |
| `fp.m61.pow_mod.d64_k2` | 3.866 ms | 534.0 us | 3.800 ms | 3.510 ms | 599.8 us | 1.017x | ok `394fb320787700fe` |
| `fp.m61.is_irreducible.d24` | 661.1 us | 69.8 us | 660.1 us | 638.1 us | 93.1 us | 1.001x | ok `9047cf1d3fc29abf` |
| `fp.m61.factor.d26` | 2.160 ms | 541.1 us | 2.108 ms | 2.048 ms | 330.4 us | 1.025x | ok `931b8d0b8b91b0b4` |
| `fp.m127.elem_muladd.n200000` | 8.017 ms | 1.638 ms | 7.522 ms | 7.230 ms | 448.3 us | 1.066x | ok `16fab7b139aa9aab` |
| `fp.m127.elem_inv.n50000` | 80.600 ms | 8.909 ms | 76.807 ms | 73.699 ms | 3.130 ms | 1.049x | ok `c7d67d6018036bf1` |
| `fp.m127.elem_pow.n10000` | 47.514 ms | 2.917 ms | 46.007 ms | 45.117 ms | 911.1 us | 1.033x | ok `308896e55b40c815` |
| `fp.m127.mul_school.d64` | 686.9 us | 240.3 us | 688.2 us | 667.0 us | 40.1 us | 0.998x | ok `a9664408b0fb88e7` |
| `fp.m127.mul_school.d512` | 2.147 ms | 901.9 us | 1.920 ms | 1.864 ms | 277.9 us | 1.118x | ok `98eca357802203c0` |
| `fp.m127.mul_karatsuba.d512` | 1.167 ms | 455.0 us | 1.001 ms | 938.9 us | 162.9 us | 1.165x | ok `f773c153fde34172` |
| `fp.m127.mul_karatsuba.d2048` | 4.285 ms | 963.6 us | 4.054 ms | 3.162 ms | 503.3 us | 1.057x | ok `695e0d0f13c14e87` |
| `fp.m127.mul_auto.d2048` | 3.694 ms | 822.5 us | 3.509 ms | 3.243 ms | 382.3 us | 1.053x | ok `1bc667b834b9e559` |
| `fp.m127.divmod.d2000_1000` | 7.728 ms | 822.7 us | 7.537 ms | 7.170 ms | 704.5 us | 1.025x | ok `a8babdcd119eb404` |
| `fp.m127.gcd.d400` | 4.106 ms | 499.8 us | 4.123 ms | 3.930 ms | 474.3 us | 0.996x | ok `1b15fac38c54accc` |
| `fp.m127.gcd_ext.d300` | 6.725 ms | 824.5 us | 7.146 ms | 6.388 ms | 1.051 ms | 0.941x | ok `45a1b9001453c6ef` |
| `fp.m127.pow_mod.d64_k1` | 17.791 ms | 1.134 ms | 17.639 ms | 16.928 ms | 937.4 us | 1.009x | ok `9e9e535b3c4e4e7d` |
| `fp.m127.is_irreducible.d16` | 2.084 ms | 447.2 us | 2.058 ms | 1.973 ms | 319.6 us | 1.013x | ok `2cfb69ee2746aa9d` |
| `fp.m127.factor.d19` | 13.376 ms | 689.1 us | 13.129 ms | 11.987 ms | 1.592 ms | 1.019x | ok `7b0defab2417737c` |
| `fp.p256.elem_muladd.n200000` | 82.147 ms | 9.075 ms | 82.400 ms | 78.557 ms | 7.022 ms | 0.997x | ok `cc3cab838e9b19f1` |
| `fp.p256.elem_inv.n50000` | 168.510 ms | 12.563 ms | 170.668 ms | 161.573 ms | 8.621 ms | 0.987x | ok `89f920a8eba0c743` |
| `fp.p256.elem_pow.n2000` | 48.377 ms | 1.889 ms | 48.079 ms | 46.148 ms | 2.730 ms | 1.006x | ok `40472b256a3c877d` |
| `fp.p256.mul_school.d64` | 6.399 ms | 1.095 ms | 6.401 ms | 5.736 ms | 1.211 ms | 1.000x | ok `f5d3d87e3ed0977e` |
| `fp.p256.mul_school.d512` | 15.883 ms | 698.7 us | 16.734 ms | 15.385 ms | 1.990 ms | 0.949x | ok `6fd19bd7b3ceab69` |
| `fp.p256.mul_karatsuba.d512` | 4.593 ms | 721.1 us | 5.463 ms | 4.703 ms | 1.098 ms | 0.841x | ok `4b2f43bf0f878408` |
| `fp.p256.mul_karatsuba.d2048` | 25.520 ms | 1.154 ms | 26.020 ms | 24.633 ms | 1.894 ms | 0.981x | ok `b4b89a80c2a3aaa9` |
| `fp.p256.mul_fft.d4096` | 42.874 ms | 2.387 ms | 43.821 ms | 36.728 ms | 4.416 ms | 0.978x | ok `60ed510564513dfb` |
| `fp.p256.mul_fft.d8192` | 88.457 ms | 8.528 ms | 89.654 ms | 82.334 ms | 8.445 ms | 0.987x | ok `47812935df4de5fe` |
| `fp.p256.mul_auto.d2048` | 8.882 ms | 909.0 us | 9.248 ms | 8.153 ms | 859.1 us | 0.960x | ok `3d89402081500506` |
| `fp.p256.divmod.d1000_500` | 9.004 ms | 1.542 ms | 9.487 ms | 7.950 ms | 1.175 ms | 0.949x | ok `0bebce51d98eb3f4` |
| `fp.p256.gcd.d300` | 27.483 ms | 2.656 ms | 27.893 ms | 26.259 ms | 3.394 ms | 0.985x | ok `f69e4050cb1a5861` |
| `fp.p256.gcd_ext.d200` | 30.852 ms | 2.203 ms | 31.302 ms | 30.124 ms | 3.667 ms | 0.986x | ok `1a9df499c28b9f49` |
| `fp.p256.pow_mod.d64_k1` | 170.846 ms | 6.253 ms | 174.045 ms | 162.753 ms | 10.512 ms | 0.982x | ok `dbe2b595567fbb9c` |
| `fp.p256.is_irreducible.d12` | 17.780 ms | 298.0 us | 18.251 ms | 17.082 ms | 2.279 ms | 0.974x | ok `7b839a40c98b0ef0` |
| `fp.p256.factor.d12` | 36.780 ms | 1.546 ms | 37.117 ms | 35.461 ms | 2.226 ms | 0.991x | ok `c33e3512ed388a7e` |
| `fq.q3e30.elem_muladd.n4000` | 5.268 ms | 749.2 us | 5.194 ms | 4.909 ms | 491.5 us | 1.014x | ok `16edfcc06a08028a` |
| `fq.q3e30.elem_inv.n100` | 359.8 us | 33.8 us | 377.9 us | 349.3 us | 96.4 us | 0.952x | ok `a8167c01d1cfc5dd` |
| `fq.q3e30.mul_school.d32` | 367.0 us | 25.0 us | 377.4 us | 360.4 us | 99.3 us | 0.972x | ok `ed0ca8eb7d325de3` |
| `fq.q3e30.mul_karatsuba.d128` | 2.417 ms | 251.5 us | 2.345 ms | 1.990 ms | 375.3 us | 1.031x | ok `2464b901bd242c7d` |
| `fq.q3e30.divmod.d64_32` | 404.1 us | 25.3 us | 413.4 us | 395.8 us | 65.6 us | 0.978x | ok `a2b35f414a139741` |
| `fq.q3e30.factor.d5` | 2.957 ms | 492.1 us | 3.020 ms | 2.774 ms | 336.7 us | 0.979x | ok `417a5b753c7a2a07` |
| `fq.q2e8.elem_muladd.n50000` | 350.6 us | 10.4 us | 351.3 us | 346.6 us | 9.5 us | 0.998x | ok `c0ffe16ee54f128c` |
| `fq.q2e8.elem_inv.n5000` | 12.8 us | 1.1 us | 12.8 us | 12.5 us | 1.2 us | 1.000x | ok `81d71b5687e3e735` |
| `fq.q2e8.mul_school.d128` | 30.1 us | 1.2 us | 30.4 us | 29.4 us | 2.2 us | 0.990x | ok `c056d9d992271db9` |
| `fq.q2e8.mul_karatsuba.d512` | 224.3 us | 25.2 us | 231.9 us | 219.3 us | 37.5 us | 0.967x | ok `de1089fdce26d02c` |
| `fq.q2e8.is_irreducible.d16` | 98.7 us | 18.6 us | 95.2 us | 91.8 us | 19.2 us | 1.037x | ok `9478d590a29a9cbf` |
| `fq.q2e8.factor.d28` | 314.2 us | 29.4 us | 317.8 us | 302.4 us | 76.6 us | 0.989x | ok `d2fc95c10ececedd` |
| `fq.q7e4.elem_muladd.n50000` | 759.3 us | 28.7 us | 758.3 us | 747.0 us | 33.2 us | 1.001x | ok `a6cf18a43d40c5c0` |
| `fq.q7e4.elem_inv.n5000` | 52.2 us | 3.2 us | 53.6 us | 51.6 us | 2.7 us | 0.974x | ok `8664ac6d09409f94` |
| `fq.q7e4.mul_school.d128` | 196.8 us | 30.7 us | 196.5 us | 157.1 us | 29.1 us | 1.001x | ok `71d186560b1d8153` |
| `fq.q7e4.mul_karatsuba.d512` | 1.318 ms | 214.9 us | 1.320 ms | 1.246 ms | 177.7 us | 0.999x | ok `6a2e9b751ab48e99` |
| `fq.q7e4.gcd.d200` | 592.0 us | 110.3 us | 576.1 us | 548.4 us | 126.9 us | 1.028x | ok `8cf19ac90e512feb` |
| `fq.q7e4.factor.d20` | 352.9 us | 102.4 us | 357.1 us | 304.2 us | 105.9 us | 0.988x | ok `ce23318b1b5fd48e` |
| `fq.m61e2.elem_muladd.n100000` | 23.326 ms | 655.8 us | 23.635 ms | 21.480 ms | 2.338 ms | 0.987x | ok `1814347c960518d8` |
| `fq.m61e2.elem_inv.n1000` | 950.6 us | 98.1 us | 958.7 us | 930.4 us | 97.2 us | 0.992x | ok `45326d8095032432` |
| `fq.m61e2.mul_school.d128` | 293.7 us | 19.4 us | 290.9 us | 279.5 us | 23.1 us | 1.010x | ok `2cb38bd639886e3c` |
| `fq.m61e2.mul_karatsuba.d512` | 3.448 ms | 342.5 us | 3.453 ms | 2.756 ms | 447.6 us | 0.999x | ok `f138ad86c0a7a838` |
| `fq.m61e2.mul_fft.d1024` | 36.946 ms | 1.910 ms | 36.929 ms | 35.778 ms | 1.206 ms | 1.000x | ok `a45642d67f5a9912` |
| `fq.m61e2.mul_fft.d4096` | 94.189 ms | 3.824 ms | 93.441 ms | 89.310 ms | 6.658 ms | 1.008x | ok `51844e8f95a9c0ab` |
| `fq.m61e2.mul_auto.d2048` | 6.617 ms | 775.1 us | 6.546 ms | 5.498 ms | 796.9 us | 1.011x | ok `76e66ad22fb3fa2e` |
| `fq.m61e2.divmod.d500_250` | 1.117 ms | 216.9 us | 1.088 ms | 1.069 ms | 98.0 us | 1.027x | ok `53509a111940ee2a` |
| `fq.m61e2.gcd.d150` | 3.404 ms | 525.9 us | 3.438 ms | 3.183 ms | 615.6 us | 0.990x | ok `bd8c3c98cc56df84` |
| `fq.m61e2.factor.d7` | 6.779 ms | 927.0 us | 6.402 ms | 5.877 ms | 961.3 us | 1.059x | ok `d61430d11ee187c3` |
| `fq.build.q3e30` | 248.7 us | 24.6 us | 246.7 us | 231.6 us | 100.2 us | 1.008x | ok `4973a3e567e2fc74` |
| `fq.build.q7e20` | 245.9 us | 65.5 us | 236.9 us | 223.1 us | 58.1 us | 1.038x | ok `9c377cd554aa108e` |
| `fast.ntt.dft.k10` | 3.084 ms | 487.0 us | 3.063 ms | 2.824 ms | 844.3 us | 1.007x | ok `b61f454e63ae6d8d` |
| `fast.ntt.dft.k11` | 3.143 ms | 595.0 us | 3.144 ms | 2.872 ms | 623.6 us | 1.000x | ok `311d4fed0f42034c` |
| `fast.ntt.dft.k12` | 3.227 ms | 672.1 us | 3.110 ms | 2.901 ms | 1.050 ms | 1.038x | ok `70691e178cfc0c9d` |
| `fast.ntt.idft.k12` | 3.165 ms | 466.0 us | 3.043 ms | 2.924 ms | 428.5 us | 1.040x | ok `e0677b59c4e33615` |
| `fast.ntt.toeplitz_lower_inverse.n1024` | 105.8 us | 18.5 us | 106.2 us | 96.9 us | 26.0 us | 0.996x | ok `7749eb8ff8cac3d0` |
| `fast.ntt.toeplitz_lower_inverse.n8192` | 2.231 ms | 390.8 us | 2.281 ms | 1.727 ms | 390.1 us | 0.978x | ok `0d31554b19d8b996` |
| `fast.ntt.toeplitz_lower_vec.n4096` | 422.8 us | 107.0 us | 402.1 us | 334.8 us | 126.1 us | 1.051x | ok `86e36ad2a1eacb05` |
| `fast.ntt.toeplitz_divmod.d2000_1000` | 6.405 ms | 488.4 us | 6.404 ms | 6.137 ms | 513.7 us | 1.000x | ok `b138cdcceebcc32e` |
| `fast.p7.toeplitz_lower_inverse.n1024` | 149.7 us | 49.1 us | 140.8 us | 131.5 us | 19.9 us | 1.063x | ok `5df08c92a814155a` |
| `fast.p7.toeplitz_divmod.d2000_1000` | 3.285 ms | 821.8 us | 3.109 ms | 2.939 ms | 417.7 us | 1.056x | ok `433d13716b2f95e8` |
