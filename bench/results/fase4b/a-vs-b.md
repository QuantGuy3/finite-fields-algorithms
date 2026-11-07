# Comparativa de benchmarks alcp

Base: **f4-f4**. Speedup = mediana(base) / mediana(variante); > 1 significa más rápido.
IQR = rango intercuartílico de las muestras. Las medias geométricas excluyen escenarios con error.

- **f4-f4**: modo interleave, 2026-11-07T18:45:43+01:00 (nota: fase 4 nivel B: original vs nivel A vs nivel B)  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -DALCP_RELINK=1`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-65-g0aaa877 (con cambios sin commit), sha256 8eb181a79710  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: sí
- **f4-b**: modo interleave, 2026-11-07T18:45:56+01:00 (nota: fase 4 nivel B: original vs nivel A vs nivel B)  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -DALCP_RELINK=10`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-66-g9e8e925 (con cambios sin commit), sha256 a89c6fdf61a3  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: sí

Checksums: todos los escenarios coinciden con la base.

## Media geométrica del speedup

| ámbito | f4-b |
|---|---|
| **global** | 1.253x (n=128) |
| grupo `fast` | 1.287x (n=10) |
| grupo `fp.m127` | 1.540x (n=14) |
| grupo `fp.m61` | 1.141x (n=14) |
| grupo `fp.ntt` | 1.101x (n=16) |
| grupo `fp.p101` | 1.491x (n=14) |
| grupo `fp.p256` | 1.221x (n=16) |
| grupo `fp.p7` | 1.469x (n=14) |
| grupo `fq.build` | 0.920x (n=2) |
| grupo `fq.m61e2` | 1.230x (n=10) |
| grupo `fq.q2e8` | 1.058x (n=6) |
| grupo `fq.q3e30` | 1.093x (n=6) |
| grupo `fq.q7e4` | 0.995x (n=6) |
| operación `build` | 0.920x (n=2) |
| operación `dft` | 0.949x (n=3) |
| operación `divmod` | 1.873x (n=8) |
| operación `elem_inv` | 0.977x (n=10) |
| operación `elem_muladd` | 0.988x (n=10) |
| operación `elem_pow` | 1.016x (n=6) |
| operación `factor` | 1.044x (n=10) |
| operación `gcd` | 0.930x (n=8) |
| operación `gcd_ext` | 1.311x (n=6) |
| operación `idft` | 0.960x (n=1) |
| operación `is_irreducible` | 3.409x (n=7) |
| operación `mul_auto` | 3.731x (n=7) |
| operación `mul_fft` | 1.007x (n=6) |
| operación `mul_karatsuba` | 0.994x (n=16) |
| operación `mul_school` | 0.984x (n=16) |
| operación `pow_mod` | 1.589x (n=6) |
| operación `toeplitz_divmod` | 1.097x (n=2) |
| operación `toeplitz_lower_inverse` | 2.450x (n=3) |
| operación `toeplitz_lower_vec` | 0.861x (n=1) |

## Detalle por escenario

| escenario | base mediana | base IQR | f4-b mediana | f4-b mín | f4-b IQR | speedup | checksum |
|---|---:|---:|---:|---:|---:|---:|---|
| `fp.p7.elem_muladd.n400000` | 1.723 ms | 55.2 us | 1.715 ms | 1.674 ms | 101.2 us | 1.005x | ok `cf4797a827e7f9ec` |
| `fp.p7.elem_inv.n100000` | 343.7 us | 11.6 us | 341.7 us | 334.2 us | 15.3 us | 1.006x | ok `7beb5a3238b919cb` |
| `fp.p7.elem_pow.n100000` | 963.0 us | 40.4 us | 948.9 us | 925.2 us | 47.2 us | 1.015x | ok `a7f8c429d0d81b73` |
| `fp.p7.mul_school.d64` | 69.2 us | 4.1 us | 70.9 us | 66.8 us | 8.2 us | 0.976x | ok `559357b5bc3c57c2` |
| `fp.p7.mul_school.d512` | 166.5 us | 16.0 us | 166.2 us | 161.5 us | 16.2 us | 1.002x | ok `5936b7480b151eb3` |
| `fp.p7.mul_karatsuba.d512` | 136.2 us | 21.7 us | 135.1 us | 122.3 us | 17.9 us | 1.009x | ok `d72c59b658ac0379` |
| `fp.p7.mul_karatsuba.d4096` | 1.422 ms | 174.6 us | 1.470 ms | 1.082 ms | 307.6 us | 0.968x | ok `00f3dde93cbd25a4` |
| `fp.p7.mul_auto.d2048` | 728.5 us | 154.4 us | 76.9 us | 73.4 us | 10.0 us | 9.473x | ok `2b256c747e0070c7` |
| `fp.p7.divmod.d2000_1000` | 670.3 us | 62.9 us | 150.8 us | 143.2 us | 13.6 us | 4.446x | ok `958d383cda71b699` |
| `fp.p7.gcd.d500` | 866.1 us | 76.3 us | 954.8 us | 889.9 us | 104.8 us | 0.907x | ok `cdf18334cf9b13d7` |
| `fp.p7.gcd_ext.d400` | 1.874 ms | 210.5 us | 1.148 ms | 1.028 ms | 136.6 us | 1.632x | ok `984ac156886a3f82` |
| `fp.p7.pow_mod.d128_k48` | 4.740 ms | 349.6 us | 1.783 ms | 1.670 ms | 164.9 us | 2.658x | ok `e3cedb8421cd962f` |
| `fp.p7.is_irreducible.d64` | 1.689 ms | 130.4 us | 1.276 ms | 1.227 ms | 99.9 us | 1.323x | ok `9cb1bec87bc433d0` |
| `fp.p7.factor.d62` | 293.1 us | 36.3 us | 288.9 us | 277.6 us | 25.6 us | 1.015x | ok `3425aecb393956ba` |
| `fp.p101.elem_muladd.n400000` | 1.724 ms | 60.3 us | 1.725 ms | 1.676 ms | 120.4 us | 1.000x | ok `e51aa870ce6269eb` |
| `fp.p101.elem_inv.n100000` | 683.3 us | 40.7 us | 695.6 us | 663.6 us | 61.3 us | 0.982x | ok `69d137d626075364` |
| `fp.p101.elem_pow.n100000` | 2.176 ms | 180.4 us | 1.991 ms | 1.930 ms | 248.4 us | 1.093x | ok `09fec44561c97afa` |
| `fp.p101.mul_school.d64` | 78.9 us | 6.3 us | 76.6 us | 73.9 us | 6.9 us | 1.030x | ok `5a0fa852855ea8b0` |
| `fp.p101.mul_school.d512` | 201.5 us | 12.6 us | 198.8 us | 190.0 us | 24.1 us | 1.014x | ok `8289633c899317a1` |
| `fp.p101.mul_karatsuba.d512` | 130.9 us | 11.7 us | 129.3 us | 120.2 us | 16.5 us | 1.012x | ok `fea453dad54c0188` |
| `fp.p101.mul_karatsuba.d4096` | 1.460 ms | 193.1 us | 1.491 ms | 1.197 ms | 288.0 us | 0.979x | ok `82f83f0b7fe7713b` |
| `fp.p101.mul_auto.d2048` | 763.8 us | 165.0 us | 109.5 us | 99.6 us | 13.1 us | 6.979x | ok `3631fff5520b84d9` |
| `fp.p101.divmod.d2000_1000` | 804.6 us | 95.1 us | 206.2 us | 184.5 us | 25.3 us | 3.901x | ok `ce32beed1dda91cd` |
| `fp.p101.gcd.d500` | 1.129 ms | 104.0 us | 1.239 ms | 1.075 ms | 191.4 us | 0.911x | ok `5cad0c7c5f168edb` |
| `fp.p101.gcd_ext.d400` | 2.364 ms | 183.4 us | 1.336 ms | 1.219 ms | 198.5 us | 1.770x | ok `4fe0534b04a103a1` |
| `fp.p101.pow_mod.d128_k20` | 5.619 ms | 319.1 us | 2.117 ms | 1.940 ms | 172.0 us | 2.655x | ok `bc5334d07b025488` |
| `fp.p101.is_irreducible.d48` | 1.584 ms | 92.4 us | 876.0 us | 806.0 us | 102.1 us | 1.808x | ok `5146dc1467eee05b` |
| `fp.p101.factor.d100` | 1.525 ms | 175.3 us | 1.333 ms | 1.228 ms | 197.2 us | 1.145x | ok `48a7604aa4ae63ee` |
| `fp.ntt.elem_muladd.n400000` | 2.550 ms | 102.0 us | 2.522 ms | 2.452 ms | 199.7 us | 1.011x | ok `682b833fe4edd17b` |
| `fp.ntt.elem_inv.n100000` | 3.127 ms | 217.6 us | 3.024 ms | 2.932 ms | 216.7 us | 1.034x | ok `d275b741e6b3050e` |
| `fp.ntt.elem_pow.n50000` | 6.878 ms | 453.4 us | 6.783 ms | 6.670 ms | 494.4 us | 1.014x | ok `c27b88ab2935c700` |
| `fp.ntt.mul_school.d64` | 79.4 us | 6.2 us | 78.8 us | 75.8 us | 8.2 us | 1.008x | ok `67c7e1645a9cc6ae` |
| `fp.ntt.mul_school.d512` | 201.4 us | 18.6 us | 221.3 us | 194.3 us | 26.1 us | 0.910x | ok `e296393a6a637046` |
| `fp.ntt.mul_karatsuba.d512` | 134.0 us | 13.7 us | 134.8 us | 121.6 us | 15.5 us | 0.994x | ok `964e15063fb16df4` |
| `fp.ntt.mul_karatsuba.d4096` | 1.720 ms | 332.2 us | 1.687 ms | 1.381 ms | 213.7 us | 1.020x | ok `d656b678f0c2ca16` |
| `fp.ntt.mul_fft.d4096` | 1.732 ms | 211.9 us | 1.705 ms | 1.447 ms | 253.4 us | 1.016x | ok `f4c81b85220fb683` |
| `fp.ntt.mul_fft.d32768` | 16.010 ms | 1.074 ms | 15.668 ms | 14.252 ms | 2.009 ms | 1.022x | ok `8228e76950faded3` |
| `fp.ntt.mul_auto.d2048` | 761.4 us | 97.4 us | 817.3 us | 694.0 us | 114.2 us | 0.932x | ok `02e67c888996018d` |
| `fp.ntt.divmod.d2000_1000` | 799.2 us | 33.8 us | 551.4 us | 508.6 us | 63.6 us | 1.449x | ok `7adfbb2614325099` |
| `fp.ntt.gcd.d500` | 1.195 ms | 60.0 us | 1.540 ms | 1.363 ms | 233.9 us | 0.776x | ok `678962135e7a2447` |
| `fp.ntt.gcd_ext.d400` | 2.482 ms | 225.1 us | 2.580 ms | 2.314 ms | 411.6 us | 0.962x | ok `f9dc4e400a7954ee` |
| `fp.ntt.pow_mod.d128_k5` | 5.239 ms | 419.1 us | 4.176 ms | 3.845 ms | 602.1 us | 1.255x | ok `0eb2514aceb37261` |
| `fp.ntt.is_irreducible.d32` | 1.571 ms | 260.9 us | 367.8 us | 336.8 us | 44.9 us | 4.272x | ok `1f2827a0d8aeb65e` |
| `fp.ntt.factor.d40` | 1.167 ms | 110.5 us | 1.381 ms | 1.245 ms | 213.3 us | 0.845x | ok `297c76c423fe55aa` |
| `fp.m61.elem_muladd.n400000` | 2.260 ms | 78.6 us | 2.311 ms | 2.173 ms | 359.6 us | 0.978x | ok `b8e325f73ad90582` |
| `fp.m61.elem_inv.n100000` | 8.426 ms | 152.1 us | 9.308 ms | 8.371 ms | 1.070 ms | 0.905x | ok `41e7fe38de6fdcc1` |
| `fp.m61.elem_pow.n20000` | 6.189 ms | 118.2 us | 6.432 ms | 6.047 ms | 646.2 us | 0.962x | ok `a3680ece7fc07932` |
| `fp.m61.mul_school.d64` | 84.4 us | 7.7 us | 85.1 us | 78.7 us | 11.7 us | 0.992x | ok `bed6e9e86907de79` |
| `fp.m61.mul_school.d512` | 219.7 us | 25.1 us | 216.9 us | 199.8 us | 35.4 us | 1.013x | ok `43c3884710f7d468` |
| `fp.m61.mul_karatsuba.d512` | 138.0 us | 12.2 us | 144.1 us | 132.9 us | 30.9 us | 0.958x | ok `a2dcbeaa6574fcb8` |
| `fp.m61.mul_karatsuba.d4096` | 1.678 ms | 182.0 us | 1.666 ms | 1.207 ms | 340.7 us | 1.008x | ok `37d5deaa1cfeb600` |
| `fp.m61.mul_auto.d2048` | 920.1 us | 115.1 us | 595.4 us | 527.8 us | 343.3 us | 1.545x | ok `73174ffa4e90d6fa` |
| `fp.m61.divmod.d2000_1000` | 862.4 us | 87.6 us | 857.3 us | 787.4 us | 124.1 us | 1.006x | ok `a8710fc5aaa24d09` |
| `fp.m61.gcd.d500` | 1.040 ms | 86.9 us | 1.305 ms | 1.192 ms | 214.7 us | 0.797x | ok `fda394dde025c1a2` |
| `fp.m61.gcd_ext.d400` | 2.196 ms | 263.1 us | 2.420 ms | 2.189 ms | 537.2 us | 0.908x | ok `be3d81e8acc44b48` |
| `fp.m61.pow_mod.d64_k2` | 1.565 ms | 166.8 us | 1.607 ms | 1.470 ms | 216.6 us | 0.974x | ok `394fb320787700fe` |
| `fp.m61.is_irreducible.d24` | 2.520 ms | 174.0 us | 335.4 us | 311.5 us | 43.4 us | 7.512x | ok `9047cf1d3fc29abf` |
| `fp.m61.factor.d26` | 1.243 ms | 99.2 us | 1.334 ms | 1.194 ms | 165.9 us | 0.931x | ok `931b8d0b8b91b0b4` |
| `fp.m127.elem_muladd.n200000` | 24.389 ms | 883.1 us | 24.916 ms | 22.393 ms | 3.517 ms | 0.979x | ok `16fab7b139aa9aab` |
| `fp.m127.elem_inv.n50000` | 31.465 ms | 1.330 ms | 32.190 ms | 29.439 ms | 4.844 ms | 0.978x | ok `c7d67d6018036bf1` |
| `fp.m127.elem_pow.n10000` | 24.565 ms | 1.740 ms | 24.647 ms | 22.347 ms | 2.678 ms | 0.997x | ok `308896e55b40c815` |
| `fp.m127.mul_school.d64` | 1.765 ms | 180.6 us | 1.711 ms | 1.577 ms | 178.9 us | 1.032x | ok `a9664408b0fb88e7` |
| `fp.m127.mul_school.d512` | 4.063 ms | 421.9 us | 4.221 ms | 3.760 ms | 706.3 us | 0.963x | ok `98eca357802203c0` |
| `fp.m127.mul_karatsuba.d512` | 1.873 ms | 210.4 us | 1.881 ms | 1.203 ms | 304.9 us | 0.996x | ok `f773c153fde34172` |
| `fp.m127.mul_karatsuba.d2048` | 10.761 ms | 800.5 us | 10.793 ms | 8.726 ms | 1.109 ms | 0.997x | ok `695e0d0f13c14e87` |
| `fp.m127.mul_auto.d2048` | 10.659 ms | 588.0 us | 1.936 ms | 1.622 ms | 507.1 us | 5.505x | ok `1bc667b834b9e559` |
| `fp.m127.divmod.d2000_1000` | 15.288 ms | 1.083 ms | 4.508 ms | 3.977 ms | 684.5 us | 3.391x | ok `a8babdcd119eb404` |
| `fp.m127.gcd.d400` | 14.993 ms | 754.6 us | 13.266 ms | 11.408 ms | 1.707 ms | 1.130x | ok `1b15fac38c54accc` |
| `fp.m127.gcd_ext.d300` | 24.954 ms | 1.974 ms | 13.643 ms | 12.214 ms | 1.729 ms | 1.829x | ok `45a1b9001453c6ef` |
| `fp.m127.pow_mod.d64_k1` | 47.703 ms | 8.195 ms | 28.520 ms | 26.292 ms | 3.359 ms | 1.673x | ok `9e9e535b3c4e4e7d` |
| `fp.m127.is_irreducible.d16` | 40.535 ms | 5.466 ms | 5.466 ms | 4.974 ms | 328.9 us | 7.416x | ok `2cfb69ee2746aa9d` |
| `fp.m127.factor.d19` | 37.558 ms | 11.121 ms | 39.929 ms | 35.538 ms | 4.616 ms | 0.941x | ok `7b0defab2417737c` |
| `fp.p256.elem_muladd.n200000` | 33.144 ms | 5.754 ms | 34.409 ms | 31.364 ms | 1.985 ms | 0.963x | ok `cc3cab838e9b19f1` |
| `fp.p256.elem_inv.n50000` | 69.600 ms | 10.446 ms | 69.267 ms | 65.136 ms | 4.389 ms | 1.005x | ok `89f920a8eba0c743` |
| `fp.p256.elem_pow.n2000` | 20.007 ms | 3.249 ms | 19.673 ms | 18.902 ms | 1.685 ms | 1.017x | ok `40472b256a3c877d` |
| `fp.p256.mul_school.d64` | 2.539 ms | 640.1 us | 2.561 ms | 2.407 ms | 322.8 us | 0.991x | ok `f5d3d87e3ed0977e` |
| `fp.p256.mul_school.d512` | 6.618 ms | 3.266 ms | 6.939 ms | 6.058 ms | 429.5 us | 0.954x | ok `6fd19bd7b3ceab69` |
| `fp.p256.mul_karatsuba.d512` | 2.094 ms | 240.4 us | 2.105 ms | 1.488 ms | 300.9 us | 0.995x | ok `4b2f43bf0f878408` |
| `fp.p256.mul_karatsuba.d2048` | 13.506 ms | 1.275 ms | 13.812 ms | 10.502 ms | 1.470 ms | 0.978x | ok `b4b89a80c2a3aaa9` |
| `fp.p256.mul_fft.d4096` | 18.072 ms | 4.706 ms | 19.003 ms | 15.300 ms | 4.090 ms | 0.951x | ok `60ed510564513dfb` |
| `fp.p256.mul_fft.d8192` | 38.900 ms | 9.110 ms | 36.742 ms | 31.956 ms | 10.599 ms | 1.059x | ok `47812935df4de5fe` |
| `fp.p256.mul_auto.d2048` | 8.747 ms | 917.6 us | 4.082 ms | 3.217 ms | 919.3 us | 2.143x | ok `3d89402081500506` |
| `fp.p256.divmod.d1000_500` | 6.147 ms | 531.7 us | 3.437 ms | 3.134 ms | 376.4 us | 1.788x | ok `0bebce51d98eb3f4` |
| `fp.p256.gcd.d300` | 10.902 ms | 2.409 ms | 11.106 ms | 10.346 ms | 1.272 ms | 0.982x | ok `f69e4050cb1a5861` |
| `fp.p256.gcd_ext.d200` | 13.883 ms | 3.378 ms | 12.638 ms | 11.880 ms | 604.1 us | 1.099x | ok `1a9df499c28b9f49` |
| `fp.p256.pow_mod.d64_k1` | 79.660 ms | 5.504 ms | 71.379 ms | 66.040 ms | 6.654 ms | 1.116x | ok `dbe2b595567fbb9c` |
| `fp.p256.is_irreducible.d12` | 41.679 ms | 3.297 ms | 7.259 ms | 6.673 ms | 709.4 us | 5.741x | ok `7b839a40c98b0ef0` |
| `fp.p256.factor.d12` | 17.777 ms | 462.9 us | 17.662 ms | 16.748 ms | 1.222 ms | 1.007x | ok `c33e3512ed388a7e` |
| `fq.q3e30.elem_muladd.n4000` | 5.900 ms | 557.4 us | 5.857 ms | 5.505 ms | 500.7 us | 1.007x | ok `16edfcc06a08028a` |
| `fq.q3e30.elem_inv.n100` | 754.6 us | 67.1 us | 784.9 us | 717.2 us | 79.0 us | 0.961x | ok `a8167c01d1cfc5dd` |
| `fq.q3e30.mul_school.d32` | 734.1 us | 86.7 us | 767.8 us | 707.7 us | 72.4 us | 0.956x | ok `ed0ca8eb7d325de3` |
| `fq.q3e30.mul_karatsuba.d128` | 2.335 ms | 179.3 us | 2.265 ms | 1.520 ms | 369.8 us | 1.031x | ok `2464b901bd242c7d` |
| `fq.q3e30.divmod.d64_32` | 776.1 us | 68.7 us | 783.7 us | 722.9 us | 104.4 us | 0.990x | ok `a2b35f414a139741` |
| `fq.q3e30.factor.d5` | 6.980 ms | 644.3 us | 3.864 ms | 3.585 ms | 391.0 us | 1.806x | ok `417a5b753c7a2a07` |
| `fq.q2e8.elem_muladd.n50000` | 9.411 ms | 889.2 us | 9.813 ms | 8.559 ms | 1.159 ms | 0.959x | ok `c0ffe16ee54f128c` |
| `fq.q2e8.elem_inv.n5000` | 2.511 ms | 184.5 us | 2.511 ms | 2.240 ms | 285.4 us | 1.000x | ok `81d71b5687e3e735` |
| `fq.q2e8.mul_school.d128` | 618.5 us | 47.4 us | 619.2 us | 567.5 us | 105.2 us | 0.999x | ok `c056d9d992271db9` |
| `fq.q2e8.mul_karatsuba.d512` | 2.301 ms | 215.1 us | 2.377 ms | 1.643 ms | 373.7 us | 0.968x | ok `de1089fdce26d02c` |
| `fq.q2e8.is_irreducible.d16` | 2.558 ms | 257.5 us | 1.563 ms | 1.431 ms | 269.2 us | 1.636x | ok `9478d590a29a9cbf` |
| `fq.q2e8.factor.d28` | 2.536 ms | 250.5 us | 2.737 ms | 2.535 ms | 324.4 us | 0.926x | ok `d2fc95c10ececedd` |
| `fq.q7e4.elem_muladd.n50000` | 6.266 ms | 685.3 us | 6.288 ms | 5.911 ms | 714.2 us | 0.997x | ok `a6cf18a43d40c5c0` |
| `fq.q7e4.elem_inv.n5000` | 1.847 ms | 173.4 us | 1.986 ms | 1.786 ms | 223.5 us | 0.930x | ok `8664ac6d09409f94` |
| `fq.q7e4.mul_school.d128` | 294.1 us | 31.1 us | 304.6 us | 267.5 us | 58.2 us | 0.966x | ok `71d186560b1d8153` |
| `fq.q7e4.mul_karatsuba.d512` | 1.757 ms | 209.5 us | 1.807 ms | 1.214 ms | 532.3 us | 0.972x | ok `6a2e9b751ab48e99` |
| `fq.q7e4.gcd.d200` | 3.586 ms | 506.6 us | 3.612 ms | 3.189 ms | 826.8 us | 0.993x | ok `8cf19ac90e512feb` |
| `fq.q7e4.factor.d20` | 1.482 ms | 109.8 us | 1.318 ms | 1.206 ms | 236.4 us | 1.124x | ok `ce23318b1b5fd48e` |
| `fq.m61e2.elem_muladd.n100000` | 9.309 ms | 1.143 ms | 9.449 ms | 8.665 ms | 2.098 ms | 0.985x | ok `1814347c960518d8` |
| `fq.m61e2.elem_inv.n1000` | 402.1 us | 32.6 us | 411.8 us | 385.0 us | 38.2 us | 0.976x | ok `45326d8095032432` |
| `fq.m61e2.mul_school.d128` | 121.8 us | 11.0 us | 127.6 us | 112.2 us | 12.0 us | 0.955x | ok `2cb38bd639886e3c` |
| `fq.m61e2.mul_karatsuba.d512` | 1.434 ms | 240.9 us | 1.389 ms | 907.1 us | 275.1 us | 1.032x | ok `f138ad86c0a7a838` |
| `fq.m61e2.mul_fft.d1024` | 15.645 ms | 1.774 ms | 15.416 ms | 12.089 ms | 1.711 ms | 1.015x | ok `a45642d67f5a9912` |
| `fq.m61e2.mul_fft.d4096` | 42.433 ms | 3.327 ms | 43.174 ms | 41.416 ms | 4.060 ms | 0.983x | ok `51844e8f95a9c0ab` |
| `fq.m61e2.mul_auto.d2048` | 26.664 ms | 2.988 ms | 2.977 ms | 2.590 ms | 494.9 us | 8.956x | ok `76e66ad22fb3fa2e` |
| `fq.m61e2.divmod.d500_250` | 479.8 us | 54.2 us | 480.8 us | 448.3 us | 28.1 us | 0.998x | ok `53509a111940ee2a` |
| `fq.m61e2.gcd.d150` | 1.391 ms | 128.4 us | 1.403 ms | 1.302 ms | 113.1 us | 0.992x | ok `bd8c3c98cc56df84` |
| `fq.m61e2.factor.d7` | 2.938 ms | 277.3 us | 3.117 ms | 2.934 ms | 213.1 us | 0.942x | ok `d61430d11ee187c3` |
| `fq.build.q3e30` | 246.1 us | 31.0 us | 287.2 us | 257.7 us | 57.2 us | 0.857x | ok `4973a3e567e2fc74` |
| `fq.build.q7e20` | 213.6 us | 26.9 us | 216.4 us | 192.3 us | 29.4 us | 0.987x | ok `9c377cd554aa108e` |
| `fast.ntt.dft.k10` | 1.272 ms | 142.2 us | 1.359 ms | 1.225 ms | 255.5 us | 0.936x | ok `b61f454e63ae6d8d` |
| `fast.ntt.dft.k11` | 1.292 ms | 177.9 us | 1.305 ms | 1.240 ms | 342.4 us | 0.989x | ok `311d4fed0f42034c` |
| `fast.ntt.dft.k12` | 1.261 ms | 111.5 us | 1.368 ms | 1.232 ms | 177.0 us | 0.922x | ok `70691e178cfc0c9d` |
| `fast.ntt.idft.k12` | 1.287 ms | 136.2 us | 1.340 ms | 1.247 ms | 155.0 us | 0.960x | ok `e0677b59c4e33615` |
| `fast.ntt.toeplitz_lower_inverse.n1024` | 610.4 us | 41.9 us | 383.0 us | 342.7 us | 64.5 us | 1.594x | ok `7749eb8ff8cac3d0` |
| `fast.ntt.toeplitz_lower_inverse.n8192` | 6.474 ms | 639.4 us | 6.461 ms | 5.535 ms | 1.580 ms | 1.002x | ok `0d31554b19d8b996` |
| `fast.ntt.toeplitz_lower_vec.n4096` | 812.2 us | 82.8 us | 943.1 us | 746.7 us | 246.4 us | 0.861x | ok `86e36ad2a1eacb05` |
| `fast.ntt.toeplitz_divmod.d2000_1000` | 2.922 ms | 269.3 us | 2.866 ms | 2.741 ms | 417.8 us | 1.019x | ok `b138cdcceebcc32e` |
| `fast.p7.toeplitz_lower_inverse.n1024` | 1.005 ms | 112.8 us | 109.2 us | 96.6 us | 19.9 us | 9.203x | ok `5df08c92a814155a` |
| `fast.p7.toeplitz_divmod.d2000_1000` | 1.630 ms | 186.5 us | 1.382 ms | 1.246 ms | 144.3 us | 1.180x | ok `433d13716b2f95e8` |
