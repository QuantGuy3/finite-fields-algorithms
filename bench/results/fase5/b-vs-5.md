# Comparativa de benchmarks alcp

Base: **f5-b**. Speedup = mediana(base) / mediana(variante); > 1 significa más rápido.
IQR = rango intercuartílico de las muestras. Las medias geométricas excluyen escenarios con error.

- **f5-b**: modo interleave, 2026-11-07T21:21:14+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -DSUITE_R2`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\build\bench\inc-nivelB\include` @ baseline-original-68-g20f9a31, sha256 f2a9fe2a61e0  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: no (batería)
- **f5-5**: modo interleave, 2026-11-07T21:21:41+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -DSUITE_R2`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-68-g20f9a31 (con cambios sin commit), sha256 ea3a40d2b2fc  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: no (batería)

Checksums: todos los escenarios coinciden con la base.

## Media geométrica del speedup

| ámbito | f5-5 |
|---|---|
| **global** | 1.139x (n=128) |
| grupo `fast` | 1.215x (n=10) |
| grupo `fp.m127` | 1.009x (n=14) |
| grupo `fp.m61` | 1.214x (n=14) |
| grupo `fp.ntt` | 1.423x (n=16) |
| grupo `fp.p101` | 1.230x (n=14) |
| grupo `fp.p256` | 0.988x (n=16) |
| grupo `fp.p7` | 1.204x (n=14) |
| grupo `fq.build` | 1.540x (n=2) |
| grupo `fq.m61e2` | 1.006x (n=10) |
| grupo `fq.q2e8` | 1.017x (n=6) |
| grupo `fq.q3e30` | 0.983x (n=6) |
| grupo `fq.q7e4` | 1.007x (n=6) |
| operación `build` | 1.540x (n=2) |
| operación `dft` | 0.991x (n=3) |
| operación `divmod` | 1.012x (n=8) |
| operación `elem_inv` | 1.008x (n=10) |
| operación `elem_muladd` | 1.012x (n=10) |
| operación `elem_pow` | 1.015x (n=6) |
| operación `factor` | 1.036x (n=10) |
| operación `gcd` | 1.988x (n=8) |
| operación `gcd_ext` | 1.909x (n=6) |
| operación `idft` | 1.016x (n=1) |
| operación `is_irreducible` | 1.309x (n=7) |
| operación `mul_auto` | 1.104x (n=7) |
| operación `mul_fft` | 1.297x (n=6) |
| operación `mul_karatsuba` | 0.983x (n=16) |
| operación `mul_school` | 0.987x (n=16) |
| operación `pow_mod` | 1.020x (n=6) |
| operación `toeplitz_divmod` | 1.015x (n=2) |
| operación `toeplitz_lower_inverse` | 1.485x (n=3) |
| operación `toeplitz_lower_vec` | 2.108x (n=1) |

## Detalle por escenario

| escenario | base mediana | base IQR | f5-5 mediana | f5-5 mín | f5-5 IQR | speedup | checksum |
|---|---:|---:|---:|---:|---:|---:|---|
| `fp.p7.elem_muladd.n400000` | 4.086 ms | 99.5 us | 4.118 ms | 4.062 ms | 211.6 us | 0.992x | ok `cf4797a827e7f9ec` |
| `fp.p7.elem_inv.n100000` | 815.4 us | 17.9 us | 815.5 us | 808.0 us | 63.3 us | 1.000x | ok `7beb5a3238b919cb` |
| `fp.p7.elem_pow.n100000` | 2.303 ms | 222.0 us | 2.268 ms | 2.245 ms | 414.5 us | 1.015x | ok `a7f8c429d0d81b73` |
| `fp.p7.mul_school.d64` | 169.4 us | 12.6 us | 175.6 us | 167.5 us | 11.4 us | 0.965x | ok `559357b5bc3c57c2` |
| `fp.p7.mul_school.d512` | 393.2 us | 55.3 us | 408.5 us | 406.6 us | 42.1 us | 0.963x | ok `5936b7480b151eb3` |
| `fp.p7.mul_karatsuba.d512` | 346.0 us | 137.3 us | 334.5 us | 302.4 us | 72.3 us | 1.034x | ok `d72c59b658ac0379` |
| `fp.p7.mul_karatsuba.d4096` | 4.390 ms | 1.286 ms | 4.655 ms | 3.026 ms | 1.547 ms | 0.943x | ok `00f3dde93cbd25a4` |
| `fp.p7.mul_auto.d2048` | 178.2 us | 7.5 us | 178.2 us | 177.4 us | 7.0 us | 1.000x | ok `2b256c747e0070c7` |
| `fp.p7.divmod.d2000_1000` | 377.6 us | 116.9 us | 400.9 us | 343.8 us | 145.4 us | 0.942x | ok `958d383cda71b699` |
| `fp.p7.gcd.d500` | 2.350 ms | 406.7 us | 623.0 us | 611.2 us | 44.7 us | 3.772x | ok `cdf18334cf9b13d7` |
| `fp.p7.gcd_ext.d400` | 2.667 ms | 436.2 us | 1.393 ms | 1.346 ms | 181.8 us | 1.915x | ok `984ac156886a3f82` |
| `fp.p7.pow_mod.d128_k48` | 4.189 ms | 649.7 us | 4.183 ms | 3.956 ms | 754.5 us | 1.002x | ok `e3cedb8421cd962f` |
| `fp.p7.is_irreducible.d64` | 3.088 ms | 578.6 us | 1.510 ms | 1.454 ms | 186.7 us | 2.045x | ok `9cb1bec87bc433d0` |
| `fp.p7.factor.d62` | 696.5 us | 115.9 us | 659.7 us | 599.9 us | 112.7 us | 1.056x | ok `3425aecb393956ba` |
| `fp.p101.elem_muladd.n400000` | 4.150 ms | 205.4 us | 4.119 ms | 4.062 ms | 278.2 us | 1.008x | ok `e51aa870ce6269eb` |
| `fp.p101.elem_inv.n100000` | 1.621 ms | 114.9 us | 1.615 ms | 1.606 ms | 72.4 us | 1.004x | ok `69d137d626075364` |
| `fp.p101.elem_pow.n100000` | 4.985 ms | 1.032 ms | 4.794 ms | 4.680 ms | 524.0 us | 1.040x | ok `09fec44561c97afa` |
| `fp.p101.mul_school.d64` | 180.0 us | 8.3 us | 186.7 us | 186.1 us | 7.3 us | 0.964x | ok `5a0fa852855ea8b0` |
| `fp.p101.mul_school.d512` | 466.1 us | 106.2 us | 486.0 us | 478.7 us | 99.3 us | 0.959x | ok `8289633c899317a1` |
| `fp.p101.mul_karatsuba.d512` | 318.5 us | 106.5 us | 314.7 us | 295.2 us | 73.3 us | 1.012x | ok `fea453dad54c0188` |
| `fp.p101.mul_karatsuba.d4096` | 4.059 ms | 872.8 us | 3.655 ms | 2.790 ms | 1.022 ms | 1.111x | ok `82f83f0b7fe7713b` |
| `fp.p101.mul_auto.d2048` | 246.8 us | 48.0 us | 240.4 us | 238.9 us | 11.3 us | 1.027x | ok `3631fff5520b84d9` |
| `fp.p101.divmod.d2000_1000` | 469.8 us | 61.0 us | 459.5 us | 447.7 us | 30.2 us | 1.022x | ok `ce32beed1dda91cd` |
| `fp.p101.gcd.d500` | 2.844 ms | 420.3 us | 741.0 us | 723.2 us | 37.7 us | 3.838x | ok `5cad0c7c5f168edb` |
| `fp.p101.gcd_ext.d400` | 3.140 ms | 696.5 us | 1.588 ms | 1.540 ms | 136.8 us | 1.977x | ok `4fe0534b04a103a1` |
| `fp.p101.pow_mod.d128_k20` | 4.975 ms | 566.5 us | 4.944 ms | 4.696 ms | 430.6 us | 1.006x | ok `bc5334d07b025488` |
| `fp.p101.is_irreducible.d48` | 2.062 ms | 373.7 us | 1.166 ms | 1.129 ms | 104.8 us | 1.768x | ok `5146dc1467eee05b` |
| `fp.p101.factor.d100` | 3.168 ms | 614.3 us | 2.712 ms | 2.619 ms | 411.6 us | 1.168x | ok `48a7604aa4ae63ee` |
| `fp.ntt.elem_muladd.n400000` | 6.138 ms | 364.4 us | 6.065 ms | 5.953 ms | 380.5 us | 1.012x | ok `682b833fe4edd17b` |
| `fp.ntt.elem_inv.n100000` | 7.621 ms | 974.4 us | 7.583 ms | 7.201 ms | 641.3 us | 1.005x | ok `d275b741e6b3050e` |
| `fp.ntt.elem_pow.n50000` | 16.881 ms | 595.1 us | 16.885 ms | 16.287 ms | 632.1 us | 1.000x | ok `c27b88ab2935c700` |
| `fp.ntt.mul_school.d64` | 184.3 us | 7.2 us | 191.7 us | 191.1 us | 7.1 us | 0.961x | ok `67c7e1645a9cc6ae` |
| `fp.ntt.mul_school.d512` | 469.2 us | 56.2 us | 486.4 us | 484.0 us | 52.5 us | 0.965x | ok `e296393a6a637046` |
| `fp.ntt.mul_karatsuba.d512` | 314.8 us | 43.5 us | 313.1 us | 295.9 us | 48.2 us | 1.005x | ok `964e15063fb16df4` |
| `fp.ntt.mul_karatsuba.d4096` | 3.688 ms | 1.089 ms | 3.968 ms | 2.966 ms | 1.029 ms | 0.929x | ok `d656b678f0c2ca16` |
| `fp.ntt.mul_fft.d4096` | 3.780 ms | 368.7 us | 1.731 ms | 1.484 ms | 212.8 us | 2.184x | ok `f4c81b85220fb683` |
| `fp.ntt.mul_fft.d32768` | 35.542 ms | 740.9 us | 15.574 ms | 14.771 ms | 654.4 us | 2.282x | ok `8228e76950faded3` |
| `fp.ntt.mul_auto.d2048` | 1.919 ms | 208.3 us | 950.2 us | 774.4 us | 91.2 us | 2.019x | ok `02e67c888996018d` |
| `fp.ntt.divmod.d2000_1000` | 1.323 ms | 368.6 us | 1.190 ms | 1.156 ms | 136.9 us | 1.111x | ok `7adfbb2614325099` |
| `fp.ntt.gcd.d500` | 3.482 ms | 404.7 us | 786.4 us | 763.7 us | 89.8 us | 4.429x | ok `678962135e7a2447` |
| `fp.ntt.gcd_ext.d400` | 6.046 ms | 860.3 us | 1.635 ms | 1.567 ms | 202.5 us | 3.698x | ok `f9dc4e400a7954ee` |
| `fp.ntt.pow_mod.d128_k5` | 10.353 ms | 1.677 ms | 9.109 ms | 8.467 ms | 1.068 ms | 1.137x | ok `0eb2514aceb37261` |
| `fp.ntt.is_irreducible.d32` | 839.8 us | 125.8 us | 574.2 us | 557.7 us | 82.5 us | 1.462x | ok `1f2827a0d8aeb65e` |
| `fp.ntt.factor.d40` | 3.203 ms | 788.5 us | 3.052 ms | 2.931 ms | 570.7 us | 1.049x | ok `297c76c423fe55aa` |
| `fp.m61.elem_muladd.n400000` | 5.464 ms | 452.3 us | 5.298 ms | 5.237 ms | 341.9 us | 1.031x | ok `b8e325f73ad90582` |
| `fp.m61.elem_inv.n100000` | 20.927 ms | 1.095 ms | 20.686 ms | 20.098 ms | 333.8 us | 1.012x | ok `41e7fe38de6fdcc1` |
| `fp.m61.elem_pow.n20000` | 15.446 ms | 1.006 ms | 15.027 ms | 14.601 ms | 506.4 us | 1.028x | ok `a3680ece7fc07932` |
| `fp.m61.mul_school.d64` | 191.8 us | 15.3 us | 198.6 us | 198.0 us | 7.7 us | 0.966x | ok `bed6e9e86907de79` |
| `fp.m61.mul_school.d512` | 483.7 us | 23.6 us | 504.3 us | 496.7 us | 32.4 us | 0.959x | ok `43c3884710f7d468` |
| `fp.m61.mul_karatsuba.d512` | 327.1 us | 64.4 us | 325.0 us | 308.8 us | 46.3 us | 1.006x | ok `a2dcbeaa6574fcb8` |
| `fp.m61.mul_karatsuba.d4096` | 4.108 ms | 1.016 ms | 3.955 ms | 2.874 ms | 1.217 ms | 1.039x | ok `37d5deaa1cfeb600` |
| `fp.m61.mul_auto.d2048` | 1.587 ms | 368.3 us | 1.655 ms | 1.503 ms | 282.4 us | 0.959x | ok `73174ffa4e90d6fa` |
| `fp.m61.divmod.d2000_1000` | 1.929 ms | 232.3 us | 1.937 ms | 1.900 ms | 250.5 us | 0.996x | ok `a8710fc5aaa24d09` |
| `fp.m61.gcd.d500` | 3.083 ms | 579.0 us | 834.4 us | 810.8 us | 53.9 us | 3.696x | ok `fda394dde025c1a2` |
| `fp.m61.gcd_ext.d400` | 5.849 ms | 1.082 ms | 1.669 ms | 1.610 ms | 197.1 us | 3.505x | ok `be3d81e8acc44b48` |
| `fp.m61.pow_mod.d64_k2` | 3.720 ms | 862.6 us | 3.738 ms | 3.581 ms | 727.1 us | 0.995x | ok `394fb320787700fe` |
| `fp.m61.is_irreducible.d24` | 779.5 us | 106.2 us | 669.5 us | 655.4 us | 33.5 us | 1.164x | ok `9047cf1d3fc29abf` |
| `fp.m61.factor.d26` | 3.001 ms | 441.6 us | 2.958 ms | 2.848 ms | 374.7 us | 1.015x | ok `931b8d0b8b91b0b4` |
| `fp.m127.elem_muladd.n200000` | 61.799 ms | 3.534 ms | 59.233 ms | 56.466 ms | 2.770 ms | 1.043x | ok `16fab7b139aa9aab` |
| `fp.m127.elem_inv.n50000` | 81.560 ms | 7.206 ms | 80.073 ms | 74.378 ms | 4.368 ms | 1.019x | ok `c7d67d6018036bf1` |
| `fp.m127.elem_pow.n10000` | 60.308 ms | 5.714 ms | 57.659 ms | 55.619 ms | 4.207 ms | 1.046x | ok `308896e55b40c815` |
| `fp.m127.mul_school.d64` | 4.115 ms | 964.9 us | 3.950 ms | 3.739 ms | 714.1 us | 1.042x | ok `a9664408b0fb88e7` |
| `fp.m127.mul_school.d512` | 10.034 ms | 1.313 ms | 9.864 ms | 8.917 ms | 1.448 ms | 1.017x | ok `98eca357802203c0` |
| `fp.m127.mul_karatsuba.d512` | 4.162 ms | 895.6 us | 4.385 ms | 3.365 ms | 1.105 ms | 0.949x | ok `f773c153fde34172` |
| `fp.m127.mul_karatsuba.d2048` | 20.682 ms | 2.308 ms | 21.645 ms | 18.955 ms | 2.373 ms | 0.955x | ok `695e0d0f13c14e87` |
| `fp.m127.mul_auto.d2048` | 4.098 ms | 472.6 us | 4.052 ms | 3.768 ms | 293.2 us | 1.011x | ok `1bc667b834b9e559` |
| `fp.m127.divmod.d2000_1000` | 10.855 ms | 1.265 ms | 10.388 ms | 9.687 ms | 783.5 us | 1.045x | ok `a8babdcd119eb404` |
| `fp.m127.gcd.d400` | 29.688 ms | 2.555 ms | 29.494 ms | 27.806 ms | 1.603 ms | 1.007x | ok `1b15fac38c54accc` |
| `fp.m127.gcd_ext.d300` | 31.706 ms | 2.181 ms | 32.044 ms | 29.858 ms | 1.657 ms | 0.989x | ok `45a1b9001453c6ef` |
| `fp.m127.pow_mod.d64_k1` | 69.030 ms | 2.390 ms | 70.116 ms | 66.526 ms | 4.521 ms | 0.985x | ok `9e9e535b3c4e4e7d` |
| `fp.m127.is_irreducible.d16` | 13.979 ms | 1.393 ms | 13.667 ms | 12.313 ms | 734.6 us | 1.023x | ok `2cfb69ee2746aa9d` |
| `fp.m127.factor.d19` | 91.609 ms | 4.677 ms | 91.592 ms | 85.510 ms | 7.074 ms | 1.000x | ok `7b0defab2417737c` |
| `fp.p256.elem_muladd.n200000` | 81.436 ms | 1.935 ms | 80.970 ms | 75.967 ms | 4.329 ms | 1.006x | ok `cc3cab838e9b19f1` |
| `fp.p256.elem_inv.n50000` | 167.642 ms | 5.829 ms | 167.631 ms | 153.947 ms | 5.899 ms | 1.000x | ok `89f920a8eba0c743` |
| `fp.p256.elem_pow.n2000` | 47.979 ms | 1.665 ms | 49.911 ms | 46.338 ms | 3.681 ms | 0.961x | ok `40472b256a3c877d` |
| `fp.p256.mul_school.d64` | 6.137 ms | 938.5 us | 6.128 ms | 5.721 ms | 706.4 us | 1.002x | ok `f5d3d87e3ed0977e` |
| `fp.p256.mul_school.d512` | 15.915 ms | 783.2 us | 15.596 ms | 14.554 ms | 1.137 ms | 1.020x | ok `6fd19bd7b3ceab69` |
| `fp.p256.mul_karatsuba.d512` | 5.481 ms | 1.052 ms | 5.567 ms | 3.869 ms | 1.685 ms | 0.984x | ok `4b2f43bf0f878408` |
| `fp.p256.mul_karatsuba.d2048` | 25.047 ms | 8.127 ms | 28.973 ms | 23.478 ms | 8.897 ms | 0.865x | ok `b4b89a80c2a3aaa9` |
| `fp.p256.mul_fft.d4096` | 36.195 ms | 2.037 ms | 36.321 ms | 34.668 ms | 1.994 ms | 0.997x | ok `60ed510564513dfb` |
| `fp.p256.mul_fft.d8192` | 81.147 ms | 8.203 ms | 84.240 ms | 74.967 ms | 8.809 ms | 0.963x | ok `47812935df4de5fe` |
| `fp.p256.mul_auto.d2048` | 8.769 ms | 1.161 ms | 8.777 ms | 7.887 ms | 826.8 us | 0.999x | ok `3d89402081500506` |
| `fp.p256.divmod.d1000_500` | 8.644 ms | 812.9 us | 8.762 ms | 7.887 ms | 1.053 ms | 0.987x | ok `0bebce51d98eb3f4` |
| `fp.p256.gcd.d300` | 26.808 ms | 936.5 us | 27.017 ms | 24.869 ms | 1.400 ms | 0.992x | ok `f69e4050cb1a5861` |
| `fp.p256.gcd_ext.d200` | 30.828 ms | 1.874 ms | 30.903 ms | 27.977 ms | 1.830 ms | 0.998x | ok `1a9df499c28b9f49` |
| `fp.p256.pow_mod.d64_k1` | 169.075 ms | 9.196 ms | 168.258 ms | 156.385 ms | 8.396 ms | 1.005x | ok `dbe2b595567fbb9c` |
| `fp.p256.is_irreducible.d12` | 17.842 ms | 1.444 ms | 17.534 ms | 16.064 ms | 1.433 ms | 1.018x | ok `7b839a40c98b0ef0` |
| `fp.p256.factor.d12` | 44.359 ms | 1.859 ms | 43.381 ms | 41.073 ms | 2.185 ms | 1.023x | ok `c33e3512ed388a7e` |
| `fq.q3e30.elem_muladd.n4000` | 14.454 ms | 627.6 us | 13.941 ms | 13.157 ms | 642.7 us | 1.037x | ok `16edfcc06a08028a` |
| `fq.q3e30.elem_inv.n100` | 1.757 ms | 105.0 us | 1.780 ms | 1.736 ms | 228.3 us | 0.987x | ok `a8167c01d1cfc5dd` |
| `fq.q3e30.mul_school.d32` | 1.719 ms | 79.0 us | 1.725 ms | 1.703 ms | 195.0 us | 0.997x | ok `ed0ca8eb7d325de3` |
| `fq.q3e30.mul_karatsuba.d128` | 6.007 ms | 1.662 ms | 6.712 ms | 5.195 ms | 832.3 us | 0.895x | ok `2464b901bd242c7d` |
| `fq.q3e30.divmod.d64_32` | 1.769 ms | 79.9 us | 1.780 ms | 1.746 ms | 107.3 us | 0.994x | ok `a2b35f414a139741` |
| `fq.q3e30.factor.d5` | 8.865 ms | 425.2 us | 8.935 ms | 8.620 ms | 414.5 us | 0.992x | ok `417a5b753c7a2a07` |
| `fq.q2e8.elem_muladd.n50000` | 21.712 ms | 1.587 ms | 21.464 ms | 20.510 ms | 1.070 ms | 1.012x | ok `c0ffe16ee54f128c` |
| `fq.q2e8.elem_inv.n5000` | 6.502 ms | 1.483 ms | 6.179 ms | 5.543 ms | 976.4 us | 1.052x | ok `81d71b5687e3e735` |
| `fq.q2e8.mul_school.d128` | 1.454 ms | 333.5 us | 1.433 ms | 1.369 ms | 244.0 us | 1.014x | ok `c056d9d992271db9` |
| `fq.q2e8.mul_karatsuba.d512` | 6.412 ms | 859.5 us | 6.452 ms | 5.457 ms | 929.5 us | 0.994x | ok `de1089fdce26d02c` |
| `fq.q2e8.is_irreducible.d16` | 3.497 ms | 247.2 us | 3.411 ms | 3.302 ms | 197.2 us | 1.025x | ok `9478d590a29a9cbf` |
| `fq.q2e8.factor.d28` | 6.367 ms | 521.4 us | 6.314 ms | 6.006 ms | 411.0 us | 1.009x | ok `d2fc95c10ececedd` |
| `fq.q7e4.elem_muladd.n50000` | 14.374 ms | 641.8 us | 14.813 ms | 13.542 ms | 665.6 us | 0.970x | ok `a6cf18a43d40c5c0` |
| `fq.q7e4.elem_inv.n5000` | 4.634 ms | 841.7 us | 4.606 ms | 4.243 ms | 911.0 us | 1.006x | ok `8664ac6d09409f94` |
| `fq.q7e4.mul_school.d128` | 680.1 us | 161.3 us | 669.1 us | 645.8 us | 70.3 us | 1.016x | ok `71d186560b1d8153` |
| `fq.q7e4.mul_karatsuba.d512` | 4.487 ms | 748.2 us | 4.413 ms | 3.692 ms | 881.6 us | 1.017x | ok `6a2e9b751ab48e99` |
| `fq.q7e4.gcd.d200` | 8.035 ms | 298.5 us | 7.849 ms | 7.645 ms | 463.8 us | 1.024x | ok `8cf19ac90e512feb` |
| `fq.q7e4.factor.d20` | 2.937 ms | 147.1 us | 2.914 ms | 2.853 ms | 148.0 us | 1.008x | ok `ce23318b1b5fd48e` |
| `fq.m61e2.elem_muladd.n100000` | 21.798 ms | 1.568 ms | 21.625 ms | 20.313 ms | 1.943 ms | 1.008x | ok `1814347c960518d8` |
| `fq.m61e2.elem_inv.n1000` | 952.7 us | 81.9 us | 955.5 us | 925.8 us | 88.5 us | 0.997x | ok `45326d8095032432` |
| `fq.m61e2.mul_school.d128` | 283.7 us | 21.5 us | 286.8 us | 265.3 us | 17.5 us | 0.989x | ok `2cb38bd639886e3c` |
| `fq.m61e2.mul_karatsuba.d512` | 3.461 ms | 550.3 us | 3.411 ms | 2.658 ms | 706.4 us | 1.015x | ok `f138ad86c0a7a838` |
| `fq.m61e2.mul_fft.d1024` | 35.192 ms | 1.804 ms | 35.046 ms | 32.803 ms | 1.441 ms | 1.004x | ok `a45642d67f5a9912` |
| `fq.m61e2.mul_fft.d4096` | 80.880 ms | 4.434 ms | 81.648 ms | 76.510 ms | 4.433 ms | 0.991x | ok `51844e8f95a9c0ab` |
| `fq.m61e2.mul_auto.d2048` | 6.290 ms | 482.0 us | 6.332 ms | 5.282 ms | 435.8 us | 0.993x | ok `76e66ad22fb3fa2e` |
| `fq.m61e2.divmod.d500_250` | 1.086 ms | 64.2 us | 1.079 ms | 1.057 ms | 34.6 us | 1.006x | ok `53509a111940ee2a` |
| `fq.m61e2.gcd.d150` | 3.172 ms | 315.3 us | 3.147 ms | 3.060 ms | 258.7 us | 1.008x | ok `bd8c3c98cc56df84` |
| `fq.m61e2.factor.d7` | 7.203 ms | 728.7 us | 6.878 ms | 6.655 ms | 1.078 ms | 1.047x | ok `d61430d11ee187c3` |
| `fq.build.q3e30` | 647.6 us | 100.1 us | 380.5 us | 364.2 us | 58.0 us | 1.702x | ok `4973a3e567e2fc74` |
| `fq.build.q7e20` | 477.1 us | 61.3 us | 342.5 us | 328.4 us | 35.9 us | 1.393x | ok `9c377cd554aa108e` |
| `fast.ntt.dft.k10` | 2.859 ms | 254.6 us | 2.913 ms | 2.829 ms | 388.4 us | 0.981x | ok `b61f454e63ae6d8d` |
| `fast.ntt.dft.k11` | 2.933 ms | 355.3 us | 2.948 ms | 2.854 ms | 461.3 us | 0.995x | ok `311d4fed0f42034c` |
| `fast.ntt.dft.k12` | 3.035 ms | 480.8 us | 3.046 ms | 2.904 ms | 489.7 us | 0.997x | ok `70691e178cfc0c9d` |
| `fast.ntt.idft.k12` | 3.101 ms | 632.7 us | 3.052 ms | 2.936 ms | 391.6 us | 1.016x | ok `e0677b59c4e33615` |
| `fast.ntt.toeplitz_lower_inverse.n1024` | 862.0 us | 103.2 us | 556.2 us | 541.9 us | 24.8 us | 1.550x | ok `7749eb8ff8cac3d0` |
| `fast.ntt.toeplitz_lower_inverse.n8192` | 14.556 ms | 616.7 us | 6.874 ms | 5.999 ms | 581.0 us | 2.118x | ok `0d31554b19d8b996` |
| `fast.ntt.toeplitz_lower_vec.n4096` | 2.002 ms | 153.2 us | 950.0 us | 808.5 us | 114.9 us | 2.108x | ok `86e36ad2a1eacb05` |
| `fast.ntt.toeplitz_divmod.d2000_1000` | 6.619 ms | 485.5 us | 6.501 ms | 6.262 ms | 484.2 us | 1.018x | ok `b138cdcceebcc32e` |
| `fast.p7.toeplitz_lower_inverse.n1024` | 237.8 us | 19.9 us | 238.3 us | 231.3 us | 23.2 us | 0.998x | ok `5df08c92a814155a` |
| `fast.p7.toeplitz_divmod.d2000_1000` | 3.075 ms | 461.0 us | 3.041 ms | 2.967 ms | 359.4 us | 1.011x | ok `433d13716b2f95e8` |
