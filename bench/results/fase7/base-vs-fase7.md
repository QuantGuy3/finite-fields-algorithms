# Comparativa de benchmarks alcp

Base: **f6-desd3**. Speedup = mediana(base) / mediana(variante); > 1 significa más rápido.
IQR = rango intercuartílico de las muestras. Las medias geométricas excluyen escenarios con error.

- **f6-desd3**: modo interleave, 2026-11-09T11:37:52+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -fno-ident -DV3`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-81-g6221fa8 (con cambios sin commit), sha256 01b563a80c52  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: sí
- **f7i534**: modo interleave, 2026-11-09T11:38:10+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -fno-ident -Wl,--no-insert-timestamp`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-93-g0877c55, sha256 1d935c5161da  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: sí

Checksums: todos los escenarios coinciden con la base.

## Media geométrica del speedup

| ámbito | f7i534 |
|---|---|
| **global** | 1.721x (n=128) |
| grupo `fast` | 1.475x (n=10) |
| grupo `fp.m127` | 3.589x (n=14) |
| grupo `fp.m61` | 0.997x (n=14) |
| grupo `fp.ntt` | 1.339x (n=16) |
| grupo `fp.p101` | 1.031x (n=14) |
| grupo `fp.p256` | 0.980x (n=16) |
| grupo `fp.p7` | 1.047x (n=14) |
| grupo `fq.build` | 0.990x (n=2) |
| grupo `fq.m61e2` | 0.996x (n=10) |
| grupo `fq.q2e8` | 45.392x (n=6) |
| grupo `fq.q3e30` | 2.536x (n=6) |
| grupo `fq.q7e4` | 10.185x (n=6) |
| operación `build` | 0.990x (n=2) |
| operación `dft` | 1.000x (n=3) |
| operación `divmod` | 1.175x (n=8) |
| operación `elem_inv` | 3.317x (n=10) |
| operación `elem_muladd` | 2.559x (n=10) |
| operación `elem_pow` | 1.041x (n=6) |
| operación `factor` | 2.088x (n=10) |
| operación `gcd` | 1.756x (n=8) |
| operación `gcd_ext` | 1.726x (n=6) |
| operación `idft` | 0.993x (n=1) |
| operación `is_irreducible` | 2.115x (n=7) |
| operación `mul_auto` | 1.181x (n=7) |
| operación `mul_fft` | 1.597x (n=6) |
| operación `mul_karatsuba` | 1.658x (n=16) |
| operación `mul_school` | 1.726x (n=16) |
| operación `pow_mod` | 1.292x (n=6) |
| operación `toeplitz_divmod` | 1.001x (n=2) |
| operación `toeplitz_lower_inverse` | 2.677x (n=3) |
| operación `toeplitz_lower_vec` | 2.564x (n=1) |

## Detalle por escenario

| escenario | base mediana | base IQR | f7i534 mediana | f7i534 mín | f7i534 IQR | speedup | checksum |
|---|---:|---:|---:|---:|---:|---:|---|
| `fp.p7.elem_muladd.n400000` | 1.704 ms | 30.6 us | 1.705 ms | 1.674 ms | 34.1 us | 0.999x | ok `cf4797a827e7f9ec` |
| `fp.p7.elem_inv.n100000` | 335.8 us | 7.0 us | 336.1 us | 333.9 us | 7.5 us | 0.999x | ok `7beb5a3238b919cb` |
| `fp.p7.elem_pow.n100000` | 940.8 us | 20.5 us | 935.9 us | 927.6 us | 21.1 us | 1.005x | ok `a7f8c429d0d81b73` |
| `fp.p7.mul_school.d64` | 22.0 us | 0.5 us | 21.9 us | 21.4 us | 0.9 us | 1.005x | ok `559357b5bc3c57c2` |
| `fp.p7.mul_school.d512` | 53.3 us | 1.8 us | 53.9 us | 52.4 us | 2.9 us | 0.989x | ok `5936b7480b151eb3` |
| `fp.p7.mul_karatsuba.d512` | 75.8 us | 3.2 us | 77.3 us | 73.8 us | 3.2 us | 0.981x | ok `d72c59b658ac0379` |
| `fp.p7.mul_karatsuba.d4096` | 1.029 ms | 202.7 us | 1.011 ms | 736.7 us | 174.0 us | 1.018x | ok `00f3dde93cbd25a4` |
| `fp.p7.mul_auto.d2048` | 51.0 us | 0.9 us | 50.9 us | 50.3 us | 1.1 us | 1.002x | ok `2b256c747e0070c7` |
| `fp.p7.divmod.d2000_1000` | 87.6 us | 3.2 us | 87.1 us | 84.4 us | 2.3 us | 1.006x | ok `958d383cda71b699` |
| `fp.p7.gcd.d500` | 263.8 us | 8.9 us | 259.1 us | 252.2 us | 8.7 us | 1.018x | ok `cdf18334cf9b13d7` |
| `fp.p7.gcd_ext.d400` | 574.2 us | 15.6 us | 310.2 us | 298.4 us | 14.7 us | 1.851x | ok `984ac156886a3f82` |
| `fp.p7.pow_mod.d128_k48` | 817.4 us | 24.7 us | 804.6 us | 782.9 us | 26.3 us | 1.016x | ok `e3cedb8421cd962f` |
| `fp.p7.is_irreducible.d64` | 377.2 us | 10.2 us | 375.5 us | 368.1 us | 9.8 us | 1.005x | ok `9cb1bec87bc433d0` |
| `fp.p7.factor.d62` | 144.2 us | 4.9 us | 145.7 us | 140.3 us | 5.2 us | 0.990x | ok `3425aecb393956ba` |
| `fp.p101.elem_muladd.n400000` | 1.714 ms | 91.2 us | 1.683 ms | 1.674 ms | 30.3 us | 1.018x | ok `e51aa870ce6269eb` |
| `fp.p101.elem_inv.n100000` | 670.0 us | 14.9 us | 678.6 us | 663.6 us | 14.0 us | 0.987x | ok `69d137d626075364` |
| `fp.p101.elem_pow.n100000` | 2.049 ms | 52.9 us | 1.955 ms | 1.930 ms | 43.7 us | 1.048x | ok `09fec44561c97afa` |
| `fp.p101.mul_school.d64` | 23.7 us | 0.5 us | 23.7 us | 23.3 us | 0.6 us | 1.000x | ok `5a0fa852855ea8b0` |
| `fp.p101.mul_school.d512` | 61.7 us | 1.6 us | 63.0 us | 60.9 us | 2.6 us | 0.979x | ok `8289633c899317a1` |
| `fp.p101.mul_karatsuba.d512` | 77.1 us | 2.2 us | 81.6 us | 77.6 us | 3.0 us | 0.945x | ok `fea453dad54c0188` |
| `fp.p101.mul_karatsuba.d4096` | 1.271 ms | 272.8 us | 1.282 ms | 1.040 ms | 284.6 us | 0.991x | ok `82f83f0b7fe7713b` |
| `fp.p101.mul_auto.d2048` | 75.1 us | 2.1 us | 76.7 us | 74.3 us | 2.5 us | 0.979x | ok `3631fff5520b84d9` |
| `fp.p101.divmod.d2000_1000` | 124.0 us | 5.0 us | 130.6 us | 124.9 us | 12.7 us | 0.949x | ok `ce32beed1dda91cd` |
| `fp.p101.gcd.d500` | 309.0 us | 12.6 us | 312.9 us | 298.8 us | 13.1 us | 0.988x | ok `5cad0c7c5f168edb` |
| `fp.p101.gcd_ext.d400` | 655.9 us | 21.7 us | 354.4 us | 339.4 us | 16.5 us | 1.851x | ok `4fe0534b04a103a1` |
| `fp.p101.pow_mod.d128_k20` | 1.027 ms | 36.5 us | 1.060 ms | 1.023 ms | 32.7 us | 0.969x | ok `bc5334d07b025488` |
| `fp.p101.is_irreducible.d48` | 281.2 us | 10.0 us | 284.7 us | 272.7 us | 44.8 us | 0.988x | ok `5146dc1467eee05b` |
| `fp.p101.factor.d100` | 630.7 us | 27.2 us | 642.3 us | 616.3 us | 22.6 us | 0.982x | ok `48a7604aa4ae63ee` |
| `fp.ntt.elem_muladd.n400000` | 2.463 ms | 62.3 us | 2.537 ms | 2.460 ms | 65.5 us | 0.971x | ok `682b833fe4edd17b` |
| `fp.ntt.elem_inv.n100000` | 2.963 ms | 76.3 us | 3.012 ms | 2.943 ms | 119.7 us | 0.984x | ok `d275b741e6b3050e` |
| `fp.ntt.elem_pow.n50000` | 6.756 ms | 79.2 us | 6.774 ms | 6.671 ms | 134.8 us | 0.997x | ok `c27b88ab2935c700` |
| `fp.ntt.mul_school.d64` | 52.4 us | 1.5 us | 53.9 us | 52.1 us | 2.1 us | 0.972x | ok `67c7e1645a9cc6ae` |
| `fp.ntt.mul_school.d512` | 143.7 us | 4.0 us | 147.2 us | 142.7 us | 4.0 us | 0.976x | ok `e296393a6a637046` |
| `fp.ntt.mul_karatsuba.d512` | 99.6 us | 6.9 us | 99.4 us | 95.6 us | 3.3 us | 1.002x | ok `964e15063fb16df4` |
| `fp.ntt.mul_karatsuba.d4096` | 1.509 ms | 294.4 us | 1.476 ms | 1.201 ms | 258.3 us | 1.022x | ok `d656b678f0c2ca16` |
| `fp.ntt.mul_fft.d4096` | 730.6 us | 39.5 us | 193.2 us | 164.3 us | 33.6 us | 3.781x | ok `f4c81b85220fb683` |
| `fp.ntt.mul_fft.d32768` | 7.279 ms | 482.5 us | 1.222 ms | 1.010 ms | 138.1 us | 5.955x | ok `8228e76950faded3` |
| `fp.ntt.mul_auto.d2048` | 384.6 us | 26.2 us | 145.7 us | 127.6 us | 16.0 us | 2.640x | ok `02e67c888996018d` |
| `fp.ntt.divmod.d2000_1000` | 432.2 us | 21.9 us | 433.5 us | 421.6 us | 17.5 us | 0.997x | ok `7adfbb2614325099` |
| `fp.ntt.gcd.d500` | 327.0 us | 10.4 us | 329.5 us | 315.6 us | 15.4 us | 0.992x | ok `678962135e7a2447` |
| `fp.ntt.gcd_ext.d400` | 664.7 us | 30.4 us | 409.1 us | 387.3 us | 24.3 us | 1.625x | ok `f9dc4e400a7954ee` |
| `fp.ntt.pow_mod.d128_k5` | 3.662 ms | 123.8 us | 2.978 ms | 2.867 ms | 134.8 us | 1.230x | ok `0eb2514aceb37261` |
| `fp.ntt.is_irreducible.d32` | 207.7 us | 7.5 us | 208.4 us | 201.9 us | 6.8 us | 0.997x | ok `1f2827a0d8aeb65e` |
| `fp.ntt.factor.d40` | 657.8 us | 30.4 us | 667.2 us | 636.8 us | 23.9 us | 0.986x | ok `297c76c423fe55aa` |
| `fp.m61.elem_muladd.n400000` | 2.178 ms | 39.8 us | 2.233 ms | 2.167 ms | 58.1 us | 0.975x | ok `b8e325f73ad90582` |
| `fp.m61.elem_inv.n100000` | 8.336 ms | 88.5 us | 8.418 ms | 8.220 ms | 184.2 us | 0.990x | ok `41e7fe38de6fdcc1` |
| `fp.m61.elem_pow.n20000` | 6.050 ms | 89.8 us | 6.150 ms | 5.992 ms | 152.5 us | 0.984x | ok `a3680ece7fc07932` |
| `fp.m61.mul_school.d64` | 77.8 us | 2.1 us | 79.5 us | 76.7 us | 2.8 us | 0.979x | ok `bed6e9e86907de79` |
| `fp.m61.mul_school.d512` | 201.5 us | 5.0 us | 207.1 us | 200.1 us | 3.0 us | 0.973x | ok `43c3884710f7d468` |
| `fp.m61.mul_karatsuba.d512` | 128.7 us | 6.5 us | 128.2 us | 123.4 us | 5.1 us | 1.004x | ok `a2dcbeaa6574fcb8` |
| `fp.m61.mul_karatsuba.d4096` | 1.743 ms | 272.1 us | 1.765 ms | 1.465 ms | 229.3 us | 0.987x | ok `37d5deaa1cfeb600` |
| `fp.m61.mul_auto.d2048` | 584.3 us | 126.0 us | 549.6 us | 527.4 us | 63.9 us | 1.063x | ok `73174ffa4e90d6fa` |
| `fp.m61.divmod.d2000_1000` | 877.0 us | 77.2 us | 852.9 us | 825.8 us | 30.7 us | 1.028x | ok `a8710fc5aaa24d09` |
| `fp.m61.gcd.d500` | 349.6 us | 12.9 us | 349.2 us | 335.1 us | 16.4 us | 1.001x | ok `fda394dde025c1a2` |
| `fp.m61.gcd_ext.d400` | 679.6 us | 23.4 us | 659.6 us | 629.3 us | 21.9 us | 1.030x | ok `be3d81e8acc44b48` |
| `fp.m61.pow_mod.d64_k2` | 1.486 ms | 62.7 us | 1.512 ms | 1.457 ms | 54.0 us | 0.982x | ok `394fb320787700fe` |
| `fp.m61.is_irreducible.d24` | 267.1 us | 9.3 us | 273.1 us | 264.9 us | 9.0 us | 0.978x | ok `9047cf1d3fc29abf` |
| `fp.m61.factor.d26` | 864.2 us | 38.2 us | 879.6 us | 840.4 us | 53.7 us | 0.982x | ok `931b8d0b8b91b0b4` |
| `fp.m127.elem_muladd.n200000` | 23.118 ms | 555.6 us | 3.074 ms | 2.982 ms | 97.8 us | 7.521x | ok `16fab7b139aa9aab` |
| `fp.m127.elem_inv.n50000` | 29.725 ms | 1.801 ms | 28.367 ms | 27.015 ms | 1.979 ms | 1.048x | ok `c7d67d6018036bf1` |
| `fp.m127.elem_pow.n10000` | 22.892 ms | 638.6 us | 18.585 ms | 18.150 ms | 403.3 us | 1.232x | ok `308896e55b40c815` |
| `fp.m127.mul_school.d64` | 1.615 ms | 57.4 us | 284.9 us | 275.5 us | 9.8 us | 5.667x | ok `a9664408b0fb88e7` |
| `fp.m127.mul_school.d512` | 3.779 ms | 178.4 us | 790.9 us | 773.0 us | 22.5 us | 4.779x | ok `98eca357802203c0` |
| `fp.m127.mul_karatsuba.d512` | 1.901 ms | 234.8 us | 400.3 us | 380.5 us | 25.6 us | 4.748x | ok `f773c153fde34172` |
| `fp.m127.mul_karatsuba.d2048` | 10.435 ms | 359.7 us | 1.577 ms | 1.258 ms | 121.9 us | 6.617x | ok `695e0d0f13c14e87` |
| `fp.m127.mul_auto.d2048` | 1.688 ms | 115.3 us | 1.452 ms | 1.351 ms | 87.0 us | 1.163x | ok `1bc667b834b9e559` |
| `fp.m127.divmod.d2000_1000` | 4.064 ms | 321.0 us | 3.035 ms | 2.941 ms | 96.5 us | 1.339x | ok `a8babdcd119eb404` |
| `fp.m127.gcd.d400` | 11.836 ms | 274.2 us | 1.686 ms | 1.633 ms | 63.8 us | 7.022x | ok `1b15fac38c54accc` |
| `fp.m127.gcd_ext.d300` | 12.212 ms | 254.1 us | 2.697 ms | 2.609 ms | 110.5 us | 4.527x | ok `45a1b9001453c6ef` |
| `fp.m127.pow_mod.d64_k1` | 26.140 ms | 542.4 us | 6.794 ms | 6.644 ms | 207.7 us | 3.848x | ok `9e9e535b3c4e4e7d` |
| `fp.m127.is_irreducible.d16` | 5.108 ms | 139.9 us | 827.7 us | 800.5 us | 45.9 us | 6.172x | ok `2cfb69ee2746aa9d` |
| `fp.m127.factor.d19` | 30.555 ms | 246.8 us | 5.050 ms | 4.902 ms | 200.7 us | 6.050x | ok `7b0defab2417737c` |
| `fp.p256.elem_muladd.n200000` | 31.430 ms | 757.2 us | 31.774 ms | 31.407 ms | 367.6 us | 0.989x | ok `cc3cab838e9b19f1` |
| `fp.p256.elem_inv.n50000` | 64.193 ms | 609.5 us | 64.697 ms | 63.789 ms | 2.135 ms | 0.992x | ok `89f920a8eba0c743` |
| `fp.p256.elem_pow.n2000` | 18.832 ms | 268.5 us | 18.828 ms | 18.657 ms | 762.4 us | 1.000x | ok `40472b256a3c877d` |
| `fp.p256.mul_school.d64` | 2.409 ms | 66.5 us | 2.445 ms | 2.378 ms | 138.0 us | 0.985x | ok `f5d3d87e3ed0977e` |
| `fp.p256.mul_school.d512` | 6.177 ms | 247.0 us | 6.221 ms | 6.021 ms | 198.4 us | 0.993x | ok `6fd19bd7b3ceab69` |
| `fp.p256.mul_karatsuba.d512` | 2.040 ms | 197.6 us | 2.104 ms | 1.728 ms | 220.3 us | 0.970x | ok `4b2f43bf0f878408` |
| `fp.p256.mul_karatsuba.d2048` | 12.392 ms | 436.0 us | 12.475 ms | 12.106 ms | 306.1 us | 0.993x | ok `b4b89a80c2a3aaa9` |
| `fp.p256.mul_fft.d4096` | 16.068 ms | 3.374 ms | 18.206 ms | 14.577 ms | 3.564 ms | 0.883x | ok `60ed510564513dfb` |
| `fp.p256.mul_fft.d8192` | 33.367 ms | 6.371 ms | 39.152 ms | 31.723 ms | 5.954 ms | 0.852x | ok `47812935df4de5fe` |
| `fp.p256.mul_auto.d2048` | 3.466 ms | 224.3 us | 3.501 ms | 3.322 ms | 295.6 us | 0.990x | ok `3d89402081500506` |
| `fp.p256.divmod.d1000_500` | 3.281 ms | 311.4 us | 3.209 ms | 3.060 ms | 275.4 us | 1.022x | ok `0bebce51d98eb3f4` |
| `fp.p256.gcd.d300` | 10.394 ms | 383.5 us | 10.213 ms | 10.052 ms | 277.4 us | 1.018x | ok `f69e4050cb1a5861` |
| `fp.p256.gcd_ext.d200` | 11.720 ms | 599.6 us | 11.515 ms | 11.335 ms | 155.9 us | 1.018x | ok `1a9df499c28b9f49` |
| `fp.p256.pow_mod.d64_k1` | 65.413 ms | 1.381 ms | 64.242 ms | 63.376 ms | 1.482 ms | 1.018x | ok `dbe2b595567fbb9c` |
| `fp.p256.is_irreducible.d12` | 6.789 ms | 177.7 us | 6.801 ms | 6.527 ms | 198.0 us | 0.998x | ok `7b839a40c98b0ef0` |
| `fp.p256.factor.d12` | 14.331 ms | 468.6 us | 14.632 ms | 14.028 ms | 620.0 us | 0.979x | ok `c33e3512ed388a7e` |
| `fq.q3e30.elem_muladd.n4000` | 3.886 ms | 178.0 us | 2.075 ms | 2.036 ms | 64.5 us | 1.873x | ok `16edfcc06a08028a` |
| `fq.q3e30.elem_inv.n100` | 721.8 us | 29.2 us | 147.9 us | 144.5 us | 5.8 us | 4.880x | ok `a8167c01d1cfc5dd` |
| `fq.q3e30.mul_school.d32` | 439.4 us | 17.2 us | 154.4 us | 149.5 us | 5.5 us | 2.846x | ok `ed0ca8eb7d325de3` |
| `fq.q3e30.mul_karatsuba.d128` | 1.760 ms | 251.1 us | 1.019 ms | 782.0 us | 261.7 us | 1.726x | ok `2464b901bd242c7d` |
| `fq.q3e30.divmod.d64_32` | 444.6 us | 13.5 us | 165.2 us | 163.3 us | 4.4 us | 2.691x | ok `a2b35f414a139741` |
| `fq.q3e30.factor.d5` | 2.530 ms | 91.9 us | 1.150 ms | 1.137 ms | 29.5 us | 2.200x | ok `417a5b753c7a2a07` |
| `fq.q2e8.elem_muladd.n50000` | 7.726 ms | 298.3 us | 147.5 us | 143.2 us | 4.7 us | 52.380x | ok `c0ffe16ee54f128c` |
| `fq.q2e8.elem_inv.n5000` | 2.305 ms | 73.3 us | 5.7 us | 5.6 us | 0.2 us | 404.316x | ok `81d71b5687e3e735` |
| `fq.q2e8.mul_school.d128` | 431.4 us | 14.2 us | 12.5 us | 12.1 us | 0.4 us | 34.512x | ok `c056d9d992271db9` |
| `fq.q2e8.mul_karatsuba.d512` | 2.142 ms | 368.3 us | 93.7 us | 89.9 us | 3.2 us | 22.864x | ok `de1089fdce26d02c` |
| `fq.q2e8.is_irreducible.d16` | 1.235 ms | 56.0 us | 38.9 us | 37.3 us | 1.8 us | 31.751x | ok `9478d590a29a9cbf` |
| `fq.q2e8.factor.d28` | 2.102 ms | 79.0 us | 127.5 us | 122.6 us | 4.7 us | 16.485x | ok `d2fc95c10ececedd` |
| `fq.q7e4.elem_muladd.n50000` | 5.505 ms | 117.9 us | 314.1 us | 308.5 us | 10.7 us | 17.526x | ok `a6cf18a43d40c5c0` |
| `fq.q7e4.elem_inv.n5000` | 1.786 ms | 73.6 us | 21.6 us | 21.3 us | 0.6 us | 82.667x | ok `8664ac6d09409f94` |
| `fq.q7e4.mul_school.d128` | 203.2 us | 7.6 us | 75.2 us | 62.0 us | 5.1 us | 2.702x | ok `71d186560b1d8153` |
| `fq.q7e4.mul_karatsuba.d512` | 1.545 ms | 223.0 us | 534.5 us | 518.1 us | 18.4 us | 2.891x | ok `6a2e9b751ab48e99` |
| `fq.q7e4.gcd.d200` | 3.012 ms | 152.0 us | 235.1 us | 226.7 us | 8.4 us | 12.810x | ok `8cf19ac90e512feb` |
| `fq.q7e4.factor.d20` | 986.9 us | 40.1 us | 128.2 us | 121.0 us | 10.9 us | 7.699x | ok `ce23318b1b5fd48e` |
| `fq.m61e2.elem_muladd.n100000` | 9.063 ms | 187.4 us | 9.293 ms | 8.716 ms | 460.5 us | 0.975x | ok `1814347c960518d8` |
| `fq.m61e2.elem_inv.n1000` | 392.2 us | 12.0 us | 396.8 us | 381.7 us | 13.8 us | 0.988x | ok `45326d8095032432` |
| `fq.m61e2.mul_school.d128` | 121.9 us | 4.2 us | 121.4 us | 116.3 us | 6.0 us | 1.004x | ok `2cb38bd639886e3c` |
| `fq.m61e2.mul_karatsuba.d512` | 1.427 ms | 403.3 us | 1.428 ms | 1.076 ms | 222.5 us | 0.999x | ok `f138ad86c0a7a838` |
| `fq.m61e2.mul_fft.d1024` | 14.612 ms | 1.168 ms | 14.900 ms | 13.468 ms | 1.151 ms | 0.981x | ok `a45642d67f5a9912` |
| `fq.m61e2.mul_fft.d4096` | 40.347 ms | 2.196 ms | 40.409 ms | 39.635 ms | 1.585 ms | 0.998x | ok `51844e8f95a9c0ab` |
| `fq.m61e2.mul_auto.d2048` | 2.619 ms | 180.1 us | 2.591 ms | 2.192 ms | 146.6 us | 1.011x | ok `76e66ad22fb3fa2e` |
| `fq.m61e2.divmod.d500_250` | 458.1 us | 16.1 us | 456.1 us | 443.1 us | 17.9 us | 1.004x | ok `53509a111940ee2a` |
| `fq.m61e2.gcd.d150` | 1.357 ms | 50.0 us | 1.370 ms | 1.341 ms | 36.7 us | 0.991x | ok `bd8c3c98cc56df84` |
| `fq.m61e2.factor.d7` | 2.628 ms | 219.9 us | 2.597 ms | 2.493 ms | 227.2 us | 1.012x | ok `d61430d11ee187c3` |
| `fq.build.q3e30` | 111.6 us | 4.2 us | 112.4 us | 108.5 us | 3.6 us | 0.993x | ok `4973a3e567e2fc74` |
| `fq.build.q7e20` | 106.5 us | 4.1 us | 107.9 us | 104.7 us | 2.7 us | 0.987x | ok `9c377cd554aa108e` |
| `fast.ntt.dft.k10` | 1.238 ms | 235.8 us | 1.211 ms | 1.179 ms | 23.6 us | 1.022x | ok `b61f454e63ae6d8d` |
| `fast.ntt.dft.k11` | 1.229 ms | 43.3 us | 1.228 ms | 1.206 ms | 38.8 us | 1.001x | ok `311d4fed0f42034c` |
| `fast.ntt.dft.k12` | 1.227 ms | 36.8 us | 1.257 ms | 1.222 ms | 64.2 us | 0.977x | ok `70691e178cfc0c9d` |
| `fast.ntt.idft.k12` | 1.236 ms | 49.4 us | 1.245 ms | 1.219 ms | 39.1 us | 0.993x | ok `e0677b59c4e33615` |
| `fast.ntt.toeplitz_lower_inverse.n1024` | 229.4 us | 8.2 us | 42.0 us | 39.5 us | 2.1 us | 5.462x | ok `7749eb8ff8cac3d0` |
| `fast.ntt.toeplitz_lower_inverse.n8192` | 2.863 ms | 220.2 us | 810.9 us | 738.4 us | 88.0 us | 3.531x | ok `0d31554b19d8b996` |
| `fast.ntt.toeplitz_lower_vec.n4096` | 397.7 us | 29.4 us | 155.1 us | 140.6 us | 13.3 us | 2.564x | ok `86e36ad2a1eacb05` |
| `fast.ntt.toeplitz_divmod.d2000_1000` | 2.685 ms | 103.4 us | 2.655 ms | 2.588 ms | 86.4 us | 1.011x | ok `b138cdcceebcc32e` |
| `fast.p7.toeplitz_lower_inverse.n1024` | 56.7 us | 2.5 us | 57.0 us | 54.2 us | 4.1 us | 0.995x | ok `5df08c92a814155a` |
| `fast.p7.toeplitz_divmod.d2000_1000` | 1.168 ms | 34.8 us | 1.180 ms | 1.155 ms | 40.9 us | 0.990x | ok `433d13716b2f95e8` |
