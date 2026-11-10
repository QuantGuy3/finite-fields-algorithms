# Comparativa de benchmarks alcp

Base: **f7fin**. Speedup = mediana(base) / mediana(variante); > 1 significa más rápido.
IQR = rango intercuartílico de las muestras. Las medias geométricas excluyen escenarios con error.

- **f7fin**: modo interleave, 2026-11-10T17:21:13+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -fno-ident -Wl,--build-id=sha1`  
  include `~\AppData\Local\Temp\claude\204af1c0d38157dc\2e18e96a-361a-47be-ad99-f0341eb0cc55\scratchpad\base_f7fin\include` @ ?, sha256 cabe17191f10  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: sí
- **f8fin**: modo interleave, 2026-11-10T17:21:32+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -fno-ident -Wl,--build-id=sha1`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-113-gba76a78, sha256 b943c18962ce  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: sí

Checksums: todos los escenarios coinciden con la base.

## Media geométrica del speedup

| ámbito | f8fin |
|---|---|
| **global** | 1.076x (n=128) |
| grupo `fast` | 0.986x (n=10) |
| grupo `fp.m127` | 0.995x (n=14) |
| grupo `fp.m61` | 1.010x (n=14) |
| grupo `fp.ntt` | 1.083x (n=16) |
| grupo `fp.p101` | 1.273x (n=14) |
| grupo `fp.p256` | 1.009x (n=16) |
| grupo `fp.p7` | 1.281x (n=14) |
| grupo `fq.build` | 1.230x (n=2) |
| grupo `fq.m61e2` | 1.014x (n=10) |
| grupo `fq.q2e8` | 0.982x (n=6) |
| grupo `fq.q3e30` | 1.050x (n=6) |
| grupo `fq.q7e4` | 1.071x (n=6) |
| operación `build` | 1.230x (n=2) |
| operación `dft` | 0.987x (n=3) |
| operación `divmod` | 1.010x (n=8) |
| operación `elem_inv` | 0.998x (n=10) |
| operación `elem_muladd` | 1.020x (n=10) |
| operación `elem_pow` | 0.999x (n=6) |
| operación `factor` | 1.086x (n=10) |
| operación `gcd` | 1.414x (n=8) |
| operación `gcd_ext` | 1.354x (n=6) |
| operación `idft` | 1.002x (n=1) |
| operación `is_irreducible` | 1.187x (n=7) |
| operación `mul_auto` | 1.022x (n=7) |
| operación `mul_fft` | 1.013x (n=6) |
| operación `mul_karatsuba` | 1.020x (n=16) |
| operación `mul_school` | 0.997x (n=16) |
| operación `pow_mod` | 1.333x (n=6) |
| operación `toeplitz_divmod` | 1.014x (n=2) |
| operación `toeplitz_lower_inverse` | 0.960x (n=3) |
| operación `toeplitz_lower_vec` | 0.993x (n=1) |

## Detalle por escenario

| escenario | base mediana | base IQR | f8fin mediana | f8fin mín | f8fin IQR | speedup | checksum |
|---|---:|---:|---:|---:|---:|---:|---|
| `fp.p7.elem_muladd.n400000` | 1.724 ms | 29.4 us | 1.722 ms | 1.674 ms | 41.8 us | 1.001x | ok `cf4797a827e7f9ec` |
| `fp.p7.elem_inv.n100000` | 343.4 us | 6.5 us | 343.5 us | 333.8 us | 8.0 us | 1.000x | ok `7beb5a3238b919cb` |
| `fp.p7.elem_pow.n100000` | 953.4 us | 20.4 us | 951.7 us | 925.2 us | 28.9 us | 1.002x | ok `a7f8c429d0d81b73` |
| `fp.p7.mul_school.d64` | 22.2 us | 0.6 us | 22.2 us | 21.4 us | 0.5 us | 1.000x | ok `559357b5bc3c57c2` |
| `fp.p7.mul_school.d512` | 54.1 us | 1.6 us | 54.2 us | 52.4 us | 1.6 us | 0.998x | ok `5936b7480b151eb3` |
| `fp.p7.mul_karatsuba.d512` | 76.6 us | 3.8 us | 78.3 us | 73.9 us | 3.7 us | 0.978x | ok `d72c59b658ac0379` |
| `fp.p7.mul_karatsuba.d4096` | 1.028 ms | 240.3 us | 1.024 ms | 803.9 us | 206.9 us | 1.004x | ok `00f3dde93cbd25a4` |
| `fp.p7.mul_auto.d2048` | 51.7 us | 1.4 us | 52.0 us | 50.4 us | 1.1 us | 0.994x | ok `2b256c747e0070c7` |
| `fp.p7.divmod.d2000_1000` | 88.6 us | 3.2 us | 89.2 us | 85.1 us | 3.2 us | 0.993x | ok `958d383cda71b699` |
| `fp.p7.gcd.d500` | 173.1 us | 2.6 us | 59.1 us | 55.6 us | 2.2 us | 2.929x | ok `cdf18334cf9b13d7` |
| `fp.p7.gcd_ext.d400` | 279.8 us | 13.4 us | 144.1 us | 128.6 us | 17.1 us | 1.942x | ok `984ac156886a3f82` |
| `fp.p7.pow_mod.d128_k48` | 827.5 us | 30.3 us | 359.0 us | 333.8 us | 18.7 us | 2.305x | ok `e3cedb8421cd962f` |
| `fp.p7.is_irreducible.d64` | 297.8 us | 6.3 us | 162.6 us | 155.0 us | 6.3 us | 1.831x | ok `9cb1bec87bc433d0` |
| `fp.p7.factor.d62` | 138.1 us | 4.4 us | 100.0 us | 95.6 us | 4.9 us | 1.381x | ok `3425aecb393956ba` |
| `fp.p101.elem_muladd.n400000` | 1.714 ms | 34.2 us | 1.716 ms | 1.676 ms | 38.6 us | 0.999x | ok `e51aa870ce6269eb` |
| `fp.p101.elem_inv.n100000` | 678.5 us | 9.0 us | 678.8 us | 663.6 us | 15.6 us | 1.000x | ok `69d137d626075364` |
| `fp.p101.elem_pow.n100000` | 1.976 ms | 35.9 us | 1.981 ms | 1.932 ms | 46.4 us | 0.998x | ok `09fec44561c97afa` |
| `fp.p101.mul_school.d64` | 24.1 us | 0.6 us | 24.1 us | 23.2 us | 0.6 us | 1.000x | ok `5a0fa852855ea8b0` |
| `fp.p101.mul_school.d512` | 62.6 us | 1.5 us | 63.0 us | 61.2 us | 1.6 us | 0.994x | ok `8289633c899317a1` |
| `fp.p101.mul_karatsuba.d512` | 70.9 us | 7.3 us | 75.8 us | 67.3 us | 10.4 us | 0.935x | ok `fea453dad54c0188` |
| `fp.p101.mul_karatsuba.d4096` | 982.5 us | 176.1 us | 988.4 us | 776.7 us | 233.2 us | 0.994x | ok `82f83f0b7fe7713b` |
| `fp.p101.mul_auto.d2048` | 76.5 us | 1.8 us | 76.6 us | 74.3 us | 1.9 us | 0.999x | ok `3631fff5520b84d9` |
| `fp.p101.divmod.d2000_1000` | 129.5 us | 2.1 us | 129.3 us | 124.6 us | 5.2 us | 1.002x | ok `ce32beed1dda91cd` |
| `fp.p101.gcd.d500` | 168.5 us | 3.6 us | 66.5 us | 61.6 us | 2.2 us | 2.534x | ok `5cad0c7c5f168edb` |
| `fp.p101.gcd_ext.d400` | 293.1 us | 9.8 us | 154.2 us | 144.3 us | 10.1 us | 1.900x | ok `4fe0534b04a103a1` |
| `fp.p101.pow_mod.d128_k20` | 1.059 ms | 42.8 us | 430.1 us | 407.7 us | 21.4 us | 2.463x | ok `bc5334d07b025488` |
| `fp.p101.is_irreducible.d48` | 218.5 us | 6.1 us | 122.6 us | 117.2 us | 11.8 us | 1.782x | ok `5146dc1467eee05b` |
| `fp.p101.factor.d100` | 616.5 us | 22.1 us | 407.5 us | 380.0 us | 44.3 us | 1.513x | ok `48a7604aa4ae63ee` |
| `fp.ntt.elem_muladd.n400000` | 2.527 ms | 43.1 us | 2.538 ms | 2.470 ms | 488.1 us | 0.995x | ok `682b833fe4edd17b` |
| `fp.ntt.elem_inv.n100000` | 3.046 ms | 119.6 us | 3.041 ms | 2.965 ms | 173.2 us | 1.002x | ok `d275b741e6b3050e` |
| `fp.ntt.elem_pow.n50000` | 6.766 ms | 61.4 us | 6.829 ms | 6.688 ms | 138.8 us | 0.991x | ok `c27b88ab2935c700` |
| `fp.ntt.mul_school.d64` | 53.7 us | 1.4 us | 53.8 us | 52.0 us | 1.4 us | 0.998x | ok `67c7e1645a9cc6ae` |
| `fp.ntt.mul_school.d512` | 146.2 us | 4.2 us | 147.1 us | 142.8 us | 3.8 us | 0.994x | ok `e296393a6a637046` |
| `fp.ntt.mul_karatsuba.d512` | 97.9 us | 7.5 us | 100.2 us | 95.6 us | 8.5 us | 0.977x | ok `964e15063fb16df4` |
| `fp.ntt.mul_karatsuba.d4096` | 1.177 ms | 245.1 us | 1.204 ms | 903.8 us | 244.7 us | 0.978x | ok `d656b678f0c2ca16` |
| `fp.ntt.mul_fft.d4096` | 182.0 us | 25.8 us | 179.8 us | 154.1 us | 23.2 us | 1.012x | ok `f4c81b85220fb683` |
| `fp.ntt.mul_fft.d32768` | 1.157 ms | 137.9 us | 1.145 ms | 950.9 us | 169.1 us | 1.010x | ok `8228e76950faded3` |
| `fp.ntt.mul_auto.d2048` | 137.5 us | 15.2 us | 139.9 us | 115.4 us | 14.5 us | 0.983x | ok `02e67c888996018d` |
| `fp.ntt.divmod.d2000_1000` | 187.1 us | 5.3 us | 186.1 us | 180.3 us | 6.0 us | 1.005x | ok `7adfbb2614325099` |
| `fp.ntt.gcd.d500` | 182.4 us | 5.0 us | 82.2 us | 78.8 us | 2.6 us | 2.219x | ok `678962135e7a2447` |
| `fp.ntt.gcd_ext.d400` | 339.9 us | 14.0 us | 203.0 us | 192.9 us | 11.3 us | 1.674x | ok `f9dc4e400a7954ee` |
| `fp.ntt.pow_mod.d128_k5` | 1.370 ms | 47.6 us | 1.400 ms | 1.362 ms | 37.6 us | 0.979x | ok `0eb2514aceb37261` |
| `fp.ntt.is_irreducible.d32` | 181.5 us | 5.4 us | 174.8 us | 169.5 us | 5.2 us | 1.038x | ok `1f2827a0d8aeb65e` |
| `fp.ntt.factor.d40` | 665.5 us | 27.0 us | 659.2 us | 623.4 us | 25.9 us | 1.010x | ok `297c76c423fe55aa` |
| `fp.m61.elem_muladd.n400000` | 2.234 ms | 18.0 us | 2.216 ms | 2.169 ms | 53.5 us | 1.008x | ok `b8e325f73ad90582` |
| `fp.m61.elem_inv.n100000` | 8.415 ms | 45.4 us | 8.356 ms | 8.226 ms | 64.6 us | 1.007x | ok `41e7fe38de6fdcc1` |
| `fp.m61.elem_pow.n20000` | 6.112 ms | 34.5 us | 6.096 ms | 5.997 ms | 66.9 us | 1.003x | ok `a3680ece7fc07932` |
| `fp.m61.mul_school.d64` | 79.4 us | 1.9 us | 79.2 us | 76.6 us | 2.3 us | 1.003x | ok `bed6e9e86907de79` |
| `fp.m61.mul_school.d512` | 206.1 us | 5.7 us | 206.9 us | 200.2 us | 6.2 us | 0.996x | ok `43c3884710f7d468` |
| `fp.m61.mul_karatsuba.d512` | 126.8 us | 3.7 us | 130.5 us | 124.0 us | 9.7 us | 0.972x | ok `a2dcbeaa6574fcb8` |
| `fp.m61.mul_karatsuba.d4096` | 1.418 ms | 194.2 us | 1.387 ms | 1.106 ms | 215.6 us | 1.022x | ok `37d5deaa1cfeb600` |
| `fp.m61.mul_auto.d2048` | 657.1 us | 143.9 us | 544.2 us | 524.5 us | 24.3 us | 1.207x | ok `73174ffa4e90d6fa` |
| `fp.m61.divmod.d2000_1000` | 852.8 us | 20.2 us | 859.1 us | 829.0 us | 24.0 us | 0.993x | ok `a8710fc5aaa24d09` |
| `fp.m61.gcd.d500` | 348.4 us | 15.3 us | 358.1 us | 346.2 us | 12.2 us | 0.973x | ok `fda394dde025c1a2` |
| `fp.m61.gcd_ext.d400` | 686.2 us | 15.4 us | 684.8 us | 665.8 us | 24.5 us | 1.002x | ok `be3d81e8acc44b48` |
| `fp.m61.pow_mod.d64_k2` | 1.506 ms | 33.3 us | 1.515 ms | 1.465 ms | 52.0 us | 0.994x | ok `394fb320787700fe` |
| `fp.m61.is_irreducible.d24` | 274.2 us | 10.2 us | 278.8 us | 268.4 us | 16.7 us | 0.983x | ok `9047cf1d3fc29abf` |
| `fp.m61.factor.d26` | 864.9 us | 30.4 us | 865.5 us | 840.0 us | 30.8 us | 0.999x | ok `931b8d0b8b91b0b4` |
| `fp.m127.elem_muladd.n200000` | 3.054 ms | 51.2 us | 3.067 ms | 2.998 ms | 56.6 us | 0.996x | ok `16fab7b139aa9aab` |
| `fp.m127.elem_inv.n50000` | 28.039 ms | 1.089 ms | 26.735 ms | 26.039 ms | 1.744 ms | 1.049x | ok `c7d67d6018036bf1` |
| `fp.m127.elem_pow.n10000` | 18.270 ms | 206.0 us | 18.338 ms | 18.183 ms | 157.1 us | 0.996x | ok `308896e55b40c815` |
| `fp.m127.mul_school.d64` | 283.5 us | 3.6 us | 288.0 us | 275.4 us | 8.7 us | 0.984x | ok `a9664408b0fb88e7` |
| `fp.m127.mul_school.d512` | 795.1 us | 20.1 us | 798.9 us | 773.2 us | 24.2 us | 0.995x | ok `98eca357802203c0` |
| `fp.m127.mul_karatsuba.d512` | 393.3 us | 17.1 us | 398.8 us | 381.5 us | 14.0 us | 0.986x | ok `f773c153fde34172` |
| `fp.m127.mul_karatsuba.d2048` | 1.300 ms | 230.4 us | 1.331 ms | 1.007 ms | 198.3 us | 0.977x | ok `695e0d0f13c14e87` |
| `fp.m127.mul_auto.d2048` | 1.379 ms | 70.6 us | 1.438 ms | 1.373 ms | 79.2 us | 0.959x | ok `1bc667b834b9e559` |
| `fp.m127.divmod.d2000_1000` | 3.005 ms | 63.2 us | 3.027 ms | 2.950 ms | 64.3 us | 0.993x | ok `a8babdcd119eb404` |
| `fp.m127.gcd.d400` | 1.672 ms | 42.5 us | 1.683 ms | 1.631 ms | 47.6 us | 0.994x | ok `1b15fac38c54accc` |
| `fp.m127.gcd_ext.d300` | 2.703 ms | 107.0 us | 2.715 ms | 2.625 ms | 53.1 us | 0.996x | ok `45a1b9001453c6ef` |
| `fp.m127.pow_mod.d64_k1` | 6.867 ms | 129.1 us | 6.851 ms | 6.672 ms | 93.1 us | 1.002x | ok `9e9e535b3c4e4e7d` |
| `fp.m127.is_irreducible.d16` | 834.2 us | 157.7 us | 830.3 us | 799.0 us | 22.5 us | 1.005x | ok `2cfb69ee2746aa9d` |
| `fp.m127.factor.d19` | 5.134 ms | 426.7 us | 5.137 ms | 4.956 ms | 155.5 us | 0.999x | ok `7b0defab2417737c` |
| `fp.p256.elem_muladd.n200000` | 31.756 ms | 568.1 us | 31.117 ms | 30.661 ms | 825.8 us | 1.021x | ok `cc3cab838e9b19f1` |
| `fp.p256.elem_inv.n50000` | 65.304 ms | 1.485 ms | 64.525 ms | 64.082 ms | 441.7 us | 1.012x | ok `89f920a8eba0c743` |
| `fp.p256.elem_pow.n2000` | 19.233 ms | 1.156 ms | 19.133 ms | 18.849 ms | 149.3 us | 1.005x | ok `40472b256a3c877d` |
| `fp.p256.mul_school.d64` | 2.458 ms | 220.7 us | 2.481 ms | 2.410 ms | 32.2 us | 0.990x | ok `f5d3d87e3ed0977e` |
| `fp.p256.mul_school.d512` | 6.249 ms | 369.5 us | 6.255 ms | 6.090 ms | 130.6 us | 0.999x | ok `6fd19bd7b3ceab69` |
| `fp.p256.mul_karatsuba.d512` | 1.739 ms | 282.0 us | 1.686 ms | 1.391 ms | 280.3 us | 1.032x | ok `4b2f43bf0f878408` |
| `fp.p256.mul_karatsuba.d2048` | 11.168 ms | 1.819 ms | 11.537 ms | 10.548 ms | 1.191 ms | 0.968x | ok `b4b89a80c2a3aaa9` |
| `fp.p256.mul_fft.d4096` | 15.235 ms | 485.5 us | 15.141 ms | 14.735 ms | 1.359 ms | 1.006x | ok `60ed510564513dfb` |
| `fp.p256.mul_fft.d8192` | 32.342 ms | 4.241 ms | 31.930 ms | 31.285 ms | 1.088 ms | 1.013x | ok `47812935df4de5fe` |
| `fp.p256.mul_auto.d2048` | 3.584 ms | 905.7 us | 3.496 ms | 3.087 ms | 192.1 us | 1.025x | ok `3d89402081500506` |
| `fp.p256.divmod.d1000_500` | 3.390 ms | 326.7 us | 3.211 ms | 3.065 ms | 235.5 us | 1.056x | ok `0bebce51d98eb3f4` |
| `fp.p256.gcd.d300` | 10.457 ms | 243.6 us | 10.386 ms | 10.191 ms | 277.5 us | 1.007x | ok `f69e4050cb1a5861` |
| `fp.p256.gcd_ext.d200` | 11.807 ms | 303.1 us | 11.802 ms | 11.577 ms | 265.1 us | 1.000x | ok `1a9df499c28b9f49` |
| `fp.p256.pow_mod.d64_k1` | 66.133 ms | 579.0 us | 65.272 ms | 64.685 ms | 957.9 us | 1.013x | ok `dbe2b595567fbb9c` |
| `fp.p256.is_irreducible.d12` | 6.882 ms | 225.2 us | 6.839 ms | 6.561 ms | 154.3 us | 1.006x | ok `7b839a40c98b0ef0` |
| `fp.p256.factor.d12` | 14.372 ms | 320.4 us | 14.399 ms | 13.821 ms | 779.4 us | 0.998x | ok `c33e3512ed388a7e` |
| `fq.q3e30.elem_muladd.n4000` | 2.093 ms | 63.2 us | 1.804 ms | 1.751 ms | 42.8 us | 1.160x | ok `16edfcc06a08028a` |
| `fq.q3e30.elem_inv.n100` | 149.3 us | 5.1 us | 149.7 us | 144.7 us | 3.5 us | 0.997x | ok `a8167c01d1cfc5dd` |
| `fq.q3e30.mul_school.d32` | 154.3 us | 4.1 us | 151.5 us | 146.8 us | 3.9 us | 1.018x | ok `ed0ca8eb7d325de3` |
| `fq.q3e30.mul_karatsuba.d128` | 753.0 us | 94.9 us | 744.0 us | 633.7 us | 117.0 us | 1.012x | ok `2464b901bd242c7d` |
| `fq.q3e30.divmod.d64_32` | 169.9 us | 5.4 us | 164.1 us | 159.2 us | 3.5 us | 1.035x | ok `a2b35f414a139741` |
| `fq.q3e30.factor.d5` | 1.180 ms | 38.9 us | 1.089 ms | 1.059 ms | 33.5 us | 1.084x | ok `417a5b753c7a2a07` |
| `fq.q2e8.elem_muladd.n50000` | 148.2 us | 3.6 us | 147.5 us | 144.0 us | 4.0 us | 1.005x | ok `c0ffe16ee54f128c` |
| `fq.q2e8.elem_inv.n5000` | 5.4 us | 0.5 us | 5.8 us | 5.6 us | 0.1 us | 0.931x | ok `81d71b5687e3e735` |
| `fq.q2e8.mul_school.d128` | 12.7 us | 0.3 us | 12.8 us | 12.1 us | 0.3 us | 0.992x | ok `c056d9d992271db9` |
| `fq.q2e8.mul_karatsuba.d512` | 92.5 us | 5.4 us | 93.1 us | 89.9 us | 3.0 us | 0.994x | ok `de1089fdce26d02c` |
| `fq.q2e8.is_irreducible.d16` | 38.8 us | 1.0 us | 39.3 us | 37.9 us | 1.0 us | 0.987x | ok `9478d590a29a9cbf` |
| `fq.q2e8.factor.d28` | 127.6 us | 5.1 us | 129.5 us | 124.7 us | 4.6 us | 0.985x | ok `d2fc95c10ececedd` |
| `fq.q7e4.elem_muladd.n50000` | 316.2 us | 9.5 us | 315.9 us | 308.5 us | 7.7 us | 1.001x | ok `a6cf18a43d40c5c0` |
| `fq.q7e4.elem_inv.n5000` | 22.2 us | 1.3 us | 22.5 us | 21.4 us | 1.2 us | 0.987x | ok `8664ac6d09409f94` |
| `fq.q7e4.mul_school.d128` | 76.4 us | 2.8 us | 76.6 us | 65.3 us | 2.4 us | 0.997x | ok `71d186560b1d8153` |
| `fq.q7e4.mul_karatsuba.d512` | 535.9 us | 15.6 us | 351.0 us | 294.9 us | 68.0 us | 1.527x | ok `6a2e9b751ab48e99` |
| `fq.q7e4.gcd.d200` | 235.2 us | 8.6 us | 235.2 us | 226.7 us | 7.6 us | 1.000x | ok `8cf19ac90e512feb` |
| `fq.q7e4.factor.d20` | 130.0 us | 9.3 us | 129.6 us | 122.4 us | 6.7 us | 1.003x | ok `ce23318b1b5fd48e` |
| `fq.m61e2.elem_muladd.n100000` | 9.532 ms | 767.8 us | 9.291 ms | 9.038 ms | 231.7 us | 1.026x | ok `1814347c960518d8` |
| `fq.m61e2.elem_inv.n1000` | 398.9 us | 37.6 us | 399.8 us | 384.4 us | 15.1 us | 0.998x | ok `45326d8095032432` |
| `fq.m61e2.mul_school.d128` | 121.1 us | 9.6 us | 122.2 us | 117.8 us | 3.7 us | 0.991x | ok `2cb38bd639886e3c` |
| `fq.m61e2.mul_karatsuba.d512` | 1.234 ms | 309.1 us | 1.149 ms | 895.6 us | 223.1 us | 1.073x | ok `f138ad86c0a7a838` |
| `fq.m61e2.mul_fft.d1024` | 12.954 ms | 1.020 ms | 12.547 ms | 11.906 ms | 507.7 us | 1.032x | ok `a45642d67f5a9912` |
| `fq.m61e2.mul_fft.d4096` | 35.797 ms | 2.692 ms | 35.702 ms | 34.169 ms | 743.4 us | 1.003x | ok `51844e8f95a9c0ab` |
| `fq.m61e2.mul_auto.d2048` | 2.710 ms | 198.6 us | 2.690 ms | 2.249 ms | 346.9 us | 1.008x | ok `76e66ad22fb3fa2e` |
| `fq.m61e2.divmod.d500_250` | 468.2 us | 10.9 us | 466.6 us | 448.6 us | 10.5 us | 1.003x | ok `53509a111940ee2a` |
| `fq.m61e2.gcd.d150` | 1.398 ms | 39.1 us | 1.403 ms | 1.351 ms | 32.2 us | 0.997x | ok `bd8c3c98cc56df84` |
| `fq.m61e2.factor.d7` | 2.556 ms | 64.2 us | 2.534 ms | 2.461 ms | 223.5 us | 1.009x | ok `d61430d11ee187c3` |
| `fq.build.q3e30` | 100.2 us | 6.1 us | 76.9 us | 74.1 us | 2.9 us | 1.303x | ok `4973a3e567e2fc74` |
| `fq.build.q7e20` | 97.6 us | 2.9 us | 84.1 us | 80.0 us | 3.6 us | 1.161x | ok `9c377cd554aa108e` |
| `fast.ntt.dft.k10` | 1.220 ms | 33.6 us | 1.241 ms | 1.189 ms | 49.3 us | 0.983x | ok `b61f454e63ae6d8d` |
| `fast.ntt.dft.k11` | 1.229 ms | 28.6 us | 1.240 ms | 1.197 ms | 45.2 us | 0.991x | ok `311d4fed0f42034c` |
| `fast.ntt.dft.k12` | 1.256 ms | 35.8 us | 1.272 ms | 1.215 ms | 36.0 us | 0.987x | ok `70691e178cfc0c9d` |
| `fast.ntt.idft.k12` | 1.266 ms | 47.2 us | 1.263 ms | 1.228 ms | 26.7 us | 1.002x | ok `e0677b59c4e33615` |
| `fast.ntt.toeplitz_lower_inverse.n1024` | 42.3 us | 2.5 us | 45.7 us | 42.3 us | 2.1 us | 0.926x | ok `7749eb8ff8cac3d0` |
| `fast.ntt.toeplitz_lower_inverse.n8192` | 777.9 us | 64.9 us | 809.1 us | 712.5 us | 99.4 us | 0.961x | ok `0d31554b19d8b996` |
| `fast.ntt.toeplitz_lower_vec.n4096` | 149.3 us | 18.3 us | 150.4 us | 125.7 us | 21.6 us | 0.993x | ok `86e36ad2a1eacb05` |
| `fast.ntt.toeplitz_divmod.d2000_1000` | 2.620 ms | 61.7 us | 2.668 ms | 2.599 ms | 41.6 us | 0.982x | ok `b138cdcceebcc32e` |
| `fast.p7.toeplitz_lower_inverse.n1024` | 57.1 us | 3.5 us | 57.5 us | 54.6 us | 3.4 us | 0.993x | ok `5df08c92a814155a` |
| `fast.p7.toeplitz_divmod.d2000_1000` | 1.245 ms | 29.0 us | 1.189 ms | 1.153 ms | 31.9 us | 1.047x | ok `433d13716b2f95e8` |
