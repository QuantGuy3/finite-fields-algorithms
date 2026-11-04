# Comparativa de benchmarks alcp

Base: **final-ref-O3**. Speedup = mediana(base) / mediana(variante); > 1 significa más rápido.
IQR = rango intercuartílico de las muestras. Las medias geométricas excluyen escenarios con error.

- **final-ref-O3**: modo interleave, 2026-11-04T18:49:25+01:00 (nota: fase3 medición final: baseline-original -O3 frente a HEAD 32d6922 (-O3 portable y -O3 -march=native))  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3`  
  include `~\AppData\Local\Temp\claude\C--Users-<usuario>-Documents-Archivos-Proyectos\2e18e96a-361a-47be-ad99-f0341eb0cc55\scratchpad\fin\include` @ ?, sha256 9b306a5d8af5  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: sí
- **final-O3**: modo interleave, 2026-11-04T18:49:39+01:00 (nota: fase3 medición final: baseline-original -O3 frente a HEAD 32d6922 (-O3 portable y -O3 -march=native))  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-62-g32d6922, sha256 af0c61d3df27  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: sí
- **final-native**: modo interleave, 2026-11-04T18:49:53+01:00 (nota: fase3 medición final: baseline-original -O3 frente a HEAD 32d6922 (-O3 portable y -O3 -march=native))  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -march=native`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-62-g32d6922, sha256 af0c61d3df27  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: sí

Checksums: todos los escenarios coinciden con la base.

## Media geométrica del speedup

| ámbito | final-O3 | final-native |
|---|---|---|
| **global** | 24.362x (n=128) | 24.570x (n=128) |
| grupo `fast` | 26.341x (n=10) | 26.492x (n=10) |
| grupo `fp.m127` | 3.432x (n=14) | 3.482x (n=14) |
| grupo `fp.m61` | 43.427x (n=14) | 44.401x (n=14) |
| grupo `fp.ntt` | 52.327x (n=16) | 53.105x (n=16) |
| grupo `fp.p101` | 60.352x (n=14) | 61.194x (n=14) |
| grupo `fp.p256` | 3.809x (n=16) | 3.722x (n=16) |
| grupo `fp.p7` | 61.648x (n=14) | 62.602x (n=14) |
| grupo `fq.build` | 46.152x (n=2) | 45.390x (n=2) |
| grupo `fq.m61e2` | 17.670x (n=10) | 17.945x (n=10) |
| grupo `fq.q2e8` | 48.594x (n=6) | 49.585x (n=6) |
| grupo `fq.q3e30` | 58.622x (n=6) | 59.152x (n=6) |
| grupo `fq.q7e4` | 39.140x (n=6) | 38.796x (n=6) |
| operación `build` | 46.152x (n=2) | 45.390x (n=2) |
| operación `dft` | 12.794x (n=3) | 12.977x (n=3) |
| operación `divmod` | 58.553x (n=8) | 60.423x (n=8) |
| operación `elem_inv` | 8.044x (n=10) | 8.115x (n=10) |
| operación `elem_muladd` | 10.453x (n=10) | 10.544x (n=10) |
| operación `elem_pow` | 2.615x (n=6) | 2.645x (n=6) |
| operación `factor` | 19.201x (n=10) | 19.178x (n=10) |
| operación `gcd` | 11.986x (n=8) | 12.016x (n=8) |
| operación `gcd_ext` | 19.048x (n=6) | 19.209x (n=6) |
| operación `idft` | 13.954x (n=1) | 13.869x (n=1) |
| operación `is_irreducible` | 29.943x (n=7) | 30.104x (n=7) |
| operación `mul_auto` | 34.698x (n=7) | 35.025x (n=7) |
| operación `mul_fft` | 14.150x (n=6) | 13.745x (n=6) |
| operación `mul_karatsuba` | 69.588x (n=16) | 70.229x (n=16) |
| operación `mul_school` | 57.409x (n=16) | 59.082x (n=16) |
| operación `pow_mod` | 35.282x (n=6) | 35.156x (n=6) |
| operación `toeplitz_divmod` | 31.044x (n=2) | 32.616x (n=2) |
| operación `toeplitz_lower_inverse` | 50.540x (n=3) | 49.293x (n=3) |
| operación `toeplitz_lower_vec` | 44.234x (n=1) | 44.090x (n=1) |

## Detalle por escenario

| escenario | base mediana | base IQR | final-O3 mediana | final-O3 mín | final-O3 IQR | speedup | final-native mediana | final-native mín | final-native IQR | speedup | checksum |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `fp.p7.elem_muladd.n400000` | 40.406 ms | 2.587 ms | 1.721 ms | 1.677 ms | 23.0 us | 23.480x | 1.716 ms | 1.676 ms | 42.2 us | 23.552x | ok `cf4797a827e7f9ec` |
| `fp.p7.elem_inv.n100000` | 12.543 ms | 734.5 us | 344.8 us | 334.2 us | 6.2 us | 36.373x | 343.7 us | 334.2 us | 9.5 us | 36.495x | ok `7beb5a3238b919cb` |
| `fp.p7.elem_pow.n100000` | 14.054 ms | 867.0 us | 964.7 us | 932.3 us | 46.1 us | 14.568x | 967.3 us | 908.3 us | 147.9 us | 14.529x | ok `a7f8c429d0d81b73` |
| `fp.p7.mul_school.d64` | 7.466 ms | 400.1 us | 69.6 us | 67.5 us | 0.9 us | 107.264x | 67.0 us | 65.0 us | 1.9 us | 111.427x | ok `559357b5bc3c57c2` |
| `fp.p7.mul_school.d512` | 23.225 ms | 1.399 ms | 167.6 us | 162.7 us | 4.0 us | 138.576x | 157.4 us | 153.0 us | 4.3 us | 147.557x | ok `5936b7480b151eb3` |
| `fp.p7.mul_karatsuba.d512` | 10.840 ms | 542.0 us | 129.7 us | 122.1 us | 15.5 us | 83.576x | 138.7 us | 123.2 us | 21.3 us | 78.152x | ok `d72c59b658ac0379` |
| `fp.p7.mul_karatsuba.d4096` | 303.047 ms | 8.794 ms | 1.374 ms | 1.054 ms | 195.5 us | 220.590x | 1.326 ms | 1.054 ms | 148.3 us | 228.525x | ok `00f3dde93cbd25a4` |
| `fp.p7.mul_auto.d2048` | 98.294 ms | 2.182 ms | 720.0 us | 552.7 us | 154.9 us | 136.519x | 702.9 us | 510.6 us | 151.3 us | 139.850x | ok `2b256c747e0070c7` |
| `fp.p7.divmod.d2000_1000` | 85.479 ms | 4.407 ms | 665.9 us | 640.2 us | 31.3 us | 128.366x | 630.7 us | 614.7 us | 15.9 us | 135.530x | ok `958d383cda71b699` |
| `fp.p7.gcd.d500` | 24.332 ms | 848.6 us | 867.4 us | 826.2 us | 53.3 us | 28.051x | 852.1 us | 829.1 us | 38.8 us | 28.555x | ok `cdf18334cf9b13d7` |
| `fp.p7.gcd_ext.d400` | 49.499 ms | 1.554 ms | 1.890 ms | 1.775 ms | 165.9 us | 26.187x | 1.857 ms | 1.771 ms | 95.8 us | 26.657x | ok `984ac156886a3f82` |
| `fp.p7.pow_mod.d128_k48` | 461.880 ms | 26.076 ms | 4.849 ms | 4.550 ms | 250.7 us | 95.257x | 4.784 ms | 4.639 ms | 240.2 us | 96.549x | ok `e3cedb8421cd962f` |
| `fp.p7.is_irreducible.d64` | 266.619 ms | 11.094 ms | 3.696 ms | 3.496 ms | 226.4 us | 72.141x | 3.683 ms | 3.547 ms | 169.4 us | 72.391x | ok `9cb1bec87bc433d0` |
| `fp.p7.factor.d62` | 10.882 ms | 1.630 ms | 286.6 us | 263.4 us | 16.7 us | 37.969x | 283.4 us | 270.3 us | 14.4 us | 38.398x | ok `3425aecb393956ba` |
| `fp.p101.elem_muladd.n400000` | 38.021 ms | 2.572 ms | 1.734 ms | 1.676 ms | 65.1 us | 21.933x | 1.722 ms | 1.676 ms | 42.7 us | 22.076x | ok `e51aa870ce6269eb` |
| `fp.p101.elem_inv.n100000` | 13.596 ms | 560.7 us | 686.0 us | 663.3 us | 26.2 us | 19.819x | 683.5 us | 663.5 us | 18.2 us | 19.891x | ok `69d137d626075364` |
| `fp.p101.elem_pow.n100000` | 17.155 ms | 479.3 us | 2.113 ms | 2.016 ms | 140.3 us | 8.120x | 2.027 ms | 1.968 ms | 54.3 us | 8.461x | ok `09fec44561c97afa` |
| `fp.p101.mul_school.d64` | 8.790 ms | 490.9 us | 77.2 us | 74.5 us | 3.0 us | 113.861x | 74.1 us | 71.8 us | 1.8 us | 118.624x | ok `5a0fa852855ea8b0` |
| `fp.p101.mul_school.d512` | 28.155 ms | 715.3 us | 198.1 us | 191.5 us | 3.4 us | 142.126x | 186.0 us | 180.2 us | 6.7 us | 151.372x | ok `8289633c899317a1` |
| `fp.p101.mul_karatsuba.d512` | 12.058 ms | 735.4 us | 124.7 us | 119.3 us | 3.9 us | 96.697x | 120.9 us | 116.4 us | 4.5 us | 99.736x | ok `fea453dad54c0188` |
| `fp.p101.mul_karatsuba.d4096` | 339.460 ms | 6.619 ms | 1.451 ms | 1.176 ms | 259.4 us | 233.997x | 1.367 ms | 1.081 ms | 217.5 us | 248.243x | ok `82f83f0b7fe7713b` |
| `fp.p101.mul_auto.d2048` | 110.976 ms | 3.162 ms | 708.8 us | 558.9 us | 167.4 us | 156.579x | 726.6 us | 534.4 us | 149.2 us | 152.743x | ok `3631fff5520b84d9` |
| `fp.p101.divmod.d2000_1000` | 102.765 ms | 2.192 ms | 772.7 us | 748.0 us | 14.6 us | 132.995x | 739.2 us | 714.7 us | 19.1 us | 139.012x | ok `ce32beed1dda91cd` |
| `fp.p101.gcd.d500` | 29.470 ms | 1.305 ms | 1.098 ms | 1.056 ms | 51.1 us | 26.850x | 1.097 ms | 1.053 ms | 73.9 us | 26.863x | ok `5cad0c7c5f168edb` |
| `fp.p101.gcd_ext.d400` | 58.965 ms | 3.786 ms | 2.302 ms | 2.211 ms | 204.0 us | 25.618x | 2.279 ms | 2.176 ms | 138.9 us | 25.870x | ok `4fe0534b04a103a1` |
| `fp.p101.pow_mod.d128_k20` | 562.083 ms | 15.165 ms | 5.360 ms | 5.196 ms | 157.0 us | 104.860x | 5.478 ms | 5.281 ms | 253.9 us | 102.606x | ok `bc5334d07b025488` |
| `fp.p101.is_irreducible.d48` | 319.544 ms | 9.193 ms | 3.943 ms | 3.826 ms | 114.7 us | 81.039x | 4.086 ms | 3.896 ms | 134.6 us | 78.212x | ok `5146dc1467eee05b` |
| `fp.p101.factor.d100` | 79.335 ms | 1.962 ms | 1.468 ms | 1.413 ms | 63.1 us | 54.061x | 1.503 ms | 1.423 ms | 95.4 us | 52.788x | ok `48a7604aa4ae63ee` |
| `fp.ntt.elem_muladd.n400000` | 41.541 ms | 976.1 us | 2.529 ms | 2.467 ms | 66.2 us | 16.427x | 2.514 ms | 2.444 ms | 79.3 us | 16.525x | ok `682b833fe4edd17b` |
| `fp.ntt.elem_inv.n100000` | 16.898 ms | 606.8 us | 3.035 ms | 2.918 ms | 159.7 us | 5.568x | 3.027 ms | 2.942 ms | 74.7 us | 5.583x | ok `d275b741e6b3050e` |
| `fp.ntt.elem_pow.n50000` | 13.926 ms | 309.9 us | 6.805 ms | 6.705 ms | 189.9 us | 2.046x | 6.811 ms | 6.679 ms | 97.1 us | 2.045x | ok `c27b88ab2935c700` |
| `fp.ntt.mul_school.d64` | 9.476 ms | 1.559 ms | 79.2 us | 76.6 us | 2.2 us | 119.649x | 75.1 us | 73.6 us | 1.8 us | 126.181x | ok `67c7e1645a9cc6ae` |
| `fp.ntt.mul_school.d512` | 28.689 ms | 1.133 ms | 199.4 us | 193.6 us | 3.5 us | 143.878x | 187.5 us | 182.0 us | 5.4 us | 153.009x | ok `e296393a6a637046` |
| `fp.ntt.mul_karatsuba.d512` | 11.893 ms | 283.3 us | 125.2 us | 119.6 us | 10.3 us | 94.989x | 123.3 us | 117.4 us | 8.7 us | 96.453x | ok `964e15063fb16df4` |
| `fp.ntt.mul_karatsuba.d4096` | 345.066 ms | 8.693 ms | 1.675 ms | 1.282 ms | 275.1 us | 206.047x | 1.592 ms | 1.205 ms | 307.1 us | 216.805x | ok `d656b678f0c2ca16` |
| `fp.ntt.mul_fft.d4096` | 86.099 ms | 3.414 ms | 1.649 ms | 1.481 ms | 168.1 us | 52.200x | 1.627 ms | 1.456 ms | 144.4 us | 52.911x | ok `f4c81b85220fb683` |
| `fp.ntt.mul_fft.d32768` | 1.020 s | 45.502 ms | 15.171 ms | 14.553 ms | 358.4 us | 67.218x | 15.572 ms | 14.684 ms | 1.063 ms | 65.487x | ok `8228e76950faded3` |
| `fp.ntt.mul_auto.d2048` | 33.280 ms | 2.196 ms | 736.2 us | 690.7 us | 61.5 us | 45.202x | 748.7 us | 694.1 us | 22.4 us | 44.451x | ok `02e67c888996018d` |
| `fp.ntt.divmod.d2000_1000` | 111.025 ms | 2.835 ms | 791.7 us | 764.7 us | 35.0 us | 140.236x | 758.7 us | 732.1 us | 22.7 us | 146.336x | ok `7adfbb2614325099` |
| `fp.ntt.gcd.d500` | 31.526 ms | 795.3 us | 1.161 ms | 1.113 ms | 50.6 us | 27.144x | 1.173 ms | 1.124 ms | 44.3 us | 26.881x | ok `678962135e7a2447` |
| `fp.ntt.gcd_ext.d400` | 441.073 ms | 11.528 ms | 2.420 ms | 2.314 ms | 136.1 us | 182.284x | 2.369 ms | 2.290 ms | 143.5 us | 186.209x | ok `f9dc4e400a7954ee` |
| `fp.ntt.pow_mod.d128_k5` | 437.299 ms | 17.976 ms | 5.075 ms | 4.920 ms | 109.7 us | 86.171x | 5.107 ms | 4.899 ms | 207.8 us | 85.634x | ok `0eb2514aceb37261` |
| `fp.ntt.is_irreducible.d32` | 438.192 ms | 17.607 ms | 4.920 ms | 4.759 ms | 478.3 us | 89.060x | 4.775 ms | 4.651 ms | 127.6 us | 91.775x | ok `1f2827a0d8aeb65e` |
| `fp.ntt.factor.d40` | 65.231 ms | 1.825 ms | 1.100 ms | 1.058 ms | 33.6 us | 59.307x | 1.101 ms | 1.061 ms | 29.5 us | 59.253x | ok `297c76c423fe55aa` |
| `fp.m61.elem_muladd.n400000` | 38.740 ms | 764.3 us | 2.236 ms | 2.160 ms | 52.4 us | 17.329x | 2.205 ms | 2.155 ms | 56.8 us | 17.567x | ok `b8e325f73ad90582` |
| `fp.m61.elem_inv.n100000` | 20.767 ms | 502.7 us | 8.522 ms | 8.262 ms | 191.3 us | 2.437x | 8.452 ms | 8.300 ms | 256.2 us | 2.457x | ok `41e7fe38de6fdcc1` |
| `fp.m61.elem_pow.n20000` | 8.161 ms | 241.3 us | 6.178 ms | 6.016 ms | 99.9 us | 1.321x | 6.109 ms | 5.997 ms | 119.4 us | 1.336x | ok `a3680ece7fc07932` |
| `fp.m61.mul_school.d64` | 9.191 ms | 465.0 us | 82.1 us | 79.3 us | 2.3 us | 111.952x | 79.3 us | 76.6 us | 1.5 us | 115.905x | ok `bed6e9e86907de79` |
| `fp.m61.mul_school.d512` | 28.693 ms | 723.7 us | 206.7 us | 198.7 us | 3.6 us | 138.815x | 194.5 us | 188.4 us | 6.7 us | 147.522x | ok `43c3884710f7d468` |
| `fp.m61.mul_karatsuba.d512` | 11.852 ms | 483.1 us | 130.6 us | 125.3 us | 3.8 us | 90.747x | 128.3 us | 122.9 us | 9.8 us | 92.374x | ok `a2dcbeaa6574fcb8` |
| `fp.m61.mul_karatsuba.d4096` | 342.592 ms | 5.924 ms | 1.651 ms | 1.398 ms | 212.5 us | 207.468x | 1.591 ms | 1.207 ms | 214.7 us | 215.304x | ok `37d5deaa1cfeb600` |
| `fp.m61.mul_auto.d2048` | 115.703 ms | 2.599 ms | 884.1 us | 680.4 us | 204.2 us | 130.870x | 823.5 us | 637.0 us | 138.7 us | 140.510x | ok `73174ffa4e90d6fa` |
| `fp.m61.divmod.d2000_1000` | 111.933 ms | 2.871 ms | 815.4 us | 786.9 us | 20.2 us | 137.274x | 784.2 us | 755.6 us | 21.5 us | 142.736x | ok `a8710fc5aaa24d09` |
| `fp.m61.gcd.d500` | 30.952 ms | 1.232 ms | 1.016 ms | 969.4 us | 37.8 us | 30.470x | 1.003 ms | 966.4 us | 36.1 us | 30.851x | ok `fda394dde025c1a2` |
| `fp.m61.gcd_ext.d400` | 60.673 ms | 1.599 ms | 2.104 ms | 2.023 ms | 91.7 us | 28.837x | 2.075 ms | 2.004 ms | 73.5 us | 29.236x | ok `be3d81e8acc44b48` |
| `fp.m61.pow_mod.d64_k2` | 144.598 ms | 9.549 ms | 1.501 ms | 1.432 ms | 55.3 us | 96.367x | 1.500 ms | 1.458 ms | 49.0 us | 96.392x | ok `394fb320787700fe` |
| `fp.m61.is_irreducible.d24` | 464.373 ms | 14.595 ms | 5.870 ms | 5.725 ms | 175.4 us | 79.107x | 5.879 ms | 5.731 ms | 197.3 us | 78.988x | ok `9047cf1d3fc29abf` |
| `fp.m61.factor.d26` | 50.410 ms | 1.571 ms | 1.167 ms | 1.129 ms | 41.0 us | 43.200x | 1.176 ms | 1.139 ms | 57.8 us | 42.851x | ok `931b8d0b8b91b0b4` |
| `fp.m127.elem_muladd.n200000` | 22.649 ms | 640.6 us | 23.868 ms | 22.873 ms | 1.314 ms | 0.949x | 23.011 ms | 22.548 ms | 865.5 us | 0.984x | ok `16fab7b139aa9aab` |
| `fp.m127.elem_inv.n50000` | 31.349 ms | 1.239 ms | 31.541 ms | 28.600 ms | 1.187 ms | 0.994x | 30.103 ms | 27.735 ms | 1.342 ms | 1.041x | ok `c7d67d6018036bf1` |
| `fp.m127.elem_pow.n10000` | 23.126 ms | 523.4 us | 23.017 ms | 22.513 ms | 574.1 us | 1.005x | 22.890 ms | 22.510 ms | 382.6 us | 1.010x | ok `308896e55b40c815` |
| `fp.m127.mul_school.d64` | 9.800 ms | 309.9 us | 1.664 ms | 1.587 ms | 74.6 us | 5.888x | 1.639 ms | 1.564 ms | 91.3 us | 5.979x | ok `a9664408b0fb88e7` |
| `fp.m127.mul_school.d512` | 30.611 ms | 1.125 ms | 3.962 ms | 3.797 ms | 179.3 us | 7.726x | 3.940 ms | 3.744 ms | 164.9 us | 7.769x | ok `98eca357802203c0` |
| `fp.m127.mul_karatsuba.d512` | 13.878 ms | 297.6 us | 1.883 ms | 1.483 ms | 272.8 us | 7.370x | 1.862 ms | 1.527 ms | 234.1 us | 7.454x | ok `f773c153fde34172` |
| `fp.m127.mul_karatsuba.d2048` | 133.990 ms | 5.077 ms | 10.215 ms | 9.755 ms | 229.8 us | 13.117x | 10.208 ms | 9.747 ms | 332.2 us | 13.125x | ok `695e0d0f13c14e87` |
| `fp.m127.mul_auto.d2048` | 134.225 ms | 4.254 ms | 10.377 ms | 9.782 ms | 396.8 us | 12.935x | 10.309 ms | 9.892 ms | 495.5 us | 13.020x | ok `1bc667b834b9e559` |
| `fp.m127.divmod.d2000_1000` | 114.183 ms | 5.082 ms | 14.878 ms | 14.120 ms | 569.9 us | 7.675x | 14.909 ms | 14.006 ms | 652.3 us | 7.659x | ok `a8babdcd119eb404` |
| `fp.m127.gcd.d400` | 20.845 ms | 630.0 us | 14.644 ms | 14.376 ms | 595.1 us | 1.423x | 14.495 ms | 14.018 ms | 524.1 us | 1.438x | ok `1b15fac38c54accc` |
| `fp.m127.gcd_ext.d300` | 38.173 ms | 2.251 ms | 24.770 ms | 23.712 ms | 1.704 ms | 1.541x | 24.244 ms | 23.121 ms | 1.097 ms | 1.575x | ok `45a1b9001453c6ef` |
| `fp.m127.pow_mod.d64_k1` | 216.862 ms | 11.435 ms | 46.121 ms | 44.753 ms | 997.8 us | 4.702x | 45.593 ms | 44.021 ms | 2.156 ms | 4.756x | ok `9e9e535b3c4e4e7d` |
| `fp.m127.is_irreducible.d16` | 324.543 ms | 28.493 ms | 98.568 ms | 95.805 ms | 4.621 ms | 3.293x | 95.470 ms | 93.016 ms | 3.037 ms | 3.399x | ok `2cfb69ee2746aa9d` |
| `fp.m127.factor.d19` | 81.982 ms | 3.033 ms | 36.634 ms | 35.452 ms | 1.036 ms | 2.238x | 36.534 ms | 35.405 ms | 4.072 ms | 2.244x | ok `7b0defab2417737c` |
| `fp.p256.elem_muladd.n200000` | 31.647 ms | 1.390 ms | 32.507 ms | 30.831 ms | 2.582 ms | 0.974x | 32.436 ms | 31.262 ms | 508.7 us | 0.976x | ok `cc3cab838e9b19f1` |
| `fp.p256.elem_inv.n50000` | 65.694 ms | 1.209 ms | 65.912 ms | 64.000 ms | 4.178 ms | 0.997x | 65.592 ms | 64.392 ms | 2.466 ms | 1.002x | ok `89f920a8eba0c743` |
| `fp.p256.elem_pow.n2000` | 19.268 ms | 440.0 us | 19.379 ms | 18.890 ms | 758.4 us | 0.994x | 19.097 ms | 18.810 ms | 444.7 us | 1.009x | ok `40472b256a3c877d` |
| `fp.p256.mul_school.d64` | 13.579 ms | 487.7 us | 2.449 ms | 2.395 ms | 111.1 us | 5.545x | 2.485 ms | 2.395 ms | 74.1 us | 5.464x | ok `f5d3d87e3ed0977e` |
| `fp.p256.mul_school.d512` | 42.851 ms | 1.611 ms | 6.237 ms | 6.023 ms | 311.8 us | 6.871x | 6.282 ms | 6.062 ms | 259.0 us | 6.821x | ok `6fd19bd7b3ceab69` |
| `fp.p256.mul_karatsuba.d512` | 18.006 ms | 693.8 us | 2.025 ms | 1.708 ms | 237.2 us | 8.891x | 2.101 ms | 1.745 ms | 253.5 us | 8.568x | ok `4b2f43bf0f878408` |
| `fp.p256.mul_karatsuba.d2048` | 169.346 ms | 5.200 ms | 12.362 ms | 11.894 ms | 346.5 us | 13.699x | 12.608 ms | 12.026 ms | 1.038 ms | 13.432x | ok `b4b89a80c2a3aaa9` |
| `fp.p256.mul_fft.d4096` | 102.178 ms | 6.096 ms | 17.667 ms | 15.056 ms | 2.378 ms | 5.783x | 18.348 ms | 14.844 ms | 1.581 ms | 5.569x | ok `60ed510564513dfb` |
| `fp.p256.mul_fft.d8192` | 221.459 ms | 11.420 ms | 34.617 ms | 30.881 ms | 7.201 ms | 6.397x | 38.807 ms | 31.656 ms | 7.808 ms | 5.707x | ok `47812935df4de5fe` |
| `fp.p256.mul_auto.d2048` | 42.971 ms | 1.764 ms | 8.446 ms | 6.891 ms | 308.7 us | 5.088x | 8.465 ms | 7.161 ms | 772.2 us | 5.076x | ok `3d89402081500506` |
| `fp.p256.divmod.d1000_500` | 39.697 ms | 3.676 ms | 6.029 ms | 5.807 ms | 286.8 us | 6.585x | 6.147 ms | 5.837 ms | 395.6 us | 6.458x | ok `0bebce51d98eb3f4` |
| `fp.p256.gcd.d300` | 15.937 ms | 932.2 us | 10.623 ms | 10.310 ms | 270.5 us | 1.500x | 11.157 ms | 10.444 ms | 708.6 us | 1.428x | ok `f69e4050cb1a5861` |
| `fp.p256.gcd_ext.d200` | 117.905 ms | 6.405 ms | 13.414 ms | 13.086 ms | 627.8 us | 8.789x | 13.872 ms | 13.084 ms | 648.6 us | 8.499x | ok `1a9df499c28b9f49` |
| `fp.p256.pow_mod.d64_k1` | 390.662 ms | 15.978 ms | 78.987 ms | 77.653 ms | 1.242 ms | 4.946x | 80.487 ms | 77.999 ms | 3.344 ms | 4.854x | ok `dbe2b595567fbb9c` |
| `fp.p256.is_irreducible.d12` | 335.484 ms | 14.037 ms | 118.995 ms | 116.397 ms | 2.065 ms | 2.819x | 121.084 ms | 115.983 ms | 6.510 ms | 2.771x | ok `7b839a40c98b0ef0` |
| `fp.p256.factor.d12` | 34.762 ms | 1.331 ms | 18.040 ms | 17.229 ms | 1.292 ms | 1.927x | 18.407 ms | 17.432 ms | 970.1 us | 1.889x | ok `c33e3512ed388a7e` |
| `fq.q3e30.elem_muladd.n4000` | 195.733 ms | 8.494 ms | 5.657 ms | 5.463 ms | 416.2 us | 34.601x | 5.679 ms | 5.377 ms | 616.9 us | 34.464x | ok `16edfcc06a08028a` |
| `fq.q3e30.elem_inv.n100` | 329.537 ms | 5.378 ms | 9.928 ms | 9.535 ms | 850.7 us | 33.194x | 9.652 ms | 9.370 ms | 457.9 us | 34.140x | ok `a8167c01d1cfc5dd` |
| `fq.q3e30.mul_school.d32` | 51.708 ms | 1.784 ms | 728.2 us | 711.2 us | 20.0 us | 71.007x | 733.4 us | 707.7 us | 29.7 us | 70.504x | ok `ed0ca8eb7d325de3` |
| `fq.q3e30.mul_karatsuba.d128` | 414.091 ms | 13.251 ms | 2.287 ms | 2.053 ms | 135.5 us | 181.027x | 2.272 ms | 2.001 ms | 173.6 us | 182.238x | ok `2464b901bd242c7d` |
| `fq.q3e30.divmod.d64_32` | 56.404 ms | 1.770 ms | 803.8 us | 780.0 us | 11.9 us | 70.172x | 780.0 us | 749.9 us | 50.2 us | 72.313x | ok `a2b35f414a139741` |
| `fq.q3e30.factor.d5` | 424.297 ms | 15.472 ms | 10.831 ms | 10.474 ms | 497.5 us | 39.176x | 10.829 ms | 10.283 ms | 997.3 us | 39.183x | ok `417a5b753c7a2a07` |
| `fq.q2e8.elem_muladd.n50000` | 214.448 ms | 5.037 ms | 9.697 ms | 8.587 ms | 1.473 ms | 22.116x | 9.039 ms | 8.473 ms | 305.0 us | 23.724x | ok `c0ffe16ee54f128c` |
| `fq.q2e8.elem_inv.n5000` | 238.362 ms | 15.256 ms | 11.749 ms | 11.296 ms | 245.2 us | 20.288x | 11.679 ms | 11.224 ms | 602.1 us | 20.410x | ok `81d71b5687e3e735` |
| `fq.q2e8.mul_school.d128` | 64.960 ms | 2.082 ms | 596.1 us | 568.6 us | 32.7 us | 108.965x | 593.6 us | 571.0 us | 32.9 us | 109.433x | ok `c056d9d992271db9` |
| `fq.q2e8.mul_karatsuba.d512` | 347.915 ms | 15.739 ms | 2.343 ms | 1.972 ms | 324.2 us | 148.516x | 2.349 ms | 1.896 ms | 289.9 us | 148.083x | ok `de1089fdce26d02c` |
| `fq.q2e8.is_irreducible.d16` | 387.433 ms | 16.401 ms | 6.864 ms | 6.641 ms | 611.2 us | 56.440x | 6.684 ms | 6.402 ms | 511.5 us | 57.961x | ok `9478d590a29a9cbf` |
| `fq.q2e8.factor.d28` | 94.381 ms | 4.136 ms | 2.937 ms | 2.836 ms | 191.6 us | 32.131x | 2.888 ms | 2.774 ms | 255.9 us | 32.680x | ok `d2fc95c10ececedd` |
| `fq.q7e4.elem_muladd.n50000` | 95.928 ms | 3.021 ms | 6.015 ms | 5.877 ms | 89.1 us | 15.948x | 6.220 ms | 5.735 ms | 815.7 us | 15.423x | ok `a6cf18a43d40c5c0` |
| `fq.q7e4.elem_inv.n5000` | 190.181 ms | 4.515 ms | 10.614 ms | 10.216 ms | 326.2 us | 17.917x | 10.699 ms | 10.084 ms | 811.0 us | 17.775x | ok `8664ac6d09409f94` |
| `fq.q7e4.mul_school.d128` | 37.994 ms | 1.110 ms | 278.9 us | 267.4 us | 12.7 us | 136.230x | 281.5 us | 270.4 us | 8.2 us | 134.971x | ok `71d186560b1d8153` |
| `fq.q7e4.mul_karatsuba.d512` | 203.321 ms | 4.475 ms | 1.753 ms | 1.264 ms | 323.1 us | 115.991x | 1.816 ms | 1.556 ms | 324.7 us | 111.992x | ok `6a2e9b751ab48e99` |
| `fq.q7e4.gcd.d200` | 101.832 ms | 2.776 ms | 3.767 ms | 3.498 ms | 780.4 us | 27.035x | 3.672 ms | 3.439 ms | 195.7 us | 27.729x | ok `8cf19ac90e512feb` |
| `fq.q7e4.factor.d20` | 55.945 ms | 992.9 us | 1.900 ms | 1.833 ms | 61.5 us | 29.452x | 1.885 ms | 1.793 ms | 96.8 us | 29.673x | ok `ce23318b1b5fd48e` |
| `fq.m61e2.elem_muladd.n100000` | 83.801 ms | 3.844 ms | 8.893 ms | 8.567 ms | 294.5 us | 9.423x | 9.016 ms | 8.487 ms | 578.9 us | 9.294x | ok `1814347c960518d8` |
| `fq.m61e2.elem_inv.n1000` | 141.572 ms | 3.483 ms | 14.600 ms | 14.096 ms | 456.3 us | 9.696x | 14.713 ms | 14.254 ms | 441.5 us | 9.622x | ok `45326d8095032432` |
| `fq.m61e2.mul_school.d128` | 14.260 ms | 1.121 ms | 120.4 us | 112.2 us | 4.0 us | 118.437x | 113.7 us | 108.8 us | 4.5 us | 125.416x | ok `2cb38bd639886e3c` |
| `fq.m61e2.mul_karatsuba.d512` | 78.097 ms | 3.225 ms | 1.458 ms | 1.047 ms | 315.5 us | 53.579x | 1.396 ms | 1.087 ms | 211.8 us | 55.947x | ok `f138ad86c0a7a838` |
| `fq.m61e2.mul_fft.d1024` | 82.650 ms | 1.766 ms | 14.730 ms | 14.049 ms | 1.274 ms | 5.611x | 14.890 ms | 13.799 ms | 2.363 ms | 5.551x | ok `a45642d67f5a9912` |
| `fq.m61e2.mul_fft.d4096` | 448.872 ms | 32.476 ms | 40.730 ms | 38.860 ms | 2.015 ms | 11.021x | 40.690 ms | 39.116 ms | 688.8 us | 11.031x | ok `51844e8f95a9c0ab` |
| `fq.m61e2.mul_auto.d2048` | 183.834 ms | 7.080 ms | 25.267 ms | 24.280 ms | 1.112 ms | 7.276x | 25.069 ms | 23.484 ms | 1.024 ms | 7.333x | ok `76e66ad22fb3fa2e` |
| `fq.m61e2.divmod.d500_250` | 55.390 ms | 4.313 ms | 467.2 us | 428.7 us | 195.9 us | 118.544x | 438.8 us | 417.2 us | 18.4 us | 126.215x | ok `53509a111940ee2a` |
| `fq.m61e2.gcd.d150` | 38.532 ms | 4.168 ms | 3.252 ms | 3.098 ms | 341.3 us | 11.848x | 3.213 ms | 3.087 ms | 136.0 us | 11.994x | ok `bd8c3c98cc56df84` |
| `fq.m61e2.factor.d7` | 156.930 ms | 10.606 ms | 19.370 ms | 18.999 ms | 397.6 us | 8.102x | 19.328 ms | 18.855 ms | 930.7 us | 8.119x | ok `d61430d11ee187c3` |
| `fq.build.q3e30` | 21.544 ms | 1.177 ms | 472.6 us | 453.6 us | 22.9 us | 45.585x | 482.4 us | 466.3 us | 24.4 us | 44.659x | ok `4973a3e567e2fc74` |
| `fq.build.q7e20` | 17.811 ms | 1.528 ms | 381.2 us | 368.1 us | 35.1 us | 46.725x | 386.1 us | 370.9 us | 15.9 us | 46.132x | ok `9c377cd554aa108e` |
| `fast.ntt.dft.k10` | 15.034 ms | 650.2 us | 1.236 ms | 1.186 ms | 202.4 us | 12.164x | 1.233 ms | 1.193 ms | 37.7 us | 12.191x | ok `b61f454e63ae6d8d` |
| `fast.ntt.dft.k11` | 15.829 ms | 799.0 us | 1.276 ms | 1.196 ms | 193.2 us | 12.401x | 1.237 ms | 1.200 ms | 24.6 us | 12.794x | ok `311d4fed0f42034c` |
| `fast.ntt.dft.k12` | 17.568 ms | 508.0 us | 1.266 ms | 1.213 ms | 58.5 us | 13.881x | 1.254 ms | 1.211 ms | 45.6 us | 14.012x | ok `70691e178cfc0c9d` |
| `fast.ntt.idft.k12` | 17.768 ms | 815.8 us | 1.273 ms | 1.223 ms | 72.3 us | 13.954x | 1.281 ms | 1.230 ms | 36.1 us | 13.869x | ok `e0677b59c4e33615` |
| `fast.ntt.toeplitz_lower_inverse.n1024` | 25.240 ms | 1.263 ms | 605.3 us | 583.4 us | 22.0 us | 41.698x | 619.2 us | 596.5 us | 26.4 us | 40.762x | ok `7749eb8ff8cac3d0` |
| `fast.ntt.toeplitz_lower_inverse.n8192` | 274.538 ms | 16.440 ms | 6.394 ms | 6.032 ms | 322.7 us | 42.939x | 6.833 ms | 6.245 ms | 909.7 us | 40.176x | ok `0d31554b19d8b996` |
| `fast.ntt.toeplitz_lower_vec.n4096` | 34.690 ms | 1.177 ms | 784.2 us | 720.7 us | 67.4 us | 44.234x | 786.8 us | 737.0 us | 26.3 us | 44.090x | ok `86e36ad2a1eacb05` |
| `fast.ntt.toeplitz_divmod.d2000_1000` | 59.339 ms | 1.280 ms | 2.824 ms | 2.735 ms | 64.8 us | 21.012x | 2.778 ms | 2.684 ms | 41.3 us | 21.359x | ok `b138cdcceebcc32e` |
| `fast.p7.toeplitz_lower_inverse.n1024` | 69.466 ms | 5.936 ms | 963.5 us | 925.2 us | 31.0 us | 72.101x | 949.8 us | 909.3 us | 51.0 us | 73.138x | ok `5df08c92a814155a` |
| `fast.p7.toeplitz_divmod.d2000_1000` | 75.448 ms | 5.653 ms | 1.645 ms | 1.584 ms | 77.0 us | 45.865x | 1.515 ms | 1.458 ms | 39.0 us | 49.804x | ok `433d13716b2f95e8` |

## Antes/después: casos históricos

Casos de `bench/bench.cpp` medidos con `bench/bench_historico.cpp` (mediana de 5 repeticiones por ronda, 2 rondas alternadas, `-O3`; datos en `bench/results/fase3/final-historico.json`). «Antes» son las cifras históricas que cita el encargo, que no se han vuelto a medir; la comparación válida es con «referencia hoy».

| caso | antes (histórico, ms) | referencia hoy (ms) | final GMP (ms) | final 64 bits (ms) | speedup ref → 64 bits |
|---|---|---|---|---|---|
| factor grado 5, q = 3^30 | 823 | 439.5 | 92.97 | 11.61 | 37.9x |
| producto grado 4000, escolar | 3752 | 1750.4 | 160.58 | 12.82 | 136.5x |
| producto grado 4000, Karatsuba | 606 | 331.8 | 20.51 | 1.54 | 214.9x |
| producto grado 4000, FFT | 62 | 33.8 | 1.12 | 0.72 | 46.8x |
| producto grado 4000, Auto | — | 33.5 | 1.07 | 0.73 | 45.9x |
