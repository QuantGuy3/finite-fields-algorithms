# Comparativa de benchmarks alcp

Base: **f4-ref**. Speedup = mediana(base) / mediana(variante); > 1 significa más rápido.
IQR = rango intercuartílico de las muestras. Las medias geométricas excluyen escenarios con error.

- **f4-ref**: modo interleave, 2026-11-07T18:45:30+01:00 (nota: fase 4 nivel B: original vs nivel A vs nivel B)  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3`  
  include `~\AppData\Local\Temp\claude\204af1c0d38157dc\2e18e96a-361a-47be-ad99-f0341eb0cc55\scratchpad\inc-ref\include` @ ?, sha256 9b306a5d8af5  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: sí
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

| ámbito | f4-f4 | f4-b |
|---|---|---|
| **global** | 29.378x (n=128) | 36.811x (n=128) |
| grupo `fast` | 26.826x (n=10) | 34.529x (n=10) |
| grupo `fp.m127` | 3.804x (n=14) | 5.859x (n=14) |
| grupo `fp.m61` | 47.370x (n=14) | 54.044x (n=14) |
| grupo `fp.ntt` | 58.339x (n=16) | 64.236x (n=16) |
| grupo `fp.p101` | 67.433x (n=14) | 100.527x (n=14) |
| grupo `fp.p256` | 4.103x (n=16) | 5.008x (n=16) |
| grupo `fp.p7` | 67.498x (n=14) | 99.152x (n=14) |
| grupo `fq.build` | 88.280x (n=2) | 81.206x (n=2) |
| grupo `fq.m61e2` | 33.025x (n=10) | 40.605x (n=10) |
| grupo `fq.q2e8` | 76.107x (n=6) | 80.553x (n=6) |
| grupo `fq.q3e30` | 98.580x (n=6) | 107.770x (n=6) |
| grupo `fq.q7e4` | 54.772x (n=6) | 54.507x (n=6) |
| operación `build` | 88.280x (n=2) | 81.206x (n=2) |
| operación `dft` | 12.953x (n=3) | 12.288x (n=3) |
| operación `divmod` | 61.102x (n=8) | 114.457x (n=8) |
| operación `elem_inv` | 21.055x (n=10) | 20.573x (n=10) |
| operación `elem_muladd` | 10.679x (n=10) | 10.554x (n=10) |
| operación `elem_pow` | 2.643x (n=6) | 2.684x (n=6) |
| operación `factor` | 25.640x (n=10) | 26.761x (n=10) |
| operación `gcd` | 13.998x (n=8) | 13.013x (n=8) |
| operación `gcd_ext` | 19.874x (n=6) | 26.049x (n=6) |
| operación `idft` | 14.371x (n=1) | 13.799x (n=1) |
| operación `is_irreducible` | 81.485x (n=7) | 277.781x (n=7) |
| operación `mul_auto` | 35.786x (n=7) | 133.506x (n=7) |
| operación `mul_fft` | 13.963x (n=6) | 14.060x (n=6) |
| operación `mul_karatsuba` | 70.313x (n=16) | 69.924x (n=16) |
| operación `mul_school` | 58.446x (n=16) | 57.529x (n=16) |
| operación `pow_mod` | 36.841x (n=6) | 58.543x (n=6) |
| operación `toeplitz_divmod` | 32.205x (n=2) | 35.320x (n=2) |
| operación `toeplitz_lower_inverse` | 51.736x (n=3) | 126.731x (n=3) |
| operación `toeplitz_lower_vec` | 43.028x (n=1) | 37.053x (n=1) |

## Detalle por escenario

| escenario | base mediana | base IQR | f4-f4 mediana | f4-f4 mín | f4-f4 IQR | speedup | f4-b mediana | f4-b mín | f4-b IQR | speedup | checksum |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `fp.p7.elem_muladd.n400000` | 40.591 ms | 3.249 ms | 1.723 ms | 1.676 ms | 55.2 us | 23.557x | 1.715 ms | 1.674 ms | 101.2 us | 23.670x | ok `cf4797a827e7f9ec` |
| `fp.p7.elem_inv.n100000` | 12.710 ms | 1.286 ms | 343.7 us | 334.1 us | 11.6 us | 36.980x | 341.7 us | 334.2 us | 15.3 us | 37.197x | ok `7beb5a3238b919cb` |
| `fp.p7.elem_pow.n100000` | 14.433 ms | 2.881 ms | 963.0 us | 924.4 us | 40.4 us | 14.987x | 948.9 us | 925.2 us | 47.2 us | 15.210x | ok `a7f8c429d0d81b73` |
| `fp.p7.mul_school.d64` | 7.872 ms | 865.8 us | 69.2 us | 67.1 us | 4.1 us | 113.763x | 70.9 us | 66.8 us | 8.2 us | 111.035x | ok `559357b5bc3c57c2` |
| `fp.p7.mul_school.d512` | 24.244 ms | 2.488 ms | 166.5 us | 161.8 us | 16.0 us | 145.607x | 166.2 us | 161.5 us | 16.2 us | 145.870x | ok `5936b7480b151eb3` |
| `fp.p7.mul_karatsuba.d512` | 10.876 ms | 1.029 ms | 136.2 us | 121.1 us | 21.7 us | 79.822x | 135.1 us | 122.3 us | 17.9 us | 80.501x | ok `d72c59b658ac0379` |
| `fp.p7.mul_karatsuba.d4096` | 322.072 ms | 36.550 ms | 1.422 ms | 1.076 ms | 174.6 us | 226.429x | 1.470 ms | 1.082 ms | 307.6 us | 219.104x | ok `00f3dde93cbd25a4` |
| `fp.p7.mul_auto.d2048` | 104.255 ms | 9.286 ms | 728.5 us | 528.6 us | 154.4 us | 143.110x | 76.9 us | 73.4 us | 10.0 us | 1355.726x | ok `2b256c747e0070c7` |
| `fp.p7.divmod.d2000_1000` | 91.553 ms | 11.128 ms | 670.3 us | 643.4 us | 62.9 us | 136.586x | 150.8 us | 143.2 us | 13.6 us | 607.319x | ok `958d383cda71b699` |
| `fp.p7.gcd.d500` | 26.439 ms | 2.050 ms | 866.1 us | 827.8 us | 76.3 us | 30.527x | 954.8 us | 889.9 us | 104.8 us | 27.691x | ok `cdf18334cf9b13d7` |
| `fp.p7.gcd_ext.d400` | 51.268 ms | 3.287 ms | 1.874 ms | 1.797 ms | 210.5 us | 27.360x | 1.148 ms | 1.028 ms | 136.6 us | 44.658x | ok `984ac156886a3f82` |
| `fp.p7.pow_mod.d128_k48` | 474.926 ms | 46.096 ms | 4.740 ms | 4.614 ms | 349.6 us | 100.193x | 1.783 ms | 1.670 ms | 164.9 us | 266.341x | ok `e3cedb8421cd962f` |
| `fp.p7.is_irreducible.d64` | 281.532 ms | 21.242 ms | 1.689 ms | 1.618 ms | 130.4 us | 166.715x | 1.276 ms | 1.227 ms | 99.9 us | 220.584x | ok `9cb1bec87bc433d0` |
| `fp.p7.factor.d62` | 11.204 ms | 694.6 us | 293.1 us | 273.1 us | 36.3 us | 38.227x | 288.9 us | 277.6 us | 25.6 us | 38.783x | ok `3425aecb393956ba` |
| `fp.p101.elem_muladd.n400000` | 40.127 ms | 1.872 ms | 1.724 ms | 1.705 ms | 60.3 us | 23.274x | 1.725 ms | 1.676 ms | 120.4 us | 23.263x | ok `e51aa870ce6269eb` |
| `fp.p101.elem_inv.n100000` | 14.209 ms | 943.1 us | 683.3 us | 663.6 us | 40.7 us | 20.794x | 695.6 us | 663.6 us | 61.3 us | 20.425x | ok `69d137d626075364` |
| `fp.p101.elem_pow.n100000` | 17.852 ms | 1.103 ms | 2.176 ms | 2.017 ms | 180.4 us | 8.203x | 1.991 ms | 1.930 ms | 248.4 us | 8.966x | ok `09fec44561c97afa` |
| `fp.p101.mul_school.d64` | 9.588 ms | 865.8 us | 78.9 us | 74.4 us | 6.3 us | 121.521x | 76.6 us | 73.9 us | 6.9 us | 125.170x | ok `5a0fa852855ea8b0` |
| `fp.p101.mul_school.d512` | 29.713 ms | 2.531 ms | 201.5 us | 191.6 us | 12.6 us | 147.457x | 198.8 us | 190.0 us | 24.1 us | 149.460x | ok `8289633c899317a1` |
| `fp.p101.mul_karatsuba.d512` | 12.921 ms | 1.176 ms | 130.9 us | 123.5 us | 11.7 us | 98.711x | 129.3 us | 120.2 us | 16.5 us | 99.932x | ok `fea453dad54c0188` |
| `fp.p101.mul_karatsuba.d4096` | 368.458 ms | 24.423 ms | 1.460 ms | 1.154 ms | 193.1 us | 252.334x | 1.491 ms | 1.197 ms | 288.0 us | 247.088x | ok `82f83f0b7fe7713b` |
| `fp.p101.mul_auto.d2048` | 118.133 ms | 12.027 ms | 763.8 us | 565.7 us | 165.0 us | 154.655x | 109.5 us | 99.6 us | 13.1 us | 1079.337x | ok `3631fff5520b84d9` |
| `fp.p101.divmod.d2000_1000` | 112.089 ms | 10.178 ms | 804.6 us | 749.2 us | 95.1 us | 139.310x | 206.2 us | 184.5 us | 25.3 us | 543.463x | ok `ce32beed1dda91cd` |
| `fp.p101.gcd.d500` | 31.625 ms | 3.297 ms | 1.129 ms | 1.050 ms | 104.0 us | 28.017x | 1.239 ms | 1.075 ms | 191.4 us | 25.514x | ok `5cad0c7c5f168edb` |
| `fp.p101.gcd_ext.d400` | 63.657 ms | 5.662 ms | 2.364 ms | 2.220 ms | 183.4 us | 26.924x | 1.336 ms | 1.219 ms | 198.5 us | 47.655x | ok `4fe0534b04a103a1` |
| `fp.p101.pow_mod.d128_k20` | 627.358 ms | 38.646 ms | 5.619 ms | 5.316 ms | 319.1 us | 111.643x | 2.117 ms | 1.940 ms | 172.0 us | 296.413x | ok `bc5334d07b025488` |
| `fp.p101.is_irreducible.d48` | 353.192 ms | 96.214 ms | 1.584 ms | 1.500 ms | 92.4 us | 222.975x | 876.0 us | 806.0 us | 102.1 us | 403.187x | ok `5146dc1467eee05b` |
| `fp.p101.factor.d100` | 85.523 ms | 3.347 ms | 1.525 ms | 1.439 ms | 175.3 us | 56.066x | 1.333 ms | 1.228 ms | 197.2 us | 64.173x | ok `48a7604aa4ae63ee` |
| `fp.ntt.elem_muladd.n400000` | 44.694 ms | 3.360 ms | 2.550 ms | 2.469 ms | 102.0 us | 17.524x | 2.522 ms | 2.452 ms | 199.7 us | 17.725x | ok `682b833fe4edd17b` |
| `fp.ntt.elem_inv.n100000` | 17.360 ms | 1.058 ms | 3.127 ms | 2.981 ms | 217.6 us | 5.552x | 3.024 ms | 2.932 ms | 216.7 us | 5.741x | ok `d275b741e6b3050e` |
| `fp.ntt.elem_pow.n50000` | 14.560 ms | 3.019 ms | 6.878 ms | 6.683 ms | 453.4 us | 2.117x | 6.783 ms | 6.670 ms | 494.4 us | 2.146x | ok `c27b88ab2935c700` |
| `fp.ntt.mul_school.d64` | 9.979 ms | 882.7 us | 79.4 us | 76.5 us | 6.2 us | 125.685x | 78.8 us | 75.8 us | 8.2 us | 126.641x | ok `67c7e1645a9cc6ae` |
| `fp.ntt.mul_school.d512` | 30.477 ms | 6.101 ms | 201.4 us | 193.6 us | 18.6 us | 151.327x | 221.3 us | 194.3 us | 26.1 us | 137.688x | ok `e296393a6a637046` |
| `fp.ntt.mul_karatsuba.d512` | 12.605 ms | 1.218 ms | 134.0 us | 120.8 us | 13.7 us | 94.067x | 134.8 us | 121.6 us | 15.5 us | 93.509x | ok `964e15063fb16df4` |
| `fp.ntt.mul_karatsuba.d4096` | 366.833 ms | 15.392 ms | 1.720 ms | 1.444 ms | 332.2 us | 213.312x | 1.687 ms | 1.381 ms | 213.7 us | 217.512x | ok `d656b678f0c2ca16` |
| `fp.ntt.mul_fft.d4096` | 83.670 ms | 20.972 ms | 1.732 ms | 1.465 ms | 211.9 us | 48.317x | 1.705 ms | 1.447 ms | 253.4 us | 49.075x | ok `f4c81b85220fb683` |
| `fp.ntt.mul_fft.d32768` | 1.175 s | 284.020 ms | 16.010 ms | 14.797 ms | 1.074 ms | 73.410x | 15.668 ms | 14.252 ms | 2.009 ms | 75.015x | ok `8228e76950faded3` |
| `fp.ntt.mul_auto.d2048` | 36.461 ms | 4.412 ms | 761.4 us | 685.8 us | 97.4 us | 47.886x | 817.3 us | 694.0 us | 114.2 us | 44.611x | ok `02e67c888996018d` |
| `fp.ntt.divmod.d2000_1000` | 123.448 ms | 7.790 ms | 799.2 us | 767.6 us | 33.8 us | 154.455x | 551.4 us | 508.6 us | 63.6 us | 223.881x | ok `7adfbb2614325099` |
| `fp.ntt.gcd.d500` | 33.346 ms | 2.711 ms | 1.195 ms | 1.122 ms | 60.0 us | 27.917x | 1.540 ms | 1.363 ms | 233.9 us | 21.652x | ok `678962135e7a2447` |
| `fp.ntt.gcd_ext.d400` | 485.236 ms | 69.305 ms | 2.482 ms | 2.329 ms | 225.1 us | 195.486x | 2.580 ms | 2.314 ms | 411.6 us | 188.091x | ok `f9dc4e400a7954ee` |
| `fp.ntt.pow_mod.d128_k5` | 475.165 ms | 26.343 ms | 5.239 ms | 4.993 ms | 419.1 us | 90.692x | 4.176 ms | 3.845 ms | 602.1 us | 113.776x | ok `0eb2514aceb37261` |
| `fp.ntt.is_irreducible.d32` | 481.551 ms | 17.970 ms | 1.571 ms | 1.466 ms | 260.9 us | 306.427x | 367.8 us | 336.8 us | 44.9 us | 1309.096x | ok `1f2827a0d8aeb65e` |
| `fp.ntt.factor.d40` | 67.161 ms | 8.262 ms | 1.167 ms | 1.081 ms | 110.5 us | 57.565x | 1.381 ms | 1.245 ms | 213.3 us | 48.646x | ok `297c76c423fe55aa` |
| `fp.m61.elem_muladd.n400000` | 40.601 ms | 3.053 ms | 2.260 ms | 2.164 ms | 78.6 us | 17.965x | 2.311 ms | 2.173 ms | 359.6 us | 17.572x | ok `b8e325f73ad90582` |
| `fp.m61.elem_inv.n100000` | 20.870 ms | 669.2 us | 8.426 ms | 8.310 ms | 152.1 us | 2.477x | 9.308 ms | 8.371 ms | 1.070 ms | 2.242x | ok `41e7fe38de6fdcc1` |
| `fp.m61.elem_pow.n20000` | 8.240 ms | 269.6 us | 6.189 ms | 6.066 ms | 118.2 us | 1.332x | 6.432 ms | 6.047 ms | 646.2 us | 1.281x | ok `a3680ece7fc07932` |
| `fp.m61.mul_school.d64` | 9.500 ms | 328.1 us | 84.4 us | 79.2 us | 7.7 us | 112.555x | 85.1 us | 78.7 us | 11.7 us | 111.629x | ok `bed6e9e86907de79` |
| `fp.m61.mul_school.d512` | 30.037 ms | 1.856 ms | 219.7 us | 198.6 us | 25.1 us | 136.749x | 216.9 us | 199.8 us | 35.4 us | 138.483x | ok `43c3884710f7d468` |
| `fp.m61.mul_karatsuba.d512` | 12.215 ms | 239.9 us | 138.0 us | 126.0 us | 12.2 us | 88.512x | 144.1 us | 132.9 us | 30.9 us | 84.765x | ok `a2dcbeaa6574fcb8` |
| `fp.m61.mul_karatsuba.d4096` | 365.817 ms | 30.208 ms | 1.678 ms | 1.284 ms | 182.0 us | 217.956x | 1.666 ms | 1.207 ms | 340.7 us | 219.591x | ok `37d5deaa1cfeb600` |
| `fp.m61.mul_auto.d2048` | 127.410 ms | 15.345 ms | 920.1 us | 703.6 us | 115.1 us | 138.466x | 595.4 us | 527.8 us | 343.3 us | 213.990x | ok `73174ffa4e90d6fa` |
| `fp.m61.divmod.d2000_1000` | 121.389 ms | 11.191 ms | 862.4 us | 792.4 us | 87.6 us | 140.758x | 857.3 us | 787.4 us | 124.1 us | 141.595x | ok `a8710fc5aaa24d09` |
| `fp.m61.gcd.d500` | 33.352 ms | 2.558 ms | 1.040 ms | 975.7 us | 86.9 us | 32.057x | 1.305 ms | 1.192 ms | 214.7 us | 25.553x | ok `fda394dde025c1a2` |
| `fp.m61.gcd_ext.d400` | 64.561 ms | 6.991 ms | 2.196 ms | 2.053 ms | 263.1 us | 29.393x | 2.420 ms | 2.189 ms | 537.2 us | 26.681x | ok `be3d81e8acc44b48` |
| `fp.m61.pow_mod.d64_k2` | 158.080 ms | 24.562 ms | 1.565 ms | 1.456 ms | 166.8 us | 100.997x | 1.607 ms | 1.470 ms | 216.6 us | 98.394x | ok `394fb320787700fe` |
| `fp.m61.is_irreducible.d24` | 502.138 ms | 25.194 ms | 2.520 ms | 2.341 ms | 174.0 us | 199.293x | 335.4 us | 311.5 us | 43.4 us | 1497.131x | ok `9047cf1d3fc29abf` |
| `fp.m61.factor.d26` | 54.745 ms | 3.592 ms | 1.243 ms | 1.155 ms | 99.2 us | 44.051x | 1.334 ms | 1.194 ms | 165.9 us | 41.029x | ok `931b8d0b8b91b0b4` |
| `fp.m127.elem_muladd.n200000` | 24.390 ms | 2.497 ms | 24.389 ms | 22.923 ms | 883.1 us | 1.000x | 24.916 ms | 22.393 ms | 3.517 ms | 0.979x | ok `16fab7b139aa9aab` |
| `fp.m127.elem_inv.n50000` | 33.188 ms | 2.750 ms | 31.465 ms | 29.290 ms | 1.330 ms | 1.055x | 32.190 ms | 29.439 ms | 4.844 ms | 1.031x | ok `c7d67d6018036bf1` |
| `fp.m127.elem_pow.n10000` | 24.339 ms | 1.420 ms | 24.565 ms | 22.555 ms | 1.740 ms | 0.991x | 24.647 ms | 22.347 ms | 2.678 ms | 0.988x | ok `308896e55b40c815` |
| `fp.m127.mul_school.d64` | 10.682 ms | 395.7 us | 1.765 ms | 1.576 ms | 180.6 us | 6.053x | 1.711 ms | 1.577 ms | 178.9 us | 6.244x | ok `a9664408b0fb88e7` |
| `fp.m127.mul_school.d512` | 33.126 ms | 1.641 ms | 4.063 ms | 3.697 ms | 421.9 us | 8.153x | 4.221 ms | 3.760 ms | 706.3 us | 7.848x | ok `98eca357802203c0` |
| `fp.m127.mul_karatsuba.d512` | 14.804 ms | 1.165 ms | 1.873 ms | 1.534 ms | 210.4 us | 7.904x | 1.881 ms | 1.203 ms | 304.9 us | 7.872x | ok `f773c153fde34172` |
| `fp.m127.mul_karatsuba.d2048` | 143.352 ms | 13.446 ms | 10.761 ms | 9.706 ms | 800.5 us | 13.321x | 10.793 ms | 8.726 ms | 1.109 ms | 13.281x | ok `695e0d0f13c14e87` |
| `fp.m127.mul_auto.d2048` | 144.144 ms | 18.820 ms | 10.659 ms | 9.963 ms | 588.0 us | 13.524x | 1.936 ms | 1.622 ms | 507.1 us | 74.445x | ok `1bc667b834b9e559` |
| `fp.m127.divmod.d2000_1000` | 126.455 ms | 20.209 ms | 15.288 ms | 14.153 ms | 1.083 ms | 8.272x | 4.508 ms | 3.977 ms | 684.5 us | 28.051x | ok `a8babdcd119eb404` |
| `fp.m127.gcd.d400` | 22.352 ms | 2.776 ms | 14.993 ms | 13.928 ms | 754.6 us | 1.491x | 13.266 ms | 11.408 ms | 1.707 ms | 1.685x | ok `1b15fac38c54accc` |
| `fp.m127.gcd_ext.d300` | 39.933 ms | 4.309 ms | 24.954 ms | 22.898 ms | 1.974 ms | 1.600x | 13.643 ms | 12.214 ms | 1.729 ms | 2.927x | ok `45a1b9001453c6ef` |
| `fp.m127.pow_mod.d64_k1` | 230.995 ms | 20.814 ms | 47.703 ms | 43.808 ms | 8.195 ms | 4.842x | 28.520 ms | 26.292 ms | 3.359 ms | 8.099x | ok `9e9e535b3c4e4e7d` |
| `fp.m127.is_irreducible.d16` | 337.766 ms | 29.594 ms | 40.535 ms | 37.172 ms | 5.466 ms | 8.333x | 5.466 ms | 4.974 ms | 328.9 us | 61.796x | ok `2cfb69ee2746aa9d` |
| `fp.m127.factor.d19` | 85.242 ms | 8.945 ms | 37.558 ms | 35.425 ms | 11.121 ms | 2.270x | 39.929 ms | 35.538 ms | 4.616 ms | 2.135x | ok `7b0defab2417737c` |
| `fp.p256.elem_muladd.n200000` | 32.723 ms | 3.643 ms | 33.144 ms | 30.601 ms | 5.754 ms | 0.987x | 34.409 ms | 31.364 ms | 1.985 ms | 0.951x | ok `cc3cab838e9b19f1` |
| `fp.p256.elem_inv.n50000` | 68.130 ms | 6.294 ms | 69.600 ms | 63.773 ms | 10.446 ms | 0.979x | 69.267 ms | 65.136 ms | 4.389 ms | 0.984x | ok `89f920a8eba0c743` |
| `fp.p256.elem_pow.n2000` | 19.872 ms | 2.906 ms | 20.007 ms | 18.954 ms | 3.249 ms | 0.993x | 19.673 ms | 18.902 ms | 1.685 ms | 1.010x | ok `40472b256a3c877d` |
| `fp.p256.mul_school.d64` | 14.156 ms | 1.264 ms | 2.539 ms | 2.422 ms | 640.1 us | 5.576x | 2.561 ms | 2.407 ms | 322.8 us | 5.527x | ok `f5d3d87e3ed0977e` |
| `fp.p256.mul_school.d512` | 44.060 ms | 4.733 ms | 6.618 ms | 6.169 ms | 3.266 ms | 6.658x | 6.939 ms | 6.058 ms | 429.5 us | 6.350x | ok `6fd19bd7b3ceab69` |
| `fp.p256.mul_karatsuba.d512` | 18.362 ms | 1.727 ms | 2.094 ms | 1.751 ms | 240.4 us | 8.769x | 2.105 ms | 1.488 ms | 300.9 us | 8.724x | ok `4b2f43bf0f878408` |
| `fp.p256.mul_karatsuba.d2048` | 170.875 ms | 23.492 ms | 13.506 ms | 11.830 ms | 1.275 ms | 12.651x | 13.812 ms | 10.502 ms | 1.470 ms | 12.372x | ok `b4b89a80c2a3aaa9` |
| `fp.p256.mul_fft.d4096` | 103.964 ms | 22.852 ms | 18.072 ms | 15.280 ms | 4.706 ms | 5.753x | 19.003 ms | 15.300 ms | 4.090 ms | 5.471x | ok `60ed510564513dfb` |
| `fp.p256.mul_fft.d8192` | 245.531 ms | 60.048 ms | 38.900 ms | 31.756 ms | 9.110 ms | 6.312x | 36.742 ms | 31.956 ms | 10.599 ms | 6.683x | ok `47812935df4de5fe` |
| `fp.p256.mul_auto.d2048` | 45.487 ms | 4.753 ms | 8.747 ms | 7.430 ms | 917.6 us | 5.200x | 4.082 ms | 3.217 ms | 919.3 us | 11.144x | ok `3d89402081500506` |
| `fp.p256.divmod.d1000_500` | 42.534 ms | 6.860 ms | 6.147 ms | 5.914 ms | 531.7 us | 6.920x | 3.437 ms | 3.134 ms | 376.4 us | 12.375x | ok `0bebce51d98eb3f4` |
| `fp.p256.gcd.d300` | 17.190 ms | 1.427 ms | 10.902 ms | 10.348 ms | 2.409 ms | 1.577x | 11.106 ms | 10.346 ms | 1.272 ms | 1.548x | ok `f69e4050cb1a5861` |
| `fp.p256.gcd_ext.d200` | 126.287 ms | 11.863 ms | 13.883 ms | 13.160 ms | 3.378 ms | 9.096x | 12.638 ms | 11.880 ms | 604.1 us | 9.993x | ok `1a9df499c28b9f49` |
| `fp.p256.pow_mod.d64_k1` | 401.416 ms | 43.627 ms | 79.660 ms | 77.234 ms | 5.504 ms | 5.039x | 71.379 ms | 66.040 ms | 6.654 ms | 5.624x | ok `dbe2b595567fbb9c` |
| `fp.p256.is_irreducible.d12` | 349.066 ms | 36.663 ms | 41.679 ms | 39.047 ms | 3.297 ms | 8.375x | 7.259 ms | 6.673 ms | 709.4 us | 48.084x | ok `7b839a40c98b0ef0` |
| `fp.p256.factor.d12` | 36.630 ms | 4.283 ms | 17.777 ms | 17.230 ms | 462.9 us | 2.061x | 17.662 ms | 16.748 ms | 1.222 ms | 2.074x | ok `c33e3512ed388a7e` |
| `fq.q3e30.elem_muladd.n4000` | 208.788 ms | 25.666 ms | 5.900 ms | 5.425 ms | 557.4 us | 35.387x | 5.857 ms | 5.505 ms | 500.7 us | 35.646x | ok `16edfcc06a08028a` |
| `fq.q3e30.elem_inv.n100` | 336.271 ms | 28.802 ms | 754.6 us | 712.7 us | 67.1 us | 445.628x | 784.9 us | 717.2 us | 79.0 us | 428.398x | ok `a8167c01d1cfc5dd` |
| `fq.q3e30.mul_school.d32` | 54.268 ms | 5.858 ms | 734.1 us | 708.1 us | 86.7 us | 73.924x | 767.8 us | 707.7 us | 72.4 us | 70.680x | ok `ed0ca8eb7d325de3` |
| `fq.q3e30.mul_karatsuba.d128` | 420.173 ms | 50.658 ms | 2.335 ms | 2.023 ms | 179.3 us | 179.969x | 2.265 ms | 1.520 ms | 369.8 us | 185.482x | ok `2464b901bd242c7d` |
| `fq.q3e30.divmod.d64_32` | 56.044 ms | 8.040 ms | 776.1 us | 726.5 us | 68.7 us | 72.208x | 783.7 us | 722.9 us | 104.4 us | 71.513x | ok `a2b35f414a139741` |
| `fq.q3e30.factor.d5` | 422.848 ms | 30.597 ms | 6.980 ms | 6.644 ms | 644.3 us | 60.584x | 3.864 ms | 3.585 ms | 391.0 us | 109.434x | ok `417a5b753c7a2a07` |
| `fq.q2e8.elem_muladd.n50000` | 211.448 ms | 20.546 ms | 9.411 ms | 8.609 ms | 889.2 us | 22.467x | 9.813 ms | 8.559 ms | 1.159 ms | 21.547x | ok `c0ffe16ee54f128c` |
| `fq.q2e8.elem_inv.n5000` | 240.700 ms | 20.415 ms | 2.511 ms | 2.283 ms | 184.5 us | 95.849x | 2.511 ms | 2.240 ms | 285.4 us | 95.875x | ok `81d71b5687e3e735` |
| `fq.q2e8.mul_school.d128` | 65.324 ms | 3.687 ms | 618.5 us | 571.0 us | 47.4 us | 105.617x | 619.2 us | 567.5 us | 105.2 us | 105.497x | ok `c056d9d992271db9` |
| `fq.q2e8.mul_karatsuba.d512` | 351.070 ms | 65.204 ms | 2.301 ms | 1.873 ms | 215.1 us | 152.596x | 2.377 ms | 1.643 ms | 373.7 us | 147.670x | ok `de1089fdce26d02c` |
| `fq.q2e8.is_irreducible.d16` | 385.210 ms | 51.728 ms | 2.558 ms | 2.343 ms | 257.5 us | 150.570x | 1.563 ms | 1.431 ms | 269.2 us | 246.401x | ok `9478d590a29a9cbf` |
| `fq.q2e8.factor.d28` | 94.306 ms | 4.862 ms | 2.536 ms | 2.413 ms | 250.5 us | 37.186x | 2.737 ms | 2.535 ms | 324.4 us | 34.452x | ok `d2fc95c10ececedd` |
| `fq.q7e4.elem_muladd.n50000` | 96.061 ms | 6.360 ms | 6.266 ms | 5.884 ms | 685.3 us | 15.330x | 6.288 ms | 5.911 ms | 714.2 us | 15.276x | ok `a6cf18a43d40c5c0` |
| `fq.q7e4.elem_inv.n5000` | 187.750 ms | 13.076 ms | 1.847 ms | 1.738 ms | 173.4 us | 101.632x | 1.986 ms | 1.786 ms | 223.5 us | 94.523x | ok `8664ac6d09409f94` |
| `fq.q7e4.mul_school.d128` | 37.406 ms | 4.287 ms | 294.1 us | 268.7 us | 31.1 us | 127.187x | 304.6 us | 267.5 us | 58.2 us | 122.803x | ok `71d186560b1d8153` |
| `fq.q7e4.mul_karatsuba.d512` | 207.789 ms | 20.506 ms | 1.757 ms | 1.495 ms | 209.5 us | 118.230x | 1.807 ms | 1.214 ms | 532.3 us | 114.972x | ok `6a2e9b751ab48e99` |
| `fq.q7e4.gcd.d200` | 103.416 ms | 8.105 ms | 3.586 ms | 3.195 ms | 506.6 us | 28.836x | 3.612 ms | 3.189 ms | 826.8 us | 28.632x | ok `8cf19ac90e512feb` |
| `fq.q7e4.factor.d20` | 59.220 ms | 3.752 ms | 1.482 ms | 1.369 ms | 109.8 us | 39.965x | 1.318 ms | 1.206 ms | 236.4 us | 44.928x | ok `ce23318b1b5fd48e` |
| `fq.m61e2.elem_muladd.n100000` | 86.492 ms | 7.013 ms | 9.309 ms | 8.639 ms | 1.143 ms | 9.291x | 9.449 ms | 8.665 ms | 2.098 ms | 9.154x | ok `1814347c960518d8` |
| `fq.m61e2.elem_inv.n1000` | 145.270 ms | 7.758 ms | 402.1 us | 385.3 us | 32.6 us | 361.277x | 411.8 us | 385.0 us | 38.2 us | 352.767x | ok `45326d8095032432` |
| `fq.m61e2.mul_school.d128` | 14.282 ms | 868.1 us | 121.8 us | 113.2 us | 11.0 us | 117.261x | 127.6 us | 112.2 us | 12.0 us | 111.931x | ok `2cb38bd639886e3c` |
| `fq.m61e2.mul_karatsuba.d512` | 77.449 ms | 3.365 ms | 1.434 ms | 1.105 ms | 240.9 us | 54.003x | 1.389 ms | 907.1 us | 275.1 us | 55.743x | ok `f138ad86c0a7a838` |
| `fq.m61e2.mul_fft.d1024` | 82.440 ms | 5.655 ms | 15.645 ms | 13.679 ms | 1.774 ms | 5.269x | 15.416 ms | 12.089 ms | 1.711 ms | 5.348x | ok `a45642d67f5a9912` |
| `fq.m61e2.mul_fft.d4096` | 463.403 ms | 72.292 ms | 42.433 ms | 39.799 ms | 3.327 ms | 10.921x | 43.174 ms | 41.416 ms | 4.060 ms | 10.733x | ok `51844e8f95a9c0ab` |
| `fq.m61e2.mul_auto.d2048` | 194.194 ms | 13.809 ms | 26.664 ms | 23.212 ms | 2.988 ms | 7.283x | 2.977 ms | 2.590 ms | 494.9 us | 65.230x | ok `76e66ad22fb3fa2e` |
| `fq.m61e2.divmod.d500_250` | 54.523 ms | 3.025 ms | 479.8 us | 444.2 us | 54.2 us | 113.638x | 480.8 us | 448.3 us | 28.1 us | 113.401x | ok `53509a111940ee2a` |
| `fq.m61e2.gcd.d150` | 39.513 ms | 2.027 ms | 1.391 ms | 1.281 ms | 128.4 us | 28.410x | 1.403 ms | 1.302 ms | 113.1 us | 28.169x | ok `bd8c3c98cc56df84` |
| `fq.m61e2.factor.d7` | 157.637 ms | 2.701 ms | 2.938 ms | 2.750 ms | 277.3 us | 53.663x | 3.117 ms | 2.934 ms | 213.1 us | 50.568x | ok `d61430d11ee187c3` |
| `fq.build.q3e30` | 21.771 ms | 1.995 ms | 246.1 us | 228.5 us | 31.0 us | 88.466x | 287.2 us | 257.7 us | 57.2 us | 75.819x | ok `4973a3e567e2fc74` |
| `fq.build.q7e20` | 18.822 ms | 1.325 ms | 213.6 us | 196.0 us | 26.9 us | 88.095x | 216.4 us | 192.3 us | 29.4 us | 86.976x | ok `9c377cd554aa108e` |
| `fast.ntt.dft.k10` | 15.402 ms | 1.970 ms | 1.272 ms | 1.176 ms | 142.2 us | 12.108x | 1.359 ms | 1.225 ms | 255.5 us | 11.331x | ok `b61f454e63ae6d8d` |
| `fast.ntt.dft.k11` | 16.386 ms | 721.5 us | 1.292 ms | 1.180 ms | 177.9 us | 12.686x | 1.305 ms | 1.240 ms | 342.4 us | 12.552x | ok `311d4fed0f42034c` |
| `fast.ntt.dft.k12` | 17.847 ms | 2.143 ms | 1.261 ms | 1.185 ms | 111.5 us | 14.149x | 1.368 ms | 1.232 ms | 177.0 us | 13.046x | ok `70691e178cfc0c9d` |
| `fast.ntt.idft.k12` | 18.494 ms | 1.801 ms | 1.287 ms | 1.196 ms | 136.2 us | 14.371x | 1.340 ms | 1.247 ms | 155.0 us | 13.799x | ok `e0677b59c4e33615` |
| `fast.ntt.toeplitz_lower_inverse.n1024` | 25.911 ms | 3.418 ms | 610.4 us | 584.2 us | 41.9 us | 42.450x | 383.0 us | 342.7 us | 64.5 us | 67.654x | ok `7749eb8ff8cac3d0` |
| `fast.ntt.toeplitz_lower_inverse.n8192` | 291.507 ms | 28.639 ms | 6.474 ms | 6.136 ms | 639.4 us | 45.025x | 6.461 ms | 5.535 ms | 1.580 ms | 45.120x | ok `0d31554b19d8b996` |
| `fast.ntt.toeplitz_lower_vec.n4096` | 34.947 ms | 5.522 ms | 812.2 us | 708.0 us | 82.8 us | 43.028x | 943.1 us | 746.7 us | 246.4 us | 37.053x | ok `86e36ad2a1eacb05` |
| `fast.ntt.toeplitz_divmod.d2000_1000` | 61.528 ms | 7.648 ms | 2.922 ms | 2.749 ms | 269.3 us | 21.056x | 2.866 ms | 2.741 ms | 417.8 us | 21.466x | ok `b138cdcceebcc32e` |
| `fast.p7.toeplitz_lower_inverse.n1024` | 72.812 ms | 9.044 ms | 1.005 ms | 934.1 us | 112.8 us | 72.450x | 109.2 us | 96.6 us | 19.9 us | 666.781x | ok `5df08c92a814155a` |
| `fast.p7.toeplitz_divmod.d2000_1000` | 80.302 ms | 9.742 ms | 1.630 ms | 1.519 ms | 186.5 us | 49.256x | 1.382 ms | 1.246 ms | 144.3 us | 58.114x | ok `433d13716b2f95e8` |
