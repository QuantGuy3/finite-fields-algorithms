# Comparativa de benchmarks alcp

Base: **f4-ref**. Speedup = mediana(base) / mediana(variante); > 1 significa más rápido.
IQR = rango intercuartílico de las muestras. Las medias geométricas excluyen escenarios con error.

- **f4-ref**: modo interleave, 2026-11-07T17:10:12+01:00 (nota: fase 4 nivel A: original vs fase 3 vs fase 4)  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3`  
  include `~\AppData\Local\Temp\claude\204af1c0d38157dc\2e18e96a-361a-47be-ad99-f0341eb0cc55\scratchpad\inc-ref\include` @ ?, sha256 9b306a5d8af5  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: sí
- **f4-f3**: modo interleave, 2026-11-07T17:10:26+01:00 (nota: fase 4 nivel A: original vs fase 3 vs fase 4)  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3`  
  include `~\AppData\Local\Temp\claude\204af1c0d38157dc\2e18e96a-361a-47be-ad99-f0341eb0cc55\scratchpad\inc-f3\include` @ ?, sha256 af0c61d3df27  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: sí
- **f4-f4**: modo interleave, 2026-11-07T17:10:40+01:00 (nota: fase 4 nivel A: original vs fase 3 vs fase 4)  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -DALCP_RELINK=1`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-65-g0aaa877 (con cambios sin commit), sha256 8eb181a79710  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: sí

Checksums: todos los escenarios coinciden con la base.

## Media geométrica del speedup

| ámbito | f4-f3 | f4-f4 |
|---|---|---|
| **global** | 23.267x (n=128) | 27.736x (n=128) |
| grupo `fast` | 26.539x (n=10) | 27.547x (n=10) |
| grupo `fp.m127` | 3.223x (n=14) | 3.492x (n=14) |
| grupo `fp.m61` | 42.265x (n=14) | 43.879x (n=14) |
| grupo `fp.ntt` | 50.336x (n=16) | 53.081x (n=16) |
| grupo `fp.p101` | 56.525x (n=14) | 58.677x (n=14) |
| grupo `fp.p256` | 3.508x (n=16) | 3.883x (n=16) |
| grupo `fp.p7` | 57.613x (n=14) | 61.423x (n=14) |
| grupo `fq.build` | 44.700x (n=2) | 88.853x (n=2) |
| grupo `fq.m61e2` | 17.097x (n=10) | 33.780x (n=10) |
| grupo `fq.q2e8` | 49.075x (n=6) | 79.703x (n=6) |
| grupo `fq.q3e30` | 56.057x (n=6) | 96.077x (n=6) |
| grupo `fq.q7e4` | 37.262x (n=6) | 54.197x (n=6) |
| operación `build` | 44.700x (n=2) | 88.853x (n=2) |
| operación `dft` | 11.836x (n=3) | 12.269x (n=3) |
| operación `divmod` | 54.517x (n=8) | 56.710x (n=8) |
| operación `elem_inv` | 7.637x (n=10) | 20.211x (n=10) |
| operación `elem_muladd` | 10.123x (n=10) | 10.294x (n=10) |
| operación `elem_pow` | 2.538x (n=6) | 2.520x (n=6) |
| operación `factor` | 18.177x (n=10) | 24.124x (n=10) |
| operación `gcd` | 11.267x (n=8) | 13.203x (n=8) |
| operación `gcd_ext` | 18.075x (n=6) | 18.121x (n=6) |
| operación `idft` | 14.756x (n=1) | 15.118x (n=1) |
| operación `is_irreducible` | 27.784x (n=7) | 74.707x (n=7) |
| operación `mul_auto` | 31.340x (n=7) | 31.444x (n=7) |
| operación `mul_fft` | 13.590x (n=6) | 13.533x (n=6) |
| operación `mul_karatsuba` | 67.902x (n=16) | 65.992x (n=16) |
| operación `mul_school` | 54.552x (n=16) | 54.333x (n=16) |
| operación `pow_mod` | 33.770x (n=6) | 34.079x (n=6) |
| operación `toeplitz_divmod` | 33.549x (n=2) | 35.127x (n=2) |
| operación `toeplitz_lower_inverse` | 52.905x (n=3) | 54.053x (n=3) |
| operación `toeplitz_lower_vec` | 42.514x (n=1) | 46.239x (n=1) |

## Detalle por escenario

| escenario | base mediana | base IQR | f4-f3 mediana | f4-f3 mín | f4-f3 IQR | speedup | f4-f4 mediana | f4-f4 mín | f4-f4 IQR | speedup | checksum |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `fp.p7.elem_muladd.n400000` | 45.395 ms | 3.836 ms | 2.168 ms | 1.704 ms | 638.0 us | 20.934x | 2.026 ms | 1.706 ms | 406.7 us | 22.407x | ok `cf4797a827e7f9ec` |
| `fp.p7.elem_inv.n100000` | 13.969 ms | 1.939 ms | 392.5 us | 336.0 us | 102.6 us | 35.590x | 394.8 us | 334.1 us | 59.8 us | 35.387x | ok `7beb5a3238b919cb` |
| `fp.p7.elem_pow.n100000` | 16.142 ms | 1.810 ms | 1.086 ms | 939.6 us | 266.8 us | 14.860x | 1.143 ms | 975.9 us | 236.2 us | 14.117x | ok `a7f8c429d0d81b73` |
| `fp.p7.mul_school.d64` | 8.251 ms | 1.020 ms | 81.3 us | 67.9 us | 19.9 us | 101.423x | 87.1 us | 69.8 us | 18.6 us | 94.727x | ok `559357b5bc3c57c2` |
| `fp.p7.mul_school.d512` | 26.481 ms | 3.554 ms | 204.5 us | 178.6 us | 35.3 us | 129.490x | 206.0 us | 165.5 us | 36.2 us | 128.547x | ok `5936b7480b151eb3` |
| `fp.p7.mul_karatsuba.d512` | 12.612 ms | 1.683 ms | 158.9 us | 133.6 us | 38.2 us | 79.344x | 160.2 us | 124.3 us | 35.8 us | 78.724x | ok `d72c59b658ac0379` |
| `fp.p7.mul_karatsuba.d4096` | 357.284 ms | 45.744 ms | 1.668 ms | 1.252 ms | 393.6 us | 214.225x | 1.730 ms | 1.291 ms | 385.3 us | 206.570x | ok `00f3dde93cbd25a4` |
| `fp.p7.mul_auto.d2048` | 113.514 ms | 13.933 ms | 979.2 us | 711.1 us | 285.7 us | 115.925x | 918.0 us | 679.6 us | 328.9 us | 123.653x | ok `2b256c747e0070c7` |
| `fp.p7.divmod.d2000_1000` | 97.192 ms | 11.218 ms | 874.5 us | 705.4 us | 167.3 us | 111.140x | 816.4 us | 676.7 us | 215.3 us | 119.049x | ok `958d383cda71b699` |
| `fp.p7.gcd.d500` | 26.748 ms | 3.383 ms | 1.070 ms | 900.1 us | 209.9 us | 24.993x | 998.5 us | 856.9 us | 207.1 us | 26.788x | ok `cdf18334cf9b13d7` |
| `fp.p7.gcd_ext.d400` | 57.563 ms | 7.121 ms | 2.349 ms | 1.979 ms | 482.3 us | 24.508x | 2.336 ms | 1.902 ms | 469.6 us | 24.646x | ok `984ac156886a3f82` |
| `fp.p7.pow_mod.d128_k48` | 521.049 ms | 102.070 ms | 5.649 ms | 5.045 ms | 353.7 us | 92.234x | 5.614 ms | 4.795 ms | 1.184 ms | 92.809x | ok `e3cedb8421cd962f` |
| `fp.p7.is_irreducible.d64` | 301.902 ms | 39.728 ms | 4.322 ms | 3.750 ms | 1.201 ms | 69.847x | 1.936 ms | 1.659 ms | 480.5 us | 155.917x | ok `9cb1bec87bc433d0` |
| `fp.p7.factor.d62` | 11.947 ms | 1.196 ms | 338.2 us | 284.4 us | 93.9 us | 35.325x | 342.2 us | 279.3 us | 99.3 us | 34.907x | ok `3425aecb393956ba` |
| `fp.p101.elem_muladd.n400000` | 41.795 ms | 4.800 ms | 1.911 ms | 1.716 ms | 415.0 us | 21.865x | 1.979 ms | 1.686 ms | 484.8 us | 21.114x | ok `e51aa870ce6269eb` |
| `fp.p101.elem_inv.n100000` | 14.941 ms | 1.760 ms | 799.9 us | 678.6 us | 167.8 us | 18.679x | 762.0 us | 667.3 us | 188.1 us | 19.606x | ok `69d137d626075364` |
| `fp.p101.elem_pow.n100000` | 19.379 ms | 7.218 ms | 2.477 ms | 2.181 ms | 411.7 us | 7.823x | 2.392 ms | 2.027 ms | 642.7 us | 8.103x | ok `09fec44561c97afa` |
| `fp.p101.mul_school.d64` | 9.877 ms | 1.863 ms | 87.9 us | 74.8 us | 14.3 us | 112.363x | 87.4 us | 74.5 us | 25.7 us | 113.006x | ok `5a0fa852855ea8b0` |
| `fp.p101.mul_school.d512` | 32.014 ms | 7.917 ms | 227.3 us | 194.8 us | 60.9 us | 140.814x | 225.3 us | 195.1 us | 61.6 us | 142.095x | ok `8289633c899317a1` |
| `fp.p101.mul_karatsuba.d512` | 13.449 ms | 2.079 ms | 143.6 us | 122.4 us | 44.0 us | 93.654x | 158.5 us | 121.8 us | 35.1 us | 84.850x | ok `fea453dad54c0188` |
| `fp.p101.mul_karatsuba.d4096` | 384.447 ms | 47.013 ms | 1.735 ms | 1.249 ms | 507.2 us | 221.571x | 1.755 ms | 1.393 ms | 396.4 us | 219.096x | ok `82f83f0b7fe7713b` |
| `fp.p101.mul_auto.d2048` | 127.492 ms | 9.010 ms | 1.007 ms | 727.1 us | 290.6 us | 126.568x | 1.072 ms | 774.5 us | 289.2 us | 118.896x | ok `3631fff5520b84d9` |
| `fp.p101.divmod.d2000_1000` | 118.105 ms | 11.124 ms | 1.034 ms | 842.1 us | 210.1 us | 114.266x | 1.105 ms | 910.1 us | 219.2 us | 106.902x | ok `ce32beed1dda91cd` |
| `fp.p101.gcd.d500` | 33.143 ms | 5.295 ms | 1.357 ms | 1.199 ms | 332.9 us | 24.431x | 1.406 ms | 1.112 ms | 238.5 us | 23.569x | ok `5cad0c7c5f168edb` |
| `fp.p101.gcd_ext.d400` | 66.147 ms | 12.673 ms | 2.779 ms | 2.474 ms | 501.0 us | 23.802x | 2.885 ms | 2.326 ms | 664.3 us | 22.924x | ok `4fe0534b04a103a1` |
| `fp.p101.pow_mod.d128_k20` | 632.875 ms | 49.858 ms | 6.123 ms | 5.670 ms | 1.190 ms | 103.357x | 6.562 ms | 5.443 ms | 1.603 ms | 96.442x | ok `bc5334d07b025488` |
| `fp.p101.is_irreducible.d48` | 360.535 ms | 44.702 ms | 4.814 ms | 4.080 ms | 1.005 ms | 74.901x | 1.965 ms | 1.518 ms | 492.4 us | 183.450x | ok `5146dc1467eee05b` |
| `fp.p101.factor.d100` | 89.374 ms | 11.389 ms | 1.795 ms | 1.503 ms | 280.5 us | 49.780x | 1.896 ms | 1.512 ms | 447.6 us | 47.146x | ok `48a7604aa4ae63ee` |
| `fp.ntt.elem_muladd.n400000` | 46.925 ms | 3.285 ms | 2.787 ms | 2.597 ms | 388.1 us | 16.839x | 2.864 ms | 2.531 ms | 712.8 us | 16.387x | ok `682b833fe4edd17b` |
| `fp.ntt.elem_inv.n100000` | 18.277 ms | 1.975 ms | 3.336 ms | 2.986 ms | 442.9 us | 5.478x | 3.538 ms | 2.970 ms | 681.4 us | 5.166x | ok `d275b741e6b3050e` |
| `fp.ntt.elem_pow.n50000` | 15.714 ms | 2.930 ms | 7.266 ms | 6.799 ms | 1.157 ms | 2.163x | 7.541 ms | 6.797 ms | 1.619 ms | 2.084x | ok `c27b88ab2935c700` |
| `fp.ntt.mul_school.d64` | 10.289 ms | 1.380 ms | 93.5 us | 76.9 us | 16.9 us | 109.989x | 91.6 us | 76.9 us | 25.7 us | 112.330x | ok `67c7e1645a9cc6ae` |
| `fp.ntt.mul_school.d512` | 32.427 ms | 2.804 ms | 235.7 us | 206.2 us | 32.0 us | 137.608x | 232.3 us | 195.7 us | 50.6 us | 139.593x | ok `e296393a6a637046` |
| `fp.ntt.mul_karatsuba.d512` | 13.873 ms | 1.087 ms | 151.9 us | 132.6 us | 43.0 us | 91.331x | 156.5 us | 124.9 us | 35.2 us | 88.647x | ok `964e15063fb16df4` |
| `fp.ntt.mul_karatsuba.d4096` | 385.661 ms | 39.526 ms | 1.946 ms | 1.572 ms | 382.6 us | 198.161x | 2.014 ms | 1.659 ms | 354.9 us | 191.500x | ok `d656b678f0c2ca16` |
| `fp.ntt.mul_fft.d4096` | 86.591 ms | 14.900 ms | 2.040 ms | 1.631 ms | 487.1 us | 42.453x | 2.103 ms | 1.579 ms | 450.9 us | 41.185x | ok `f4c81b85220fb683` |
| `fp.ntt.mul_fft.d32768` | 1.360 s | 162.145 ms | 17.586 ms | 15.808 ms | 3.424 ms | 77.322x | 18.911 ms | 15.772 ms | 5.257 ms | 71.904x | ok `8228e76950faded3` |
| `fp.ntt.mul_auto.d2048` | 37.191 ms | 7.881 ms | 931.1 us | 723.7 us | 269.2 us | 39.941x | 925.7 us | 727.5 us | 260.0 us | 40.176x | ok `02e67c888996018d` |
| `fp.ntt.divmod.d2000_1000` | 129.413 ms | 17.295 ms | 968.5 us | 812.2 us | 176.9 us | 133.615x | 919.2 us | 807.6 us | 91.5 us | 140.781x | ok `7adfbb2614325099` |
| `fp.ntt.gcd.d500` | 35.572 ms | 4.882 ms | 1.402 ms | 1.181 ms | 230.2 us | 25.376x | 1.398 ms | 1.160 ms | 140.9 us | 25.440x | ok `678962135e7a2447` |
| `fp.ntt.gcd_ext.d400` | 499.782 ms | 66.808 ms | 2.874 ms | 2.470 ms | 516.3 us | 173.910x | 2.987 ms | 2.426 ms | 475.9 us | 167.302x | ok `f9dc4e400a7954ee` |
| `fp.ntt.pow_mod.d128_k5` | 491.457 ms | 62.188 ms | 6.100 ms | 5.198 ms | 688.2 us | 80.563x | 6.069 ms | 5.652 ms | 455.0 us | 80.985x | ok `0eb2514aceb37261` |
| `fp.ntt.is_irreducible.d32` | 493.736 ms | 66.223 ms | 6.066 ms | 5.026 ms | 1.424 ms | 81.399x | 1.805 ms | 1.634 ms | 192.0 us | 273.492x | ok `1f2827a0d8aeb65e` |
| `fp.ntt.factor.d40` | 76.367 ms | 10.751 ms | 1.246 ms | 1.117 ms | 388.8 us | 61.305x | 1.425 ms | 1.160 ms | 454.3 us | 53.578x | ok `297c76c423fe55aa` |
| `fp.m61.elem_muladd.n400000` | 44.279 ms | 5.203 ms | 2.628 ms | 2.257 ms | 593.1 us | 16.847x | 2.614 ms | 2.403 ms | 744.4 us | 16.942x | ok `b8e325f73ad90582` |
| `fp.m61.elem_inv.n100000` | 23.086 ms | 4.542 ms | 9.324 ms | 8.471 ms | 1.545 ms | 2.476x | 9.446 ms | 8.302 ms | 2.185 ms | 2.444x | ok `41e7fe38de6fdcc1` |
| `fp.m61.elem_pow.n20000` | 9.104 ms | 556.6 us | 6.637 ms | 6.087 ms | 1.038 ms | 1.372x | 7.147 ms | 6.431 ms | 1.673 ms | 1.274x | ok `a3680ece7fc07932` |
| `fp.m61.mul_school.d64` | 10.533 ms | 1.507 ms | 98.4 us | 79.7 us | 23.8 us | 107.038x | 110.3 us | 85.9 us | 18.5 us | 95.490x | ok `bed6e9e86907de79` |
| `fp.m61.mul_school.d512` | 32.791 ms | 4.693 ms | 239.1 us | 207.9 us | 68.7 us | 137.116x | 260.7 us | 216.7 us | 50.7 us | 125.782x | ok `43c3884710f7d468` |
| `fp.m61.mul_karatsuba.d512` | 13.404 ms | 3.069 ms | 155.0 us | 128.5 us | 50.5 us | 86.480x | 161.5 us | 136.9 us | 33.6 us | 82.999x | ok `a2dcbeaa6574fcb8` |
| `fp.m61.mul_karatsuba.d4096` | 389.649 ms | 82.561 ms | 1.925 ms | 1.431 ms | 413.8 us | 202.394x | 1.997 ms | 1.582 ms | 381.0 us | 195.102x | ok `37d5deaa1cfeb600` |
| `fp.m61.mul_auto.d2048` | 129.616 ms | 24.322 ms | 1.073 ms | 711.2 us | 350.5 us | 120.747x | 1.114 ms | 823.5 us | 278.9 us | 116.394x | ok `73174ffa4e90d6fa` |
| `fp.m61.divmod.d2000_1000` | 143.245 ms | 54.413 ms | 987.9 us | 813.8 us | 249.1 us | 145.000x | 1.025 ms | 911.2 us | 201.6 us | 139.697x | ok `a8710fc5aaa24d09` |
| `fp.m61.gcd.d500` | 35.640 ms | 12.277 ms | 1.206 ms | 1.009 ms | 261.4 us | 29.550x | 1.191 ms | 1.082 ms | 111.6 us | 29.917x | ok `fda394dde025c1a2` |
| `fp.m61.gcd_ext.d400` | 70.118 ms | 12.450 ms | 2.529 ms | 2.141 ms | 475.2 us | 27.726x | 2.468 ms | 2.232 ms | 343.5 us | 28.415x | ok `be3d81e8acc44b48` |
| `fp.m61.pow_mod.d64_k2` | 163.777 ms | 28.157 ms | 1.755 ms | 1.479 ms | 264.2 us | 93.307x | 1.706 ms | 1.587 ms | 633.1 us | 96.021x | ok `394fb320787700fe` |
| `fp.m61.is_irreducible.d24` | 520.920 ms | 95.956 ms | 7.224 ms | 6.401 ms | 995.2 us | 72.108x | 2.904 ms | 2.637 ms | 671.4 us | 179.355x | ok `9047cf1d3fc29abf` |
| `fp.m61.factor.d26` | 58.131 ms | 11.659 ms | 1.421 ms | 1.185 ms | 241.1 us | 40.911x | 1.447 ms | 1.273 ms | 352.8 us | 40.171x | ok `931b8d0b8b91b0b4` |
| `fp.m127.elem_muladd.n200000` | 29.344 ms | 6.165 ms | 29.311 ms | 25.162 ms | 3.530 ms | 1.001x | 30.507 ms | 26.129 ms | 11.961 ms | 0.962x | ok `16fab7b139aa9aab` |
| `fp.m127.elem_inv.n50000` | 34.996 ms | 7.539 ms | 37.353 ms | 29.968 ms | 11.020 ms | 0.937x | 37.782 ms | 31.586 ms | 8.811 ms | 0.926x | ok `c7d67d6018036bf1` |
| `fp.m127.elem_pow.n10000` | 26.622 ms | 5.126 ms | 30.903 ms | 24.298 ms | 14.961 ms | 0.861x | 28.794 ms | 24.823 ms | 7.309 ms | 0.925x | ok `308896e55b40c815` |
| `fp.m127.mul_school.d64` | 11.218 ms | 1.747 ms | 1.976 ms | 1.679 ms | 265.4 us | 5.676x | 1.953 ms | 1.706 ms | 287.3 us | 5.744x | ok `a9664408b0fb88e7` |
| `fp.m127.mul_school.d512` | 34.200 ms | 5.913 ms | 4.709 ms | 3.893 ms | 593.1 us | 7.263x | 4.602 ms | 4.037 ms | 875.8 us | 7.431x | ok `98eca357802203c0` |
| `fp.m127.mul_karatsuba.d512` | 16.073 ms | 1.946 ms | 2.082 ms | 1.555 ms | 454.6 us | 7.721x | 2.162 ms | 1.722 ms | 444.3 us | 7.433x | ok `f773c153fde34172` |
| `fp.m127.mul_karatsuba.d2048` | 154.834 ms | 22.410 ms | 13.711 ms | 11.000 ms | 4.234 ms | 11.293x | 13.795 ms | 11.056 ms | 3.940 ms | 11.224x | ok `695e0d0f13c14e87` |
| `fp.m127.mul_auto.d2048` | 152.004 ms | 27.912 ms | 13.694 ms | 10.980 ms | 3.560 ms | 11.100x | 13.943 ms | 11.277 ms | 4.286 ms | 10.902x | ok `1bc667b834b9e559` |
| `fp.m127.divmod.d2000_1000` | 141.670 ms | 41.837 ms | 18.702 ms | 15.066 ms | 5.922 ms | 7.575x | 17.670 ms | 15.415 ms | 3.688 ms | 8.017x | ok `a8babdcd119eb404` |
| `fp.m127.gcd.d400` | 25.413 ms | 5.861 ms | 18.157 ms | 14.849 ms | 5.694 ms | 1.400x | 17.575 ms | 15.097 ms | 5.505 ms | 1.446x | ok `1b15fac38c54accc` |
| `fp.m127.gcd_ext.d300` | 44.188 ms | 8.048 ms | 29.521 ms | 24.521 ms | 4.870 ms | 1.497x | 28.626 ms | 25.387 ms | 6.746 ms | 1.544x | ok `45a1b9001453c6ef` |
| `fp.m127.pow_mod.d64_k1` | 248.101 ms | 32.895 ms | 56.121 ms | 47.925 ms | 9.646 ms | 4.421x | 54.557 ms | 47.427 ms | 12.430 ms | 4.548x | ok `9e9e535b3c4e4e7d` |
| `fp.m127.is_irreducible.d16` | 362.519 ms | 58.601 ms | 122.106 ms | 100.122 ms | 35.919 ms | 2.969x | 49.226 ms | 41.368 ms | 8.547 ms | 7.364x | ok `2cfb69ee2746aa9d` |
| `fp.m127.factor.d19` | 93.483 ms | 23.189 ms | 48.057 ms | 38.104 ms | 14.602 ms | 1.945x | 44.604 ms | 37.807 ms | 10.614 ms | 2.096x | ok `7b0defab2417737c` |
| `fp.p256.elem_muladd.n200000` | 35.600 ms | 5.368 ms | 46.840 ms | 33.455 ms | 12.876 ms | 0.760x | 39.671 ms | 33.221 ms | 5.817 ms | 0.897x | ok `cc3cab838e9b19f1` |
| `fp.p256.elem_inv.n50000` | 72.783 ms | 10.893 ms | 79.333 ms | 71.322 ms | 22.491 ms | 0.917x | 78.426 ms | 66.334 ms | 11.163 ms | 0.928x | ok `89f920a8eba0c743` |
| `fp.p256.elem_pow.n2000` | 21.776 ms | 3.851 ms | 24.192 ms | 21.005 ms | 4.165 ms | 0.900x | 23.873 ms | 20.564 ms | 1.954 ms | 0.912x | ok `40472b256a3c877d` |
| `fp.p256.mul_school.d64` | 15.184 ms | 2.046 ms | 3.249 ms | 2.558 ms | 662.8 us | 4.674x | 3.150 ms | 2.594 ms | 637.5 us | 4.821x | ok `f5d3d87e3ed0977e` |
| `fp.p256.mul_school.d512` | 49.226 ms | 6.833 ms | 7.663 ms | 6.522 ms | 1.552 ms | 6.424x | 7.774 ms | 6.624 ms | 1.306 ms | 6.332x | ok `6fd19bd7b3ceab69` |
| `fp.p256.mul_karatsuba.d512` | 20.589 ms | 5.236 ms | 2.446 ms | 1.809 ms | 496.1 us | 8.418x | 2.508 ms | 1.794 ms | 484.2 us | 8.208x | ok `4b2f43bf0f878408` |
| `fp.p256.mul_karatsuba.d2048` | 200.443 ms | 32.049 ms | 18.175 ms | 13.259 ms | 5.726 ms | 11.029x | 17.429 ms | 14.333 ms | 3.519 ms | 11.500x | ok `b4b89a80c2a3aaa9` |
| `fp.p256.mul_fft.d4096` | 114.672 ms | 16.903 ms | 21.965 ms | 16.921 ms | 2.914 ms | 5.221x | 21.885 ms | 17.133 ms | 3.876 ms | 5.240x | ok `60ed510564513dfb` |
| `fp.p256.mul_fft.d8192` | 276.483 ms | 58.211 ms | 44.782 ms | 35.902 ms | 8.568 ms | 6.174x | 47.310 ms | 35.794 ms | 7.919 ms | 5.844x | ok `47812935df4de5fe` |
| `fp.p256.mul_auto.d2048` | 52.935 ms | 9.844 ms | 9.877 ms | 7.838 ms | 2.569 ms | 5.360x | 9.538 ms | 8.153 ms | 1.835 ms | 5.550x | ok `3d89402081500506` |
| `fp.p256.divmod.d1000_500` | 46.785 ms | 8.880 ms | 7.614 ms | 6.143 ms | 1.819 ms | 6.145x | 6.957 ms | 6.310 ms | 2.912 ms | 6.725x | ok `0bebce51d98eb3f4` |
| `fp.p256.gcd.d300` | 19.665 ms | 2.949 ms | 13.572 ms | 10.743 ms | 3.086 ms | 1.449x | 12.611 ms | 11.132 ms | 2.425 ms | 1.559x | ok `f69e4050cb1a5861` |
| `fp.p256.gcd_ext.d200` | 135.187 ms | 30.525 ms | 16.321 ms | 13.954 ms | 3.544 ms | 8.283x | 15.832 ms | 14.150 ms | 1.537 ms | 8.539x | ok `1a9df499c28b9f49` |
| `fp.p256.pow_mod.d64_k1` | 454.745 ms | 67.453 ms | 97.128 ms | 83.122 ms | 14.391 ms | 4.682x | 91.880 ms | 85.909 ms | 16.190 ms | 4.949x | ok `dbe2b595567fbb9c` |
| `fp.p256.is_irreducible.d12` | 384.215 ms | 57.374 ms | 147.863 ms | 128.439 ms | 38.262 ms | 2.598x | 46.576 ms | 41.665 ms | 7.858 ms | 8.249x | ok `7b839a40c98b0ef0` |
| `fp.p256.factor.d12` | 41.220 ms | 6.521 ms | 21.275 ms | 18.553 ms | 6.609 ms | 1.937x | 21.025 ms | 19.030 ms | 6.002 ms | 1.960x | ok `c33e3512ed388a7e` |
| `fq.q3e30.elem_muladd.n4000` | 224.674 ms | 53.834 ms | 6.924 ms | 5.927 ms | 925.1 us | 32.450x | 6.851 ms | 5.878 ms | 800.8 us | 32.792x | ok `16edfcc06a08028a` |
| `fq.q3e30.elem_inv.n100` | 381.140 ms | 48.138 ms | 11.741 ms | 10.584 ms | 2.249 ms | 32.462x | 848.8 us | 748.1 us | 105.7 us | 449.008x | ok `a8167c01d1cfc5dd` |
| `fq.q3e30.mul_school.d32` | 60.350 ms | 11.591 ms | 872.0 us | 744.8 us | 132.8 us | 69.209x | 872.4 us | 751.9 us | 180.7 us | 69.177x | ok `ed0ca8eb7d325de3` |
| `fq.q3e30.mul_karatsuba.d128` | 484.070 ms | 116.234 ms | 2.593 ms | 2.024 ms | 705.3 us | 186.705x | 2.671 ms | 2.198 ms | 414.7 us | 181.215x | ok `2464b901bd242c7d` |
| `fq.q3e30.divmod.d64_32` | 63.587 ms | 10.469 ms | 1.033 ms | 817.6 us | 246.8 us | 61.544x | 907.1 us | 796.8 us | 204.0 us | 70.099x | ok `a2b35f414a139741` |
| `fq.q3e30.factor.d5` | 482.621 ms | 101.797 ms | 13.029 ms | 11.237 ms | 2.814 ms | 37.042x | 7.939 ms | 7.361 ms | 1.801 ms | 60.790x | ok `417a5b753c7a2a07` |
| `fq.q2e8.elem_muladd.n50000` | 241.496 ms | 60.485 ms | 10.733 ms | 9.820 ms | 1.107 ms | 22.500x | 10.472 ms | 9.656 ms | 1.921 ms | 23.061x | ok `c0ffe16ee54f128c` |
| `fq.q2e8.elem_inv.n5000` | 273.032 ms | 56.228 ms | 14.234 ms | 12.345 ms | 2.102 ms | 19.182x | 2.691 ms | 2.469 ms | 367.2 us | 101.448x | ok `81d71b5687e3e735` |
| `fq.q2e8.mul_school.d128` | 80.773 ms | 21.615 ms | 719.1 us | 595.5 us | 174.4 us | 112.325x | 701.0 us | 626.3 us | 194.8 us | 115.225x | ok `c056d9d992271db9` |
| `fq.q2e8.mul_karatsuba.d512` | 446.636 ms | 115.381 ms | 2.567 ms | 2.063 ms | 630.1 us | 173.981x | 2.663 ms | 2.157 ms | 441.2 us | 167.744x | ok `de1089fdce26d02c` |
| `fq.q2e8.is_irreducible.d16` | 449.714 ms | 112.900 ms | 8.336 ms | 7.556 ms | 1.895 ms | 53.951x | 2.952 ms | 2.578 ms | 583.3 us | 152.368x | ok `9478d590a29a9cbf` |
| `fq.q2e8.factor.d28` | 108.890 ms | 28.117 ms | 3.547 ms | 3.003 ms | 931.1 us | 30.698x | 2.927 ms | 2.695 ms | 545.1 us | 37.207x | ok `d2fc95c10ececedd` |
| `fq.q7e4.elem_muladd.n50000` | 116.325 ms | 22.085 ms | 6.922 ms | 6.218 ms | 2.133 ms | 16.806x | 7.104 ms | 6.429 ms | 1.754 ms | 16.374x | ok `a6cf18a43d40c5c0` |
| `fq.q7e4.elem_inv.n5000` | 221.477 ms | 45.053 ms | 13.823 ms | 11.419 ms | 3.015 ms | 16.022x | 2.436 ms | 1.978 ms | 447.6 us | 90.915x | ok `8664ac6d09409f94` |
| `fq.q7e4.mul_school.d128` | 44.158 ms | 8.985 ms | 363.1 us | 278.9 us | 76.1 us | 121.615x | 339.3 us | 285.1 us | 67.3 us | 130.145x | ok `71d186560b1d8153` |
| `fq.q7e4.mul_karatsuba.d512` | 241.649 ms | 67.309 ms | 2.024 ms | 1.566 ms | 468.7 us | 119.421x | 2.072 ms | 1.649 ms | 372.3 us | 116.632x | ok `6a2e9b751ab48e99` |
| `fq.q7e4.gcd.d200` | 118.474 ms | 19.637 ms | 4.712 ms | 3.744 ms | 974.9 us | 25.146x | 4.027 ms | 3.530 ms | 895.8 us | 29.416x | ok `8cf19ac90e512feb` |
| `fq.q7e4.factor.d20` | 63.374 ms | 9.558 ms | 2.328 ms | 1.922 ms | 363.0 us | 27.220x | 1.662 ms | 1.477 ms | 168.8 us | 38.125x | ok `ce23318b1b5fd48e` |
| `fq.m61e2.elem_muladd.n100000` | 99.557 ms | 21.805 ms | 10.677 ms | 9.284 ms | 1.495 ms | 9.324x | 10.457 ms | 9.184 ms | 1.687 ms | 9.520x | ok `1814347c960518d8` |
| `fq.m61e2.elem_inv.n1000` | 165.917 ms | 34.966 ms | 19.018 ms | 15.170 ms | 3.350 ms | 8.724x | 454.9 us | 412.2 us | 64.9 us | 364.733x | ok `45326d8095032432` |
| `fq.m61e2.mul_school.d128` | 16.530 ms | 4.030 ms | 145.2 us | 120.3 us | 30.5 us | 113.885x | 142.7 us | 122.3 us | 28.0 us | 115.881x | ok `2cb38bd639886e3c` |
| `fq.m61e2.mul_karatsuba.d512` | 91.228 ms | 19.866 ms | 1.567 ms | 1.218 ms | 349.5 us | 58.205x | 1.636 ms | 1.253 ms | 272.2 us | 55.749x | ok `f138ad86c0a7a838` |
| `fq.m61e2.mul_fft.d1024` | 105.071 ms | 26.610 ms | 19.167 ms | 15.526 ms | 3.959 ms | 5.482x | 17.668 ms | 15.955 ms | 2.208 ms | 5.947x | ok `a45642d67f5a9912` |
| `fq.m61e2.mul_fft.d4096` | 581.071 ms | 49.551 ms | 53.501 ms | 43.996 ms | 13.356 ms | 10.861x | 51.003 ms | 45.524 ms | 8.188 ms | 11.393x | ok `51844e8f95a9c0ab` |
| `fq.m61e2.mul_auto.d2048` | 230.512 ms | 45.718 ms | 32.680 ms | 26.075 ms | 7.139 ms | 7.054x | 31.546 ms | 29.233 ms | 4.762 ms | 7.307x | ok `76e66ad22fb3fa2e` |
| `fq.m61e2.divmod.d500_250` | 62.296 ms | 18.318 ms | 562.7 us | 440.2 us | 118.4 us | 110.708x | 550.9 us | 492.0 us | 123.7 us | 113.090x | ok `53509a111940ee2a` |
| `fq.m61e2.gcd.d150` | 43.848 ms | 9.090 ms | 3.943 ms | 3.196 ms | 416.2 us | 11.120x | 1.514 ms | 1.402 ms | 138.2 us | 28.964x | ok `bd8c3c98cc56df84` |
| `fq.m61e2.factor.d7` | 177.981 ms | 37.062 ms | 23.252 ms | 19.834 ms | 4.742 ms | 7.655x | 3.347 ms | 3.159 ms | 974.3 us | 53.183x | ok `d61430d11ee187c3` |
| `fq.build.q3e30` | 25.288 ms | 4.417 ms | 581.6 us | 488.1 us | 130.0 us | 43.476x | 284.0 us | 257.4 us | 64.8 us | 89.042x | ok `4973a3e567e2fc74` |
| `fq.build.q7e20` | 21.297 ms | 4.539 ms | 463.4 us | 393.2 us | 106.2 us | 45.959x | 240.2 us | 218.8 us | 42.8 us | 88.664x | ok `9c377cd554aa108e` |
| `fast.ntt.dft.k10` | 16.755 ms | 5.564 ms | 1.498 ms | 1.288 ms | 394.6 us | 11.181x | 1.430 ms | 1.329 ms | 240.4 us | 11.720x | ok `b61f454e63ae6d8d` |
| `fast.ntt.dft.k11` | 17.934 ms | 4.170 ms | 1.503 ms | 1.292 ms | 523.3 us | 11.935x | 1.462 ms | 1.338 ms | 226.7 us | 12.263x | ok `311d4fed0f42034c` |
| `fast.ntt.dft.k12` | 19.016 ms | 2.439 ms | 1.531 ms | 1.302 ms | 215.0 us | 12.424x | 1.480 ms | 1.373 ms | 197.7 us | 12.850x | ok `70691e178cfc0c9d` |
| `fast.ntt.idft.k12` | 22.751 ms | 6.322 ms | 1.542 ms | 1.324 ms | 296.8 us | 14.756x | 1.505 ms | 1.393 ms | 394.0 us | 15.118x | ok `e0677b59c4e33615` |
| `fast.ntt.toeplitz_lower_inverse.n1024` | 31.057 ms | 28.428 ms | 713.7 us | 611.9 us | 121.2 us | 43.515x | 690.7 us | 636.5 us | 121.2 us | 44.964x | ok `7749eb8ff8cac3d0` |
| `fast.ntt.toeplitz_lower_inverse.n8192` | 342.472 ms | 62.472 ms | 7.576 ms | 6.436 ms | 1.351 ms | 45.207x | 7.714 ms | 6.869 ms | 1.096 ms | 44.397x | ok `0d31554b19d8b996` |
| `fast.ntt.toeplitz_lower_vec.n4096` | 42.040 ms | 8.226 ms | 988.8 us | 804.8 us | 179.6 us | 42.514x | 909.2 us | 824.8 us | 157.8 us | 46.239x | ok `86e36ad2a1eacb05` |
| `fast.ntt.toeplitz_divmod.d2000_1000` | 73.337 ms | 14.590 ms | 3.258 ms | 2.896 ms | 422.0 us | 22.508x | 3.197 ms | 3.080 ms | 175.9 us | 22.943x | ok `b138cdcceebcc32e` |
| `fast.p7.toeplitz_lower_inverse.n1024` | 87.326 ms | 17.663 ms | 1.160 ms | 967.8 us | 215.0 us | 75.275x | 1.104 ms | 1.046 ms | 179.3 us | 79.111x | ok `5df08c92a814155a` |
| `fast.p7.toeplitz_divmod.d2000_1000` | 95.295 ms | 20.591 ms | 1.906 ms | 1.688 ms | 211.7 us | 50.005x | 1.772 ms | 1.681 ms | 131.0 us | 53.781x | ok `433d13716b2f95e8` |
