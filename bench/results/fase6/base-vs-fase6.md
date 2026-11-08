# Comparativa de benchmarks alcp

Base: **f6-base**. Speedup = mediana(base) / mediana(variante); > 1 significa más rápido.
IQR = rango intercuartílico de las muestras. Las medias geométricas excluyen escenarios con error.

- **f6-base**: modo interleave, 2026-11-08T14:44:20+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3`  
  include `~\AppData\Local\Temp\claude\204af1c0d38157dc\2e18e96a-361a-47be-ad99-f0341eb0cc55\scratchpad\inc-2b36cde\include` @ ?, sha256 7bc05e5e10a5  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: no (batería)
- **f6-int**: modo interleave, 2026-11-08T14:44:41+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-81-g6221fa8, sha256 75031a6ad260  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: no (batería)
- **f6-desd3**: modo interleave, 2026-11-08T14:45:01+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -fno-ident -DV3`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-81-g6221fa8 (con cambios sin commit), sha256 01b563a80c52  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: no (batería)

Checksums: todos los escenarios coinciden con la base.

## Media geométrica del speedup

| ámbito | f6-int | f6-desd3 |
|---|---|---|
| **global** | 1.164x (n=128) | 1.169x (n=128) |
| grupo `fast` | 1.038x (n=10) | 1.066x (n=10) |
| grupo `fp.m127` | 1.046x (n=14) | 1.072x (n=14) |
| grupo `fp.m61` | 1.035x (n=14) | 1.014x (n=14) |
| grupo `fp.ntt` | 1.155x (n=16) | 1.141x (n=16) |
| grupo `fp.p101` | 1.483x (n=14) | 1.488x (n=14) |
| grupo `fp.p256` | 1.012x (n=16) | 1.001x (n=16) |
| grupo `fp.p7` | 1.521x (n=14) | 1.543x (n=14) |
| grupo `fq.build` | 1.322x (n=2) | 1.335x (n=2) |
| grupo `fq.m61e2` | 0.958x (n=10) | 1.004x (n=10) |
| grupo `fq.q2e8` | 1.224x (n=6) | 1.187x (n=6) |
| grupo `fq.q3e30` | 1.388x (n=6) | 1.399x (n=6) |
| grupo `fq.q7e4` | 1.140x (n=6) | 1.141x (n=6) |
| operación `build` | 1.322x (n=2) | 1.335x (n=2) |
| operación `dft` | 0.939x (n=3) | 0.999x (n=3) |
| operación `divmod` | 1.218x (n=8) | 1.220x (n=8) |
| operación `elem_inv` | 1.000x (n=10) | 0.980x (n=10) |
| operación `elem_muladd` | 1.046x (n=10) | 1.051x (n=10) |
| operación `elem_pow` | 0.996x (n=6) | 0.998x (n=6) |
| operación `factor` | 1.408x (n=10) | 1.410x (n=10) |
| operación `gcd` | 1.006x (n=8) | 1.012x (n=8) |
| operación `gcd_ext` | 1.009x (n=6) | 1.006x (n=6) |
| operación `idft` | 0.998x (n=1) | 1.029x (n=1) |
| operación `is_irreducible` | 1.215x (n=7) | 1.209x (n=7) |
| operación `mul_auto` | 1.033x (n=7) | 1.092x (n=7) |
| operación `mul_fft` | 1.007x (n=6) | 0.994x (n=6) |
| operación `mul_karatsuba` | 1.197x (n=16) | 1.212x (n=16) |
| operación `mul_school` | 1.527x (n=16) | 1.540x (n=16) |
| operación `pow_mod` | 1.335x (n=6) | 1.279x (n=6) |
| operación `toeplitz_divmod` | 1.008x (n=2) | 1.019x (n=2) |
| operación `toeplitz_lower_inverse` | 1.183x (n=3) | 1.201x (n=3) |
| operación `toeplitz_lower_vec` | 1.049x (n=1) | 1.025x (n=1) |

## Detalle por escenario

| escenario | base mediana | base IQR | f6-int mediana | f6-int mín | f6-int IQR | speedup | f6-desd3 mediana | f6-desd3 mín | f6-desd3 IQR | speedup | checksum |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `fp.p7.elem_muladd.n400000` | 4.175 ms | 415.1 us | 4.074 ms | 4.062 ms | 212.7 us | 1.025x | 4.113 ms | 4.062 ms | 174.9 us | 1.015x | ok `cf4797a827e7f9ec` |
| `fp.p7.elem_inv.n100000` | 817.2 us | 99.8 us | 815.5 us | 808.1 us | 38.8 us | 1.002x | 814.6 us | 808.1 us | 16.1 us | 1.003x | ok `7beb5a3238b919cb` |
| `fp.p7.elem_pow.n100000` | 2.252 ms | 48.2 us | 2.270 ms | 2.256 ms | 98.2 us | 0.992x | 2.276 ms | 2.257 ms | 125.3 us | 0.989x | ok `a7f8c429d0d81b73` |
| `fp.p7.mul_school.d64` | 171.1 us | 6.4 us | 55.3 us | 54.7 us | 0.3 us | 3.094x | 52.5 us | 51.9 us | 0.5 us | 3.259x | ok `559357b5bc3c57c2` |
| `fp.p7.mul_school.d512` | 407.1 us | 7.9 us | 127.7 us | 127.0 us | 4.9 us | 3.188x | 127.4 us | 127.0 us | 0.6 us | 3.195x | ok `5936b7480b151eb3` |
| `fp.p7.mul_karatsuba.d512` | 347.3 us | 93.8 us | 209.3 us | 186.6 us | 53.2 us | 1.660x | 210.2 us | 178.2 us | 115.2 us | 1.652x | ok `d72c59b658ac0379` |
| `fp.p7.mul_karatsuba.d4096` | 4.418 ms | 1.990 ms | 3.129 ms | 2.617 ms | 695.8 us | 1.412x | 2.672 ms | 2.239 ms | 423.1 us | 1.654x | ok `00f3dde93cbd25a4` |
| `fp.p7.mul_auto.d2048` | 178.2 us | 4.9 us | 122.0 us | 121.3 us | 2.0 us | 1.461x | 122.3 us | 121.5 us | 8.6 us | 1.457x | ok `2b256c747e0070c7` |
| `fp.p7.divmod.d2000_1000` | 353.8 us | 24.0 us | 214.2 us | 208.8 us | 25.8 us | 1.651x | 212.8 us | 207.8 us | 26.6 us | 1.663x | ok `958d383cda71b699` |
| `fp.p7.gcd.d500` | 626.7 us | 23.4 us | 625.1 us | 614.0 us | 13.0 us | 1.003x | 624.9 us | 612.6 us | 17.5 us | 1.003x | ok `cdf18334cf9b13d7` |
| `fp.p7.gcd_ext.d400` | 1.397 ms | 57.7 us | 1.387 ms | 1.362 ms | 59.3 us | 1.007x | 1.411 ms | 1.362 ms | 70.5 us | 0.990x | ok `984ac156886a3f82` |
| `fp.p7.pow_mod.d128_k48` | 4.175 ms | 598.5 us | 1.980 ms | 1.909 ms | 324.4 us | 2.108x | 1.932 ms | 1.901 ms | 109.5 us | 2.160x | ok `e3cedb8421cd962f` |
| `fp.p7.is_irreducible.d64` | 1.518 ms | 108.4 us | 949.2 us | 904.2 us | 156.7 us | 1.599x | 950.1 us | 903.3 us | 195.5 us | 1.598x | ok `9cb1bec87bc433d0` |
| `fp.p7.factor.d62` | 663.1 us | 93.8 us | 362.9 us | 338.5 us | 46.6 us | 1.827x | 362.7 us | 340.5 us | 41.9 us | 1.828x | ok `3425aecb393956ba` |
| `fp.p101.elem_muladd.n400000` | 4.076 ms | 98.9 us | 4.102 ms | 4.062 ms | 181.1 us | 0.994x | 4.116 ms | 4.062 ms | 497.9 us | 0.990x | ok `e51aa870ce6269eb` |
| `fp.p101.elem_inv.n100000` | 1.616 ms | 36.8 us | 1.613 ms | 1.603 ms | 30.2 us | 1.002x | 1.614 ms | 1.604 ms | 53.4 us | 1.001x | ok `69d137d626075364` |
| `fp.p101.elem_pow.n100000` | 4.919 ms | 644.3 us | 4.951 ms | 4.852 ms | 635.1 us | 0.994x | 4.985 ms | 4.769 ms | 302.5 us | 0.987x | ok `09fec44561c97afa` |
| `fp.p101.mul_school.d64` | 186.6 us | 6.3 us | 59.8 us | 59.3 us | 0.6 us | 3.120x | 56.9 us | 56.4 us | 0.6 us | 3.279x | ok `5a0fa852855ea8b0` |
| `fp.p101.mul_school.d512` | 481.4 us | 23.1 us | 149.4 us | 147.9 us | 3.2 us | 3.222x | 147.8 us | 147.4 us | 7.0 us | 3.257x | ok `8289633c899317a1` |
| `fp.p101.mul_karatsuba.d512` | 324.9 us | 19.4 us | 177.8 us | 168.2 us | 43.7 us | 1.827x | 175.3 us | 166.5 us | 31.6 us | 1.853x | ok `fea453dad54c0188` |
| `fp.p101.mul_karatsuba.d4096` | 3.940 ms | 669.6 us | 2.918 ms | 2.589 ms | 367.1 us | 1.351x | 2.966 ms | 2.488 ms | 294.0 us | 1.329x | ok `82f83f0b7fe7713b` |
| `fp.p101.mul_auto.d2048` | 240.2 us | 7.3 us | 180.4 us | 178.8 us | 9.2 us | 1.331x | 180.6 us | 178.7 us | 14.9 us | 1.330x | ok `3631fff5520b84d9` |
| `fp.p101.divmod.d2000_1000` | 457.6 us | 27.7 us | 306.2 us | 294.3 us | 33.7 us | 1.494x | 317.1 us | 295.7 us | 43.4 us | 1.443x | ok `ce32beed1dda91cd` |
| `fp.p101.gcd.d500` | 738.2 us | 39.6 us | 747.2 us | 726.9 us | 35.8 us | 0.988x | 751.7 us | 723.4 us | 82.8 us | 0.982x | ok `5cad0c7c5f168edb` |
| `fp.p101.gcd_ext.d400` | 1.601 ms | 102.5 us | 1.608 ms | 1.565 ms | 100.3 us | 0.995x | 1.613 ms | 1.560 ms | 137.1 us | 0.993x | ok `4fe0534b04a103a1` |
| `fp.p101.pow_mod.d128_k20` | 4.931 ms | 386.9 us | 2.668 ms | 2.471 ms | 1.224 ms | 1.848x | 2.592 ms | 2.468 ms | 257.0 us | 1.903x | ok `bc5334d07b025488` |
| `fp.p101.is_irreducible.d48` | 1.173 ms | 113.5 us | 698.0 us | 667.1 us | 196.9 us | 1.681x | 686.1 us | 666.0 us | 77.9 us | 1.710x | ok `5146dc1467eee05b` |
| `fp.p101.factor.d100` | 2.715 ms | 487.2 us | 1.629 ms | 1.523 ms | 190.6 us | 1.667x | 1.628 ms | 1.506 ms | 527.1 us | 1.668x | ok `48a7604aa4ae63ee` |
| `fp.ntt.elem_muladd.n400000` | 6.125 ms | 490.8 us | 6.613 ms | 5.935 ms | 1.460 ms | 0.926x | 6.079 ms | 5.930 ms | 703.6 us | 1.008x | ok `682b833fe4edd17b` |
| `fp.ntt.elem_inv.n100000` | 7.418 ms | 457.5 us | 7.451 ms | 7.111 ms | 1.149 ms | 0.996x | 7.693 ms | 7.200 ms | 466.7 us | 0.964x | ok `d275b741e6b3050e` |
| `fp.ntt.elem_pow.n50000` | 16.867 ms | 605.2 us | 16.563 ms | 16.364 ms | 694.3 us | 1.018x | 17.079 ms | 16.633 ms | 435.4 us | 0.988x | ok `c27b88ab2935c700` |
| `fp.ntt.mul_school.d64` | 192.3 us | 72.2 us | 129.2 us | 125.4 us | 16.9 us | 1.488x | 128.9 us | 125.4 us | 10.4 us | 1.492x | ok `67c7e1645a9cc6ae` |
| `fp.ntt.mul_school.d512` | 510.7 us | 42.5 us | 346.3 us | 343.7 us | 18.4 us | 1.475x | 347.3 us | 343.7 us | 19.7 us | 1.470x | ok `e296393a6a637046` |
| `fp.ntt.mul_karatsuba.d512` | 306.6 us | 36.8 us | 250.0 us | 235.8 us | 21.6 us | 1.226x | 240.4 us | 231.2 us | 17.5 us | 1.275x | ok `964e15063fb16df4` |
| `fp.ntt.mul_karatsuba.d4096` | 4.196 ms | 657.5 us | 3.741 ms | 2.898 ms | 1.035 ms | 1.121x | 3.251 ms | 2.752 ms | 863.5 us | 1.291x | ok `d656b678f0c2ca16` |
| `fp.ntt.mul_fft.d4096` | 1.777 ms | 155.5 us | 1.850 ms | 1.483 ms | 272.4 us | 0.961x | 1.744 ms | 1.508 ms | 158.8 us | 1.019x | ok `f4c81b85220fb683` |
| `fp.ntt.mul_fft.d32768` | 16.757 ms | 1.952 ms | 15.939 ms | 14.423 ms | 2.538 ms | 1.051x | 17.482 ms | 14.798 ms | 2.829 ms | 0.958x | ok `8228e76950faded3` |
| `fp.ntt.mul_auto.d2048` | 930.8 us | 103.2 us | 964.5 us | 825.3 us | 133.4 us | 0.965x | 921.0 us | 809.2 us | 112.1 us | 1.011x | ok `02e67c888996018d` |
| `fp.ntt.divmod.d2000_1000` | 1.196 ms | 90.7 us | 1.048 ms | 1.019 ms | 92.2 us | 1.142x | 1.053 ms | 1.020 ms | 123.6 us | 1.136x | ok `7adfbb2614325099` |
| `fp.ntt.gcd.d500` | 788.3 us | 44.5 us | 796.8 us | 767.1 us | 37.2 us | 0.989x | 804.4 us | 766.5 us | 46.7 us | 0.980x | ok `678962135e7a2447` |
| `fp.ntt.gcd_ext.d400` | 1.643 ms | 143.4 us | 1.631 ms | 1.587 ms | 61.8 us | 1.007x | 1.624 ms | 1.580 ms | 105.0 us | 1.012x | ok `f9dc4e400a7954ee` |
| `fp.ntt.pow_mod.d128_k5` | 9.788 ms | 2.837 ms | 7.170 ms | 6.554 ms | 797.9 us | 1.365x | 9.405 ms | 8.640 ms | 1.397 ms | 1.041x | ok `0eb2514aceb37261` |
| `fp.ntt.is_irreducible.d32` | 568.8 us | 26.7 us | 462.9 us | 450.1 us | 34.0 us | 1.229x | 501.6 us | 488.8 us | 27.0 us | 1.134x | ok `1f2827a0d8aeb65e` |
| `fp.ntt.factor.d40` | 3.061 ms | 399.1 us | 1.616 ms | 1.538 ms | 207.7 us | 1.894x | 1.684 ms | 1.566 ms | 389.1 us | 1.818x | ok `297c76c423fe55aa` |
| `fp.m61.elem_muladd.n400000` | 5.378 ms | 252.2 us | 5.307 ms | 5.227 ms | 154.4 us | 1.013x | 5.389 ms | 5.223 ms | 573.7 us | 0.998x | ok `b8e325f73ad90582` |
| `fp.m61.elem_inv.n100000` | 20.691 ms | 455.9 us | 20.559 ms | 20.150 ms | 265.8 us | 1.006x | 20.650 ms | 20.100 ms | 462.9 us | 1.002x | ok `41e7fe38de6fdcc1` |
| `fp.m61.elem_pow.n20000` | 15.005 ms | 640.6 us | 15.123 ms | 14.602 ms | 343.9 us | 0.992x | 15.096 ms | 14.573 ms | 690.2 us | 0.994x | ok `a3680ece7fc07932` |
| `fp.m61.mul_school.d64` | 199.4 us | 26.7 us | 187.4 us | 184.1 us | 95.0 us | 1.064x | 187.4 us | 185.6 us | 29.2 us | 1.064x | ok `bed6e9e86907de79` |
| `fp.m61.mul_school.d512` | 499.5 us | 21.5 us | 489.1 us | 474.5 us | 188.5 us | 1.021x | 489.3 us | 481.8 us | 32.9 us | 1.021x | ok `43c3884710f7d468` |
| `fp.m61.mul_karatsuba.d512` | 317.2 us | 27.0 us | 324.4 us | 298.2 us | 35.7 us | 0.978x | 332.8 us | 300.1 us | 61.5 us | 0.953x | ok `a2dcbeaa6574fcb8` |
| `fp.m61.mul_karatsuba.d4096` | 4.427 ms | 1.271 ms | 3.358 ms | 2.952 ms | 370.2 us | 1.319x | 4.529 ms | 3.058 ms | 1.233 ms | 0.978x | ok `37d5deaa1cfeb600` |
| `fp.m61.mul_auto.d2048` | 1.657 ms | 261.0 us | 1.677 ms | 1.532 ms | 174.6 us | 0.988x | 1.849 ms | 1.581 ms | 713.3 us | 0.896x | ok `73174ffa4e90d6fa` |
| `fp.m61.divmod.d2000_1000` | 1.941 ms | 118.9 us | 2.276 ms | 2.242 ms | 323.8 us | 0.853x | 2.045 ms | 1.992 ms | 434.0 us | 0.949x | ok `a8710fc5aaa24d09` |
| `fp.m61.gcd.d500` | 841.0 us | 57.1 us | 845.4 us | 812.4 us | 38.8 us | 0.995x | 857.6 us | 817.6 us | 70.8 us | 0.981x | ok `fda394dde025c1a2` |
| `fp.m61.gcd_ext.d400` | 1.678 ms | 167.5 us | 1.663 ms | 1.624 ms | 71.1 us | 1.009x | 1.691 ms | 1.634 ms | 263.9 us | 0.993x | ok `be3d81e8acc44b48` |
| `fp.m61.pow_mod.d64_k2` | 3.682 ms | 223.0 us | 3.820 ms | 3.659 ms | 432.4 us | 0.964x | 3.731 ms | 3.518 ms | 706.7 us | 0.987x | ok `394fb320787700fe` |
| `fp.m61.is_irreducible.d24` | 678.0 us | 62.8 us | 659.3 us | 643.2 us | 30.9 us | 1.028x | 652.9 us | 635.2 us | 36.1 us | 1.038x | ok `9047cf1d3fc29abf` |
| `fp.m61.factor.d26` | 2.995 ms | 486.5 us | 2.191 ms | 2.048 ms | 572.6 us | 1.367x | 2.115 ms | 2.056 ms | 193.7 us | 1.416x | ok `931b8d0b8b91b0b4` |
| `fp.m127.elem_muladd.n200000` | 57.666 ms | 2.470 ms | 59.785 ms | 54.902 ms | 4.208 ms | 0.965x | 56.221 ms | 53.865 ms | 5.932 ms | 1.026x | ok `16fab7b139aa9aab` |
| `fp.m127.elem_inv.n50000` | 83.183 ms | 9.606 ms | 80.727 ms | 74.604 ms | 7.999 ms | 1.030x | 82.783 ms | 73.733 ms | 14.105 ms | 1.005x | ok `c7d67d6018036bf1` |
| `fp.m127.elem_pow.n10000` | 57.858 ms | 4.533 ms | 60.527 ms | 54.637 ms | 7.422 ms | 0.956x | 57.835 ms | 55.303 ms | 5.387 ms | 1.000x | ok `308896e55b40c815` |
| `fp.m127.mul_school.d64` | 3.895 ms | 730.4 us | 3.961 ms | 3.784 ms | 708.2 us | 0.983x | 3.868 ms | 3.775 ms | 567.2 us | 1.007x | ok `a9664408b0fb88e7` |
| `fp.m127.mul_school.d512` | 11.269 ms | 2.628 ms | 9.397 ms | 8.936 ms | 998.7 us | 1.199x | 9.337 ms | 8.900 ms | 1.183 ms | 1.207x | ok `98eca357802203c0` |
| `fp.m127.mul_karatsuba.d512` | 5.157 ms | 1.207 ms | 4.767 ms | 3.633 ms | 787.1 us | 1.082x | 4.361 ms | 3.590 ms | 680.4 us | 1.182x | ok `f773c153fde34172` |
| `fp.m127.mul_karatsuba.d2048` | 22.170 ms | 3.994 ms | 23.104 ms | 20.530 ms | 4.612 ms | 0.960x | 21.024 ms | 20.203 ms | 1.373 ms | 1.054x | ok `695e0d0f13c14e87` |
| `fp.m127.mul_auto.d2048` | 4.124 ms | 303.2 us | 4.115 ms | 3.780 ms | 594.0 us | 1.002x | 4.199 ms | 3.889 ms | 485.6 us | 0.982x | ok `1bc667b834b9e559` |
| `fp.m127.divmod.d2000_1000` | 11.897 ms | 2.181 ms | 10.216 ms | 9.118 ms | 1.097 ms | 1.164x | 10.363 ms | 9.482 ms | 1.301 ms | 1.148x | ok `a8babdcd119eb404` |
| `fp.m127.gcd.d400` | 30.942 ms | 3.096 ms | 30.376 ms | 28.675 ms | 2.692 ms | 1.019x | 29.547 ms | 28.415 ms | 1.798 ms | 1.047x | ok `1b15fac38c54accc` |
| `fp.m127.gcd_ext.d300` | 30.919 ms | 1.682 ms | 32.850 ms | 30.975 ms | 4.300 ms | 0.941x | 31.333 ms | 29.807 ms | 5.008 ms | 0.987x | ok `45a1b9001453c6ef` |
| `fp.m127.pow_mod.d64_k1` | 76.939 ms | 14.561 ms | 68.196 ms | 64.518 ms | 3.182 ms | 1.128x | 70.053 ms | 65.862 ms | 13.587 ms | 1.098x | ok `9e9e535b3c4e4e7d` |
| `fp.m127.is_irreducible.d16` | 13.448 ms | 1.963 ms | 13.691 ms | 12.306 ms | 1.395 ms | 0.982x | 13.267 ms | 12.225 ms | 1.751 ms | 1.014x | ok `2cfb69ee2746aa9d` |
| `fp.m127.factor.d19` | 104.751 ms | 18.271 ms | 80.200 ms | 77.655 ms | 5.362 ms | 1.306x | 80.132 ms | 74.807 ms | 4.652 ms | 1.307x | ok `7b0defab2417737c` |
| `fp.p256.elem_muladd.n200000` | 80.432 ms | 4.528 ms | 81.064 ms | 75.246 ms | 8.352 ms | 0.992x | 80.513 ms | 77.456 ms | 3.340 ms | 0.999x | ok `cc3cab838e9b19f1` |
| `fp.p256.elem_inv.n50000` | 171.105 ms | 18.505 ms | 168.316 ms | 158.211 ms | 4.617 ms | 1.017x | 171.150 ms | 160.770 ms | 8.196 ms | 1.000x | ok `89f920a8eba0c743` |
| `fp.p256.elem_pow.n2000` | 49.667 ms | 7.635 ms | 48.553 ms | 45.792 ms | 5.059 ms | 1.023x | 48.203 ms | 45.408 ms | 3.039 ms | 1.030x | ok `40472b256a3c877d` |
| `fp.p256.mul_school.d64` | 5.963 ms | 726.1 us | 6.294 ms | 5.784 ms | 956.2 us | 0.947x | 6.015 ms | 5.734 ms | 480.7 us | 0.991x | ok `f5d3d87e3ed0977e` |
| `fp.p256.mul_school.d512` | 15.332 ms | 484.6 us | 16.112 ms | 14.845 ms | 1.236 ms | 0.952x | 16.141 ms | 14.459 ms | 1.534 ms | 0.950x | ok `6fd19bd7b3ceab69` |
| `fp.p256.mul_karatsuba.d512` | 5.677 ms | 1.020 ms | 5.623 ms | 4.746 ms | 1.024 ms | 1.010x | 5.812 ms | 4.965 ms | 774.5 us | 0.977x | ok `4b2f43bf0f878408` |
| `fp.p256.mul_karatsuba.d2048` | 24.511 ms | 1.093 ms | 25.997 ms | 23.737 ms | 4.510 ms | 0.943x | 24.668 ms | 23.329 ms | 3.923 ms | 0.994x | ok `b4b89a80c2a3aaa9` |
| `fp.p256.mul_fft.d4096` | 39.815 ms | 7.749 ms | 42.941 ms | 34.525 ms | 5.714 ms | 0.927x | 43.004 ms | 36.551 ms | 6.728 ms | 0.926x | ok `60ed510564513dfb` |
| `fp.p256.mul_fft.d8192` | 87.136 ms | 12.168 ms | 82.286 ms | 76.811 ms | 11.629 ms | 1.059x | 86.897 ms | 83.166 ms | 5.808 ms | 1.003x | ok `47812935df4de5fe` |
| `fp.p256.mul_auto.d2048` | 8.836 ms | 1.731 ms | 8.704 ms | 8.068 ms | 1.967 ms | 1.015x | 8.773 ms | 8.125 ms | 756.4 us | 1.007x | ok `3d89402081500506` |
| `fp.p256.divmod.d1000_500` | 8.946 ms | 1.193 ms | 8.479 ms | 7.649 ms | 980.9 us | 1.055x | 8.652 ms | 7.864 ms | 1.082 ms | 1.034x | ok `0bebce51d98eb3f4` |
| `fp.p256.gcd.d300` | 27.344 ms | 1.860 ms | 27.478 ms | 25.351 ms | 1.707 ms | 0.995x | 26.455 ms | 25.102 ms | 2.208 ms | 1.034x | ok `f69e4050cb1a5861` |
| `fp.p256.gcd_ext.d200` | 33.416 ms | 6.518 ms | 30.363 ms | 28.613 ms | 1.829 ms | 1.101x | 31.453 ms | 28.562 ms | 4.685 ms | 1.062x | ok `1a9df499c28b9f49` |
| `fp.p256.pow_mod.d64_k1` | 168.868 ms | 25.515 ms | 172.406 ms | 163.102 ms | 13.061 ms | 0.979x | 178.536 ms | 166.349 ms | 21.259 ms | 0.946x | ok `dbe2b595567fbb9c` |
| `fp.p256.is_irreducible.d12` | 17.638 ms | 718.4 us | 17.334 ms | 16.442 ms | 780.8 us | 1.018x | 17.736 ms | 16.623 ms | 7.312 ms | 0.994x | ok `7b839a40c98b0ef0` |
| `fp.p256.factor.d12` | 44.775 ms | 2.391 ms | 37.454 ms | 34.296 ms | 1.937 ms | 1.195x | 41.192 ms | 35.113 ms | 6.238 ms | 1.087x | ok `c33e3512ed388a7e` |
| `fq.q3e30.elem_muladd.n4000` | 14.001 ms | 428.6 us | 10.486 ms | 9.223 ms | 722.6 us | 1.335x | 10.920 ms | 9.308 ms | 2.936 ms | 1.282x | ok `16edfcc06a08028a` |
| `fq.q3e30.elem_inv.n100` | 1.822 ms | 190.9 us | 1.840 ms | 1.729 ms | 278.3 us | 0.990x | 2.010 ms | 1.728 ms | 495.2 us | 0.907x | ok `a8167c01d1cfc5dd` |
| `fq.q3e30.mul_school.d32` | 1.756 ms | 187.9 us | 1.023 ms | 1.003 ms | 93.6 us | 1.717x | 1.031 ms | 1.006 ms | 213.1 us | 1.704x | ok `ed0ca8eb7d325de3` |
| `fq.q3e30.mul_karatsuba.d128` | 5.936 ms | 1.333 ms | 4.430 ms | 3.672 ms | 709.5 us | 1.340x | 3.777 ms | 2.988 ms | 1.187 ms | 1.572x | ok `2464b901bd242c7d` |
| `fq.q3e30.divmod.d64_32` | 1.779 ms | 161.0 us | 1.065 ms | 1.042 ms | 48.1 us | 1.671x | 1.075 ms | 1.040 ms | 509.3 us | 1.655x | ok `a2b35f414a139741` |
| `fq.q3e30.factor.d5` | 8.932 ms | 451.3 us | 6.341 ms | 5.891 ms | 541.8 us | 1.409x | 6.148 ms | 5.963 ms | 676.3 us | 1.453x | ok `417a5b753c7a2a07` |
| `fq.q2e8.elem_muladd.n50000` | 22.198 ms | 1.418 ms | 19.238 ms | 18.226 ms | 2.466 ms | 1.154x | 19.888 ms | 18.047 ms | 2.971 ms | 1.116x | ok `c0ffe16ee54f128c` |
| `fq.q2e8.elem_inv.n5000` | 6.102 ms | 859.5 us | 6.132 ms | 5.411 ms | 972.6 us | 0.995x | 6.330 ms | 5.403 ms | 1.034 ms | 0.964x | ok `81d71b5687e3e735` |
| `fq.q2e8.mul_school.d128` | 1.453 ms | 153.6 us | 1.005 ms | 983.0 us | 29.9 us | 1.446x | 1.021 ms | 979.5 us | 134.9 us | 1.423x | ok `c056d9d992271db9` |
| `fq.q2e8.mul_karatsuba.d512` | 6.036 ms | 1.777 ms | 4.468 ms | 3.539 ms | 1.234 ms | 1.351x | 5.038 ms | 3.784 ms | 1.011 ms | 1.198x | ok `de1089fdce26d02c` |
| `fq.q2e8.is_irreducible.d16` | 3.409 ms | 364.5 us | 2.954 ms | 2.821 ms | 267.0 us | 1.154x | 2.926 ms | 2.810 ms | 629.7 us | 1.165x | ok `9478d590a29a9cbf` |
| `fq.q2e8.factor.d28` | 6.779 ms | 863.5 us | 5.224 ms | 5.011 ms | 405.9 us | 1.298x | 5.172 ms | 4.943 ms | 646.1 us | 1.311x | ok `d2fc95c10ececedd` |
| `fq.q7e4.elem_muladd.n50000` | 14.852 ms | 625.4 us | 13.467 ms | 13.055 ms | 550.5 us | 1.103x | 13.235 ms | 12.673 ms | 743.4 us | 1.122x | ok `a6cf18a43d40c5c0` |
| `fq.q7e4.elem_inv.n5000` | 4.433 ms | 280.7 us | 4.545 ms | 4.244 ms | 638.9 us | 0.975x | 4.595 ms | 4.266 ms | 663.3 us | 0.965x | ok `8664ac6d09409f94` |
| `fq.q7e4.mul_school.d128` | 662.4 us | 51.5 us | 468.9 us | 451.3 us | 104.6 us | 1.413x | 466.3 us | 451.6 us | 39.4 us | 1.421x | ok `71d186560b1d8153` |
| `fq.q7e4.mul_karatsuba.d512` | 4.166 ms | 847.7 us | 3.867 ms | 3.284 ms | 715.1 us | 1.077x | 3.915 ms | 3.122 ms | 492.0 us | 1.064x | ok `6a2e9b751ab48e99` |
| `fq.q7e4.gcd.d200` | 8.215 ms | 619.1 us | 7.469 ms | 7.126 ms | 482.7 us | 1.100x | 7.561 ms | 7.089 ms | 525.7 us | 1.086x | ok `8cf19ac90e512feb` |
| `fq.q7e4.factor.d20` | 2.986 ms | 817.4 us | 2.452 ms | 2.335 ms | 313.4 us | 1.218x | 2.407 ms | 2.344 ms | 195.6 us | 1.241x | ok `ce23318b1b5fd48e` |
| `fq.m61e2.elem_muladd.n100000` | 22.515 ms | 4.703 ms | 22.265 ms | 21.365 ms | 534.2 us | 1.011x | 22.766 ms | 21.615 ms | 6.890 ms | 0.989x | ok `1814347c960518d8` |
| `fq.m61e2.elem_inv.n1000` | 946.3 us | 42.3 us | 953.5 us | 928.4 us | 98.6 us | 0.992x | 953.6 us | 934.9 us | 170.8 us | 0.992x | ok `45326d8095032432` |
| `fq.m61e2.mul_school.d128` | 282.2 us | 12.0 us | 286.0 us | 271.1 us | 21.0 us | 0.987x | 293.7 us | 278.1 us | 22.0 us | 0.961x | ok `2cb38bd639886e3c` |
| `fq.m61e2.mul_karatsuba.d512` | 3.079 ms | 500.1 us | 3.444 ms | 2.827 ms | 337.1 us | 0.894x | 3.539 ms | 3.025 ms | 732.8 us | 0.870x | ok `f138ad86c0a7a838` |
| `fq.m61e2.mul_fft.d1024` | 36.804 ms | 1.319 ms | 35.902 ms | 33.142 ms | 1.322 ms | 1.025x | 35.809 ms | 33.287 ms | 1.991 ms | 1.028x | ok `a45642d67f5a9912` |
| `fq.m61e2.mul_fft.d4096` | 92.539 ms | 4.720 ms | 90.452 ms | 84.523 ms | 3.953 ms | 1.023x | 89.358 ms | 85.628 ms | 4.121 ms | 1.036x | ok `51844e8f95a9c0ab` |
| `fq.m61e2.mul_auto.d2048` | 6.398 ms | 836.9 us | 9.647 ms | 5.376 ms | 4.022 ms | 0.663x | 6.005 ms | 5.402 ms | 437.0 us | 1.065x | ok `76e66ad22fb3fa2e` |
| `fq.m61e2.divmod.d500_250` | 1.071 ms | 60.8 us | 1.093 ms | 1.037 ms | 74.6 us | 0.980x | 1.107 ms | 1.061 ms | 296.6 us | 0.968x | ok `53509a111940ee2a` |
| `fq.m61e2.gcd.d150` | 3.179 ms | 567.9 us | 3.290 ms | 3.174 ms | 391.2 us | 0.966x | 3.224 ms | 3.165 ms | 307.0 us | 0.986x | ok `bd8c3c98cc56df84` |
| `fq.m61e2.factor.d7` | 7.288 ms | 763.7 us | 6.546 ms | 5.828 ms | 968.0 us | 1.113x | 6.193 ms | 5.836 ms | 778.4 us | 1.177x | ok `d61430d11ee187c3` |
| `fq.build.q3e30` | 381.3 us | 38.8 us | 278.4 us | 262.6 us | 67.9 us | 1.370x | 277.2 us | 261.7 us | 58.3 us | 1.376x | ok `4973a3e567e2fc74` |
| `fq.build.q7e20` | 343.4 us | 22.1 us | 269.1 us | 255.4 us | 105.3 us | 1.276x | 265.1 us | 253.4 us | 54.2 us | 1.295x | ok `9c377cd554aa108e` |
| `fast.ntt.dft.k10` | 2.962 ms | 561.3 us | 3.166 ms | 2.924 ms | 983.4 us | 0.936x | 2.942 ms | 2.840 ms | 355.1 us | 1.007x | ok `b61f454e63ae6d8d` |
| `fast.ntt.dft.k11` | 2.942 ms | 351.1 us | 3.151 ms | 2.963 ms | 1.624 ms | 0.934x | 2.936 ms | 2.882 ms | 505.2 us | 1.002x | ok `311d4fed0f42034c` |
| `fast.ntt.dft.k12` | 2.977 ms | 355.1 us | 3.141 ms | 2.975 ms | 536.3 us | 0.948x | 3.009 ms | 2.920 ms | 569.6 us | 0.989x | ok `70691e178cfc0c9d` |
| `fast.ntt.idft.k12` | 3.107 ms | 494.6 us | 3.112 ms | 2.996 ms | 811.8 us | 0.998x | 3.018 ms | 2.956 ms | 992.0 us | 1.029x | ok `e0677b59c4e33615` |
| `fast.ntt.toeplitz_lower_inverse.n1024` | 569.8 us | 139.5 us | 571.8 us | 543.0 us | 150.4 us | 0.997x | 564.2 us | 541.4 us | 88.3 us | 1.010x | ok `7749eb8ff8cac3d0` |
| `fast.ntt.toeplitz_lower_inverse.n8192` | 7.058 ms | 729.3 us | 7.141 ms | 6.294 ms | 525.7 us | 0.988x | 6.981 ms | 6.436 ms | 756.3 us | 1.011x | ok `0d31554b19d8b996` |
| `fast.ntt.toeplitz_lower_vec.n4096` | 1.016 ms | 157.7 us | 968.2 us | 864.7 us | 89.4 us | 1.049x | 990.9 us | 876.2 us | 104.5 us | 1.025x | ok `86e36ad2a1eacb05` |
| `fast.ntt.toeplitz_divmod.d2000_1000` | 6.452 ms | 206.2 us | 6.644 ms | 6.416 ms | 569.1 us | 0.971x | 6.654 ms | 6.396 ms | 373.1 us | 0.970x | ok `b138cdcceebcc32e` |
| `fast.p7.toeplitz_lower_inverse.n1024` | 239.1 us | 27.2 us | 142.3 us | 135.3 us | 25.3 us | 1.680x | 141.0 us | 133.6 us | 20.0 us | 1.696x | ok `5df08c92a814155a` |
| `fast.p7.toeplitz_divmod.d2000_1000` | 3.036 ms | 589.4 us | 2.901 ms | 2.795 ms | 464.5 us | 1.046x | 2.836 ms | 2.795 ms | 198.1 us | 1.070x | ok `433d13716b2f95e8` |
