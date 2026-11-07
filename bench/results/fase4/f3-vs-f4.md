# Comparativa de benchmarks alcp

Base: **f4-f3**. Speedup = mediana(base) / mediana(variante); > 1 significa más rápido.
IQR = rango intercuartílico de las muestras. Las medias geométricas excluyen escenarios con error.

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

| ámbito | f4-f4 |
|---|---|
| **global** | 1.192x (n=128) |
| grupo `fast` | 1.038x (n=10) |
| grupo `fp.m127` | 1.084x (n=14) |
| grupo `fp.m61` | 1.038x (n=14) |
| grupo `fp.ntt` | 1.055x (n=16) |
| grupo `fp.p101` | 1.038x (n=14) |
| grupo `fp.p256` | 1.107x (n=16) |
| grupo `fp.p7` | 1.066x (n=14) |
| grupo `fq.build` | 1.988x (n=2) |
| grupo `fq.m61e2` | 1.976x (n=10) |
| grupo `fq.q2e8` | 1.624x (n=6) |
| grupo `fq.q3e30` | 1.714x (n=6) |
| grupo `fq.q7e4` | 1.454x (n=6) |
| operación `build` | 1.988x (n=2) |
| operación `dft` | 1.037x (n=3) |
| operación `divmod` | 1.040x (n=8) |
| operación `elem_inv` | 2.647x (n=10) |
| operación `elem_muladd` | 1.017x (n=10) |
| operación `elem_pow` | 0.993x (n=6) |
| operación `factor` | 1.327x (n=10) |
| operación `gcd` | 1.172x (n=8) |
| operación `gcd_ext` | 1.003x (n=6) |
| operación `idft` | 1.025x (n=1) |
| operación `is_irreducible` | 2.689x (n=7) |
| operación `mul_auto` | 1.003x (n=7) |
| operación `mul_fft` | 0.996x (n=6) |
| operación `mul_karatsuba` | 0.972x (n=16) |
| operación `mul_school` | 0.996x (n=16) |
| operación `pow_mod` | 1.009x (n=6) |
| operación `toeplitz_divmod` | 1.047x (n=2) |
| operación `toeplitz_lower_inverse` | 1.022x (n=3) |
| operación `toeplitz_lower_vec` | 1.088x (n=1) |

## Detalle por escenario

| escenario | base mediana | base IQR | f4-f4 mediana | f4-f4 mín | f4-f4 IQR | speedup | checksum |
|---|---:|---:|---:|---:|---:|---:|---|
| `fp.p7.elem_muladd.n400000` | 2.168 ms | 638.0 us | 2.026 ms | 1.706 ms | 406.7 us | 1.070x | ok `cf4797a827e7f9ec` |
| `fp.p7.elem_inv.n100000` | 392.5 us | 102.6 us | 394.8 us | 334.1 us | 59.8 us | 0.994x | ok `7beb5a3238b919cb` |
| `fp.p7.elem_pow.n100000` | 1.086 ms | 266.8 us | 1.143 ms | 975.9 us | 236.2 us | 0.950x | ok `a7f8c429d0d81b73` |
| `fp.p7.mul_school.d64` | 81.3 us | 19.9 us | 87.1 us | 69.8 us | 18.6 us | 0.934x | ok `559357b5bc3c57c2` |
| `fp.p7.mul_school.d512` | 204.5 us | 35.3 us | 206.0 us | 165.5 us | 36.2 us | 0.993x | ok `5936b7480b151eb3` |
| `fp.p7.mul_karatsuba.d512` | 158.9 us | 38.2 us | 160.2 us | 124.3 us | 35.8 us | 0.992x | ok `d72c59b658ac0379` |
| `fp.p7.mul_karatsuba.d4096` | 1.668 ms | 393.6 us | 1.730 ms | 1.291 ms | 385.3 us | 0.964x | ok `00f3dde93cbd25a4` |
| `fp.p7.mul_auto.d2048` | 979.2 us | 285.7 us | 918.0 us | 679.6 us | 328.9 us | 1.067x | ok `2b256c747e0070c7` |
| `fp.p7.divmod.d2000_1000` | 874.5 us | 167.3 us | 816.4 us | 676.7 us | 215.3 us | 1.071x | ok `958d383cda71b699` |
| `fp.p7.gcd.d500` | 1.070 ms | 209.9 us | 998.5 us | 856.9 us | 207.1 us | 1.072x | ok `cdf18334cf9b13d7` |
| `fp.p7.gcd_ext.d400` | 2.349 ms | 482.3 us | 2.336 ms | 1.902 ms | 469.6 us | 1.006x | ok `984ac156886a3f82` |
| `fp.p7.pow_mod.d128_k48` | 5.649 ms | 353.7 us | 5.614 ms | 4.795 ms | 1.184 ms | 1.006x | ok `e3cedb8421cd962f` |
| `fp.p7.is_irreducible.d64` | 4.322 ms | 1.201 ms | 1.936 ms | 1.659 ms | 480.5 us | 2.232x | ok `9cb1bec87bc433d0` |
| `fp.p7.factor.d62` | 338.2 us | 93.9 us | 342.2 us | 279.3 us | 99.3 us | 0.988x | ok `3425aecb393956ba` |
| `fp.p101.elem_muladd.n400000` | 1.911 ms | 415.0 us | 1.979 ms | 1.686 ms | 484.8 us | 0.966x | ok `e51aa870ce6269eb` |
| `fp.p101.elem_inv.n100000` | 799.9 us | 167.8 us | 762.0 us | 667.3 us | 188.1 us | 1.050x | ok `69d137d626075364` |
| `fp.p101.elem_pow.n100000` | 2.477 ms | 411.7 us | 2.392 ms | 2.027 ms | 642.7 us | 1.036x | ok `09fec44561c97afa` |
| `fp.p101.mul_school.d64` | 87.9 us | 14.3 us | 87.4 us | 74.5 us | 25.7 us | 1.006x | ok `5a0fa852855ea8b0` |
| `fp.p101.mul_school.d512` | 227.3 us | 60.9 us | 225.3 us | 195.1 us | 61.6 us | 1.009x | ok `8289633c899317a1` |
| `fp.p101.mul_karatsuba.d512` | 143.6 us | 44.0 us | 158.5 us | 121.8 us | 35.1 us | 0.906x | ok `fea453dad54c0188` |
| `fp.p101.mul_karatsuba.d4096` | 1.735 ms | 507.2 us | 1.755 ms | 1.393 ms | 396.4 us | 0.989x | ok `82f83f0b7fe7713b` |
| `fp.p101.mul_auto.d2048` | 1.007 ms | 290.6 us | 1.072 ms | 774.5 us | 289.2 us | 0.939x | ok `3631fff5520b84d9` |
| `fp.p101.divmod.d2000_1000` | 1.034 ms | 210.1 us | 1.105 ms | 910.1 us | 219.2 us | 0.936x | ok `ce32beed1dda91cd` |
| `fp.p101.gcd.d500` | 1.357 ms | 332.9 us | 1.406 ms | 1.112 ms | 238.5 us | 0.965x | ok `5cad0c7c5f168edb` |
| `fp.p101.gcd_ext.d400` | 2.779 ms | 501.0 us | 2.885 ms | 2.326 ms | 664.3 us | 0.963x | ok `4fe0534b04a103a1` |
| `fp.p101.pow_mod.d128_k20` | 6.123 ms | 1.190 ms | 6.562 ms | 5.443 ms | 1.603 ms | 0.933x | ok `bc5334d07b025488` |
| `fp.p101.is_irreducible.d48` | 4.814 ms | 1.005 ms | 1.965 ms | 1.518 ms | 492.4 us | 2.449x | ok `5146dc1467eee05b` |
| `fp.p101.factor.d100` | 1.795 ms | 280.5 us | 1.896 ms | 1.512 ms | 447.6 us | 0.947x | ok `48a7604aa4ae63ee` |
| `fp.ntt.elem_muladd.n400000` | 2.787 ms | 388.1 us | 2.864 ms | 2.531 ms | 712.8 us | 0.973x | ok `682b833fe4edd17b` |
| `fp.ntt.elem_inv.n100000` | 3.336 ms | 442.9 us | 3.538 ms | 2.970 ms | 681.4 us | 0.943x | ok `d275b741e6b3050e` |
| `fp.ntt.elem_pow.n50000` | 7.266 ms | 1.157 ms | 7.541 ms | 6.797 ms | 1.619 ms | 0.964x | ok `c27b88ab2935c700` |
| `fp.ntt.mul_school.d64` | 93.5 us | 16.9 us | 91.6 us | 76.9 us | 25.7 us | 1.021x | ok `67c7e1645a9cc6ae` |
| `fp.ntt.mul_school.d512` | 235.7 us | 32.0 us | 232.3 us | 195.7 us | 50.6 us | 1.014x | ok `e296393a6a637046` |
| `fp.ntt.mul_karatsuba.d512` | 151.9 us | 43.0 us | 156.5 us | 124.9 us | 35.2 us | 0.971x | ok `964e15063fb16df4` |
| `fp.ntt.mul_karatsuba.d4096` | 1.946 ms | 382.6 us | 2.014 ms | 1.659 ms | 354.9 us | 0.966x | ok `d656b678f0c2ca16` |
| `fp.ntt.mul_fft.d4096` | 2.040 ms | 487.1 us | 2.103 ms | 1.579 ms | 450.9 us | 0.970x | ok `f4c81b85220fb683` |
| `fp.ntt.mul_fft.d32768` | 17.586 ms | 3.424 ms | 18.911 ms | 15.772 ms | 5.257 ms | 0.930x | ok `8228e76950faded3` |
| `fp.ntt.mul_auto.d2048` | 931.1 us | 269.2 us | 925.7 us | 727.5 us | 260.0 us | 1.006x | ok `02e67c888996018d` |
| `fp.ntt.divmod.d2000_1000` | 968.5 us | 176.9 us | 919.2 us | 807.6 us | 91.5 us | 1.054x | ok `7adfbb2614325099` |
| `fp.ntt.gcd.d500` | 1.402 ms | 230.2 us | 1.398 ms | 1.160 ms | 140.9 us | 1.003x | ok `678962135e7a2447` |
| `fp.ntt.gcd_ext.d400` | 2.874 ms | 516.3 us | 2.987 ms | 2.426 ms | 475.9 us | 0.962x | ok `f9dc4e400a7954ee` |
| `fp.ntt.pow_mod.d128_k5` | 6.100 ms | 688.2 us | 6.069 ms | 5.652 ms | 455.0 us | 1.005x | ok `0eb2514aceb37261` |
| `fp.ntt.is_irreducible.d32` | 6.066 ms | 1.424 ms | 1.805 ms | 1.634 ms | 192.0 us | 3.360x | ok `1f2827a0d8aeb65e` |
| `fp.ntt.factor.d40` | 1.246 ms | 388.8 us | 1.425 ms | 1.160 ms | 454.3 us | 0.874x | ok `297c76c423fe55aa` |
| `fp.m61.elem_muladd.n400000` | 2.628 ms | 593.1 us | 2.614 ms | 2.403 ms | 744.4 us | 1.006x | ok `b8e325f73ad90582` |
| `fp.m61.elem_inv.n100000` | 9.324 ms | 1.545 ms | 9.446 ms | 8.302 ms | 2.185 ms | 0.987x | ok `41e7fe38de6fdcc1` |
| `fp.m61.elem_pow.n20000` | 6.637 ms | 1.038 ms | 7.147 ms | 6.431 ms | 1.673 ms | 0.929x | ok `a3680ece7fc07932` |
| `fp.m61.mul_school.d64` | 98.4 us | 23.8 us | 110.3 us | 85.9 us | 18.5 us | 0.892x | ok `bed6e9e86907de79` |
| `fp.m61.mul_school.d512` | 239.1 us | 68.7 us | 260.7 us | 216.7 us | 50.7 us | 0.917x | ok `43c3884710f7d468` |
| `fp.m61.mul_karatsuba.d512` | 155.0 us | 50.5 us | 161.5 us | 136.9 us | 33.6 us | 0.960x | ok `a2dcbeaa6574fcb8` |
| `fp.m61.mul_karatsuba.d4096` | 1.925 ms | 413.8 us | 1.997 ms | 1.582 ms | 381.0 us | 0.964x | ok `37d5deaa1cfeb600` |
| `fp.m61.mul_auto.d2048` | 1.073 ms | 350.5 us | 1.114 ms | 823.5 us | 278.9 us | 0.964x | ok `73174ffa4e90d6fa` |
| `fp.m61.divmod.d2000_1000` | 987.9 us | 249.1 us | 1.025 ms | 911.2 us | 201.6 us | 0.963x | ok `a8710fc5aaa24d09` |
| `fp.m61.gcd.d500` | 1.206 ms | 261.4 us | 1.191 ms | 1.082 ms | 111.6 us | 1.012x | ok `fda394dde025c1a2` |
| `fp.m61.gcd_ext.d400` | 2.529 ms | 475.2 us | 2.468 ms | 2.232 ms | 343.5 us | 1.025x | ok `be3d81e8acc44b48` |
| `fp.m61.pow_mod.d64_k2` | 1.755 ms | 264.2 us | 1.706 ms | 1.587 ms | 633.1 us | 1.029x | ok `394fb320787700fe` |
| `fp.m61.is_irreducible.d24` | 7.224 ms | 995.2 us | 2.904 ms | 2.637 ms | 671.4 us | 2.487x | ok `9047cf1d3fc29abf` |
| `fp.m61.factor.d26` | 1.421 ms | 241.1 us | 1.447 ms | 1.273 ms | 352.8 us | 0.982x | ok `931b8d0b8b91b0b4` |
| `fp.m127.elem_muladd.n200000` | 29.311 ms | 3.530 ms | 30.507 ms | 26.129 ms | 11.961 ms | 0.961x | ok `16fab7b139aa9aab` |
| `fp.m127.elem_inv.n50000` | 37.353 ms | 11.020 ms | 37.782 ms | 31.586 ms | 8.811 ms | 0.989x | ok `c7d67d6018036bf1` |
| `fp.m127.elem_pow.n10000` | 30.903 ms | 14.961 ms | 28.794 ms | 24.823 ms | 7.309 ms | 1.073x | ok `308896e55b40c815` |
| `fp.m127.mul_school.d64` | 1.976 ms | 265.4 us | 1.953 ms | 1.706 ms | 287.3 us | 1.012x | ok `a9664408b0fb88e7` |
| `fp.m127.mul_school.d512` | 4.709 ms | 593.1 us | 4.602 ms | 4.037 ms | 875.8 us | 1.023x | ok `98eca357802203c0` |
| `fp.m127.mul_karatsuba.d512` | 2.082 ms | 454.6 us | 2.162 ms | 1.722 ms | 444.3 us | 0.963x | ok `f773c153fde34172` |
| `fp.m127.mul_karatsuba.d2048` | 13.711 ms | 4.234 ms | 13.795 ms | 11.056 ms | 3.940 ms | 0.994x | ok `695e0d0f13c14e87` |
| `fp.m127.mul_auto.d2048` | 13.694 ms | 3.560 ms | 13.943 ms | 11.277 ms | 4.286 ms | 0.982x | ok `1bc667b834b9e559` |
| `fp.m127.divmod.d2000_1000` | 18.702 ms | 5.922 ms | 17.670 ms | 15.415 ms | 3.688 ms | 1.058x | ok `a8babdcd119eb404` |
| `fp.m127.gcd.d400` | 18.157 ms | 5.694 ms | 17.575 ms | 15.097 ms | 5.505 ms | 1.033x | ok `1b15fac38c54accc` |
| `fp.m127.gcd_ext.d300` | 29.521 ms | 4.870 ms | 28.626 ms | 25.387 ms | 6.746 ms | 1.031x | ok `45a1b9001453c6ef` |
| `fp.m127.pow_mod.d64_k1` | 56.121 ms | 9.646 ms | 54.557 ms | 47.427 ms | 12.430 ms | 1.029x | ok `9e9e535b3c4e4e7d` |
| `fp.m127.is_irreducible.d16` | 122.106 ms | 35.919 ms | 49.226 ms | 41.368 ms | 8.547 ms | 2.481x | ok `2cfb69ee2746aa9d` |
| `fp.m127.factor.d19` | 48.057 ms | 14.602 ms | 44.604 ms | 37.807 ms | 10.614 ms | 1.077x | ok `7b0defab2417737c` |
| `fp.p256.elem_muladd.n200000` | 46.840 ms | 12.876 ms | 39.671 ms | 33.221 ms | 5.817 ms | 1.181x | ok `cc3cab838e9b19f1` |
| `fp.p256.elem_inv.n50000` | 79.333 ms | 22.491 ms | 78.426 ms | 66.334 ms | 11.163 ms | 1.012x | ok `89f920a8eba0c743` |
| `fp.p256.elem_pow.n2000` | 24.192 ms | 4.165 ms | 23.873 ms | 20.564 ms | 1.954 ms | 1.013x | ok `40472b256a3c877d` |
| `fp.p256.mul_school.d64` | 3.249 ms | 662.8 us | 3.150 ms | 2.594 ms | 637.5 us | 1.031x | ok `f5d3d87e3ed0977e` |
| `fp.p256.mul_school.d512` | 7.663 ms | 1.552 ms | 7.774 ms | 6.624 ms | 1.306 ms | 0.986x | ok `6fd19bd7b3ceab69` |
| `fp.p256.mul_karatsuba.d512` | 2.446 ms | 496.1 us | 2.508 ms | 1.794 ms | 484.2 us | 0.975x | ok `4b2f43bf0f878408` |
| `fp.p256.mul_karatsuba.d2048` | 18.175 ms | 5.726 ms | 17.429 ms | 14.333 ms | 3.519 ms | 1.043x | ok `b4b89a80c2a3aaa9` |
| `fp.p256.mul_fft.d4096` | 21.965 ms | 2.914 ms | 21.885 ms | 17.133 ms | 3.876 ms | 1.004x | ok `60ed510564513dfb` |
| `fp.p256.mul_fft.d8192` | 44.782 ms | 8.568 ms | 47.310 ms | 35.794 ms | 7.919 ms | 0.947x | ok `47812935df4de5fe` |
| `fp.p256.mul_auto.d2048` | 9.877 ms | 2.569 ms | 9.538 ms | 8.153 ms | 1.835 ms | 1.035x | ok `3d89402081500506` |
| `fp.p256.divmod.d1000_500` | 7.614 ms | 1.819 ms | 6.957 ms | 6.310 ms | 2.912 ms | 1.094x | ok `0bebce51d98eb3f4` |
| `fp.p256.gcd.d300` | 13.572 ms | 3.086 ms | 12.611 ms | 11.132 ms | 2.425 ms | 1.076x | ok `f69e4050cb1a5861` |
| `fp.p256.gcd_ext.d200` | 16.321 ms | 3.544 ms | 15.832 ms | 14.150 ms | 1.537 ms | 1.031x | ok `1a9df499c28b9f49` |
| `fp.p256.pow_mod.d64_k1` | 97.128 ms | 14.391 ms | 91.880 ms | 85.909 ms | 16.190 ms | 1.057x | ok `dbe2b595567fbb9c` |
| `fp.p256.is_irreducible.d12` | 147.863 ms | 38.262 ms | 46.576 ms | 41.665 ms | 7.858 ms | 3.175x | ok `7b839a40c98b0ef0` |
| `fp.p256.factor.d12` | 21.275 ms | 6.609 ms | 21.025 ms | 19.030 ms | 6.002 ms | 1.012x | ok `c33e3512ed388a7e` |
| `fq.q3e30.elem_muladd.n4000` | 6.924 ms | 925.1 us | 6.851 ms | 5.878 ms | 800.8 us | 1.011x | ok `16edfcc06a08028a` |
| `fq.q3e30.elem_inv.n100` | 11.741 ms | 2.249 ms | 848.8 us | 748.1 us | 105.7 us | 13.832x | ok `a8167c01d1cfc5dd` |
| `fq.q3e30.mul_school.d32` | 872.0 us | 132.8 us | 872.4 us | 751.9 us | 180.7 us | 1.000x | ok `ed0ca8eb7d325de3` |
| `fq.q3e30.mul_karatsuba.d128` | 2.593 ms | 705.3 us | 2.671 ms | 2.198 ms | 414.7 us | 0.971x | ok `2464b901bd242c7d` |
| `fq.q3e30.divmod.d64_32` | 1.033 ms | 246.8 us | 907.1 us | 796.8 us | 204.0 us | 1.139x | ok `a2b35f414a139741` |
| `fq.q3e30.factor.d5` | 13.029 ms | 2.814 ms | 7.939 ms | 7.361 ms | 1.801 ms | 1.641x | ok `417a5b753c7a2a07` |
| `fq.q2e8.elem_muladd.n50000` | 10.733 ms | 1.107 ms | 10.472 ms | 9.656 ms | 1.921 ms | 1.025x | ok `c0ffe16ee54f128c` |
| `fq.q2e8.elem_inv.n5000` | 14.234 ms | 2.102 ms | 2.691 ms | 2.469 ms | 367.2 us | 5.289x | ok `81d71b5687e3e735` |
| `fq.q2e8.mul_school.d128` | 719.1 us | 174.4 us | 701.0 us | 626.3 us | 194.8 us | 1.026x | ok `c056d9d992271db9` |
| `fq.q2e8.mul_karatsuba.d512` | 2.567 ms | 630.1 us | 2.663 ms | 2.157 ms | 441.2 us | 0.964x | ok `de1089fdce26d02c` |
| `fq.q2e8.is_irreducible.d16` | 8.336 ms | 1.895 ms | 2.952 ms | 2.578 ms | 583.3 us | 2.824x | ok `9478d590a29a9cbf` |
| `fq.q2e8.factor.d28` | 3.547 ms | 931.1 us | 2.927 ms | 2.695 ms | 545.1 us | 1.212x | ok `d2fc95c10ececedd` |
| `fq.q7e4.elem_muladd.n50000` | 6.922 ms | 2.133 ms | 7.104 ms | 6.429 ms | 1.754 ms | 0.974x | ok `a6cf18a43d40c5c0` |
| `fq.q7e4.elem_inv.n5000` | 13.823 ms | 3.015 ms | 2.436 ms | 1.978 ms | 447.6 us | 5.674x | ok `8664ac6d09409f94` |
| `fq.q7e4.mul_school.d128` | 363.1 us | 76.1 us | 339.3 us | 285.1 us | 67.3 us | 1.070x | ok `71d186560b1d8153` |
| `fq.q7e4.mul_karatsuba.d512` | 2.024 ms | 468.7 us | 2.072 ms | 1.649 ms | 372.3 us | 0.977x | ok `6a2e9b751ab48e99` |
| `fq.q7e4.gcd.d200` | 4.712 ms | 974.9 us | 4.027 ms | 3.530 ms | 895.8 us | 1.170x | ok `8cf19ac90e512feb` |
| `fq.q7e4.factor.d20` | 2.328 ms | 363.0 us | 1.662 ms | 1.477 ms | 168.8 us | 1.401x | ok `ce23318b1b5fd48e` |
| `fq.m61e2.elem_muladd.n100000` | 10.677 ms | 1.495 ms | 10.457 ms | 9.184 ms | 1.687 ms | 1.021x | ok `1814347c960518d8` |
| `fq.m61e2.elem_inv.n1000` | 19.018 ms | 3.350 ms | 454.9 us | 412.2 us | 64.9 us | 41.807x | ok `45326d8095032432` |
| `fq.m61e2.mul_school.d128` | 145.2 us | 30.5 us | 142.7 us | 122.3 us | 28.0 us | 1.018x | ok `2cb38bd639886e3c` |
| `fq.m61e2.mul_karatsuba.d512` | 1.567 ms | 349.5 us | 1.636 ms | 1.253 ms | 272.2 us | 0.958x | ok `f138ad86c0a7a838` |
| `fq.m61e2.mul_fft.d1024` | 19.167 ms | 3.959 ms | 17.668 ms | 15.955 ms | 2.208 ms | 1.085x | ok `a45642d67f5a9912` |
| `fq.m61e2.mul_fft.d4096` | 53.501 ms | 13.356 ms | 51.003 ms | 45.524 ms | 8.188 ms | 1.049x | ok `51844e8f95a9c0ab` |
| `fq.m61e2.mul_auto.d2048` | 32.680 ms | 7.139 ms | 31.546 ms | 29.233 ms | 4.762 ms | 1.036x | ok `76e66ad22fb3fa2e` |
| `fq.m61e2.divmod.d500_250` | 562.7 us | 118.4 us | 550.9 us | 492.0 us | 123.7 us | 1.022x | ok `53509a111940ee2a` |
| `fq.m61e2.gcd.d150` | 3.943 ms | 416.2 us | 1.514 ms | 1.402 ms | 138.2 us | 2.605x | ok `bd8c3c98cc56df84` |
| `fq.m61e2.factor.d7` | 23.252 ms | 4.742 ms | 3.347 ms | 3.159 ms | 974.3 us | 6.948x | ok `d61430d11ee187c3` |
| `fq.build.q3e30` | 581.6 us | 130.0 us | 284.0 us | 257.4 us | 64.8 us | 2.048x | ok `4973a3e567e2fc74` |
| `fq.build.q7e20` | 463.4 us | 106.2 us | 240.2 us | 218.8 us | 42.8 us | 1.929x | ok `9c377cd554aa108e` |
| `fast.ntt.dft.k10` | 1.498 ms | 394.6 us | 1.430 ms | 1.329 ms | 240.4 us | 1.048x | ok `b61f454e63ae6d8d` |
| `fast.ntt.dft.k11` | 1.503 ms | 523.3 us | 1.462 ms | 1.338 ms | 226.7 us | 1.027x | ok `311d4fed0f42034c` |
| `fast.ntt.dft.k12` | 1.531 ms | 215.0 us | 1.480 ms | 1.373 ms | 197.7 us | 1.034x | ok `70691e178cfc0c9d` |
| `fast.ntt.idft.k12` | 1.542 ms | 296.8 us | 1.505 ms | 1.393 ms | 394.0 us | 1.025x | ok `e0677b59c4e33615` |
| `fast.ntt.toeplitz_lower_inverse.n1024` | 713.7 us | 121.2 us | 690.7 us | 636.5 us | 121.2 us | 1.033x | ok `7749eb8ff8cac3d0` |
| `fast.ntt.toeplitz_lower_inverse.n8192` | 7.576 ms | 1.351 ms | 7.714 ms | 6.869 ms | 1.096 ms | 0.982x | ok `0d31554b19d8b996` |
| `fast.ntt.toeplitz_lower_vec.n4096` | 988.8 us | 179.6 us | 909.2 us | 824.8 us | 157.8 us | 1.088x | ok `86e36ad2a1eacb05` |
| `fast.ntt.toeplitz_divmod.d2000_1000` | 3.258 ms | 422.0 us | 3.197 ms | 3.080 ms | 175.9 us | 1.019x | ok `b138cdcceebcc32e` |
| `fast.p7.toeplitz_lower_inverse.n1024` | 1.160 ms | 215.0 us | 1.104 ms | 1.046 ms | 179.3 us | 1.051x | ok `5df08c92a814155a` |
| `fast.p7.toeplitz_divmod.d2000_1000` | 1.906 ms | 211.7 us | 1.772 ms | 1.681 ms | 131.0 us | 1.076x | ok `433d13716b2f95e8` |
