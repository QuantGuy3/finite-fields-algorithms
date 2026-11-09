# Comparativa de benchmarks alcp

Base: **f7i534**. Speedup = mediana(base) / mediana(variante); > 1 significa más rápido.
IQR = rango intercuartílico de las muestras. Las medias geométricas excluyen escenarios con error.

- **f7i534**: modo interleave, 2026-11-09T15:45:45+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -fno-ident -Wl,--no-insert-timestamp`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-93-g0877c55, sha256 1d935c5161da  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: no (batería)
- **f7resto3**: modo interleave, 2026-11-09T15:45:46+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ 15.2.0, GMP 6.3.0, flags `-O3 -fno-ident -Wl,--no-insert-timestamp -Wa,-mbig-obj -Wl,--build-id`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-99-gd3c455e (con cambios sin commit), sha256 ce0110a2d72b  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: no (batería)

Checksums: todos los escenarios coinciden con la base.

## Media geométrica del speedup

| ámbito | f7resto3 |
|---|---|
| **global** | 1.017x (n=25) |
| grupo `fast` | 0.986x (n=1) |
| grupo `fp.m127` | 1.001x (n=4) |
| grupo `fp.m61` | 0.992x (n=3) |
| grupo `fp.ntt` | 1.045x (n=3) |
| grupo `fp.p101` | 1.079x (n=3) |
| grupo `fp.p256` | 1.014x (n=3) |
| grupo `fp.p7` | 1.030x (n=3) |
| grupo `fq.m61e2` | 0.984x (n=2) |
| grupo `fq.q2e8` | 1.003x (n=1) |
| grupo `fq.q3e30` | 0.986x (n=1) |
| grupo `fq.q7e4` | 0.997x (n=1) |
| operación `dft` | 0.986x (n=1) |
| operación `factor` | 0.981x (n=2) |
| operación `gcd_ext` | 1.049x (n=6) |
| operación `mul_karatsuba` | 1.011x (n=16) |

## Detalle por escenario

| escenario | base mediana | base IQR | f7resto3 mediana | f7resto3 mín | f7resto3 IQR | speedup | checksum |
|---|---:|---:|---:|---:|---:|---:|---|
| `fp.p7.mul_karatsuba.d512` | 187.4 us | 27.4 us | 184.2 us | 167.4 us | 34.9 us | 1.017x | ok `d72c59b658ac0379` |
| `fp.p7.mul_karatsuba.d4096` | 2.622 ms | 299.0 us | 2.607 ms | 2.170 ms | 302.9 us | 1.006x | ok `00f3dde93cbd25a4` |
| `fp.p7.gcd_ext.d400` | 842.3 us | 89.2 us | 787.8 us | 669.8 us | 94.2 us | 1.069x | ok `984ac156886a3f82` |
| `fp.p101.mul_karatsuba.d512` | 201.3 us | 39.0 us | 185.6 us | 169.5 us | 36.6 us | 1.085x | ok `fea453dad54c0188` |
| `fp.p101.mul_karatsuba.d4096` | 2.625 ms | 277.2 us | 2.583 ms | 2.118 ms | 272.2 us | 1.016x | ok `82f83f0b7fe7713b` |
| `fp.p101.gcd_ext.d400` | 965.8 us | 78.8 us | 847.9 us | 726.6 us | 91.3 us | 1.139x | ok `4fe0534b04a103a1` |
| `fp.ntt.mul_karatsuba.d512` | 251.7 us | 37.2 us | 247.7 us | 237.3 us | 31.7 us | 1.016x | ok `964e15063fb16df4` |
| `fp.ntt.mul_karatsuba.d4096` | 2.900 ms | 378.0 us | 2.875 ms | 2.446 ms | 261.0 us | 1.009x | ok `d656b678f0c2ca16` |
| `fp.ntt.gcd_ext.d400` | 1.083 ms | 94.6 us | 973.4 us | 861.0 us | 85.7 us | 1.113x | ok `f9dc4e400a7954ee` |
| `fp.m61.mul_karatsuba.d512` | 330.6 us | 43.0 us | 328.5 us | 308.8 us | 39.7 us | 1.006x | ok `a2dcbeaa6574fcb8` |
| `fp.m61.mul_karatsuba.d4096` | 3.285 ms | 505.2 us | 3.250 ms | 2.788 ms | 295.9 us | 1.010x | ok `37d5deaa1cfeb600` |
| `fp.m61.gcd_ext.d400` | 1.657 ms | 124.3 us | 1.726 ms | 1.654 ms | 102.5 us | 0.960x | ok `be3d81e8acc44b48` |
| `fp.m127.mul_karatsuba.d512` | 1.013 ms | 106.3 us | 1.021 ms | 921.2 us | 103.0 us | 0.992x | ok `f773c153fde34172` |
| `fp.m127.mul_karatsuba.d2048` | 3.087 ms | 230.3 us | 3.074 ms | 2.694 ms | 207.6 us | 1.004x | ok `695e0d0f13c14e87` |
| `fp.m127.gcd_ext.d300` | 7.049 ms | 840.6 us | 6.922 ms | 6.441 ms | 871.3 us | 1.018x | ok `45a1b9001453c6ef` |
| `fp.m127.factor.d19` | 12.806 ms | 1.048 ms | 12.927 ms | 11.962 ms | 961.3 us | 0.991x | ok `7b0defab2417737c` |
| `fp.p256.mul_karatsuba.d512` | 4.365 ms | 674.2 us | 4.186 ms | 3.728 ms | 280.5 us | 1.043x | ok `4b2f43bf0f878408` |
| `fp.p256.mul_karatsuba.d2048` | 23.202 ms | 1.716 ms | 23.314 ms | 20.085 ms | 1.631 ms | 0.995x | ok `b4b89a80c2a3aaa9` |
| `fp.p256.gcd_ext.d200` | 30.643 ms | 3.162 ms | 30.459 ms | 29.066 ms | 1.906 ms | 1.006x | ok `1a9df499c28b9f49` |
| `fq.q3e30.mul_karatsuba.d128` | 2.158 ms | 192.3 us | 2.188 ms | 1.851 ms | 227.5 us | 0.986x | ok `2464b901bd242c7d` |
| `fq.q2e8.mul_karatsuba.d512` | 229.6 us | 13.1 us | 228.9 us | 224.1 us | 15.5 us | 1.003x | ok `de1089fdce26d02c` |
| `fq.q7e4.mul_karatsuba.d512` | 1.313 ms | 55.9 us | 1.318 ms | 1.280 ms | 80.0 us | 0.997x | ok `6a2e9b751ab48e99` |
| `fq.m61e2.mul_karatsuba.d512` | 3.005 ms | 341.0 us | 3.011 ms | 2.526 ms | 288.4 us | 0.998x | ok `f138ad86c0a7a838` |
| `fq.m61e2.factor.d7` | 6.998 ms | 844.5 us | 7.210 ms | 6.655 ms | 875.1 us | 0.971x | ok `d61430d11ee187c3` |
| `fast.ntt.dft.k10` | 3.151 ms | 383.5 us | 3.195 ms | 3.023 ms | 469.9 us | 0.986x | ok `b61f454e63ae6d8d` |
