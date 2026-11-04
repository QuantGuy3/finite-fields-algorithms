# Comparativa de benchmarks alcp

Base: **py-final-ref**. Speedup = mediana(base) / mediana(variante); > 1 significa más rápido.
IQR = rango intercuartílico de las muestras. Las medias geométricas excluyen escenarios con error.

- **py-final-ref**: modo interleave, 2026-11-04T19:01:26+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ (Rev8, Built by MSYS2 project) 15.2.0, GMP ?, flags `-O3`  
  include `~\AppData\Local\Temp\claude\C--Users-<usuario>-Documents-Archivos-Proyectos\2e18e96a-361a-47be-ad99-f0341eb0cc55\scratchpad\fin\include` @ ?, sha256 9b306a5d8af5  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: sí
- **py-final**: modo interleave, 2026-11-04T19:01:38+01:00  
  CPU AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0  
  compilador g++ (Rev8, Built by MSYS2 project) 15.2.0, GMP ?, flags `-O3`  
  include `~\Documents\Archivos\Proyectos\finite-fields-algorithms\include` @ baseline-original-62-g32d6922, sha256 af0c61d3df27  
  energía: Power Scheme GUID: 381b4222-f694-41f0-9685-ff5bb260df2e  (Balanced); corriente: sí

Checksums: todos los escenarios coinciden con la base.

## Media geométrica del speedup

| ámbito | py-final |
|---|---|
| **global** | 9.701x (n=19) |
| grupo `py.fq` | 12.920x (n=3) |
| grupo `py.m521` | 2.593x (n=4) |
| grupo `py.misc` | 1.103x (n=1) |
| grupo `py.ntt` | 16.624x (n=10) |
| grupo `py.p101` | 32.382x (n=1) |
| operación `coeffs` | 1.982x (n=2) |
| operación `dft` | 8.752x (n=1) |
| operación `factor` | 32.289x (n=2) |
| operación `fp_mul_loop` | 2.939x (n=2) |
| operación `is_prime` | 1.103x (n=1) |
| operación `mul_auto` | 17.878x (n=3) |
| operación `mul_fft` | 42.810x (n=1) |
| operación `mul_karatsuba` | 87.193x (n=1) |
| operación `mul_school` | 144.929x (n=1) |
| operación `poly_from_list` | 3.366x (n=3) |
| operación `roundtrip` | 4.249x (n=1) |
| operación `toeplitz_lower_inverse` | 32.813x (n=1) |

## Detalle por escenario

| escenario | base mediana | base IQR | py-final mediana | py-final mín | py-final IQR | speedup | checksum |
|---|---:|---:|---:|---:|---:|---:|---|
| `py.ntt.poly_from_list.d1000` | 5.395 ms | 338.9 us | 950.0 us | 898.1 us | 105.0 us | 5.679x | ok `8627fe77db0d9194` |
| `py.ntt.coeffs.d1000` | 6.059 ms | 256.8 us | 2.067 ms | 1.972 ms | 141.6 us | 2.932x | ok `0fc88e478ee1e3c7` |
| `py.ntt.roundtrip.d1000` | 8.457 ms | 353.1 us | 1.990 ms | 1.902 ms | 85.8 us | 4.249x | ok `aa3b36944fe8dc0a` |
| `py.ntt.mul_auto.d1000` | 34.234 ms | 1.912 ms | 782.4 us | 756.0 us | 22.0 us | 43.756x | ok `0d4a33cb5f0e28e9` |
| `py.ntt.mul_fft.d1000` | 33.561 ms | 1.183 ms | 784.0 us | 756.6 us | 23.8 us | 42.810x | ok `368464de3541086c` |
| `py.ntt.mul_karatsuba.d1000` | 70.757 ms | 2.375 ms | 811.5 us | 780.4 us | 30.4 us | 87.193x | ok `6776e3c212b85eaa` |
| `py.ntt.mul_school.d1000` | 108.335 ms | 3.927 ms | 747.5 us | 722.3 us | 31.8 us | 144.929x | ok `eebc7f1a8dc65835` |
| `py.p101.factor.d40` | 10.935 ms | 3.586 ms | 337.7 us | 248.2 us | 90.4 us | 32.382x | ok `19d2f8b8f63708d8` |
| `py.ntt.fp_mul_loop.n20000` | 19.895 ms | 790.0 us | 5.935 ms | 5.660 ms | 399.6 us | 3.352x | ok `9bb1099e4af5fb76` |
| `py.m521.fp_mul_loop.n5000` | 28.024 ms | 1.328 ms | 10.876 ms | 10.465 ms | 881.3 us | 2.577x | ok `a33d46746e9ee96f` |
| `py.m521.poly_from_list.d200` | 8.857 ms | 360.1 us | 2.236 ms | 2.109 ms | 138.0 us | 3.962x | ok `387f4c5743d4f4dd` |
| `py.m521.coeffs.d200` | 5.650 ms | 185.7 us | 4.216 ms | 3.996 ms | 215.3 us | 1.340x | ok `d070fe4d50d6532f` |
| `py.m521.mul_auto.d200` | 5.651 ms | 231.5 us | 1.710 ms | 1.644 ms | 112.3 us | 3.305x | ok `fac612a2571258f9` |
| `py.is_prime.n2000` | 10.731 ms | 226.2 us | 9.731 ms | 9.305 ms | 443.9 us | 1.103x | ok `fac62a4469352205` |
| `py.fq9.poly_from_list.d500` | 6.769 ms | 340.6 us | 3.992 ms | 3.796 ms | 240.4 us | 1.696x | ok `439133366a0b2497` |
| `py.fq9.mul_auto.d500` | 55.913 ms | 1.342 ms | 1.415 ms | 1.094 ms | 183.4 us | 39.512x | ok `c8af7f9ed94034c4` |
| `py.fq.q3e30.factor.d5` | 295.800 ms | 66.510 ms | 9.188 ms | 7.487 ms | 3.330 ms | 32.196x | ok `92e0b4a372f8a109` |
| `py.ntt.dft.k10` | 11.096 ms | 381.4 us | 1.268 ms | 1.202 ms | 74.4 us | 8.752x | ok `4534ca8416019033` |
| `py.ntt.toeplitz_lower_inverse.n1024` | 25.088 ms | 688.6 us | 764.6 us | 726.9 us | 36.5 us | 32.813x | ok `2eafed6cdc613c27` |
