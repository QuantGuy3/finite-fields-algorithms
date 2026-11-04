# Comparativa con otras librerías: tabla completa

Generado por `bench/informe_comparativa.py` desde `bench/results/comparativa/comparativa.json`. Tiempo mínimo por operación (la construcción del cuerpo y de las entradas no se cronometra).

| Escenario | alcp original | alcp final | alcp final, 1 hilo | alcp paquete | python-flint | galois | SymPy puro | SymPy + FLINT | Resultados iguales |
|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| `fp30.mul.d100` | 561.5 µs | 6.3 µs | 6.3 µs | 6.5 µs | 3.6 µs | 56.9 µs | 5.31 ms | 6.6 µs | sí |
| `fp30.mul.d1000` | 6.47 ms | 120.5 µs | 120.9 µs | 121.3 µs | 68.8 µs | 2.47 ms | 232.89 ms | 72.2 µs | sí |
| `fp30.mul.d10000` | 172.08 ms | 2.97 ms | 4.13 ms | 2.96 ms | 1.69 ms | 251.96 ms | 10.04 s | 1.56 ms | sí |
| `fp30.divmod.d4000` | 431.27 ms | 3.03 ms | 3.03 ms | 3.09 ms | 688.0 µs | 6.42 ms | 3.56 s | 686.8 µs | sí |
| `fp30.gcd.d1000` | 123.68 ms | 4.30 ms | 4.27 ms | 4.29 ms | 582.5 µs | 87.25 ms | 2.79 s | 586.5 µs | sí |
| `fp101.factor.d40` | 8.35 ms | 183.7 µs | 181.1 µs | 186.9 µs | 70.5 µs | 53.33 ms | 8.08 ms | 106.5 µs | sí |
| `fp101.factor.d100` | 479.95 ms | 6.84 ms | 6.66 ms | 6.64 ms | 449.9 µs | 420.29 ms | 242.55 ms | 484.1 µs | sí |
| `fp101.irred.d200` | 20.13 s | 213.63 ms | 207.23 ms | 204.49 ms | 7.67 ms | 286.63 ms | 920.61 ms | 2.50 ms | sí |
| `fp127.mul.d500` | 13.47 ms | 1.41 ms | 3.06 ms | 1.38 ms | 326.2 µs | 16.76 ms | 97.65 ms | 328.9 µs | sí |
| `fp127.factor.d20` | 108.30 ms | 29.79 ms | 30.27 ms | 29.85 ms | 2.56 ms | 2.28 s | 79.15 ms | error | sí |
| `f2e8.mul.d1000` | 1.02 s | 4.39 ms | 15.23 ms | 4.33 ms | 765.4 µs | 2.10 ms | n/s | n/s | sí |
| `f2e8.factor.d30` | 892.00 ms | 13.78 ms | 13.74 ms | 13.78 ms | 8.44 ms | 50.99 ms | n/s | n/s | sí |
| `f3e30.mul.d200` | 854.81 ms | 3.30 ms | 11.93 ms | 3.25 ms | 648.7 µs | 6.85 s | n/s | n/s | sí |
| `f3e30.factor.d5` | 335.73 ms | 9.07 ms | 9.29 ms | 9.25 ms | 2.14 ms | 3.19 s | n/s | n/s | sí |

n/s: la librería no ofrece la operación (SymPy no tiene polinomios sobre F_(p^n)). error: SymPy 1.14 sobre python-flint 0.9 lanza `TypeError` al factorizar con p de 127 bits.

## Cociente frente a python-flint (media geométrica; > 1 = más lento que FLINT)

| Librería | Escenarios | Media geométrica | Mejor | Peor |
|---|---:|---:|---:|---:|
| alcp original | 14 | 248.06x | 41.31x (`fp127.mul.d500`) | 2625.2x (`fp101.irred.d200`) |
| alcp final | 14 | 4.61x | 1.63x (`f2e8.factor.d30`) | 27.9x (`fp101.irred.d200`) |
| alcp final, 1 hilo | 14 | 5.96x | 1.63x (`f2e8.factor.d30`) | 27.0x (`fp101.irred.d200`) |
| alcp paquete | 14 | 4.60x | 1.63x (`f2e8.factor.d30`) | 26.7x (`fp101.irred.d200`) |
| galois | 14 | 109.86x | 2.75x (`f2e8.mul.d1000`) | 10566.7x (`f3e30.mul.d200`) |
| SymPy puro | 10 | 741.76x | 30.89x (`fp127.factor.d20`) | 5930.8x (`fp30.mul.d10000`) |
| SymPy + FLINT | 9 | 0.99x | 0.33x (`fp101.irred.d200`) | 1.8x (`fp30.mul.d100`) |

## Versiones y máquina

- alcp original: baseline-original (67331a2)
- alcp final: rama local
- alcp final, 1 hilo: rama local (ALCP_THREADS=1)
- alcp paquete: rama local
- python-flint: 0.9.0
- galois: 0.4.11
- SymPy puro: 1.14.0
- SymPy + FLINT: 1.14.0 + python-flint 0.9.0
- Máquina: AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0, Python 3.15.1, enchufado, plan de energía «Balanced»
