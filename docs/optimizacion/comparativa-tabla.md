# Comparativa con otras librerías: tabla completa

Generado por `bench/informe_comparativa.py` desde `bench/results/comparativa/comparativa.json`. Tiempo mínimo por operación (la construcción del cuerpo y de las entradas no se cronometra).

| Escenario | alcp original | alcp fase 3, 1 hilo | alcp fase 4 | alcp fase 4, 1 hilo | alcp paquete | python-flint | galois | SymPy puro | SymPy + FLINT | Resultados iguales |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| `fp30.mul.d100` | 565.8 µs | 7.4 µs | 6.9 µs | 6.8 µs | 7.1 µs | 4.0 µs | 67.0 µs | 5.37 ms | 7.4 µs | sí |
| `fp30.mul.d1000` | 6.97 ms | 143.9 µs | 135.2 µs | 139.3 µs | 134.5 µs | 83.0 µs | 2.55 ms | 261.29 ms | 81.3 µs | sí |
| `fp30.mul.d10000` | 170.11 ms | 4.32 ms | 3.43 ms | 4.61 ms | 3.29 ms | 1.95 ms | 263.32 ms | 11.39 s | 1.73 ms | sí |
| `fp30.divmod.d4000` | 534.97 ms | 3.87 ms | 3.38 ms | 3.23 ms | 3.42 ms | 729.9 µs | 7.82 ms | 4.74 s | 765.6 µs | sí |
| `fp30.gcd.d1000` | 157.55 ms | 5.84 ms | 4.83 ms | 4.89 ms | 4.73 ms | 734.9 µs | 106.25 ms | 3.70 s | 611.6 µs | sí |
| `fp30.irred.d400` | > límite | 5.34 s | 1.75 s | 1.66 s | 1.65 s | 247.34 ms | 20.85 s | 23.62 s | 62.55 ms | sí |
| `fp30.factor.d300` | 38.50 s | 451.02 ms | 397.39 ms | 390.98 ms | 391.63 ms | 25.02 ms | 5.58 s | 25.02 s | 25.94 ms | sí |
| `fp101.factor.d40` | 9.05 ms | 197.0 µs | 200.7 µs | 205.3 µs | 216.9 µs | 76.9 µs | 61.42 ms | 9.20 ms | 117.7 µs | sí |
| `fp101.factor.d100` | 610.60 ms | 8.61 ms | 7.86 ms | 7.91 ms | 8.09 ms | 493.8 µs | 479.68 ms | 267.49 ms | 552.6 µs | sí |
| `fp101.irred.d200` | 22.31 s | 222.20 ms | 93.99 ms | 94.38 ms | 89.46 ms | 8.53 ms | 295.77 ms | 979.66 ms | 2.60 ms | sí |
| `fp127.mul.d500` | 14.55 ms | 3.38 ms | 1.75 ms | 3.31 ms | 1.66 ms | 343.1 µs | 18.25 ms | 109.64 ms | 346.9 µs | sí |
| `fp127.factor.d20` | 115.51 ms | 32.62 ms | 32.75 ms | 31.74 ms | 31.37 ms | 2.69 ms | 2.50 s | 88.34 ms | error | sí |
| `f2e8.mul.d1000` | 1.10 s | 16.52 ms | 5.31 ms | 16.11 ms | 5.32 ms | 867.0 µs | 2.27 ms | n/s | n/s | sí |
| `f2e8.factor.d30` | 961.10 ms | 14.90 ms | 13.53 ms | 13.41 ms | 13.39 ms | 8.55 ms | 61.04 ms | n/s | n/s | sí |
| `f3e30.mul.d200` | 915.10 ms | 12.73 ms | 3.73 ms | 12.69 ms | 3.81 ms | 560.8 µs | 7.57 s | n/s | n/s | sí |
| `f3e30.factor.d5` | 372.45 ms | 9.94 ms | 6.44 ms | 6.01 ms | 5.93 ms | 2.21 ms | 3.40 s | n/s | n/s | sí |

n/s: la librería no ofrece la operación (SymPy no tiene polinomios sobre F_(p^n)). error: SymPy 1.14 sobre python-flint 0.9 lanza `TypeError` al factorizar con p de 127 bits.

## Cociente frente a python-flint (media geométrica; > 1 = más lento que FLINT)

| Librería | Escenarios | Media geométrica | Mejor | Peor |
|---|---:|---:|---:|---:|
| alcp original | 15 | 285.49x | 42.40x (`fp127.mul.d500`) | 2617.2x (`fp101.irred.d200`) |
| alcp fase 3, 1 hilo | 16 | 7.22x | 1.73x (`fp30.mul.d1000`) | 26.1x (`fp101.irred.d200`) |
| alcp fase 4 | 16 | 4.84x | 1.58x (`f2e8.factor.d30`) | 15.9x (`fp101.factor.d100`) |
| alcp fase 4, 1 hilo | 16 | 5.87x | 1.57x (`f2e8.factor.d30`) | 22.6x (`f3e30.mul.d200`) |
| alcp paquete | 16 | 4.78x | 1.57x (`f2e8.factor.d30`) | 16.4x (`fp101.factor.d100`) |
| galois | 16 | 115.79x | 2.61x (`f2e8.mul.d1000`) | 13496.8x (`f3e30.mul.d200`) |
| SymPy puro | 12 | 653.13x | 32.85x (`fp127.factor.d20`) | 6491.5x (`fp30.divmod.d4000`) |
| SymPy + FLINT | 11 | 0.86x | 0.25x (`fp30.irred.d400`) | 1.9x (`fp30.mul.d100`) |

## Versiones y máquina

- alcp original: baseline-original (67331a2)
- alcp fase 3, 1 hilo: mejoras-cpp-pybind antes de la fase 4 (ALCP_THREADS=1)
- alcp fase 4: rama local
- alcp fase 4, 1 hilo: rama local (ALCP_THREADS=1)
- alcp paquete: rama local
- python-flint: 0.9.0
- galois: 0.4.11
- SymPy puro: 1.14.0
- SymPy + FLINT: 1.14.0 + python-flint 0.9.0
- Máquina: AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0, Python 3.15.1, enchufado, plan de energía «Balanced»
