# Comparativa con otras librerías: tabla completa

Generado por `bench/informe_comparativa.py` desde `bench/results/comparativa/comparativa.json`. Tiempo mínimo por operación (la construcción del cuerpo y de las entradas no se cronometra).

| Escenario | alcp original | alcp fase 3, 1 hilo | alcp nivel A, 1 hilo | alcp nivel B | alcp nivel B, 1 hilo | alcp paquete | python-flint | galois | SymPy puro | SymPy + FLINT | Resultados iguales |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| `fp30.mul.d100` | 565.8 µs | 6.0 µs | 6.0 µs | 6.4 µs | 6.1 µs | 6.7 µs | 3.7 µs | 67.0 µs | 5.37 ms | 6.9 µs | sí |
| `fp30.mul.d1000` | 6.97 ms | 122.3 µs | 122.8 µs | 107.9 µs | 108.0 µs | 108.8 µs | 70.6 µs | 2.55 ms | 261.29 ms | 74.1 µs | sí |
| `fp30.mul.d10000` | 170.11 ms | 4.26 ms | 4.23 ms | 3.21 ms | 4.33 ms | 3.21 ms | 1.82 ms | 263.32 ms | 11.39 s | 1.62 ms | sí |
| `fp30.divmod.d4000` | 534.97 ms | 3.21 ms | 3.27 ms | 2.20 ms | 3.35 ms | 2.21 ms | 706.0 µs | 7.82 ms | 4.74 s | 708.8 µs | sí |
| `fp30.gcd.d1000` | 157.55 ms | 4.53 ms | 4.70 ms | 4.75 ms | 4.74 ms | 4.66 ms | 594.2 µs | 106.25 ms | 3.70 s | 606.1 µs | sí |
| `fp30.irred.d400` | > límite | 5.06 s | 1.52 s | 272.84 ms | 269.40 ms | 259.41 ms | 216.15 ms | 20.85 s | 23.62 s | 58.09 ms | sí |
| `fp30.factor.d300` | 38.50 s | 343.85 ms | 345.21 ms | 133.30 ms | 132.93 ms | 133.70 ms | 23.51 ms | 5.58 s | 25.02 s | 23.56 ms | sí |
| `fp101.factor.d40` | 9.05 ms | 179.9 µs | 183.9 µs | 208.9 µs | 209.8 µs | 214.1 µs | 69.2 µs | 61.42 ms | 9.20 ms | 103.8 µs | sí |
| `fp101.factor.d100` | 610.60 ms | 6.56 ms | 6.68 ms | 4.46 ms | 4.38 ms | 4.38 ms | 450.0 µs | 479.68 ms | 267.49 ms | 481.1 µs | sí |
| `fp101.irred.d200` | 22.31 s | 200.98 ms | 77.47 ms | 27.70 ms | 27.50 ms | 28.12 ms | 7.63 ms | 295.77 ms | 979.66 ms | 2.49 ms | sí |
| `fp127.mul.d500` | 14.55 ms | 3.05 ms | 3.13 ms | 308.3 µs | 303.2 µs | 303.2 µs | 330.9 µs | 18.25 ms | 109.64 ms | 328.9 µs | sí |
| `fp127.factor.d20` | 115.51 ms | 29.60 ms | 29.51 ms | 17.37 ms | 17.53 ms | 17.70 ms | 2.58 ms | 2.50 s | 88.34 ms | error | sí |
| `f2e8.mul.d1000` | 1.10 s | 14.99 ms | 15.11 ms | 1.03 ms | 1.02 ms | 1.03 ms | 762.7 µs | 2.27 ms | n/s | n/s | sí |
| `f2e8.factor.d30` | 961.10 ms | 13.60 ms | 12.63 ms | 7.84 ms | 8.02 ms | 7.84 ms | 8.31 ms | 61.04 ms | n/s | n/s | sí |
| `f3e30.mul.d200` | 915.10 ms | 11.90 ms | 11.91 ms | 949.4 µs | 945.1 µs | 948.8 µs | 521.3 µs | 7.57 s | n/s | n/s | sí |
| `f3e30.factor.d5` | 372.45 ms | 9.01 ms | 5.61 ms | 5.81 ms | 5.78 ms | 5.78 ms | 2.11 ms | 3.40 s | n/s | n/s | sí |

n/s: la librería no ofrece la operación (SymPy no tiene polinomios sobre F_(p^n)). error: SymPy 1.14 sobre python-flint 0.9 lanza `TypeError` al factorizar con p de 127 bits.

## Cociente frente a python-flint (media geométrica; > 1 = más lento que FLINT)

| Librería | Escenarios | Media geométrica | Mejor | Peor |
|---|---:|---:|---:|---:|
| alcp original | 15 | 310.92x | 43.96x (`fp127.mul.d500`) | 2925.5x (`fp101.irred.d200`) |
| alcp fase 3, 1 hilo | 16 | 6.89x | 1.62x (`fp30.mul.d100`) | 26.3x (`fp101.irred.d200`) |
| alcp nivel A, 1 hilo | 16 | 5.87x | 1.52x (`f2e8.factor.d30`) | 22.8x (`f3e30.mul.d200`) |
| alcp nivel B | 16 | 2.57x | 0.93x (`fp127.mul.d500`) | 9.9x (`fp101.factor.d100`) |
| alcp nivel B, 1 hilo | 16 | 2.67x | 0.92x (`fp127.mul.d500`) | 9.7x (`fp101.factor.d100`) |
| alcp paquete | 16 | 2.57x | 0.92x (`fp127.mul.d500`) | 9.7x (`fp101.factor.d100`) |
| galois | 16 | 126.50x | 2.97x (`f2e8.mul.d1000`) | 14519.5x (`f3e30.mul.d200`) |
| SymPy puro | 12 | 718.22x | 34.30x (`fp127.factor.d20`) | 6711.3x (`fp30.divmod.d4000`) |
| SymPy + FLINT | 11 | 0.88x | 0.27x (`fp30.irred.d400`) | 1.9x (`fp30.mul.d100`) |

## Versiones y máquina

- alcp original: baseline-original (67331a2)
- alcp fase 3, 1 hilo: mejoras-cpp-pybind antes de la fase 4 (ALCP_THREADS=1)
- alcp nivel A, 1 hilo: fase4-nivel-a (ALCP_THREADS=1)
- alcp nivel B: rama local
- alcp nivel B, 1 hilo: rama local (ALCP_THREADS=1)
- alcp paquete: rama local
- python-flint: 0.9.0
- galois: 0.4.11
- SymPy puro: 1.14.0
- SymPy + FLINT: 1.14.0 + python-flint 0.9.0
- Máquina: AMD Ryzen 7 7730U with Radeon Graphics (8 núcleos / 16 hilos), Windows-11-10.0.26300-SP0, Python 3.15.1, enchufado, plan de energía «Balanced»
