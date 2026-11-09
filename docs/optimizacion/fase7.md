# Fase 7: half-GCD V3e, `FpN<2>`, NTT radix 4 con AVX2 y FqTab/FqSwar

Fecha: 9 de noviembre de 2026. Rama `fase7`, sobre `mejoras-cpp-pybind` en `80b983d` (fase 6).
Implementa las propuestas 2, 3, 4 y 7 de [investigacion-fase6.md](investigacion-fase6.md),
que la fase 6 había dejado pendientes.

## Resumen

| # | Cambio | Dónde | Medido (suite intercalada) |
|---|---|---|---|
| 2 | Half-GCD V3e: resto y aplicación de la matriz en palabras, producto de matrices con NTT cíclica (con recuento de acarreos), segunda llamada recursiva sin copias y nuevos umbrales (960/480/256/256 con NTT) | `polynomial.hpp` (`hgcd_rem_w`, `aplicar_w`, `mul_mat_w`, `segunda_w`, `umbral_sigue_gcd`, `hgcd_con_ntt`) | `gcd_ext` 1,6-1,85x en F_7, F_101 y F_998244353 (grado 400); el agente midió 1,4-2,5x en `gcd` y `gcd_ext` desde grado 2048 |
| 3 | `FpN<2>`: F_p con 2^63 ≤ p < 2^128 en dos palabras con Montgomery, sin GMP; backend `'u128'` en el binding | `include/alcp/fpn.hpp`, `polynomial.hpp` (ramas de umbral), `python-bindings/bindings.cpp` (`Var = variant<Fp, Fp64, FpN2>`) | grupo `fp.m127` **3,59x**; `gcd` 7,0x, `gcd_ext` 4,5x |
| 4 | Convolución NTT radix 4 (`dif_r4`/`dit_r4`) con producto punto a punto de Montgomery fusionado con n^-1, y camino AVX2 de 32 bits para p < 2^30, n ≥ 16 | `fft.hpp`, `fft_tables.hpp` (`Word32NTTTables`, `build_tables32`) | grupo `fp.ntt` 1,34x; `toeplitz_lower_inverse` 3,5-5,5x |
| 7 | Extensiones pequeñas: `FqTab` (tablas de log/antilog, q < 2^16) y `FqSwar` (coeficientes empaquetados en una palabra) | `include/alcp/fq_tab.hpp`, `fq_swar.hpp`, `polynomial.hpp` (`kRepresentaFq64`), binding (`FqW::Var`) | F_(2^8) **45x**; F_(7^4) 10x; F_(3^30) 2,5x |

**Suite C++** (128 escenarios, intercalada con la fase 6, 3 rondas × 5 repeticiones): **1,72x**
de media geométrica global y **todos los checksums iguales** a la base. La base es el binario
`f6-desd3`, que tiene exactamente las cabeceras de `80b983d`, recompilado con otros flags del
enlazador porque Smart App Control bloqueaba el ejecutable original.

## Método

Se siguió el mismo de la fase 6:

1. **Implementación en paralelo.** Un agente por propuesta, cada uno en su git worktree.
2. **Revisión adversaria** independiente de cada rama. Las propuestas 2, 4 y 7 se aprobaron con
   observaciones menores. La 3 no, por el fallo descrito abajo; se corrigió en la integración.
3. **Integración en serie** en `fase7`: fusiones de `fase7-p2-hgcd`, `fase7-p4-ntt`,
   `fase7-p3-fpn` y `fase7-p7-fq`, resolución de conflictos (`bench_suite.cpp`,
   `prefer_kronecker`, `visitar2` genérico del binding) y `build.py all`.
4. **Medición** con `bench/run_bench.py interleave` en la máquina sin otros procesos pesados.

## Dos correcciones de la integración

- **FpN2 y el generador aleatorio (revisión).** Un producto pequeño (por ejemplo 1×1) con
  `algo='fft'` en un cuerpo `FpN2` recién creado no pedía las tablas FFT, así que la raíz se
  sacaba del generador en otro momento que con GMP y las secuencias aleatorias posteriores
  divergían. Ahora `MulAlgo::FFT` llama a `fft::tables` antes de `convolve` cuando el cuerpo es
  `FpN2`. Prueba: `test_fft_corto_frio` en `tests/test_fase7_fpn.cpp`.
- **`backend` de los cuerpos de extensión (fusión).** Tras añadir `FpN2`, la propiedad
  `backend` de `Fq` devolvía `'u128'`: fallaban 11 pruebas de pytest. Ahora se decide por el
  tipo del cuerpo (`gmp` para `Fp`/`Fq`, `u128` para `FpN2`, `u64` para el resto).

## Resultados

Tabla completa en [bench/results/fase7/base-vs-fase7.md](../../bench/results/fase7/base-vs-fase7.md);
JSON en `bench/results/fase7/suite-base.json` y `suite-fase7.json`.

| Grupo | Speedup | | Operación | Speedup |
|---|---|---|---|---|
| `fq.q2e8` | 45,39x | | `elem_inv` | 3,32x |
| `fq.q7e4` | 10,19x | | `elem_muladd` | 2,56x |
| `fp.m127` | 3,59x | | `is_irreducible` | 2,12x |
| `fq.q3e30` | 2,54x | | `factor` | 2,09x |
| `fast` | 1,48x | | `gcd` / `gcd_ext` | 1,76x / 1,73x |
| `fp.ntt` | 1,34x | | `mul_school` / `mul_karatsuba` | 1,73x / 1,66x |
| `fp.p7` / `fp.p101` | 1,05x / 1,03x | | `mul_fft` | 1,60x |
| `fp.m61`, `fq.m61e2`, `fq.build` | 0,99-1,00x | | `divmod` | 1,18x |
| `fp.p256` | 0,98x | | global | **1,72x** |

**Regresiones aparentes descartadas.** En la suite salieron por debajo de 0,95x
`fp.p256.mul_fft.d4096` (0,88x), `fp.p256.mul_fft.d8192` (0,85x, con un 20 % de desviación
típica), `fp.p101.mul_karatsuba.d512` (0,945x) y `fp.p101.divmod.d2000_1000` (0,949x). Al
repetirlos escenario a escenario (7 rondas, al menos 15 repeticiones y 0,3 s cada uno) quedan
los cuatro en **0,99x**, con los checksums iguales: era ruido
([repeticion_regresiones.txt](../../bench/results/fase7/repeticion_regresiones.txt)).

## Comparativa con FLINT

`python bench/compare_libs.py --libs alcp_native,flint` con python-flint 0.9.0, binding de la
fase 7. Los tiempos son medianas y los checksums coinciden en todos los escenarios. La columna
«alcp antes» es la de [comparativa.md](comparativa.md) (7 de noviembre, antes de las fases 6 y 7).
JSON en `bench/results/fase7/comparativa_flint.json`. Una cifra mayor que 1 en la última
columna significa que alcp es más rápido.

| Escenario | alcp antes | alcp fase 7 | FLINT | FLINT / alcp fase 7 |
|---|---|---|---|---|
| `fp30.mul.d100` | 6.5 µs | 1.5 µs | 4.5 µs | 3.00x |
| `fp30.mul.d1000` | 110.1 µs | 9.0 µs | 70.3 µs | 7.81x |
| `fp30.mul.d10000` | 3.63 ms | 272.7 µs | 1.77 ms | 6.50x |
| `fp30.divmod.d4000` | 2.66 ms | 1.21 ms | 695.2 µs | 0.58x |
| `fp30.gcd.d1000` | 5.36 ms | 1.28 ms | 592.0 µs | 0.46x |
| `fp30.irred.d400` | 290.26 ms | 96.75 ms | 218.05 ms | 2.25x |
| `fp30.factor.d300` | 133.82 ms | 22.48 ms | 23.50 ms | 1.05x |
| `fp101.factor.d40` | 222.4 µs | 99.2 µs | 70.9 µs | 0.71x |
| `fp101.factor.d100` | 4.55 ms | 792.0 µs | 454.9 µs | 0.57x |
| `fp101.irred.d200` | 28.30 ms | 10.48 ms | 7.98 ms | 0.76x |
| `fp127.mul.d500` | 311.7 µs | 255.3 µs | 333.2 µs | 1.31x |
| `fp127.factor.d20` | 17.56 ms | 1.52 ms | 2.69 ms | 1.77x |
| `f2e8.mul.d1000` | 1.07 ms | 250.5 µs | 807.1 µs | 3.22x |
| `f2e8.factor.d30` | 8.03 ms | 187.5 µs | 8.62 ms | 45.96x |
| `f3e30.mul.d200` | 971.9 µs | 593.5 µs | 533.2 µs | 0.90x |
| `f3e30.factor.d5` | 6.22 ms | 1.17 ms | 2.15 ms | 1.84x |

alcp supera a FLINT en el producto, la irreducibilidad en F_998244353, F_(2^127−1) y las
extensiones pequeñas. Sigue por detrás en `divmod` (0,58x) y `gcd` (0,46x) en
F_998244353 y en `factor`/`irred` con p = 101 (0,57-0,76x). Esas son las candidatas para la
recalibración de Newton y de los umbrales de half-GCD que quedan pendientes. Es una sola
pasada sin intercalar, así que diferencias menores del 10 % no son significativas.

## Pruebas

- `tests/test_fase7_*.cpp` (hgcd, NTT, FpN2 con 3558 comprobaciones, FqTab/FqSwar) y las
  pruebas Python correspondientes.
- `build.py all`: «TODO CORRECTO» (16 tests C++, oráculo y 837 pruebas de pytest). Los avisos
  `FPX_FACTOR`/`FQX_FACTOR` del oráculo vienen de la referencia Python y ya existían antes.

## Observaciones de la revisión no resueltas

- `FqTab` reutiliza `kFrobPresupuesto` y `kParKaratsubaMin` de `Fp64` sin calibrar.
- Con q cercano a 2^16, `FqTab` reserva unos 768 KiB de tablas por cuerpo.
- `fft::convolve` no comprueba que las entradas estén en forma canónica.
- El cambio de `repr` de los cuerpos con p de 64-128 bits (ahora `FpN2`) cuando se pide
  `backend='gmp'` debe documentarse en la referencia del binding.
- La convolución AVX2 se desactiva a -O0, así que la cobertura C++ no la ejecuta.

## Pendiente

- Caché de raíz determinista (propuesta 3d) y recalibración de `prefer_newton` y `medio::kMin`
  con las nuevas NTT.
- Acercarse a FLINT en `divmod`, `gcd` y `factor` con p pequeño (ver la comparativa).
- Cobertura C++ (`build.py cov-cpp`; con las
  variantes nuevas puede acercarse al límite de secciones COFF de MinGW).
