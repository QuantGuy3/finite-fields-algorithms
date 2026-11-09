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

## Propuesta 3d: tablas NTT con raíz determinista

**Problema.** En `Fp64`, las rutas que no deben consumir el generador aleatorio
(`mul_no_rng`, Newton, `Reducer`, `Composer`, half-GCD) solo podían usar la NTT si la caché
FFT del cuerpo ya tenía tablas, es decir, si la raíz base ya se había sacado del generador.
Desde Python cada cuerpo nace con la caché vacía, así que `divmod` o `gcd` iban por Kronecker
o Karatsuba.

**Solución** (`fft.hpp`: `raiz_det`, `tablas_sin_rng`; `fp64.hpp`: `fft_cache_det_`):

- Esas rutas piden las tablas con `fft::tablas_sin_rng`. Si la caché normal tiene tablas,
  devuelve esas. Si está vacía, construye tablas con una raíz determinista (el primer
  a = 2, 3, … que no es resto cuadrático, elevado a (p−1)/2^cap) y las guarda en una caché
  aparte.
- El producto calculado con la NTT es la convolución exacta módulo p y no depende de la raíz.
  `fft_base_root()`, `tables()` y `fft_cache_` no cambian, así que el generador se consume
  exactamente igual que antes.
- **Riesgo de mezclar raíces.** Una transformada solo puede combinarse con otra hecha con las
  mismas tablas. Un `Reducer` vive durante todo un `factor()`, y entre dos reducciones la caché
  puede llenarse con la raíz del generador. Para evitarlo, `medio::Fija` guarda el puntero a
  sus tablas, y `Reducer` y `Composer` combinan con esas. Las demás rutas piden las tablas una
  vez por operación.

**Medido** con `compare_libs.py --filter fp30` en una sola pasada. La máquina estaba más lenta
que en la comparativa anterior: FLINT tardó unas 2,3 veces más. Por eso se comparan cocientes y
no tiempos absolutos.

| Escenario | tiempo FLINT / tiempo alcp, antes de 3d | ídem, con 3d |
|---|---|---|
| `fp30.divmod.d4000` | 0,58x | **1,73x** |
| `fp30.gcd.d1000` | 0,46x | 0,53x |
| `fp30.factor.d300` | 1,05x | 1,06x |
| `fp30.irred.d400` | 2,25x | 2,29x |

Una cifra mayor que 1 significa que alcp es más rápido. `gcd` apenas cambia: con grado 1000 el
coste está en el Euclides clásico de la cola, no en los productos. JSON en
`bench/results/fase7/comparativa_flint_raizdet.json`.

**Pruebas.** `tests/test_fase7_raizdet.cpp` tiene 43 comprobaciones:

- `raiz_det` tiene orden exactamente 2^cap;
- en frío, `divmod`, `inverse_series`, `gcd` y `Composer` coinciden con las versiones escolares,
  no consumen el generador y dejan `fft_cache_` vacía;
- un `Reducer` y un `Composer` creados en frío siguen siendo correctos después de llenar la
  caché con la raíz del generador.

Además pasan todas las pruebas previas de consumo del generador en frío
(`test_fase6_newton`, `test_fase7_hgcd`, `test_fase6_kronecker`…).

**Revisión adversaria: aprobada.**

- Un programa diferencial compilado contra las cabeceras de antes y de después dio la misma
  salida byte a byte en 117 combinaciones: 13 primos (de 2 a 2^62), 3 semillas y modos frío,
  caliente y mixto. Comparaba resultados, `fft_cache_.empty()` y el estado del generador tras
  cada operación.
- Prueba de mezcla de raíces: 560 comprobaciones, 0 fallos. Su control negativo (una copia
  que combina con las tablas de la caché en vez de `Fija::tab`) da 123 fallos, así que la
  prueba detecta la mezcla.
- Observaciones menores corregidas: `raiz_det` lanza una excepción si la capacidad FFT es 0
  (p = 2), los mensajes de `valida_k` nombran las dos funciones, y se actualizaron los
  comentarios que aún decían «solo con la caché llena».

## Euclides: resto con cociente de grado 1 en una pasada

En el Euclides clásico casi todos los pasos tienen cociente de grado 1. La división escolar
aplicaba dos filas, `a += c1·x·b` y `a += c0·b`, cada una con su reducción módulo p. Con
p < 2^31, `euclides::resto` las aplica ahora en una sola pasada: suma
`a[j] + c1·b[j−1] + c0·b[j]` en 64 bits y reduce una vez con Barrett. Los factores son los
mismos, y cociente y resto son únicos, así que el resultado es idéntico.

Con p ≥ 2^31 se dejó el bucle general. La versión fusionada medía 0,94x en `fp.m61.gcd`, por
la reducción de 128 bits adicional del coeficiente a[db].

Pruebas: sección 4 de `test_fase7_raizdet.cpp` (49 comprobaciones en total). Compara con un
Euclides de referencia hecho con `divmod_school` para p = 3, 101, 998244353, 2^31 − 1,
2^31 + 11 y 2^61 − 1, con coeficientes p − 1 y con b[db−1] = 0.

## Medición de 3d + resto fusionado

Suite intercalada frente a la fase 7 inicial (`f7i534`), 3 rondas × 5 repeticiones, **todos
los checksums iguales**: **1,025x** global. El portátil estaba más lento, probablemente con
batería, pero la comparación intercalada es justa.

| Escenario | Speedup |
|---|---|
| `fp.ntt.divmod.d2000_1000` | 2,38x |
| `fp.p101.gcd.d500` / `fp.ntt.gcd.d500` / `fp.p7.gcd.d500` | 1,59x / 1,60x / 1,36x |
| `gcd_ext` con p < 2^31 | 1,12-1,18x |
| `fp.m61`, `fp.m127`, `fp.p256`, F_q | 0,96-1,04x (sin cambios) |

Algunas multiplicaciones de Karatsuba salieron a 0,82-0,93x, aunque ese código no cambió. Al
repetirlas escenario a escenario (7 rondas, 11 repeticiones, 0,2 s) quedan entre 0,99x y
1,09x: era ruido. Resultados en `bench/results/fase7/fase7-vs-3d-resto.md` y
`fase7-vs-3d-resto_repeticion.md`.

El binario de la suite necesita ahora `-Wa,-mbig-obj` para enlazar, igual que la cobertura (`run_bench.py build` lo añade en Windows):
sin él aparecen referencias indefinidas a símbolos de `std::async`.

## Recalibración de la división de Newton

Desde la última calibración (fase 5) cambiaron tres cosas: la división escolar es más rápida
con p < 2^32 (acumulador de 64 bits de la fase 6), la NTT es más rápida (radix 4) y Newton ya
no necesita la caché llena (3d). Se volvió a ejecutar `bench/calibra_newton.cpp` (Fp64, con
entradas variadas y el mínimo de 5 tandas) con tamaños intermedios. Resultados en
`bench/results/fase7/calibra_newton.txt` y `calibra_newton_fino.txt`.

Umbrales de `prefer_newton` (amortizado / división aislada): el primer tamaño en que Newton
gana a la escolar.

| Cuerpo | Antes | Ahora | Por qué |
|---|---|---|---|
| p < 2^16 (F_101) | 64 / 192 | **80 / 384** | La aislada a 192 era 0,65x frente a la escolar. Con 384 da 1,07x |
| p < 2^32 con NTT (F_998244353) | 128 / 512 | **112 / 224** | La aislada ya gana 1,07x a 224 y 1,44x a 256. La amortizada, 1,32x a 112 |
| p > 2^32 (F_(2^61−1)) | 768 / 1536 | **512 / 1792** | La amortizada gana 1,10x a 512. La aislada era 0,94x a 1536 y da 1,05x a 1792 |

En p < 2^32 sin NTT no se volvió a medir; se deja 192 / 768. `medio::kMin` (64) no cambia:
solo gobierna la recursión de `inverse_series`, y los nuevos umbrales están por encima.

Suite intercalada frente a la versión anterior: **1,011x** global, **checksums iguales**.
`fp.ntt.pow_mod.d128_k5` va **2,24x** más rápido (con 112 el `Reducer` ya usa Newton en
grado 128); el resto queda en ±5 % (`bench/results/fase7/newton-recalibrado.md`).

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

- Recalibrar Newton con 2^16 < p < 2^32 sin NTT (p. ej. 2^31 − 1).
- Acercarse a FLINT en `gcd` con F_998244353 y en `factor`/`irred` con p = 101.
- Informe de cobertura C++. `build.py cov-cpp` ya compila todo (se añadió `-Wa,-mbig-obj`
  a los flags de cobertura: sin él, `bindings.o` a -O0 no enlazaba) y **todos los tests
  instrumentados y pytest pasan**. Pero gcovr 8.7 aborta porque gcov da cuentas negativas en
  los tests largos (`test_fase6_newton`, `test_fase7_fq`; por ejemplo `fp64.hpp:228` y
  `field.hpp:116`), y `--gcov-ignore-parse-errors` no lo evita. Hay que averiguar si son
  contadores desbordados o el error 68080 de gcov.
