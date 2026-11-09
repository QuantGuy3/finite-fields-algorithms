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

**Posibles regresiones** (por debajo de 0,95x):

- `fp.p256.mul_fft.d4096` 0,88x y `d8192` 0,85x. La desviación típica de esos escenarios es
  del 20 % (3-6 ms sobre 16-39 ms) y la fase 7 no toca el camino multi-limb de F_p con
  p ≥ 2^62. Probablemente es ruido, **sin confirmar**: conviene repetir esos dos escenarios
  con más rondas.
- `fp.p101.mul_karatsuba.d512` 0,945x y `fp.p101.divmod.d2000_1000` 0,949x, dentro del margen
  de ruido habitual.

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
- Comparativa con FLINT (`compare_libs.py`) y cobertura C++ (`build.py cov-cpp`; con las
  variantes nuevas puede acercarse al límite de secciones COFF de MinGW).
- Repetir `fp.p256.mul_fft` para descartar la regresión.
