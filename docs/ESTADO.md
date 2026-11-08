# Estado del trabajo

Documento para retomar el trabajo si se interrumpe. El contrato de diseño está en [diseno.md](diseno.md).

## Encargo

Sin tocar los `.py` originales de `python/`, trabajar sobre el port C++ y pybind11:

1. Mejoras propias.
2. Tests: pytest con cobertura >= 99 % y rápidos, más tests C++ donde corresponda.
3. `ejemplos.py`, didáctico y con salida estética.
4. Conversión de tipos generosa y pretty printers.
5. Optimización de velocidad sin cambiar la lógica de los algoritmos, documentada y cuantificada en un documento estético.
6. Git local y en la nube.

## Hecho (rama `mejoras-cpp-pybind`, en GitHub)

| Fase | Contenido | Verificación |
|---|---|---|
| 1. Núcleo C++ | Validaciones, cuerpo común, `division_by_zero`, `prime_power_factor` con `mpz_root`, `factor()` canónico, `evaluate`/`compose`/`is_squarefree`, concepto `FiniteField` | `test_all` 171.162, `test_robustez` 410, `smoke` 15 + 242, `xcheck` 924: 0 fallos |
| 1. Build | `build.py` (tests, oracle, module, pytest, cov-cpp, bench, all, ci), `build.sh`, CMake/CTest, CI de GitHub (Linux + Windows MSYS2), oráculo corregido | `build.py all`: todo correcto |
| 1. Bench | `bench/bench_suite.cpp` (128 escenarios, checksums), `run_bench.py` (build/run/interleave/compare/module), `bench_python.py` | Compila contra las cabeceras originales y las nuevas |
| 2. Binding | `bindings.cpp` compatible y ampliado, conversión de enteros vía bytes, traducción de excepciones | 160 tests en `test_native.py`; `bindings.cpp` cubierto al 100 % |
| 2. Paquete | `python-bindings/alcp/`: `GF`, `K["x"]`, conversión generosa, analizador, pretty print Unicode/ASCII/LaTeX | 692 tests en unos 11 s; cobertura del 100 % en líneas y ramas |
| 2. Ejemplos | `python-bindings/ejemplos.py`: 10 secciones, color ANSI, unas 14 comprobaciones | Unos 4 s |
| 3. Optimización | GMP en sitio, reducción perezosa, backend `Fp64`/`Fq64`, NTT en palabras y multi-limb, Karatsuba sobre vistas, umbrales, paralelismo, capa Python con cachés. Informe en [optimizacion/README.md](optimizacion/README.md) | Suite C++ **24,4x** frente a `baseline-original` (128 checksums iguales); módulo Python 9,7x; `build.py all` correcto, 709 tests de pytest con cobertura del 100 % |
| 3. Comparativa | `bench/compare_libs.py` frente a python-flint, galois y SymPy (14 escenarios, resultados verificados). Informe en [optimizacion/comparativa.md](optimizacion/comparativa.md) | alcp pasa de 248x a 6,0x el tiempo de FLINT (1 hilo, media geométrica); más rápido que galois y SymPy puro en 13 de 14 y 10 de 10 escenarios |

Cobertura C++ (tests C++ + pytest), tras la fase 3: 94,4 % de líneas, 98,4 % de funciones y `bindings.cpp` al 100 %. Antes de la fase 3 era el 98,3 % de líneas; la diferencia está en instanciaciones de los ayudantes paralelos y multi-limb que los tests no recorren.

| 4. Nivel A | División de Newton (`divmod`, `Reducer` en `pow_mod`/`is_irreducible`/`distinct_degree`) con umbrales calibrados, inverso en F_q por Euclides extendido, Ben-Or sin la comprobación final redundante. Informe en [optimizacion/fase4.md](optimizacion/fase4.md) | Suite 1,19x sobre la fase 3 (27,7x frente al original, 128 checksums iguales); `test_fase4` 610 comprobaciones; distancia a FLINT 7,2x → 5,9x |
| 4. Nivel B | Kronecker sobre Z (`kronecker.hpp`), half-GCD en `gcd`/`gcd_ext`, Frobenius por composición modular; umbrales de Newton recalibrados. Informe en [optimizacion/fase4.md](optimizacion/fase4.md) | Suite 1,25x sobre el nivel A (36,8x frente al original, 128 checksums iguales); `test_nivel_b` 31 bloques; distancia a FLINT 5,9x → 2,7x |
| 5 | NTT perezosa sin saltos, Euclides en sitio en Fp64 (también en el caso base de half-GCD) y DDF de Kaltofen–Shoup; despacho, Newton y half-GCD recalibrados. Informe en [optimizacion/fase5.md](optimizacion/fase5.md) | Suite 1,14x sobre el nivel B (128 checksums iguales); `test_fase5` 26 comprobaciones; distancia a FLINT 2,7x → 1,66x (medido con batería) |
| 5. Corrección | `gcd_ext` en Fp64 por encima del umbral de half-GCD: sigue con `hgcd` y termina con el Euclides extendido en sitio (`cola_ext_w`), en lugar de `paso` sobre `Poly`. Detalle en [optimizacion/fase5.md](optimizacion/fase5.md) | De 0,39-0,74x a 1,4-5,1x frente al clásico (3,5-8,5x sobre la versión anterior) en grados 1024-4096; `test_fase5` 27 comprobaciones |

## Pendiente

- Fase 5 hecha (propuestas 1 a 3 de [optimizacion/investigacion-fase5.md](optimizacion/investigacion-fase5.md)); quedan las propuestas 4 a 10.
- Fase 6: propuestas en [optimizacion/investigacion-fase6.md](optimizacion/investigacion-fase6.md). La regresión de `gcd_ext` (punto 0) ya está corregida.
- Half-GCD: acelerar los productos de matrices; con el Euclides en sitio solo gana desde grado 2048-4096 en Fp64.
- Repetir la suite y la comparativa de la fase 5 enchufado (se midieron con batería).
- Investigar el `gcd` de grado 500 en Fp64, a 0,78x en la suite, que no se reproduce aislado.
- Volver a medir la cobertura C++ (`build.py cov-cpp`) tras el nivel B y la fase 5.

- Terminar y probar el asignador rápido de la rama `exp/toolchain-wip` (1,23-1,25x preliminar).
- Recalibrar en otra máquina y en reposo los umbrales de despacho y de paralelismo.
- Revisar `fp.m127.elem_muladd` (~5 % más lento que el original con GMP).
- PGO y AVX2 en la NTT de `Fp64`, sin medir.
- Borrar a mano los restos de `.git/worktrees/med-*` y `wf_627ab9d6-2d7-*` (permisos) y recrear el worktree `scratchpad/baseline`, cuyo árbol está vacío.

## Resueltas en la fase 3

- `fp.ntt.gcd_ext.d400` (el despacho Auto elegía la FFT con un factor de grado 1): 441 → 2,4 ms.
- La FFT a grados pequeños: NTT en palabras de máquina y umbrales de despacho por backend.
