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

Cobertura C++ (tests C++ + pytest), tras la fase 3: 94,4 % de líneas, 98,4 % de funciones y `bindings.cpp` al 100 %. Antes de la fase 3 era el 98,3 % de líneas; la diferencia está en instanciaciones de los ayudantes paralelos y multi-limb que los tests no recorren.

## Pendiente

- Terminar y probar el asignador rápido de la rama `exp/toolchain-wip` (1,23-1,25x preliminar).
- Recalibrar en otra máquina y en reposo los umbrales de despacho y de paralelismo.
- Revisar `fp.m127.elem_muladd` (~5 % más lento que el original con GMP).
- PGO y AVX2 en la NTT de `Fp64`, sin medir.
- Borrar a mano los restos de `.git/worktrees/med-*` y `wf_627ab9d6-2d7-*` (permisos) y recrear el worktree `scratchpad/baseline`, cuyo árbol está vacío.

## Ideas apuntadas por los agentes para la fase 3

- `fp.ntt.gcd_ext.d400` es unas 10 veces más lento que en F_7: el despacho Auto elige FFT aunque un factor sea de grado 1.
- La FFT apenas gana a grados pequeños (×2,4 en grado 256), por un sobrecoste fijo.
- En GF(2^127) una multiplicación tarda unos 0,7 ms y `minimal_polynomial` unos 8,7 s.
- El inverso en F_q por Fermat cuesta Θ(n³ log p).
