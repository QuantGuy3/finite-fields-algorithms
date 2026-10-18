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

Cobertura C++ (tests C++ + pytest): 98,3 % de líneas y 100 % de funciones.

## Pendiente

- Fase 3: optimización de velocidad, medida con `bench/run_bench.py interleave` contra `baseline-original` (-O3):
  - flags (`-march`, LTO, PGO);
  - GMP en sitio y reducción perezosa;
  - backend de 64 bits para p pequeño;
  - FFT (cuadrado, tablas), umbrales de despacho, asignador de memoria, paralelismo.
- Informe final estético y README actualizado.

## Ideas apuntadas por los agentes para la fase 3

- `fp.ntt.gcd_ext.d400` es unas 10 veces más lento que en F_7: el despacho Auto elige FFT aunque un factor sea de grado 1.
- La FFT apenas gana a grados pequeños (×2,4 en grado 256), por un sobrecoste fijo.
- En GF(2^127) una multiplicación tarda unos 0,7 ms y `minimal_polynomial` unos 8,7 s.
- El inverso en F_q por Fermat cuesta Θ(n³ log p).
