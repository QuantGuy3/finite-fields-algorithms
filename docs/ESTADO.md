# Estado del trabajo (pausado el 2026-10-16)

Documento para retomar el trabajo. El contrato de diseño está en [diseno.md](diseno.md).

## Encargo

Sin tocar los `.py` originales de `python/`, trabajar sobre el port C++ y pybind11:

1. Mejoras propias (corrección de fallos y robustez).
2. Tests: pytest con cobertura >= 99 % y rápidos, más tests C++ donde corresponda.
3. `ejemplos.py` didáctico y con salida estética.
4. Conversión de tipos generosa en la capa Python (enteros, cadenas, listas...) y pretty printers.
5. Optimización de velocidad sin cambiar la lógica de los algoritmos, documentada y
   cuantificada (qué funcionó, qué no y cuánto) en un documento estético.
6. Git local y en la nube (GitHub).

## Hecho

- Rama `mejoras-cpp-pybind` con:
  - `docs/diseno.md`: contrato de API y reparto de tests;
  - `.gitignore` ampliado;
  - `.gitattributes`: LF en el repositorio, excepto `python/`, que se conserva byte a byte.
- Etiqueta `baseline-original` en el commit original (`67331a2`), que es la línea base de las mediciones.
- Entorno local:
  - `.venv` con pybind11 3.1, pytest 9, pytest-cov, pytest-xdist e hypothesis;
  - g++ 15.2 (MSYS2 ucrt64) y GMP 6.3.
- Línea base medida con el `bench/bench.cpp` original (`-O3 -march=native`, Ryzen 7 7730U). Es una sola ejecución y orientativa:

| Escenario | Tiempo |
|---|---|
| Construir F_(3^30) | 67 ms |
| Factorizar un polinomio de grado 5 sobre F_(3^30) | 823 ms |
| Multiplicar dos polinomios de grado 4000 en F_998244353, escolar | 3752 ms |
| ídem, Karatsuba | 606 ms |
| ídem, FFT | 62 ms |
| Factorizar un polinomio de grado 48 sobre F_7 | 12 ms |

- Tests originales: `test_all` (99.634 comprobaciones) y `crosscheck` (704) pasan.

## Trabajo parcial de la fase 1 (ramas `wip/`, SIN VERIFICAR)

La fase 1 (3 agentes en paralelo) se detuvo a medias para ahorrar créditos:

| Rama | Contenido | Estado |
|---|---|---|
| `wip/fase1-nucleo-cpp` | `errors.hpp` (`division_by_zero`) y primeros cambios en `bigint.hpp` / `fp.hpp` | Empezado. No compilado ni probado |
| `wip/fase1-build` | `oracle.py` (rutas a `../python`, `py_ok` normalizado) y `showcase.py` (imports) | Empezado. Sin `build.py` ni CI |
| bench | Nada escrito todavía | Por hacer |

## Próximos pasos

1. Terminar la fase 1 según `docs/diseno.md` §4 y el script del workflow de la fase 1:
   - núcleo C++ y tests C++ (`include/`, `tests/`);
   - `build.py`, `build.sh`, CMake, `pyproject.toml`, CI y oráculo;
   - harness de benchmarks (`bench/bench_suite.cpp`, `run_bench.py`, `bench_python.py`), compilable también contra las cabeceras de `baseline-original`.
2. Integrar las ramas en `mejoras-cpp-pybind`, ejecutar todos los tests y medir la línea base en reposo.
3. Fase 2:
   - `bindings.cpp`: API compatible y ampliada, conversión rápida de enteros, traducción de excepciones;
   - paquete `python-bindings/alcp/` con conversión de tipos y pretty print;
   - pytest con >= 99 % de cobertura;
   - `ejemplos.py`.
4. Fase 3: optimizaciones medidas, por ejemplo:
   - flags `-O3`, `-march`, LTO y PGO;
   - operaciones GMP en sitio y reducción perezosa;
   - backend de 64 bits (Montgomery) para p pequeño;
   - caché y estructura de la FFT, paralelismo, asignador de memoria;
   - conversión de enteros en el binding;
   - umbrales de despacho.
5. Informe final estético y push.

## Defectos conocidos (de la revisión inicial) que quedan por corregir

- `build.sh`: `run_exe` siempre devuelve 0; rutas a `python/` incorrectas.
- `bindings.cpp:64`: `from_coeffs` mueve el `shared_ptr` y lo desreferencia en la misma llamada.
- Sin comprobación de cuerpo común en las operaciones de `Poly` (comportamiento indefinido con `Fq` distintos).
- Bucles o recursiones infinitas y accesos fuera de rango:
  - `toeplitz_*` sin validar (n = 0);
  - `multiplicity` con divisor constante;
  - `radical` con constantes;
  - `Fq(p, 0)` y `Fq` con módulo vacío.
- `prime_power_factor` con división por tentativa: no termina para p grande.
- `crosscheck` no elimina `\r`; con `-DNDEBUG` los smoke tests no comprueban nada.
