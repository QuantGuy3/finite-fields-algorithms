# Diseño de la mejora del port C++ / pybind11

Documento de referencia (contrato) para el trabajo sobre el port C++23 y su
binding Python. **Los `.py` originales de `python/` (las tres prácticas) no se
modifican bajo ningún concepto**: se usan solo como referencia/oráculo.

## 1. Estructura del repositorio

```
include/alcp/            Núcleo C++23 header-only (GMP)
tests/                   Tests C++ (arnés CHECK de tests/check.hpp; sin assert)
bench/                   Benchmarks C++ y Python + resultados (bench/results/)
python/                  Prácticas Python ORIGINALES (intocables)
python-bindings/
  bindings.cpp           Binding pybind11 -> python-bindings/alcp/alcp_native<EXT_SUFFIX>
  alcp/                  Paquete Python (capa de alto nivel sobre alcp_native)
  tests/                 Suite pytest (cobertura >= 99 % del paquete alcp)
  ejemplos.py            Recorrido didáctico por la librería
  oracle.py              Genera vectores de referencia con las prácticas originales
  showcase.py, compare_nemo.jl
build.py                 Driver de compilación multiplataforma (Windows MinGW + Linux)
build.sh                 Envoltorio: exec python build.py "$@"
CMakeLists.txt           Tests C++ con CTest
pyproject.toml           Configuración de pytest y coverage
.github/workflows/ci.yml Integración continua (GitHub Actions)
docs/                    Documentación (este archivo, informe de optimización)
```

Todo lo generado va a `build/` (ignorado por git), salvo el módulo compilado,
que se coloca dentro del paquete (`python-bindings/alcp/alcp_native*.pyd|so`,
también ignorado).

## 2. Reparto de tests

| Qué se prueba | Dónde | Herramienta |
|---|---|---|
| Núcleo C++ (aritmética, polinomios, factorización, FFT, Toeplitz, validaciones) | `tests/*.cpp` | arnés `CHECK`, CTest |
| Paridad C++ vs prácticas Python originales | `tests/crosscheck.cpp` + `oracle.py` | vectores generados en `build/` |
| Binding y capa Python (conversión de tipos, parser, pretty print, errores) | `python-bindings/tests/` | pytest + pytest-cov (>= 99 % líneas y ramas) |
| Cobertura C++ alcanzada por pytest y por los tests C++ | `build.py cov-cpp` | gcov + gcovr (informativo) |

Los tests deben ser rápidos: la suite pytest completa en < 15 s y los tests C++
en < 10 s en el portátil de desarrollo.

## 3. Convenciones de código

- Comentarios nuevos **en español**. Se conserva el estilo del código existente.
- C++: no se rompe la API existente (solo se añade). Las cabeceras siguen siendo
  header-only y compilan con `-std=c++23 -Wall -Wextra` sin avisos.
- Errores C++:
  - `alcp::division_by_zero` (deriva de `std::domain_error`): invertir 0,
    dividir por el polinomio cero. pybind11 lo traduce a `ZeroDivisionError`.
  - `std::invalid_argument`: argumentos inválidos (p no primo, módulo
    reducible, longitudes incorrectas, operandos de cuerpos distintos...).
    pybind11 lo traduce a `ValueError`.
  - `std::domain_error` (otros): p. ej. polinomio no invertible módulo m -> `ValueError`.
- Python: docstrings y mensajes de error en español; `__slots__` en las clases de valor.

## 4. Mejoras del núcleo C++ (fase 1)

Corrección de los defectos de la revisión y robustez, **sin cambiar la lógica de
los algoritmos**:

1. Comprobación de cuerpo común en todas las operaciones binarias de `Poly`
   (`same_field`: mismo puntero o mismos parámetros; si no, `invalid_argument`).
2. `prime_power_factor` por detección de potencia perfecta con GMP (`mpz_root`),
   O(log q) raíces en lugar de división por tentativa.
3. Validación previa en `Fq` (módulo vacío, `n = 0`), en `fast.hpp` (n, longitudes,
   orden exacto de la raíz en `dft`), en `fft.hpp` (precondiciones de tamaño),
   en `multiplicity` (divisor nulo o constante), `radical` (constantes),
   `distinct_degree` / `equal_degree` (precondiciones).
4. `factor()` devuelve los factores en orden canónico (grado y coeficientes).
5. Añadidos: `Poly::evaluate` (Horner), `Poly::compose`, `Poly::is_squarefree`,
   `Fq::from_vector`, concepto C++20 `FiniteField` en `field.hpp`,
   `Poly(const Field&&) = delete` para evitar punteros colgantes.
6. Comentarios corregidos (el test de irreducibilidad es el de Ben-Or, etc.).

## 5. Módulo nativo `alcp_native` (fase 2)

Superconjunto **compatible** con la API original (los scripts antiguos siguen
funcionando), más:

- Conversión entero Python <-> `mpz` rápida y sin límite de dígitos
  (camino rápido para enteros de 64 bits; bytes con `mpz_import`/`mpz_export`
  para el resto; nada de cadenas decimales).
- `set_seed(s)`; traducción de excepciones (`ZeroDivisionError`, `ValueError`).
- `Fp` y `Fq`: igualdad y hash por parámetros; operaciones de elemento
  completas (`add, sub, neg, mul, div, inv, pow, random`), en `Fp` sobre enteros
  y en `Fq` sobre **índices en base p**; `Fq.from_coeffs(list) -> índice`,
  `Fq.to_coeffs(índice) -> list`; `fft_capacity()`, `fft_base_root()`.
- `PolyFp` / `PolyFq`: lo existente + `evaluate`, `compose`, `radical`,
  `distinct_degree`, `equal_degree`, `is_squarefree`, `leading`, `coeff`,
  `scalar_mul`; los operandos de cuerpos distintos lanzan `ValueError`.
- Algoritmos rápidos (`dft, idft, toeplitz_*`, `divmod_toeplitz`) para `Fp` **y** `Fq`.

## 6. Paquete Python `alcp` (fase 2)

Capa de alto nivel, con conversión de tipos generosa y representación cuidada:

```python
import alcp
F = alcp.GF(7)                      # cuerpo primo
K = alcp.GF(9, modulus="x^2+1")     # extensión (sin módulo: el primer irreducible en un orden fijo, determinista; modulus="random" para uno aleatorio)
a = K.gen                           # clase de x, se imprime 'a'
R = K["x"]                          # anillo de polinomios (también alcp.PolynomialRing(K, "x"))
x = R.gen
f = R("(a+1)x^3 + 2x + a")          # cadenas, listas, dicts, enteros, elementos...
g = f + 2                           # 2 se interpreta como 2·1
h = f * "x + a"                     # cadenas también
print(f.factor())                   # (a + 1) · (x + …)² · …
```

Reglas de conversión (de más a menos específica):

| Entrada | En un cuerpo | En un anillo de polinomios |
|---|---|---|
| `int`, `numpy.integer`, `bool` (vía `operator.index`) | n·1 | constante n·1 |
| `fractions.Fraction` | num · den⁻¹ | constante |
| `str` | se analiza (`"a^2+1"`, `"2*a"`, `"a²"`) | se analiza (`"x^3 + (a+1)x - 2"`) |
| `list` / `tuple` | coeficientes sobre F_p (orden creciente) | coeficientes (cada uno convertido) |
| `dict {exp: coef}` | — | polinomio disperso |
| elemento del subcuerpo primo | se sumerge | constante |
| `float`, otros | `TypeError` con mensaje claro | `TypeError` |

Representación: `str` en Unicode (`x³ + 2x + 1`), `repr` inequívoco,
`_repr_latex_` y `_repr_pretty_` para Jupyter/IPython, especificadores de
formato (`f"{f:ascii}"`, `f"{f:latex}"`, `f"{f:unicode}"`) y estilo global
configurable (`alcp.set_print_style`).

## 7. Optimización de velocidad (fase 3)

Sin tocar la lógica de los algoritmos. Cada experimento se mide contra la línea
base (código original, `-O3`) con `bench/run_bench.py` (mediana de varias
repeticiones, ejecuciones intercaladas, metadatos de máquina). Se documenta lo
que funcionó, lo que no y cuánto, en `docs/` y en el informe final.
