# Álgebra Computacional sobre Cuerpos Finitos

Implementación de algoritmos de álgebra computacional sobre cuerpos finitos, desarrollada en dos versiones: Python original y port a C++23 con aritmética de precisión arbitraria.

## Contenido

### Python (Prácticas originales)
Tres prácticas que cubren progresivamente la teoría y algoritmia de cuerpos finitos:

- `cuerpos_finitos.py` — Práctica 1: aritmética en F_p, F_p[x] y F_q = F_p[x]/(g); test de irreducibilidad de Rabin; factorización en cuerpos finitos
- `factorizacion.py` — Práctica 2: factorización de Cantor–Zassenhaus (square-free, distinct-degree, equal-degree)
- `algoritmos_rapidos.py` — Práctica 3: multiplicación de Karatsuba, matrices de Toeplitz, FFT/IFFT de Cooley–Tukey

### C++23 (Port con BigInt + templates)
Re-implementación completa del código Python en C++23, con enteros de precisión arbitraria (GMP) y diseño genérico mediante conceptos y plantillas.

- Un único conjunto de algoritmos genéricos que funciona sobre cualquier campo que modele el concepto `Field`
- `Fp` (cuerpo primo) y `Fq` (extensión) comparten exactamente el mismo código de polinomios
- `Poly<Field>` incluye: aritmética, MCD, inverso modular, potencia modular, test de irreducibilidad, factorización completa y FFT
- La multiplicación de polinomios selecciona automáticamente el algoritmo óptimo (schoolbook, Karatsuba o FFT) según el tamaño de los operandos
- Binding Python mediante `pybind11`

## Tecnologías

- Python 3
- C++23 (GMP/mpz_class, pybind11, CMake)

## Mejoras de la rama `mejoras-cpp-pybind`

Trabajo sobre el port C++ y el binding pybind11. Las prácticas originales de `python/` no se han modificado.

| Área | Qué hay | Documento |
|---|---|---|
| Núcleo C++ | Validaciones, errores tipados, `factor()` canónico, `evaluate`/`compose`/`is_squarefree`, concepto `FiniteField` | [docs/diseno.md](docs/diseno.md) |
| Paquete Python `alcp` | `GF(q)`, `K["x"]`, conversión de tipos generosa (enteros, cadenas, listas, dicts, `Fraction`), impresión Unicode/ASCII/LaTeX | [docs/diseno.md](docs/diseno.md) §6 |
| Tests | Tests C++ (más de 170 000 comprobaciones), paridad con las prácticas originales y 709 tests de pytest con cobertura del 100 % | `build.py all` |
| Ejemplos | Recorrido didáctico por la librería | [python-bindings/ejemplos.py](python-bindings/ejemplos.py) |
| Velocidad | Mismos algoritmos, 24,4x más rápido de media en la suite C++ (128 escenarios) y 9,7x en el módulo Python | [docs/optimizacion/README.md](docs/optimizacion/README.md) |
| Comparativa | Frente a python-flint, galois y SymPy | [docs/optimizacion/comparativa.md](docs/optimizacion/comparativa.md) |
| Fase 4 (nivel A) | División de Newton, inverso en F_q por Euclides y test de irreducibilidad sin la comprobación redundante: 27,7x frente al original y distancia media a FLINT de 7,2x a 5,9x | [docs/optimizacion/fase4.md](docs/optimizacion/fase4.md) |
| Fase 4 (nivel B) | Producto de Kronecker sobre Z con GMP, half-GCD y Frobenius por composición modular: 36,8x frente al original y distancia media a FLINT de 5,9x a 2,7x | [docs/optimizacion/fase4.md](docs/optimizacion/fase4.md) |
| Fase 5 | NTT con mariposas perezosas (3,5x por mariposa), Euclides en sitio en Fp64 (2,4-3,6x) y DDF de Kaltofen–Shoup: suite 1,14x sobre el nivel B y distancia media a FLINT de 2,7x a 1,66x | [docs/optimizacion/fase5.md](docs/optimizacion/fase5.md) |
| Fase 6 | Newton con producto medio y transformadas precalculadas, F_p con p pequeño (acumulador de 64 bits, Kronecker sin divisiones), DDF sin `pow_mod` repetido y conversiones directas en el binding: suite 1,17x sobre la fase 5 (con batería); binding 4-11x en conversiones | [docs/optimizacion/fase6.md](docs/optimizacion/fase6.md) |
| Fase 7 | Half-GCD V3e, F_p con p de 64-128 bits sin GMP (`FpN<2>`), NTT radix 4 con AVX2 y extensiones pequeñas con tablas o empaquetadas: suite 1,72x sobre la fase 6; F_(2^8) 45x, F_(2^127−1) 3,6x | [docs/optimizacion/fase7.md](docs/optimizacion/fase7.md) |

### Uso rápido

Se necesitan g++ con C++23 y GMP (en Windows, MSYS2 ucrt64), Python 3 con `pybind11` y, para los tests, `pytest` y `pytest-cov`.

```bash
python build.py all
```

```bash
python python-bindings/ejemplos.py
```

```python
import alcp
K = alcp.GF(9, modulus="x^2+1")
R = K["x"]
f = R("x^4 + a*x + 2")
print(f.factor())
```
