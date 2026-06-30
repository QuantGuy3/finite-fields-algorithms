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
