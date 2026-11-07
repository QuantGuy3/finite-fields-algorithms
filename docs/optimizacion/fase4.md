# Fase 4 (nivel A): algoritmos rápidos que ya estaban en la librería

Rama `fase4-nivel-a`, noviembre de 2026. Esta fase **sí cambia algoritmos**, pero solo
aprovecha piezas que ya existían en la librería o en las prácticas originales. No
añade algoritmos nuevos como half-GCD, Kronecker o composición modular, que serían el
nivel B. El algoritmo original de cada operación se conserva y sigue siendo accesible.

## 1. Qué cambia

| Cambio | Antes | Ahora | Algoritmo original disponible como |
|---|---|---|---|
| División con resto | Escolar, O(n·d) | Por **inversión de Newton** (`rev(q) = rev(f)·rev(g)⁻¹ mod x^(d+1)`), O(M(d) + M(n,d)), cuando el divisor y el cociente son grandes | `divmod_school` |
| `pow_mod`, `is_irreducible`, `distinct_degree` | Cada reducción mod m con la división escolar | Clase `Reducer`: calcula rev(m)⁻¹ **una sola vez** y reduce cada producto con dos multiplicaciones | `divmod_school` en cada paso |
| Inverso en F_q | Fermat, a^(q−2): unas 2·log₂ q multiplicaciones de F_q | **Euclides extendido** en F_p[t] sobre (g, a): O(e²) | `FqT::inv_fermat` |
| Test de irreducibilidad | Ben-Or y, además, x^(q^n) ≡ x (mod f) calculado desde cero | Ben-Or sin esa comprobación final, que es redundante | `ben_or_original` en `tests/test_fase4.cpp` |

La inversión de series es la iteración de Newton de `fp_toep_inf_inv` (Práctica 3,
`algoritmos_rapidos.py`), con la misma recursión: k = ⌈n/2⌉ y x_n = x_k·(2 − a·x_k) mod x^n.
La práctica la tenía implementada, pero su `fp_x_divmod` no la usaba: resolvía el sistema
triangular por sustitución, en O(n²).

**Generador aleatorio.** Los productos nuevos (`mul_no_rng`) nunca calculan la raíz de la
FFT, que consume números aleatorios. Así, la secuencia aleatoria es la misma que antes y
`Fq(p, n)` con módulo aleatorio, `factor()` y el resto dan exactamente los mismos
resultados. `test_fase4` lo comprueba, y los 128 checksums de la suite lo confirman.

## 2. Cuándo compensa Newton: calibración

Con la reducción perezosa de la fase 3, la división escolar ya es muy rápida (en Fp64,
unos 0,8 ns por operación de cuerpo). Newton solo gana con divisores y cocientes grandes
y productos rápidos. El umbral (`prefer_newton` en `polynomial.hpp`) depende de **si hay
FFT** y de **lo que cuesta cada producto de elementos**. Medidas de
`bench/calibra_newton.cpp` (resultados en `bench/results/fase4/calibra_newton.txt`),
dividiendo un polinomio de grado 2n − 1 entre uno de grado n:

| Cuerpo | Escolar / Newton aislada | Escolar / Newton con inverso amortizado | Umbral elegido (aislada / amortizado) |
|---|---|---|---|
| Fp64, F_998244353 (NTT) | 0,80x a 1024; 1,48x a 2048 | 0,93x a 256; **1,41x a 512**; 2,32x a 1024 | 1536 / 512 |
| Fp64, F_101 y F_(2^61−1) (sin NTT) | < 1 hasta 2048 | 1,0x a 1024; 2,0x a 2048 | nunca / 1536 |
| GMP, p de 255 bits (NTT multi-limb) | 1,05x a 512; 2,0x a 1024 | 1,51x a 256; 10,9x a 2048 | 768 / 256 |
| GMP, 2^127 − 1 (sin NTT) | 1,24x a 2048 | 0,96x a 512; 1,88x a 1024 | 1536 / 768 |
| F_256 (e = 8) | 1,43x a 1024 | 1,03x a 256; 1,49x a 384 | 750 / 437 (6000/e, 3500/e) |
| F_(3^30) (e = 30) | 1,25x a 192 | 1,09x a 96; 2,17x a 128 | 200 / 116 |
| F_9 (e = 2) | < 1 | 1,37x a 2048 | 3000 / 1750 |
| F_((2^61−1)^2) (FFT genérica) | < 1 | 0,92x a 2048 | umbral duplicado: no se usa en la práctica |

### Lección: la primera calibración estaba sesgada

La primera calibración repetía siempre las mismas entradas y daba Newton ganador desde
grado 128 en F_998244353. En la comparativa, sin embargo, `factor()` de grado 300 salió un
29 % más lento que en la fase 3. La causa estaba en la medida. Con una entrada repetida,
el **predictor de saltos** del procesador aprende las restas condicionales de la NTT, y
la FFT parece mucho más rápida. Por ejemplo, el producto `a*a` de grado 300 tarda 23 µs
repitiendo la misma entrada y 77 µs rotando entre 64 distintas. La calibración buena
rota entre 16 entradas. Con ella, el umbral amortizado de Fp64 con NTT subió de 128 a 512
y la regresión desapareció. La medida sesgada se conserva en
`bench/results/fase4/calibra_newton_entrada_fija.txt`.

## 3. Resultados

### Suite C++ (128 escenarios, 4 rondas intercaladas, enchufado)

| | Media geométrica |
|---|---:|
| Fase 4 frente a fase 3 | **1,19x** |
| Fase 3 frente al original | 23,3x |
| Fase 4 frente al original | **27,7x** |

Los 128 checksums coinciden en las tres versiones. Por operación, frente a la fase 3:

| Operación | Speedup | Motivo |
|---|---:|---|
| `is_irreducible` (7) | **2,69x** (2,2x – 3,4x) | Sin la comprobación final redundante |
| `elem_inv` (10) | **2,65x**; en F_q de 5,3x a 41,8x | Euclides extendido. F_p no cambia (ya usaba GMP o Euclides en 64 bits) |
| `build` de `Fq` aleatorio (2) | **1,99x** | Busca el módulo con `is_irreducible` |
| `factor` (10) | **1,33x**; en F_((2^61−1)^2) 6,9x | Inverso en F_q (`monic`, divisiones) e irreducibilidad |
| Resto | entre 0,97x y 1,10x | Sin cambios; dentro del ruido |

La división de Newton apenas aparece en la suite porque sus grados son pequeños: `pow_mod`
llega como mucho a grado 128, por debajo de los umbrales. Hay escenarios que bajan hasta un
13 % frente a la fase 3, por ejemplo `fp.ntt.factor.d40`. No siguen el camino nuevo y sus
medianas por ronda se solapan, así que es ruido de una medición que se fue calentando a lo
largo de las rondas.

### Comparativa con python-flint (16 operaciones, un hilo)

`bench/compare_libs.py`, con la variante `alcp_fase3` (módulo de la fase 3 con un hilo)
añadida. Se añadieron dos escenarios de grado alto: irreducibilidad de grado 400 y
factorización de grado 300 en F_998244353. Las 16 operaciones dan el mismo resultado en
todas las librerías. Distancia a FLINT como media geométrica: **fase 3 7,22x → fase 4
5,87x**.

| Operación | Fase 3 | Fase 4 | Mejora | FLINT | Distancia a FLINT |
|---|---:|---:|---:|---:|---:|
| Irreducibilidad grado 400, F_998244353 | 5,34 s | 1,66 s | **3,21x** | 247 ms | 21,6x → 6,7x |
| Irreducibilidad grado 200, F_101 | 222 ms | 94,4 ms | **2,35x** | 8,5 ms | 26,1x → 11,1x |
| Factorizar grado 5, F_(3^30) | 9,9 ms | 6,0 ms | **1,65x** | 2,2 ms | 4,5x → 2,7x |
| División 4000/2000, F_998244353 | 3,9 ms | 3,2 ms | 1,20x | 0,73 ms | 5,3x → 4,4x |
| Factorizar grado 300, F_998244353 | 451 ms | 391 ms | 1,15x | 25 ms | 18,0x → 15,6x |
| Resto (11 operaciones) | | | entre 0,94x y 1,19x | | sin cambio claro |

Las variantes se midieron en dos tandas del mismo día (ver `tandas` en
`bench/results/comparativa/comparativa.json`). Entre tandas, operaciones que no cambian,
como el MCD, varían hasta un ±20 %. Solo las tres primeras filas son mejoras claras.

## 4. Verificación

- `build.py all`: 7 tests C++. `tests/test_fase4.cpp` aporta 610 comprobaciones nuevas:
  - Newton frente a la división escolar en 7 tipos de cuerpo, con y sin tablas FFT;
  - `pow_mod` frente a la potencia con división escolar;
  - el inverso por Euclides frente al de Fermat en **todos** los elementos de F_9, F_256 y
    F_(7^4), y en muestras de F_(3^30) y F_((2^61−1)^2);
  - la irreducibilidad frente al test original y frente al número de irreducibles
    mónicos de Gauss, (1/d)·Σ μ(e)·q^(d/e), para F_2 hasta grado 8, F_3 hasta 5 y F_9
    hasta 3;
  - el generador aleatorio queda intacto.

  Además pasan el oráculo de paridad con las prácticas, `xcheck` (924) y 709 tests de
  pytest con un 100 % de cobertura del paquete.
- Suite C++: 128/128 checksums iguales al original y a la fase 3.
- `python/` intacto.

## 5. Limitaciones

- Una sola máquina (Zen 3, Windows): los umbrales pueden ser distintos en otra CPU.
- La comparativa de fase 3 frente a fase 4 mezcla dos tandas. Las mejoras de menos de un
  ~20 % no son concluyentes.
- Los umbrales de F_q se ajustan con una fórmula en función de e, medida en tres
  extensiones (F_9, F_256, F_(3^30)) más una con FFT.
- Smart App Control bloqueó varias veces python-flint, galois y ejecutables recién
  compilados. Se volvieron a medir y no afecta a los resultados, pero alarga las medidas.

## 6. Conclusiones y recomendaciones

- El nivel A reduce la distancia media con FLINT de 7,2x a 5,9x con un hilo, y de forma
  muy marcada en la irreducibilidad (de 21,6x a 6,7x). La división de Newton aporta poco
  en los tamaños habituales, porque la división escolar con reducción perezosa ya es muy
  eficiente. Solo compensa a partir de unos cientos de coeficientes, y antes con
  elementos caros (GMP, extensiones grandes).
- **Propuesta (fuera del nivel A): hacer sin saltos las restas condicionales de la NTT en
  palabras.** La diferencia de 3x entre entradas repetidas y variadas indica que la NTT
  sufre fallos de predicción de saltos con datos reales. Es una hipótesis: no se ha
  comprobado en el ensamblador si el compilador genera saltos o `cmov`. Si se confirma,
  escribir esas restas con máscaras aceleraría todos los productos FFT, y con ellos
  `pow_mod`, la factorización y Newton, sin cambiar ningún algoritmo.
- Lo que queda frente a FLINT es sobre todo algorítmico: half-GCD, Kronecker en F_q,
  composición modular para Frobenius y NTT con varios primos. Es el nivel B.

---

# Fase 4, nivel B: algoritmos nuevos

Rama `fase4-nivel-b`, construida sobre el nivel A. Añade tres algoritmos que no estaban en
la librería. Cada uno da **exactamente el mismo resultado** que el original, que se conserva:

| Algoritmo | Sustituye a | Original disponible como |
|---|---|---|
| **Producto por sustitución de Kronecker sobre Z** (`include/alcp/kronecker.hpp`, `MulAlgo::Kronecker`, `"kronecker"` en Python): los coeficientes se empaquetan en huecos de b bits de un entero, se multiplica con `mpn_mul` de GMP (Toom-Cook y Schönhage–Strassen) y se desempaqueta reduciendo módulo p. En F_q es Kronecker en dos niveles: cada elemento ocupa un bloque de 2e − 1 huecos y cada bloque se reduce módulo g | Escolar, Karatsuba o NTT, según el cuerpo | `MulAlgo::Schoolbook`, `Karatsuba`, `FFT` |
| **Half-GCD** (Thull–Yap) en `gcd` y `gcd_ext`: la recursión sobre la mitad alta de los coeficientes devuelve la matriz 2×2 de cocientes. Los cocientes son los del Euclides clásico, así que el mcd **y los coeficientes de Bézout** son idénticos | Euclides clásico, O(n²) | `gcd_classic`, `gcd_ext_classic` |
| **Frobenius por composición modular** (Brent–Kung): H_i = H_{i−1}(H_1) mod f con las potencias de H_1 precalculadas, en `is_irreducible` y `distinct_degree` | H_i = H_{i−1}^q con `pow_mod` | `FrobeniusAlgo::Potencia` |

Se descartó la **NTT con varios primos y CRT** que se había propuesto. Kronecker con GMP
cubre el mismo objetivo, productos rápidos para cualquier p, también p ≥ 2^62 y
extensiones, con un único mecanismo y sin raíces de la unidad.

**Generador aleatorio.** Ni Kronecker ni las matrices de half-GCD ni la composición
calculan la raíz de la FFT. Cuando el despacho elige Kronecker, la raíz se sigue pidiendo
exactamente cuando la habría pedido el original. El camino rápido de `gcd_ext` solo se usa
si esa raíz ya existe, porque el `gcd_ext` original la podía pedir a mitad del algoritmo.
`test_fp64` y `test_nivel_b` lo comprueban.

## Calibración (entradas variadas, Zen 3)

Los ficheros completos están en `bench/results/fase4/`:

- **Kronecker** (`calibra_kronecker_hilos.txt` y `calibra_kronecker_1hilo.txt`, frente al
  despacho anterior). En F_q y en F_p con GMP gana desde 32 coeficientes, de 1,5x a 25x y
  hasta 59x en F_9 con 4096 coeficientes. En Fp64 sin NTT gana desde 32 si p < 2^32 (hasta
  20x con p = 101) y desde 256 si p es mayor. En Fp64 con NTT solo gana entre 128 y 1024
  coeficientes (1,5x a 2,5x); fuera de ese rango, la NTT en palabras con hilos empata o gana.
  **Nota:** en esas tablas, las filas de p = 2^61 − 1 se midieron antes de corregir un fallo
  que truncaba p, porque `mpz_get_ui` es de 32 bits en Windows. El trabajo medido es el
  mismo, pero esos productos daban un resultado erróneo, que los tests detectaron.
- **Half-GCD** (`calibra_hgcd.txt`). `gcd_ext` gana desde grado 128-512 según el cuerpo, con
  hasta 8,7x a grado 4096. `gcd` gana desde 384 (F_p con GMP, F_q y p < 2^16), desde 1024
  (Fp64 con p < 2^32) y desde 2048 (p mayor). Acumular cada paso clásico con 2 productos en
  lugar de un producto de matrices completo adelantó estos puntos de corte.
- **Frobenius** (`calibra_frobenius.txt`, con entradas irreducibles). La composición gana
  siempre, de 1,2x con q = 2 a 18x con q = 2^127 − 1. `Auto` aplica la regla de alquilar o
  comprar: cambia a composición cuando lo gastado en potencias iguala el coste de
  precalcular, de modo que los reducibles, que salen pronto del bucle, no pagan nada.
- **Newton, otra vez** (`calibra_newton.txt`; la del nivel A está en
  `calibra_newton_nivel_a.txt`). Con productos de Kronecker los umbrales bajan mucho. Con
  el inverso amortizado: F_101 desde 64, F_(2^127−1) desde 96, F_(3^30) desde 32.

## Resultados

### Suite C++ (128 escenarios, 4 rondas intercaladas, enchufado)

| | Media geométrica |
|---|---:|
| Nivel B frente al nivel A | **1,25x** |
| Nivel A frente al original | 29,4x |
| Nivel B frente al original | **36,8x** |

Los 128 checksums son iguales en las tres versiones (`bench/results/fase4b/`). Por
operación, frente al nivel A:

| Operación | Speedup |
|---|---:|
| `mul_auto` | 3,73x |
| `is_irreducible` | 3,41x |
| `toeplitz_lower_inverse` | 2,45x |
| `divmod` | 1,87x |
| `pow_mod` | 1,59x |
| `gcd_ext` | 1,31x |
| `factor` | 1,04x |
| `gcd` | 0,93x |

Por cuerpo: F_(2^127−1) 1,54x, F_101 1,49x, F_7 1,47x, `fast` 1,29x, F_((2^61−1)^2) 1,23x.

**Anomalía sin resolver.** En la suite, el `gcd` de grado 500 en Fp64 sale a 0,78x en
F_998244353 y a 0,80x en F_(2^61−1), de forma constante en las 4 rondas. A ese grado
half-GCD no se activa y se ejecuta el mismo Euclides clásico. Aislado, con las mismas
entradas (factor común de grado n/8) y las cabeceras de los dos niveles, tarda lo mismo:
1211 µs con el nivel A y 1200 µs con el B. Puede ser un efecto del binario completo de la
suite, como la disposición del código, pero no está confirmado. También bajan
`fq.build.q3e30` (0,86x) y `fast.ntt.toeplitz_lower_vec` (0,86x), con medianas por ronda
que se solapan.

### Comparativa con python-flint (16 operaciones, un hilo, una misma tanda)

Distancia media a FLINT (media geométrica): **fase 3 6,9x → nivel A 5,9x → nivel B 2,7x**.
Las 16 operaciones dan el mismo resultado en todas las librerías.

| Operación | Fase 3 | Nivel A | Nivel B | FLINT | Distancia a FLINT |
|---|---:|---:|---:|---:|---:|
| Producto grado 500, F_(2^127−1) | 3,0 ms | 3,1 ms | **0,30 ms** | 0,33 ms | 9,2x → 9,5x → **0,92x** |
| Producto grado 1000, F_256 | 15,0 ms | 15,1 ms | **1,0 ms** | 0,76 ms | 19,7x → 19,8x → 1,3x |
| Producto grado 200, F_(3^30) | 11,9 ms | 11,9 ms | **0,95 ms** | 0,52 ms | 22,8x → 22,8x → 1,8x |
| Irreducibilidad grado 400, F_998244353 | 5,06 s | 1,52 s | **269 ms** | 216 ms | 23,4x → 7,0x → 1,2x |
| Irreducibilidad grado 200, F_101 | 201 ms | 77,5 ms | **27,5 ms** | 7,6 ms | 26,3x → 10,2x → 3,6x |
| Factorizar grado 300, F_998244353 | 344 ms | 345 ms | **133 ms** | 23,5 ms | 14,6x → 14,7x → 5,7x |
| Factorizar grado 100, F_101 | 6,6 ms | 6,7 ms | 4,4 ms | 0,45 ms | 14,6x → 14,9x → 9,7x |
| Factorizar grado 20, F_(2^127−1) | 29,6 ms | 29,5 ms | 17,5 ms | 2,6 ms | 11,5x → 11,5x → 6,8x |
| Factorizar grado 30, F_256 | 13,6 ms | 12,6 ms | 8,0 ms | 8,3 ms | 1,6x → 1,5x → 1,0x |
| División 4000/2000, F_998244353 | 3,2 ms | 3,3 ms | 3,4 ms | 0,71 ms | 4,5x → 4,6x → 4,7x |
| MCD grado 1000, F_998244353 | 4,5 ms | 4,7 ms | 4,7 ms | 0,59 ms | 7,6x → 7,9x → 8,0x |

Los dos últimos no cambian: a esos tamaños, en F_998244353, Newton y half-GCD quedan por
debajo de sus umbrales.

## Verificación del nivel B

- `build.py all`: 8 tests C++, oráculo, módulo y 710 tests de pytest con un 100 % de
  cobertura del paquete. `tests/test_nivel_b.cpp` comprueba:
  - Kronecker frente al escolar en 11 cuerpos, incluidos el peor caso de huecos con
    coeficientes p − 1, p = 2^62 − 57 y p = 2^521 − 1;
  - half-GCD frente al clásico (mcd y Bézout idénticos) en 6 cuerpos, con F_2, mcd grandes,
    ceros y grados iguales;
  - el generador aleatorio;
  - la composición frente a `compose` seguido de `%`;
  - Frobenius e `is_irreducible` / `distinct_degree` frente al algoritmo original.
- Suite C++: 128/128 checksums iguales al original.
- `python/` intacto.
- No se ha vuelto a medir la cobertura de líneas C++ (`build.py cov-cpp`).

## Conclusiones del nivel B

- Kronecker con GMP es la pieza de más impacto: lleva los productos de F_q y de p grande al
  nivel de FLINT, e incluso por delante con p = 2^127 − 1. Además arrastra a Newton,
  half-GCD y `pow_mod`.
- La irreducibilidad queda a 1,2x de FLINT en F_998244353. La factorización sigue a 3-10x:
  la fase por grados iguales (Cantor–Zassenhaus) y el número de MCD siguen siendo los del
  algoritmo original.
- **Siguientes pasos posibles:**
  - restas sin saltos en la NTT en palabras (el sesgo de 3x del predictor de saltos);
  - DDF de Kaltofen–Shoup (paso pequeño y paso gigante) para reducir el número de MCD;
  - una calibración de los umbrales en otra máquina.
