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
