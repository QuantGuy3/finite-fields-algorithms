# Investigación: cómo acelerar más la implementación (propuesta de fase 5)

> **Estado:** las propuestas 1 a 3 están implementadas en la rama `fase5`; resultados en [fase5.md](fase5.md).

Fecha: 7 de noviembre de 2026. Punto de partida: rama `mejoras-cpp-pybind` en `7880a88`, con
las fases 3 y 4 (niveles A y B) ya fusionadas.

Este documento separa lo **medido** (experimentos en `bench/experimentos/`, salidas en
`bench/results/investigacion/mediciones.txt`) de lo **estimado** (cálculos de coste o datos
de la bibliografía, sin implementar aún). Las estimaciones son orientativas y hay que
confirmarlas con la suite antes de darlas por buenas.

## Resumen

Hay dos hallazgos medidos y uno algorítmico, por este orden de prioridad:

1. **La NTT en palabras es entre 3 y 10 veces más lenta de lo que debería.** Cada mariposa
   tarda unos **14 ns**, cuando una mariposa de Shoup escalar suele costar 1–2 ns. La causa
   está en las tres comparaciones por mariposa: con datos aleatorios el predictor de saltos
   falla la mitad de las veces. Una versión sin saltos baja a 5,2 ns (2,7x) y la mariposa
   perezosa de Harvey (2014) a **4,0 ns (3,5x)**, con resultados idénticos bit a bit. Esto
   también explica por qué Kronecker con GMP le ganaba a la NTT entre 128 y 1024
   coeficientes.
2. **El Euclides clásico gasta unos 11 ns por operación de coeficiente.** Al reservar y
   reducir acumuladores de 128 bits en cada `a % b`, cada paso recorre el vector varias
   veces. Un Euclides en sitio sobre `u64`, con multiplicador de Shoup y sin saltos, es
   **2,4x a 3,6x** más rápido y obtiene el mismo mcd. Ese bucle es la base de `gcd`, del caso
   base de half-GCD y de la DDF.
3. **En `factor`, la DDF consume el 96–100 % del tiempo.** Calcula un paso de Frobenius y un
   MCD por grado, unos n/2 de cada. La DDF de **Kaltofen–Shoup** (paso pequeño y paso
   gigante), que es la que usan FLINT y NTL, reduce eso a unas 2·√(n/2) composiciones y
   √(n/2) MCD. Para grado 300 en F_998244353 se estima unas **4x** solo por el algoritmo.

Las tres propuestas dan **exactamente el mismo resultado** que la versión actual. Las dos
primeras son de implementación: no cambian ni el algoritmo ni el orden de las operaciones
aritméticas.

## 1. Dónde se va el tiempo (medido)

Medido con `bench/experimentos/perfil_factor.cpp` y `perfil_mul.cpp`, en Zen 3 con g++ 15.2,
`-O3 -march=native`, un hilo y entradas variadas.

**Factorización** (ms por llamada, polinomios mónicos aleatorios):

| Caso | Total | Radical | DDF: Frobenius | DDF: MCD | EDF | Multiplicidad |
|---|---:|---:|---:|---:|---:|---:|
| F_998244353, grado 300 | 380 | 1,4 | **278 (73 %)** | **99 (26 %)**, 103 MCD | 0,2 | 0,3 |
| F_101, grado 100 | 10,0 | 0,2 | 4,7 (47 %) | 4,9 (49 %), 35 MCD | 0,1 | 0,1 |
| F_101, grado 300 | 141 | 1,3 | 46 (33 %) | **93 (66 %)**, 104 MCD | 0,03 | 0,2 |

La fase por grados iguales (Cantor–Zassenhaus) es despreciable con entradas aleatorias,
porque tienen pocos factores. Todo el coste está en la DDF.

**Primitivas en F_998244353** (µs; la NTT existe para estos tamaños):

| n | `mul` Auto | NTT | Kronecker | `divmod` 2n/n | `gcd` n |
|---:|---:|---:|---:|---:|---:|
| 300 | 59 | 188 | 61 | 180 | 1 261 |
| 1000 | 303 | 485 | 284 | 1 568 | 11 643 |
| 2000 | 1 039 | 1 025 | 886 | 4 321 | 32 944 |

Una NTT de 2^11 puntos tarda 160 µs, es decir, 14 ns por cada una de sus 11 264
mariposas. El producto de grado 1000 hace tres transformadas, que suman unos 480 µs de los
485 µs totales. La NTT es, por tanto, todo el coste del producto.

## 2. Propuestas priorizadas

| # | Propuesta | Tipo | Evidencia | Ganancia esperada | Esfuerzo |
|---|---|---|---|---|---|
| 1 | Mariposas sin saltos y perezosas (Harvey) en la NTT en palabras | Implementación | **Medido**: 3,5x en la transformada | Producto NTT 2,5–3,5x; arrastra a divmod, Newton, pow_mod, Frobenius y half-GCD | Bajo |
| 2 | Euclides y división escolar en sitio con multiplicador de Shoup | Implementación | **Medido**: 2,4–3,6x en `gcd_classic` | MCD de la DDF 2,5–3,5x; adelanta los umbrales de half-GCD | Bajo |
| 3 | DDF de Kaltofen–Shoup (paso pequeño y paso gigante) | Algoritmo | Estimado (§3.3) | `factor` unas 4x a grado 300; crece con n | Medio |
| 4 | Más NTT: sin permutación de bits, radix 4 y AVX2 para p < 2^32 | Implementación | Estimado y bibliografía | Otras 1,5–3x en la NTT | Medio-alto |
| 5 | Producto medio y transformadas precalculadas del módulo en `Reducer` | Algoritmo | Estimado y bibliografía | pow_mod y Frobenius 1,3–1,6x | Medio |
| 6 | FFT truncada (van der Hoeven) | Algoritmo | Bibliografía | Hasta 2x justo después de cada potencia de 2; 1,2–1,4x de media | Medio |
| 7 | Kronecker en varios puntos (KS2/KS4 de Harvey, 2009) | Algoritmo | Bibliografía | Productos en F_q y p grande 1,3–1,8x | Medio |
| 8 | Cuerpo F_2[x] empaquetado en bits con `PCLMULQDQ` | Nuevo backend | Bibliografía (gf2x) | 10x a 100x en F_2 y F_2^k | Alto |
| 9 | Tablas de logaritmos de Zech para F_q pequeño (q ≤ 2^16) | Nuevo backend | Estimado | Producto de elementos de F_q de ~e² operaciones a 3 accesos a tabla | Medio |
| 10 | Grupo de hilos persistente en lugar de `std::async` | Implementación | Medido indirectamente | Pequeña; solo con n ≥ 8192 | Bajo |

Recomendación: las propuestas 1 y 2 primero, porque son baratas, están medidas y aceleran
todo lo demás. Después la 3, que es la que más acerca `factor` a FLINT. Tras cada una hay
que recalibrar los umbrales (Kronecker frente a NTT, Newton, half-GCD), porque cambian los
costes relativos.

## 3. Detalle

### 3.1 NTT sin saltos y perezosa (propuesta 1)

La mariposa actual (`include/alcp/fft.hpp`, `word::transform`) hace tres comparaciones con
datos impredecibles. Una está en `mul_shoup` y dos en la suma y la resta, todas con salto
condicional. En `bench/experimentos/ntt_sin_saltos.cpp` se compara con dos variantes:

- **Sin saltos.** Cada corrección se escribe con máscaras (`s += p & -(s >> 63)`), y el
  compilador ya no puede emitir saltos.
- **Perezosa** (Harvey, 2014). Los valores viven en [0, 4p), así que cada mariposa hace una
  sola corrección condicional, de la entrada u a [0, 2p), y el producto de Shoup se queda
  en [0, 2p) sin corregir. Al final se normaliza una vez. Hace falta 4p < 2^64, que se
  cumple porque `kMaxBits = 62`.

| Tamaño | Actual | Sin saltos | Perezosa |
|---|---:|---:|---:|
| 2^11 | 14,0 ns | 5,2 ns | **4,0 ns** |
| 2^13 | 14,2 ns | 5,5 ns | **4,4 ns** |
| 2^16 | 14,9 ns | 6,1 ns | **5,3 ns** |

Las tres dan el mismo vector. La aritmética es la misma módulo p y solo cambia el
representante intermedio, que se normaliza al final. El producto punto a punto
(`mul_mod`, con una división de 128 bits) y la transformada inversa admiten el mismo trato.

Esto también explica la «calibración sesgada» de la fase 4. Con entradas repetidas, el
predictor aprendía los saltos y la NTT parecía hasta 3 veces más rápida.

### 3.2 Euclides en sitio (propuesta 2)

`gcd_classic` llama a `a % b`. En Fp64, `divmod_school` usa la reducción perezosa: reserva
un vector de acumuladores de 128 bits, lo inicializa, acumula y reduce todos los
coeficientes. En la sucesión de restos normal cada cociente tiene grado 1, así que esa
preparación cuesta más que el trabajo útil. El experimento
(`bench/experimentos/euclides_en_sitio.cpp`) hace los **mismos cocientes en el mismo orden**
sobre dos vectores `u64` que se intercambian. Cada fila es un `axpy` con el multiplicador de
Shoup de −cociente, precalculado una vez por fila, y sin saltos:

| Cuerpo | n = 300 | n = 1000 |
|---|---:|---:|
| F_998244353 | 1 318 → 367 µs (**3,6x**) | 11 632 → 3 455 µs (**3,4x**) |
| F_101 | 1 184 → 342 µs (3,5x) | 10 754 → 3 318 µs (3,2x) |
| F_(2^61−1) | 1 044 → 426 µs (2,4x) | 9 791 → 3 615 µs (2,7x) |

El mcd es idéntico en todos los casos. Lo mismo se aplica a `gcd_ext_classic`, que además
actualiza los cofactores en sitio, y al caso base de `hgcd`. Ese `axpy` se puede vectorizar
con AVX2 cuando p < 2^32 (propuesta 4).

### 3.3 DDF de Kaltofen–Shoup (propuesta 3)

La DDF actual calcula, para d = 1, 2, …, H_d = x^(q^d) mod f y gcd(H_d − x, f). Son unos
n/2 pasos de Frobenius y n/2 MCD. La variante de paso pequeño y paso gigante (Kaltofen y
Shoup, 1998; von zur Gathen y Shoup, 1992) toma l ≈ √(n/2) y hace lo siguiente:

1. Pasos pequeños: H_1, …, H_l, con l composiciones modulares.
2. Pasos gigantes: G_j = H_(lj) = G_(j−1)(H_l), con unas n/(2l) composiciones.
3. Para cada j, I_j = ∏_(i<l) (G_j − H_i) mod f. Son l productos modulares, que no son
   composiciones y salen baratos, y un **único** MCD gcd(I_j, f) por intervalo.
4. Solo en los intervalos con factor se separan los grados exactos, con MCD sobre el trozo
   pequeño.

La salida (pares grado–producto) es única, así que el resultado es idéntico. Además la DDF
no usa el generador aleatorio.

Coste estimado para grado 300 en F_998244353, con los tiempos de §1 (1,85 ms por
composición, unos 0,13 ms por producto modular y 0,96 ms por MCD):

| | Actual | Kaltofen–Shoup (l ≈ 13) |
|---|---:|---:|
| Composiciones | ~150 → 278 ms | ~25 → 46 ms |
| Productos modulares | — | ~150 → 20 ms |
| MCD | 103 → 99 ms | ~12 (+ separación) → 15 ms |
| **Total DDF** | **378 ms** | **~80 ms (≈ 4,7x)** |

Con las propuestas 1 y 2 aplicadas, las composiciones, los productos y los MCD se abaratan
a la vez. Comparada con FLINT (23,5 ms en este caso), la distancia pasaría de 5,7x a una
cifra estimada entre 1x y 2x, pendiente de medir.

### 3.4 Otras propuestas

- **Más NTT (4).**
  - Quitar la permutación de inversión de bits: la directa se hace en DIF y la inversa en
    DIT, y el producto punto a punto no depende del orden.
  - Fusionar las etapas pequeñas, con `half` igual a 1 o 2, que hoy son bucles casi vacíos.
  - Radix 4.
  - Con p < 2^31, vectorizar con AVX2: 8 carriles de 32 bits con Montgomery o Shoup. Es lo
    que hacen NTL y la `fft_small` de FLINT 3, que trabaja en coma flotante con FMA y
    primos de 50 bits.
  - La ganancia es incierta y depende de la máquina. Requiere detectar la CPU en tiempo de
    ejecución o compilar con `-march`.
- **Producto medio y transformadas del módulo (5).** En la iteración de Newton,
  x(2 − a·x) solo necesita los coeficientes centrales del producto. El producto medio
  (Hanrot, Quercia y Zimmermann, 2004) lo calcula con una transformada del tamaño de un
  producto normal y ahorra cerca de 1/3 de la inversión. En `Reducer`, además, las
  transformadas de rev(m)⁻¹ y de m se pueden guardar una vez, igual que el `FFTRep` de NTL,
  con lo que cada reducción hace 2 directas y 2 inversas en lugar de 3 y 3.
- **FFT truncada (6).** Hoy un producto de longitud 2049 usa una NTT de 4096. La TFT
  (van der Hoeven, 2004) hace que el coste crezca de forma suave con la longitud.
  **Restricción:** `tables(K, k)` tiene que seguir llamándose en el mismo punto y con el
  mismo k para no alterar el generador aleatorio.
- **Kronecker en varios puntos (7).** En lugar de evaluar en 2^b con huecos de b bits, se
  evalúa en 2^b y 2^(−b), o en ±2^b (KS2/KS4, Harvey 2009), y se reconstruye. Los enteros
  quedan de la mitad de bits. Afecta a F_q y a p grande, donde Kronecker ya es el algoritmo
  elegido.
- **F_2[x] empaquetado (8).** 64 coeficientes por palabra y producto sin acarreo con
  `PCLMULQDQ`, como gf2x (Brent, Gaudry, Thomé y Zimmermann, 2008). Es un tipo de cuerpo
  nuevo y requiere bastante trabajo, pero es la mayor ganancia posible para F_2 y F_2^k.
- **Logaritmos de Zech (9).** Para F_q con q ≤ 2^16, cada elemento se guarda como un
  logaritmo discreto: el producto es una suma y la suma usa una tabla. Afecta sobre todo a
  F_256 y F_9 en tamaños pequeños, por debajo del umbral de Kronecker.

### 3.5 Descartadas

- **Composición modular de Kedlaya–Umans.** Su complejidad asintótica es mejor, pero sus
  constantes no la hacen práctica a estos tamaños. Ninguna librería de referencia la usa.
- **Schönhage–Strassen propio.** GMP ya lo usa dentro de `mpn_mul` para los productos de
  Kronecker grandes.
- **Berlekamp para q pequeño.** Solo compensa con muchos factores de grado pequeño. Con la
  DDF de Kaltofen–Shoup y entradas aleatorias no se espera ganancia; se revisaría si el
  caso de uso lo pidiera.

## 4. Restricciones que hay que mantener

- `python/` no se toca.
- Resultados idénticos: 128 checksums de la suite, oráculo y tests de cada fase.
- **Generador aleatorio.** Ninguna propuesta puede cambiar cuándo se pide la raíz de la FFT
  (`tables()`). La propuesta 3 no usa aleatoriedad. Si se cambiara la EDF (traza por
  composición), cambiaría el consumo del generador y por tanto la secuencia posterior;
  queda fuera de este plan.
- Cada umbral se calibra con entradas variadas (rotando 8 o 16), nunca con la misma
  entrada repetida.

## 5. Limitaciones de esta investigación

- Todo se ha medido en una sola máquina (Zen 3, Windows, MinGW). En otras CPU el efecto del
  predictor de saltos y de AVX2 puede ser distinto.
- Los experimentos de §3.1 y §3.2 son prototipos fuera de la librería. Integrados puede
  haber costes de conversión o de reserva que reduzcan la ganancia.
- Las cifras de las propuestas 3 a 9 son estimaciones de coste o vienen de la bibliografía,
  no de mediciones propias.
- Las referencias se citan de memoria. Conviene comprobar los datos bibliográficos antes de
  usarlos en un trabajo formal.

## Referencias

- D. Harvey, «Faster arithmetic for number-theoretic transforms», *Journal of Symbolic
  Computation* 60 (2014), 113–119.
- D. Harvey, «Faster polynomial multiplication via multipoint Kronecker substitution»,
  *Journal of Symbolic Computation* 44 (2009), 1502–1510.
- E. Kaltofen y V. Shoup, «Subquadratic-time factoring of polynomials over finite fields»,
  *Mathematics of Computation* 67 (1998), 1179–1197.
- J. von zur Gathen y V. Shoup, «Computing Frobenius maps and factoring polynomials»,
  *Computational Complexity* 2 (1992), 187–224.
- G. Hanrot, M. Quercia y P. Zimmermann, «The middle product algorithm I», *Applicable
  Algebra in Engineering, Communication and Computing* 14 (2004), 415–438.
- J. van der Hoeven, «The truncated Fourier transform and applications», *ISSAC 2004*.
- R. P. Brent, P. Gaudry, E. Thomé y P. Zimmermann, «Faster multiplication in GF(2)[x]»,
  *ANTS VIII* (2008).
- K. S. Kedlaya y C. Umans, «Fast polynomial factorization and modular composition», *SIAM
  Journal on Computing* 40 (2011).
- J. von zur Gathen y J. Gerhard, *Modern Computer Algebra*, 3.ª ed., Cambridge University
  Press (2013), capítulos 8, 9, 11 y 14.
- V. Shoup, NTL (`ZZ_pX`, `FFTRep`); documentación de FLINT (`nmod_poly`, `fft_small`).
