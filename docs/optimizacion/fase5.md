# Fase 5: NTT perezosa, Euclides en sitio y DDF de Kaltofen–Shoup

Rama `fase5`, construida sobre `mejoras-cpp-pybind` (fase 4, nivel B, `7880a88`). Aplica las
propuestas 1 a 3 de [investigacion-fase5.md](investigacion-fase5.md). Las tres dan
**exactamente el mismo resultado** que antes, y el generador aleatorio se consume igual.

| # | Cambio | Dónde | Tipo |
|---|---|---|---|
| 1 | Mariposas perezosas de Harvey, sin saltos, en la NTT en palabras (valores en [0, 4p) entre etapas y normalización al final) | `include/alcp/fft.hpp`, `word::transform` | Implementación |
| 2 | Euclides clásico **en sitio** en Fp64 sobre dos vectores `u64`, con multiplicador de Shoup por fila y sin saltos. Lo usan `gcd_classic`, `gcd_ext_classic` y el caso base de `hgcd` | `include/alcp/polynomial.hpp`, `namespace euclides`, `gcd_classic_w`, `gcd_ext_classic_w`, `hgcd_base_w` | Implementación |
| 3 | DDF de **Kaltofen–Shoup** (paso pequeño y paso gigante) en `distinct_degree` y, por tanto, en `factor`. La de la fase 4 queda como `distinct_degree_classic()` | `include/alcp/polynomial.hpp`, `ddf_ks` | Algoritmo |

## 1. NTT perezosa

La mariposa anterior hacía tres comparaciones con datos impredecibles y gastaba unos
**14 ns** por mariposa por los fallos del predictor de saltos. Ahora cada mariposa hace una
sola corrección, escrita con máscara, y el producto de Shoup se deja en [0, 2p). Como
p < 2^62, se cumple 4p < 2^64. Dentro de la librería:

| Transformada | Antes | Fase 5 |
|---|---:|---:|
| 2^11 puntos | 160 µs (14,2 ns por mariposa) | 51 µs (4,5 ns) |
| 2^13 puntos | 744 µs (14,0 ns) | 214 µs (4,0 ns) |
| 2^16 puntos | 8,5 ms (16,2 ns) | 3,0 ms (5,8 ns) |

**Recalibración del despacho** (`bench/calibra_ntt_kronecker.cpp`, en
`bench/results/fase5/calibra_ntt_kronecker_*.txt`):

- Con NTT disponible, la NTT gana ya a Kronecker desde 64 coeficientes. Por ejemplo, con
  1024 tarda 154 µs frente a 292 µs, y con 2048, 336 µs frente a 976 µs. `prefer_kronecker`
  ya no elige Kronecker si la NTT se puede usar.
- `prefer_fft` en Fp64 pasa de `small >= 48` a `small >= 64` en el caso equilibrado.
- `mul_no_rng` (productos de Newton y de half-GCD, que no pueden pedir la raíz de la FFT)
  razona como si no hubiera FFT cuando la caché está vacía, y así elige Kronecker. Una
  primera versión no lo hacía: caía a Karatsuba y `fp.ntt.divmod` salía a 0,35x en la suite.
- Newton en Fp64 con p < 2^32 y NTT: 128 con el inverso amortizado y 512 aislada (antes
  192 y 768).

## 2. Euclides en sitio

En cada paso, `a % b` reservaba e inicializaba un vector de acumuladores de 128 bits.
Con cocientes de grado 1, que es el caso normal, eso costaba más que el propio trabajo.
Ahora el paso actualiza en sitio con un `axpy` sin saltos. Los cocientes son los mismos y
salen en el mismo orden:

- Si el cociente y el divisor son grandes (`prefer_newton`), se usa la división de Newton
  de siempre.
- Los cofactores de `gcd_ext_classic` se actualizan en sitio. Si el cociente tiene más de
  64 coeficientes, se usa el producto rápido (`mul_no_rng`).
- **Generador aleatorio.** El original calculaba q·x1 con la regla Auto, que pide las
  tablas FFT, y con ellas la raíz aleatoria, la primera vez que el producto tiene 64 o más
  coeficientes. El camino en sitio repite esa llamada a `fft::tables` exactamente en el
  mismo punto. `test_fase5` lo comprueba.

Con el clásico entre 2,4x y 3,6x más rápido, half-GCD solo compensa con grados mucho
mayores. Se recalibró con `bench/calibra_hgcd.cpp -DALCP_HGCD_UMBRAL=128 -DALCP_HGCD_MIN=…`
(salidas en `bench/results/fase5/calibra_hgcd_min*.txt`):

| Fp64 | Umbral `gcd` / `gcd_ext` antes | Fase 5 |
|---|---|---|
| p < 2^16 | 384 / 128 | 2048 / 1024 |
| p < 2^32 | 1024 / 384 | 4096 / 2048 |
| p mayor | 2048 / 512 | 8192 / 4096 |

En Fp64, `kHgcdMin`, el grado del caso base de la recursión, pasa de 128 a 256.

## 3. DDF de Kaltofen–Shoup

Con l = ⌈√(n/2)⌉, el algoritmo:

1. Calcula los pasos pequeños h_i = x^(q^i) y los gigantes G_j = x^(q^(lj)).
2. Para cada intervalo de grados (l(j−1), lj], hace un único MCD de ∏(G_j − h_i) con lo
   que queda por factorizar.
3. Separa el producto encontrado por grados crecientes.

Cada paso, ya sea a → a^q (pequeño) o a → a^(q^l) (gigante), se hace elevando con
`pow_mod` o componiendo con Brent–Kung (k = ⌈√m⌉ potencias guardadas). Se elige lo que
cueste menos productos modulares, de modo que en F_2 se eleva y con q grande se compone.
Cuando lo que queda baja a la mitad del módulo de trabajo, todo se reduce módulo ese resto.

- **Resultado idéntico.** La salida (grado, producto de los irreducibles de ese grado),
  ordenada por grado, es única.
- **Generador aleatorio.** El primer paso de Frobenius es común a los dos caminos. Si
  después de él la raíz de la FFT aún no existe y el cuerpo admite una FFT de 64 o más, se
  sigue por el camino de la fase 4. Pasa, por ejemplo, con p = 193 o p = 257, y el test lo
  cubre.
- **Cuándo se usa** (`usar_ddf_ks`). Con grado ≥ 40 y bitlen(q) + popcount(q) − 1 ≥ 4,
  es decir, q ≥ 5. Con q = 2, 3 o 4, desde grado 256: con q diminuto cada paso del
  original es casi un solo cuadrado.

Calibración (`bench/calibra_ddf.cpp`, un hilo, partes libres de cuadrados de polinomios
aleatorios; salida en `bench/results/fase5/calibra_ddf.txt`). Las salidas son idénticas en
todos los casos:

| Cuerpo | Grado 48 | 64 | 96 | 150 | 300 |
|---|---:|---:|---:|---:|---:|
| F_101 | 1,60x | 1,82x | 2,10x | 2,34x | 2,84x |
| F_998244353 | 2,29x | 2,44x | 2,37x | 1,83x | **4,14x** |
| F_(2^61−1) | 2,33x | 2,17x | 2,16x | **5,87x** | — |
| F_9 | 1,65x | 1,73x | 1,99x | 2,65x | — |
| F_(2^127−1) | 2,71x | 2,63x | 2,63x | — | — |
| F_2 | 0,59x | 0,59x | 0,71x | 0,89x | 1,13x (600: 2,65x) |

## Resultados

**Advertencia de medición.** El portátil funcionó con batería durante todas las medidas
(plan equilibrado), así que los tiempos absolutos son más altos que en fases anteriores.
Por ejemplo, `fp127.factor.d20` del nivel B tarda 47,8 ms en esta tanda y 17,5 ms en la
del nivel B. Las comparaciones son válidas porque cada una se midió en la misma tanda
(suite intercalada; comparativa con todas las librerías seguidas), pero conviene repetirlas
enchufado.

### Suite C++ (128 escenarios, 6 rondas intercaladas)

Nivel B frente a la fase 5: **1,14x de media geométrica** y 128 de 128 checksums iguales
(`bench/results/fase5/b-vs-5.md`). Por operación:

| Operación | Speedup |
|---|---:|
| `toeplitz_lower_vec` | 2,11x |
| `gcd` | 1,99x |
| `gcd_ext` | 1,91x |
| `fq.build` | 1,54x |
| `toeplitz_lower_inverse` | 1,49x |
| `is_irreducible` | 1,31x |
| `mul_fft` | 1,30x |
| `mul_auto` | 1,10x |
| `factor` | 1,04x |

Los escenarios de `gcd` y `gcd_ext` de grado 400-500 en Fp64 van de 3,4x a 4,3x, y
`mul_fft` de grado 4096 y 32 768 a 2,2x. En `factor` la mejora es pequeña porque los
escenarios de la suite son de grado 5 a 30, por debajo del umbral de la DDF nueva. Lo que
sale por debajo de 1 (hasta 0,87x) está en código que no ha cambiado, sobre todo
Karatsuba con hilos y aritmética de F_q. En la primera tanda de 4 rondas esos mismos
escenarios salieron distintos, así que se atribuye al ruido de medir con batería.

### Comparativa con python-flint (16 operaciones, un hilo, una misma tanda)

Distancia media a FLINT (media geométrica): **nivel A 6,0x → nivel B 2,7x → fase 5 1,66x**.
Las 16 operaciones dan el mismo resultado en todas las librerías. Datos en
`bench/results/fase5/comparativa_fase5.json`.

| Operación | Nivel B | Fase 5 | FLINT | Distancia a FLINT (B → 5) | Fase 5 frente a B |
|---|---:|---:|---:|---:|---:|
| Producto grado 1000, F_998244353 | 266 µs | **146 µs** | 167 µs | 1,6x → **0,87x** | 1,82x |
| Producto grado 10 000, F_998244353 | 10,64 ms | **3,81 ms** | 4,36 ms | 2,4x → **0,87x** | 2,79x |
| División 4000/2000, F_998244353 | 9,10 ms | 3,28 ms | 1,67 ms | 5,5x → 1,97x | 2,77x |
| MCD grado 1000, F_998244353 | 11,41 ms | 3,00 ms | 1,50 ms | 7,6x → 2,00x | 3,80x |
| Irreducibilidad grado 400, F_998244353 | 704 ms | **366 ms** | 586 ms | 1,2x → **0,62x** | 1,92x |
| Factorizar grado 300, F_998244353 | 352 ms | 111 ms | 61,8 ms | 5,7x → 1,79x | 3,18x |
| Factorizar grado 100, F_101 | 12,36 ms | 3,32 ms | 1,11 ms | 11,1x → 2,98x | 3,72x |
| Irreducibilidad grado 200, F_101 | 79,4 ms | 35,2 ms | 20,3 ms | 3,9x → 1,73x | 2,25x |
| Factorizar grado 40, F_101 | 556 µs | 466 µs | 173 µs | 3,2x → 2,69x | 1,19x |
| Producto grado 500, F_(2^127−1) | 741 µs | 742 µs | 801 µs | 0,93x → 0,93x | 1,00x |
| Factorizar grado 20, F_(2^127−1) | 47,8 ms | 47,7 ms | 6,97 ms | 6,9x → 6,85x | 1,00x |
| Producto grado 1000, F_256 | 2,57 ms | 2,61 ms | 1,95 ms | 1,3x → 1,34x | 0,98x |
| Factorizar grado 30, F_256 | 21,1 ms | 20,9 ms | 24,5 ms | 0,86x → 0,85x | 1,01x |
| Producto grado 200, F_(3^30) | 2,36 ms | 2,38 ms | 1,28 ms | 1,8x → 1,85x | 0,99x |
| Factorizar grado 5, F_(3^30) | 15,4 ms | 15,5 ms | 5,33 ms | 2,9x → 2,91x | 0,99x |
| Producto grado 100, F_998244353 | 15 µs | 16 µs | 9 µs | 1,7x → 1,81x | 0,91x |

Las filas que no cambian son extensiones, p grande o grados por debajo de los umbrales:
ninguna de las tres mejoras les afecta. El producto de grado 100 (0,91x) es una medida de
16 µs que incluye el coste de la llamada desde Python. En la suite C++, el mismo tamaño
sale igual.

## Verificación

- `build.py all`, con el intérprete del entorno `.venv`:
  - 9 tests C++, incluido el nuevo `tests/test_fase5.cpp` con 26 comprobaciones: NTT frente
    al escolar también con el primo 29·2^57 + 1, cercano a 2^62; Euclides en sitio frente
    a una copia del original en 6 cuerpos; half-GCD en Fp64 por encima de los nuevos
    umbrales; generador aleatorio en 5 primos; DDF nueva frente a la anterior en 7 cuerpos;
    `factor` reconstruye la entrada;
  - el oráculo frente a Python;
  - el módulo;
  - 710 tests de pytest con un 100 % de cobertura del paquete.
- Suite: 128 de 128 checksums iguales al nivel B.
- `python/` intacto.
- No se ha vuelto a medir la cobertura de líneas C++ (`build.py cov-cpp`).

## Pendiente y siguientes pasos

- **Half-GCD.** Su coste está en los productos de matrices (`aplicar`, `mul_mat`). Con el
  Euclides en sitio solo gana desde grado 2048-4096, y por eso el MCD de grado 1000 sigue a
  2x de FLINT. Se podría:
  - guardar las transformadas de las entradas de la matriz;
  - aplicar la matriz solo a la parte alta;
  - usar un caso base en sitio más grande.
- **Factorización en F_q y con p grande.** Grado 20 en F_(2^127−1) va a 6,9x de FLINT y
  grado 5 en F_(3^30) a 2,9x. A esos grados la DDF nueva no se activa. El coste está en la
  fase de grados iguales y en la aritmética de los elementos.
- Las propuestas 4 a 10 de la investigación: más NTT (sin permutación, radix 4, AVX2),
  producto medio, FFT truncada, Kronecker en varios puntos, F_2[x] empaquetado y
  logaritmos de Zech.
- Repetir la suite y la comparativa enchufado.
