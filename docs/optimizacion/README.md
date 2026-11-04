# Fase 3: optimización de velocidad — informe final

Rama `mejoras-cpp-pybind`, medición final sobre el commit `32d6922` (4 de noviembre de 2026).
Referencia: etiqueta `baseline-original` (67331a2), es decir, el código original
compilado con `-O3`. Todas las cifras salen de los JSON de `bench/results/fase3/final-*.json`.
Las tablas completas por escenario están en [final.md](final.md) y [final.html](final.html)
(suite C++) y en [final-python.md](final-python.md) (módulo Python). La medición técnica a
técnica, hecha antes de la integración, está en [mediciones.md](mediciones.md).

Regla del encargo: optimizar **sin cambiar la lógica de los algoritmos**. Escolar,
Karatsuba, NTT radix-2, Euclides, Ben-Or, Cantor–Zassenhaus e inverso de Fermat en F_q
siguen siendo los mismos algoritmos, con las mismas operaciones matemáticas. Solo cambian
la representación de los datos, la reducción modular, la gestión de memoria, los umbrales
de despacho entre algoritmos que ya existían y el reparto del trabajo entre hilos.

## 1. Resumen ejecutivo

- **Corrección.** En los 128 escenarios de `bench_suite`, el checksum de la versión final
  coincide con el de la referencia en todas las variantes medidas: `-O3`, `-march=native`,
  `--backend gmp` y `ALCP_THREADS=1`. En el módulo Python coinciden los 19 escenarios.
  `build.py all` termina sin errores.
- **Suite C++ completa: 24,4x** de media geométrica (128 escenarios, 4 rondas
  intercaladas). La suma de las medianas de los 128 escenarios baja de **16,0 s a 1,15 s**.
- **Sin el backend de 64 bits** (`--backend gmp`, todo con GMP): **4,7x**. **Con un solo
  hilo** (`ALCP_THREADS=1`): **21,5x**. El paralelismo aporta, por tanto, un 13,5 % de la
  media global, aunque mucho más en Karatsuba y en la FFT grande.
- **Módulo Python** (`bench_python.py`, 19 escenarios, módulo original frente al final):
  **9,7x**.
- **Caso histórico** (factorización de grado 5 sobre F_(3^30)): 439 ms → **11,6 ms** con el
  backend de 64 bits (37,9x) y 93 ms todo con GMP (4,7x). Ver la sección 3.
- **`-march=native`**: 24,6x frente a 24,4x (+0,9 %, dentro del ruido). Se mantiene el
  perfil **portable** por defecto (ver la sección 6).

### Speedup por grupo (media geométrica, `final-O3` frente a la referencia)

| Grupo | Cuerpo | Final (`-O3`) | Solo GMP | 1 hilo | `-march=native` |
|---|---|---:|---:|---:|---:|
| `fp.p7` | F_7 | 61,6x | 3,9x | 52,2x | 62,6x |
| `fp.p101` | F_101 | 60,4x | 4,1x | 52,9x | 61,2x |
| `fq.q3e30` | F_(3^30) | 58,6x | 7,0x | 49,3x | 59,2x |
| `fp.ntt` | F_998244353 | 52,3x | 7,2x | 46,5x | 53,1x |
| `fq.q2e8` | F_(2^8) | 48,6x | 6,3x | 43,0x | 49,6x |
| `fq.build` | `Fq(p, n)` aleatorio | 46,2x | 3,2x | 46,6x | 45,4x |
| `fp.m61` | F_(2^61−1) | 43,4x | 3,7x | 39,0x | 44,4x |
| `fq.q7e4` | F_(7^4) | 39,1x | 6,2x | 35,8x | 38,8x |
| `fast` | `fast.hpp` (NTT y F_7) | 26,3x | 7,7x | 25,7x | 26,5x |
| `fq.m61e2` | F_((2^61−1)^2) | 17,7x | 4,7x | 16,7x | 17,9x |
| `fp.p256` | F_r, r de 255 bits | 3,8x | 3,8x | 3,2x | 3,7x |
| `fp.m127` | F_(2^127−1) | 3,4x | 3,5x | 2,8x | 3,5x |
| **global** | 128 escenarios | **24,4x** | **4,7x** | **21,5x** | **24,6x** |

Los grupos `fp.m127` y `fp.p256` tienen p ≥ 2^62, así que no usan el backend de 64 bits:
su mejora (3,4-3,8x) viene de GMP en sitio, la reducción perezosa, la NTT multi-limb y el
paralelismo. Las columnas «Solo GMP» y «1 hilo» se midieron justo después del intercalado
principal, pero no intercaladas con la referencia (ver la sección 7).

### Por operación (final `-O3`)

| Operación | Speedup | Operación | Speedup |
|---|---:|---|---:|
| `mul_karatsuba` (16) | 69,6x | `factor` (10) | 19,2x |
| `divmod` (8) | 58,6x | `gcd_ext` (6) | 19,0x |
| `mul_school` (16) | 57,4x | `mul_fft` (6) | 14,2x |
| `toeplitz_lower_inverse` (3) | 50,5x | `idft` / `dft` | 14,0x / 12,8x |
| `toeplitz_lower_vec` (1) | 44,2x | `gcd` (8) | 12,0x |
| `pow_mod` (6) | 35,3x | `elem_muladd` (10) | 10,5x |
| `mul_auto` (7) | 34,7x | `elem_inv` (10) | 8,0x |
| `toeplitz_divmod` (2) | 31,0x | `elem_pow` (6) | 2,6x |
| `is_irreducible` (7) | 29,9x | | |

### Los 10 escenarios que más mejoran

| Escenario | Referencia | Final | Speedup |
|---|---:|---:|---:|
| `fp.p101.mul_karatsuba.d4096` | 339,5 ms | 1,45 ms | 234x |
| `fp.p7.mul_karatsuba.d4096` | 303,1 ms | 1,37 ms | 221x |
| `fp.m61.mul_karatsuba.d4096` | 342,6 ms | 1,65 ms | 208x |
| `fp.ntt.mul_karatsuba.d4096` | 345,1 ms | 1,67 ms | 206x |
| `fp.ntt.gcd_ext.d400` | 441,1 ms | 2,42 ms | 182x |
| `fq.q3e30.mul_karatsuba.d128` | 414,1 ms | 2,29 ms | 181x |
| `fp.p101.mul_auto.d2048` | 111,0 ms | 0,71 ms | 157x |
| `fq.q2e8.mul_karatsuba.d512` | 347,9 ms | 2,34 ms | 149x |
| `fp.ntt.mul_school.d512` | 28,7 ms | 0,20 ms | 144x |
| `fp.p101.mul_school.d512` | 28,2 ms | 0,20 ms | 142x |

### Los que menos mejoran

| Escenario | Referencia | Final | Speedup | Comentario |
|---|---:|---:|---:|---|
| `fp.m127.elem_muladd.n200000` | 22,7 ms | 23,9 ms | 0,95x | Las 4 medianas por ronda del final (23,4-25,3 ms) quedan todas por encima de las de la referencia (22,3-23,1 ms). Probablemente es una pérdida real de un ~5 %, pequeña, en `Fp` (GMP) con un solo elemento |
| `fp.p256.elem_muladd.n200000` | 31,7 ms | 32,5 ms | 0,97x | Rangos por ronda solapados: ruido |
| `fp.m127.elem_inv.n50000` | 31,4 ms | 31,5 ms | 0,99x | El inverso de GMP (`mpz_invert`) no cambia |
| `fp.p256.elem_pow.n2000` | 19,3 ms | 19,4 ms | 0,99x | Ídem `mpz_powm` |
| `fp.p256.elem_inv.n50000` | 65,7 ms | 65,9 ms | 1,00x | |
| `fp.m127.elem_pow.n10000` | 23,1 ms | 23,0 ms | 1,00x | |
| `fp.m61.elem_pow.n20000` | 8,2 ms | 6,2 ms | 1,32x | `elem_pow` en 64 bits está limitado por la cadena secuencial de cuadrados |
| `fp.m127.gcd.d400` | 20,9 ms | 14,6 ms | 1,42x | |

En las operaciones de un solo elemento con p ≥ 2^62, el coste ya estaba en GMP
(`mpz_mul`, `mpz_invert`, `mpz_powm`) y la versión final no lo cambia.

## 2. Qué se integró

Orden de integración: toolchain → multiplicación → aritmética → capa Python, y después un
commit propio de integración. En la columna «medida», «serie» indica la medición técnica
a técnica de [mediciones.md](mediciones.md), con corriente y cada técnica sobre la
anterior de su rama. Las cifras de esa columna no son acumulativas entre técnicas.

### 2.1 Aritmética de cuerpos (rama `worktree-wf_e9ad1fab-e9f-1`)

| Commit | Técnica | Medida (serie) |
|---|---|---|
| `489728b` | GMP en sitio en `Fp`/`Fq` (API C de GMP, sin temporales de gmpxx; búfer reutilizado en `Fq::mul`) | 2,01x en los escenarios afectados (`fq.*` 2-4x) |
| `0de8cf9` | Reducción perezosa en `mul_school` y `divmod` (concepto opcional `LazyAccumulation`: se acumula sin reducir y se reduce una vez por coeficiente) | 4,68x (hasta 13x en `mul_school`/`divmod`) |
| `dcb905c` | Suma y resta de `Poly` sin relleno de ceros ni copias | 1,17x focalizada (Karatsuba 1,06-1,50x; neutra en `gcd`) |
| `65a5f1c` + `5c444c0` | Backend de palabra de máquina `Fp64` (Barrett, p < 2^62) y `FqT<Base>` (`Fq64 = FqT<Fp64>`); selección `--backend` en la suite | 8,69x; 34x en los afectados |
| `4fb4f1a` | Binding: el módulo elige `u64` o GMP automáticamente (`std::variant`) | Módulo 5,83x frente al de 5f49217 |
| `2297387` | Barrett de `Fp64` sin desplazamientos variables de 128 bits | 1,07x focalizada (solo `elem_*`: 1,07-1,23x) |
| `8f1d896` | `Fp64::pow` con intermedios en forma de Montgomery (misma cadena de cuadrados) | `elem_pow` 1,60x |

### 2.2 Multiplicación (rama `worktree-wf_e9ad1fab-e9f-2`)

| Commit | Técnica | Medida (serie) |
|---|---|---|
| `d4580af` | FFT: tablas de acceso directo, inversión de bits precalculada, twiddles por etapa contiguos | 1,21x |
| `f68a18b` | Cuadrado con una sola transformada directa | Neutra en la suite (1,02x); −31 % aislada. Se mantiene: ver la sección 4 |
| `9bfc1e7` | Karatsuba sobre vistas (puntero + longitud), sin `slice`/`shift` | 1,18x focalizada |
| `1ad8c02` | NTT en `uint64_t` con Shoup para p < 2^62 | 4,51x (hasta 44x) |
| `70e3453` | NTT multi-limb (`mpn`, Montgomery) para p ≥ 2^62 | 2,45x |
| `3660189` | Umbrales de despacho Auto por backend y por operando menor | 1,02x de media; `fp.p256.gcd_ext` 2,0x y `fp.ntt.factor` 1,6x |
| `06e46b5` | `fast::dft/idft` en palabras de máquina | 11,6x |
| `e8e07e6` | Karatsuba en paralelo (`std::async`) | 1,83x en la cadena; `mul_auto` de F_p 2,9-3,1x |
| `b0a9fe7` | NTT paralela (mariposas genéricas y transformadas directas concurrentes) | 1,27x en la cadena |

### 2.3 Integración (commits `f30bafd` y `32d6922`)

Las dos ramas C++ se desarrollaron y midieron por separado. Al combinarlas aparecieron
interacciones que no se veían en ninguna de las dos por separado:

| Cambio | Por qué | Efecto medido (1 ejecución de 5 repeticiones, mediana) |
|---|---|---|
| `Fp64` usa la NTT en palabras (concepto `NativeWordField` en `fft.hpp` y `fast.hpp`) | Sin esto, `Fp64` caía en la NTT genérica elemento a elemento | Necesario para que `fp.ntt.mul_fft` y `fast.ntt.*` conserven las ganancias de `1ad8c02`/`06e46b5` con el backend `u64` |
| Umbral Auto propio de `Fp64`: FFT si small ≥ 128, o small ≥ 48 y large ≤ 4·small | Con el escolar de 128 bits, la regla «Word: L ≥ 4» de `Fp` elegía la FFT para productos diminutos | `fp.ntt.gcd_ext.d400` 11,7 → 2,47 ms; `fp.ntt.factor.d40` 2,05 → 1,11 ms |
| `kParKaratsubaMin = 1024` para elementos de 64 bits | Con 64 se creaban hilos para productos de microsegundos | `fp.p7.pow_mod` 82,9 → 4,7 ms; `fp.p7.is_irreducible` 53,2 → 3,7 ms |
| `kParChunk = 1024` para elementos de `Fq64` | La reserva de memoria de cada elemento escala mal con 16 hilos | `fq.m61e2.mul_fft.d1024` 32,6 → 11,9 ms |
| El búfer de `FqT::mul` deja de ser `thread_local` | Con MinGW-w64 y enlace estático, un `thread_local` con destructor corrompe el montón al terminar los hilos de `std::async` (0xC0000374, reproducido con un programa mínimo). `fq.m61e2.mul_fft` **abortaba** | Ahora hay un búfer global con indicador atómico, y los hilos que no lo consiguen usan un acumulador local |
| `-fprofile-update=atomic` en `cov-cpp` | Los contadores de cobertura se corrompían en las regiones paralelas (gcovr: «negative hit value») | — |

Hay tests nuevos en `tests/test_fp64.cpp`, sección 7: NTT de `Fp64` alrededor de los
umbrales, `fast::dft/idft`, NTT multi-limb paralela y `Fq64` con FFT y Karatsuba
paralelos, siempre frente a GMP.

### 2.4 Capa Python (rama `worktree-wf_e9ad1fab-e9f-4`)

| Commit | Técnica | Medida (serie, mismo `.pyd`) |
|---|---|---|
| `abc2f4d` | Camino rápido para operandos `int` en `+ - *` | 1,15x |
| `27873ea` | Analizador: `x^k` como monomio en lugar de potenciar (no toca ningún algoritmo del núcleo) | 1,69x; `parse_long` 5,1x |
| `29d5628` | Caché de cadenas ya analizadas | 94x (solo con cadenas repetidas) |
| `0ab679a` | `_pretty`: cada coeficiente distinto se formatea una vez | 1,47x |
| `a3399d9` | Caché de textos de elementos de una extensión | 5,19x |
| `cab2c2a` | Pickle de elementos más directo | 1,22x |
| `e683a88` | Listas de coeficientes con caminos rápidos | 1,28x |
| `c70bed5` | `_pretty`: potencias en línea | 1,12x |
| `8fcc151` | Métodos nativos ligados en `__slots__` | Neutra (1,03x). Se mantiene porque `abc2f4d` depende de ella |

La capa completa, de 5f49217 a c70bed5, da **2,05x** (41 escenarios `alcp.*`, mismo
`.pyd`; ver [mediciones.md](mediciones.md)). Con el módulo final, `bench_alcp.py`
(`final-alcp.json`) da estos cocientes capa/nativo: operaciones de elementos de F_p de
1,55x a 1,9x (la llamada pybind11 cuesta unos 1,3 µs), F_9 y GF(2^8) de 1,1x a 1,3x,
construcción desde lista 1,6-1,8x, producto de grado 100 1,04x, `factor`/`roots` 1,03-1,05x
e impresión 5-8x. El analizador es 0,08-0,6x con cadenas repetidas (gracias a la caché) y
147-403x con cadenas nuevas.

### 2.5 Toolchain

Solo se integró `2e5cc5c`: la opción `--perfil portable|v3|nativo` de `build.py`, con
`portable` por defecto. El resto de la rama no se integró (sección 4).

## 3. Casos históricos: antes y después

Se usan los casos de `bench/bench.cpp`: factorizar un polinomio aleatorio de grado 5 sobre
F_(3^30) y multiplicar dos polinomios de grado 4000 en F_998244353. Se miden con
`bench/bench_historico.cpp`, que hace las mismas llamadas con la misma semilla, pero
repite cada medida (mediana de 5 repeticiones por ronda, 2 rondas alternadas, `-O3`) y
añade el backend de 64 bits. Los resultados están en `bench/results/fase3/final-historico.json`.

| Caso | Antes (cifras históricas) | Referencia hoy | Final, GMP (`Fp`/`Fq`) | Final, 64 bits (`Fp64`/`Fq64`) | Speedup (ref. → 64 bits) |
|---|---:|---:|---:|---:|---:|
| `factor` grado 5, q = 3^30 | 823 ms | 439,5 ms | 93,0 ms | **11,6 ms** | 37,9x |
| Producto grado 4000, escolar | 3752 ms | 1750,4 ms | 160,6 ms | **12,8 ms** | 137x |
| Producto grado 4000, Karatsuba | 606 ms | 331,8 ms | 20,5 ms (1 hilo: 81,8) | **1,54 ms** (1 hilo: 4,07) | 215x |
| Producto grado 4000, FFT | 62 ms | 33,8 ms | 1,12 ms | **0,72 ms** | 47x |
| Producto grado 4000, Auto (`*`) | — | 33,5 ms | 1,07 ms | **0,73 ms** | 46x |

- La columna «Antes» son las cifras históricas que cita el encargo. **No se han vuelto a
  medir** y sus condiciones (batería, flags, carga) no constan. La columna «Referencia
  hoy» es el código original medido ahora en las mismas condiciones que la versión final.
  Es la comparación válida.
- En la suite, el escenario equivalente con módulo fijo, `fq.q3e30.factor.d5`, pasa de
  424,3 ms a **10,83 ms** (39x).
- `bench.cpp` sin cambios (una sola medida, tipos GMP) da 440 → 92 ms (factor) y
  1748/328/36 → 159/21/2,7 ms (escolar/Karatsuba/FFT).

## 4. Qué no funcionó, o no se integró, y por qué

| Técnica | Resultado (serie) | Decisión |
|---|---|---|
| `0438b74` Resultado creado en línea en `+ - * /` (capa Python) | 1,02x con control 0,97x: ruido. Domina la llamada pybind11 | **Revertida** (`5ce76db`): duplicaba código |
| `06d922a` Caminos rápidos para `// % divmod` | 1,04x con grado 100 | **Revertida** (`87f4a27`) |
| `a833bcb` Caminos rápidos int/elemento en la conversión y en `f(v)` | 0,98x | **Revertida** (`5d0556f`) |
| `8fcc151` Métodos nativos ligados | 1,03x (control 1,03x) | Se mantiene: es requisito de `abc2f4d`, que funciona (revertirla entra en conflicto) |
| `f68a18b` Cuadrado con una sola transformada | Neutra en la suite (domina `divmod`); −31 % en un microbenchmark de cuadrados | Se mantiene: los caminos en palabras y multi-limb posteriores (`1ad8c02`, `70e3453`, `b0a9fe7`) se apoyan en esa rama del código, y evita trabajo idéntico sin coste |
| `-march=x86-64-v3` / `-march=native` | 1,01-1,02x en la serie; **+0,9 %** en la medición final (24,57x frente a 24,36x; por escenario entre 0,89x y 1,09x) | No por defecto (sección 6) |
| LTO (`06ebf59`) | 0,99-1,02x: hay una única unidad de traducción | No integrado |
| Asignador rápido `--fast-alloc` (`584b9aa`, wip) | 1,23-1,25x en 19 escenarios, pero el código está sin terminar y sin tests | **No integrado**; pendiente |
| PGO | No se midió: no había un flujo utilizable | Pendiente |
| Tokenizador con regex (`exp/analizador-tokenizador-regex`, `ff9ad21`) | Neto ≈ +3 %, dentro del ruido: el coste está en las operaciones con polinomios de cada término | Descartado (rama `exp/`) |
| AVX2 en la NTT de `Fp64` | No implementado | Pendiente |
| GCD binario para `elem_inv` en 64 bits | Cambiaría el algoritmo (Euclides extendido) | Descartado por la regla del encargo |

## 5. Coste de las mejoras de robustez (fases 1 y 2)

5f49217 (final de la fase 2, con validaciones, excepciones, `factor()` canónico, etc.)
frente a `baseline-original`, ambas con `-O3`, en la suite completa: **1,04x** con
corriente (3 rondas) y 1,08x con batería. Las mejoras de robustez **no costaron
velocidad**. Nueve escenarios quedan entre 0,85x y 0,95x, todos con un IQR del mismo orden
que la diferencia (el peor, `fp.ntt.elem_muladd`, 0,85x con un IQR del 44 %). Ver
[mediciones.md](mediciones.md).

## 6. Flags por defecto

- **Perfil por defecto de `build.py`: `portable`** (x86-64 genérico, `-O3` en el módulo y
  `-O2` en los tests, sin LTO). `-march=native` da +0,9 % de media geométrica en la
  medición final, dentro del ruido, y además el binario resultante solo funciona en CPU
  compatibles con la que compila. El `.pyd` portable funciona en cualquier PC x86-64. El
  único perfil con ganancia clara (el asignador wip) no está terminado.
- Quien compile solo para su máquina puede usar `build.py --perfil nativo` (o `v3`):
  los checksums no cambian.
- `cov-cpp` añade `-fprofile-update=atomic` (sección 2.3).

## 7. Metodología y condiciones

| Aspecto | Valor |
|---|---|
| Equipo | AMD Ryzen 7 7730U (Zen 3, 8 núcleos/16 hilos), Windows 11 Home 10.0.26300 |
| Compilador / GMP | g++ 15.2 (MSYS2 ucrt64), GMP 6.3, `-std=c++23 -Wall -Wextra`, `-static -lgmpxx -lgmp` |
| Energía | **Enchufado** (batería al 99 %, `metadata.power` en los JSON), plan **Equilibrado** (no se cambió) |
| Carga | Ningún otro proceso pesado durante la medición (este proceso era el único que compilaba o medía) |
| Suite C++ | `run_bench.py interleave --labels final-ref-O3 final-O3 final-native --rounds 4 --reps 5 --min-time 0.05`: 4 rondas A B C intercaladas, ≥ 5 repeticiones y 1 calentamiento por ronda. Se agregan todas las muestras y se usa la mediana |
| Mismo arnés | Las tres variantes usan el mismo `bench_suite.cpp` (sha256 en `build_info.json`). La referencia se compiló contra `git archive baseline-original include` (directorio extraído sin metadatos git; hay un sha256 de las cabeceras en el JSON) |
| Variantes extra | `--backend gmp` y `ALCP_THREADS=1`: 3 rondas alternadas entre sí justo después del intercalado principal, unidas con `une_rondas.py`. **No** se intercalaron con la referencia |
| Módulo Python | `bench_python.py`, módulo original (`bindings.cpp` y cabeceras de `baseline-original`) frente al final, 4 rondas alternadas a mano y unidas por escenario |
| Ruido | IQR relativo mediano por escenario: 4,3 % en la referencia, 4,7 % en la final y 4,8 % en la nativa. 23 escenarios de la final superan el 10 % (casi todos de menos de 1 ms, donde se nota el planificador del SO). Las medias geométricas son robustas a eso; las diferencias por escenario por debajo del ~5 % no se consideran significativas |
| Checksums | FNV-1a de 64 bits por escenario; `compare` termina con código 0 en todos los conjuntos |

Archivos: `bench/results/fase3/final-ref-O3.json`, `final-O3.json`, `final-native.json`,
`final-O3-gmp.json`, `final-O3-t1.json`, `final-py-ref.json`, `final-py.json`,
`final-alcp.json` y `final-historico.json`.

### Verificación

- `build.py all`: 6 tests C++ (`test_all` 171 122 comprobaciones, `test_fp64` 868,
  `test_robustez` 410, `smoke` 15 + 242, `xcheck` 924), oráculo, módulo y 709 tests de
  pytest con cobertura del paquete del **100 %**.
- `build.py cov-cpp`: **94,4 % de líneas** (5345/5664), 98,4 % de funciones y 87,6 % de
  ramas. `bindings.cpp` está al 100 %. Antes de la fase 3 era el 98,3 % de líneas: gcovr
  cuenta las líneas por instanciación de plantilla, y los ayudantes de paralelismo y
  multi-limb (`fft.hpp`, 78 %) tienen instanciaciones que los tests no recorren (por
  ejemplo, la lectura de `ALCP_THREADS` o la NTT paralela para tipos de elemento que no
  llegan a ese tamaño).
- `python-bindings/ejemplos.py --sin-color`: correcto. `python/` intacto
  (`git diff --quiet baseline-original HEAD -- python/`).

## 8. Limitaciones

- **Una sola máquina.** Los umbrales de despacho (`prefer_fft`, el umbral de `Fp64`) y los
  de paralelismo (`kParKaratsubaMin`, `kParChunk`, `kParWordMin`, `kParLimbMin`) se
  calibraron en este Zen 3 con 16 hilos. En otra CPU, con otro número de núcleos o en
  Linux, el punto óptimo puede ser distinto. Los umbrales de la integración se fijaron con
  medidas de 1 ejecución.
- **El paralelismo aumenta el tiempo de CPU total.** El speedup de Karatsuba y de la FFT
  grande es de tiempo de pared. Con `ALCP_THREADS=1` la media global es 21,5x.
- **Las cifras de la sección 2 no son aditivas.** Se midieron técnica a técnica en cada
  rama, antes de integrar. La cifra global válida es la de la sección 1.
- **Plan de energía Equilibrado.** No se cambió la configuración del sistema. El
  intercalado compensa la deriva térmica, pero los tiempos absolutos pueden variar unos
  puntos con otro plan.
- **Columnas «Solo GMP» y «1 hilo».** No están intercaladas con la referencia: tienen un
  error algo mayor que la columna principal.
- **Pequeña pérdida probable** en `fp.m127.elem_muladd` (~5 %, `Fp` con GMP). El resto de
  operaciones de un elemento con p ≥ 2^62 quedan igual.
- **Smart App Control** bloquea a veces binarios recién compilados. Se vuelven a enlazar
  y no afecta a los tiempos.
- El análisis es **descriptivo**: atribuir la mejora de un escenario a una técnica
  concreta se apoya en la medición técnica a técnica, no en la final.

## 9. Recomendaciones

1. Terminar y probar el asignador rápido (`--fast-alloc`): es la única técnica de
   toolchain con ganancia clara (1,23-1,25x en el subconjunto medido).
2. Volver a calibrar los umbrales en reposo con un plan de alto rendimiento y en otra
   máquina (Linux, otro número de núcleos), sobre todo `kParKaratsubaMin` y el umbral
   FFT de `Fp64`.
3. Revisar `Fp::mul` con GMP en un solo elemento (`fp.m127.elem_muladd`, ~5 %).
4. Evitar `thread_local` con destructores no triviales en todo el código que pueda
   ejecutarse en las regiones paralelas, mientras se use MinGW con enlace estático.
5. Opcional: AVX2 en la NTT de `Fp64` y un flujo PGO reproducible.

Todas las cifras son de este equipo y de esta sesión. Antes de reutilizarlas conviene
revisarlas y, si se cambia de entorno, repetir la medición con `bench/run_bench.py`.
