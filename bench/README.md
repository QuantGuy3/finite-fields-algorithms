# Benchmarks de alcp

Harness reproducible para medir el rendimiento del port C++ (`include/alcp`) y del
módulo Python `alcp_native`, y para cuantificar cada optimización **contra el código
original** (etiqueta git `baseline-original`).

| Archivo | Para qué sirve |
|---|---|
| `bench_suite.cpp` | Suite C++ de 128 escenarios. Compila sin cambios contra las cabeceras originales y contra las mejoradas: solo usa API que ya existía en la versión original. |
| `run_bench.py` | Orquestador: `build`, `run`, `interleave`, `compare` y `module`. |
| `bench_python.py` | Benchmarks a nivel Python del módulo `alcp_native`, usando solo la API del binding original. |
| `bench.cpp` | Benchmark histórico original (se conserva tal cual). |
| `results/` | Resultados JSON versionados. Por ahora solo hay un ejemplo pequeño `smoke-*`, que **no es una medición**. |

Los ejecutables y módulos compilados van a `build/bench/<etiqueta>/`, que git ignora.

## Requisitos

- g++ con C++23 (probado con g++ 15.2 de MSYS2 ucrt64) y GMP/gmpxx (probado con 6.3).
- El intérprete del entorno virtual del proyecto (`.venv`), con `pybind11` para el subcomando `module`.

En los ejemplos, `PY` es el intérprete del entorno virtual:

```bash
PY=.venv/Scripts/python.exe        # Windows (Git Bash)
PY=.venv/bin/python                # Linux
```

## Metodología

### Entradas deterministas e idénticas entre versiones

- Las entradas de cada escenario se generan con un PRNG propio (splitmix64) con una
  semilla fija, que se deriva del nombre del escenario. Este PRNG es independiente
  del RNG de la librería: si una versión optimizada consume números aleatorios de
  otra forma, las entradas no cambian.
- También se llama a `set_seed(semilla fija)` antes de preparar el escenario y antes
  de cada repetición. Así la aleatoriedad interna de la librería (Cantor–Zassenhaus,
  raíz de la FFT, `Fq(p, n)`) sigue el mismo camino en todas las repeticiones de una
  misma versión.
- Los cuerpos F_q usan **módulos irreducibles fijos** escritos en el código. Se
  calcularon una vez con la versión original; el de F_(2^8) es el polinomio de AES.
- Las entradas se generan **antes** de cronometrar. Solo se mide la operación; el
  resultado de una repetición se destruye fuera del cronómetro.
- Cada escenario crea su propio cuerpo, con su propia caché de tablas FFT, para que
  los escenarios no dependan unos de otros.

### Repeticiones y estadísticos

- Se hace 1 calentamiento (`--warmup`). Después se repite hasta cumplir dos
  condiciones a la vez: al menos `--reps` repeticiones (7 por defecto) y al menos
  `--min-time` segundos acumulados (0,3 s por defecto). Hay un tope de `--max-reps`.
- Se guardan **todas** las muestras.
- El estimador principal es la **mediana**, porque es robusta frente a
  interrupciones del sistema operativo. También se informan el **mínimo** y el
  **IQR** (rango intercuartílico) como medida de dispersión.
- Speedup = mediana(base) / mediana(variante). Un valor mayor que 1 significa que la
  variante es más rápida.
- Se calcula la **media geométrica** de los speedups por grupo (cuerpo), por
  operación y global. Los escenarios con error de corrección se excluyen de las medias.
- Los tamaños están elegidos para que cada repetición tarde entre ~5 ms y ~1 s con la
  versión original compilada con `-O3`. Algunos escenarios ligeros repiten la
  operación varias veces dentro de cada repetición (campo `inner`).

### Checksums: misma salida en todas las variantes

- Cada escenario calcula un checksum FNV-1a de 64 bits de su resultado. Los
  coeficientes se vuelcan byte a byte en orden explícito, sin depender de la
  representación interna. `compare` exige que el checksum de cada variante coincida
  con el de la base; si no coincide, lo marca en **rojo** como **ERROR de corrección**
  y termina con código 1.
- **Factorizaciones**: el checksum no depende del orden de los factores. Se ordenan
  los hashes de cada par (factor, multiplicidad) antes de combinarlos, porque una
  versión puede devolver los factores en otro orden (p. ej. en orden canónico).
  Además se reconstruye `leading · Π uᵉ`; si no coincide con la entrada, el checksum
  lleva el prefijo `INVALIDO-` y se cuenta como error.
- **Elemento**: además del valor final se guarda una "foto" del acumulador cada
  1024 iteraciones, para que el checksum distinga resultados incluso en F_7.
- **`is_irreducible`**: el checksum incluye el polinomio de entrada. El irreducible se
  busca en la preparación, así que, si una versión defectuosa eligiera otro, se
  detectaría.
- **`fq.build.*`** (`Fq(p, n)` con módulo aleatorio): el checksum solo comprueba p, n
  y que el módulo sea mónico, de grado n e irreducible. Otra versión puede elegir otro
  irreducible igual de válido.
- Si el checksum de un escenario cambia entre rondas de `interleave`, se guarda como
  `VARIABLE:...` y se trata como error, porque indica no determinismo.

### Intercalado

En un portátil la frecuencia de la CPU varía con la temperatura. Medir primero todas
las variantes A y después todas las B introduce un sesgo sistemático. Por eso
`interleave` ejecuta A B C A B C… durante `--rounds` rondas y agrega las muestras de
todas ellas. Hay dos modos:

- **Por defecto**: cada turno ejecuta la suite completa (o la filtrada).
- **`--per-scenario`**: el intercalado es escenario a escenario. Es más fino, pero
  más lento, porque cada proceso repite la preparación del escenario.

### Recomendaciones para medir

- Portátil **enchufado a la corriente**: con batería, Windows limita la CPU. El JSON
  registra si hay corriente (`metadata.power`).
- Plan de energía de **alto rendimiento** o "máximo rendimiento". El JSON guarda la
  salida de `powercfg /getactivescheme`. Este harness **no cambia** el plan: se elige
  a mano en Configuración > Sistema > Energía, o con `powercfg /setactive <GUID>`.
- **Cerrar programas**: navegador, IDE indexando, OneDrive, actualizaciones y escaneos
  del antivirus. No compilar a la vez, ni siquiera otros agentes.
- Dejar enfriar el equipo unos minutos antes de empezar, y usar `interleave` con
  varias rondas (p. ej. `--rounds 5`) en lugar de una sola ejecución larga.
- Comparar siempre variantes compiladas con el **mismo compilador y la misma GMP**. El
  JSON guarda ambas versiones.
- Desconfiar de diferencias menores que el IQR, o de alrededor de un 3-5 %, sin
  repetir la medición.

## Uso

### `build`: compilar la suite contra unas cabeceras

```bash
$PY bench/run_bench.py build --include <DIR> --label <NOMBRE> [--cxxflags "-O3 ..."] [--ldflags "..."] [--cxx g++]
```

- Compila con `g++ -std=c++23 -Wall -Wextra <cxxflags> -I<DIR> bench/bench_suite.cpp <ldflags>`.
  El resultado es `build/bench/<NOMBRE>/bench_suite.exe`.
- `--ldflags` vale por defecto `-static -lgmpxx -lgmp` en Windows (ejecutable
  autocontenido) y `-lgmpxx -lgmp` en Linux.
- Guarda `build_info.json` con lo siguiente: flags, versión del compilador, ruta del
  include, commit git del include (`describe`, `dirty`), sha256 de las cabeceras y
  sha256 de la suite.
- `--cxxflags -O3` funciona sin comillas especiales: el script lo reescribe como
  `--cxxflags=-O3`.
- LTO: `--cxxflags "-O3 -flto" --ldflags "-flto -static -lgmpxx -lgmp"`.
- PGO (todavía no se ha probado con este harness): compilar con `-fprofile-generate`
  y ejecutar `run`; después, recompilar con la misma etiqueta y
  `-fprofile-use -fprofile-partial-training`. Puede requerir `-fprofile-dir=<dir>`
  para que GCC encuentre los `.gcda`.

### `run`: ejecutar una variante

```bash
$PY bench/run_bench.py run --label <NOMBRE> [--filter RE] [--reps N] [--min-time S] [--warmup N] [--out RUTA] [--note TEXTO]
```

Escribe `bench/results/<NOMBRE>.json`, o la ruta de `--out`. El JSON incluye las
muestras, los checksums y estos metadatos: CPU, núcleos físicos y lógicos, SO,
Python, plan de energía, corriente/batería, fecha, compilador, GMP, flags y commit
del include. Las rutas personales se anonimizan (`~`, `<usuario>`).

### `interleave`: variantes alternadas

```bash
$PY bench/run_bench.py interleave --labels A B C [--rounds 5] [--filter RE] [--reps 3] [--min-time 0] [--per-scenario] [--out-dir DIR]
```

Escribe un JSON por etiqueta con las muestras de todas las rondas
(`round_medians_s` guarda la mediana de cada ronda).

### `compare`: tabla comparativa

```bash
$PY bench/run_bench.py compare --base A --labels B C [--md salida.md] [--html salida.html] [--results-dir DIR]
```

- Muestra, por escenario, la mediana, el mínimo, el IQR y el speedup. También
  muestra las medias geométricas por grupo, por operación y global, y la verificación
  de checksums.
- `--base` y `--labels` aceptan etiquetas (que se buscan en `bench/results/`) o rutas
  a ficheros `.json`. Esto incluye la salida directa de `bench_suite --json` y los
  JSON de `bench_python.py`.
- `--html` genera un informe autocontenido con los errores en rojo.
- Termina con código 1 si hay algún error de corrección.

### `module`: compilar el binding pybind11

```bash
$PY bench/run_bench.py module --source <bindings.cpp> --include <DIR> --label <NOMBRE> [--cxxflags "-O3"]
```

- Genera `build/bench/<NOMBRE>/module/alcp_native<EXT_SUFFIX>` (p. ej.
  `.cp312-win_amd64.pyd`).
- En Windows enlaza g++ (MinGW) contra el CPython de MSVC con
  `-L<python>/libs -lpython312 -static -lgmpxx -lgmp`. El `.pyd` solo depende de
  DLL del sistema y de `python312.dll`.
- En Linux usa `-fPIC` y GMP dinámica.
- Al terminar, importa el módulo en un proceso aparte para comprobar que carga.

### `bench_python.py`: coste visto desde Python

```bash
$PY bench/bench_python.py --module-dir build/bench/<NOMBRE>/module --json bench/results/py-<NOMBRE>.json [--filter RE] [--reps N] [--min-time S] [--list]
```

Mide lo siguiente:

- Construcción de polinomios de grado 1000 desde una lista, `coeffs()` y la ida y
  vuelta.
- Productos de grado 1000 en F_998244353 con `*` y `.mul(b, "fft"|"karatsuba"|"school")`.
- `factor()` de grado 40 en F_101 y el caso histórico de grado 5 en F_(3^30).
- 20 000 llamadas a `Fp.mul` en un bucle Python (sobrecarga por llamada).
- Conversión de enteros grandes con p = 2^521 − 1 (grado 200).
- `is_prime`.
- `Fq(3, [1,0,1]).poly(...)`, `dft` y `toeplitz_lower_inverse`.

Funciona con el módulo compilado desde el `bindings.cpp` original. Desactiva el
recolector de basura durante cada medición. El JSON tiene el mismo formato que los
de `run_bench.py`, así que se compara con
`run_bench.py compare --base py-A.json --labels py-B.json`.

## Flujo típico: optimización frente a la línea base

```bash
# 1. Línea base: worktree con el código original (solo lectura)
git worktree add ../alcp-baseline baseline-original

# 2. Compilar las variantes
$PY bench/run_bench.py build --include ../alcp-baseline/include --label base-O3 --cxxflags -O3
$PY bench/run_bench.py build --include include --label nuevo-O3 --cxxflags -O3
$PY bench/run_bench.py build --include include --label nuevo-native --cxxflags "-O3 -march=native"

# 3. Medir intercalando (equipo en reposo, enchufado, alto rendimiento)
$PY bench/run_bench.py interleave --labels base-O3 nuevo-O3 nuevo-native --rounds 5

# 4. Comparar
$PY bench/run_bench.py compare --base base-O3 --labels nuevo-O3 nuevo-native \
    --md bench/results/comparativa.md --html bench/results/comparativa.html

# 5. Nivel Python
$PY bench/run_bench.py module --source ../alcp-baseline/python-bindings/bindings.cpp \
    --include ../alcp-baseline/include --label base-O3
$PY bench/bench_python.py --module-dir build/bench/base-O3/module --json bench/results/py-base-O3.json
```

## Catálogo de escenarios

Nombre: `<grupo>.<operación>.<tamaño>`. Usa `--list` para ver la descripción completa,
y `--filter` (expresión regular ECMAScript) para elegir un subconjunto.

| Grupo | Cuerpo | Escenarios (tamaños) |
|---|---|---|
| `fp.p7` | F_7 | `elem_muladd` n=4·10⁵, `elem_inv` 10⁵, `elem_pow` 10⁵; `mul_school` d64 (x20), d512; `mul_karatsuba` d512, d4096; `mul_auto` d2048; `divmod` 2000/1000; `gcd` d500; `gcd_ext` d400; `pow_mod` x^(7^48) mod f, f de grado 128; `is_irreducible` d64; `factor` d62 |
| `fp.p101` | F_101 | ídem; `pow_mod` 101^20; `is_irreducible` d48; `factor` d100 |
| `fp.ntt` | F_998244353 (NTT, 2^23 \| p−1) | ídem; además `mul_fft` d4096, d32768; `pow_mod` p^5; `is_irreducible` d32; `factor` d40 |
| `fp.m61` | F_(2^61−1) (sin NTT) | ídem; `pow_mod` d64, p^2; `is_irreducible` d24; `factor` d26 |
| `fp.m127` | F_(2^127−1) | `mul_karatsuba` d512, d2048; `pow_mod` d64, p; `is_irreducible` d16; `factor` d19 |
| `fp.p256` | F_r, r de 255 bits (BLS12-381, 2^32 \| r−1) | `mul_fft` d4096, d8192; `divmod` 1000/500; `pow_mod` d64, p; `is_irreducible` d12; `factor` d12 |
| `fq.q3e30` | F_(3^30), módulo fijo | `elem_muladd`, `elem_inv`, `mul_school` d32, `mul_karatsuba` d128, `divmod` 64/32, `factor` d5 (caso histórico) |
| `fq.q2e8` | F_(2^8), módulo AES | elemento, `mul_school` d128, `mul_karatsuba` d512, `is_irreducible` d16, `factor` d28 (característica 2: traza) |
| `fq.q7e4` | F_(7^4), módulo fijo | elemento, `mul_school` d128, `mul_karatsuba` d512, `gcd` d200, `factor` d20 |
| `fq.m61e2` | F_(p^2), p = 2^61−1, x^2+1 | elemento, `mul_school` d128, `mul_karatsuba` d512, `mul_fft` d1024, d4096, `mul_auto` d2048, `divmod`, `gcd`, `factor` d7 |
| `fq.build` | `Fq(p, n)` con módulo aleatorio | F_(3^30) y F_(7^20) |
| `fast` | `fast.hpp` | `dft` k=10, 11, 12; `idft` k=12, con la raíz fija 3^((p−1)/2^k); `toeplitz_lower_inverse` n=1024, 8192 (NTT) y n=1024 (F_7); `toeplitz_lower_vec` n=4096; `toeplitz_divmod` 2000/1000 (NTT y F_7) |

Las entradas de `factor` son productos `c · Π uᵢ^{mᵢ}` de polinomios mónicos
aleatorios, con factores repetidos, de modo que se ejercitan `radical`,
`distinct_degree`, `equal_degree` y las multiplicidades. Las de `gcd` y `gcd_ext`
comparten un factor común de grado ~d/8, para que el mcd no sea trivial.

## Formato de los resultados

`bench/results/<etiqueta>.json` tiene estos campos:

- `format`, `label`, `kind` (`cpp` o `python`), `mode` (`single` o `interleave`),
  `note`.
- `metadata`: máquina, energía, `build` (el `build_info.json` del build) y `suite`
  (compilador y GMP que informa el propio ejecutable).
- `runs`: una entrada por ejecución o ronda.
- `scenarios[]`, con estos campos para cada escenario:
  - `name`, `group`, `op`, `desc`, `inner`;
  - `samples_s` (todas las muestras, en segundos), `median_s`, `min_s`;
  - `round_medians_s`, `setup_s`;
  - `checksum`, `error`.

## Limitaciones

- La medición depende del estado de la máquina. Los resultados `smoke-*` se
  obtuvieron con el equipo ocupado y solo validan que todo funciona; no sirven para
  sacar conclusiones de rendimiento.
- En `factor`, `fq.build.*` y la primera FFT de cada cuerpo, el tiempo depende del
  camino aleatorio de la librería. Es idéntico entre repeticiones de una misma versión,
  pero una versión que consuma el RNG de otra forma recorrerá otro camino, igual de
  válido. Esa parte de la diferencia de tiempo no se debe a la optimización.
- La búsqueda del irreducible de los escenarios `is_irreducible` se hace en la
  preparación (no cronometrada) y puede tardar algunos segundos por escenario.
- El checksum FNV-1a de 64 bits detecta diferencias accidentales, pero no es
  criptográfico.
- Las medias geométricas resumen escenarios heterogéneos. Para conclusiones
  concretas, mire la tabla por escenario y el IQR.
