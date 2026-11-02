# Fase 3: mediciones en serie de las optimizaciones de velocidad

Etapa de medición serial (4 de noviembre de 2026, 15:24-18:11). Solo se compilaron y
midieron variantes; no se modificó código fuente. Los JSON con todas las muestras están en
`bench/results/fase3/` y los scripts de la campaña, con las rutas anonimizadas, en
`bench/results/fase3/scripts/`.

## Resumen

- **Corrección.** Todas las variantes medidas (27 compilaciones de `bench_suite`, 4 módulos
  `alcp_native` y 14 versiones de la capa Python) dan **los mismos checksums** que la línea
  base original. `run_bench.py compare` termina con código 0 en los 20 conjuntos de
  resultados. No hay ninguna técnica con veredicto «error».
- **Coste de robustez (5f49217 frente a `baseline-original`, ambas con `-O3`).** Media
  geométrica de **1,04x** en los 128 escenarios con corriente (1,08x con batería). Es decir,
  las fases 1-2 **no** costaron velocidad. Hay 9 escenarios entre 0,85x y 0,95x, todos
  con un IQR del mismo orden que la diferencia.
- **Rama aritmética (HEAD 8f1d896), suite completa frente a la referencia:** **18,8x** con
  el backend automático de 64 bits y **2,84x** con `--backend gmp`.
- **Rama multiplicación (HEAD b0a9fe7):** **1,70x** con todos los hilos y **1,40x** con
  `ALCP_THREADS=1`. Las ganancias se concentran en `mul_fft`, `mul_auto`, `mul_karatsuba`,
  `gcd_ext`, `dft` y Toeplitz.
- **Módulo Python** (`bench_python.py`, 19 escenarios): binding de la rama aritmética
  **8,97x** frente al módulo original y **5,83x** frente al de 5f49217.
- **Capa Python** (`bench_alcp.py`, 41 escenarios `alcp.*`, mismo `.pyd`): **2,05x** de
  media geométrica de 5f49217 a c70bed5. Lo que más pesa son las cachés y el analizador;
  en las operaciones de elementos la ganancia es del 4 al 38 %.
- **No funcionaron o son neutras:** `f68a18b` (cuadrado con una transformada: invisible en
  la suite), `8fcc151`, `a833bcb`, `06d922a` y **`0438b74`** (capa Python; este último
  conviene revertirlo), la rama experimental `exp/analizador-tokenizador-regex`
  (`ff9ad21`) y los perfiles de toolchain `-march=x86-64-v3`, `-march=native` y LTO
  (±2 %). La técnica `2297387` (Barrett sin desplazamientos variables) solo gana en las
  operaciones de elementos (+7 % de media en la medición focalizada).
- **El asignador rápido** de la rama wip de toolchain (`--fast-alloc`, sin terminar) da
  **1,23-1,25x** en el subconjunto medido. Es la única técnica de toolchain con efecto
  claro, pero no tiene tests: está pendiente.

## Condiciones de medición

| Aspecto | Valor |
|---|---|
| Equipo | AMD Ryzen 7 7730U (Zen 3, 8 núcleos/16 hilos), Windows 11 Home 10.0.26300 |
| Compilador / GMP | g++ 15.2 (MSYS2 ucrt64), GMP 6.3, `-std=c++23 -Wall -Wextra -O3`, `-static -lgmpxx -lgmp` |
| Plan de energía | **Equilibrado** (`381b4222-…`, «Balanced»). No se cambió: modificar la configuración del sistema queda fuera del alcance de esta etapa |
| Corriente | **Batería hasta las 17:06:29** (evento Kernel-Power 105, «Power source change»); **enchufado** a partir de entonces. Con batería, la CPU rendía ~2,3 veces menos: la suite de la referencia tardaba 199 s y con corriente 85 s |
| Suspensiones | Entre las 15:47 y las 16:32 la tapa se cerró varias veces y el equipo entró en **Modern Standby** (Kernel-Power 506/507, «Reason: Lid»). Las rondas afectadas se descartaron (ver abajo) |
| Carga | Ningún otro proceso pesado (los agentes de desarrollo habían terminado). Carga de base del 10-18 % (navegador, Google Drive, la propia aplicación) |
| Rondas | Extremos con batería: 3 rondas (1-2 válidas por variante); extremos con corriente: 3 rondas. Cadenas: 3 rondas con batería y 3 (aritmética) o 2 (multiplicación) con corriente. Focalizadas: 6 rondas × ≥ 5 repeticiones (`--min-time 0.2`). Módulo: 6 rondas. Capa Python: 3 rondas |
| Repeticiones | `--reps 3 --min-time 0` (1 calentamiento) por ronda en las cadenas y los extremos |
| Ruido | IQR relativo mediano por escenario: extremos con corriente 11,7 %, con batería 15,4 %; cadena aritmética con corriente 18,8 %; cadena de multiplicación 12,1 %; focalizadas 7,9-10,8 %; toolchain 7,7 %; capa Python 13,3 %. **El ruido supera el 10 % de la mediana en muchos escenarios**, por eso se repitieron las técnicas dudosas en mediciones focalizadas de 6 rondas |

### Cómo se trató el ruido

1. **Detección de suspensiones.** La primera tanda de extremos (15:24-16:38) coincidió con
   cinco intervalos de Modern Standby: una ejecución de 82 s pasó a durar 910 s. El script
   `limpia.py` cruza la hora de cada ronda (`runs[].finished` y `total_s`) con los eventos
   506/507 del registro del sistema, y descarta las rondas que se solapan con un intervalo
   de suspensión. El resultado filtrado está en `extremos/limpio/` y es lo que se usa en
   las tablas «con batería». Desde las 16:38, `inter2.py` comprueba los eventos después
   de cada ejecución y la repite si hubo suspensión. No hubo ninguna más: las 267
   ejecuciones de `scripts/ejecuciones.log` están marcadas como `ok`.
2. **Batería y corriente por separado.** Las tablas principales (secciones 1-6) usan los
   datos **con corriente**. Las de batería se conservan en el apéndice A, como
   comprobación de coherencia: los speedups relativos son parecidos porque cada
   comparación se intercala, pero los tiempos absolutos no.
3. **Intercalado.** Todas las comparaciones C++ alternan las variantes con
   `run_bench.py interleave`, mediante el envoltorio `inter2.py`. Este añade las
   variantes virtuales `etiqueta=build:backend=gmp` y `:threads=1`, y detecta las
   suspensiones.
4. **Mediciones focalizadas.** Las técnicas dudosas, o con resultados contradictorios entre
   batería y corriente, se midieron otra vez con 6 rondas y ≥ 5 repeticiones, filtrando a
   sus escenarios (sección 6). **Para esas técnicas, el veredicto sale de la medición
   focalizada.**

### Línea base y compilaciones

- **Referencia** = `baseline-original` (67331a2) con `-O3`. El worktree
  `scratchpad/baseline` existe, pero su árbol de trabajo está vacío (archivos borrados
  sin confirmar), así que no se tocó. Se usó otro worktree limpio en el mismo commit,
  `scratchpad/med-b0` (git worktree add --detach … baseline-original). Su `include/`
  coincide con `git archive baseline-original`.
- **Cada commit** se compiló contra el `include/` de su propio worktree temporal
  (`scratchpad/med-<x>`, git worktree add --detach) con **el mismo `bench_suite.cpp`**: el
  de 8f1d896. Es la única versión que admite `--backend auto|gmp`; compila sin avisos
  contra las cabeceras originales y conserva los nombres y checksums de los 128
  escenarios. Así, todas las variantes ejecutan exactamente el mismo código de medida.
- Las etiquetas `a2`, `a5` y `th0-O3` se volvieron a enlazar con `-DALCP_RELINK=1`
  porque Smart App Control bloqueó el primer ejecutable. Ese define no afecta al código.
- Los perfiles de toolchain (`tb0-*` sobre la referencia, `th0-*` sobre 5f49217) se
  compilaron con el `run_bench.py` de la rama wip 584b9aa (`--fast-alloc` hace
  `-include alcp/alloc.hpp -DALCP_FAST_ALLOC`).
- Módulos: `run_bench.py module -O3` con el `bindings.cpp` y el `include/` de cada
  commit: `modb0` = referencia, `modh0` = 5f49217, `moda6` = 8f1d896 (incluye el binding
  4fb4f1a) y `modm9` = b0a9fe7.
- Capa Python: `git archive <commit> python-bindings/alcp` con **el mismo** `.pyd`
  (`modh0`, de 5f49217, enlazado como vínculo duro porque Smart App Control bloquea las
  copias). Se midió con `bench_alcp.py` de c70bed5.

### Cómo leer las tablas

- **mediana**: mediana de todas las muestras de todas las rondas válidas, en ms. **IQR %**:
  rango intercuartílico entre la mediana. `(!)` indica IQR > 10 %.
- **paso** = mediana(commit k−1) / mediana(commit k) de la misma rama, en el mismo
  intercalado. Un valor > 1 significa que el commit k es más rápido.
- **vs ref** (acumulado) = (ref/5f49217)\_extremos × (5f49217/commit k)\_cadena. Se
  encadena porque la referencia no está en el intercalado de las cadenas, por coste: un
  único paso de la referencia tarda 85-200 s.
- **media geom. paso (afectados)**: media geométrica sobre los escenarios que casan con el
  campo `escenarios_afectados` de la técnica, dentro del subconjunto filtrado de la cadena.
  **(control)**: media sobre los demás escenarios de ese subconjunto. Debería estar en
  ≈1,00x.
- **Criterio del veredicto**: «funciona» si la mejora en los escenarios afectados supera
  claramente el ruido (≥ 5 % de media geométrica, confirmada en las dos condiciones o en la
  medición focalizada) sin regresiones fuera del IQR; «neutra» si queda dentro de ±5 % o
  del IQR; «no funciona» si empeora de forma reproducible; «error» si cambia algún
  checksum.

## Veredictos por técnica

### Rama aritmética (`worktree-wf_e9ad1fab-e9f-1`)

| commit | técnica | paso (afectados) con corriente | con batería | focalizada | acumulado vs ref (afectados) | veredicto |
|---|---|---:|---:|---:|---:|---|
| 489728b | GMP en sitio en Fp y Fq | 2,01x (0,95-3,96) | 1,84x | — | 2,07x | **funciona** |
| 0de8cf9 | Reducción perezosa (LazyAccumulation) en `mul_school`/`divmod` | 4,68x (1,24-13,2) | 4,44x | — | 6,88x | **funciona** |
| dcb905c | Suma y resta de Poly sin relleno ni copias | 1,01x (0,91-1,12) | 1,06x | **1,17x** (Karatsuba 1,06-1,50; gcd 0,91-1,05) | 2,94x | **funciona** en Karatsuba y `gcd_ext`; neutra en `gcd` |
| 65a5f1c (+5c444c0, 4fb4f1a) | Backend Fp64/Fq64 (auto) | 8,69x (0,73-28,1) | 8,68x | — | 34,0x | **funciona**. Única regresión: `fp.m61.elem_pow` 0,73x, que 8f1d896 corrige |
| 65a5f1c con `--backend gmp` | Plantilla FqT<Base> con GMP | 1,00x (0,81-1,18) | 1,01x | — | 4,02x | **neutra**: la plantilla no cuesta nada |
| 2297387 | Fp64: Barrett sin desplazamientos variables | 0,96x (0,80-1,22) | 1,05x | **1,07x** (elem_muladd 1,14-1,22; elem_pow 1,07-1,23; resto 0,98-1,03) | 18,6x | **funciona** solo en operaciones de elementos; en polinomios no se distingue del ruido |
| 8f1d896 | Fp64::pow con intermedios de Montgomery | 1,60x (1,39-1,81) | 1,61x | — | 3,13x (elem_pow) | **funciona**: `fp.m61.elem_pow` pasa de 0,73x a 1,23x frente a la referencia |
| 5c444c0 | bench_suite `--backend` | infraestructura | | | | n/a |
| 4fb4f1a | Binding con std::variant (u64/GMP) | módulo 5,83x frente a 5f49217 | | | 8,97x frente al módulo original | **funciona** (sección 4) |

### Rama multiplicación (`worktree-wf_e9ad1fab-e9f-2`)

| commit | técnica | paso (afectados) con corriente | con batería | focalizada | acumulado vs ref (afectados) | veredicto |
|---|---|---:|---:|---:|---:|---|
| d4580af | FFT: tablas directas, bit-reversal precalculado, mariposa sin temporal | 1,21x (0,85-1,58) | 1,10x | — | 1,22x | **funciona** |
| f68a18b | Cuadrado con una sola transformada | 1,02x (1,00-1,04) | 1,05x | — | 1,15x | **neutra en la suite**: en `pow_mod`/`is_irreducible` domina `divmod`. El desarrollador midió -31 % en a*a aislado; aquí no se volvió a medir |
| 9bfc1e7 | Karatsuba sobre vistas | 0,94x (0,67-1,28; IQR de hasta el 144 %) | 1,11x | **1,18x** (F_p 1,16-1,34; F_q 1,01-1,06) | 1,09x | **funciona**. Las pérdidas en F_q de la etapa de desarrollo y de la cadena con corriente no se reproducen con 6 rondas |
| 1ad8c02 | NTT en uint64_t + Shoup (p < 2^62) | 4,51x (0,87-44,3) | 4,06x | — | 5,21x | **funciona** |
| 70e3453 | NTT multi-limb mpn/Montgomery (p ≥ 2^62) | 2,45x (1,28-3,98) | 2,60x | — | 3,40x | **funciona** |
| 3660189 | Umbrales de despacho Auto | 1,04x (0,57-2,09) | 1,10x | **1,02x** (`fp.p256.gcd_ext` 2,01x; `fp.ntt.factor` 1,58x en la cadena; `fp.ntt.gcd*` 0,94x dentro del IQR) | 2,59x | **funciona** de forma localizada (p256 y ntt.factor); en el resto es neutra |
| 06e46b5 | fast::dft/idft en palabras de máquina | 11,58x (10,3-12,2) | 8,89x | — | 12,18x | **funciona** |
| e8e07e6 | Karatsuba paralelo (std::async) | 1,83x (0,46-4,22) | 2,17x | **1,23x** (`mul_auto` de F_p sin FFT 2,94-3,13x; gcd 0,96-1,05) | 2,78x | **funciona**. Aumenta el tiempo de CPU total |
| e8e07e6 con `ALCP_THREADS=1` | Coste de la infraestructura paralela | 1,15x | 1,02x | **1,01x** (0,98-1,04) | 1,75x | sin coste con un hilo |
| b0a9fe7 | NTT paralela | 1,27x (0,85-2,95) | 1,40x | **1,08x** (`fq.m61e2.mul_fft.d1024` 1,38x; Toeplitz n1024 1,02x) | 4,55x | **funciona**. No hay regresión en las FFT pequeñas |
| b0a9fe7 con `ALCP_THREADS=1` | Coste de la infraestructura paralela | 0,95x (0,75-1,01) | 1,03x | — | 3,47x | neutra, dentro del ruido |

Parte del paralelismo frente al resto (suite completa con corriente): m9 da 1,70x y m9t1
1,40x frente a la referencia. `mul_karatsuba` pasa de 1,19x con un hilo a 3,53x con todos,
y `mul_auto` de 2,38x a 4,92x. El paralelismo **aumenta el tiempo de CPU total**: aquí
solo se mide el tiempo de reloj.

### Toolchain (`worktree-wf_e9ad1fab-e9f-3`, rama wip 584b9aa; el agente no entregó informe)

| perfil | sobre la referencia | sobre 5f49217 | veredicto |
|---|---:|---:|---|
| `-march=x86-64-v3` | 1,02x (0,95-1,08) | 1,01x (0,93-1,13) | **neutra** |
| `-march=native` | 1,02x (0,97-1,08) | 1,02x (0,98-1,17) | **neutra** |
| LTO (`-flto -flto-partition=one`) | 0,99x (0,94-1,04) | 1,02x (0,96-1,19) | **neutra** (una sola unidad de traducción; GMP y libstdc++ sin LTO) |
| asignador `--fast-alloc` (wip) | **1,23x** (1,11-1,37) | **1,25x** (1,15-1,48) | **funciona** en esta medición (checksums iguales), pero el código está sin terminar y sin tests |
| PGO | — | — | no medido: la rama wip no llegó a tener un flujo de PGO utilizable |

Subconjunto de 19 escenarios, 3 rondas con corriente, IQR relativo mediano del 7,7 %.

### Capa Python (`worktree-wf_e9ad1fab-e9f-4`)

| commit | técnica | paso (afectados) | resto `alcp.*` | veredicto |
|---|---|---:|---:|---|
| 8fcc151 | Métodos nativos ligados en `__slots__` | 1,03x (0,87-1,13) | 1,03x | **neutra**: igual que el control |
| abc2f4d | Camino rápido para int en + - * | 1,15x (0,98-1,36) | 0,98x | **funciona** |
| 27873ea | Analizador: x^k como monomio | 1,69x (`parse_long` 5,07x) | 1,02x | **funciona**. Toca el analizador, no un algoritmo del núcleo |
| 29d5628 | Caché de cadenas analizadas | 94x (1,02-1670) | 1,02x | **funciona** con cadenas repetidas; con cadenas nuevas no cambia nada |
| 0ab679a | `_pretty`: cada coeficiente se formatea una vez | 1,47x (1,06-3,12) | 1,01x | **funciona** |
| a3399d9 | Caché de textos de elementos | 5,19x (1,14-23,5) | 1,00x | **funciona** |
| cab2c2a | Pickle de elementos más directo | 1,22x | 0,98x | **funciona** |
| a833bcb | Caminos rápidos int/elemento en la conversión y en f(v) | 0,98x (0,95-1,00) | 1,01x | **neutra** |
| e683a88 | Listas de coeficientes con caminos rápidos | 1,28x (1,26-1,30) | 1,02x | **funciona** |
| 06d922a | Caminos rápidos para // % divmod | 1,04x (1,02-1,06) | 1,01x | **neutra** con grado 100: el beneficio esperado está en grados pequeños, que la batería no cubre |
| 0438b74 | Resultado creado en línea (DUDOSA) | 1,02x (0,97-1,06) | 0,97x | **neutra**: no llega al 3-5 % que pedía su autor. Se recomienda revertirlo |
| c70bed5 | `_pretty`: potencias en línea | 1,12x (1,05-1,21) | 0,99x | **funciona** |

## Qué no funcionó

| técnica | rama / commit | medida | comentario |
|---|---|---|---|
| Tokenizador con regex en el analizador | `exp/analizador-tokenizador-regex` (ff9ad21, sobre 27873ea) | 1,05x en los escenarios del analizador (1,02-1,14) frente a 1,02x en el resto: efecto neto ≈ +3 %, dentro del ruido | Confirma el descarte: el coste está en las operaciones con polinomios de cada término, no en el tokenizador |
| Cuadrado con una sola transformada | f68a18b | 1,02-1,05x en `pow_mod`/`is_irreducible` | Correcto, pero no se nota en la suite (domina `divmod`) |
| Métodos ligados en `__slots__` | 8fcc151 | 1,03x, igual que el control | Sin efecto propio apreciable |
| Caminos rápidos en conversión y f(v) | a833bcb | 0,98x | Sin efecto |
| Caminos rápidos en // % | 06d922a | 1,04x con grado 100 | Sin efecto donde manda el nativo |
| Resultado en línea en + - * / | 0438b74 | 1,02x (control 0,97x) | **Revertir**: duplica código sin ganancia medible |
| `-march=x86-64-v3` / `-march=native` | rama wip toolchain | 1,01-1,02x | La GMP de MSYS2 ya despacha por CPU; el código propio no se vectoriza de forma apreciable |
| LTO | rama wip toolchain (06ebf59) | 0,99-1,02x | Una sola unidad de traducción: no hay nada que ganar |
| Barrett sin desplazamientos variables (efecto global) | 2297387 | 0,96x / 1,05x en la cadena; 1,07x focalizada | Solo mejora `elem_*`. En las operaciones de polinomios queda en el ruido |
| Infraestructura paralela con 1 hilo | e8e07e6 y b0a9fe7 con `ALCP_THREADS=1` | 1,01x y 0,95-1,03x | No empeora, pero sin hilos tampoco aporta |

Ninguna técnica empeora de forma reproducible. Las regresiones puntuales de las cadenas
(p. ej. `fq.m61e2.gcd.d150` 0,46-0,57x o `fp.p101.mul_auto.d2048` 0,67x) tenían IQR del
90-144 % y no se reprodujeron en las mediciones focalizadas (0,99-1,05x).

## Pendientes

- Repetir las cadenas completas **con corriente y plan de alto rendimiento**, con el equipo
  sin suspensiones: la cadena de multiplicación con corriente solo tiene 2 rondas.
- La rama de toolchain está sin terminar. El asignador `--fast-alloc` merece completarse
  y probarse (tests y build.py all); PGO no se pudo medir.
- `f68a18b`: medir con un escenario de cuadrados puros si se quiere cuantificar el -31 %
  aislado que midió el desarrollador.
- No se volvió a medir la duración de pytest en reposo (pendiente de la capa Python): esta
  etapa no ejecutó `build.py all`, porque no tocaba código.
- Recalibrar en reposo los umbrales de paralelismo y de despacho (kParWordMin,
  kParLimbMin, kParKaratsubaMin, prefer_fft).
- Revertir `0438b74` al integrar la capa Python.

## Organización de `bench/results/fase3/`

| directorio | contenido |
|---|---|
| `extremos/` (+ `limpio/`) | Suite completa, 6 variantes (b0, h0, a6, a6gmp, m9, m9t1), 3 rondas con batería. `limpio/` excluye las rondas con suspensión |
| `extremos_ac/`, `extremos_ac_extra/` | Las mismas variantes, 2+1 rondas con corriente (datos principales) |
| `aritmetica/`, `aritmetica_ac/` | Cadena h0→a1→…→a6 (+ a4 con `--backend gmp`), 37 escenarios |
| `multiplicacion/`, `multiplicacion_ac/` | Cadena h0→m1→…→m9 (+ m8/m9 con `ALCP_THREADS=1`), 35 escenarios |
| `focal/` | Mediciones focalizadas de 6 rondas: `barrett`, `sumaresta`, `karatsuba_vistas`, `ntt_par`, `umbrales_par` |
| `toolchain/` | 10 perfiles de flags (5 sobre la referencia y 5 sobre 5f49217), 19 escenarios |
| `modulo/` | `bench_python.py`, 4 módulos × 6 rondas (`py-<módulo>-r<k>.json`) |
| `capa_python/` | `bench_alcp.py --interleave`, q00 = 5f49217 … q12 = c70bed5, q13 = ff9ad21 (exp) |
| `scripts/` | Scripts de la campaña (rutas sustituidas por `<SCRATCH>`/`<REPO>`) y registro de ejecuciones |

Etiquetas: `b0` = baseline-original, `h0` = 5f49217, `a1`…`a6` = 489728b, 0de8cf9, dcb905c,
4fb4f1a, 2297387, 8f1d896; `m1`…`m9` = d4580af, f68a18b, 9bfc1e7, 1ad8c02, 70e3453,
3660189, 06e46b5, e8e07e6, b0a9fe7.

---

# Tablas detalladas (con corriente)

<!-- datos: extremos_ac + extremos_ac_extra (3 rondas con corriente) -->
## 1. Extremos: suite completa (128 escenarios)

| variante | escenarios | checksums = ref | errores | media geom. vs ref | escenarios con IQR > 10 % |
|---|---:|---:|---:|---:|---:|
| `b0` | 128 | 128 | 0 | 1.00x (ref) | 74 |
| `h0` | 128 | 128 | 0 | 1.04x | 61 |
| `a6` | 128 | 128 | 0 | 18.76x | 87 |
| `a6gmp` | 128 | 128 | 0 | 2.84x | 64 |
| `m9` | 128 | 128 | 0 | 1.70x | 74 |
| `m9t1` | 128 | 128 | 0 | 1.40x | 74 |

### Media geométrica del speedup frente a la referencia, por grupo

| grupo | n | h0 | a6 | a6gmp | m9 | m9t1 |
|---|---:|---:|---:|---:|---:|---:|
| `fast.ntt` | 8 | 1.04x | 17.20x | 1.06x | 10.09x | 10.19x |
| `fast.p7` | 2 | 1.00x | 50.58x | 1.72x | 1.67x | 0.99x |
| `fp.m127` | 14 | 1.03x | 2.63x | 2.64x | 1.28x | 1.08x |
| `fp.m61` | 14 | 1.03x | 36.04x | 2.71x | 1.24x | 1.04x |
| `fp.ntt` | 16 | 1.00x | 38.30x | 2.40x | 2.92x | 2.56x |
| `fp.p101` | 14 | 1.04x | 52.13x | 3.09x | 1.32x | 1.06x |
| `fp.p256` | 16 | 1.06x | 2.23x | 2.39x | 1.95x | 1.60x |
| `fp.p7` | 14 | 1.09x | 50.57x | 2.83x | 1.31x | 1.06x |
| `fq.build` | 2 | 1.00x | 44.08x | 3.50x | 0.97x | 0.81x |
| `fq.m61e2` | 10 | 1.02x | 13.66x | 3.72x | 1.37x | 0.98x |
| `fq.q2e8` | 6 | 1.02x | 37.67x | 5.16x | 1.21x | 0.90x |
| `fq.q3e30` | 6 | 1.04x | 50.29x | 5.72x | 1.18x | 1.02x |
| `fq.q7e4` | 6 | 1.03x | 30.43x | 5.33x | 1.18x | 0.95x |

### Por operación

| operación | n | h0 | a6 | a6gmp | m9 | m9t1 |
|---|---:|---:|---:|---:|---:|---:|
| `build` | 2 | 1.00x | 44.08x | 3.50x | 0.97x | 0.81x |
| `dft` | 3 | 1.08x | 12.11x | 1.09x | 10.86x | 10.77x |
| `divmod` | 8 | 1.05x | 60.21x | 10.15x | 0.99x | 0.99x |
| `elem_inv` | 10 | 1.02x | 7.56x | 1.52x | 0.97x | 0.97x |
| `elem_muladd` | 10 | 1.01x | 9.65x | 1.54x | 0.98x | 0.93x |
| `elem_pow` | 6 | 1.06x | 2.75x | 0.96x | 0.96x | 1.04x |
| `factor` | 10 | 1.02x | 17.63x | 3.14x | 1.04x | 1.04x |
| `gcd` | 8 | 1.07x | 12.20x | 1.94x | 1.04x | 1.04x |
| `gcd_ext` | 6 | 1.06x | 9.66x | 1.39x | 2.01x | 2.03x |
| `idft` | 1 | 1.06x | 11.93x | 1.09x | 8.39x | 8.84x |
| `is_irreducible` | 7 | 1.04x | 30.16x | 4.59x | 1.01x | 1.02x |
| `mul_auto` | 7 | 1.04x | 14.39x | 2.08x | 4.92x | 2.38x |
| `mul_fft` | 6 | 0.98x | 5.24x | 1.21x | 7.72x | 5.57x |
| `mul_karatsuba` | 16 | 1.04x | 28.61x | 3.88x | 3.53x | 1.19x |
| `mul_school` | 16 | 1.06x | 63.26x | 9.68x | 0.96x | 1.00x |
| `pow_mod` | 6 | 1.04x | 30.57x | 4.44x | 1.21x | 1.24x |
| `toeplitz_divmod` | 2 | 0.98x | 29.40x | 1.13x | 1.26x | 1.04x |
| `toeplitz_lower_inverse` | 3 | 1.01x | 34.22x | 1.35x | 9.25x | 7.47x |
| `toeplitz_lower_vec` | 1 | 1.01x | 26.70x | 1.06x | 22.33x | 22.65x |

### Coste de robustez (5f49217 frente a la referencia): los 10 escenarios más lentos y los 5 más rápidos

| escenario | ref mediana ms (IQR %) | 5f49217 mediana ms (IQR %) | speedup |
|---|---:|---:|---:|
| `fp.ntt.elem_muladd.n400000` | 45.55 (5) | 53.45 (44) | 0.85x |
| `fp.ntt.mul_fft.d32768` | 1060 (13) | 1162 (23) | 0.91x |
| `fp.ntt.mul_school.d64` | 9.89 (5) | 10.79 (7) | 0.92x |
| `fp.ntt.mul_fft.d4096` | 79.74 (17) | 85.23 (9) | 0.94x |
| `fp.m127.elem_inv.n50000` | 32.12 (13) | 34.31 (10) | 0.94x |
| `fq.q3e30.divmod.d64_32` | 61.57 (16) | 65.39 (21) | 0.94x |
| `fq.m61e2.mul_fft.d1024` | 87.37 (13) | 92.65 (3) | 0.94x |
| `fp.p101.factor.d100` | 86.27 (6) | 91.20 (22) | 0.95x |
| `fp.ntt.elem_pow.n50000` | 14.81 (5) | 15.65 (19) | 0.95x |
| `fq.build.q7e20` | 19.25 (18) | 20.16 (6) | 0.95x |
| `fp.p7.elem_pow.n100000` | 16.19 (29) | 13.79 (12) | 1.17x |
| `fp.p256.gcd_ext.d200` | 149 (24) | 126 (4) | 1.18x |
| `fp.p7.mul_school.d64` | 8.42 (29) | 7.01 (12) | 1.20x |
| `fp.p7.mul_school.d512` | 26.93 (34) | 22.09 (17) | 1.22x |
| `fp.p256.gcd.d300` | 25.60 (29) | 17.88 (11) | 1.43x |

<!-- datos: aritmetica_ac -->
## 2.1. Cadena aritmetica

| variante | escenarios | checksums = ref (b0) | = 5f49217 | errores | IQR > 10 % |
|---|---:|---:|---:|---:|---:|
| `a1` | 37 | 37 | 37 | 0 | 23 |
| `a2` | 37 | 37 | 37 | 0 | 15 |
| `a3` | 37 | 37 | 37 | 0 | 33 |
| `a4` | 37 | 37 | 37 | 0 | 29 |
| `a4gmp` | 37 | 37 | 37 | 0 | 31 |
| `a5` | 37 | 37 | 37 | 0 | 37 |
| `a6` | 37 | 37 | 37 | 0 | 33 |
| `h0` | 37 | 37 | 37 | 0 | 26 |

| commit | técnica | n afect. | media geom. paso (afectados) | rango paso | media geom. paso (control) | media geom. acumulada vs ref (afectados) |
|---|---|---:|---:|---:|---:|---:|
| 489728b | GMP en sitio en Fp y Fq | 18 | 2.01x | 0.95-3.96 | 1.00x | 2.07x |
| 0de8cf9 | Reducción perezosa en mul_school/divmod (LazyAccumulation) | 25 | 4.68x | 1.24-13.24 | 1.09x | 6.88x |
| dcb905c | Suma/resta de Poly sin relleno ni copias | 5 | 1.01x | 0.91-1.12 | 0.95x | 2.94x |
| 65a5f1c (+5c444c0, 4fb4f1a) | Backend Fp64/Fq64 (auto) | 31 | 8.69x | 0.73-28.14 | 0.99x | 34.00x |
| 65a5f1c, --backend gmp | FqT<Base> con backend GMP forzado (coste de la plantilla) | 37 | 1.00x | 0.81-1.18 | — | 4.02x |
| 2297387 | Fp64: Barrett sin desplazamientos variables | 23 | 0.96x | 0.80-1.22 | 0.98x | 18.61x |
| 8f1d896 | Fp64::pow en forma de Montgomery | 3 | 1.60x | 1.39-1.81 | 0.99x | 3.13x |

<details><summary>489728b — GMP en sitio en Fp y Fq: detalle por escenario</summary>

| escenario | h0 mediana ms (IQR %) | a1 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.m127.elem_muladd.n200000` | 28.04 (6) | 26.55 (9) | 1.06x | 1.09x |
| `fp.m127.mul_school.d512` | 34.66 (8) | 34.08 (6) | 1.02x | 1.06x |
| `fp.m127.divmod.d2000_1000` | 136 (30) | 136 (14) (!) | 1.00x | 1.04x |
| `fp.p256.elem_muladd.n200000` | 52.21 (35) | 35.90 (14) (!) | 1.45x | 1.45x |
| `fp.p256.mul_school.d512` | 49.24 (39) | 51.71 (19) (!) | 0.95x | 0.94x |
| `fp.p256.divmod.d1000_500` | 85.22 (54) | 86.93 (16) (!) | 0.98x | 1.00x |
| `fq.q3e30.elem_muladd.n4000` | 245 (20) | 61.88 (12) (!) | 3.96x | 4.10x |
| `fq.q3e30.mul_school.d32` | 62.59 (12) | 17.35 (10) (!) | 3.61x | 3.65x |
| `fq.q3e30.divmod.d64_32` | 64.72 (26) | 18.00 (9) (!) | 3.60x | 3.39x |
| `fq.q2e8.elem_muladd.n50000` | 252 (8) | 83.46 (28) (!) | 3.02x | 3.06x |
| `fq.q2e8.mul_school.d128` | 75.99 (9) | 30.26 (24) (!) | 2.51x | 2.59x |
| `fq.q2e8.factor.d28` | 114 (19) | 43.72 (19) (!) | 2.62x | 2.62x |
| `fq.q7e4.mul_school.d128` | 42.51 (6) | 16.30 (10) | 2.61x | 2.96x |
| `fq.q7e4.gcd.d200` | 119 (17) | 45.24 (6) (!) | 2.62x | 2.77x |
| `fq.m61e2.elem_muladd.n100000` | 98.40 (15) | 40.56 (22) (!) | 2.43x | 2.37x |
| `fq.m61e2.mul_school.d128` | 16.12 (12) | 7.03 (3) (!) | 2.29x | 2.52x |
| `fq.m61e2.mul_karatsuba.d512` | 96.08 (7) | 44.77 (17) (!) | 2.15x | 2.32x |
| `fq.m61e2.divmod.d500_250` | 62.11 (27) | 27.06 (30) (!) | 2.30x | 2.49x |

</details>

<details><summary>0de8cf9 — Reducción perezosa en mul_school/divmod (LazyAccumulation): detalle por escenario</summary>

| escenario | a1 mediana ms (IQR %) | a2 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p101.mul_school.d512` | 33.18 (24) | 2.71 (14) (!) | 12.24x | 13.30x |
| `fp.p101.mul_karatsuba.d512` | 14.13 (10) | 5.32 (1) | 2.65x | 2.83x |
| `fp.p101.divmod.d2000_1000` | 124 (12) | 10.15 (2) (!) | 12.26x | 13.64x |
| `fp.p101.gcd.d500` | 34.13 (30) | 19.75 (4) (!) | 1.73x | 1.66x |
| `fp.p101.gcd_ext.d400` | 66.30 (10) | 42.15 (4) (!) | 1.57x | 1.85x |
| `fp.p101.factor.d100` | 90.22 (7) | 26.60 (4) | 3.39x | 3.57x |
| `fp.ntt.mul_school.d512` | 32.34 (8) | 2.85 (6) | 11.34x | 11.51x |
| `fp.ntt.factor.d40` | 77.16 (13) | 24.30 (38) (!) | 3.18x | 2.99x |
| `fp.m61.mul_school.d512` | 36.86 (41) | 2.90 (7) (!) | 12.69x | 12.16x |
| `fp.m61.divmod.d2000_1000` | 136 (5) | 11.01 (7) | 12.36x | 12.28x |
| `fp.m61.factor.d26` | 58.74 (5) | 19.03 (9) | 3.09x | 3.07x |
| `fp.m127.mul_school.d512` | 34.08 (6) | 4.22 (3) | 8.07x | 8.55x |
| `fp.m127.divmod.d2000_1000` | 136 (14) | 15.69 (3) (!) | 8.66x | 9.02x |
| `fp.p256.mul_school.d512` | 51.71 (19) | 6.83 (12) (!) | 7.57x | 7.12x |
| `fp.p256.divmod.d1000_500` | 86.93 (16) | 6.56 (34) (!) | 13.24x | 13.21x |
| `fq.q3e30.mul_school.d32` | 17.35 (10) | 6.25 (6) (!) | 2.78x | 10.12x |
| `fq.q3e30.divmod.d64_32` | 18.00 (9) | 7.50 (13) (!) | 2.40x | 8.13x |
| `fq.q2e8.mul_school.d128` | 30.26 (24) | 5.39 (2) (!) | 5.62x | 14.56x |
| `fq.q2e8.factor.d28` | 43.72 (19) | 23.95 (8) (!) | 1.83x | 4.78x |
| `fq.q7e4.mul_school.d128` | 16.30 (10) | 2.58 (3) | 6.33x | 18.75x |
| `fq.q7e4.gcd.d200` | 45.24 (6) | 31.05 (27) (!) | 1.46x | 4.04x |
| `fq.m61e2.mul_school.d128` | 7.03 (3) | 0.896 (26) (!) | 7.84x | 19.73x |
| `fq.m61e2.mul_karatsuba.d512` | 44.77 (17) | 17.78 (32) (!) | 2.52x | 5.83x |
| `fq.m61e2.divmod.d500_250` | 27.06 (30) | 3.25 (8) (!) | 8.34x | 20.76x |
| `fast.p7.toeplitz_divmod.d2000_1000` | 85.31 (29) | 68.63 (13) (!) | 1.24x | 1.28x |

</details>

<details><summary>dcb905c — Suma/resta de Poly sin relleno ni copias: detalle por escenario</summary>

| escenario | a2 mediana ms (IQR %) | a3 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p101.mul_karatsuba.d512` | 5.32 (1) | 4.74 (26) (!) | 1.12x | 3.18x |
| `fp.p101.gcd.d500` | 19.75 (4) | 21.28 (27) (!) | 0.93x | 1.54x |
| `fp.p101.gcd_ext.d400` | 42.15 (4) | 41.62 (19) (!) | 1.01x | 1.87x |
| `fq.q7e4.gcd.d200` | 31.05 (27) | 34.25 (32) (!) | 0.91x | 3.66x |
| `fq.m61e2.mul_karatsuba.d512` | 17.78 (32) | 15.82 (38) (!) | 1.12x | 6.56x |

</details>

<details><summary>65a5f1c (+5c444c0, 4fb4f1a) — Backend Fp64/Fq64 (auto): detalle por escenario</summary>

| escenario | a3 mediana ms (IQR %) | a4 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p101.elem_muladd.n400000` | 47.99 (26) | 2.69 (31) (!) | 17.83x | 17.84x |
| `fp.p101.elem_pow.n100000` | 22.56 (37) | 4.45 (28) (!) | 5.07x | 5.14x |
| `fp.p101.mul_school.d512` | 3.12 (53) | 0.220 (24) (!) | 14.21x | 164x |
| `fp.p101.mul_karatsuba.d512` | 4.74 (26) | 0.281 (21) (!) | 16.88x | 53.70x |
| `fp.p101.divmod.d2000_1000` | 11.66 (17) | 0.870 (15) (!) | 13.40x | 159x |
| `fp.p101.gcd.d500` | 21.28 (27) | 1.33 (14) (!) | 15.98x | 24.63x |
| `fp.p101.gcd_ext.d400` | 41.62 (19) | 2.67 (11) (!) | 15.58x | 29.20x |
| `fp.p101.factor.d100` | 27.82 (31) | 1.84 (12) (!) | 15.15x | 51.73x |
| `fp.ntt.elem_pow.n50000` | 15.57 (22) | 12.81 (9) (!) | 1.22x | 1.18x |
| `fp.ntt.mul_school.d512` | 3.07 (37) | 0.206 (14) (!) | 14.91x | 159x |
| `fp.ntt.mul_fft.d4096` | 90.97 (22) | 3.23 (18) (!) | 28.14x | 24.43x |
| `fp.ntt.factor.d40` | 27.19 (12) | 1.55 (21) (!) | 17.52x | 46.82x |
| `fp.m61.elem_muladd.n400000` | 42.30 (19) | 2.89 (28) (!) | 14.62x | 15.79x |
| `fp.m61.elem_inv.n100000` | 21.50 (9) | 10.42 (28) (!) | 2.06x | 2.31x |
| `fp.m61.elem_pow.n20000` | 8.67 (4) | 11.88 (3) | 0.73x | 0.92x |
| `fp.m61.mul_school.d512` | 3.18 (16) | 0.222 (7) (!) | 14.34x | 159x |
| `fp.m61.divmod.d2000_1000` | 10.94 (24) | 0.852 (7) (!) | 12.83x | 159x |
| `fp.m61.factor.d26` | 21.37 (30) | 1.30 (12) (!) | 16.50x | 45.15x |
| `fq.q3e30.elem_muladd.n4000` | 58.89 (13) | 6.59 (6) (!) | 8.93x | 38.52x |
| `fq.q3e30.mul_school.d32` | 6.40 (22) | 0.842 (17) (!) | 7.60x | 75.13x |
| `fq.q3e30.divmod.d64_32` | 7.64 (11) | 1.04 (14) (!) | 7.35x | 58.62x |
| `fq.q2e8.elem_muladd.n50000` | 83.59 (19) | 10.75 (30) (!) | 7.78x | 23.73x |
| `fq.q2e8.mul_school.d128` | 5.44 (35) | 0.675 (19) (!) | 8.07x | 116x |
| `fq.q2e8.factor.d28` | 25.22 (20) | 3.41 (11) (!) | 7.40x | 33.62x |
| `fq.q7e4.mul_school.d128` | 2.73 (24) | 0.325 (4) (!) | 8.40x | 149x |
| `fq.q7e4.gcd.d200` | 34.25 (32) | 4.32 (15) (!) | 7.93x | 29.04x |
| `fq.m61e2.elem_muladd.n100000` | 39.25 (15) | 11.91 (7) (!) | 3.30x | 8.08x |
| `fq.m61e2.mul_school.d128` | 0.928 (97) | 0.148 (22) (!) | 6.26x | 119x |
| `fq.m61e2.mul_karatsuba.d512` | 15.82 (38) | 3.65 (21) (!) | 4.34x | 28.44x |
| `fq.m61e2.divmod.d500_250` | 3.28 (4) | 0.617 (16) (!) | 5.31x | 109x |
| `fast.p7.toeplitz_divmod.d2000_1000` | 62.57 (18) | 2.36 (18) (!) | 26.54x | 37.29x |

</details>

<details><summary>65a5f1c, --backend gmp — FqT<Base> con backend GMP forzado (coste de la plantilla): detalle por escenario</summary>

| escenario | a3 mediana ms (IQR %) | a4gmp mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p101.elem_muladd.n400000` | 47.99 (26) | 45.55 (10) (!) | 1.05x | 1.05x |
| `fp.p101.elem_pow.n100000` | 22.56 (37) | 19.16 (16) (!) | 1.18x | 1.19x |
| `fp.p101.mul_school.d512` | 3.12 (53) | 2.79 (13) (!) | 1.12x | 12.91x |
| `fp.p101.mul_karatsuba.d512` | 4.74 (26) | 4.18 (22) (!) | 1.13x | 3.61x |
| `fp.p101.divmod.d2000_1000` | 11.66 (17) | 10.40 (19) (!) | 1.12x | 13.31x |
| `fp.p101.gcd.d500` | 21.28 (27) | 19.39 (10) (!) | 1.10x | 1.69x |
| `fp.p101.gcd_ext.d400` | 41.62 (19) | 39.40 (12) (!) | 1.06x | 1.98x |
| `fp.p101.factor.d100` | 27.82 (31) | 24.86 (11) (!) | 1.12x | 3.82x |
| `fp.ntt.elem_pow.n50000` | 15.57 (22) | 14.63 (5) (!) | 1.06x | 1.04x |
| `fp.ntt.mul_school.d512` | 3.07 (37) | 2.90 (9) (!) | 1.06x | 11.30x |
| `fp.ntt.mul_fft.d4096` | 90.97 (22) | 83.26 (16) (!) | 1.09x | 0.95x |
| `fp.ntt.factor.d40` | 27.19 (12) | 25.11 (14) (!) | 1.08x | 2.89x |
| `fp.m61.elem_muladd.n400000` | 42.30 (19) | 44.98 (8) (!) | 0.94x | 1.02x |
| `fp.m61.elem_inv.n100000` | 21.50 (9) | 21.72 (11) (!) | 0.99x | 1.11x |
| `fp.m61.elem_pow.n20000` | 8.67 (4) | 8.52 (4) | 1.02x | 1.28x |
| `fp.m61.mul_school.d512` | 3.18 (16) | 2.91 (12) (!) | 1.09x | 12.13x |
| `fp.m61.divmod.d2000_1000` | 10.94 (24) | 11.31 (21) (!) | 0.97x | 11.96x |
| `fp.m61.factor.d26` | 21.37 (30) | 21.64 (20) (!) | 0.99x | 2.70x |
| `fp.m127.elem_muladd.n200000` | 26.13 (29) | 26.60 (44) (!) | 0.98x | 1.09x |
| `fp.m127.mul_school.d512` | 4.52 (32) | 4.55 (152) (!) | 0.99x | 7.94x |
| `fp.m127.divmod.d2000_1000` | 16.78 (14) | 16.91 (83) (!) | 0.99x | 8.37x |
| `fp.p256.elem_muladd.n200000` | 35.24 (5) | 37.86 (204) (!) | 0.93x | 1.38x |
| `fp.p256.mul_school.d512` | 7.47 (16) | 6.97 (264) (!) | 1.07x | 6.98x |
| `fp.p256.divmod.d1000_500` | 6.56 (21) | 7.12 (353) (!) | 0.92x | 12.17x |
| `fq.q3e30.elem_muladd.n4000` | 58.89 (13) | 64.16 (50) (!) | 0.92x | 3.96x |
| `fq.q3e30.mul_school.d32` | 6.40 (22) | 6.47 (8) (!) | 0.99x | 9.78x |
| `fq.q3e30.divmod.d64_32` | 7.64 (11) | 8.01 (13) (!) | 0.95x | 7.61x |
| `fq.q2e8.elem_muladd.n50000` | 83.59 (19) | 91.76 (10) (!) | 0.91x | 2.78x |
| `fq.q2e8.mul_school.d128` | 5.44 (35) | 6.41 (24) (!) | 0.85x | 12.23x |
| `fq.q2e8.factor.d28` | 25.22 (20) | 27.20 (54) (!) | 0.93x | 4.21x |
| `fq.q7e4.mul_school.d128` | 2.73 (24) | 2.70 (47) (!) | 1.01x | 17.86x |
| `fq.q7e4.gcd.d200` | 34.25 (32) | 38.41 (26) (!) | 0.89x | 3.26x |
| `fq.m61e2.elem_muladd.n100000` | 39.25 (15) | 46.62 (17) (!) | 0.84x | 2.06x |
| `fq.m61e2.mul_school.d128` | 0.928 (97) | 0.941 (8) (!) | 0.99x | 18.79x |
| `fq.m61e2.mul_karatsuba.d512` | 15.82 (38) | 14.80 (17) (!) | 1.07x | 7.01x |
| `fq.m61e2.divmod.d500_250` | 3.28 (4) | 4.05 (42) (!) | 0.81x | 16.63x |
| `fast.p7.toeplitz_divmod.d2000_1000` | 62.57 (18) | 66.73 (17) (!) | 0.94x | 1.32x |

</details>

<details><summary>2297387 — Fp64: Barrett sin desplazamientos variables: detalle por escenario</summary>

| escenario | a4 mediana ms (IQR %) | a5 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p101.elem_muladd.n400000` | 2.69 (31) | 2.21 (14) (!) | 1.22x | 21.70x |
| `fp.p101.elem_pow.n100000` | 4.45 (28) | 3.86 (48) (!) | 1.15x | 5.92x |
| `fp.p101.mul_school.d512` | 0.220 (24) | 0.233 (47) (!) | 0.94x | 155x |
| `fp.p101.mul_karatsuba.d512` | 0.281 (21) | 0.307 (43) (!) | 0.91x | 49.11x |
| `fp.p101.divmod.d2000_1000` | 0.870 (15) | 0.918 (35) (!) | 0.95x | 151x |
| `fp.p101.gcd.d500` | 1.33 (14) | 1.47 (22) (!) | 0.91x | 22.31x |
| `fp.p101.gcd_ext.d400` | 2.67 (11) | 3.28 (21) (!) | 0.82x | 23.81x |
| `fp.p101.factor.d100` | 1.84 (12) | 2.27 (50) (!) | 0.81x | 41.91x |
| `fp.ntt.elem_pow.n50000` | 12.81 (9) | 11.09 (87) (!) | 1.16x | 1.37x |
| `fp.ntt.mul_school.d512` | 0.206 (14) | 0.204 (283) (!) | 1.01x | 161x |
| `fp.ntt.mul_fft.d4096` | 3.23 (18) | 3.30 (574) (!) | 0.98x | 23.94x |
| `fp.ntt.factor.d40` | 1.55 (21) | 1.84 (224) (!) | 0.84x | 39.55x |
| `fp.m61.elem_muladd.n400000` | 2.89 (28) | 2.89 (134) (!) | 1.00x | 15.79x |
| `fp.m61.elem_inv.n100000` | 10.42 (28) | 9.80 (83) (!) | 1.06x | 2.45x |
| `fp.m61.elem_pow.n20000` | 11.88 (3) | 11.84 (19) (!) | 1.00x | 0.92x |
| `fp.m61.mul_school.d512` | 0.222 (7) | 0.278 (21) (!) | 0.80x | 127x |
| `fp.m61.divmod.d2000_1000` | 0.852 (7) | 1.01 (21) (!) | 0.84x | 134x |
| `fp.m61.factor.d26` | 1.30 (12) | 1.52 (23) (!) | 0.85x | 38.58x |
| `fp.m127.elem_muladd.n200000` | 26.99 (29) | 30.71 (21) (!) | 0.88x | 0.94x |
| `fp.p256.elem_muladd.n200000` | 36.28 (7) | 37.63 (73) (!) | 0.96x | 1.39x |
| `fq.q3e30.elem_muladd.n4000` | 6.59 (6) | 6.34 (56) (!) | 1.04x | 40.06x |
| `fq.q2e8.elem_muladd.n50000` | 10.75 (30) | 10.75 (88) (!) | 1.00x | 23.72x |
| `fq.m61e2.elem_muladd.n100000` | 11.91 (7) | 11.93 (28) (!) | 1.00x | 8.06x |

</details>

<details><summary>8f1d896 — Fp64::pow en forma de Montgomery: detalle por escenario</summary>

| escenario | a5 mediana ms (IQR %) | a6 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p101.elem_pow.n100000` | 3.86 (48) | 2.38 (21) (!) | 1.62x | 9.61x |
| `fp.ntt.elem_pow.n50000` | 11.09 (87) | 7.97 (29) (!) | 1.39x | 1.90x |
| `fp.m61.elem_pow.n20000` | 11.84 (19) | 6.53 (18) (!) | 1.81x | 1.67x |

</details>

<!-- datos: multiplicacion_ac -->
## 2.2. Cadena multiplicacion

| variante | escenarios | checksums = ref (b0) | = 5f49217 | errores | IQR > 10 % |
|---|---:|---:|---:|---:|---:|
| `h0` | 35 | 35 | 35 | 0 | 25 |
| `m1` | 35 | 35 | 35 | 0 | 18 |
| `m2` | 35 | 35 | 35 | 0 | 12 |
| `m3` | 35 | 35 | 35 | 0 | 24 |
| `m4` | 35 | 35 | 35 | 0 | 17 |
| `m5` | 35 | 35 | 35 | 0 | 17 |
| `m6` | 35 | 35 | 35 | 0 | 20 |
| `m7` | 35 | 35 | 35 | 0 | 18 |
| `m8` | 35 | 35 | 35 | 0 | 27 |
| `m8t1` | 35 | 35 | 35 | 0 | 5 |
| `m9` | 35 | 35 | 35 | 0 | 27 |
| `m9t1` | 35 | 35 | 35 | 0 | 22 |

| commit | técnica | n afect. | media geom. paso (afectados) | rango paso | media geom. paso (control) | media geom. acumulada vs ref (afectados) |
|---|---|---:|---:|---:|---:|---:|
| d4580af | FFT: tablas directas, bit-reversal precalculado, mariposa sin temporal | 17 | 1.21x | 0.85-1.58 | 1.00x | 1.22x |
| f68a18b | Cuadrado con una sola transformada | 2 | 1.02x | 1.00-1.04 | 1.06x | 1.15x |
| 9bfc1e7 | Karatsuba sobre vistas | 10 | 0.94x | 0.67-1.28 | 0.94x | 1.09x |
| 1ad8c02 | NTT en uint64_t + Shoup (p < 2^62) | 15 | 4.51x | 0.87-44.34 | 1.02x | 5.21x |
| 70e3453 | NTT multi-limb mpn/Montgomery (p >= 2^62) | 4 | 2.45x | 1.28-3.98 | 1.02x | 3.40x |
| 3660189 | Umbrales de despacho Auto | 11 | 1.04x | 0.57-2.09 | 0.94x | 2.59x |
| 06e46b5 | fast::dft/idft en palabras de máquina | 4 | 11.58x | 10.28-12.21 | 0.98x | 12.18x |
| e8e07e6 | Karatsuba paralelo (std::async) | 13 | 1.83x | 0.46-4.22 | 0.86x | 2.78x |
| b0a9fe7 | NTT paralela | 12 | 1.27x | 0.85-2.95 | 1.10x | 4.55x |
| e8e07e6 con ALCP_THREADS=1 | Karatsuba paralelo limitado a 1 hilo (sobrecoste) | 13 | 1.15x | 0.99-1.55 | 1.06x | 1.75x |
| b0a9fe7 con ALCP_THREADS=1 | NTT paralela limitada a 1 hilo (sobrecoste) | 12 | 0.95x | 0.75-1.01 | 0.95x | 3.47x |

<details><summary>d4580af — FFT: tablas directas, bit-reversal precalculado, mariposa sin temporal: detalle por escenario</summary>

| escenario | h0 mediana ms (IQR %) | m1 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p101.mul_auto.d2048` | 126 (9) | 149 (39) (!) | 0.85x | 0.89x |
| `fp.ntt.mul_fft.d4096` | 84.95 (22) | 83.78 (69) (!) | 1.01x | 0.95x |
| `fp.ntt.mul_fft.d32768` | 1267 (30) | 947 (31) (!) | 1.34x | 1.22x |
| `fp.ntt.mul_auto.d2048` | 40.08 (29) | 34.47 (67) (!) | 1.16x | 1.24x |
| `fp.p256.mul_fft.d4096` | 115 (11) | 85.30 (18) (!) | 1.35x | 1.38x |
| `fp.p256.mul_auto.d2048` | 49.81 (6) | 36.93 (5) | 1.35x | 1.38x |
| `fq.m61e2.mul_fft.d1024` | 88.46 (8) | 77.93 (1) | 1.14x | 1.07x |
| `fq.m61e2.mul_fft.d4096` | 487 (22) | 400 (8) (!) | 1.22x | 1.16x |
| `fq.m61e2.mul_auto.d2048` | 233 (30) | 181 (29) (!) | 1.29x | 1.29x |
| `fast.ntt.dft.k10` | 17.80 (30) | 16.04 (18) (!) | 1.11x | 1.14x |
| `fast.ntt.dft.k11` | 20.05 (33) | 17.58 (2) (!) | 1.14x | 1.28x |
| `fast.ntt.dft.k12` | 21.21 (33) | 19.45 (3) (!) | 1.09x | 1.21x |
| `fast.ntt.idft.k12` | 23.65 (32) | 19.53 (4) (!) | 1.21x | 1.28x |
| `fast.ntt.toeplitz_lower_inverse.n1024` | 31.92 (35) | 22.40 (1) (!) | 1.42x | 1.45x |
| `fast.ntt.toeplitz_lower_inverse.n8192` | 388 (63) | 246 (7) (!) | 1.58x | 1.55x |
| `fast.ntt.toeplitz_lower_vec.n4096` | 42.15 (16) | 30.28 (2) (!) | 1.39x | 1.40x |
| `fast.ntt.toeplitz_divmod.d2000_1000` | 69.63 (35) | 61.94 (3) (!) | 1.12x | 1.10x |

</details>

<details><summary>f68a18b — Cuadrado con una sola transformada: detalle por escenario</summary>

| escenario | m1 mediana ms (IQR %) | m2 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.ntt.pow_mod.d128_k5` | 502 (7) | 484 (22) (!) | 1.04x | 1.08x |
| `fp.p256.pow_mod.d64_k1` | 401 (5) | 399 (7) | 1.00x | 1.21x |

</details>

<details><summary>9bfc1e7 — Karatsuba sobre vistas: detalle por escenario</summary>

| escenario | m2 mediana ms (IQR %) | m3 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p101.mul_karatsuba.d512` | 13.83 (7) | 10.77 (2) | 1.28x | 1.29x |
| `fp.p101.mul_karatsuba.d4096` | 361 (4) | 319 (9) | 1.13x | 1.32x |
| `fp.p101.mul_auto.d2048` | 118 (1) | 175 (144) (!) | 0.67x | 0.76x |
| `fp.ntt.mul_karatsuba.d512` | 13.45 (9) | 12.55 (17) (!) | 1.07x | 1.20x |
| `fp.ntt.mul_auto.d2048` | 28.86 (6) | 40.48 (32) (!) | 0.71x | 1.05x |
| `fp.p256.mul_auto.d2048` | 38.54 (33) | 42.27 (20) (!) | 0.91x | 1.20x |
| `fq.q3e30.mul_karatsuba.d128` | 464 (10) | 551 (112) (!) | 0.84x | 0.92x |
| `fq.q2e8.mul_karatsuba.d512` | 398 (16) | 381 (29) (!) | 1.04x | 1.00x |
| `fq.m61e2.mul_karatsuba.d512` | 84.81 (7) | 83.03 (32) (!) | 1.02x | 1.13x |
| `fq.m61e2.mul_auto.d2048` | 173 (5) | 200 (20) (!) | 0.87x | 1.17x |

</details>

<details><summary>1ad8c02 — NTT en uint64_t + Shoup (p < 2^62): detalle por escenario</summary>

| escenario | m3 mediana ms (IQR %) | m4 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.ntt.mul_karatsuba.d512` | 12.55 (17) | 10.92 (9) (!) | 1.15x | 1.38x |
| `fp.ntt.mul_fft.d4096` | 89.74 (115) | 2.33 (6) (!) | 38.53x | 34.12x |
| `fp.ntt.mul_fft.d32768` | 1060 (60) | 23.91 (2) (!) | 44.34x | 48.35x |
| `fp.ntt.mul_auto.d2048` | 40.48 (32) | 1.000 (1) (!) | 40.49x | 42.61x |
| `fp.ntt.gcd_ext.d400` | 451 (14) | 48.81 (3) (!) | 9.25x | 11.64x |
| `fp.ntt.pow_mod.d128_k5` | 469 (19) | 405 (21) (!) | 1.16x | 1.29x |
| `fp.ntt.factor.d40` | 87.03 (71) | 65.91 (37) (!) | 1.32x | 1.09x |
| `fast.ntt.dft.k10` | 16.83 (4) | 15.95 (3) | 1.06x | 1.14x |
| `fast.ntt.dft.k11` | 18.21 (3) | 20.81 (72) (!) | 0.87x | 1.08x |
| `fast.ntt.dft.k12` | 19.55 (2) | 18.61 (44) (!) | 1.05x | 1.27x |
| `fast.ntt.idft.k12` | 20.23 (5) | 19.35 (5) | 1.05x | 1.30x |
| `fast.ntt.toeplitz_lower_inverse.n1024` | 23.22 (5) | 1.51 (1) | 15.33x | 21.45x |
| `fast.ntt.toeplitz_lower_inverse.n8192` | 263 (14) | 13.87 (6) (!) | 18.96x | 27.62x |
| `fast.ntt.toeplitz_lower_vec.n4096` | 29.91 (12) | 1.71 (3) (!) | 17.46x | 24.80x |
| `fast.ntt.toeplitz_divmod.d2000_1000` | 63.18 (5) | 56.58 (5) | 1.12x | 1.21x |

</details>

<details><summary>70e3453 — NTT multi-limb mpn/Montgomery (p >= 2^62): detalle por escenario</summary>

| escenario | m4 mediana ms (IQR %) | m5 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p256.mul_fft.d4096` | 94.04 (25) | 23.61 (22) (!) | 3.98x | 4.99x |
| `fp.p256.mul_auto.d2048` | 37.72 (7) | 11.62 (22) (!) | 3.25x | 4.37x |
| `fp.p256.gcd_ext.d200` | 106 (2) | 48.77 (21) (!) | 2.18x | 3.70x |
| `fp.p256.pow_mod.d64_k1` | 375 (9) | 293 (4) | 1.28x | 1.65x |

</details>

<details><summary>3660189 — Umbrales de despacho Auto: detalle por escenario</summary>

| escenario | m5 mediana ms (IQR %) | m6 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p101.mul_auto.d2048` | 98.68 (6) | 102 (20) (!) | 0.97x | 1.30x |
| `fp.p101.gcd_ext.d400` | 62.42 (9) | 66.36 (4) | 0.94x | 1.08x |
| `fp.ntt.mul_auto.d2048` | 1.03 (5) | 1.02 (4) | 1.00x | 41.69x |
| `fp.ntt.gcd_ext.d400` | 47.72 (1) | 52.94 (32) (!) | 0.90x | 10.73x |
| `fp.ntt.factor.d40` | 61.90 (3) | 39.27 (14) (!) | 1.58x | 1.82x |
| `fp.p256.mul_auto.d2048` | 11.62 (22) | 10.41 (2) (!) | 1.12x | 4.88x |
| `fp.p256.gcd_ext.d200` | 48.77 (21) | 23.35 (1) (!) | 2.09x | 7.73x |
| `fq.q7e4.gcd.d200` | 110 (11) | 116 (23) (!) | 0.94x | 1.04x |
| `fq.q7e4.factor.d20` | 60.32 (42) | 62.88 (3) (!) | 0.96x | 1.00x |
| `fq.m61e2.mul_auto.d2048` | 175 (14) | 183 (18) (!) | 0.96x | 1.27x |
| `fq.m61e2.gcd.d150` | 40.73 (3) | 71.45 (92) (!) | 0.57x | 0.60x |

</details>

<details><summary>06e46b5 — fast::dft/idft en palabras de máquina: detalle por escenario</summary>

| escenario | m6 mediana ms (IQR %) | m7 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fast.ntt.dft.k10` | 17.79 (82) | 1.73 (6) (!) | 10.28x | 10.53x |
| `fast.ntt.dft.k11` | 19.75 (18) | 1.67 (5) (!) | 11.83x | 13.44x |
| `fast.ntt.dft.k12` | 20.99 (19) | 1.72 (11) (!) | 12.21x | 13.71x |
| `fast.ntt.idft.k12` | 26.69 (26) | 2.21 (5) (!) | 12.08x | 11.36x |

</details>

<details><summary>e8e07e6 — Karatsuba paralelo (std::async): detalle por escenario</summary>

| escenario | m7 mediana ms (IQR %) | m8 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p101.mul_karatsuba.d512` | 12.07 (9) | 3.49 (16) (!) | 3.46x | 3.99x |
| `fp.p101.mul_karatsuba.d4096` | 337 (7) | 96.70 (9) | 3.48x | 4.35x |
| `fp.p101.mul_auto.d2048` | 109 (14) | 35.78 (9) (!) | 3.05x | 3.72x |
| `fp.ntt.mul_karatsuba.d512` | 11.13 (4) | 3.70 (14) (!) | 3.01x | 4.06x |
| `fp.ntt.mul_auto.d2048` | 1.07 (2) | 1.32 (63) (!) | 0.81x | 32.23x |
| `fp.p256.mul_auto.d2048` | 10.83 (7) | 11.01 (10) | 0.98x | 4.62x |
| `fq.q3e30.mul_karatsuba.d128` | 665 (74) | 165 (6) (!) | 4.04x | 3.09x |
| `fq.q2e8.mul_karatsuba.d512` | 544 (73) | 129 (18) (!) | 4.22x | 2.95x |
| `fq.q7e4.gcd.d200` | 153 (124) | 126 (27) (!) | 1.21x | 0.96x |
| `fq.q7e4.factor.d20` | 59.98 (12) | 64.22 (16) (!) | 0.93x | 0.98x |
| `fq.m61e2.mul_karatsuba.d512` | 76.89 (26) | 29.61 (14) (!) | 2.60x | 3.17x |
| `fq.m61e2.mul_auto.d2048` | 263 (63) | 202 (44) (!) | 1.31x | 1.16x |
| `fq.m61e2.gcd.d150` | 41.91 (19) | 90.90 (76) (!) | 0.46x | 0.47x |

</details>

<details><summary>b0a9fe7 — NTT paralela: detalle por escenario</summary>

| escenario | m8 mediana ms (IQR %) | m9 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p101.mul_auto.d2048` | 35.78 (9) | 41.10 (66) (!) | 0.87x | 3.24x |
| `fp.ntt.mul_fft.d4096` | 2.73 (4) | 2.56 (23) (!) | 1.06x | 30.99x |
| `fp.ntt.mul_fft.d32768` | 23.73 (10) | 21.97 (54) (!) | 1.08x | 52.62x |
| `fp.ntt.mul_auto.d2048` | 1.32 (63) | 1.14 (25) (!) | 1.16x | 37.46x |
| `fp.p256.mul_fft.d4096` | 24.42 (6) | 28.70 (72) (!) | 0.85x | 4.11x |
| `fp.p256.mul_auto.d2048` | 11.01 (10) | 12.82 (62) (!) | 0.86x | 3.96x |
| `fq.m61e2.mul_karatsuba.d512` | 29.61 (14) | 31.91 (25) (!) | 0.93x | 2.94x |
| `fq.m61e2.mul_fft.d1024` | 88.36 (31) | 67.79 (22) (!) | 1.30x | 1.23x |
| `fq.m61e2.mul_fft.d4096` | 556 (90) | 189 (60) (!) | 2.95x | 2.47x |
| `fq.m61e2.mul_auto.d2048` | 202 (44) | 102 (10) (!) | 1.98x | 2.30x |
| `fq.m61e2.divmod.d500_250` | 77.08 (44) | 58.64 (16) (!) | 1.31x | 1.18x |
| `fq.m61e2.gcd.d150` | 90.90 (76) | 42.43 (1) (!) | 2.14x | 1.01x |

</details>

<details><summary>e8e07e6 con ALCP_THREADS=1 — Karatsuba paralelo limitado a 1 hilo (sobrecoste): detalle por escenario</summary>

| escenario | m7 mediana ms (IQR %) | m8t1 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p101.mul_karatsuba.d512` | 12.07 (9) | 11.34 (6) | 1.06x | 1.23x |
| `fp.p101.mul_karatsuba.d4096` | 337 (7) | 305 (5) | 1.10x | 1.38x |
| `fp.p101.mul_auto.d2048` | 109 (14) | 99.76 (26) (!) | 1.09x | 1.33x |
| `fp.ntt.mul_karatsuba.d512` | 11.13 (4) | 10.77 (5) | 1.03x | 1.40x |
| `fp.ntt.mul_auto.d2048` | 1.07 (2) | 1.07 (12) (!) | 1.00x | 39.68x |
| `fp.p256.mul_auto.d2048` | 10.83 (7) | 10.24 (2) | 1.06x | 4.96x |
| `fq.q3e30.mul_karatsuba.d128` | 665 (74) | 433 (5) (!) | 1.53x | 1.17x |
| `fq.q2e8.mul_karatsuba.d512` | 544 (73) | 352 (10) (!) | 1.55x | 1.08x |
| `fq.q7e4.gcd.d200` | 153 (124) | 110 (6) (!) | 1.38x | 1.10x |
| `fq.q7e4.factor.d20` | 59.98 (12) | 60.70 (8) (!) | 0.99x | 1.04x |
| `fq.m61e2.mul_karatsuba.d512` | 76.89 (26) | 77.23 (1) (!) | 1.00x | 1.22x |
| `fq.m61e2.mul_auto.d2048` | 263 (63) | 184 (18) (!) | 1.43x | 1.27x |
| `fq.m61e2.gcd.d150` | 41.91 (19) | 42.53 (5) (!) | 0.99x | 1.01x |

</details>

<details><summary>b0a9fe7 con ALCP_THREADS=1 — NTT paralela limitada a 1 hilo (sobrecoste): detalle por escenario</summary>

| escenario | m8t1 mediana ms (IQR %) | m9t1 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p101.mul_auto.d2048` | 99.76 (26) | 133 (62) (!) | 0.75x | 1.00x |
| `fp.ntt.mul_fft.d4096` | 2.24 (4) | 2.37 (4) | 0.95x | 33.50x |
| `fp.ntt.mul_fft.d32768` | 24.93 (5) | 24.88 (6) | 1.00x | 46.45x |
| `fp.ntt.mul_auto.d2048` | 1.07 (12) | 1.09 (18) (!) | 0.98x | 39.05x |
| `fp.p256.mul_fft.d4096` | 22.51 (3) | 22.61 (4) | 1.00x | 5.22x |
| `fp.p256.mul_auto.d2048` | 10.24 (2) | 10.43 (3) | 0.98x | 4.87x |
| `fq.m61e2.mul_karatsuba.d512` | 77.23 (1) | 79.43 (11) (!) | 0.97x | 1.18x |
| `fq.m61e2.mul_fft.d1024` | 78.60 (3) | 82.08 (14) (!) | 0.96x | 1.02x |
| `fq.m61e2.mul_fft.d4096` | 390 (5) | 413 (17) (!) | 0.94x | 1.13x |
| `fq.m61e2.mul_auto.d2048` | 184 (18) | 182 (11) (!) | 1.01x | 1.28x |
| `fq.m61e2.divmod.d500_250` | 58.08 (11) | 59.79 (12) (!) | 0.97x | 1.16x |
| `fq.m61e2.gcd.d150` | 42.53 (5) | 44.16 (33) (!) | 0.96x | 0.97x |

</details>

## 3. Toolchain (rama wip, perfiles de flags)

| cabeceras | perfil | media geom. vs -O3 | rango | checksums = ref | IQR > 10 % |
|---|---|---:|---:|---:|---:|
| baseline-original | v3 | 1.02x | 0.95-1.08 | 19/19 | 5 |
| baseline-original | nat | 1.02x | 0.97-1.08 | 19/19 | 5 |
| baseline-original | lto | 0.99x | 0.94-1.04 | 19/19 | 7 |
| baseline-original | alloc | 1.23x | 1.11-1.37 | 19/19 | 9 |
| 5f49217 | v3 | 1.01x | 0.93-1.13 | 19/19 | 8 |
| 5f49217 | nat | 1.02x | 0.98-1.17 | 19/19 | 5 |
| 5f49217 | lto | 1.02x | 0.96-1.19 | 19/19 | 7 |
| 5f49217 | alloc | 1.25x | 1.15-1.48 | 19/19 | 6 |

## 4. Módulo nativo visto desde Python (bench_python.py)

| escenario | modb0 mediana ms (IQR %) | modh0 mediana ms (IQR %) | moda6 mediana ms (IQR %) | modm9 mediana ms (IQR %) | 5f49217/ref | a6/ref | a6/5f49217 | m9/5f49217 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `py.ntt.poly_from_list.d1000` | 6.38 (24) | 2.14 (12) | 1.02 (14) | 1.89 (12) | 2.98x | 6.24x | 2.09x | 1.13x |
| `py.ntt.coeffs.d1000` | 7.15 (51) | 2.23 (18) | 2.32 (13) | 2.07 (13) | 3.21x | 3.08x | 0.96x | 1.08x |
| `py.ntt.roundtrip.d1000` | 10.13 (23) | 3.25 (12) | 2.22 (14) | 3.00 (13) | 3.12x | 4.57x | 1.46x | 1.08x |
| `py.ntt.mul_auto.d1000` | 41.84 (57) | 36.92 (12) | 1.37 (14) | 1.20 (12) | 1.13x | 30.48x | 26.89x | 30.77x |
| `py.ntt.mul_fft.d1000` | 38.50 (81) | 36.34 (16) | 1.39 (14) | 1.21 (17) | 1.06x | 27.67x | 26.12x | 29.91x |
| `py.ntt.mul_karatsuba.d1000` | 83.37 (71) | 76.59 (12) | 1.02 (15) | 18.83 (41) | 1.09x | 81.71x | 75.06x | 4.07x |
| `py.ntt.mul_school.d1000` | 120 (14) | 118 (16) | 0.798 (16) | 120 (11) | 1.02x | 150x | 148x | 0.98x |
| `py.p101.factor.d40` | 11.40 (41) | 11.05 (37) | 0.383 (31) | 10.87 (33) | 1.03x | 29.76x | 28.86x | 1.02x |
| `py.ntt.fp_mul_loop.n20000` | 22.06 (10) | 10.13 (18) | 6.61 (14) | 10.88 (21) | 2.18x | 3.34x | 1.53x | 0.93x |
| `py.m521.fp_mul_loop.n5000` | 30.22 (20) | 11.26 (16) | 12.01 (9) | 12.34 (15) | 2.68x | 2.52x | 0.94x | 0.91x |
| `py.m521.poly_from_list.d200` | 9.49 (9) | 2.41 (12) | 2.47 (11) | 2.50 (16) | 3.93x | 3.84x | 0.98x | 0.96x |
| `py.m521.coeffs.d200` | 6.18 (9) | 3.91 (23) | 4.30 (14) | 4.34 (25) | 1.58x | 1.44x | 0.91x | 0.90x |
| `py.m521.mul_auto.d200` | 6.38 (12) | 6.17 (31) | 2.16 (13) | 5.41 (8) | 1.03x | 2.95x | 2.85x | 1.14x |
| `py.is_prime.n2000` | 12.23 (18) | 10.13 (16) | 10.60 (15) | 10.66 (5) | 1.21x | 1.15x | 0.96x | 0.95x |
| `py.fq9.poly_from_list.d500` | 7.70 (28) | 5.18 (19) | 4.49 (8) | 5.24 (10) | 1.49x | 1.72x | 1.15x | 0.99x |
| `py.fq9.mul_auto.d500` | 63.36 (14) | 67.01 (29) | 2.83 (21) | 16.40 (37) | 0.95x | 22.41x | 23.70x | 4.09x |
| `py.fq.q3e30.factor.d5` | 336 (33) | 309 (26) | 10.00 (79) | 346 (31) | 1.09x | 33.62x | 30.92x | 0.89x |
| `py.ntt.dft.k10` | 13.27 (8) | 10.80 (15) | 1.40 (19) | 1.87 (19) | 1.23x | 9.49x | 7.73x | 5.78x |
| `py.ntt.toeplitz_lower_inverse.n1024` | 28.30 (8) | 27.53 (14) | 1.16 (14) | 1.63 (14) | 1.03x | 24.31x | 23.64x | 16.88x |

| variante | media geom. vs módulo original | checksums distintos | errores |
|---|---:|---:|---:|
| `py-modh0` | 1.54x | 0/19 | 0 |
| `py-moda6` | 8.97x | 0/19 | 0 |
| `py-modm9` | 3.24x | 0/19 | 0 |
| `py-moda6` frente a `py-modh0` | 5.83x | | |

## 5. Capa Python (bench_alcp.py, mismo .pyd de 5f49217)

| commit | técnica | n afect. | media geom. paso (afectados) | rango | media geom. resto alcp.* |
|---|---|---:|---:|---:|---:|
| 8fcc151 | métodos nativos ligados en __slots__ | 12 | 1.03x | 0.87-1.13 | 1.03x |
| abc2f4d | camino rápido int en + - * | 3 | 1.15x | 0.98-1.36 | 0.98x |
| 27873ea | analizador: x^k como monomio | 5 | 1.69x | 0.96-5.07 | 1.02x |
| 29d5628 | caché de cadenas analizadas | 4 | 94.40x | 1.02-1670.05 | 1.02x |
| 0ab679a | _pretty: cada coeficiente formateado una vez | 4 | 1.47x | 1.06-3.12 | 1.01x |
| a3399d9 | caché de textos de elementos | 2 | 5.19x | 1.14-23.54 | 1.00x |
| cab2c2a | pickle de elementos más directo | 1 | 1.22x | 1.22-1.22 | 0.98x |
| a833bcb | caminos rápidos int/elemento en conversión y f(v) | 5 | 0.98x | 0.95-1.00 | 1.01x |
| e683a88 | listas de coeficientes con caminos rápidos | 2 | 1.28x | 1.26-1.30 | 1.02x |
| 06d922a | caminos rápidos // % divmod | 2 | 1.04x | 1.02-1.06 | 1.01x |
| 0438b74 | resultado en línea en + - * / (DUDOSA) | 12 | 1.02x | 0.97-1.06 | 0.97x |
| c70bed5 | _pretty: potencias en línea | 4 | 1.12x | 1.05-1.21 | 0.99x |
| ff9ad21 (exp/) | EXPERIMENTO tokenizador regex | 5 | 1.05x | 1.02-1.14 | 1.02x |

Checksums: 0 discrepancias frente a 5f49217 y 0 errores en 14 variantes.

### Antes y después (5f49217 frente a c70bed5): tiempo y cociente paquete/nativo

| escenario | 5f49217 ms | c70bed5 ms | speedup | cociente alcp/nat antes | después |
|---|---:|---:|---:|---:|---:|
| `alcp.fp.add_loop.n20000` | 14.91 | 13.86 | 1.08x | 1.40x | 1.44x |
| `alcp.fp.mul_loop.n20000` | 16.30 | 12.90 | 1.26x | 1.49x | 1.32x |
| `alcp.fp.sub_loop.n20000` | 16.30 | 13.00 | 1.25x | 1.49x | 1.28x |
| `alcp.fp.div_loop.n5000` | 5.22 | 4.33 | 1.21x | 1.44x | 1.29x |
| `alcp.fp.mixed_int.n10000` | 18.57 | 13.48 | 1.38x | 1.76x | 1.32x |
| `alcp.fp.pow_loop.n5000` | 4.77 | 4.59 | 1.04x | 1.73x | 1.67x |
| `alcp.fp.eq_loop.n20000` | 2.35 | 2.33 | 1.00x | 2.59x | 2.50x |
| `alcp.f9.mul_add.n10000` | 37.55 | 35.66 | 1.05x | 1.17x | 1.19x |
| `alcp.gf256.mul_loop.n20000` | 97.68 | 88.17 | 1.11x | 1.11x | 1.02x |
| `alcp.gf256.add_int.n10000` | 25.74 | 22.48 | 1.15x | 1.24x | 1.08x |
| `alcp.fp.convert_int.n20000` | 6.16 | 5.82 | 1.06x | 0.74x | 0.72x |
| `alcp.fp.sum_elems.n20000` | 15.65 | 13.69 | 1.14x | 1.53x | 1.28x |
| `alcp.poly.from_list.d200` | 0.742 | 0.614 | 1.21x | 1.55x | 1.36x |
| `alcp.poly.from_elems.d200` | 0.832 | 0.663 | 1.26x | 1.82x | 1.44x |
| `alcp.poly.from_dict.d1000` | 1.58 | 1.56 | 1.02x | 1.18x | 1.21x |
| `alcp.poly.from_str.d40` | 16.81 | 0.013 | 1283x | 153x | 0.13x |
| `alcp.poly.from_str_distinct.d40` | 16.67 | 11.58 | 1.44x | 149x | 114x |
| `alcp.poly.from_str_f9.d40` | 20.95 | 0.008 | 2794x | 170x | 0.06x |
| `alcp.poly.parse_long.d500` | 211 | 38.39 | 5.51x | 3512x | 676x |
| `alcp.poly.from_str_repeat.n200` | 4.40 | 0.105 | 41.84x | 19.23x | 0.48x |
| `alcp.poly.add.d100` | 1.62 | 1.46 | 1.11x | 1.13x | 1.06x |
| `alcp.poly.sub.d100` | 1.48 | 1.44 | 1.03x | 1.05x | 1.05x |
| `alcp.poly.mul.d100` | 62.87 | 58.66 | 1.07x | 0.98x | 0.91x |
| `alcp.poly.floordiv.d100` | 31.24 | 31.45 | 0.99x | 1.00x | 1.03x |
| `alcp.poly.mod.d100` | 32.99 | 34.49 | 0.96x | 1.02x | 1.14x |
| `alcp.poly.small_mul.d5` | 2.00 | 2.02 | 0.99x | 1.12x | 1.19x |
| `alcp.poly.add_int.d100` | 1.22 | 1.13 | 1.08x | 1.18x | 1.10x |
| `alcp.poly.rmul_int.d100` | 3.43 | 3.15 | 1.09x | 1.16x | 1.10x |
| `alcp.poly.mul_elem_f9.d100` | 17.87 | 17.47 | 1.02x | 1.01x | 1.03x |
| `alcp.poly.eq.d100` | 0.299 | 0.301 | 0.99x | 1.17x | 1.27x |
| `alcp.poly.eval.d100.n500` | 6.47 | 5.79 | 1.12x | 1.04x | 0.93x |
| `alcp.poly.eval_elem.d100.n500` | 6.51 | 5.96 | 1.09x | 1.08x | 1.09x |
| `alcp.poly.factor.d40` | 10.98 | 9.98 | 1.10x | 0.95x | 0.98x |
| `alcp.poly.roots.d40` | 3.32 | 3.02 | 1.10x | 1.00x | 1.01x |
| `alcp.poly.str_unicode.d500` | 2.14 | 1.53 | 1.39x | 4.98x | 3.80x |
| `alcp.poly.str_ascii.d500` | 1.20 | 0.972 | 1.23x | 2.76x | 2.42x |
| `alcp.poly.str_latex.d500` | 1.40 | 0.959 | 1.46x | 3.41x | 2.43x |
| `alcp.poly.str_f9.d200` | 2.72 | 0.656 | 4.14x | 13.76x | 3.65x |
| `alcp.elem.str_gf256.n2000` | 7.97 | 0.340 | 23.45x | — | — |
| `alcp.pickle.poly.d500` | 2.12 | 1.98 | 1.07x | 1.14x | 1.15x |
| `alcp.pickle.elems.n1000` | 3.10 | 2.51 | 1.23x | 21.66x | 19.25x |

Media geométrica del speedup de la capa (41 escenarios alcp.*): 2.05x

<details><summary>8fcc151 — métodos nativos ligados en __slots__</summary>

| escenario | q00 mediana ms (IQR %) | q01 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `alcp.fp.add_loop.n20000` | 14.91 (18) | 13.55 (11) (!) | 1.10x | — |
| `alcp.fp.mul_loop.n20000` | 16.30 (13) | 15.81 (21) (!) | 1.03x | — |
| `alcp.fp.sub_loop.n20000` | 16.30 (21) | 14.38 (18) (!) | 1.13x | — |
| `alcp.fp.div_loop.n5000` | 5.22 (17) | 4.75 (21) (!) | 1.10x | — |
| `alcp.fp.mixed_int.n10000` | 18.57 (14) | 21.35 (15) (!) | 0.87x | — |
| `alcp.fp.pow_loop.n5000` | 4.77 (17) | 5.02 (18) (!) | 0.95x | — |
| `alcp.fp.eq_loop.n20000` | 2.35 (17) | 2.46 (19) (!) | 0.95x | — |
| `alcp.f9.mul_add.n10000` | 37.55 (17) | 36.81 (20) (!) | 1.02x | — |
| `alcp.gf256.mul_loop.n20000` | 97.68 (23) | 89.68 (7) (!) | 1.09x | — |
| `alcp.gf256.add_int.n10000` | 25.74 (16) | 25.12 (21) (!) | 1.02x | — |
| `alcp.fp.convert_int.n20000` | 6.16 (23) | 5.68 (6) (!) | 1.09x | — |
| `alcp.fp.sum_elems.n20000` | 15.65 (11) | 14.47 (5) (!) | 1.08x | — |

</details>

<details><summary>abc2f4d — camino rápido int en + - *</summary>

| escenario | q01 mediana ms (IQR %) | q02 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `alcp.fp.mixed_int.n10000` | 21.35 (15) | 15.66 (34) (!) | 1.36x | — |
| `alcp.gf256.add_int.n10000` | 25.12 (21) | 22.26 (11) (!) | 1.13x | — |
| `alcp.fp.sum_elems.n20000` | 14.47 (5) | 14.71 (27) (!) | 0.98x | — |

</details>

<details><summary>27873ea — analizador: x^k como monomio</summary>

| escenario | q02 mediana ms (IQR %) | q03 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `alcp.poly.from_str.d40` | 18.57 (44) | 13.97 (87) (!) | 1.33x | — |
| `alcp.poly.from_str_distinct.d40` | 16.56 (8) | 12.68 (19) (!) | 1.31x | — |
| `alcp.poly.from_str_f9.d40` | 19.85 (7) | 12.36 (20) (!) | 1.61x | — |
| `alcp.poly.parse_long.d500` | 197 (10) | 38.81 (15) (!) | 5.07x | — |
| `alcp.poly.from_str_repeat.n200` | 4.05 (11) | 4.22 (22) (!) | 0.96x | — |

</details>

<details><summary>29d5628 — caché de cadenas analizadas</summary>

| escenario | q03 mediana ms (IQR %) | q04 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `alcp.poly.from_str.d40` | 13.97 (87) | 0.013 (29) (!) | 1083x | — |
| `alcp.poly.from_str_distinct.d40` | 12.68 (19) | 12.38 (17) (!) | 1.02x | — |
| `alcp.poly.from_str_f9.d40` | 12.36 (20) | 0.007 (28) (!) | 1670x | — |
| `alcp.poly.from_str_repeat.n200` | 4.22 (22) | 0.098 (18) (!) | 42.88x | — |

</details>

<details><summary>0ab679a — _pretty: cada coeficiente formateado una vez</summary>

| escenario | q04 mediana ms (IQR %) | q05 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `alcp.poly.str_unicode.d500` | 2.15 (24) | 1.63 (18) (!) | 1.32x | — |
| `alcp.poly.str_ascii.d500` | 1.13 (20) | 1.07 (9) (!) | 1.06x | — |
| `alcp.poly.str_latex.d500` | 1.21 (16) | 1.13 (18) (!) | 1.07x | — |
| `alcp.poly.str_f9.d200` | 2.53 (38) | 0.811 (12) (!) | 3.12x | — |

</details>

<details><summary>a3399d9 — caché de textos de elementos</summary>

| escenario | q05 mediana ms (IQR %) | q06 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `alcp.poly.str_f9.d200` | 0.811 (12) | 0.708 (4) (!) | 1.14x | — |
| `alcp.elem.str_gf256.n2000` | 7.83 (19) | 0.332 (15) (!) | 23.54x | — |

</details>

<details><summary>cab2c2a — pickle de elementos más directo</summary>

| escenario | q06 mediana ms (IQR %) | q07 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `alcp.pickle.elems.n1000` | 3.03 (11) | 2.48 (11) (!) | 1.22x | — |

</details>

<details><summary>a833bcb — caminos rápidos int/elemento en conversión y f(v)</summary>

| escenario | q07 mediana ms (IQR %) | q08 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `alcp.poly.add_int.d100` | 1.23 (13) | 1.26 (22) (!) | 0.98x | — |
| `alcp.poly.rmul_int.d100` | 3.28 (12) | 3.37 (14) (!) | 0.97x | — |
| `alcp.poly.mul_elem_f9.d100` | 17.32 (32) | 18.14 (92) (!) | 0.95x | — |
| `alcp.poly.eval.d100.n500` | 6.08 (15) | 6.10 (11) (!) | 1.00x | — |
| `alcp.poly.eval_elem.d100.n500` | 5.99 (10) | 6.06 (9) | 0.99x | — |

</details>

<details><summary>e683a88 — listas de coeficientes con caminos rápidos</summary>

| escenario | q08 mediana ms (IQR %) | q09 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `alcp.poly.from_list.d200` | 0.752 (18) | 0.581 (23) (!) | 1.30x | — |
| `alcp.poly.from_elems.d200` | 0.792 (14) | 0.629 (24) (!) | 1.26x | — |

</details>

<details><summary>06d922a — caminos rápidos // % divmod</summary>

| escenario | q09 mediana ms (IQR %) | q10 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `alcp.poly.floordiv.d100` | 30.08 (15) | 28.48 (13) (!) | 1.06x | — |
| `alcp.poly.mod.d100` | 30.04 (9) | 29.58 (18) (!) | 1.02x | — |

</details>

<details><summary>0438b74 — resultado en línea en + - * / (DUDOSA)</summary>

| escenario | q10 mediana ms (IQR %) | q11 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `alcp.fp.add_loop.n20000` | 13.57 (12) | 13.03 (16) (!) | 1.04x | — |
| `alcp.fp.mul_loop.n20000` | 13.60 (9) | 13.50 (19) (!) | 1.01x | — |
| `alcp.fp.sub_loop.n20000` | 13.22 (10) | 13.60 (13) (!) | 0.97x | — |
| `alcp.fp.div_loop.n5000` | 4.14 (11) | 4.07 (13) (!) | 1.02x | — |
| `alcp.fp.mixed_int.n10000` | 13.82 (8) | 13.33 (12) (!) | 1.04x | — |
| `alcp.fp.pow_loop.n5000` | 4.38 (16) | 4.40 (15) (!) | 0.99x | — |
| `alcp.fp.eq_loop.n20000` | 2.33 (6) | 2.21 (22) (!) | 1.05x | — |
| `alcp.f9.mul_add.n10000` | 34.27 (14) | 33.03 (19) (!) | 1.04x | — |
| `alcp.gf256.mul_loop.n20000` | 86.29 (13) | 84.91 (12) (!) | 1.02x | — |
| `alcp.gf256.add_int.n10000` | 21.29 (11) | 21.71 (11) (!) | 0.98x | — |
| `alcp.fp.convert_int.n20000` | 5.69 (22) | 5.68 (13) (!) | 1.00x | — |
| `alcp.fp.sum_elems.n20000` | 14.07 (6) | 13.33 (12) (!) | 1.06x | — |

</details>

<details><summary>c70bed5 — _pretty: potencias en línea</summary>

| escenario | q11 mediana ms (IQR %) | q12 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `alcp.poly.str_unicode.d500` | 1.61 (24) | 1.53 (10) (!) | 1.05x | — |
| `alcp.poly.str_ascii.d500` | 1.17 (11) | 0.972 (8) (!) | 1.21x | — |
| `alcp.poly.str_latex.d500` | 1.12 (15) | 0.959 (15) (!) | 1.17x | — |
| `alcp.poly.str_f9.d200` | 0.703 (8) | 0.656 (14) (!) | 1.07x | — |

</details>

<details><summary>ff9ad21 (exp/) — EXPERIMENTO tokenizador regex</summary>

| escenario | q03 mediana ms (IQR %) | q13 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `alcp.poly.from_str.d40` | 13.97 (87) | 13.29 (28) (!) | 1.05x | — |
| `alcp.poly.from_str_distinct.d40` | 12.68 (19) | 12.43 (106) (!) | 1.02x | — |
| `alcp.poly.from_str_f9.d40` | 12.36 (20) | 11.88 (7) (!) | 1.04x | — |
| `alcp.poly.parse_long.d500` | 38.81 (15) | 37.72 (9) (!) | 1.03x | — |
| `alcp.poly.from_str_repeat.n200` | 4.22 (22) | 3.70 (12) (!) | 1.14x | — |

</details>

## 6. Mediciones focalizadas de las técnicas dudosas (6 rondas x >= 5 repeticiones, con corriente)

### 2297387 Barrett sin desplazamientos variables (a4 -> a5)

Media geométrica 1.07x en 24 escenarios (rango 0.98-1.23); checksums distintos de la referencia: 0; errores: 0.

| escenario | a4 mediana ms (IQR %) | a5 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p7.elem_muladd.n400000` | 2.24 (37) | 1.84 (12) (!) | 1.22x | — |
| `fp.p7.elem_inv.n100000` | 0.375 (28) | 0.370 (12) (!) | 1.01x | — |
| `fp.p7.elem_pow.n100000` | 1.43 (31) | 1.34 (30) (!) | 1.07x | — |
| `fp.p7.mul_school.d512` | 0.204 (22) | 0.188 (20) (!) | 1.09x | — |
| `fp.p7.divmod.d2000_1000` | 0.781 (34) | 0.766 (13) (!) | 1.02x | — |
| `fp.p7.factor.d62` | 0.351 (17) | 0.346 (14) (!) | 1.02x | — |
| `fp.p101.elem_muladd.n400000` | 2.24 (8) | 1.85 (13) (!) | 1.21x | — |
| `fp.p101.elem_inv.n100000` | 0.709 (10) | 0.697 (9) (!) | 1.02x | — |
| `fp.p101.elem_pow.n100000` | 3.95 (21) | 3.22 (8) (!) | 1.23x | — |
| `fp.p101.mul_school.d512` | 0.228 (13) | 0.213 (5) (!) | 1.07x | — |
| `fp.p101.divmod.d2000_1000` | 0.874 (9) | 0.855 (9) | 1.02x | — |
| `fp.p101.factor.d100` | 1.85 (10) | 1.81 (8) (!) | 1.02x | — |
| `fp.ntt.elem_muladd.n400000` | 3.06 (9) | 2.67 (7) | 1.15x | — |
| `fp.ntt.elem_inv.n100000` | 3.02 (6) | 3.10 (13) (!) | 0.98x | — |
| `fp.ntt.elem_pow.n50000` | 12.95 (8) | 11.00 (52) (!) | 1.18x | — |
| `fp.ntt.mul_school.d512` | 0.229 (7) | 0.228 (24) (!) | 1.01x | — |
| `fp.ntt.divmod.d2000_1000` | 0.898 (10) | 0.910 (32) (!) | 0.99x | — |
| `fp.ntt.factor.d40` | 1.53 (16) | 1.49 (20) (!) | 1.03x | — |
| `fp.m61.elem_muladd.n400000` | 2.72 (22) | 2.38 (8) (!) | 1.14x | — |
| `fp.m61.elem_inv.n100000` | 8.35 (5) | 8.35 (4) | 1.00x | — |
| `fp.m61.elem_pow.n20000` | 11.69 (5) | 10.57 (9) | 1.11x | — |
| `fp.m61.mul_school.d512` | 0.249 (8) | 0.225 (9) | 1.11x | — |
| `fp.m61.divmod.d2000_1000` | 0.931 (7) | 0.930 (11) (!) | 1.00x | — |
| `fp.m61.factor.d26` | 1.33 (7) | 1.35 (12) (!) | 0.98x | — |

### dcb905c suma/resta de Poly (a2 -> a3)

Media geométrica 1.17x en 30 escenarios (rango 0.91-1.50); checksums distintos de la referencia: 0; errores: 0.

| escenario | a2 mediana ms (IQR %) | a3 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p7.mul_karatsuba.d512` | 5.28 (10) | 4.15 (8) (!) | 1.27x | — |
| `fp.p7.mul_karatsuba.d4096` | 151 (7) | 113 (7) | 1.33x | — |
| `fp.p7.gcd.d500` | 16.87 (11) | 16.69 (8) (!) | 1.01x | — |
| `fp.p7.gcd_ext.d400` | 38.44 (21) | 33.14 (4) (!) | 1.16x | — |
| `fp.p101.mul_karatsuba.d512` | 5.44 (8) | 4.04 (12) (!) | 1.34x | — |
| `fp.p101.mul_karatsuba.d4096` | 151 (6) | 118 (20) (!) | 1.28x | — |
| `fp.p101.gcd.d500` | 19.78 (8) | 21.25 (12) (!) | 0.93x | — |
| `fp.p101.gcd_ext.d400` | 45.29 (16) | 41.11 (9) (!) | 1.10x | — |
| `fp.ntt.mul_karatsuba.d512` | 5.17 (21) | 4.16 (13) (!) | 1.24x | — |
| `fp.ntt.mul_karatsuba.d4096` | 155 (7) | 128 (47) (!) | 1.21x | — |
| `fp.ntt.gcd.d500` | 20.45 (12) | 22.42 (10) (!) | 0.91x | — |
| `fp.ntt.gcd_ext.d400` | 526 (18) | 494 (5) (!) | 1.07x | — |
| `fp.m61.mul_karatsuba.d512` | 5.36 (6) | 4.36 (6) | 1.23x | — |
| `fp.m61.mul_karatsuba.d4096` | 165 (12) | 126 (12) (!) | 1.31x | — |
| `fp.m61.gcd.d500` | 25.80 (7) | 27.31 (20) (!) | 0.94x | — |
| `fp.m61.gcd_ext.d400` | 53.37 (7) | 48.49 (8) | 1.10x | — |
| `fp.m127.mul_karatsuba.d512` | 7.75 (9) | 5.18 (7) | 1.50x | — |
| `fp.m127.mul_karatsuba.d2048` | 72.11 (14) | 50.33 (20) (!) | 1.43x | — |
| `fp.m127.gcd.d400` | 15.74 (8) | 15.35 (7) | 1.03x | — |
| `fp.m127.gcd_ext.d300` | 32.06 (7) | 25.34 (8) | 1.27x | — |
| `fp.p256.mul_karatsuba.d512` | 8.91 (8) | 6.23 (7) | 1.43x | — |
| `fp.p256.mul_karatsuba.d2048` | 85.89 (11) | 59.19 (10) (!) | 1.45x | — |
| `fp.p256.gcd.d300` | 11.91 (23) | 11.46 (7) (!) | 1.04x | — |
| `fp.p256.gcd_ext.d200` | 131 (11) | 131 (15) (!) | 1.00x | — |
| `fq.q3e30.mul_karatsuba.d128` | 65.47 (14) | 61.58 (12) (!) | 1.06x | — |
| `fq.q2e8.mul_karatsuba.d512` | 74.48 (14) | 64.86 (22) (!) | 1.15x | — |
| `fq.q7e4.mul_karatsuba.d512` | 36.94 (9) | 30.18 (13) (!) | 1.22x | — |
| `fq.q7e4.gcd.d200` | 33.28 (8) | 32.77 (18) (!) | 1.02x | — |
| `fq.m61e2.mul_karatsuba.d512` | 18.35 (12) | 13.95 (10) (!) | 1.32x | — |
| `fq.m61e2.gcd.d150` | 16.69 (9) | 15.96 (6) | 1.05x | — |

### 9bfc1e7 Karatsuba sobre vistas (m2 -> m3)

Media geométrica 1.18x en 14 escenarios (rango 1.01-1.34); checksums distintos de la referencia: 0; errores: 0.

| escenario | m2 mediana ms (IQR %) | m3 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p7.mul_karatsuba.d512` | 12.25 (15) | 9.35 (6) (!) | 1.31x | — |
| `fp.p7.mul_karatsuba.d4096` | 325 (5) | 262 (6) | 1.24x | — |
| `fp.p101.mul_karatsuba.d512` | 12.82 (3) | 11.07 (5) | 1.16x | — |
| `fp.p101.mul_karatsuba.d4096` | 380 (11) | 319 (11) (!) | 1.19x | — |
| `fp.m61.mul_karatsuba.d512` | 13.22 (8) | 11.07 (6) | 1.19x | — |
| `fp.m61.mul_karatsuba.d4096` | 393 (8) | 318 (11) (!) | 1.23x | — |
| `fp.m127.mul_karatsuba.d512` | 15.67 (8) | 11.73 (7) | 1.34x | — |
| `fp.m127.mul_karatsuba.d2048` | 147 (8) | 110 (10) (!) | 1.33x | — |
| `fp.p256.mul_karatsuba.d512` | 19.53 (11) | 15.72 (13) (!) | 1.24x | — |
| `fp.p256.mul_karatsuba.d2048` | 191 (11) | 150 (7) (!) | 1.28x | — |
| `fq.q3e30.mul_karatsuba.d128` | 458 (6) | 455 (18) (!) | 1.01x | — |
| `fq.q2e8.mul_karatsuba.d512` | 393 (7) | 370 (10) | 1.06x | — |
| `fq.q7e4.mul_karatsuba.d512` | 225 (7) | 218 (10) | 1.03x | — |
| `fq.m61e2.mul_karatsuba.d512` | 85.69 (3) | 82.34 (9) | 1.04x | — |

### b0a9fe7 NTT paralela (m8 -> m9)

Media geométrica 1.08x en 8 escenarios (rango 1.00-1.38); checksums distintos de la referencia: 0; errores: 0.

| escenario | m8 mediana ms (IQR %) | m9 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.ntt.mul_fft.d4096` | 2.36 (5) | 2.12 (13) (!) | 1.11x | — |
| `fq.m61e2.mul_karatsuba.d512` | 30.24 (7) | 30.27 (9) | 1.00x | — |
| `fq.m61e2.mul_fft.d1024` | 85.09 (10) | 61.74 (8) (!) | 1.38x | — |
| `fq.m61e2.gcd.d150` | 45.29 (7) | 43.36 (8) | 1.04x | — |
| `fast.ntt.toeplitz_lower_inverse.n1024` | 1.50 (8) | 1.47 (10) | 1.02x | — |
| `fast.ntt.toeplitz_lower_inverse.n8192` | 14.97 (8) | 14.37 (8) | 1.04x | — |
| `fast.ntt.toeplitz_lower_vec.n4096` | 1.77 (10) | 1.69 (18) (!) | 1.05x | — |
| `fast.ntt.toeplitz_divmod.d2000_1000` | 60.29 (8) | 58.74 (6) | 1.03x | — |

### 3660189 umbrales de despacho Auto (m5 -> m6)

Media geométrica 1.02x en 22 escenarios (rango 0.94-2.01); checksums distintos de la referencia: 0; errores: 0.

| escenario | m5 mediana ms (IQR %) | m6 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p7.mul_auto.d2048` | 79.78 (12) | 79.14 (6) (!) | 1.01x | — |
| `fp.p7.gcd.d500` | 24.13 (2) | 25.24 (20) (!) | 0.96x | — |
| `fp.p7.gcd_ext.d400` | 48.86 (4) | 51.18 (12) (!) | 0.95x | — |
| `fp.p101.mul_auto.d2048` | 92.13 (4) | 91.83 (15) (!) | 1.00x | — |
| `fp.p101.gcd.d500` | 30.00 (3) | 30.32 (15) (!) | 0.99x | — |
| `fp.p101.gcd_ext.d400` | 57.95 (5) | 59.55 (13) (!) | 0.97x | — |
| `fp.ntt.mul_auto.d2048` | 0.909 (4) | 0.972 (13) (!) | 0.94x | — |
| `fp.ntt.gcd.d500` | 31.31 (3) | 33.22 (12) (!) | 0.94x | — |
| `fp.ntt.gcd_ext.d400` | 44.15 (8) | 47.07 (39) (!) | 0.94x | — |
| `fp.m61.mul_auto.d2048` | 93.73 (6) | 92.87 (11) (!) | 1.01x | — |
| `fp.m61.gcd.d500` | 31.18 (9) | 31.76 (15) (!) | 0.98x | — |
| `fp.m61.gcd_ext.d400` | 60.52 (13) | 60.72 (8) (!) | 1.00x | — |
| `fp.m127.mul_auto.d2048` | 103 (10) | 101 (12) (!) | 1.02x | — |
| `fp.m127.gcd.d400` | 21.30 (17) | 21.18 (4) (!) | 1.01x | — |
| `fp.m127.gcd_ext.d300` | 40.86 (15) | 38.75 (11) (!) | 1.05x | — |
| `fp.p256.mul_auto.d2048` | 9.89 (10) | 9.76 (9) (!) | 1.01x | — |
| `fp.p256.gcd.d300` | 16.12 (11) | 16.68 (14) (!) | 0.97x | — |
| `fp.p256.gcd_ext.d200` | 45.95 (24) | 22.84 (7) (!) | 2.01x | — |
| `fq.q7e4.gcd.d200` | 105 (10) | 105 (15) (!) | 0.99x | — |
| `fq.q7e4.factor.d20` | 56.33 (16) | 57.20 (10) (!) | 0.98x | — |
| `fq.m61e2.mul_auto.d2048` | 166 (19) | 164 (13) (!) | 1.01x | — |
| `fq.m61e2.gcd.d150` | 39.29 (10) | 39.22 (8) (!) | 1.00x | — |

### e8e07e6 Karatsuba paralelo (m7 -> m8)

Media geométrica 1.23x en 22 escenarios (rango 0.96-3.13); checksums distintos de la referencia: 0; errores: 0.

| escenario | m7 mediana ms (IQR %) | m8 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p7.mul_auto.d2048` | 78.67 (3) | 25.13 (14) (!) | 3.13x | — |
| `fp.p7.gcd.d500` | 24.53 (6) | 24.27 (25) (!) | 1.01x | — |
| `fp.p7.gcd_ext.d400` | 50.02 (4) | 50.96 (10) | 0.98x | — |
| `fp.p101.mul_auto.d2048` | 92.98 (7) | 30.88 (9) | 3.01x | — |
| `fp.p101.gcd.d500` | 30.97 (7) | 30.11 (18) (!) | 1.03x | — |
| `fp.p101.gcd_ext.d400` | 58.98 (10) | 61.65 (18) (!) | 0.96x | — |
| `fp.ntt.mul_auto.d2048` | 0.944 (7) | 0.941 (8) | 1.00x | — |
| `fp.ntt.gcd.d500` | 31.75 (7) | 31.77 (10) (!) | 1.00x | — |
| `fp.ntt.gcd_ext.d400` | 45.13 (8) | 44.64 (24) (!) | 1.01x | — |
| `fp.m61.mul_auto.d2048` | 97.38 (17) | 31.41 (14) (!) | 3.10x | — |
| `fp.m61.gcd.d500` | 32.10 (11) | 30.60 (12) (!) | 1.05x | — |
| `fp.m61.gcd_ext.d400` | 61.34 (6) | 61.05 (6) | 1.00x | — |
| `fp.m127.mul_auto.d2048` | 102 (11) | 34.56 (10) (!) | 2.94x | — |
| `fp.m127.gcd.d400` | 21.79 (8) | 22.11 (5) | 0.99x | — |
| `fp.m127.gcd_ext.d300` | 40.29 (13) | 40.65 (7) (!) | 0.99x | — |
| `fp.p256.mul_auto.d2048` | 9.73 (5) | 9.54 (8) | 1.02x | — |
| `fp.p256.gcd.d300` | 16.46 (5) | 16.35 (6) | 1.01x | — |
| `fp.p256.gcd_ext.d200` | 22.14 (12) | 22.32 (4) (!) | 0.99x | — |
| `fq.q7e4.gcd.d200` | 107 (9) | 105 (6) | 1.02x | — |
| `fq.q7e4.factor.d20` | 56.98 (5) | 57.53 (5) | 0.99x | — |
| `fq.m61e2.mul_auto.d2048` | 167 (14) | 169 (11) (!) | 0.99x | — |
| `fq.m61e2.gcd.d150` | 39.40 (10) | 39.83 (6) | 0.99x | — |

### e8e07e6 con ALCP_THREADS=1 (coste de la infraestructura paralela) (m7 -> m8t1)

Media geométrica 1.01x en 22 escenarios (rango 0.98-1.04); checksums distintos de la referencia: 0; errores: 0.

| escenario | m7 mediana ms (IQR %) | m8t1 mediana ms (IQR %) | paso | vs ref |
|---|---:|---:|---:|---:|
| `fp.p7.mul_auto.d2048` | 78.67 (3) | 79.79 (6) | 0.99x | — |
| `fp.p7.gcd.d500` | 24.53 (6) | 24.65 (12) (!) | 1.00x | — |
| `fp.p7.gcd_ext.d400` | 50.02 (4) | 50.98 (10) | 0.98x | — |
| `fp.p101.mul_auto.d2048` | 92.98 (7) | 93.14 (7) | 1.00x | — |
| `fp.p101.gcd.d500` | 30.97 (7) | 29.98 (6) | 1.03x | — |
| `fp.p101.gcd_ext.d400` | 58.98 (10) | 59.05 (7) (!) | 1.00x | — |
| `fp.ntt.mul_auto.d2048` | 0.944 (7) | 0.935 (3) | 1.01x | — |
| `fp.ntt.gcd.d500` | 31.75 (7) | 31.06 (4) | 1.02x | — |
| `fp.ntt.gcd_ext.d400` | 45.13 (8) | 44.17 (3) | 1.02x | — |
| `fp.m61.mul_auto.d2048` | 97.38 (17) | 94.47 (9) (!) | 1.03x | — |
| `fp.m61.gcd.d500` | 32.10 (11) | 31.37 (5) (!) | 1.02x | — |
| `fp.m61.gcd_ext.d400` | 61.34 (6) | 60.82 (7) | 1.01x | — |
| `fp.m127.mul_auto.d2048` | 102 (11) | 99.51 (9) (!) | 1.02x | — |
| `fp.m127.gcd.d400` | 21.79 (8) | 21.75 (6) | 1.00x | — |
| `fp.m127.gcd_ext.d300` | 40.29 (13) | 38.79 (5) (!) | 1.04x | — |
| `fp.p256.mul_auto.d2048` | 9.73 (5) | 9.58 (5) | 1.02x | — |
| `fp.p256.gcd.d300` | 16.46 (5) | 16.15 (4) | 1.02x | — |
| `fp.p256.gcd_ext.d200` | 22.14 (12) | 22.10 (4) (!) | 1.00x | — |
| `fq.q7e4.gcd.d200` | 107 (9) | 103 (5) | 1.04x | — |
| `fq.q7e4.factor.d20` | 56.98 (5) | 56.15 (7) | 1.01x | — |
| `fq.m61e2.mul_auto.d2048` | 167 (14) | 162 (6) (!) | 1.03x | — |
| `fq.m61e2.gcd.d150` | 39.40 (10) | 38.38 (4) | 1.03x | — |


---

## Apéndice A. Mismas tablas con batería (comprobación de coherencia)

Datos de `extremos/limpio/` (rondas sin suspensión), `aritmetica/` y `multiplicacion/`, medidos con batería (CPU limitada a ~40 % del rendimiento con corriente). Los tiempos absolutos no son comparables con los de las tablas principales; los speedups relativos sí.

<!-- datos: extremos -->
### 1. Extremos: suite completa (128 escenarios)

| variante | escenarios | checksums = ref | errores | media geom. vs ref | escenarios con IQR > 10 % |
|---|---:|---:|---:|---:|---:|
| `b0` | 128 | 128 | 0 | 1.00x (ref) | 98 |
| `h0` | 128 | 128 | 0 | 1.08x | 82 |
| `a6` | 128 | 128 | 0 | 19.35x | 94 |
| `a6gmp` | 128 | 128 | 0 | 2.77x | 107 |
| `m9` | 128 | 128 | 0 | 1.92x | 85 |
| `m9t1` | 128 | 128 | 0 | 1.50x | 56 |

#### Media geométrica del speedup frente a la referencia, por grupo

| grupo | n | h0 | a6 | a6gmp | m9 | m9t1 |
|---|---:|---:|---:|---:|---:|---:|
| `fast.ntt` | 8 | 0.94x | 16.12x | 0.92x | 10.01x | 10.00x |
| `fast.p7` | 2 | 1.05x | 56.63x | 1.96x | 1.79x | 1.21x |
| `fp.m127` | 14 | 1.18x | 2.90x | 2.90x | 1.58x | 1.19x |
| `fp.m61` | 14 | 1.09x | 40.57x | 2.97x | 1.47x | 1.12x |
| `fp.ntt` | 16 | 1.03x | 37.22x | 2.07x | 2.93x | 2.38x |
| `fp.p101` | 14 | 1.14x | 57.74x | 3.18x | 1.58x | 1.21x |
| `fp.p256` | 16 | 1.09x | 2.07x | 2.15x | 2.14x | 1.70x |
| `fp.p7` | 14 | 1.10x | 52.94x | 3.06x | 1.46x | 1.09x |
| `fq.build` | 2 | 1.05x | 45.11x | 3.34x | 1.00x | 1.03x |
| `fq.m61e2` | 10 | 1.05x | 14.01x | 3.71x | 1.76x | 1.13x |
| `fq.q2e8` | 6 | 1.10x | 44.54x | 4.67x | 1.41x | 1.14x |
| `fq.q3e30` | 6 | 1.07x | 46.12x | 4.63x | 1.32x | 1.06x |
| `fq.q7e4` | 6 | 1.09x | 30.23x | 5.13x | 1.34x | 1.14x |

#### Por operación

| operación | n | h0 | a6 | a6gmp | m9 | m9t1 |
|---|---:|---:|---:|---:|---:|---:|
| `build` | 2 | 1.05x | 45.11x | 3.34x | 1.00x | 1.03x |
| `dft` | 3 | 0.90x | 11.00x | 0.87x | 10.05x | 10.44x |
| `divmod` | 8 | 1.03x | 65.04x | 9.87x | 1.02x | 1.05x |
| `elem_inv` | 10 | 1.14x | 8.16x | 1.56x | 1.08x | 1.13x |
| `elem_muladd` | 10 | 1.16x | 11.24x | 1.60x | 1.09x | 1.14x |
| `elem_pow` | 6 | 1.11x | 2.73x | 0.96x | 1.08x | 1.09x |
| `factor` | 10 | 1.07x | 17.60x | 2.77x | 1.07x | 1.11x |
| `gcd` | 8 | 1.06x | 12.25x | 1.80x | 0.94x | 1.03x |
| `gcd_ext` | 6 | 1.13x | 10.48x | 1.47x | 2.17x | 2.19x |
| `idft` | 1 | 0.99x | 10.37x | 0.90x | 8.38x | 8.77x |
| `is_irreducible` | 7 | 1.05x | 30.31x | 3.97x | 0.98x | 1.05x |
| `mul_auto` | 7 | 1.12x | 15.17x | 2.12x | 8.02x | 2.54x |
| `mul_fft` | 6 | 1.06x | 5.74x | 1.15x | 10.57x | 6.24x |
| `mul_karatsuba` | 16 | 1.08x | 29.19x | 3.93x | 5.40x | 1.24x |
| `mul_school` | 16 | 1.13x | 60.74x | 9.76x | 1.01x | 1.08x |
| `pow_mod` | 6 | 1.07x | 31.72x | 4.25x | 1.23x | 1.25x |
| `toeplitz_divmod` | 2 | 1.04x | 31.15x | 1.18x | 1.16x | 1.10x |
| `toeplitz_lower_inverse` | 3 | 1.00x | 35.90x | 1.35x | 10.82x | 8.53x |
| `toeplitz_lower_vec` | 1 | 0.89x | 23.58x | 0.99x | 22.10x | 19.53x |

#### Coste de robustez (5f49217 frente a la referencia): los 10 escenarios más lentos y los 5 más rápidos

| escenario | ref mediana ms (IQR %) | 5f49217 mediana ms (IQR %) | speedup |
|---|---:|---:|---:|
| `fast.ntt.dft.k11` | 38.53 (5) | 47.32 (26) | 0.81x |
| `fp.p7.elem_pow.n100000` | 35.12 (72) | 40.81 (20) | 0.86x |
| `fast.ntt.toeplitz_lower_vec.n4096` | 82.22 (9) | 92.04 (16) | 0.89x |
| `fp.m61.divmod.d2000_1000` | 299 (5) | 326 (23) | 0.92x |
| `fast.ntt.dft.k10` | 34.69 (14) | 37.51 (4) | 0.92x |
| `fq.m61e2.divmod.d500_250` | 138 (27) | 148 (19) | 0.93x |
| `fp.p101.mul_karatsuba.d4096` | 888 (11) | 938 (8) | 0.95x |
| `fq.build.q7e20` | 42.71 (18) | 45.00 (8) | 0.95x |
| `fq.m61e2.factor.d7` | 389 (9) | 408 (24) | 0.95x |
| `fast.ntt.toeplitz_lower_inverse.n8192` | 747 (23) | 783 (9) | 0.95x |
| `fp.p101.mul_school.d64` | 28.66 (41) | 21.43 (24) | 1.34x |
| `fp.m127.elem_pow.n10000` | 75.05 (37) | 55.23 (5) | 1.36x |
| `fp.p7.elem_inv.n100000` | 40.60 (27) | 29.08 (32) | 1.40x |
| `fp.p101.elem_pow.n100000` | 60.44 (32) | 43.12 (20) | 1.40x |
| `fp.m127.mul_school.d64` | 37.90 (14) | 22.69 (4) | 1.67x |

<!-- datos: aritmetica -->
### 2.1. Cadena aritmetica

| variante | escenarios | checksums = ref (b0) | = 5f49217 | errores | IQR > 10 % |
|---|---:|---:|---:|---:|---:|
| `a1` | 37 | 37 | 37 | 0 | 12 |
| `a2` | 37 | 37 | 37 | 0 | 18 |
| `a3` | 37 | 37 | 37 | 0 | 20 |
| `a4` | 37 | 37 | 37 | 0 | 17 |
| `a4gmp` | 37 | 37 | 37 | 0 | 8 |
| `a5` | 37 | 37 | 37 | 0 | 14 |
| `a6` | 37 | 37 | 37 | 0 | 18 |
| `h0` | 37 | 37 | 37 | 0 | 5 |

| commit | técnica | n afect. | media geom. paso (afectados) | rango paso | media geom. paso (control) | media geom. acumulada vs ref (afectados) |
|---|---|---:|---:|---:|---:|---:|
| 489728b | GMP en sitio en Fp y Fq | 18 | 1.84x | 0.98-3.46 | 1.00x | 2.02x |
| 0de8cf9 | Reducción perezosa en mul_school/divmod (LazyAccumulation) | 25 | 4.44x | 1.22-11.52 | 1.01x | 6.53x |
| dcb905c | Suma/resta de Poly sin relleno ni copias | 5 | 1.06x | 0.94-1.27 | 1.00x | 2.99x |
| 65a5f1c (+5c444c0, 4fb4f1a) | Backend Fp64/Fq64 (auto) | 31 | 8.68x | 0.71-29.08 | 0.98x | 34.75x |
| 65a5f1c, --backend gmp | FqT<Base> con backend GMP forzado (coste de la plantilla) | 37 | 1.01x | 0.86-1.24 | — | 4.13x |
| 2297387 | Fp64: Barrett sin desplazamientos variables | 23 | 1.05x | 0.85-1.28 | 1.03x | 20.65x |
| 8f1d896 | Fp64::pow en forma de Montgomery | 3 | 1.61x | 1.58-1.66 | 1.02x | 3.39x |








<!-- datos: multiplicacion -->
### 2.2. Cadena multiplicacion

| variante | escenarios | checksums = ref (b0) | = 5f49217 | errores | IQR > 10 % |
|---|---:|---:|---:|---:|---:|
| `h0` | 35 | 35 | 35 | 0 | 22 |
| `m1` | 35 | 35 | 35 | 0 | 15 |
| `m2` | 35 | 35 | 35 | 0 | 8 |
| `m3` | 35 | 35 | 35 | 0 | 19 |
| `m4` | 35 | 35 | 35 | 0 | 19 |
| `m5` | 35 | 35 | 35 | 0 | 18 |
| `m6` | 35 | 35 | 35 | 0 | 19 |
| `m7` | 35 | 35 | 35 | 0 | 26 |
| `m8` | 35 | 35 | 35 | 0 | 24 |
| `m8t1` | 35 | 35 | 35 | 0 | 24 |
| `m9` | 35 | 35 | 35 | 0 | 27 |
| `m9t1` | 35 | 35 | 35 | 0 | 35 |

| commit | técnica | n afect. | media geom. paso (afectados) | rango paso | media geom. paso (control) | media geom. acumulada vs ref (afectados) |
|---|---|---:|---:|---:|---:|---:|
| d4580af | FFT: tablas directas, bit-reversal precalculado, mariposa sin temporal | 17 | 1.10x | 0.99-1.24 | 1.01x | 1.12x |
| f68a18b | Cuadrado con una sola transformada | 2 | 1.05x | 1.04-1.06 | 1.00x | 1.15x |
| 9bfc1e7 | Karatsuba sobre vistas | 10 | 1.11x | 0.95-1.25 | 1.00x | 1.24x |
| 1ad8c02 | NTT en uint64_t + Shoup (p < 2^62) | 15 | 4.06x | 0.94-32.83 | 0.99x | 4.42x |
| 70e3453 | NTT multi-limb mpn/Montgomery (p >= 2^62) | 4 | 2.60x | 1.23-3.95 | 1.01x | 3.45x |
| 3660189 | Umbrales de despacho Auto | 11 | 1.10x | 0.90-1.93 | 0.97x | 2.68x |
| 06e46b5 | fast::dft/idft en palabras de máquina | 4 | 8.89x | 5.99-10.45 | 0.98x | 7.92x |
| e8e07e6 | Karatsuba paralelo (std::async) | 13 | 2.17x | 0.93-5.31 | 1.05x | 3.54x |
| b0a9fe7 | NTT paralela | 12 | 1.40x | 0.77-4.78 | 0.96x | 6.14x |
| e8e07e6 con ALCP_THREADS=1 | Karatsuba paralelo limitado a 1 hilo (sobrecoste) | 13 | 1.02x | 0.84-1.20 | 1.08x | 1.66x |
| b0a9fe7 con ALCP_THREADS=1 | NTT paralela limitada a 1 hilo (sobrecoste) | 12 | 1.03x | 0.96-1.09 | 1.03x | 3.60x |












