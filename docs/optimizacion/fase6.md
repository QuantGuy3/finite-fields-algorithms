# Fase 6: Newton con producto medio, p pequeño, DDF sin pow_mod repetido y binding directo

Fecha: 8 de noviembre de 2026. Rama `fase6`, sobre `mejoras-cpp-pybind` en `2b36cde` (la
fase 5 con la corrección de `gcd_ext`). Implementa las propuestas 1, 5, 6 y 8 de
[investigacion-fase6.md](investigacion-fase6.md); las demás (half-GCD V3e, `FpN<2>`, NTT
radix 4/AVX2, FqTab/FqSwar) quedan pendientes.

## Resumen

| # | Cambio | Dónde | Medido (en serie) |
|---|---|---|---|
| 1 | División de Newton con producto medio; `Reducer` y `Composer` con transformadas NTT precalculadas | `polynomial.hpp` (`namespace medio`, `inverse_series`, `divmod_with_inverse`, `Reducer`, `Composer`) | `fp.ntt.divmod` 1,14x; `pow_mod` y `factor` en F_998244353 hasta 1,8x |
| 5 | Fp64 con p < 2^32: acumulador de 64 bits y producto de 64 bits sin vigilar el acarreo; Kronecker sin `u128 % p` y con desempaquetado directo en F_q | `fp64.hpp`, `kronecker.hpp`, `polynomial.hpp`, `fq.hpp`, `field.hpp` | F_7 y F_101: grupos 1,54x y 1,49x; `mul_school` 3,2-3,3x; F_(3^30) 1,40x |
| 6 | DDF de la fase 4: tras cada factor encontrado compone de inmediato si compensa, en lugar de rehacer un `pow_mod` | `polynomial.hpp` (Frobenius: `gastado_`, `compone_ya`, `rng_intacto`) | `fp.m127.factor.d19` 1,31x; `fp.m61.factor.d26` 1,42x |
| 8 | Conversión directa `int` ↔ `u64` en el binding para Fp64 | `python-bindings/bindings.cpp` | `K.poly(lista)` 6,1x / 11x; `coeffs()` 4,4x / 3,6x; elementos 1,6x |

**Suite C++** (128 escenarios, intercalada con la base, 3 rondas): **1,17x** de media
geométrica y **todos los checksums iguales** a la base. Medida con el portátil **con
batería**; la comparación es justa porque las variantes se intercalan, pero las cifras
absolutas no son comparables con las de fases medidas enchufado.

## Método

1. **Implementación en paralelo.** Un agente por propuesta, cada uno en su git worktree,
   con un fichero de pruebas propio que compara con la versión anterior y con el
   generador aleatorio en frío.
2. **Revisión adversaria.** Un revisor independiente por propuesta intentó refutar la
   implementación: programas diferenciales contra la cabecera de `21e9fe6`, cotas de
   desbordamiento, consumo del generador y `build.py all`. Las propuestas 1, 5 y 8 se
   aprobaron. La 6 no: sin la salvaguarda `rng_intacto` cambiaba el consumo del
   generador en 5 de 6552 casos reales (por ejemplo F_577 de grado 80 en frío, o cuerpos
   con v2(q−1) = 6), y las pruebas no lo cubrían. Se corrigió en una segunda ronda con
   una prueba que lo habría detectado.
3. **Integración y medición en serie.** Las cuatro ramas se fusionaron sin conflictos.
   La suite intercalada detectó dos efectos de la integración, corregidos aquí (abajo).

## Detalle

### 1. Newton con producto medio y transformadas precalculadas

- `inverse_series` hace cada paso de Newton con una convolución cíclica de tamaño
  2^⌈log2 n⌉ y reutiliza la transformada de x_k: 5 transformadas en lugar de 6 que
  podían ser del doble de tamaño.
- El resto r = f − g·q se calcula módulo x^N − 1 con N = 2^⌈log2 deg g⌉.
- `Reducer` guarda la transformada de m (plegado) y las de rev(m)⁻¹, con cocientes de
  Shoup, en un `shared_ptr` perezoso que comparten sus copias. `reduce` hace 4
  transformadas en lugar de 6 de tamaño ~2n. El inverso se trunca para no doblar la NTT
  con n = 2^j + 1 (reserva del revisor de la investigación).
- `Composer` transforma h y el paso gigante una sola vez.
- **Generador aleatorio.** Todo esto solo actúa con la caché FFT ya llena
  (`ntt_sin_rng`), la misma condición que `mul_no_rng`.
- **Límites.** Solo Fp64. La división desde Python sigue por el camino frío (Kronecker),
  porque el binding crea cuerpos nuevos: falta la caché de raíz determinista
  (propuesta 3d de la investigación). Los umbrales de Newton no se han recalibrado.
  `Reducer` no es seguro entre hilos si sus copias se usan a la vez (hoy no ocurre).

### 5. F_p con p pequeño y Kronecker de F_q

- **Acumulador de 64 bits** en `mul_school`, `divmod_school` y `Composer::combina` de
  `Poly<Fp64>` cuando `terminos·(p−1)^2 < 2^64` (`Fp64::acc64_cabe`, con la cota en un
  comentario).
- **Producto de 64 bits sin vigilar el acarreo** (`acc_fma_peq`) con p < 2^32. Cota:
  harían falta más de 2^64 términos para desbordar 2^128.
- **Barrett de 64 bits** (`detail::Barrett64`) en `reduce128` y en el desempaquetado de
  Kronecker, en lugar de `u128 % p`.
- **Kronecker**: en Fp64 se empaqueta palabra a palabra (`poner`) y, con huecos de hasta
  64 bits, se lee cada hueco de una vez (`tomar`). En Fq64 los huecos van directos al
  acumulador sin reducir.

### 6. DDF sin el `pow_mod` repetido

- La sucesión de Frobenius cuenta `gastado_` (los `pow_mod` hechos sumando todos los
  módulos). Tras un cambio de módulo compone de inmediato si la regla de alquilar o
  comprar con ese gasto lo justifica, y solo si un `pow_mod` cuesta más del doble que un
  paso compuesto. Con ese margen no cambia nada en F_2 a F_5.
- `rng_intacto()` impide adelantar la composición cuando eso podría pedir tablas FFT
  antes que el original. Es imprescindible: ver la revisión.

### 8. Binding

- Nuevos `type_caster` (`EntInt`, `SalInt`, `ListaInt`, `ListaSal`). El camino rápido
  solo admite `int` exactos que caben en 64 bits con signo (o `list`/`tuple` exactas de
  ellos), detectados con `PyLong_AsLongLongAndOverflow` (API estable, no
  `PyUnstable_Long_IsCompact`). Todo lo demás (enteros grandes, `bool`, subclases,
  `numpy.int64`, `__index__`, `str`, `None`...) va por el camino anterior, así que las
  excepciones y sus mensajes no cambian.
- `python-bindings/tests/test_fase6_conversiones.py`; cobertura de pytest al 100 %.

## Dos correcciones de la integración

Medidas con la suite intercalada (`bench/results/fase6/`):

1. **La comprobación `p < 2^32` dentro de `Fp64::acc_fma`** (propuesta 5) se evaluaba en
   cada término y frenaba los bucles de p grande: Karatsuba de grado 4096 y `divmod` en
   F_(2^61−1) hasta 0,84-0,93x. Ahora `acc_fma` es la de siempre y la variante rápida se
   elige una vez por bucle (`detail::con_acc_fma`, en `field.hpp`) en `mul_school` y
   `combina`. En `FqT` se mantiene la comprobación por término (`acc_fma_rama`), porque
   ahí el despacho por producto medía un 3-8 % peor.
2. **`divmod_school` con la lambda** de `con_acc_fma` seguía a 0,85x en F_(2^61−1); con
   el bucle directo, 0,95-0,96x. Con p pequeño ese bucle no se usa (actúa el acumulador
   de 64 bits).

Se probó también volver al `reduce128` anterior: empeora F_(2^61−1) (`divmod` 0,83x,
`mul_school` 0,86x), así que se queda el de Barrett.

## Resultados

Suite intercalada, base `2b36cde` frente a la fase 6 (columna `f6-desd3`, el mismo
código que el commit final). Media geométrica por grupo:

| Grupo | Speedup | Grupo | Speedup |
|---|---:|---|---:|
| `fp.p7` | 1,54x | `fq.q3e30` | 1,40x |
| `fp.p101` | 1,49x | `fq.build` | 1,34x |
| `fp.ntt` | 1,14x | `fq.q2e8` | 1,19x |
| `fp.m127` | 1,07x | `fq.q7e4` | 1,14x |
| `fp.m61` | 1,01x | `fq.m61e2` | 1,00x |
| `fp.p256` | 1,00x | `fast` | 1,07x |
| **global (128)** | **1,17x** | | |

Mayores ganancias: `mul_school` en F_7/F_101 (3,2-3,3x), `pow_mod` en F_7 (2,2x) y F_101
(1,9x), `factor` en F_7 (1,8x), F_998244353 (1,8x), F_101 (1,7x) y F_(3^30) (1,45x).

**Regresiones.** La única estable es `fp.m61.divmod.d2000_1000`, a **0,95-0,96x** (8
rondas). Otros escenarios bajo 0,95x en la tanda completa (`fp.m61.mul_auto`,
`fp.p256.mul_fft.d4096`, `fq.q3e30.elem_inv`, `fq.m61e2.mul_karatsuba`) dieron
0,99-1,07x al repetirlos con 8 rondas: era ruido de la batería.

## Pruebas

- Nuevos: `test_fase6_newton` (30), `test_fase6_kronecker` (63), `test_fase6_ddf` (28),
  `test_fase6_conversiones.py`.
- `build.py all`: 12 tests C++, oráculo y 759 pytest, «TODO CORRECTO». Los avisos del
  oráculo (`FPX_FACTOR`, `FQX_FACTOR`) son de la referencia Python original y ya salían
  antes.
- No se ha vuelto a medir la cobertura C++ (`build.py cov-cpp`).

## Pendiente

- Repetir la suite enchufado y la comparativa con FLINT (`compare_libs.py`).
- Caché de raíz determinista para que `divmod` desde Python use la NTT (propuesta 3d).
- Recalibrar `prefer_newton` y `medio::kMin` con el Newton nuevo.
- Propuestas 2, 3, 4 y 7 de la investigación.
