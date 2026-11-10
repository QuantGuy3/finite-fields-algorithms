# Fase 8: Reducer por tabla, Euclides con AVX2, FqSwar con NTT y cabos sueltos

Fecha: 9-10 de noviembre de 2026. Rama `fase8`, sobre el final de la fase 7 (`04c1909`). Ataca
las diferencias con FLINT que quedaban tras la fase 7 ([fase7.md](fase7.md), «Comparativa
final» y «Perfil de factor»). Implementa las propuestas 1, 2, 4 (en parte), 5 y 6. La 3, una
DDF con matriz de Frobenius, se descartó por decisión del usuario.

## Resumen

| # | Cambio | Dónde | Medido (en serie) |
|---|---|---|---|
| 1 | **Reducer por tabla.** Para p < 2^16 y 8 ≤ n ≤ 384, con AVX2, el resto de a (grado ≤ 2n−1) módulo m fijo se calcula como bajos(a) + Σ a_{n+i}·(x^(n+i) mod m). Las filas están precalculadas en 32 bits y la suma es exacta en 8 carriles de 32 bits o 4 de 64. La tabla se construye tras 1 + n/128 reducciones con el mismo módulo y se comparte entre copias con `std::call_once`. El `Composer` sin NTT usa la misma suma | `polynomial.hpp` (`namespace tabla`, `Reducer`, `Composer`) | `reduce` de grado 100 en F_101 unas 4,7x más rápido. F_101 1,27x y F_7 1,28x en la suite. `factor` grado 100 en F_101 pasa de 0,61x a 1,34x frente a FLINT |
| 2 | **Euclides clásico con AVX2** (p < 2^31). Las filas de la división (`fila2` para cociente de grado 1 y `axpy`) usan productos de Shoup de 32 bits con `_mm256_mul_epu32`. Mismo despacho que la NTT; el camino escalar queda intacto | `polynomial.hpp` (`namespace euclides::avx2`) | `gcd` 1,41x y `gcd_ext` 1,35x en la suite. `gcd` grado 1000 en F_998244353 pasa de 0,65x a 2,22x frente a FLINT |
| 5 | **FqSwar.** Los productos grandes convolucionan los huecos de Kronecker con la NTT módulo 998244353, mediante un `Fp64` auxiliar con tablas sin generador. Además: huecos de 16 bits usados directamente como acumulador SWAR, reducción módulo g con tabla t^k mod g, y Kronecker desde 8 coeficientes | `fq_swar.hpp`, `kronecker.hpp` | Producto grado 200 en F_(3^30): de 0,91x a 2,45x frente a FLINT |
| 4 | **Constantes de FqTab.** `par_karatsuba_min(K)` vale 512 con p impar y 1024 con p = 2. `kFrobPresupuesto` no cambia: no hubo diferencia clara | `polynomial.hpp`, `bench/calibra_fqtab.cpp` | F_(3^7) y F_(7^4), Karatsuba de grado 512: de 1,37x a 1,74x (medido con carga) |
| 6 | **Cabos sueltos de la revisión de la fase 7.** `fft::convolve` rechaza coeficientes no canónicos antes de pedir tablas. Las tablas deterministas se liberan cuando la caché normal ya tiene tablas. Las docstrings de `Fp` y `backend` explican cuándo el `repr` incluye `backend='gmp'`. Se limpiaron las ramas y los worktrees antiguos de agentes | `fft.hpp`, `bindings.cpp` | — |

**Suite C++** (128 escenarios, intercalada con el final de la fase 7, 3 rondas × 5
repeticiones): **1,076x** de media geométrica y **todos los checksums iguales**.

| Grupo / operación | Speedup |
|---|---|
| `fp.p7` / `fp.p101` | 1,28x / 1,27x |
| `fp.ntt` | 1,08x |
| `gcd` / `gcd_ext` | 1,41x / 1,35x |
| `pow_mod` / `is_irreducible` / `factor` | 1,33x / 1,19x / 1,09x |
| `fp.m61`, `fp.m127`, `fp.p256`, `fq.m61e2` | 0,99-1,01x (sin cambios) |

Tres escenarios de código que no se tocó salieron a 0,93x: `fp.p101.mul_karatsuba.d512`,
`fq.q2e8.elem_inv.n5000` y `fast.ntt.toeplitz_lower_inverse.n1024`. En las fases anteriores
casos así resultaron ser ruido al repetirlos, pero estos **no se han repetido**.

## Comparativa con FLINT

`compare_libs.py --libs alcp_native,flint` (python-flint 0.9.0), una pasada, con los checksums
iguales. La tabla muestra tiempo FLINT / tiempo alcp; mayor que 1 significa que alcp es más
rápido.

| Escenario | Final de la fase 7 | Fase 8 |
|---|---|---|
| `fp30.gcd.d1000` | 0,65x | **2,22x** |
| `fp101.factor.d100` | 0,61x | **1,34x** |
| `fp101.factor.d40` | 0,78x | 0,86x |
| `fp101.irred.d200` | 1,10x | **2,32x** |
| `f3e30.mul.d200` | 0,91x | **2,45x** |
| `fp30.irred.d400` | 3,02x | 3,48x |
| `fp30.divmod.d4000` | 1,82x | 1,74x |
| resto (producto, F_(2^127−1), F_256) | 1,3-50x | sin cambios significativos |

Solo `fp101.factor.d40` sigue por detrás de FLINT. Es un producto de factores conocidos de
grado 40, demasiado pequeño para que se construya la tabla del `Reducer`.

## Método

Se siguió el de las fases anteriores: un agente por propuesta en su worktree, una revisión
adversaria independiente por rama, integración en serie y medición sin otros procesos
pesados. Las tres revisiones **aprobaron**:

- **Propuesta 1:** 101.828 comprobaciones de las primitivas, incluido el peor caso de las
  cotas de 32/64 bits; diferencial byte a byte en 946 líneas y 5 configuraciones de
  compilación; ensamblador revisado. La revisión dejó cuatro observaciones, corregidas en la
  rama antes de integrar: sincronización entre hilos, versiones escalares en lugar de ramas
  que devolvían ceros sin AVX2, límite de grado del `Composer` y prueba de la frontera 32/64.
- **Propuesta 2:** 1.169.176 casos de las primitivas con centinelas de memoria; diferencial
  byte a byte con 13 primos.
- **Propuesta 5:** diferencial de 53.281 líneas byte a byte, todos los tamaños de 1 a 4000,
  bordes con coeficientes máximos (hasta 2^23 huecos) y prueba con hilos. La errata del
  comentario que señaló se corrigió en la integración.

**Pila y AVX2 en MinGW.** GCC en MinGW-w64 no alinea la pila a 32 bytes. Las revisiones
comprobaron en el ensamblador (de -O1 a -O3, con y sin `-march=native`) que no hay volcados
alineados de registros ymm a la pila: solo `vmovdqu`, que no exige alineación. Eso lo decide
el compilador; si se cambia de versión de GCC, conviene repetir la comprobación.

## Pruebas

- Nuevos tests: `test_fase8_reducer.cpp` (219 comprobaciones), `test_fase8_euclides.cpp`
  (43), `test_fase8_fq.cpp` (78) y `test_fase8_cabos.cpp` (18).
- `build.py all`: «TODO CORRECTO» (21 tests C++, oráculo y 837 pruebas de pytest con
  cobertura del 100 %).

## Pendiente

- Recalibrar Newton con 2^16 < p < 2^32 sin NTT (p. ej. 2^31 − 1).
- Umbrales del half-GCD tras el Euclides con AVX2. Según el agente de la propuesta 2, el
  clásico gana hasta unos 6144 con NTT (umbral actual 1792) y hasta unos 5120 en F_101
  (actual 768). Hay que confirmarlo con un barrido en serie.
- Repetir los tres escenarios que salieron a 0,93x.
- `fp101.factor.d40`: el producto de Kronecker (unos 3 µs en grado 100) pesa ya más que la
  reducción; un producto escolar vectorizado como el del `Reducer` podría bajarlo.
