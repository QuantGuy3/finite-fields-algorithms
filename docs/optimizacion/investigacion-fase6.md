# Investigación: siguientes optimizaciones (propuesta de fase 6)

Fecha: 7 de noviembre de 2026. Punto de partida: `mejoras-cpp-pybind` en `8a5a28e` (fase 5).
Rama de la investigación: `investigacion-fase6`. **No cambia el código de la librería.**

## Método

Seis temas, investigados con un workflow de 18 agentes en tres fases:

1. **Prototipos.** Un agente por tema, todos en paralelo. Cada uno escribió perfiles y
   prototipos fuera de la librería, en `bench/experimentos/fase6/<tema>/`, y comprobó que
   el resultado coincidía exactamente con el de la librería.
2. **Medición limpia.** Los binarios se volvieron a ejecutar uno a uno, con el portátil
   **enchufado** y sin otros procesos midiendo. Salidas en `bench/results/fase6/<tema>.txt`.
3. **Verificación adversaria.** Un revisor por tema intentó refutar las conclusiones:
   comparación justa, comprobación de corrección, generador aleatorio y estimaciones.
   Ninguna conclusión principal quedó refutada. Los seis temas salieron «con reservas» y
   las correcciones se incorporan abajo.

Las cifras son de la medición limpia. Todo es Zen 3 (Ryzen 7 7730U), Windows, MinGW
g++ 15.2, `-O3` portable y un hilo. FLINT es python-flint 0.9.0, medido desde Python. El
detalle completo, con las reservas de cada revisor, está en
`bench/results/fase6/workflow_resultados.json`.

**Medido** significa prototipo ejecutado frente a la librería, con salida idéntica.
**Estimado** significa extrapolación de piezas medidas.

## Dos fallos encontrados en el código actual

1. **Regresión de `gcd_ext` en Fp64 desde el umbral de half-GCD (fase 5).** Medido:
   `P::gcd_ext` de grado 2048 tarda 45,5 ms y el Euclides clásico 17,5 ms, es decir, 2,6x
   más lento. En grado 4096: 171 ms frente a 68 ms.
   - **Causa.** Tras el primer `hgcd`, `r0` queda por debajo del umbral y el bucle sigue
     con `paso` sobre `Poly`, que no es en sitio.
   - **Corrección.** Seguir con `hgcd` hasta el caso base (umbral solo de entrada), o
     terminar la cola con el Euclides en sitio arrastrando los cofactores de la matriz.
   - **Ganancia medida con el prototipo V3e:** 7,1x en grado 2048 y 11x en grado 4096.
   - **Prioridad máxima.** El resultado es correcto; solo es lento.
   - **Corregido** en mejoras-cpp-pybind (commit `2b36cde`): cola en sitio; ver [fase5.md](fase5.md).
2. **La DDF de la fase 4 rehace un `pow_mod` tras cada factor encontrado.**
   - **Causa.** `cambiar_modulo` pone `pasos_ = 0`, y `quiere_componer` exige al menos un
     `pow_mod` antes de volver a componer.
   - **Efecto.** Afecta a grados por debajo de 40, donde no actúa Kaltofen–Shoup.
   - **Ganancia medida:** componer de inmediato hace la DDF 1,9-2,3x más rápida en
     F_(2^127−1) de grado 20, con la misma salida.
   - **Matiz del revisor.** Componer solo si compensa, aplicando la regla de alquilar o
     comprar con lo ya gastado, y medirlo también con q pequeño.

## Resultados por tema

### 1. `fpn`: F_p de tamaño fijo para p de 65 a 128 bits

- **Perfil** de `factor()` de grado 20 en F_(2^127−1): el 92 % es `pow_mod` de la DDF.
  Hay unas 195 000 reservas de GMP por llamada. Solo con un asignador de listas libres
  `factor()` ya va 1,7x más rápido, así que cerca del 40 % del tiempo es `malloc`.
- **Prototipo `FpN2`.** Dos palabras, Montgomery con doble REDC sobre el representante
  canónico, y acumulador perezoso de 5 palabras. Cumple `FiniteField` y
  `LazyAccumulation`, así que `Poly<FpN2>` ejecuta los algoritmos de la librería sin
  cambios.

  | Medido (p = 2^127 − 1) | Fp (mpz) | FpN2 | Ganancia |
  |---|---:|---:|---:|
  | producto de elementos | 66,7 ns | 11,0 ns | 6,1x |
  | producto 19x19 / resto 38/20 / gcd de grado 20 | 10,9 / 13,5 / 61,7 µs | 1,76 / 2,00 / 19,2 µs | 6,2x / 6,8x / 3,2x |
  | `pow_mod` de grado 20, exponente q | 6,30 ms | 0,98 ms | 6,4x |
  | `factor()`, 16 entradas de grado 20 | 15,8 ms | 2,76 ms | 5,7x |
  | Escenario `fp127.factor.d20` (FLINT 2,59 ms) | 18,1 ms (7,0x de FLINT) | 3,08 ms (**1,19x**) | 5,9x |

- **Estimado:** junto con la corrección de la DDF, unas 0,5-0,6x de FLINT. No se ha
  medido de forma conjunta.
- **Reservas:**
  - No se probaron primos de más de 127 bits; hay ramas de acarreo sin ejecutar.
  - Con p de 63-64 bits, `inv` necesita otro código.
  - La FFT con `FpN2` y su generador aleatorio no se han probado, porque 2^127 − 1 no
    tiene NTT.
  - En Linux con glibc las reservas son más baratas, así que la ganancia puede ser menor.
- **Esfuerzo:** medio, de 2 a 4 días. Un `fpn.hpp` de unas 500 líneas que siga el
  modelo de `fp64.hpp`, más el despacho en el binding.

### 2. `hgcd`: half-GCD en Fp64

- **Perfil** de `hgcd` en F_998244353, grado 1024-4096:
  - caso base en sitio, 24-35 %;
  - `aplicar`, 44-52 %;
  - `mul_mat`, 16-20 %;
  - fuera de `hgcd`, el `gcd` hace además un `aplicar` final de tamaño completo.
- **Prototipo V3e.** Mantiene los mismos cocientes, y el mcd y Bézout son idénticos.
  Combina tres cambios:
  - `hgcd` devuelve los restos y aplica la matriz solo a la parte baja;
  - la matriz 2x2 se multiplica con transformadas NTT reutilizadas;
  - la **envoltura cíclica** de tamaño 2^(k−1) se usa cuando la longitud se pasa poco de
    una potencia de 2.

  | Medido (F_998244353) | Antes | V3e | FLINT |
  |---|---:|---:|---:|
  | `hgcd` n = 1024 / 2048 / 4096 | — | 1,91x / 2,02x / 2,00x | — |
  | `mul_mat` n = 1024 | 119 µs | 51 µs (2,3x) | — |
  | `gcd` n = 2048 (clásico, umbral 512) | 4944 µs | 3121 µs (1,58x) | 2789 µs |
  | `gcd` n = 4096 / 8192 | 15,7 / 45,3 ms | 7,8 / 18,9 ms (2,0x / 2,4x) | 8,2 / 22,2 ms |
  | `gcd` n = 1000 | 1243 µs | 1129 µs (1,10x, en el ruido) | 644 µs |
  | `gcd_ext` n = 1000 | 4384 µs | 2495 µs (1,76x) | 1733 µs |

- **Conclusión:** con V3e, half-GCD iguala o supera a FLINT desde grado 4096. Strassen–
  Winograd 2x2 es 1,2x **más lento** y se descarta.
- **El MCD de grado 1000 sigue a 1,9x de FLINT.** Ese hueco está en el Euclides clásico
  en sitio, que es 1,6-1,9x más lento que el de FLINT, no en half-GCD.
- **Reservas del revisor:**
  - La cifra de `gcd_ext` es la de «`hgcd` hasta el caso base», no la de la cola en sitio.
  - La envoltura acumula en u128 sin reducir: con primos NTT cercanos a 2^62 y
    coeficientes p − 1 **desborda**. Hay que limitarla a p < 2^32 o reducir cada ~8
    términos.
  - El nuevo umbral de `gcd` debe exigir que la caché FFT ya tenga tablas.
- **Esfuerzo:** medio, unas 300 líneas.

### 3. `newton`: división de Newton, `Reducer` y `pow_mod`

- **La división desde Python va por el camino frío.** El binding crea un `Fp64` nuevo y
  `divmod` solo usa `mul_no_rng`, que con la caché vacía elige Kronecker. Medido en C++:
  1,33 ms en frío frente a 0,78 ms en caliente (1,7x), sin tocar el algoritmo. FLINT tarda
  0,70 ms.
- **Producto medio en Newton y resto módulo x^N − 1.** Medido: `divmod` 4000/2000 pasa de
  778 a 386-403 µs, 1,94x frente a caliente y 3,3x frente a frío, con lo que quedaría a
  0,55x de FLINT. La inversión de la serie gana 2,6x.
  - **Reserva:** 4000/2000 es casi el mejor caso de potencias de 2. Con otros tamaños se
    espera de 1,2x a 2,6x, y hay que medir uno desfavorable.
- **`Reducer` con las transformadas del módulo y de rev(m)⁻¹ precalculadas** (el `FFTRep`
  de NTL). Medido:
  - `reduce` va 2,0x más rápido;
  - `pow_mod` va 1,6x más rápido: de 2,73 a 1,71 ms en grado 300 (FLINT 1,54 ms) y de
    6,0 a 3,75 ms en grado 1000 (FLINT 10,1 ms);
  - `Composer` (Brent–Kung) va 1,7x más rápido.
  - **Reserva:** con n = 2^j + 1 conviene truncar el inverso para no doblar la NTT.
- **Ya estaba hecho:** el cuadrado usa una sola transformada directa (0,68x del producto).
- **Propuesta d), sin medir:** tablas FFT con una raíz **determinista** en una caché
  aparte, para `mul_no_rng`. No cambia el generador aleatorio, porque la raíz no influye
  en el producto (0 diferencias en 8 pruebas). Falta implementar la búsqueda de la raíz,
  demostrar que termina y medir su efecto en el MCD y en Euclides en frío.
- **Estimado:** `factor` de grado 300 e `is_irreducible` mejorarían, pero el revisor
  señaló que el recuento de productos usado no cuadra con el tiempo total. Hay que
  instrumentarlo antes de dar una cifra.

### 4. `ntt`: estructura de la NTT y AVX2

Con el portátil enchufado, la mariposa actual cuesta **1,5-2,0 ns**, no 4-5: aquella cifra
era con batería.

**Perfil.** Las mariposas son el 75-85 % del tiempo. La permutación de bits pesa poco con
n = 2048, pero ronda el 18 % con n = 32 768 por fallos de caché.

| Variante (convolución completa; n = 256 / 2048 / 32 768) | Medido |
|---|---:|
| a) Directa DIF sin permutación + inversa DIT | 1,07x / 1,05x / 1,36x |
| c) Producto punto a punto de Montgomery fusionado con n⁻¹ | ~1,01x |
| b) Radix 4, etapa len = 4 especializada y Montgomery (portable) | **1,57x / 1,47x / 1,79x** |
| d) AVX2, 8 carriles de 32 bits, p < 2^30, despacho en tiempo de ejecución | **7,8x / 8,3x / 9,9x** |
| Salto de potencia de 2 (longitud 2049 frente a 2048), en la librería | 2,1x más caro |

- **Frente a FLINT, desde Python** (esta rueda de python-flint probablemente no incluye
  `fft_small`): alcp ya va a 0,89x en grado 1000 y a 0,99x en grado 10 000.
- **Reservas del revisor:**
  - No hay que cambiar `word::transform`, porque `fft::transform` es API pública y espera
    el orden natural. Hay que añadir rutinas propias de la convolución y conservar
    `t.swaps`.
  - Falta prototipar el camino del cuadrado.
  - AVX2 necesita un respaldo para n < 16.
  - En MinGW, GCC no garantiza la alineación de 32 bytes de la pila (bug 54412), así que
    hay que probarlo con varias compilaciones.
- **Esfuerzo:**
  - escalar: 1-2 días;
  - AVX2: 3-5 días, con tablas de 32 bits derivadas de los mismos twiddles, sin generador
    aleatorio.

### 5. `factor`: q pequeño y extensiones

**Perfil.** En F_(3^30) y F_256 casi todo es aritmética de elementos de F_q: por ejemplo,
5 891 `acc_fma` y 236 `inv` de 5 µs cada uno por `factor()` de grado 5. En `mul_extension`,
el `mpn_mul` es solo el 35-41 %; el resto es desempaquetar y reducir módulo g. En F_101
pesan el desempaquetado `u128 % p` de Kronecker y el acumulador de 128 bits.

| Prototipo (medido; salida idéntica) | Antes → después | FLINT |
|---|---|---:|
| **FqTab** (F_(2^e), e ≤ 16, elemento en una palabra con tablas exp/log): `factor` de F_256 de grado 30 | 6,7 ms → 0,18 ms (**37x**) | 5,3 ms |
| FqTab: producto de grado 1000 en F_256 | 1080 → 297 µs (3,6x) | 665 µs |
| F_101: `mul_primo` sin `u128 % p` + acumulador de 64 bits: `factor` de grado 40 / 100 | 211 → 136 µs (1,56x) / 1383 → 802 µs (1,72x) | 78 / 475 µs |
| F_101: producto de grado 16-1000 | 1,45-1,99x | — |
| **FqSwar** (F_(3^30) en 64 bits, 4 huecos de 16 bits): `factor` de grado 5 | 5,3 ms → 2,06 ms (2,5x) | 3,35 ms |
| Kronecker en F_q con desempaquetado directo: producto de grado 200 en F_(3^30) | 985 → 738 µs (1,34x); con aritmética de 64 bits, 584 µs (1,7x) | 543 µs |
| KS2 de Harvey (evaluar en 2^b y 2^−b) | ≤ 1,07x: **se descarta** | — |
| `pow_mod` con ventana deslizante | 1,17x con exponente de 48 bits; 0,64x con 7 bits: solo para exponentes grandes | — |

**Reservas del revisor (importantes):**

- El medidor reiniciaba el índice de entrada en cada tanda. En los casos lentos solo
  medía las primeras entradas, así que las cifras de F_(3^30) de grado 5 y F_256 de grado
  30 hay que repetirlas con el medidor corregido. Por recuento de operaciones, FqSwar
  rondaría 2,9x.
- La comprobación del generador se hizo con la caché ya caliente. Falta una prueba en frío
  y en un cuerpo con FFT.
- El producto de grado 200 en F_(3^30) queda a 1,08x de FLINT, no por debajo, salvo con
  `-march=native`.

### 6. `python`: coste fijo de las llamadas

- **El binding no explica la distancia a FLINT en el producto de grado 100.** alcp_native
  tarda 6,4 µs, C++ puro 5,85 µs y FLINT 3,7 µs. Eliminar todo pybind11 ganaría como mucho
  1,08x: el hueco está en el núcleo C++ a tamaños pequeños.
- **Donde sí se gana: conversión directa entre `int` y `u64`** en Fp64 (medido):
  - `K.poly(lista)`: 7,6x en grado 100 y 13,9x en grado 1000;
  - `coeffs()`: 4,4x;
  - operaciones de elementos: 1,8x (335 → 184 ns).
  Hoy cada coeficiente pasa por `std::vector<BigInt>` y un `mpz`.
- **Con tipos CPython propios (`METH_FASTCALL`):** elementos de 335 a 28 ns, al nivel de
  FLINT. El producto de grado 10 bajaría de 0,8 a 0,25 µs.
  - Cuesta de 3 a 5 días.
  - Solo mueve una fila de la comparativa, un 6-8 %.
- **Reservas del revisor:**
  - Hay que conservar los mensajes `TypeError` de `str` y `None`.
  - `PyUnstable_Long_IsCompact` requiere CPython ≥ 3.12: hace falta otro camino para
    3.11 o anteriores.
  - Falta probar `numpy.int64`.

## Prioridades recomendadas

| # | Cambio | Evidencia | Ganancia | Esfuerzo |
|---|---|---|---|---|
| 0 | **Corregir la regresión de `gcd_ext`** (seguir con `hgcd` hasta el caso base) | Medido | 2,6x frente a hoy en grado ≥ 2048 (7-11x con V3e) | Bajo |
| 1 | Producto medio en Newton + `Reducer`/`Composer` con transformadas precalculadas | Medido | `divmod` 1,9x (3,3x en frío); `pow_mod` y composición 1,6-1,7x | Medio |
| 2 | Half-GCD V3e (restos, NTT reutilizada, envoltura limitada a p < 2^32) | Medido | `gcd` 2,0-2,4x desde grado 2048; `gcd_ext` 1,8x en grado 1000 | Medio |
| 3 | `FpN<2>` para p de 65-128 bits | Medido | 5,7-6,4x con p grande; `factor` de grado 20 de 7,0x a 1,2x de FLINT | Medio |
| 4 | NTT radix 4 (portable) y luego AVX2 con despacho | Medido | 1,5-1,8x; con AVX2, 8-10x en la convolución | Medio-alto |
| 5 | F_101: `mul_primo` sin `%` y acumulador de 64 bits; Kronecker de F_q con desempaquetado directo | Medido | 1,3-1,7x | Bajo |
| 6 | DDF de la fase 4 sin el `pow_mod` repetido | Medido | DDF 1,9-2,3x por debajo de grado 40 | Bajo |
| 7 | FqTab para F_(2^e) y FqSwar para F_(p^e) con p pequeño | Medido, pendiente de repetir | F_256: 3,6-37x; F_(3^30): ~2,5-2,9x | Medio-alto |
| 8 | Conversión directa `int`/`u64` en el binding | Medido | Conversiones 4-14x; elementos 1,8x | Bajo-medio |
| — | KS2, Strassen–Winograd 2x2, Montgomery fusionado, ventana con exponentes pequeños | Medido | ≤ 1,07x o pérdidas | Descartados |

Todas las propuestas mantienen los resultados idénticos y el consumo del generador
aleatorio. Las que tocan productos «sin generador» solo pueden usar la NTT si la caché ya
tiene tablas, salvo que se adopte la caché determinista de la propuesta d).

## Limitaciones

- **Prototipos fuera de la librería.** Integrados puede haber costes de conversión y
  reservas, y los umbrales tendrán que recalibrarse.
- **Una sola máquina** (Zen 3, Windows y MinGW). Las ganancias que dependen de `malloc`
  (FpN2) o de AVX2 pueden variar en otras plataformas.
- **FLINT se midió desde Python**, con el coste de su binding. Probablemente es una rueda
  sin `fft_small`.
- **Pendiente de repetir.** Las cifras de FqSwar y FqTab en `factor` deben repetirse con
  el medidor corregido, como se indica en el tema 5.
