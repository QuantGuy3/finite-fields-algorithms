# Equivalencia de SHAs tras la reescritura del historial

El 10 de noviembre de 2026 se reescribió el historial a la estructura `main` / `develop` /
`feature/*` (ver [flujo-git.md](flujo-git.md)). Cada commit conserva su contenido, autor y
fecha, pero cambia de SHA. Los informes de `docs/optimizacion/` y algunos mensajes citan SHAs
antiguos; esta tabla da el commit nuevo de cada uno de ellos. La copia completa del historial
anterior está en el bundle de respaldo local, fuera del repositorio.

| SHA antiguo | SHA nuevo | Commit |
|---|---|---|
| `67331a2` | `73b7663` | feat(binding): subir el binding pybind11 y el oráculo Python |
| `424ad1d` | `00b24dc` | wip(fase2): borrador del paquete Python alcp sin tests |
| `5f49217` | `3f00476` | docs: estado tras la fase 2 y módulo por defecto determinista en el diseño |
| `489728b` | `e910d63` | perf(cuerpos): usar la API C de GMP en sitio en Fp y Fq |
| `d4580af` | `e6ae5fc` | perf(fft): tablas de acceso directo, bit-reversal precalculado y twiddles contiguos |
| `8fcc151` | `c3b7c61` | perf(python): cachear en __slots__ los métodos nativos ligados del cuerpo |
| `0de8cf9` | `1eca6a1` | perf(poly): reducción perezosa en mul_school y divmod con ganchos del cuerpo |
| `f68a18b` | `04ad088` | perf(fft): cuadrado con una sola transformada directa |
| `dcb905c` | `deba46a` | perf(poly): suma y resta sin relleno de ceros ni copias intermedias |
| `abc2f4d` | `26bcd07` | perf(python): camino rápido para operandos int en + - * de elementos |
| `27873ea` | `21ff29a` | perf(python): el analizador construye x^k como monomio en lugar de potenciar |
| `9bfc1e7` | `fb8339e` | perf(poly): Karatsuba sobre vistas, sin slice/shift ni Poly temporales |
| `29d5628` | `28db98e` | perf(python): caché de cadenas ya analizadas en anillos y cuerpos |
| `0ab679a` | `c4f3c76` | perf(python): _pretty.polinomio formatea cada coeficiente distinto una sola vez |
| `a3399d9` | `cfa3d8d` | perf(python): caché por cuerpo de los textos de elementos de una extensión |
| `1ad8c02` | `b6163e4` | perf(fft): NTT en palabras de máquina para F_p con p < 2^62 |
| `cab2c2a` | `5236244` | perf(python): pickle de elementos más directo con __reduce_ex__ y reconstrucción rápida |
| `a833bcb` | `f20e73d` | perf(python): caminos rápidos int/elemento en la conversión de operandos y en f(v) |
| `65a5f1c` | `6a9dd41` | perf(cuerpos): backend de 64 bits Fp64 (Barrett) y extensión genérica FqT<Base> |
| `5c444c0` | `c00d04f` | perf(bench): selección de backend --backend auto\|gmp en bench_suite |
| `e683a88` | `fce226f` | perf(python): conversión de listas de coeficientes con caminos rápidos int/elemento |
| `70e3453` | `9ccfd4f` | perf(fft): NTT multi-limb mpn con Montgomery para F_p con p >= 2^62 |
| `3660189` | `843e589` | perf(poly): umbrales de despacho Auto por backend y por operando menor |
| `06e46b5` | `99c7521` | perf(fast): dft/idft explícitas en palabras de máquina para F_p con p < 2^62 |
| `e8e07e6` | `3e0fb5e` | perf(poly): paralelizar los niveles superiores de Karatsuba con std::async |
| `06d922a` | `00c08b8` | perf(python): caminos rápidos para // % divmod entre polinomios del mismo anillo |
| `4fb4f1a` | `68fd40e` | perf(binding): elegir automáticamente el backend de 64 bits en Fp y Fq |
| `b0a9fe7` | `b892fce` | perf(fft): paralelizar la NTT (mariposas genéricas y transformadas directas) |
| `0438b74` | `932062c` | perf(python): crear el resultado en línea en los caminos rápidos de + - * / |
| `c70bed5` | `8f70927` | perf(python): _pretty.polinomio escribe las potencias en línea |
| `2297387` | `48296ab` | perf(fp64): Barrett sin desplazamientos variables de 128 bits |
| `8f1d896` | `a8d4d80` | perf(fp64): potencia con intermedios en forma de Montgomery |
| `f30bafd` | `f01bfb4` | perf(integración): NTT en palabras para Fp64, umbrales por backend y búfer de Fq seguro con hilos |
| `5ce76db` | `6de3395` | revert(python): descartar 0438b74 (resultado en línea en + - * /) |
| `87f4a27` | `be5ebd0` | revert(python): descartar 06d922a (caminos rápidos para // % divmod) |
| `5d0556f` | `394ae29` | revert(python): descartar a833bcb (caminos rápidos en conversión y en f(v)) |
| `32d6922` | `8074fb8` | build(cov): contadores de cobertura atómicos y fft sin instanciar el camino GMP para Fp64 |
| `7880a88` | `493c8cb` | perf(fase4): nivel B con Kronecker, half-GCD y Frobenius por composición |
| `8a5a28e` | `f56c5f9` | perf(fase5): NTT perezosa, Euclides en sitio y DDF de Kaltofen-Shoup |
| `2b36cde` | `c1f98a1` | fix(fase5): gcd_ext en Fp64 termina con el Euclides extendido en sitio |
| `21e9fe6` | `f28dd67` | merge(feature/fase6-investigacion): investigación medida de la fase 6 |
| `80b983d` | `7475486` | docs(fase6): informe de la fase 6 y resultados medidos |
| `04c1909` | `5b47eaa` | build(cov): gcov en formato de texto para evitar cuentas truncadas |
| `4200f5b` | `8c724e8` | perf(fase8): Reducer por tabla de x^(n+i) mod m en Fp64 sin NTT |
