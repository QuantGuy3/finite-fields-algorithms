# Comparativa de velocidad con otras librerías de cuerpos finitos

Rama `mejoras-cpp-pybind`, 4 de noviembre de 2026. Los datos están en
`bench/results/comparativa/comparativa.json` y la tabla completa, generada a partir de ellos,
en [comparativa-tabla.md](comparativa-tabla.md). Para reproducirla:

```bash
python bench/compare_libs.py --json bench/results/comparativa/comparativa.json
python bench/informe_comparativa.py bench/results/comparativa/comparativa.json docs/optimizacion/comparativa-tabla.md
```

> **Actualización (fase 4).** La tabla completa, [comparativa-tabla.md](comparativa-tabla.md), ya incluye la fase 4, la variante «alcp fase 3, 1 hilo» y dos escenarios de grado alto (16 en total). Con un hilo, la distancia media a FLINT pasa de 6,9x (fase 3) a 5,9x (nivel A) y a 2,7x (nivel B), medidos en una misma tanda. El análisis está en [fase4.md](fase4.md). El texto de abajo describe la comparativa de la fase 3.

## 1. Qué se compara

| Variante | Qué es |
|---|---|
| alcp original | El binding del código **original** (`baseline-original`, 67331a2), compilado con `-O3` |
| alcp final | El binding de la versión final (fase 3), con paralelismo |
| alcp final, 1 hilo | Lo mismo con `ALCP_THREADS=1`. Las demás librerías usan un solo hilo, así que **esta es la comparación justa** |
| alcp paquete | La capa Python de alto nivel (`alcp.GF`, `K["x"]`) sobre el binding final |
| python-flint 0.9 | Interfaz Python de [FLINT](https://flintlib.org) 3, la referencia en C del sector. Usa `nmod_poly` si p cabe en una palabra, `fmpz_mod_poly` para p grande y para el test de irreducibilidad (`nmod_poly` no lo ofrece) y `fq_default_poly` para extensiones |
| galois 0.4.11 | Librería sobre NumPy compilada con numba (JIT) |
| SymPy 1.14, Python puro | `SYMPY_GROUND_TYPES=python`: aritmética de polinomios en Python |
| SymPy 1.14 sobre FLINT | SymPy, cuando python-flint está instalado, delega los polinomios en él (`DUP_Flint`) |

No se incluyó **Nemo.jl**, aunque Julia 1.12 está instalado y el repositorio tiene
`compare_nemo.jl`: el Control de aplicaciones de Windows (Smart App Control) bloquea sus DLL
precompiladas. Es una configuración de seguridad del sistema y no se modificó.

## 2. Metodología

- **Mismas entradas.** Las entradas se generan con `random.Random(semilla fija por escenario)`,
  como coeficientes enteros. En las extensiones, cada coeficiente es el índice en base p del
  elemento, con el **mismo** polinomio de definición en todas las librerías: el de AES para
  F_256 y el módulo fijo de `bench_suite.cpp` para F_(3^30).
- **Resultados verificados.** Cada librería normaliza su resultado a una forma canónica
  (coeficientes como índices, factores mónicos ordenados con su exponente), y se compara
  un hash. **En los 14 escenarios, todas las librerías que terminaron dieron el mismo
  resultado.**
- **Solo se cronometra la operación.** La creación del cuerpo y la conversión de las
  entradas quedan fuera. Así se excluye, por ejemplo, la compilación JIT de galois (unos
  20 s la primera vez) y no se penaliza a ninguna librería por su formato de entrada.
- **Un proceso nuevo** para cada par (librería, escenario), con un calentamiento y
  repeticiones hasta reunir 5 muestras y 0,5 s (las operaciones de más de 20 s se miden una
  vez). Se informa del **mínimo**. El recolector de basura se desactiva durante la medida.
- **Máquina:** AMD Ryzen 7 7730U (8 núcleos/16 hilos), Windows 11, Python 3.15.1,
  enchufado y en reposo (no se medía nada más).

**Corrección aplicada.** galois guarda el resultado de `is_irreducible` en una
`functools.lru_cache`, así que tras el calentamiento solo se medía la consulta a la caché
(6 µs). El escenario se repitió llamando a la función sin caché (`__wrapped__`): 287 ms.

## 3. Resultados

Tiempo mínimo por operación. En negrita, el más rápido de cada fila.

| Operación | alcp original | alcp final, 1 hilo | python-flint | galois | SymPy puro |
|---|---:|---:|---:|---:|---:|
| Producto grado 100, F_998244353 | 561 µs | 6,3 µs | **3,6 µs** | 57 µs | 5,3 ms |
| Producto grado 1000, F_998244353 | 6,47 ms | 121 µs | **69 µs** | 2,47 ms | 233 ms |
| Producto grado 10 000, F_998244353 | 172 ms | 4,13 ms (2,97 ms con hilos) | **1,69 ms** | 252 ms | 10,0 s |
| `divmod` 4000 / 2000, F_998244353 | 431 ms | 3,03 ms | **0,69 ms** | 6,42 ms | 3,56 s |
| MCD grado 1000, F_998244353 | 124 ms | 4,27 ms | **0,58 ms** | 87 ms | 2,79 s |
| Factorizar grado 40 (con multiplicidades), F_101 | 8,35 ms | 181 µs | **71 µs** | 53 ms | 8,1 ms |
| Factorizar grado 100, F_101 | 480 ms | 6,66 ms | **0,45 ms** | 420 ms | 243 ms |
| Irreducibilidad grado 200, F_101 | 20,1 s | 207 ms | 7,7 ms (SymPy+FLINT: **2,5 ms**) | 287 ms | 921 ms |
| Producto grado 500, F_(2^127−1) | 13,5 ms | 3,06 ms (1,41 ms con hilos) | **0,33 ms** | 16,8 ms | 98 ms |
| Factorizar grado 20, F_(2^127−1) | 108 ms | 30,3 ms | **2,56 ms** | 2,28 s | 79 ms |
| Producto grado 1000, F_256 | 1,02 s | 15,2 ms (4,39 ms con hilos) | **0,77 ms** | 2,10 ms | n/s |
| Factorizar grado 30, F_256 | 892 ms | 13,7 ms | **8,44 ms** | 51 ms | n/s |
| Producto grado 200, F_(3^30) | 855 ms | 11,9 ms (3,30 ms con hilos) | **0,65 ms** | 6,85 s | n/s |
| Factorizar grado 5, F_(3^30) | 336 ms | 9,29 ms | **2,14 ms** | 3,19 s | n/s |

n/s: no soportado (SymPy no tiene polinomios sobre F_(p^n)).

### Resumen: cociente frente a python-flint (media geométrica; 1 = igual que FLINT)

| Librería | Media geométrica | Rango |
|---|---:|---|
| alcp original | **248x** más lento | 41x – 2625x |
| alcp final, 1 hilo | **5,96x** más lento | 1,6x – 27x |
| alcp final (con hilos) | 4,61x más lento | 1,6x – 28x |
| galois | 110x más lento | 2,8x – 10 567x |
| SymPy puro | 742x más lento | 31x – 5931x (10 escenarios) |
| SymPy sobre FLINT | 0,99x (igual) | 0,33x – 1,8x (9 escenarios) |

## 4. Lectura de los resultados

1. **La fase 3 cambió la posición de alcp.** El código original era unas 250 veces más
   lento que FLINT, e incluso más lento que SymPy en Python puro en 4 de los 10 escenarios
   comparables (las factorizaciones y la irreducibilidad). La versión final, con un hilo, está a unas 6 veces de FLINT y supera a
   galois y a SymPy puro en todos los escenarios salvo uno (el producto en F_256, donde
   galois usa tablas de logaritmos).
2. **FLINT sigue siendo más rápido en todo**, y la distancia es mayor justo donde FLINT usa
   **algoritmos asintóticamente mejores**. Igualarlos queda fuera de la regla del encargo,
   que era no cambiar los algoritmos:
   - MCD (7x): FLINT usa *half-GCD*, O(M(n) log n). alcp, el Euclides clásico, O(n²).
   - Irreducibilidad y factorización (12-27x): FLINT usa potencias de Frobenius con
     composición modular (Kaltofen–Shoup, Brent–Kung). alcp, Ben-Or y Cantor–Zassenhaus con
     `pow_mod` por cuadrados repetidos.
   - `divmod` (4,4x): FLINT usa división por inversión de Newton, O(M(n)). La división de
     alcp es la escolar, O(n²).
   - Producto (1,7-20x): misma familia de algoritmos, pero FLINT tiene años de ajuste fino
     (Kronecker substitution en extensiones, aritmética vectorizada, FFT de pequeños primos).
     Es donde alcp queda más cerca (1,7-2,4x en F_998244353).
3. **En las extensiones es donde más ganó alcp respecto al original** (36-72x en F_256 y
   F_(3^30)), gracias al backend de 64 bits y a la reducción perezosa. La mayor distancia
   restante con FLINT está en el producto (18-20x con un hilo): FLINT empaqueta el
   polinomio de F_q en un entero (Kronecker) y multiplica con una sola FFT, mientras que
   alcp aplica Karatsuba o la NTT elemento a elemento.
4. **La capa Python de alto nivel no cuesta casi nada** en operaciones de polinomios: «alcp
   paquete» queda dentro de un ±4 % del binding en los 14 escenarios.
5. **SymPy sobre FLINT** iguala a FLINT (es FLINT con una capa fina). En la irreducibilidad
   es incluso más rápido que python-flint directo: `DUP_Flint.is_irreducible` factoriza
   con `nmod_poly.factor()` y comprueba que salga un único factor (2,5 ms), mientras que el
   adaptador de python-flint usa `fmpz_mod_poly.is_irreducible` (7,7 ms), porque
   `nmod_poly` no ofrece ese test. Es decir, el mejor tiempo de FLINT para esta operación
   es 2,5 ms, y la distancia real de alcp en ella es de unas 80 veces. Con p = 2^127 − 1, SymPy 1.14 sobre python-flint 0.9
   **falla** al factorizar (`TypeError` en el constructor de `fmpz_mod_poly`): es un fallo
   de compatibilidad entre esas dos versiones.

## 5. Limitaciones

- Una sola máquina y una sola sesión. Los cocientes por escenario por debajo de un ~10 %
  no son significativos.
- Se mide el **mínimo**, que es el mejor estimador del coste sin interferencias, pero no
  refleja la variabilidad.
- Los escenarios son 14 y se eligieron para cubrir las operaciones del proyecto (producto,
  división, MCD, irreducibilidad y factorización) en cuerpos primos pequeños, de 30 bits y
  de 127 bits, y en dos extensiones. Otras operaciones o tamaños pueden dar cocientes
  distintos.
- El paralelismo de alcp solo actúa en productos grandes. La columna «1 hilo» es la
  comparación en igualdad de condiciones.
- Cada librería se usa tal como viene instalada. No se compiló FLINT con opciones
  específicas para esta CPU.

## 6. Conclusiones y recomendaciones

- **Conclusión:** tras la fase 3, alcp pasa de ser la opción más lenta (después de SymPy en
  Python puro) a ser la más rápida de las librerías que no usan FLINT, con resultados
  idénticos. Frente a FLINT queda a unas 6 veces de media con un hilo, y la mayor parte de
  esa distancia se debe a diferencias de **algoritmo**, no de implementación.
- **Para producción** sobre cuerpos finitos en Python, la opción más rápida es python-flint.
  alcp tiene sentido como librería didáctica y legible: su código sigue los algoritmos del
  curso y los resultados se comprueban contra las prácticas originales.
- **Si se quisiera acercarse a FLINT** (fuera del alcance de este encargo, porque supone
  cambiar algoritmos): half-GCD, división por inversión de Newton, Kronecker substitution
  en F_q y potencias de Frobenius por composición modular en la factorización.
