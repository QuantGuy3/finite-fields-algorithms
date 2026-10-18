// bindings.cpp — Binding pybind11 del núcleo C++ (alcp) para Python.
//
// Genera el módulo nativo `alcp_native` (python-bindings/alcp/alcp_native<EXT_SUFFIX>),
// según el contrato de docs/diseno.md, §5.  Es un superconjunto COMPATIBLE con la
// API original (Fp, Fq, PolyFp, PolyFq, is_prime, prime_power_factor, dft, idft,
// toeplitz_lower_vec, toeplitz_lower_inverse y sus métodos), de modo que los
// scripts antiguos (showcase.py, bench/bench_python.py) siguen funcionando.
//
// Representación de los elementos en Python
//   * F_p: enteros de Python.  Cualquier entero de entrada se reduce módulo p;
//     las salidas están normalizadas en [0, p).
//   * F_q = F_p[a]/(g): índices en base p.  El elemento
//     c_0 + c_1·a + ... + c_{n-1}·a^{n-1} se codifica como el entero
//     c_0 + c_1·p + ... + c_{n-1}·p^{n-1} en [0, q) (igual que
//     cuerpo_fq.elem_de_int en la práctica original).  Los índices de entrada se
//     reducen módulo q; las salidas están en [0, q).
//
// Conversión de enteros Python <-> mpz (type_caster<mpz_class>, más abajo)
//   Sin cadenas decimales (cuadráticas y limitadas a 4300 dígitos en CPython):
//   camino rápido para valores de 64 bits y, para el resto, la magnitud en bytes
//   (int.to_bytes / int.from_bytes) con mpz_import / mpz_export.  Se acepta
//   cualquier objeto con __index__ (bool, numpy.int64, ...); los float se
//   rechazan con TypeError.
//
// Excepciones del núcleo -> Python (traductor registrado en PYBIND11_MODULE)
//   alcp::division_by_zero -> ZeroDivisionError   (se comprueba ANTES que
//                                                   std::domain_error, su base)
//   std::invalid_argument, std::domain_error -> ValueError  (std::length_error
//                                   también, por el traductor por defecto)
//   std::overflow_error -> OverflowError
//
// GIL: NO se libera en ninguna llamada (no hay py::gil_scoped_release ni
// py::call_guard).  El núcleo tiene estado global mutable sin sincronizar: el
// generador aleatorio de GMP de bigint.hpp y, en cada cuerpo, la caché FFT y la
// raíz base (miembros `mutable` de Fp/Fq).  Dos hilos que lo usaran a la vez
// podrían corromperlo; mantener el GIL serializa todas las llamadas.
//
// Vida de los cuerpos: cada objeto Fp/Fq de Python se guarda en un
// std::shared_ptr (holder de pybind11) y cada polinomio comparte ese puntero,
// así que el cuerpo C++ al que apunta Poly<Field> vive al menos tanto como el
// polinomio.
#include "alcp/alcp.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <concepts>
#include <cstddef>
#include <exception>
#include <forward_list>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace py = pybind11;

using alcp::BigInt;
using alcp::Fp;
using alcp::Fq;

// ---------------------------------------------------------------------------
// Conversión de enteros Python <-> mpz
// ---------------------------------------------------------------------------
namespace alcp_py {

// Toma posesión de una referencia nueva de la API C; si es nula, propaga la
// excepción de Python pendiente.
inline py::object propio(PyObject* o) {
    if (o == nullptr) throw py::error_already_set();
    return py::reinterpret_steal<py::object>(o);
}

// int de Python (tipo exacto int) -> mpz.
inline void pylong_a_mpz(PyObject* o, mpz_ptr z) {
    int desborde = 0;
    // `o` es un int exacto: esta llamada no puede fallar, solo desbordar.
    const long long v = PyLong_AsLongLongAndOverflow(o, &desborde);
    if (desborde == 0) {
        // Camino rápido: el valor cabe en 64 bits con signo.
        if (std::in_range<long>(v)) {
            mpz_set_si(z, static_cast<long>(v));
        } else {
            // long de 32 bits (Windows): se importa la magnitud como una palabra de 64 bits.
            const unsigned long long mag = v < 0 ? 0ULL - static_cast<unsigned long long>(v)
                                                 : static_cast<unsigned long long>(v);
            mpz_import(z, 1, -1, sizeof mag, 0, 0, &mag);
            if (v < 0) mpz_neg(z, z);
        }
        return;
    }
    // Camino general (|v| >= 2^63): magnitud en bytes little-endian e importación
    // directa; coste lineal y sin límite de dígitos.
    py::object mag = desborde < 0 ? propio(PyNumber_Negative(o)) : py::reinterpret_borrow<py::object>(o);
    const auto nbytes = (mag.attr("bit_length")().cast<std::size_t>() + 7) / 8;
    py::object bytes = mag.attr("to_bytes")(nbytes, "little");
    mpz_import(z, nbytes, -1, 1, 0, 0, PyBytes_AS_STRING(bytes.ptr()));
    if (desborde < 0) mpz_neg(z, z);
}

// mpz -> int de Python (referencia nueva).
inline py::object mpz_a_pylong(mpz_srcptr z) {
    if (mpz_fits_slong_p(z)) return propio(PyLong_FromLong(mpz_get_si(z)));
    const std::size_t nbits = mpz_sizeinbase(z, 2);
    if (nbits <= 63) {
        // Cabe en long long (caso de Windows, donde long tiene 32 bits).
        unsigned long long mag = 0;
        mpz_export(&mag, nullptr, -1, sizeof mag, 0, 0, z);
        const long long v = static_cast<long long>(mag);
        return propio(PyLong_FromLongLong(mpz_sgn(z) < 0 ? -v : v));
    }
    // Camino general: se exporta la magnitud directamente a un objeto bytes y se
    // reconstruye con int.from_bytes.
    const std::size_t nbytes = (nbits + 7) / 8;
    py::object bytes = propio(PyBytes_FromStringAndSize(nullptr, static_cast<Py_ssize_t>(nbytes)));
    mpz_export(PyBytes_AS_STRING(bytes.ptr()), nullptr, -1, 1, 0, 0, z);
    py::object r = py::reinterpret_borrow<py::object>(reinterpret_cast<PyObject*>(&PyLong_Type))
                       .attr("from_bytes")(bytes, "little");
    if (mpz_sgn(z) < 0) r = propio(PyNumber_Negative(r.ptr()));
    return r;
}

}  // namespace alcp_py

// type_caster de pybind11 para alcp::BigInt (= mpz_class): cualquier parámetro
// `const BigInt&` o `std::vector<BigInt>` acepta enteros de Python, y cualquier
// BigInt devuelto se convierte en int.  OJO: las funciones deben devolver
// BigInt explícitamente, no una expresión de gmpxx (p. ej. `a + b`).
namespace pybind11::detail {
template <>
struct type_caster<mpz_class> {
    PYBIND11_TYPE_CASTER(mpz_class, const_name("int"));

    bool load(handle src, bool /*convert*/) {
        PyObject* o = src.ptr();
        if (o == nullptr) return false;
        if (PyLong_CheckExact(o)) {  // caso habitual: sin llamadas extra
            alcp_py::pylong_a_mpz(o, value.get_mpz_t());
            return true;
        }
        // Un float nunca es un entero válido (ni siquiera 2.0): error explícito.
        if (PyFloat_Check(o))
            throw type_error(std::string("se esperaba un entero (int o un objeto con __index__), "
                                         "no un valor de tipo '") + Py_TYPE(o)->tp_name + "'");
        if (!PyIndex_Check(o)) return false;  // str, list, ...: no es un entero
        // bool, subclases de int, numpy.int64, ... -> int exacto vía __index__
        // (si __index__ lanza, esa excepción se propaga tal cual).
        py::object idx = alcp_py::propio(PyNumber_Index(o));
        alcp_py::pylong_a_mpz(idx.ptr(), value.get_mpz_t());
        return true;
    }

    static handle cast(const mpz_class& v, return_value_policy /*policy*/, handle /*parent*/) {
        return alcp_py::mpz_a_pylong(v.get_mpz_t()).release();
    }
};
}  // namespace pybind11::detail

namespace alcp_py {

// Valor de un BigInt como entero de máquina T (k, n, índices, grados...).
// Negativo en un tipo sin signo -> std::invalid_argument (ValueError); fuera de
// rango -> std::overflow_error (OverflowError).
template <std::integral T>
T a_maquina(const BigInt& v, const char* nombre) {
    if constexpr (std::is_unsigned_v<T>) {
        if (sgn(v) < 0) throw std::invalid_argument(std::string(nombre) + " debe ser >= 0");
    }
    if (mpz_sizeinbase(v.get_mpz_t(), 2) <= 63) {
        unsigned long long mag = 0;
        mpz_export(&mag, nullptr, -1, sizeof mag, 0, 0, v.get_mpz_t());
        const long long x = sgn(v) < 0 ? -static_cast<long long>(mag) : static_cast<long long>(mag);
        if (std::in_range<T>(x)) return static_cast<T>(x);
    }
    throw std::overflow_error(std::string(nombre) + " es demasiado grande");
}

// ---------------------------------------------------------------------------
// Cuerpos.  El objeto Python posee el cuerpo C++ mediante std::shared_ptr.
// ---------------------------------------------------------------------------
struct FpW {
    using Field = Fp;
    Fp K;
    explicit FpW(const BigInt& p) : K(p) {}
};

struct FqW {
    using Field = Fq;
    Fq K;
    explicit FqW(const BigInt& q) : K(q) {}
    FqW(const BigInt& p, std::vector<BigInt> g) : K(Fp(p), std::move(g)) {}
};

// Textos que dependen del tipo de cuerpo (nombres de clase y docstrings).
template <class W> struct Info;
template <> struct Info<FpW> {
    static constexpr const char* cuerpo = "Fp";
    static constexpr const char* poli = "PolyFp";
    static constexpr const char* elem = "entero (se reduce módulo p)";
    static constexpr const char* sal = "entero en [0, p)";
};
template <> struct Info<FqW> {
    static constexpr const char* cuerpo = "Fq";
    static constexpr const char* poli = "PolyFq";
    static constexpr const char* elem = "índice en base p (se reduce módulo q)";
    static constexpr const char* sal = "índice en [0, q)";
};

// Docstring con marcadores: {E} = representación de un elemento de entrada,
// {S} = de salida, {K} = clase del cuerpo, {P} = clase de los polinomios.
// Devuelve un const char* con almacenamiento estático (vive lo que el módulo):
// así las llamadas a .def no crean std::string temporales.
template <class W>
const char* doc(const char* plantilla) {
    static std::forward_list<std::string> almacen;
    std::string s = plantilla;
    const std::pair<const char*, const char*> reglas[] = {
        {"{E}", Info<W>::elem}, {"{S}", Info<W>::sal}, {"{K}", Info<W>::cuerpo}, {"{P}", Info<W>::poli}};
    for (const auto& [marca, texto] : reglas)
        for (std::size_t i = s.find(marca); i != std::string::npos; i = s.find(marca, i))
            s.replace(i, 3, texto);
    return almacen.emplace_front(std::move(s)).c_str();
}

// Elemento <-> entero de Python, según el cuerpo.
inline Fp::Elem elem_in(const Fp& K, const BigInt& v) { return K.from_int(v); }
inline BigInt elem_out(const Fp& /*K*/, const Fp::Elem& e) { return e; }
inline Fq::Elem elem_in(const Fq& K, const BigInt& v) { return K.from_index(v); }
inline BigInt elem_out(const Fq& K, const Fq::Elem& e) { return K.to_index(e); }

// Listas completas.  Se construyen desde una vista (sin vector local: así no hay
// código de limpieza de excepciones en estas líneas).  No se usa
// std::ranges::to (C++23) porque la libstdc++ de GCC 13, la de la CI de Linux,
// no lo tiene.
template <class Field>
std::vector<typename Field::Elem> elems_in(const Field& K, const std::vector<BigInt>& v) {
    const auto vista = v | std::views::transform([&K](const BigInt& x) { return elem_in(K, x); });
    return std::vector<typename Field::Elem>(vista.begin(), vista.end());
}
template <class Field>
std::vector<BigInt> elems_out(const Field& K, const std::vector<typename Field::Elem>& v) {
    const auto vista = v | std::views::transform([&K](const typename Field::Elem& x) { return elem_out(K, x); });
    return std::vector<BigInt>(vista.begin(), vista.end());
}

// Identidad de un cuerpo (para __hash__): sus parámetros, como en operator==.
inline py::tuple clave(const Fp& K) { return py::make_tuple("Fp", K.characteristic()); }
inline py::tuple clave(const Fq& K) {
    return py::make_tuple("Fq", K.characteristic(), py::tuple(py::cast(K.modulus())));
}

inline std::string repr_cuerpo(const Fp& K) { return "Fp(" + K.characteristic().get_str() + ")"; }
inline std::string repr_cuerpo(const Fq& K) {
    std::string s = "Fq(" + K.characteristic().get_str() + ", [";
    for (std::size_t i = 0; i < K.modulus().size(); ++i) s += (i ? ", " : "") + K.modulus()[i].get_str();
    return s + "])";
}

// ---------------------------------------------------------------------------
// Polinomios.  PolyW comparte el cuerpo (shared_ptr) para que el Field* que
// guarda Poly no quede colgando.
// ---------------------------------------------------------------------------
template <class W>
struct PolyW {
    using Field = typename W::Field;
    using P = alcp::Poly<Field>;

    std::shared_ptr<W> F;  // mantiene vivo el cuerpo al que apunta p
    P p;

    PolyW(std::shared_ptr<W> f, P q) : F(std::move(f)), p(std::move(q)) {
        // Invariante: p apunta EXACTAMENTE al cuerpo que mantiene vivo F.  Un
        // resultado del núcleo puede apuntar al cuerpo del OTRO operando (p. ej.
        // gcd(a, b) devuelve a veces el propio b), que es igual en parámetros
        // (el núcleo lo ha comprobado) pero es otro objeto, que podría
        // destruirse antes.  En ese caso se reconstruye sobre F.
        if (&p.field() != &F->K) p = P(F->K, p.coeffs());
    }

    const Field& K() const { return F->K; }

    // Construye desde una lista de enteros.  El Poly se construye ANTES de mover
    // el shared_ptr (el código original hacía PolyW(std::move(k), Poly(*k, ...)),
    // con orden de evaluación no especificado: podía desreferenciar un puntero
    // ya movido).
    static PolyW desde_lista(std::shared_ptr<W> f, const std::vector<BigInt>& cs) {
        P poly(f->K, elems_in(f->K, cs));
        return PolyW(std::move(f), std::move(poly));
    }

    // Envuelve un resultado del núcleo sobre el mismo cuerpo que *this.
    PolyW mismo(P q) const { return PolyW(F, std::move(q)); }

    std::vector<BigInt> coeffs() const { return elems_out(K(), p.coeffs()); }
};

inline alcp::MulAlgo algoritmo(const std::string& s) {
    if (s == "auto") return alcp::MulAlgo::Auto;
    if (s == "school" || s == "schoolbook") return alcp::MulAlgo::Schoolbook;
    if (s == "karatsuba") return alcp::MulAlgo::Karatsuba;
    if (s == "fft") return alcp::MulAlgo::FFT;
    throw std::invalid_argument("algoritmo de multiplicación desconocido: '" + s +
                                "' (válidos: 'auto', 'school', 'karatsuba', 'fft')");
}

// ---------------------------------------------------------------------------
// Docstrings (en español).  Marcadores que sustituye doc<W>(): {E} = elemento
// de entrada, {S} = elemento de salida, {K} = clase del cuerpo, {P} = clase de
// los polinomios.
// ---------------------------------------------------------------------------
namespace docs {
constexpr const char* modulo =
    "Cuerpos finitos, factorización de polinomios y algoritmos rápidos (núcleo C++/GMP).\n\n"
    "Elementos: en Fp, enteros (se reducen módulo p); en Fq, índices en base p "
    "(c0 + c1·p + ... + c_{n-1}·p^(n-1) representa c0 + c1·a + ... + c_{n-1}·a^(n-1)). "
    "Los enteros pueden ser de cualquier tamaño y de cualquier tipo con __index__ "
    "(bool, numpy.int64...); los float se rechazan con TypeError.\n\n"
    "Contenido: is_prime, prime_power_factor, set_seed; los cuerpos Fp y Fq; los "
    "polinomios PolyFp y PolyFq; dft, idft, toeplitz_lower_vec, toeplitz_upper_vec, "
    "toeplitz_vec, toeplitz_lower_inverse, toeplitz_upper_inverse y divmod_toeplitz "
    "(para Fp y Fq).\n\n"
    "Errores: ZeroDivisionError (invertir 0, dividir por el polinomio cero), ValueError "
    "(argumentos inválidos, operandos de cuerpos distintos...), OverflowError (enteros "
    "que no caben donde el núcleo necesita un entero de máquina), TypeError (tipos).\n\n"
    "El GIL no se libera: el núcleo tiene estado global (generador aleatorio, cachés FFT).";

// ---- funciones del módulo ----
constexpr const char* is_prime = "¿Es n primo? Test probabilístico de GMP (Baillie-PSW + Miller-Rabin, 40 rondas).";
constexpr const char* prime_power_factor = "(p, n) con q = p^n y p primo. ValueError si q no es potencia de un primo.";
constexpr const char* int_roundtrip =
    "Uso interno (pruebas): convierte x a mpz y de vuelta a int (comprueba la conversión de enteros en "
    "ambos sentidos). No forma parte de la API pública.";
constexpr const char* set_seed =
    "Fija la semilla (entero >= 0) del generador aleatorio global del núcleo (random, random_poly, "
    "módulos aleatorios de Fq, factorización...), para resultados reproducibles.";

// ---- cuerpos ----
constexpr const char* fp = "Cuerpo primo F_p. Sus elementos son enteros en [0, p).";
constexpr const char* fp_init = "Fp(p): F_p para un primo p. ValueError si p no es primo.";
constexpr const char* fp_n = "Grado de extensión (1).";
constexpr const char* fq = "Cuerpo F_q = F_p[a]/(g), q = p^n. Sus elementos son índices en base p, en [0, q).";
constexpr const char* fq_init_q =
    "Fq(q): F_q con q = p^n y un módulo mónico irreducible aleatorio de grado n (reproducible con "
    "set_seed). ValueError si q no es potencia de un primo.";
constexpr const char* fq_init_mod =
    "Fq(p, modulus): F_p[a]/(g) con g = modulus en orden creciente (se hace mónico). ValueError si p "
    "no es primo o g no es irreducible de grado >= 1.";
constexpr const char* fq_n = "Grado de extensión n.";
constexpr const char* fq_modulus = "Módulo mónico g en orden creciente (longitud n + 1).";
constexpr const char* fq_from_index = "Índice normalizado: i mod q.";
constexpr const char* fq_from_coeffs =
    "Índice del elemento c0 + c1·a + c2·a^2 + ... (orden creciente). Conversión generosa: cada "
    "coeficiente se reduce módulo p, las listas de más de n coeficientes se reducen módulo g y [] es el 0.";
constexpr const char* fq_to_coeffs =
    "Coeficientes [c0, ..., c_{n-1}] (longitud n, en [0, p)) del elemento de índice i.";

constexpr const char* c_p = "Característica p del cuerpo.";
constexpr const char* c_q = "Cardinal q = p^n del cuerpo.";
constexpr const char* c_eq = "Igualdad por parámetros (mismos p y módulo). Con otro tipo devuelve NotImplemented.";
constexpr const char* c_hash = "Hash coherente con __eq__ (calculado a partir de los parámetros).";
constexpr const char* c_repr = "Representación evaluable: Fp(p) o Fq(p, [g0, ..., gn]).";
constexpr const char* c_copia = "Los cuerpos son inmutables: devuelve el mismo objeto.";
constexpr const char* c_fft_capacity =
    "Mayor k tal que 2^k divide a q - 1: tamaño máximo 2^k de una FFT en este cuerpo.";
constexpr const char* c_fft_base_root =
    "Raíz primitiva 2^fft_capacity()-ésima de la unidad que usa la FFT ({S}). Se elige al azar la "
    "primera vez y queda fijada en el cuerpo.";
constexpr const char* c_from_int = "Inmersión del entero a en el cuerpo: (a mod p)·1, como {S}.";
constexpr const char* c_random = "Elemento aleatorio uniforme ({S}); usa el generador global (set_seed).";
constexpr const char* c_add = "a + b. Entradas: {E}; salida: {S}.";
constexpr const char* c_sub = "a - b. Entradas: {E}; salida: {S}.";
constexpr const char* c_neg = "-a. Entrada: {E}; salida: {S}.";
constexpr const char* c_mul = "a · b. Entradas: {E}; salida: {S}.";
constexpr const char* c_div = "a / b. Entradas: {E}; salida: {S}. ZeroDivisionError si b = 0.";
constexpr const char* c_inv = "Inverso a^(-1). Entrada: {E}; salida: {S}. ZeroDivisionError si a = 0.";
constexpr const char* c_pow =
    "a^k para cualquier entero k (k < 0: potencia del inverso). Entrada: {E}; salida: {S}. "
    "ZeroDivisionError si a = 0 y k < 0.";
constexpr const char* c_poly = "Polinomio {P} con coeficientes en orden creciente [c0, c1, ...] (cada uno: {E}).";
constexpr const char* c_random_poly =
    "Polinomio {P} aleatorio de grado exactamente d (coeficiente principal no nulo); si d < 0, el "
    "polinomio cero.";

// ---- polinomios ----
constexpr const char* p_clase =
    "Polinomio inmutable con coeficientes en {K} (cada coeficiente: {E}).\n\n"
    "Se construye con {K}.poly(coeffs) o {P}(field, coeffs). Los operadores binarios exigen operandos "
    "sobre el mismo cuerpo (mismos parámetros); si no, ValueError.";
constexpr const char* p_init = "{P}(field, coeffs=[]): coeficientes en orden creciente (cada uno: {E}).";
constexpr const char* p_field =
    "Cuerpo {K} del polinomio. Devuelve el MISMO objeto Python con el que se creó mientras siga vivo; "
    "si ya no existe, un objeto nuevo que comparte el mismo cuerpo C++ (igual por ==).";
constexpr const char* p_coeffs = "Coeficientes en orden creciente, sin ceros finales ({S}); [] para el cero.";
constexpr const char* p_coeff = "Coeficiente de x^i ({S}); 0 si i > grado. ValueError si i < 0.";
constexpr const char* p_degree = "Grado (-1 para el polinomio cero).";
constexpr const char* p_is_zero = "¿Es el polinomio cero?";
constexpr const char* p_is_one = "¿Es el polinomio constante 1?";
constexpr const char* p_leading = "Coeficiente principal ({S}). ValueError si el polinomio es cero.";
constexpr const char* p_eq =
    "Igualdad: mismo cuerpo (por parámetros) y mismos coeficientes. Con otro tipo devuelve "
    "NotImplemented; sobre cuerpos distintos, False.";
constexpr const char* p_hash = "Hash coherente con __eq__ (parámetros del cuerpo y coeficientes).";
constexpr const char* p_copy = "Copia (comparte el cuerpo).";
constexpr const char* p_deepcopy = "Copia; el cuerpo, inmutable, se comparte.";
constexpr const char* p_add = "Suma.";
constexpr const char* p_sub = "Resta.";
constexpr const char* p_neg = "Opuesto.";
constexpr const char* p_mul_op = "Producto (elige el algoritmo según el tamaño, como mul(b, 'auto')).";
constexpr const char* p_floordiv = "Cociente euclídeo. ZeroDivisionError si b = 0.";
constexpr const char* p_mod = "Resto euclídeo. ZeroDivisionError si b = 0.";
constexpr const char* p_divmod_op = "divmod(a, b) -> (cociente, resto).";
constexpr const char* p_divmod =
    "(cociente, resto) con a = q·b + r y grado(r) < grado(b). ZeroDivisionError si b = 0.";
constexpr const char* p_pow = "Potencia a^k (k >= 0; ValueError si k < 0).";
constexpr const char* p_mul =
    "Producto con el algoritmo indicado: 'auto' (por defecto), 'school' (o 'schoolbook'), 'karatsuba' "
    "o 'fft' (si el cuerpo no admite ese tamaño de FFT, Karatsuba). Otro nombre: ValueError.";
constexpr const char* p_scalar_mul = "Producto por el escalar c ({E}).";
constexpr const char* p_derivative = "Derivada formal.";
constexpr const char* p_monic = "Polinomio mónico asociado (el cero se devuelve tal cual).";
constexpr const char* p_evaluate = "Valor f(x) por Horner. x: {E}; resultado: {S}.";
constexpr const char* p_compose = "Composición f(g(x)).";
constexpr const char* p_is_irreducible = "Test de irreducibilidad de Ben-Or (False para constantes y el cero).";
constexpr const char* p_is_squarefree = "¿Libre de cuadrados? (el cero no lo es; una constante no nula sí).";
constexpr const char* p_gcd = "Máximo común divisor mónico (gcd(0, 0) = 0).";
constexpr const char* p_gcd_ext = "Algoritmo de Euclides extendido: (g, x, y) con g = x·a + y·b, g mónico.";
constexpr const char* p_inv_mod = "Inverso de a módulo mod. ValueError si no es invertible.";
constexpr const char* p_pow_mod =
    "a^k módulo mod (k < 0: potencia del inverso modular). ZeroDivisionError si mod = 0.";
constexpr const char* p_multiplicity = "Mayor e tal que u^e divide a a. ValueError si grado(u) < 1 o a = 0.";
constexpr const char* p_radical =
    "Parte libre de cuadrados (producto de los factores irreducibles distintos, mónico); rad(0) = 0 y "
    "rad(c) = 1 para una constante c no nula.";
constexpr const char* p_distinct_degree =
    "Factorización en grados distintos de un polinomio libre de cuadrados: lista de (d, producto de "
    "los factores irreducibles de grado d). ValueError si grado < 1.";
constexpr const char* p_equal_degree =
    "Factorización en grado igual: lista de los factores irreducibles (mónicos) de grado r de un "
    "producto de irreducibles distintos de grado r. ValueError si no se cumple la precondición.";
constexpr const char* p_factor =
    "Factorización de Cantor-Zassenhaus: (lista de ({P} irreducible mónico, multiplicidad) en orden "
    "canónico, coeficiente principal ({S})). ValueError para el polinomio cero.";
constexpr const char* p_to_string = "Representación textual del núcleo, p. ej. '(1) + (3)*x^2'.";
constexpr const char* p_repr = "Representación compatible con la API original: {P}(<to_string()>).";

// ---- algoritmos rápidos ----
constexpr const char* r_dft =
    "DFT_{n,root}(a) en {K} con n = 2^k (Cooley-Tukey). root: raíz de la unidad de orden EXACTAMENTE "
    "2^k; a: lista de 2^k elementos ({E}). ValueError si no se cumplen las precondiciones.";
constexpr const char* r_idft =
    "Transformada inversa: n^(-1)·DFT_{n,root^(-1)}(a) en {K}; mismas precondiciones que dft.";
constexpr const char* r_lower_vec =
    "T·b con T triangular inferior de Toeplitz n×n de primera columna a[0..n-1] (len(a) >= n, "
    "len(b) == n) en {K}.";
constexpr const char* r_upper_vec =
    "T·b con T triangular superior de Toeplitz n×n de primera fila a[0..n-1] (len(a) >= n, "
    "len(b) == n) en {K}.";
constexpr const char* r_vec =
    "T·b con T de Toeplitz n×n completa en {K}: a (longitud 2n-1) es la primera fila seguida de la "
    "primera columna sin la esquina; len(b) == n.";
constexpr const char* r_lower_inverse =
    "Primera columna de T^(-1) con T triangular inferior de Toeplitz de primera columna a "
    "(len(a) >= n) en {K}. ZeroDivisionError si a[0] = 0.";
constexpr const char* r_upper_inverse =
    "Primera fila de T^(-1) con T triangular superior de Toeplitz de primera fila a (len(a) >= n) en "
    "{K}. ZeroDivisionError si a[0] = 0.";
constexpr const char* r_divmod =
    "División euclídea f = q·g + r de polinomios {P} reformulada como sistema triangular de Toeplitz: "
    "devuelve (q, r). ZeroDivisionError si g = 0; ValueError si f y g están sobre cuerpos distintos.";
}  // namespace docs

// ---------------------------------------------------------------------------
// Operaciones expuestas.  Son funciones con nombre (no lambdas dentro del
// registro) para que el registro sea una lista de una línea por función.
// Nota de cobertura: con gcov, una sentencia .def repartida en varias líneas
// deja en la primera solo el bloque de limpieza de excepciones (nunca
// ejecutado) y gcovr, que cuenta cada función por separado, la marcaría como
// no cubierta; por eso cada .def ocupa una sola línea.
// ---------------------------------------------------------------------------

// Parte común de Fp y Fq.
template <class W>
struct OpsCuerpo {
    using PW = PolyW<W>;
    using P = typename PW::P;

    static BigInt p(const W& f) { return f.K.characteristic(); }
    static BigInt q(const W& f) { return f.K.cardinality(); }
    static bool eq(const W& a, const W& b) { return a.K == b.K; }
    static py::ssize_t hash(const W& f) { return py::hash(clave(f.K)); }
    static std::string repr(const W& f) { return repr_cuerpo(f.K); }
    static std::shared_ptr<W> copia(const std::shared_ptr<W>& f) { return f; }
    static std::shared_ptr<W> copia_profunda(const std::shared_ptr<W>& f, const py::dict& /*memo*/) { return f; }
    static unsigned long fft_capacity(const W& f) { return f.K.fft_capacity(); }
    static BigInt fft_base_root(const W& f) { return elem_out(f.K, f.K.fft_base_root()); }
    static BigInt from_int(const W& f, const BigInt& a) { return elem_out(f.K, f.K.from_int(a)); }
    static BigInt random(const W& f) { return elem_out(f.K, f.K.random()); }
    static BigInt add(const W& f, const BigInt& a, const BigInt& b) {
        return elem_out(f.K, f.K.add(elem_in(f.K, a), elem_in(f.K, b)));
    }
    static BigInt sub(const W& f, const BigInt& a, const BigInt& b) {
        return elem_out(f.K, f.K.sub(elem_in(f.K, a), elem_in(f.K, b)));
    }
    static BigInt neg(const W& f, const BigInt& a) { return elem_out(f.K, f.K.neg(elem_in(f.K, a))); }
    static BigInt mul(const W& f, const BigInt& a, const BigInt& b) {
        return elem_out(f.K, f.K.mul(elem_in(f.K, a), elem_in(f.K, b)));
    }
    static BigInt div(const W& f, const BigInt& a, const BigInt& b) {
        return elem_out(f.K, f.K.div(elem_in(f.K, a), elem_in(f.K, b)));
    }
    static BigInt inv(const W& f, const BigInt& a) { return elem_out(f.K, f.K.inv(elem_in(f.K, a))); }
    static BigInt pow(const W& f, const BigInt& a, const BigInt& k) { return elem_out(f.K, f.K.pow(elem_in(f.K, a), k)); }
    static PW poly(std::shared_ptr<W> f, const std::vector<BigInt>& cs) { return PW::desde_lista(std::move(f), cs); }
    static PW random_poly(std::shared_ptr<W> f, const BigInt& d) {
        // d < 0 da el polinomio cero.  El núcleo reserva d + 1 coeficientes en
        // un int: con d = INT_MAX se desbordaría, así que se exige d + 1 <= INT_MAX.
        const int g = sgn(d) < 0 ? -1 : a_maquina<int>(d + 1, "d") - 1;
        P q = P::random(f->K, g);
        return PW(std::move(f), std::move(q));
    }
};

// Exclusivo de Fp.
struct OpsFp {
    static std::shared_ptr<FpW> crear(const BigInt& p) { return std::make_shared<FpW>(p); }
    static int n(const FpW& /*f*/) { return 1; }
};

// Exclusivo de Fq.
struct OpsFq {
    static std::shared_ptr<FqW> crear_q(const BigInt& q) { return std::make_shared<FqW>(q); }
    static std::shared_ptr<FqW> crear_mod(const BigInt& p, const std::vector<BigInt>& g) {
        return std::make_shared<FqW>(p, g);
    }
    static unsigned long n(const FqW& f) { return f.K.extension_degree(); }
    static std::vector<BigInt> modulus(const FqW& f) { return f.K.modulus(); }
    static BigInt from_index(const FqW& f, const BigInt& i) { return f.K.to_index(f.K.from_index(i)); }
    static BigInt from_coeffs(const FqW& f, const std::vector<BigInt>& v) { return f.K.to_index(f.K.from_vector(v)); }
    static std::vector<BigInt> to_coeffs(const FqW& f, const BigInt& i) { return f.K.from_index(i); }
};

// Polinomios (idéntico para PolyFp y PolyFq).
template <class W>
struct OpsPoli {
    using PW = PolyW<W>;
    using P = typename PW::P;

    static PW crear(std::shared_ptr<W> f, const std::vector<BigInt>& cs) { return PW::desde_lista(std::move(f), cs); }
    static std::shared_ptr<W> field(const PW& a) { return a.F; }
    static std::vector<BigInt> coeffs(const PW& a) { return a.coeffs(); }
    static BigInt coeff(const PW& a, const BigInt& i) { return elem_out(a.K(), a.p.coeff(a_maquina<std::size_t>(i, "i"))); }
    static int degree(const PW& a) { return a.p.degree(); }
    static bool is_zero(const PW& a) { return a.p.is_zero(); }
    static bool is_one(const PW& a) { return a.p.is_one(); }
    static BigInt leading(const PW& a) { return elem_out(a.K(), a.p.lead()); }
    static bool eq(const PW& a, const PW& b) { return a.p == b.p; }
    static py::ssize_t hash(const PW& a) {
        return py::hash(py::make_tuple(clave(a.K()), py::tuple(py::cast(a.coeffs()))));
    }
    static PW copia(const PW& a) { return a; }
    static PW copia_profunda(const PW& a, const py::dict& /*memo*/) { return a; }
    static PW add(const PW& a, const PW& b) { return a.mismo(a.p + b.p); }
    static PW sub(const PW& a, const PW& b) { return a.mismo(a.p - b.p); }
    static PW neg(const PW& a) { return a.mismo(-a.p); }
    static PW mul_op(const PW& a, const PW& b) { return a.mismo(a.p * b.p); }
    static PW floordiv(const PW& a, const PW& b) { return a.mismo(a.p / b.p); }
    static PW mod(const PW& a, const PW& b) { return a.mismo(a.p % b.p); }
    static py::tuple divmod(const PW& a, const PW& b) {
        auto [q, r] = a.p.divmod(b.p);
        return py::make_tuple(a.mismo(std::move(q)), a.mismo(std::move(r)));
    }
    static PW pow(const PW& a, const BigInt& k) { return a.mismo(a.p.pow(k)); }
    static PW mul(const PW& a, const PW& b, const std::string& algo) { return a.mismo(a.p.mul(b.p, algoritmo(algo))); }
    static PW scalar_mul(const PW& a, const BigInt& e) { return a.mismo(a.p.scalar_mul(elem_in(a.K(), e))); }
    static PW derivative(const PW& a) { return a.mismo(a.p.derivative()); }
    static PW monic(const PW& a) { return a.mismo(a.p.monic()); }
    static BigInt evaluate(const PW& a, const BigInt& x) { return elem_out(a.K(), a.p.evaluate(elem_in(a.K(), x))); }
    static PW compose(const PW& a, const PW& g) { return a.mismo(a.p.compose(g.p)); }
    static bool is_irreducible(const PW& a) { return a.p.is_irreducible(); }
    static bool is_squarefree(const PW& a) { return a.p.is_squarefree(); }
    static PW gcd(const PW& a, const PW& b) { return a.mismo(P::gcd(a.p, b.p)); }
    static py::tuple gcd_ext(const PW& a, const PW& b) {
        auto e = P::gcd_ext(a.p, b.p);
        return py::make_tuple(a.mismo(std::move(e.g)), a.mismo(std::move(e.x)), a.mismo(std::move(e.y)));
    }
    static PW inv_mod(const PW& a, const PW& m) { return a.mismo(P::inv_mod(a.p, m.p)); }
    // El módulo cero se rechaza antes de llamar al núcleo: con k < 0, P::pow_mod
    // invierte primero módulo m y lanzaría «no invertible» (ValueError) en lugar
    // del ZeroDivisionError documentado.
    static PW pow_mod(const PW& a, const BigInt& k, const PW& m) {
        if (m.p.is_zero()) throw alcp::division_by_zero("pow_mod: el módulo es el polinomio cero");
        return a.mismo(P::pow_mod(a.p, k, m.p));
    }
    static int multiplicity(const PW& a, const PW& u) { return a.p.multiplicity(u.p); }
    static PW radical(const PW& a) { return a.mismo(a.p.radical()); }
    static py::list distinct_degree(const PW& a) {
        py::list out;
        for (auto& [d, prod] : a.p.distinct_degree()) out.append(py::make_tuple(d, a.mismo(prod)));
        return out;
    }
    static py::list equal_degree(const PW& a, const BigInt& r) {
        py::list out;
        for (auto& u : a.p.equal_degree(a_maquina<int>(r, "r"))) out.append(a.mismo(u));
        return out;
    }
    static py::tuple factor(const PW& a) {
        auto f = a.p.factor();
        py::list factores;
        for (auto& [u, e] : f.factors) factores.append(py::make_tuple(a.mismo(u), e));
        return py::make_tuple(factores, elem_out(a.K(), f.leading));
    }
    static std::string to_string(const PW& a, const std::string& var) { return a.p.to_string(var); }
    static std::string repr(const PW& a) { return std::string(Info<W>::poli) + "(" + a.p.to_string() + ")"; }
};

// Algoritmos rápidos (Práctica 3): mismas funciones para Fp y Fq.
template <class W>
struct OpsRapidos {
    using PW = PolyW<W>;
    using V = std::vector<BigInt>;

    static V dft(const W& f, const BigInt& raiz, const BigInt& k, const V& a) {
        return elems_out(f.K, alcp::fast::dft(f.K, elem_in(f.K, raiz), a_maquina<unsigned long>(k, "k"), elems_in(f.K, a)));
    }
    static V idft(const W& f, const BigInt& raiz, const BigInt& k, const V& a) {
        return elems_out(f.K, alcp::fast::idft(f.K, elem_in(f.K, raiz), a_maquina<unsigned long>(k, "k"), elems_in(f.K, a)));
    }
    static V lower_vec(const W& f, const BigInt& n, const V& a, const V& b) {
        return elems_out(f.K, alcp::fast::toeplitz_lower_vec(f.K, a_maquina<std::size_t>(n, "n"), elems_in(f.K, a),
                                                             elems_in(f.K, b)));
    }
    static V upper_vec(const W& f, const BigInt& n, const V& a, const V& b) {
        return elems_out(f.K, alcp::fast::toeplitz_upper_vec(f.K, a_maquina<std::size_t>(n, "n"), elems_in(f.K, a),
                                                             elems_in(f.K, b)));
    }
    static V vec(const W& f, const BigInt& n, const V& a, const V& b) {
        return elems_out(f.K, alcp::fast::toeplitz_vec(f.K, a_maquina<std::size_t>(n, "n"), elems_in(f.K, a),
                                                       elems_in(f.K, b)));
    }
    static V lower_inverse(const W& f, const BigInt& n, const V& a) {
        return elems_out(f.K, alcp::fast::toeplitz_lower_inverse(f.K, a_maquina<std::size_t>(n, "n"), elems_in(f.K, a)));
    }
    static V upper_inverse(const W& f, const BigInt& n, const V& a) {
        return elems_out(f.K, alcp::fast::toeplitz_upper_inverse(f.K, a_maquina<std::size_t>(n, "n"), elems_in(f.K, a)));
    }
    static py::tuple divmod(const PW& f, const PW& g) {
        auto [q, r] = alcp::fast::toeplitz_divmod(f.K(), f.p, g.p);
        return py::make_tuple(f.mismo(std::move(q)), f.mismo(std::move(r)));
    }
};

// ---------------------------------------------------------------------------
// Registro (una línea por función; ver la nota de cobertura de arriba)
// ---------------------------------------------------------------------------
template <class W>
void registrar_cuerpo(py::class_<W, std::shared_ptr<W>>& c) {
    using O = OpsCuerpo<W>;
    c.def_property_readonly("p", &O::p, docs::c_p);
    c.def_property_readonly("q", &O::q, docs::c_q);
    c.def("__eq__", &O::eq, py::is_operator(), py::arg("other"), docs::c_eq);
    c.def("__hash__", &O::hash, docs::c_hash);
    c.def("__repr__", &O::repr, docs::c_repr);
    c.def("__copy__", &O::copia, docs::c_copia);
    c.def("__deepcopy__", &O::copia_profunda, py::arg("memo"), docs::c_copia);
    c.def("fft_capacity", &O::fft_capacity, docs::c_fft_capacity);
    c.def("fft_base_root", &O::fft_base_root, doc<W>(docs::c_fft_base_root));
    c.def("from_int", &O::from_int, py::arg("a"), doc<W>(docs::c_from_int));
    c.def("random", &O::random, doc<W>(docs::c_random));
    c.def("add", &O::add, py::arg("a"), py::arg("b"), doc<W>(docs::c_add));
    c.def("sub", &O::sub, py::arg("a"), py::arg("b"), doc<W>(docs::c_sub));
    c.def("neg", &O::neg, py::arg("a"), doc<W>(docs::c_neg));
    c.def("mul", &O::mul, py::arg("a"), py::arg("b"), doc<W>(docs::c_mul));
    c.def("div", &O::div, py::arg("a"), py::arg("b"), doc<W>(docs::c_div));
    c.def("inv", &O::inv, py::arg("a"), doc<W>(docs::c_inv));
    c.def("pow", &O::pow, py::arg("a"), py::arg("k"), doc<W>(docs::c_pow));
    c.def("poly", &O::poly, py::arg("coeffs"), doc<W>(docs::c_poly));
    c.def("random_poly", &O::random_poly, py::arg("d"), doc<W>(docs::c_random_poly));
}

template <class W>
void registrar_poli(py::module_& m) {
    using O = OpsPoli<W>;
    py::class_<PolyW<W>> c(m, Info<W>::poli, doc<W>(docs::p_clase));
    c.def(py::init(&O::crear), py::arg("field"), py::arg("coeffs") = std::vector<BigInt>{}, doc<W>(docs::p_init));
    c.def_property_readonly("field", &O::field, doc<W>(docs::p_field));
    c.def("coeffs", &O::coeffs, doc<W>(docs::p_coeffs));
    c.def("coeff", &O::coeff, py::arg("i"), doc<W>(docs::p_coeff));
    c.def("degree", &O::degree, docs::p_degree);
    c.def("is_zero", &O::is_zero, docs::p_is_zero);
    c.def("is_one", &O::is_one, docs::p_is_one);
    c.def("leading", &O::leading, doc<W>(docs::p_leading));
    c.def("__eq__", &O::eq, py::is_operator(), py::arg("other"), docs::p_eq);
    c.def("__hash__", &O::hash, docs::p_hash);
    c.def("__copy__", &O::copia, docs::p_copy);
    c.def("__deepcopy__", &O::copia_profunda, py::arg("memo"), docs::p_deepcopy);
    c.def("__add__", &O::add, py::is_operator(), py::arg("other"), docs::p_add);
    c.def("__sub__", &O::sub, py::is_operator(), py::arg("other"), docs::p_sub);
    c.def("__neg__", &O::neg, docs::p_neg);
    c.def("__mul__", &O::mul_op, py::is_operator(), py::arg("other"), docs::p_mul_op);
    c.def("__floordiv__", &O::floordiv, py::is_operator(), py::arg("other"), docs::p_floordiv);
    c.def("__mod__", &O::mod, py::is_operator(), py::arg("other"), docs::p_mod);
    c.def("__divmod__", &O::divmod, py::is_operator(), py::arg("other"), docs::p_divmod_op);
    c.def("divmod", &O::divmod, py::arg("other"), docs::p_divmod);
    c.def("__pow__", &O::pow, py::is_operator(), py::arg("k"), docs::p_pow);
    c.def("mul", &O::mul, py::arg("other"), py::arg("algo") = "auto", docs::p_mul);
    c.def("scalar_mul", &O::scalar_mul, py::arg("c"), doc<W>(docs::p_scalar_mul));
    c.def("derivative", &O::derivative, docs::p_derivative);
    c.def("monic", &O::monic, docs::p_monic);
    c.def("evaluate", &O::evaluate, py::arg("x"), doc<W>(docs::p_evaluate));
    c.def("compose", &O::compose, py::arg("g"), docs::p_compose);
    c.def("is_irreducible", &O::is_irreducible, docs::p_is_irreducible);
    c.def("is_squarefree", &O::is_squarefree, docs::p_is_squarefree);
    c.def("gcd", &O::gcd, py::arg("other"), docs::p_gcd);
    c.def("gcd_ext", &O::gcd_ext, py::arg("other"), docs::p_gcd_ext);
    c.def("inv_mod", &O::inv_mod, py::arg("mod"), docs::p_inv_mod);
    c.def("pow_mod", &O::pow_mod, py::arg("k"), py::arg("mod"), docs::p_pow_mod);
    c.def("multiplicity", &O::multiplicity, py::arg("u"), docs::p_multiplicity);
    c.def("radical", &O::radical, docs::p_radical);
    c.def("distinct_degree", &O::distinct_degree, docs::p_distinct_degree);
    c.def("equal_degree", &O::equal_degree, py::arg("r"), docs::p_equal_degree);
    c.def("factor", &O::factor, doc<W>(docs::p_factor));
    c.def("to_string", &O::to_string, py::arg("var") = "x", docs::p_to_string);
    c.def("__repr__", &O::repr, doc<W>(docs::p_repr));
}

template <class W>
void registrar_rapidos(py::module_& m) {
    using O = OpsRapidos<W>;
    m.def("dft", &O::dft, py::arg("field"), py::arg("root"), py::arg("k"), py::arg("a"), doc<W>(docs::r_dft));
    m.def("idft", &O::idft, py::arg("field"), py::arg("root"), py::arg("k"), py::arg("a"), doc<W>(docs::r_idft));
    m.def("toeplitz_lower_vec", &O::lower_vec, py::arg("field"), py::arg("n"), py::arg("a"), py::arg("b"), doc<W>(docs::r_lower_vec));
    m.def("toeplitz_upper_vec", &O::upper_vec, py::arg("field"), py::arg("n"), py::arg("a"), py::arg("b"), doc<W>(docs::r_upper_vec));
    m.def("toeplitz_vec", &O::vec, py::arg("field"), py::arg("n"), py::arg("a"), py::arg("b"), doc<W>(docs::r_vec));
    m.def("toeplitz_lower_inverse", &O::lower_inverse, py::arg("field"), py::arg("n"), py::arg("a"), doc<W>(docs::r_lower_inverse));
    m.def("toeplitz_upper_inverse", &O::upper_inverse, py::arg("field"), py::arg("n"), py::arg("a"), doc<W>(docs::r_upper_inverse));
    m.def("divmod_toeplitz", &O::divmod, py::arg("f"), py::arg("g"), doc<W>(docs::r_divmod));
}

// Traductor de excepciones del núcleo.  division_by_zero deriva de
// std::domain_error, así que se captura ANTES que su base.  Lo que no se
// captura aquí sigue a los traductores por defecto de pybind11 (p. ej.
// std::length_error -> ValueError, std::bad_alloc -> MemoryError).
inline void traducir_excepcion(std::exception_ptr e) {
    try {
        if (e) std::rethrow_exception(e);
    } catch (const alcp::division_by_zero& ex) {
        PyErr_SetString(PyExc_ZeroDivisionError, ex.what());
    } catch (const std::invalid_argument& ex) {
        PyErr_SetString(PyExc_ValueError, ex.what());
    } catch (const std::domain_error& ex) {
        PyErr_SetString(PyExc_ValueError, ex.what());
    } catch (const std::overflow_error& ex) {
        PyErr_SetString(PyExc_OverflowError, ex.what());
    }
}

// Funciones del módulo.
inline bool es_primo(const BigInt& n) { return alcp::is_prime(n); }
inline std::pair<BigInt, BigInt> potencia_de_primo(const BigInt& q) { return alcp::prime_power_factor(q); }
inline BigInt ida_y_vuelta(const BigInt& x) { return x; }
inline void fijar_semilla(const BigInt& s) {
    if (sgn(s) < 0) throw std::invalid_argument("set_seed: la semilla debe ser >= 0");
    if (alcp::fits_ulong(s)) alcp::set_seed(alcp::to_ulong(s));
    else alcp::detail::gmp_rng().seed(s);  // semillas que no caben en unsigned long
}

}  // namespace alcp_py

// ---------------------------------------------------------------------------
// Módulo
// ---------------------------------------------------------------------------
PYBIND11_MODULE(alcp_native, m) {
    using namespace alcp_py;
    m.doc() = docs::modulo;
    py::register_exception_translator(&traducir_excepcion);

    m.def("is_prime", &es_primo, py::arg("n"), docs::is_prime);
    m.def("prime_power_factor", &potencia_de_primo, py::arg("q"), docs::prime_power_factor);
    m.def("_int_roundtrip", &ida_y_vuelta, py::arg("x"), docs::int_roundtrip);
    m.def("set_seed", &fijar_semilla, py::arg("seed"), docs::set_seed);

    py::class_<FpW, std::shared_ptr<FpW>> fp(m, "Fp", docs::fp);
    fp.def(py::init(&OpsFp::crear), py::arg("p"), docs::fp_init);
    fp.def_property_readonly("n", &OpsFp::n, docs::fp_n);
    registrar_cuerpo(fp);

    py::class_<FqW, std::shared_ptr<FqW>> fq(m, "Fq", docs::fq);
    fq.def(py::init(&OpsFq::crear_q), py::arg("q"), docs::fq_init_q);
    fq.def(py::init(&OpsFq::crear_mod), py::arg("p"), py::arg("modulus"), docs::fq_init_mod);
    fq.def_property_readonly("n", &OpsFq::n, docs::fq_n);
    fq.def_property_readonly("modulus", &OpsFq::modulus, docs::fq_modulus);
    fq.def("from_index", &OpsFq::from_index, py::arg("i"), docs::fq_from_index);
    fq.def("from_coeffs", &OpsFq::from_coeffs, py::arg("coeffs"), docs::fq_from_coeffs);
    fq.def("to_coeffs", &OpsFq::to_coeffs, py::arg("i"), docs::fq_to_coeffs);
    registrar_cuerpo(fq);

    registrar_poli<FpW>(m);
    registrar_poli<FqW>(m);
    registrar_rapidos<FpW>(m);
    registrar_rapidos<FqW>(m);
}
