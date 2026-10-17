// test_robustez.cpp — Validaciones, errores y añadidos de la fase 1.
//
// Cubre:
//   1. Tipos EXACTOS de excepción: alcp::division_by_zero (invertir 0, dividir
//      por el polinomio cero, Toeplitz singular) frente a std::invalid_argument
//      (argumentos inválidos) y std::domain_error (resto).
//   2. Comprobación de cuerpo común en todas las operaciones binarias de Poly
//      y same_field / operator== de Fp y Fq.
//   3. prime_power_factor por detección de potencia perfecta (bases grandes).
//   4. Construcciones inválidas de Fq (módulo vacío, n = 0, n absurdo, ...).
//   5. Precondiciones de toeplitz_*, dft/idft y fft::tables/transform.
//   6. multiplicity / radical / distinct_degree / equal_degree con entradas
//      límite (constantes, cero) sin bucles infinitos.
//   7. Orden canónico de factor().
//   8. evaluate, compose, is_squarefree, coeff, Fq::from_vector.
//   9. Concepto FiniteField y prohibición de cuerpos temporales en Poly.
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <chrono>
#include <climits>
#include <cstdio>
#include <exception>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

using namespace alcp;

// ----------------------------------------------------------------------------
// Utilidades
// ----------------------------------------------------------------------------

// ¿f() lanza una excepción cuyo tipo dinámico es EXACTAMENTE Ex?  (Así se
// distingue division_by_zero de su base std::domain_error.)
template <class Ex, class Fn>
static bool throws_exactly(Fn&& f) {
    try {
        f();
    } catch (const std::exception& e) {
        if (typeid(e) == typeid(Ex)) return true;
        std::printf("    (se esperaba %s y se lanzó %s: %s)\n", typeid(Ex).name(), typeid(e).name(), e.what());
        return false;
    } catch (...) {
        std::printf("    (se lanzó una excepción que no deriva de std::exception)\n");
        return false;
    }
    std::printf("    (no se lanzó ninguna excepción; se esperaba %s)\n", typeid(Ex).name());
    return false;
}

#define CHECK_THROWS(Ex, expr) CHECK_MSG(throws_exactly<Ex>([&] { (void)(expr); }), "lanza " #Ex ": " #expr)

using P = Poly<Fp>;
using PQ = Poly<Fq>;

static double ms_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

// ----------------------------------------------------------------------------
// 1. Tipos de excepción
// ----------------------------------------------------------------------------
static void test_tipos_de_error() {
    // division_by_zero deriva de std::domain_error (compatibilidad hacia atrás).
    static_assert(std::is_base_of_v<std::domain_error, division_by_zero>);
    {
        bool capturada = false;
        try { throw division_by_zero("x"); } catch (const std::domain_error&) { capturada = true; }
        CHECK(capturada);
    }

    Fp F(7);
    CHECK_THROWS(division_by_zero, F.inv(F.zero()));
    CHECK_THROWS(division_by_zero, F.div(F.one(), F.zero()));
    CHECK_THROWS(division_by_zero, F.pow(F.zero(), BigInt(-1)));
    CHECK_THROWS(division_by_zero, F.pow(F.zero(), BigInt(-5)));
    CHECK(F.pow(F.zero(), BigInt(0)) == 1);  // 0^0 = 1 (convención de powmod)
    CHECK(F.pow(F.from_int(3), BigInt(-1)) == 5);  // 3·5 = 15 = 1 (mod 7)

    Fq F9(Fp(3), std::vector<BigInt>{1, 0, 1});
    CHECK_THROWS(division_by_zero, F9.inv(F9.zero()));
    CHECK_THROWS(division_by_zero, F9.div(F9.one(), F9.zero()));
    CHECK_THROWS(division_by_zero, F9.pow(F9.zero(), BigInt(-2)));

    P f = P::from_ints(F, {1, 2, 3});
    P z = P::zero(F);
    CHECK_THROWS(division_by_zero, f.divmod(z));
    CHECK_THROWS(division_by_zero, f / z);
    CHECK_THROWS(division_by_zero, f % z);
    CHECK_THROWS(division_by_zero, z % z);
    CHECK_THROWS(division_by_zero, fast::toeplitz_divmod(F, f, z));
    CHECK_THROWS(division_by_zero, P::pow_mod(f, BigInt(3), z));
    PQ fq = PQ::variable(F9);
    CHECK_THROWS(division_by_zero, fq.divmod(PQ::zero(F9)));
    CHECK_THROWS(division_by_zero, fast::toeplitz_divmod(F9, fq, PQ::zero(F9)));

    // Toeplitz triangular singular (a[0] == 0): invertir 0.
    std::vector<BigInt> sing{0, 1, 2};
    CHECK_THROWS(division_by_zero, fast::toeplitz_lower_inverse(F, 3, sing));
    CHECK_THROWS(division_by_zero, fast::toeplitz_upper_inverse(F, 3, sing));

    // std::domain_error "puro" (no division_by_zero).
    P x2 = P::from_ints(F, {0, 0, 1});      // x^2
    P x = P::variable(F);
    CHECK_THROWS(std::domain_error, P::inv_mod(x, x2));  // gcd(x, x^2) = x
    CHECK_THROWS(std::domain_error, P::pow_mod(x, BigInt(-1), x2));
    CHECK_THROWS(std::domain_error, z.factor());
    CHECK_THROWS(std::domain_error, f.pow(BigInt(-1)));
    CHECK_THROWS(std::domain_error, z.lead());
    CHECK_THROWS(std::domain_error, invmod(BigInt(6), BigInt(9)));

    // std::invalid_argument: argumentos inválidos.
    CHECK_THROWS(std::invalid_argument, Fp(BigInt(1)));
    CHECK_THROWS(std::invalid_argument, Fp(BigInt(9)));
    CHECK_THROWS(std::invalid_argument, Fp(BigInt(-7)));
    CHECK_THROWS(std::invalid_argument, powmod(BigInt(2), BigInt(3), BigInt(0)));
    CHECK_THROWS(std::invalid_argument, powmod(BigInt(2), BigInt(3), BigInt(-5)));
    CHECK_THROWS(std::invalid_argument, powmod(BigInt(2), BigInt(-1), BigInt(5)));
    CHECK_THROWS(std::invalid_argument, invmod(BigInt(2), BigInt(0)));
    CHECK_THROWS(std::invalid_argument, rand_below(BigInt(0)));
    CHECK_THROWS(std::invalid_argument, v2(BigInt(0)));
    CHECK_THROWS(std::overflow_error, to_ulong(BigInt(1) << 80));
}

// ----------------------------------------------------------------------------
// 2. Mismo cuerpo
// ----------------------------------------------------------------------------
static void test_mismo_cuerpo() {
    Fp A(7), B(7), C(11);
    CHECK(same_field(A, A));
    CHECK(same_field(A, B));
    CHECK(!same_field(A, C));
    CHECK(A == B);
    CHECK(!(A == C));
    CHECK(A != C);

    // Polinomios sobre objetos de cuerpo distintos pero iguales: operan.
    P fa = P::from_ints(A, {1, 2, 3});
    P fb = P::from_ints(B, {1, 2, 3});
    P gb = P::from_ints(B, {5, 1});
    CHECK(fa == fb);
    CHECK(!(fa != fb));
    CHECK(fa + gb == P::from_ints(A, {6, 3, 3}));
    CHECK(fa - gb == P::from_ints(A, {3, 1, 3}));
    CHECK((fa * gb).degree() == 3);
    CHECK(fa.mul(gb, MulAlgo::Karatsuba) == fa.mul(gb, MulAlgo::Schoolbook));
    CHECK(fa.divmod(gb).first * gb + fa.divmod(gb).second == fa);
    CHECK(P::gcd(fa, gb).degree() >= 0);
    CHECK(fa.multiplicity(gb) >= 0);

    // Cuerpos distintos: std::invalid_argument en toda operación binaria.
    P fc = P::from_ints(C, {1, 2, 3});
    P gc = P::from_ints(C, {5, 1});
    CHECK_THROWS(std::invalid_argument, fa + gc);
    CHECK_THROWS(std::invalid_argument, fa - gc);
    CHECK_THROWS(std::invalid_argument, fa * gc);
    CHECK_THROWS(std::invalid_argument, fa.mul(gc, MulAlgo::Auto));
    CHECK_THROWS(std::invalid_argument, fa.mul(gc, MulAlgo::Schoolbook));
    CHECK_THROWS(std::invalid_argument, fa.mul(gc, MulAlgo::Karatsuba));
    CHECK_THROWS(std::invalid_argument, fa.mul(gc, MulAlgo::FFT));
    CHECK_THROWS(std::invalid_argument, fa.mul_school(gc));
    CHECK_THROWS(std::invalid_argument, fa.mul_karatsuba(gc));
    CHECK_THROWS(std::invalid_argument, fa.divmod(gc));
    CHECK_THROWS(std::invalid_argument, fa / gc);
    CHECK_THROWS(std::invalid_argument, fa % gc);
    CHECK_THROWS(std::invalid_argument, P::gcd(fa, gc));
    CHECK_THROWS(std::invalid_argument, P::gcd_ext(fa, gc));
    CHECK_THROWS(std::invalid_argument, P::inv_mod(fa, gc));
    CHECK_THROWS(std::invalid_argument, P::pow_mod(fa, BigInt(5), gc));
    CHECK_THROWS(std::invalid_argument, fa.multiplicity(gc));
    CHECK_THROWS(std::invalid_argument, fa.compose(gc));
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_divmod(A, fa, gc));
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_divmod(C, fa, gc));
    // La comprobación de cuerpo va antes que la de divisor nulo.
    CHECK_THROWS(std::invalid_argument, fa.divmod(P::zero(C)));
    // operator== no lanza: cuerpos distintos -> distintos (incluso el cero).
    CHECK(!(fa == fc));
    CHECK(fa != fc);
    CHECK(!(P::zero(A) == P::zero(C)));
    CHECK(P::zero(A) == P::zero(B));

    // F_q: mismo p y mismo módulo mónico.
    Fp F3(3), F3b(3);
    Fq K1(F3, std::vector<BigInt>{1, 0, 1});      // x^2 + 1
    Fq K2(F3b, std::vector<BigInt>{2, 0, 2});     // 2x^2 + 2 -> se normaliza a x^2 + 1
    Fq K3(F3, std::vector<BigInt>{2, 1, 1});      // x^2 + x + 2 (otro módulo irreducible)
    Fq K4(Fp(5), std::vector<BigInt>{2, 0, 1});   // x^2 + 2 sobre F_5
    CHECK(K2.modulus() == K1.modulus());
    CHECK(same_field(K1, K1));
    CHECK(same_field(K1, K2));
    CHECK(K1 == K2);
    CHECK(!same_field(K1, K3));
    CHECK(!same_field(K1, K4));
    PQ q1 = PQ::variable(K1), q2 = PQ::variable(K2), q3 = PQ::variable(K3);
    CHECK(q1 == q2);
    CHECK((q1 * q2).degree() == 2);
    CHECK(!(q1 == q3));
    CHECK_THROWS(std::invalid_argument, q1 + q3);
    CHECK_THROWS(std::invalid_argument, q1 * q3);
    CHECK_THROWS(std::invalid_argument, q1.divmod(q3));
    CHECK_THROWS(std::invalid_argument, PQ::gcd(q1, q3));
    CHECK_THROWS(std::invalid_argument, q1.compose(q3));
}

// ----------------------------------------------------------------------------
// 3. prime_power_factor
// ----------------------------------------------------------------------------
static void test_prime_power_factor() {
    const BigInt m61 = (BigInt(1) << 61) - 1;  // primo de Mersenne
    const BigInt m89 = (BigInt(1) << 89) - 1;  // primo de Mersenne
    const BigInt m127 = (BigInt(1) << 127) - 1;
    CHECK(is_prime(m61) && is_prime(m89) && is_prime(m127));

    auto t0 = std::chrono::steady_clock::now();
    auto r1 = prime_power_factor(ipow(m61, 3));
    CHECK(r1.first == m61 && r1.second == 3);
    auto r2 = prime_power_factor(ipow(BigInt(1000003), 5));
    CHECK(r2.first == 1000003 && r2.second == 5);
    auto r3 = prime_power_factor(ipow(m89, 2));
    CHECK(r3.first == m89 && r3.second == 2);
    auto r4 = prime_power_factor(m127);
    CHECK(r4.first == m127 && r4.second == 1);
    auto r5 = prime_power_factor(BigInt(1) << 100);
    CHECK(r5.first == 2 && r5.second == 100);
    auto r6 = prime_power_factor(ipow(BigInt(3), 200));
    CHECK(r6.first == 3 && r6.second == 200);
    auto r7 = prime_power_factor(BigInt(4));
    CHECK(r7.first == 2 && r7.second == 2);
    auto r8 = prime_power_factor(BigInt(2));
    CHECK(r8.first == 2 && r8.second == 1);
    // No son potencias de primo: invalid_argument, también al instante.
    CHECK_THROWS(std::invalid_argument, prime_power_factor(ipow(BigInt(6), 10)));
    CHECK_THROWS(std::invalid_argument, prime_power_factor(m61 * m89));
    CHECK_THROWS(std::invalid_argument, prime_power_factor(ipow(m61, 2) * 3));
    CHECK_THROWS(std::invalid_argument, prime_power_factor(ipow(m61 * m89, 3)));
    CHECK_THROWS(std::invalid_argument, prime_power_factor(BigInt(12)));
    CHECK_THROWS(std::invalid_argument, prime_power_factor(BigInt(1)));
    CHECK_THROWS(std::invalid_argument, prime_power_factor(BigInt(0)));
    CHECK_THROWS(std::invalid_argument, prime_power_factor(BigInt(-8)));
    const double ms = ms_since(t0);
    CHECK_MSG(ms < 1000.0, "prime_power_factor con bases grandes termina al instante (< 1 s)");

    // Fq(q) con q = (2^61-1)^3: antes la división por tentativa no terminaba.
    t0 = std::chrono::steady_clock::now();
    Fq K(ipow(m61, 3));
    CHECK(K.characteristic() == m61);
    CHECK(K.extension_degree() == 3);
    auto a = K.random();
    if (!K.is_zero(a)) CHECK(K.is_one(K.mul(a, K.inv(a))));
    CHECK_MSG(ms_since(t0) < 3000.0, "Fq((2^61-1)^3) se construye en un tiempo razonable");
}

// ----------------------------------------------------------------------------
// 4. Fq inválidos y Fq::from_vector
// ----------------------------------------------------------------------------
static void test_fq_validacion() {
    Fp F3(3);
    CHECK_THROWS(std::invalid_argument, Fq(F3, std::vector<BigInt>{}));       // módulo vacío
    CHECK_THROWS(std::invalid_argument, Fq(F3, std::vector<BigInt>{1}));      // constante
    CHECK_THROWS(std::invalid_argument, Fq(F3, std::vector<BigInt>{1, 0, 0}));  // principal 0
    CHECK_THROWS(std::invalid_argument, Fq(F3, std::vector<BigInt>{1, 0, 3}));  // principal = p
    CHECK_THROWS(std::invalid_argument, Fq(F3, std::vector<BigInt>{2, 0, 1}));  // x^2-1 reducible
    CHECK_THROWS(std::invalid_argument, Fq(BigInt(3), 0UL));                   // n = 0
    CHECK_THROWS(std::invalid_argument, Fq(BigInt(3), 1UL << 20));             // n absurdo
    CHECK_THROWS(std::invalid_argument, Fq(BigInt(3), ULONG_MAX));
    CHECK_THROWS(std::invalid_argument, Fq(BigInt(4), 2UL));                   // 4 no es primo
    CHECK_THROWS(std::invalid_argument, Fq(BigInt(12)));                       // no es p^n
    CHECK_THROWS(std::invalid_argument, Fq(BigInt(1)));
    CHECK_THROWS(std::invalid_argument, Fq(BigInt(BigInt(1) << 100000)));      // n = 100000
    {
        // Módulo con un grado absurdo: se rechaza antes de calcular p^n.
        auto t0 = std::chrono::steady_clock::now();
        std::vector<BigInt> enorme(Fq::kMaxExtensionDegree + 2, BigInt(1));
        CHECK_THROWS(std::invalid_argument, Fq(F3, enorme));
        CHECK_MSG(ms_since(t0) < 500.0, "un módulo de grado absurdo se rechaza al instante");
    }
    // Casos válidos límite.
    Fq K1(BigInt(5), 1UL);  // F_5 visto como extensión de grado 1
    CHECK(K1.cardinality() == 5);
    CHECK(K1.is_one(K1.one()));
    CHECK(K1.is_one(K1.mul(K1.from_int(2), K1.from_int(3))));
    Fq K2(BigInt(2));  // F_2
    CHECK(K2.cardinality() == 2 && K2.extension_degree() == 1);

    // from_vector: reducción generosa.
    Fq F9(F3, std::vector<BigInt>{1, 0, 1});  // F_3[a]/(a^2 + 1), a^2 = -1
    using E = Fq::Elem;
    CHECK(F9.from_vector({}) == F9.zero());
    CHECK(F9.from_vector({1}) == F9.one());
    CHECK(F9.from_vector({-1, 4}) == (E{2, 1}));            // signos y valores >= p
    CHECK(F9.from_vector({0, 0, 1}) == (E{2, 0}));          // a^2 = -1 = 2
    CHECK(F9.from_vector({0, 0, 0, 1}) == (E{0, 2}));       // a^3 = -a
    CHECK(F9.from_vector({5, -7, 0, 0, 1}) == (E{0, 2}));   // 5 - 7a + a^4 = (2 + 1) + 2a = 2a
    CHECK(F9.from_vector({BigInt("123456789012345678901234567890")}) ==
          F9.from_int(BigInt("123456789012345678901234567890")));
    // Comparación con Horner en el propio cuerpo sobre vectores aleatorios.
    set_seed(777);
    Fq F125(BigInt(5), 3UL);
    const E a = F125.from_index(BigInt(5));  // la clase de x (índice p^1)
    for (int t = 0; t < 40; ++t) {
        std::vector<BigInt> v(1 + t % 9);
        for (auto& c : v) c = rand_below(BigInt(1000)) - 500;
        E esperado = F125.zero();
        for (std::size_t i = v.size(); i-- > 0;)
            esperado = F125.add(F125.mul(esperado, a), F125.from_int(v[i]));
        CHECK(F125.from_vector(v) == esperado);
    }
}

// ----------------------------------------------------------------------------
// 5. Toeplitz, DFT y FFT inválidos
// ----------------------------------------------------------------------------
static void test_rapidos_validacion() {
    Fp F(17);
    using V = std::vector<BigInt>;
    V a3{1, 2, 3}, b3{4, 5, 6}, b2{4, 5};
    // toeplitz_lower_vec / upper_vec
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_lower_vec(F, 0, V{}, V{}));
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_lower_vec(F, 3, V{1, 2}, b3));
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_lower_vec(F, 3, a3, b2));
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_lower_vec(F, 3, a3, V{1, 2, 3, 4}));
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_upper_vec(F, 0, V{}, V{}));
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_upper_vec(F, 3, V{1}, b3));
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_upper_vec(F, 3, a3, b2));
    // a más largo que n: se usan solo a[0..n-1].
    CHECK(fast::toeplitz_lower_vec(F, 3, V{1, 2, 3, 9, 9}, b3) == fast::toeplitz_lower_vec(F, 3, a3, b3));
    CHECK(fast::toeplitz_upper_vec(F, 3, V{1, 2, 3, 9}, b3) == fast::toeplitz_upper_vec(F, 3, a3, b3));
    CHECK(fast::toeplitz_lower_vec(F, 1, V{3}, V{5}) == V{15});
    // toeplitz_vec: a.size() == 2n-1
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_vec(F, 0, V{}, V{}));
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_vec(F, 3, a3, b3));
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_vec(F, 3, V{1, 2, 3, 4, 5, 6}, b3));
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_vec(F, 3, V{1, 2, 3, 4, 5}, b2));
    CHECK(fast::toeplitz_vec(F, 1, V{3}, V{5}) == V{15});
    // inversas triangulares: n = 0 (antes recursión infinita), a corto
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_lower_inverse(F, 0, V{1}));
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_lower_inverse(F, 0, V{}));
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_lower_inverse(F, 4, a3));
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_upper_inverse(F, 0, V{}));
    CHECK_THROWS(std::invalid_argument, fast::toeplitz_upper_inverse(F, 4, a3));
    CHECK(fast::toeplitz_lower_inverse(F, 1, V{2}) == V{9});  // 2·9 = 18 = 1

    // dft / idft: 3 tiene orden 16 en F_17; 2 tiene orden 8; 16 = -1 orden 2.
    const BigInt g16 = 3, g8 = 2, g2 = 16;
    V a16(16), a8(8), a1{7};
    for (int i = 0; i < 16; ++i) a16[i] = i;
    for (int i = 0; i < 8; ++i) a8[i] = i;
    CHECK_THROWS(std::invalid_argument, fast::dft(F, g16, 63, a16));        // k >= 63
    CHECK_THROWS(std::invalid_argument, fast::dft(F, g16, 200, a16));
    CHECK_THROWS(std::invalid_argument, fast::dft(F, g16, 4, a8));          // longitud != 2^k
    CHECK_THROWS(std::invalid_argument, fast::dft(F, g8, 4, a16));          // orden 8 != 16
    CHECK_THROWS(std::invalid_argument, fast::dft(F, F.one(), 2, V{1, 2, 3, 4}));
    CHECK_THROWS(std::invalid_argument, fast::dft(F, F.zero(), 1, V{1, 2}));
    CHECK_THROWS(std::invalid_argument, fast::dft(F, g16, 3, a8));          // 3 tiene orden 16, no 8
    CHECK_THROWS(std::invalid_argument, fast::dft(F, g16, 0, a1));          // k = 0 exige g = 1
    CHECK_THROWS(std::invalid_argument, fast::idft(F, g16, 63, a16));
    CHECK_THROWS(std::invalid_argument, fast::idft(F, g16, 4, a8));
    CHECK_THROWS(std::invalid_argument, fast::idft(F, g8, 4, a16));
    CHECK_THROWS(std::invalid_argument, fast::idft(F, F.zero(), 1, V{1, 2}));
    CHECK(fast::dft(F, F.one(), 0, a1) == a1);
    CHECK(fast::idft(F, F.one(), 0, a1) == a1);
    CHECK(fast::dft(F, g2, 1, V{1, 2}) == (V{3, 16}));      // (1+2, 1-2)
    CHECK(fast::idft(F, g16, 4, fast::dft(F, g16, 4, a16)) == a16);
    CHECK(fast::idft(F, g8, 3, fast::dft(F, g8, 3, a8)) == a8);
    // En característica 2 no hay raíces de orden 2^k con k >= 1.
    Fq F4(BigInt(2), 2UL);
    std::vector<Fq::Elem> c2{F4.one(), F4.one()};
    CHECK_THROWS(std::invalid_argument, fast::dft(F4, F4.one(), 1, c2));

    // fft::tables / transform: capacidad de F_7 es 1 (7 - 1 = 2·3).
    Fp F7(7);
    CHECK(F7.fft_capacity() == 1);
    CHECK_THROWS(std::invalid_argument, fft::tables(F7, 2));
    CHECK_THROWS(std::invalid_argument, fft::tables(F7, 70));
    CHECK(fft::tables(F7, 1).fwd.size() == 1);
    V v4{1, 2, 3, 4}, v3{1, 2, 3};
    V copia = v4;
    CHECK_THROWS(std::invalid_argument, fft::transform(F7, v4, false));  // 2^2 > capacidad
    CHECK(v4 == copia);  // la excepción deja el vector intacto
    CHECK_THROWS(std::invalid_argument, fft::transform(F, v3, false));   // no potencia de 2
    CHECK(v3 == (V{1, 2, 3}));
    CHECK_THROWS(std::invalid_argument, fft::convolve(F7, V{1, 2, 3}, V{1, 2}));
    CHECK(fft::fft_exponent(F7, 4) == -1);
    // F_2: capacidad 0, solo la transformada trivial.
    Fp F2(2);
    CHECK(F2.fft_capacity() == 0);
    CHECK_THROWS(std::invalid_argument, fft::tables(F2, 1));
    V uno{1};
    fft::transform(F2, uno, false);
    CHECK(uno == V{1});
    // Ida y vuelta válida en F_17 (2^4 | 16).
    V w(8);
    for (int i = 0; i < 8; ++i) w[i] = (i * 5 + 3) % 17;
    V w0 = w;
    fft::transform(F, w, false);
    fft::transform(F, w, true);
    CHECK(w == w0);
}

// ----------------------------------------------------------------------------
// 6. multiplicity, radical, distinct_degree, equal_degree en casos límite
// ----------------------------------------------------------------------------
static void test_casos_limite_factorizacion() {
    Fp F(7);
    P x = P::variable(F);
    P c5 = P::constant(F, F.from_int(5));
    P z = P::zero(F);
    P u1 = P::from_ints(F, {1, 1});     // x + 1
    P u2 = P::from_ints(F, {2, 1});     // x + 2
    P q2 = P::from_ints(F, {1, 0, 1});  // x^2 + 1 (irreducible en F_7)
    P f = u1.pow(BigInt(3)) * u2 * q2.pow(BigInt(2));

    // multiplicity
    CHECK(f.multiplicity(u1) == 3);
    CHECK(f.multiplicity(u2) == 1);
    CHECK(f.multiplicity(q2) == 2);
    CHECK(f.multiplicity(x) == 0);
    CHECK(f.multiplicity(u1.scalar_mul(F.from_int(3))) == 3);  // no hace falta que sea mónico
    CHECK(c5.multiplicity(u1) == 0);
    CHECK_THROWS(std::invalid_argument, f.multiplicity(c5));  // antes: bucle infinito
    CHECK_THROWS(std::invalid_argument, f.multiplicity(P::one(F)));
    CHECK_THROWS(std::invalid_argument, f.multiplicity(z));
    CHECK_THROWS(std::invalid_argument, z.multiplicity(u1));  // multiplicidad infinita

    // radical
    CHECK(c5.radical() == P::one(F));                         // antes: recursión infinita
    CHECK(P::one(F).radical() == P::one(F));
    CHECK(z.radical() == z);
    CHECK(f.radical() == (u1 * u2 * q2).monic());
    CHECK(x.pow(BigInt(7)).radical() == x);                   // f' = 0
    CHECK((x.pow(BigInt(14)) * u2).radical() == x * u2);      // rama del residuo p-ésimo
    CHECK(f.scalar_mul(F.from_int(4)).radical() == (u1 * u2 * q2).monic());

    // Característica que no cabe en unsigned long de 32 bits (Windows):
    // radical/factor siguen funcionando.
    const BigInt pg("4294967311");  // 2^32 + 15, primo
    CHECK(is_prime(pg));
    Fp G(pg);
    Poly<Fp> cg = Poly<Fp>::constant(G, G.from_int(12345));
    CHECK(cg.radical() == Poly<Fp>::one(G));
    Poly<Fp> lg = Poly<Fp>::from_ints(G, {3, 1});
    Poly<Fp> mg = Poly<Fp>::from_ints(G, {5, 1});
    Poly<Fp> fg = lg.pow(BigInt(2)) * mg;
    CHECK(fg.radical() == lg * mg);
    auto facg = fg.factor();
    CHECK(facg.factors.size() == 2);
    if (facg.factors.size() == 2) {
        CHECK(facg.factors[0].first == lg && facg.factors[0].second == 2);
        CHECK(facg.factors[1].first == mg && facg.factors[1].second == 1);
    }

    // distinct_degree / equal_degree: precondiciones
    CHECK_THROWS(std::invalid_argument, z.distinct_degree());
    CHECK_THROWS(std::invalid_argument, c5.distinct_degree());
    CHECK_THROWS(std::invalid_argument, z.equal_degree(1));
    CHECK_THROWS(std::invalid_argument, c5.equal_degree(1));
    CHECK_THROWS(std::invalid_argument, q2.equal_degree(0));
    CHECK_THROWS(std::invalid_argument, q2.equal_degree(-3));
    CHECK_THROWS(std::invalid_argument, (u1 * q2).equal_degree(2));  // grado 3 no múltiplo de 2
    // Irreducible de grado 2 tratado como producto de lineales: no hay bucle
    // infinito, se lanza invalid_argument tras agotar los intentos.
    CHECK_THROWS(std::invalid_argument, q2.equal_degree(1));
    CHECK_THROWS(std::invalid_argument, (q2 * P::from_ints(F, {3, 0, 1})).equal_degree(1));
    // Casos válidos.
    auto dd = (u1 * u2 * q2).distinct_degree();
    CHECK(dd.size() == 2);
    if (dd.size() == 2) {
        CHECK(dd[0].first == 1 && dd[0].second == u1 * u2);
        CHECK(dd[1].first == 2 && dd[1].second == q2);
    }
    auto ed = (u1 * u2).equal_degree(1);
    CHECK(ed.size() == 2);
    auto dd1 = u1.distinct_degree();
    CHECK(dd1.size() == 1 && dd1[0].first == 1 && dd1[0].second == u1);
    CHECK(q2.equal_degree(2).size() == 1);

    // factor de constantes y polinomios de grado 0
    auto fc = c5.factor();
    CHECK(fc.factors.empty());
    CHECK(fc.leading == 5);
}

// ----------------------------------------------------------------------------
// 7. Orden canónico de factor()
// ----------------------------------------------------------------------------
template <class Field>
static bool factores_en_orden(const typename Poly<Field>::Factorization& fac) {
    for (std::size_t i = 1; i < fac.factors.size(); ++i)
        if (!Poly<Field>::canonical_less(fac.factors[i - 1].first, fac.factors[i].first)) return false;
    return true;
}

static void test_orden_canonico() {
    Fp F(7);
    P l1 = P::from_ints(F, {1, 1}), l3 = P::from_ints(F, {3, 1}), l5 = P::from_ints(F, {5, 1});
    P q1 = P::from_ints(F, {1, 0, 1});   // x^2 + 1
    P q3 = P::from_ints(F, {3, 1, 1});   // x^2 + x + 3 (disc. 1 - 12 = -11 = 3: no residuo)
    CHECK(q1.is_irreducible() && q3.is_irreducible());
    P f = q3 * l5 * q1.pow(BigInt(2)) * l1 * l3.pow(BigInt(3));
    f = f.scalar_mul(F.from_int(4));
    // Orden esperado: grado 1 por término independiente (1, 3, 5); grado 2:
    // x^2 + 0x + 1 < x^2 + 1x + 3 (se compara primero el coeficiente de x).
    std::vector<std::pair<P, int>> esperado{{l1, 1}, {l3, 3}, {l5, 1}, {q1, 2}, {q3, 1}};
    for (unsigned long semilla : {1UL, 2UL, 3UL, 99UL}) {
        set_seed(semilla);  // el orden no depende de la aleatoriedad del algoritmo
        auto fac = f.factor();
        CHECK(fac.leading == 4);
        CHECK(fac.factors.size() == esperado.size());
        if (fac.factors.size() == esperado.size())
            for (std::size_t i = 0; i < esperado.size(); ++i) {
                CHECK(fac.factors[i].first == esperado[i].first);
                CHECK(fac.factors[i].second == esperado[i].second);
            }
        CHECK(factores_en_orden<Fp>(fac));
    }
    CHECK(P::canonical_less(l5, q1));
    CHECK(!P::canonical_less(q1, l5));
    CHECK(!P::canonical_less(l1, l1));

    // F_q: orden canónico en factorizaciones aleatorias.
    Fq F9(Fp(3), std::vector<BigInt>{1, 0, 1});
    set_seed(4242);
    for (int t = 0; t < 15; ++t) {
        PQ g = PQ::one(F9);
        for (int i = 0; i < 4; ++i) g = g * PQ::random(F9, 1 + t % 3).monic();
        auto fac = g.factor();
        CHECK(factores_en_orden<Fq>(fac));
    }
}

// ----------------------------------------------------------------------------
// 8. evaluate, compose, is_squarefree, coeff
// ----------------------------------------------------------------------------
static void test_evaluar_componer() {
    Fp F(7);
    P f = P::from_ints(F, {1, 2, 3});  // 3x^2 + 2x + 1
    CHECK(f.evaluate(F.from_int(2)) == 3);     // 12 + 4 + 1 = 17 = 3
    CHECK(f.evaluate(F.zero()) == 1);
    CHECK(P::zero(F).evaluate(F.from_int(4)) == 0);
    CHECK(P::constant(F, F.from_int(6)).evaluate(F.from_int(3)) == 6);
    // Evaluación por fuerza bruta en todo F_7.
    for (long v = 0; v < 7; ++v) {
        auto xv = F.from_int(v);
        BigInt suma = 0;
        for (std::size_t i = 0; i < f.coeffs().size(); ++i)
            suma += f.coeffs()[i] * F.pow(xv, BigInt(static_cast<unsigned long>(i)));
        CHECK(f.evaluate(xv) == F.from_int(suma));
    }

    // compose
    P x = P::variable(F);
    P g = P::from_ints(F, {4, 0, 1});  // x^2 + 4
    P fg = f.compose(g);
    CHECK(fg.degree() == 4);
    CHECK(x.compose(g) == g);
    CHECK(f.compose(x) == f);
    CHECK(P::zero(F).compose(g).is_zero());
    CHECK(f.compose(P::zero(F)) == P::constant(F, F.from_int(1)));
    CHECK(f.compose(P::constant(F, F.from_int(2))) == P::constant(F, F.from_int(3)));
    // f(g(x)) = 3(x^2+4)^2 + 2(x^2+4) + 1
    CHECK(fg == P::from_ints(F, {3 * 16 + 8 + 1, 0, 24 + 2, 0, 3}));
    for (long v = 0; v < 7; ++v) {
        auto xv = F.from_int(v);
        CHECK(fg.evaluate(xv) == f.evaluate(g.evaluate(xv)));
    }
    // compose con un cuerpo distinto pero igual (otro objeto con p = 7)
    Fp F7b(7);
    CHECK(f.compose(P::from_ints(F7b, {4, 0, 1})) == fg);

    // is_squarefree
    CHECK(!P::zero(F).is_squarefree());
    CHECK(P::constant(F, F.from_int(3)).is_squarefree());
    CHECK(x.is_squarefree());
    CHECK(P::from_ints(F, {1, 0, 1}).is_squarefree());
    CHECK(!P::from_ints(F, {1, 2, 1}).is_squarefree());  // (x+1)^2
    CHECK(!(x.pow(BigInt(7)) + P::one(F)).is_squarefree());  // (x+1)^7, f' = 0
    CHECK(!x.pow(BigInt(7)).is_squarefree());
    CHECK((x.pow(BigInt(7)) - x).is_squarefree());       // producto de todos los x - c
    CHECK(!(f * f).is_squarefree());

    // coeff
    CHECK(f.coeff(0) == 1 && f.coeff(1) == 2 && f.coeff(2) == 3);
    CHECK(f.coeff(3) == 0 && f.coeff(1000) == 0);
    CHECK(P::zero(F).coeff(0) == 0);
    Fq F9(Fp(3), std::vector<BigInt>{1, 0, 1});
    PQ h = PQ::variable(F9);
    CHECK(h.coeff(1) == F9.one());
    CHECK(h.coeff(5) == F9.zero());
    CHECK(PQ::zero(F9).coeff(0) == F9.zero());
    // evaluate / compose sobre F_9
    for (long i = 0; i < 9; ++i) {
        auto e = F9.from_index(BigInt(i));
        CHECK(h.evaluate(e) == e);
        PQ hh = h * h + PQ::constant(F9, F9.one());  // x^2 + 1
        CHECK(hh.evaluate(e) == F9.add(F9.mul(e, e), F9.one()));
        CHECK(hh.compose(hh).evaluate(e) == hh.evaluate(hh.evaluate(e)));
    }
}

// ----------------------------------------------------------------------------
// 9. Concepto FiniteField y cuerpos temporales
// ----------------------------------------------------------------------------
struct CuerpoIncompleto {  // le faltan casi todos los miembros
    using Elem = int;
    int zero() const { return 0; }
};

template <class PolyT, class Arg>
concept puede_construir_cero = requires { PolyT::zero(std::declval<Arg>()); };
template <class PolyT, class Arg>
concept puede_from_ints = requires { PolyT::from_ints(std::declval<Arg>(), std::vector<BigInt>{}); };
template <class PolyT, class Arg>
concept puede_variable = requires { PolyT::variable(std::declval<Arg>()); };

static void test_concepto_y_temporales() {
    static_assert(FiniteField<Fp>);
    static_assert(FiniteField<Fq>);
    static_assert(!FiniteField<int>);
    static_assert(!FiniteField<CuerpoIncompleto>);
    static_assert(!FiniteField<BigInt>);

    // Poly no admite un cuerpo temporal (el puntero quedaría colgando)...
    static_assert(!std::is_constructible_v<P, Fp>);
    static_assert(!std::is_constructible_v<P, Fp&&>);
    static_assert(!std::is_constructible_v<P, const Fp&&>);
    static_assert(!std::is_constructible_v<P, Fp, std::vector<BigInt>>);
    static_assert(!puede_construir_cero<P, Fp>);
    static_assert(!puede_from_ints<P, Fp>);
    static_assert(!puede_variable<PQ, Fq>);
    // ... pero sí uno con nombre (lvalue), constante o no.
    static_assert(std::is_constructible_v<P, Fp&>);
    static_assert(std::is_constructible_v<P, const Fp&>);
    static_assert(std::is_constructible_v<P, const Fp&, std::vector<BigInt>>);
    static_assert(puede_construir_cero<P, Fp&>);
    static_assert(puede_from_ints<P, const Fp&>);
    static_assert(puede_variable<PQ, Fq&>);
    CHECK(true);  // las comprobaciones anteriores son en compilación
}

int main() {
    set_seed(20261017);
    auto t0 = std::chrono::steady_clock::now();

    test_tipos_de_error();
    test_mismo_cuerpo();
    test_prime_power_factor();
    test_fq_validacion();
    test_rapidos_validacion();
    test_casos_limite_factorizacion();
    test_orden_canonico();
    test_evaluar_componer();
    test_concepto_y_temporales();

    std::printf("test_robustez: %.0f ms\n", ms_since(t0));
    return check::summary("robustez");
}
