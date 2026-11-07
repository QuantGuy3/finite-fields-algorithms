// test_fase4.cpp — Fase 4 (nivel A): algoritmos rápidos que ya estaban en la librería.
//
// Cada técnica nueva debe dar EXACTAMENTE el mismo resultado que la original:
//   1. inverse_series: a·a^{-1} = 1 (mod x^n) y error si a[0] = 0.
//   2. divmod_newton frente a divmod_school (y el despacho de divmod) en una
//      rejilla de grados, en F_p (GMP, p pequeño y grande), Fp64, F_q y Fq64,
//      con y sin tablas FFT ya construidas (los dos caminos de mul_no_rng).
//   3. pow_mod con el Reducer (Newton) frente a una potencia con la división
//      escolar.
//   4. Inverso de F_q por Euclides extendido frente a Fermat (inv_fermat):
//      todos los elementos de F_9, F_256 y F_(7^4); muestras de F_(3^30) y de
//      F_((2^61-1)^2), con los dos backends.
//   5. is_irreducible sin la comprobación final redundante frente al test
//      original (con ella) y frente al número de irreducibles mónicos de
//      Gauss: (1/d)·Σ_{e|d} μ(e)·q^{d/e}.
//   6. Generador aleatorio: la división de Newton no calcula la raíz de la
//      FFT ni consume números aleatorios (la secuencia es la del original).
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <cstdio>
#include <exception>
#include <typeinfo>
#include <vector>

using namespace alcp;

namespace {

template <class Ex, class Fn>
bool throws_exactly(Fn&& f) {
    try {
        f();
    } catch (const std::exception& e) {
        return typeid(e) == typeid(Ex);
    }
    return false;
}
#define CHECK_THROWS(Ex, expr) CHECK_MSG(throws_exactly<Ex>([&] { (void)(expr); }), "lanza " #Ex ": " #expr)

// ---------------------------------------------------------------- 1. series
template <class Field>
void test_inverse_series(const Field& K, const char* nombre) {
    using P = Poly<Field>;
    for (std::size_t n = 1; n <= 70; n += (n < 10 ? 1 : 7)) {
        P a = P::random(K, static_cast<int>(n) + 3);
        std::vector<typename Field::Elem> av = a.coeffs();
        if (K.is_zero(av[0])) av[0] = K.one();
        P inv(K, P::inverse_series(K, av, n));
        P prod = P(K, av) * inv;
        bool ok = !prod.is_zero() && K.is_one(prod.coeffs()[0]);
        for (std::size_t i = 1; i < n && ok; ++i) ok = K.is_zero(prod.coeff(i));
        CHECK_MSG(ok, nombre);
    }
    std::vector<typename Field::Elem> cero{K.zero(), K.one()};
    CHECK_THROWS(division_by_zero, P::inverse_series(K, cero, 3));
}

// ---------------------------------------------------------------- 2. divmod
template <class Field>
int test_divmod(const Field& K, const std::vector<int>& ns, const std::vector<int>& ds, const char* nombre) {
    using P = Poly<Field>;
    int casos = 0;
    for (int n : ns)
        for (int d : ds) {
            P g = P::random(K, n);
            P f = P::random(K, n + d);
            auto [qs, rs] = f.divmod_school(g);
            auto [qn, rn] = f.divmod_newton(g);
            auto [qa, ra] = f.divmod(g);
            CHECK_MSG(qs == qn && rs == rn, nombre);
            CHECK_MSG(qs == qa && rs == ra, nombre);
            CHECK_MSG(rn.degree() < g.degree(), nombre);
            ++casos;
        }
    // Dividendo de grado menor que el divisor y divisor nulo.
    P g = P::random(K, 10), f = P::random(K, 4);
    CHECK(f.divmod_newton(g).first.is_zero() && f.divmod_newton(g).second == f);
    CHECK_THROWS(division_by_zero, f.divmod_newton(P::zero(K)));
    return casos;
}

// ---------------------------------------------------------------- 3. pow_mod
template <class Field>
Poly<Field> pow_mod_escolar(Poly<Field> a, BigInt k, const Poly<Field>& m) {
    using P = Poly<Field>;
    P result = P::one(a.field()).divmod_school(m).second;
    a = a.divmod_school(m).second;
    while (k > 0) {
        if (mpz_odd_p(k.get_mpz_t())) result = (result * a).divmod_school(m).second;
        k >>= 1;
        if (k > 0) a = (a * a).divmod_school(m).second;
    }
    return result;
}

template <class Field>
void test_pow_mod(const Field& K, int grado, const BigInt& k, const char* nombre) {
    using P = Poly<Field>;
    P m = P::random(K, grado);
    P a = P::random(K, grado + 5);
    CHECK_MSG(P::pow_mod(a, k, m) == pow_mod_escolar(a, k, m), nombre);
    typename P::Reducer red(m);
    CHECK_MSG(P::pow_mod(a, k, red) == pow_mod_escolar(a, k, m), nombre);
}

// ---------------------------------------------------------------- 4. inverso F_q
template <class FqF>
void test_inv_todos(const FqF& K, const char* nombre) {
    const BigInt q = K.cardinality();
    bool ok = true;
    for (BigInt i = 1; i < q; ++i) {
        auto a = K.from_index(i);
        auto x = K.inv(a);
        ok = ok && x == K.inv_fermat(a) && K.is_one(K.mul(a, x));
    }
    CHECK_MSG(ok, nombre);
    CHECK_THROWS(division_by_zero, K.inv(K.zero()));
}

template <class FqF>
void test_inv_muestras(const FqF& K, int veces, const char* nombre) {
    bool ok = true;
    for (int t = 0; t < veces; ++t) {
        auto a = K.random();
        if (K.is_zero(a)) continue;
        auto x = K.inv(a);
        ok = ok && x == K.inv_fermat(a) && K.is_one(K.mul(a, x));
    }
    // Elementos del subcuerpo primo (a de grado 0) y el generador t.
    auto uno = K.one();
    ok = ok && K.inv(uno) == uno;
    auto t = K.from_index(K.characteristic());
    ok = ok && K.is_one(K.mul(t, K.inv(t)));
    CHECK_MSG(ok, nombre);
}

// ---------------------------------------------------------------- 5. irreducibilidad
// El test de Ben-Or ORIGINAL, con la comprobación final x^{q^n} == x (mod f).
template <class Field>
bool ben_or_original(const Poly<Field>& f) {
    using P = Poly<Field>;
    const int n = f.degree();
    if (n <= 0) return false;
    const Field& K = f.field();
    const BigInt q = K.cardinality();
    P x = P::variable(K);
    P H = x % f;
    for (int i = 1; i <= n / 2; ++i) {
        H = P::pow_mod(H, q, f);
        if (P::gcd(f, H - x).degree() != 0) return false;
    }
    P Hn = P::pow_mod(x, ipow(q, static_cast<unsigned long>(n)), f);
    return (Hn - x) % f == P::zero(K);
}

int mobius(int n) {
    int r = 1;
    for (int p = 2; p * p <= n; ++p)
        if (n % p == 0) {
            n /= p;
            if (n % p == 0) return 0;
            r = -r;
        }
    return n > 1 ? -r : r;
}

// Número de polinomios mónicos irreducibles de grado d sobre F_q (Gauss).
BigInt gauss(const BigInt& q, int d) {
    BigInt s = 0;
    for (int e = 1; e <= d; ++e)
        if (d % e == 0) s += mobius(e) * ipow(q, static_cast<unsigned long>(d / e));
    return s / d;
}

// Elemento de índice i: from_index en las extensiones, from_int en F_p.
template <class Field>
typename Field::Elem elem_de_indice(const Field& K, const BigInt& i) {
    if constexpr (requires { K.from_index(i); }) return K.from_index(i);
    else return K.from_int(i);
}

// Recorre todos los mónicos de grado d (coeficientes por índice) y cuenta.
template <class Field>
void test_irreducibles_gauss(const Field& K, int d, const char* nombre) {
    using P = Poly<Field>;
    const BigInt q = K.cardinality();
    const BigInt total = ipow(q, static_cast<unsigned long>(d));
    BigInt cuenta = 0;
    bool coinciden = true;
    for (BigInt idx = 0; idx < total; ++idx) {
        std::vector<typename Field::Elem> c;
        BigInt t = idx;
        for (int i = 0; i < d; ++i) {
            c.push_back(elem_de_indice(K, t % q));
            t /= q;
        }
        c.push_back(K.one());
        P f(K, c);
        const bool nuevo = f.is_irreducible();
        coinciden = coinciden && nuevo == ben_or_original(f);
        if (nuevo) ++cuenta;
    }
    CHECK_MSG(cuenta == gauss(q, d), nombre);
    CHECK_MSG(coinciden, nombre);
}

// ---------------------------------------------------------------- 6. RNG
template <class Field>
void test_rng_intacto(const char* nombre, Field (*crear)()) {
    using P = Poly<Field>;
    set_seed(12345);
    Field K1 = crear();
    P g1 = P::random(K1, 400), f1 = P::random(K1, 900);
    auto r1 = f1.divmod(g1);        // con 400/501 coeficientes: camino de Newton
    const bool sin_tablas = K1.fft_cache_.empty();
    const BigInt sig1 = rand_below(BigInt(1000000007));

    set_seed(12345);
    Field K2 = crear();
    P g2 = P::random(K2, 400), f2 = P::random(K2, 900);
    auto r2 = f2.divmod_school(g2);  // el original
    const BigInt sig2 = rand_below(BigInt(1000000007));

    CHECK_MSG(sin_tablas, nombre);
    CHECK_MSG(sig1 == sig2, nombre);
    CHECK_MSG(r1.first.coeffs().size() == r2.first.coeffs().size(), nombre);
}

}  // namespace

int main() {
    set_seed(20261107);

    // 1. Inverso de series de potencias
    Fp f101(101), fntt(998244353), f127((BigInt(1) << 127) - 1);
    Fp64 w101(101), wntt(998244353);
    Fq f9(Fp(BigInt(3)), std::vector<BigInt>{1, 0, 1});
    Fq64 w9(Fp64(BigInt(3)), std::vector<BigInt>{1, 0, 1});
    test_inverse_series(f101, "inverse_series F_101");
    test_inverse_series(fntt, "inverse_series F_998244353");
    test_inverse_series(w101, "inverse_series Fp64 101");
    test_inverse_series(f9, "inverse_series F_9");
    test_inverse_series(w9, "inverse_series Fq64 F_9");

    // 2. divmod_newton == divmod_school (rejilla de grados)
    const std::vector<int> ns{0, 1, 2, 7, 33, 130, 300};
    const std::vector<int> ds{0, 1, 5, 64, 257, 520};
    int casos = 0;
    casos += test_divmod(w101, ns, ds, "divmod Fp64 101");
    casos += test_divmod(wntt, ns, ds, "divmod Fp64 998244353 (sin tablas)");
    casos += test_divmod(f101, {0, 3, 40, 260}, {0, 9, 270}, "divmod F_101 GMP");
    casos += test_divmod(f127, {0, 5, 70}, {0, 3, 80}, "divmod F_(2^127-1)");
    casos += test_divmod(f9, {0, 2, 50, 140}, {0, 7, 150}, "divmod F_9");
    casos += test_divmod(w9, {0, 2, 50, 140}, {0, 7, 150}, "divmod Fq64 F_9");
    // Con las tablas FFT ya construidas (camino FFT de mul_no_rng).
    {
        using P = Poly<Fp64>;
        (void)(P::random(wntt, 600).mul(P::random(wntt, 600), MulAlgo::FFT));
        CHECK(!wntt.fft_cache_.empty());
        casos += test_divmod(wntt, {130, 300, 700}, {64, 300, 900}, "divmod Fp64 998244353 (con tablas)");
        using G = Poly<Fp>;
        (void)(G::random(fntt, 300).mul(G::random(fntt, 300), MulAlgo::FFT));
        casos += test_divmod(fntt, {40, 300}, {9, 300}, "divmod F_998244353 GMP (con tablas)");
    }
    std::printf("  divmod: %d combinaciones de grados comparadas\n", casos);

    // 3. pow_mod con Reducer frente a la potencia con división escolar
    test_pow_mod(w101, 300, BigInt(101), "pow_mod Fp64 101 grado 300");
    test_pow_mod(wntt, 600, BigInt(998244353), "pow_mod Fp64 NTT grado 600");
    test_pow_mod(f101, 280, BigInt(10007), "pow_mod F_101 GMP grado 280");
    test_pow_mod(w9, 150, BigInt(9), "pow_mod Fq64 F_9 grado 150");
    test_pow_mod(w101, 20, BigInt(12345), "pow_mod Fp64 grado 20 (sin Newton)");

    // 4. Inverso en F_q: Euclides extendido == Fermat
    Fq f256(Fp(BigInt(2)), std::vector<BigInt>{1, 1, 0, 1, 1, 0, 0, 0, 1});
    Fq64 w256(Fp64(BigInt(2)), std::vector<BigInt>{1, 1, 0, 1, 1, 0, 0, 0, 1});
    Fq64 w2401(BigInt(7), 4UL);
    test_inv_todos(f9, "inv F_9 (todos)");
    test_inv_todos(w9, "inv Fq64 F_9 (todos)");
    test_inv_todos(f256, "inv F_256 (todos)");
    test_inv_todos(w256, "inv Fq64 F_256 (todos)");
    test_inv_todos(w2401, "inv Fq64 F_(7^4) (todos)");
    const std::vector<BigInt> mod_3_30{2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0,
                                       2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1};
    Fq f3e30(Fp(BigInt(3)), mod_3_30);
    Fq64 w3e30(Fp64(BigInt(3)), mod_3_30);
    test_inv_muestras(f3e30, 200, "inv F_(3^30)");
    test_inv_muestras(w3e30, 200, "inv Fq64 F_(3^30)");
    const BigInt m61 = (BigInt(1) << 61) - 1;
    Fq fm61e2(m61, 2UL);
    test_inv_muestras(fm61e2, 200, "inv F_((2^61-1)^2)");
    Fq fq1(Fp(BigInt(5)), std::vector<BigInt>{2, 1});  // extensión de grado 1
    test_inv_todos(fq1, "inv F_5 como extensión de grado 1");

    // 5. Irreducibilidad: nuevo == original y recuento de Gauss
    Fp f2(2), f3(3);
    for (int d = 1; d <= 8; ++d) test_irreducibles_gauss(f2, d, "irreducibles F_2");
    for (int d = 1; d <= 5; ++d) test_irreducibles_gauss(f3, d, "irreducibles F_3");
    for (int d = 1; d <= 3; ++d) test_irreducibles_gauss(f9, d, "irreducibles F_9");
    for (int d = 1; d <= 3; ++d) test_irreducibles_gauss(w9, d, "irreducibles Fq64 F_9");
    {
        // Grados altos en F_101 (Newton activo): irreducible aleatorio y productos.
        using P = Poly<Fp64>;
        P u = P::random(w101, 260).monic();
        while (!u.is_irreducible()) u = P::random(w101, 260).monic();
        CHECK(ben_or_original(u));
        P v = P::random(w101, 3).monic();
        while (!v.is_irreducible()) v = P::random(w101, 3).monic();
        CHECK(!(u * v).is_irreducible());
        CHECK(!(u * u).is_irreducible());
    }

    // 6. El camino de Newton no consume el generador ni construye tablas FFT
    test_rng_intacto<Fp64>("RNG Fp64 NTT", [] { return Fp64(998244353); });
    test_rng_intacto<Fp>("RNG Fp NTT", [] { return Fp(998244353); });

    return check::summary("fase4");
}
