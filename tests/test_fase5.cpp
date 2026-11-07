// test_fase5.cpp — Fase 5: mismas salidas que antes con las tres mejoras.
//
//   1. NTT en palabras con mariposas perezosas: producto FFT frente al escolar
//      en Fp64 y en Fp (GMP) con p = 998244353 y con un primo NTT cercano a
//      2^62 (el caso límite de 4p < 2^64), con coeficientes p - 1 (peor caso)
//      y transformada directa + inversa = identidad.
//   2. Euclides en sitio (Fp64): gcd_classic y gcd_ext_classic frente a una
//      copia del algoritmo original sobre Poly (mismos mcd y cofactores de
//      Bézout), con ceros, grados iguales, a < b, mcd grande y cocientes
//      grandes (camino de Newton y producto rápido de cofactores).
//      También half-GCD en Fp64 por encima de los umbrales recalibrados.
//   3. Generador aleatorio: gcd_ext_classic y distinct_degree sobre un cuerpo
//      recién creado consumen el generador igual que el original.
//   4. DDF de Kaltofen-Shoup: distinct_degree frente a distinct_degree_classic
//      en varios cuerpos, con entradas aleatorias y con muchos factores de
//      grados pequeños (separación fina y cambio de módulo); factor() sigue
//      reconstruyendo la entrada con factores irreducibles.
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <cstdio>
#include <string>
#include <vector>

using namespace alcp;

namespace {

// ---------------------------------------------------------------- 1. NTT
template <class Field>
void test_ntt(const Field& K, const char* nombre) {
    using P = Poly<Field>;
    bool ok = true;
    for (int n : {1, 2, 31, 64, 65, 300, 1000, 2047}) {
        P a = P::random(K, n), b = P::random(K, n / 2 + 1);
        ok = ok && a.mul(b, MulAlgo::FFT) == a.mul_school(b);
        ok = ok && a.mul(a, MulAlgo::FFT) == a.mul_school(a);
    }
    // Todos los coeficientes p - 1: los valores intermedios llegan al máximo.
    P M(K, std::vector<typename Field::Elem>(1500, K.neg(K.one())));
    ok = ok && M.mul(M, MulAlgo::FFT) == M.mul_school(M);
    // Directa + inversa = identidad.
    for (std::size_t n : {std::size_t(2), std::size_t(256), std::size_t(4096)}) {
        std::vector<typename Field::Elem> v;
        for (std::size_t i = 0; i < n; ++i) v.push_back(i % 3 ? K.random() : K.neg(K.one()));
        auto w = v;
        fft::transform(K, w, false);
        fft::transform(K, w, true);
        ok = ok && w == v;
    }
    CHECK_MSG(ok, nombre);
}

// ---------------------------------------------------------------- 2. Euclides
// Copias de los algoritmos originales, con operaciones de Poly.
template <class Field>
Poly<Field> gcd_ref(Poly<Field> a, Poly<Field> b) {
    while (!b.is_zero()) {
        Poly<Field> r = a % b;
        a = std::move(b);
        b = std::move(r);
    }
    return a.is_zero() ? a : a.monic();
}

template <class Field>
typename Poly<Field>::ExtGcd gcd_ext_ref(const Poly<Field>& a, const Poly<Field>& b) {
    using P = Poly<Field>;
    const Field& K = a.field();
    P x0 = P::one(K), x1 = P::zero(K), y0 = P::zero(K), y1 = P::one(K), r0 = a, r1 = b;
    while (!r1.is_zero()) {
        auto [q, r] = r0.divmod(r1);
        r0 = std::move(r1); r1 = std::move(r);
        P tx = x0 - q * x1; x0 = std::move(x1); x1 = std::move(tx);
        P ty = y0 - q * y1; y0 = std::move(y1); y1 = std::move(ty);
    }
    if (r0.is_zero()) return {P::zero(K), P::zero(K), P::zero(K)};
    const auto inv = K.inv(r0.lead());
    return {r0.scalar_mul(inv), x0.scalar_mul(inv), y0.scalar_mul(inv)};
}

bool igual(const Poly<Fp64>& a, const Poly<Fp64>& b) {
    using P = Poly<Fp64>;
    const auto e1 = P::gcd_ext_classic(a, b), e2 = gcd_ext_ref(a, b);
    return P::gcd_classic(a, b) == gcd_ref(a, b) && P::gcd(a, b) == gcd_ref(a, b) && e1.g == e2.g &&
           e1.x == e2.x && e1.y == e2.y && P::gcd_ext(a, b).x == e2.x;
}

void test_euclides(const Fp64& K, const char* nombre) {
    using P = Poly<Fp64>;
    bool ok = true;
    const std::vector<int> grados{0, 1, 2, 7, 40, 129, 300};
    for (int da : grados)
        for (int db : grados) ok = ok && igual(P::random(K, da), P::random(K, db));
    for (int dg : {1, 30, 150}) {  // mcd grande
        P g = P::random(K, dg);
        ok = ok && igual(g * P::random(K, 200), g * P::random(K, 170));
    }
    // Cocientes grandes: división de Newton y producto rápido de cofactores.
    ok = ok && igual(P::random(K, 1500), P::random(K, 200)) && igual(P::random(K, 900), P::random(K, 3));
    P a = P::random(K, 120);
    ok = ok && igual(a, P::zero(K)) && igual(P::zero(K), a) && igual(P::zero(K), P::zero(K)) &&
         igual(a, P::one(K)) && igual(a, a);
    CHECK_MSG(ok, nombre);
}

// ---------------------------------------------------------------- 3. RNG
void test_rng() {
    using P = Poly<Fp64>;
    // gcd_ext_classic (con cofactores de >= 64 coeficientes) frente al original.
    set_seed(4242);
    Fp64 K1(998244353);
    P a1 = P::random(K1, 300), b1 = P::random(K1, 290);
    auto e1 = P::gcd_ext_classic(a1, b1);
    const BigInt s1 = rand_below(BigInt(1000000007));
    set_seed(4242);
    Fp64 K2(998244353);
    P a2 = P::random(K2, 300), b2 = P::random(K2, 290);
    auto e2 = gcd_ext_ref(a2, b2);
    const BigInt s2 = rand_below(BigInt(1000000007));
    CHECK(s1 == s2);
    CHECK(e1.x.coeffs() == e2.x.coeffs() && e1.y.coeffs() == e2.y.coeffs());
    CHECK(K1.fft_cache_.empty() == K2.fft_cache_.empty());

    // distinct_degree (KS) frente a distinct_degree_classic, sobre cuerpos nuevos.
    for (const char* ps : {"998244353", "101", "7340033", "193", "257"}) {
        set_seed(99);
        Fp64 L1{BigInt(ps)};
        P f1 = P::random(L1, 120).monic().radical();
        auto d1 = f1.distinct_degree();
        const BigInt t1 = rand_below(BigInt(1000000007));
        set_seed(99);
        Fp64 L2{BigInt(ps)};
        P f2 = P::random(L2, 120).monic().radical();
        auto d2 = f2.distinct_degree_classic();
        const BigInt t2 = rand_below(BigInt(1000000007));
        CHECK_MSG(t1 == t2 && d1.size() == d2.size() && L1.fft_cache_.empty() == L2.fft_cache_.empty(),
                  (std::string("RNG distinct_degree p = ") + ps).c_str());
    }
}

// ---------------------------------------------------------------- 4. DDF
template <class Field>
bool igual_ddf(const Poly<Field>& f) {
    const auto a = f.distinct_degree(), b = f.distinct_degree_classic();
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i].first != b[i].first || !(a[i].second == b[i].second)) return false;
    return true;
}

template <class Field>
void test_ddf(const Field& K, std::vector<int> grados, int veces, const char* nombre) {
    using P = Poly<Field>;
    bool ok = true;
    for (int n : grados)
        for (int t = 0; t < veces; ++t) ok = ok && igual_ddf(P::random(K, n).monic().radical());
    // Muchos factores pequeños: producto de polinomios aleatorios de grado
    // 1..12 (su parte libre de cuadrados), de grado total >= 100.
    for (int t = 0; t < 3; ++t) {
        P f = P::one(K);
        int d = 1;
        while (f.degree() < 100) {
            f = f * P::random(K, d).monic();
            d = d % 12 + 1;
        }
        P r = f.radical();
        ok = ok && igual_ddf(r);
        // factor() reconstruye la entrada con factores irreducibles.
        auto F = r.factor();
        P prod = P::one(K);
        for (auto& [u, e] : F.factors) {
            ok = ok && u.is_irreducible();
            for (int i = 0; i < e; ++i) prod = prod * u;
        }
        ok = ok && prod.scalar_mul(F.leading) == r;
    }
    CHECK_MSG(ok, nombre);
}

}  // namespace

int main() {
    set_seed(20261107);
    // 4179340454199820289 = 29·2^57 + 1 (primo NTT < 2^62)
    Fp64 wntt(998244353), wgrande(BigInt("4179340454199820289"));
    Fp gntt(BigInt(998244353)), ggrande(BigInt("4179340454199820289"));
    test_ntt(wntt, "NTT Fp64 p = 998244353");
    test_ntt(wgrande, "NTT Fp64 p = 29·2^57 + 1");
    test_ntt(gntt, "NTT Fp GMP p = 998244353");
    test_ntt(ggrande, "NTT Fp GMP p = 29·2^57 + 1");

    Fp64 w2(2), w3(3), w101(101), wm61((BigInt(1) << 61) - 1), w62((BigInt(1) << 62) - 57);
    test_euclides(w2, "Euclides F_2");
    test_euclides(w3, "Euclides F_3");
    test_euclides(w101, "Euclides F_101");
    test_euclides(wntt, "Euclides F_998244353");
    test_euclides(wm61, "Euclides F_(2^61-1)");
    test_euclides(w62, "Euclides F_(2^62-57)");

    // Half-GCD en Fp64 por encima de los umbrales de la fase 5 (caso base en
    // sitio, hgcd_base_w): mismos mcd y cofactores que el original.
    {
        using P = Poly<Fp64>;
        bool ok = true;
        for (const Fp64* K : {&w3, &w101}) {
            ok = ok && igual(P::random(*K, 2100), P::random(*K, 2099));
            P g = P::random(*K, 700);
            ok = ok && igual(g * P::random(*K, 1500), g * P::random(*K, 1450));
        }
        (void)(P::random(wntt, 600) * P::random(wntt, 600));  // con tablas FFT
        ok = ok && igual(P::random(wntt, 2100), P::random(wntt, 2090));
        CHECK_MSG(ok, "half-GCD Fp64 (F_3, F_101, F_998244353)");
    }

    test_rng();

    Fq64 w9(Fp64(BigInt(3)), std::vector<BigInt>{1, 0, 1});
    Fp g127((BigInt(1) << 127) - 1);
    test_ddf(w2, {40, 100, 300, 400}, 3, "DDF F_2");
    test_ddf(w3, {40, 100, 300}, 3, "DDF F_3");
    test_ddf(w101, {39, 40, 41, 64, 150, 250}, 4, "DDF F_101");
    test_ddf(wntt, {40, 100, 300}, 3, "DDF F_998244353");
    test_ddf(wm61, {40, 80}, 3, "DDF F_(2^61-1)");
    test_ddf(w9, {40, 70}, 3, "DDF F_9");
    test_ddf(g127, {40, 50}, 2, "DDF F_(2^127-1)");
    return check::summary("fase5");
}
