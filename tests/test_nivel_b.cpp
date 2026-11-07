// test_nivel_b.cpp — Fase 4, nivel B: algoritmos nuevos con el mismo resultado que los originales.
//
//   1. Producto de Kronecker sobre Z (mul_kronecker / MulAlgo::Kronecker) frente
//      al escolar: F_p pequeño, con NTT, p cercano a 2^62 (huecos de > 128 bits),
//      p con GMP (2^127 - 1, 2^521 - 1), F_q con los dos backends; operandos de
//      un coeficiente, de longitudes distintas, cuadrados.
//   2. Half-GCD: gcd == gcd_classic y gcd_ext == gcd_ext_classic (mcd Y
//      coeficientes de Bézout idénticos) en F_2 (muchas caídas de grado), F_3,
//      F_101, F_998244353 (con y sin tablas FFT), F_9, F_(2^127-1); con mcd
//      grande, grados iguales, a < b, ceros.
//   3. Generador aleatorio: gcd y gcd_ext sobre un cuerpo recién creado no
//      adelantan el cálculo de la raíz de la FFT.
//   4. Composición modular (Brent-Kung) y Frobenius: compose_mod frente a
//      compose + %, y frobenius_powers frente a pow_mod iterado.
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <cstdio>
#include <vector>

using namespace alcp;

namespace {

// ---------------------------------------------------------------- 1. Kronecker
template <class Field>
void test_kronecker(const Field& K, const std::vector<int>& grados, const char* nombre) {
    using P = Poly<Field>;
    bool ok = true;
    for (int da : grados)
        for (int db : grados) {
            P a = P::random(K, da), b = P::random(K, db);
            ok = ok && a.mul(b, MulAlgo::Kronecker) == a.mul_school(b);
        }
    for (int d : grados) {
        P a = P::random(K, d);
        ok = ok && a.mul(a, MulAlgo::Kronecker) == a.mul_school(a);  // cuadrado
    }
    // Coeficientes máximos (p - 1) en todas partes: el peor caso de los huecos.
    P M(K, std::vector<typename Field::Elem>(200, K.neg(K.one())));
    ok = ok && M.mul(M, MulAlgo::Kronecker) == M.mul_school(M);
    ok = ok && P::zero(K).mul(M, MulAlgo::Kronecker).is_zero();
    CHECK_MSG(ok, nombre);
}

// ---------------------------------------------------------------- 2. half-GCD
template <class Field>
bool igual_gcd(const Poly<Field>& a, const Poly<Field>& b) {
    using P = Poly<Field>;
    const auto e1 = P::gcd_ext(a, b), e2 = P::gcd_ext_classic(a, b);
    return P::gcd(a, b) == P::gcd_classic(a, b) && e1.g == e2.g && e1.x == e2.x && e1.y == e2.y &&
           e1.g == e1.x * a + e1.y * b;
}

template <class Field>
void test_hgcd(const Field& K, const std::vector<int>& grados, const char* nombre) {
    using P = Poly<Field>;
    bool ok = true;
    int casos = 0;
    for (int da : grados)
        for (int db : grados) {
            P a = P::random(K, da), b = P::random(K, db);
            ok = ok && igual_gcd(a, b);
            ++casos;
        }
    // mcd grande: a = g·u, b = g·v
    for (int dg : {1, 60, 200}) {
        P g = P::random(K, dg);
        P a = g * P::random(K, 300), b = g * P::random(K, 250);
        ok = ok && igual_gcd(a, b) && P::gcd(a, b).degree() >= dg;
        ++casos;
    }
    // Ceros y constantes
    P a = P::random(K, 300);
    ok = ok && igual_gcd(a, P::zero(K)) && igual_gcd(P::zero(K), a) && igual_gcd(P::zero(K), P::zero(K)) &&
         igual_gcd(a, P::one(K));
    CHECK_MSG(ok, nombre);
    std::printf("  %s: %d pares\n", nombre, casos);
}

// ---------------------------------------------------------------- 3. RNG
void test_rng_gcd() {
    using P = Poly<Fp64>;
    set_seed(777);
    Fp64 K1(998244353);
    P a1 = P::random(K1, 500), b1 = P::random(K1, 480);
    auto e1 = P::gcd_ext(a1, b1);
    auto g1 = P::gcd(a1, b1);
    const BigInt s1 = rand_below(BigInt(1000000007));

    set_seed(777);
    Fp64 K2(998244353);
    P a2 = P::random(K2, 500), b2 = P::random(K2, 480);
    auto e2 = P::gcd_ext_classic(a2, b2);
    auto g2 = P::gcd_classic(a2, b2);
    const BigInt s2 = rand_below(BigInt(1000000007));
    CHECK(s1 == s2);
    CHECK(e1.x.coeffs() == e2.x.coeffs() && g1.coeffs() == g2.coeffs());
    // gcd por sí solo no debe construir tablas FFT.
    Fp64 K3(998244353);
    P a3 = P::random(K3, 500), b3 = P::random(K3, 480);
    (void)P::gcd(a3, b3);
    CHECK(K3.fft_cache_.empty());
}

// ---------------------------------------------------------------- 4. composición y Frobenius
template <class Field>
void test_composer(const Field& K, int n, std::size_t k, const char* nombre) {
    using P = Poly<Field>;
    P f = P::random(K, n);
    typename P::Reducer red(f);
    P h = P::random(K, n - 1);
    typename P::Composer comp(h, red, k);
    bool ok = true;
    for (int d : {0, 1, n / 3, n - 1}) {
        P g = P::random(K, d);
        ok = ok && comp.compose(g) == g.compose(h) % f;
    }
    ok = ok && comp.compose(P::zero(K)).is_zero();
    CHECK_MSG(ok, nombre);
}

// distinct_degree original: H_d = H_{d-1}^q con pow_mod.
template <class Field>
std::vector<std::pair<int, Poly<Field>>> ddf_original(const Poly<Field>& f) {
    using P = Poly<Field>;
    const BigInt q = f.field().cardinality();
    std::vector<std::pair<int, P>> out;
    P g = f.monic(), x = P::variable(f.field()), H = x % g;
    int d = 1;
    while (g.degree() >= 2 * d) {
        H = P::pow_mod(H, q, g);
        P hd = P::gcd(H - x, g);
        if (hd.degree() > 0) {
            out.emplace_back(d, hd.monic());
            g = g / hd;
            H = H % g;
        }
        ++d;
    }
    if (g.degree() > 0) out.emplace_back(g.degree(), g.monic());
    return out;
}

template <class Field>
void test_frobenius(const Field& K, int n, int veces, const char* nombre) {
    using P = Poly<Field>;
    bool ok = true;
    int irreducibles = 0;
    for (int t = 0; t < veces; ++t) {
        P f = P::random(K, n).monic();
        const bool a = f.is_irreducible(FrobeniusAlgo::Potencia);
        const bool b = f.is_irreducible(FrobeniusAlgo::Composicion);
        const bool c = f.is_irreducible();
        ok = ok && a == b && a == c;
        irreducibles += a;
        // distinct_degree (Auto) frente al original, sobre la parte libre de cuadrados
        P r = f.radical();
        auto d1 = r.distinct_degree();
        auto d2 = ddf_original(r);
        ok = ok && d1.size() == d2.size();
        for (std::size_t i = 0; ok && i < d1.size(); ++i) ok = d1[i].first == d2[i].first && d1[i].second == d2[i].second;
    }
    // Un irreducible seguro: se busca y se comprueba con los dos métodos.
    P u = P::random(K, n).monic();
    while (!u.is_irreducible(FrobeniusAlgo::Potencia)) u = P::random(K, n).monic();
    ok = ok && u.is_irreducible(FrobeniusAlgo::Composicion) && u.is_irreducible();
    CHECK_MSG(ok, nombre);
    std::printf("  %s: %d de %d irreducibles\n", nombre, irreducibles, veces);
}

}  // namespace

int main() {
    set_seed(20261108);
    const std::vector<int> gk{0, 1, 2, 31, 32, 33, 100, 257, 700};
    Fp64 w2(2), w101(101), wntt(998244353), wbig((BigInt(1) << 62) - 57);  // 2^62 - 57 es primo
    Fp g7(7), g127((BigInt(1) << 127) - 1), g521((BigInt(1) << 521) - 1);
    Fq64 w9(Fp64(BigInt(3)), std::vector<BigInt>{1, 0, 1});
    Fq64 w256(Fp64(BigInt(2)), std::vector<BigInt>{1, 1, 0, 1, 1, 0, 0, 0, 1});
    Fq f9(Fp(BigInt(3)), std::vector<BigInt>{1, 0, 1});
    Fq64 wm61e2(Fp64((BigInt(1) << 61) - 1), std::vector<BigInt>{1, 0, 1});

    // 1. Kronecker
    test_kronecker(w2, gk, "Kronecker Fp64 p = 2");
    test_kronecker(w101, gk, "Kronecker Fp64 p = 101");
    test_kronecker(wntt, gk, "Kronecker Fp64 p = 998244353");
    test_kronecker(wbig, gk, "Kronecker Fp64 p = 2^62 - 57");
    test_kronecker(g7, gk, "Kronecker Fp GMP p = 7");
    test_kronecker(g127, {0, 1, 40, 150}, "Kronecker Fp GMP p = 2^127 - 1");
    test_kronecker(g521, {0, 1, 30, 70}, "Kronecker Fp GMP p = 2^521 - 1");
    test_kronecker(w9, {0, 1, 33, 120}, "Kronecker Fq64 F_9");
    test_kronecker(w256, {0, 1, 33, 120}, "Kronecker Fq64 F_256");
    test_kronecker(f9, {0, 1, 33, 120}, "Kronecker Fq GMP F_9");
    test_kronecker(wm61e2, {0, 1, 33, 90}, "Kronecker Fq64 F_((2^61-1)^2)");

    // 2. half-GCD
    const std::vector<int> gh{0, 1, 5, 127, 128, 129, 400, 1100};
    test_hgcd(w2, gh, "hgcd F_2");
    Fp64 w3(3);
    test_hgcd(w3, gh, "hgcd F_3");
    test_hgcd(w101, gh, "hgcd F_101");
    (void)(Poly<Fp64>::random(wntt, 600).mul(Poly<Fp64>::random(wntt, 600), MulAlgo::FFT));  // con tablas
    test_hgcd(wntt, gh, "hgcd F_998244353");
    test_hgcd(w9, {0, 1, 130, 400}, "hgcd Fq64 F_9");
    test_hgcd(g127, {0, 1, 130, 300}, "hgcd F_(2^127-1)");

    // 3. Generador aleatorio
    test_rng_gcd();

    // 4. Composición modular y Frobenius
    test_composer(w101, 60, 60, "Composer F_101 k = n");
    test_composer(w101, 60, 7, "Composer F_101 k = 7 (pasos gigantes)");
    test_composer(wntt, 200, 15, "Composer F_998244353");
    test_composer(w9, 40, 6, "Composer F_9");
    test_composer(g127, 30, 5, "Composer F_(2^127-1)");
    test_frobenius(w101, 40, 20, "Frobenius F_101 grado 40");
    test_frobenius(wntt, 80, 5, "Frobenius F_998244353 grado 80");
    test_frobenius(w2, 60, 30, "Frobenius F_2 grado 60");
    test_frobenius(w256, 20, 15, "Frobenius F_256 grado 20");
    test_frobenius(g127, 12, 8, "Frobenius F_(2^127-1) grado 12");
    {
        const std::vector<BigInt> mod_3_30{2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0,
                                           2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1};
        Fq64 w3e30(Fp64(BigInt(3)), mod_3_30);
        test_frobenius(w3e30, 6, 6, "Frobenius F_(3^30) grado 6");
    }

    return check::summary("nivel_b");
}
