// test_fase7_hgcd.cpp — Fase 7, propuesta 2: half-GCD «V3e» en Fp64 (hgcd con
// restos, matrices 2x2 con transformadas NTT reutilizadas y envoltura cíclica).
//
//   1. Envoltura (namespace hgcdw) con p = 29·2^57 + 1 y coeficientes p - 1:
//      la corrección de los coeficientes altos suma más de 16 productos de
//      ~2^124 (la reserva del revisor de la fase 6: un acumulador de 128 bits
//      sin reducir desborda).  Se compara con el producto escolar.
//   2. gcd y gcd_ext (g, x, y) frente a gcd_classic / gcd_ext_classic (sin
//      cambios desde 80b983d: mismos cocientes) en F_101, F_998244353,
//      F_(2^61-1), F_(29·2^57+1) y F_65537, grados 300-8200, con mcd grande,
//      a < b, b = 0, b constante, grados iguales y coeficientes p - 1.
//   3. Generador aleatorio en frío: gcd y gcd_ext sobre un cuerpo recién creado
//      consumen el generador y construyen tablas exactamente como antes.
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <cstdio>
#include <string>
#include <vector>

using namespace alcp;
using P = Poly<Fp64>;
using V = std::vector<std::uint64_t>;

namespace {

// --------------------------------------------------------------- 1. envoltura
// x1·y1 + x2·y2 con la NTT de tamaño 2^(k-1) y la corrección de final_ frente al
// producto escolar (mul_school).
bool envoltura(const Fp64& K, std::size_t l1, std::size_t l2, std::uint64_t v) {
    const std::uint64_t p = K.modulus64();
    V x1(l1, v), y1(l2, v), x2(l2, v), y2(l1, v);
    x1[0] = 1 % p;  // un coeficiente distinto para no tener un caso simétrico
    const std::size_t L = l1 + l2 - 1;
    const unsigned long k = hgcdw::elige_k(L);
    const std::size_t n = std::size_t(1) << k;
    if (n >= L) return false;  // el caso debe envolver
    const auto& t = fft::tables(K, k);
    std::vector<std::uint64_t> out;
    hgcdw::fma2(out, hgcdw::directa(t, x1, n, p), hgcdw::directa(t, y1, n, p), hgcdw::directa(t, x2, n, p),
                hgcdw::directa(t, y2, n, p), p);
    const V r = hgcdw::final_(t, out, x1, y1, x2, y2, p);
    const P ref = P(K, x1).mul_school(P(K, y1)) + P(K, x2).mul_school(P(K, y2));
    return P(K, r) == ref;
}

void test_envoltura() {
    Fp64 w62(BigInt("4179340454199820289")), w998(998244353);
    (void)(P::random(w62, 300) * P::random(w62, 300));
    (void)(P::random(w998, 300) * P::random(w998, 300));
    bool ok = hgcdw::elige_k(1088) == 10 && hgcdw::elige_k(1089) == 11 && hgcdw::elige_k(1025) == 10 &&
              hgcdw::elige_k(40) == 6;
    CHECK_MSG(ok, "elige_k: envoltura hasta 2^(k-1) + 2^(k-1)/16");
    ok = true;
    // L = 1087: los coeficientes 1024..1086 se calculan directamente, con
    // hasta 63 productos de (p - 1)^2 ~ 2^124 cada uno.
    for (const Fp64* K : {&w62, &w998})
        for (auto [l1, l2] : {std::pair<std::size_t, std::size_t>{600, 488}, {1000, 88}, {544, 545}, {70, 1018}})
            ok = ok && envoltura(*K, l1, l2, K->modulus64() - 1);
    CHECK_MSG(ok, "envoltura con coeficientes p - 1 (p = 29·2^57 + 1 y 998244353)");
    // La suma directa sin contar acarreos desbordaría: 64 productos de (p-1)^2.
    const std::uint64_t p = w62.modulus64();
    const unsigned __int128 c = static_cast<unsigned __int128>(p - 1) * (p - 1);
    CHECK_MSG(c * 64 / 64 != c, "(p - 1)^2 · 64 no cabe en 128 bits (la prueba ejercita los acarreos)");
}

// ---------------------------------------------------------- 2. gcd y gcd_ext
bool igual(const P& a, const P& b) {
    const P g = P::gcd(a, b), gc = P::gcd_classic(a, b);
    const auto e = P::gcd_ext(a, b), ec = P::gcd_ext_classic(a, b);
    return g.coeffs() == gc.coeffs() && e.g.coeffs() == ec.g.coeffs() && e.x.coeffs() == ec.x.coeffs() &&
           e.y.coeffs() == ec.y.coeffs();
}

P todos(const Fp64& K, int n, std::uint64_t v) { return P(K, V(static_cast<std::size_t>(n) + 1, v)); }

void test_gcd(const Fp64& K, const char* nombre, const std::vector<int>& grados) {
    const std::uint64_t pm1 = K.modulus64() - 1;
    bool ok = true;
    for (int n : grados) {
        ok = ok && igual(P::random(K, n), P::random(K, n - 1));
        ok = ok && igual(P::random(K, n), P::random(K, n - 300));
        ok = ok && igual(P::random(K, n / 2), P::random(K, n));  // a < b
        const P g = P::random(K, n / 3);
        ok = ok && igual(g * P::random(K, n - n / 3), g * P::random(K, n - n / 3 - 5));  // mcd grande
        const P h = P::random(K, n - 40);
        ok = ok && igual(h * P::random(K, 40), h * P::random(K, 39));  // mcd casi entero
        ok = ok && igual(P::random(K, n), P::zero(K)) && igual(P::random(K, n), P::constant(K, 7 % K.modulus64()));
        ok = ok && igual(P::zero(K), P::random(K, n));
        ok = ok && igual(P::random(K, n), P::random(K, n));  // grados iguales
        ok = ok && igual(todos(K, n, pm1), P::random(K, n - 1));
        ok = ok && igual(todos(K, n, pm1) + P::monomial(K, 3, 5 % K.modulus64()),
                         todos(K, n - 1, pm1) - P::monomial(K, static_cast<std::size_t>(n / 2), 1));
    }
    CHECK_MSG(ok, (std::string("gcd y gcd_ext = clásicos, ") + nombre).c_str());
}

// ---------------------------------------------------------------- 3. RNG
void test_rng() {
    for (const char* ps : {"998244353", "4179340454199820289", "65537", "101"}) {
        // Con gcd / gcd_ext.
        set_seed(4242);
        Fp64 K1{BigInt(ps)};
        const P a1 = P::random(K1, 3000), b1 = P::random(K1, 2999);
        const P g1 = P::gcd(a1, b1);
        const bool vacia1 = K1.fft_cache_.empty();
        const BigInt s1 = rand_below(BigInt(1000000007));
        const auto e1 = P::gcd_ext(a1, b1);
        const bool vacia1e = K1.fft_cache_.empty();
        const BigInt t1 = rand_below(BigInt(1000000007));
        // Con los clásicos (el comportamiento de siempre).
        set_seed(4242);
        Fp64 K2{BigInt(ps)};
        const P a2 = P::random(K2, 3000), b2 = P::random(K2, 2999);
        const P g2 = P::gcd_classic(a2, b2);
        const bool vacia2 = K2.fft_cache_.empty();
        const BigInt s2 = rand_below(BigInt(1000000007));
        const auto e2 = P::gcd_ext_classic(a2, b2);
        const bool vacia2e = K2.fft_cache_.empty();
        const BigInt t2 = rand_below(BigInt(1000000007));
        CHECK_MSG(s1 == s2 && t1 == t2 && vacia1 && vacia1 == vacia2 && vacia1e == vacia2e &&
                      g1.coeffs() == g2.coeffs() && e1.x.coeffs() == e2.x.coeffs() && e1.y.coeffs() == e2.y.coeffs(),
                  (std::string("RNG gcd / gcd_ext en frío p = ") + ps).c_str());
    }
}

}  // namespace

int main() {
    set_seed(20261108);
    test_envoltura();
    Fp64 w101(101), w998(998244353), wm61((BigInt(1) << 61) - 1), w62(BigInt("4179340454199820289")), w65537(65537);
    for (const Fp64* K : {&w101, &w998, &wm61, &w62, &w65537}) (void)(P::random(*K, 300) * P::random(*K, 300));
    const std::vector<int> grados{300, 700, 1000, 1500, 2100, 4200};
    test_gcd(w101, "F_101", grados);
    test_gcd(w998, "F_998244353", grados);
    test_gcd(wm61, "F_(2^61-1)", grados);
    test_gcd(w62, "F_(29·2^57+1)", grados);
    test_gcd(w65537, "F_65537", {300, 1000, 2100});
    {
        // Grado 8200 (tres niveles de recursión por encima del caso base).
        bool ok = true;
        for (const Fp64* K : {&w998, &w62, &w101}) {
            ok = ok && igual(P::random(*K, 8200), P::random(*K, 8199));
            const P g = P::random(*K, 3000);
            ok = ok && igual(g * P::random(*K, 5200), g * P::random(*K, 5150));
        }
        CHECK_MSG(ok, "gcd y gcd_ext = clásicos, grado 8200");
    }
    test_rng();
    return check::summary("test_fase7_hgcd");
}
