// test_fase6_newton.cpp — Fase 6, tema «newton»: mismas salidas que antes con
// el producto medio en la división de Newton y las transformadas precalculadas
// del Reducer y del Composer (solo Fp64 con tablas FFT ya construidas).
//
//   1. inverse_series (producto medio) frente a una copia del algoritmo
//      anterior (misma recursión con productos escolares) y frente al mismo
//      cálculo en un cuerpo «frío» (sin tablas: ruta anterior), en tamaños
//      alrededor de potencias de 2 y del umbral de 64.
//   2. divmod_newton / divmod (resto módulo x^N - 1) frente a divmod_school:
//      d1 < n, d1 > n, n potencia de 2 (plegado de g) y 2^j + 1, resto nulo,
//      coeficientes todos p - 1 y polinomios dispersos.
//   3. Reducer::reduce frente a divmod_school con deg a en {n - 1, n, n + 1,
//      2n - 2, 2n - 1} y n = 2^j, 2^j ± 1, más el cero; también un Reducer que
//      empieza en frío y sigue en caliente.
//   4. pow_mod y Composer::compose frente a copias de los algoritmos anteriores.
//   5. Cuerpos: F_998244353, F_(29·2^57 + 1) (62 bits), F_65537, F_7681
//      (capacidad FFT 9: los tamaños mayores vuelven a la ruta anterior),
//      F_101 y F_(2^61 - 1) (sin NTT útil).
//   6. Generador aleatorio en frío: divmod, inverse_series, Reducer y Composer
//      sobre un cuerpo recién creado no consumen el generador ni construyen
//      tablas; pow_mod lo consume igual que el original.
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <cstdio>
#include <string>
#include <vector>

using namespace alcp;
using P = Poly<Fp64>;
using V = std::vector<std::uint64_t>;

namespace {

// ------------------------------------------------- copias de la ruta anterior
// inverse_series anterior (fase 4) con productos escolares: x_n = x_k·(2 - a·x_k).
V inv_ref(const Fp64& K, const V& a, std::size_t n) {
    if (n == 1) return {K.inv(a[0])};
    const std::size_t k = (n + 1) / 2;
    const V xk = inv_ref(K, a, k);
    const P X(K, xk), A(K, V(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(std::min(n, a.size()))));
    V t = A.mul_school(X).coeffs();
    t.resize(n, 0);
    const P R = P::constant(K, 2) - P(K, t);
    V r = R.mul_school(X).coeffs();
    r.resize(n, 0);
    return r;
}

// pow_mod anterior con la reducción escolar.
P pow_mod_ref(P a, BigInt k, const P& m) {
    P result = P::one(m.field()).divmod_school(m).second;
    a = a.divmod_school(m).second;
    while (k > 0) {
        if (mpz_odd_p(k.get_mpz_t())) result = (result * a).divmod_school(m).second;
        k >>= 1;
        if (k > 0) a = (a * a).divmod_school(m).second;
    }
    return result;
}

// g(h) mod m por Horner (la composición es única).  Los productos usan la
// regla Auto: puede construir tablas FFT (solo se llama cuando ya da igual).
P compose_ref(const P& g, const P& h, const P& m) {
    const P hr = h.divmod_school(m).second;
    P r = P::zero(m.field());
    for (int i = g.degree(); i >= 0; --i)
        r = (r * hr + P::constant(m.field(), g.coeffs()[static_cast<std::size_t>(i)])).divmod_school(m).second;
    return r;
}

// Polinomio de grado d con todos los coeficientes p - 1.
P todo_menos_uno(const Fp64& K, int d) { return P(K, V(static_cast<std::size_t>(d) + 1, K.modulus64() - 1)); }

// Polinomio disperso de grado d (unos pocos coeficientes no nulos).
P disperso(const Fp64& K, int d) {
    V v(static_cast<std::size_t>(d) + 1, 0);
    v[0] = 1;
    v[static_cast<std::size_t>(d) / 3] = K.random();
    v[static_cast<std::size_t>(d)] = K.modulus64() - 1;
    return P(K, v);
}

// Calienta el cuerpo: un producto con la regla Auto construye las tablas FFT.
void calienta(const Fp64& K) { (void)(P::random(K, 300) * P::random(K, 300)); }

// ---------------------------------------------------------------- 1. inversión
void test_inversion(const Fp64& K, const char* nombre) {
    Fp64 Kf(K.characteristic());  // cuerpo frío: ruta anterior
    bool ok = true;
    for (std::size_t n : {1u, 2u, 31u, 63u, 64u, 65u, 127u, 128u, 129u, 255u, 256u, 257u, 300u, 513u, 1000u, 1025u,
                          2048u, 2049u}) {
        V a = P::random(K, static_cast<int>(n) + 5).coeffs();
        if (a[0] == 0) a[0] = 1;
        const V r = P::inverse_series(K, a, n);
        ok = ok && r == inv_ref(K, a, n);
        ok = ok && r == P::inverse_series(Kf, a, n);
        // a más corto que n (a mod x^n tiene menos coeficientes).
        const V corto(a.begin(), a.begin() + 3);
        ok = ok && P::inverse_series(K, corto, n) == P::inverse_series(Kf, corto, n);
    }
    ok = ok && Kf.fft_cache_.empty();
    CHECK_MSG(ok, (std::string("inverse_series ") + nombre).c_str());
}

// ---------------------------------------------------------------- 2. división
bool divide_igual(const P& f, const P& g) {
    const auto [q1, r1] = f.divmod_newton(g);
    const auto [q2, r2] = f.divmod_school(g);
    const auto [q3, r3] = f.divmod(g);
    return q1 == q2 && r1 == r2 && q3 == q2 && r3 == r2;
}

void test_division(const Fp64& K, const char* nombre) {
    bool ok = true;
    // (deg f, deg g): d1 = n, d1 < n, d1 > n, umbral 64, potencias de 2 y 2^j + 1.
    const std::pair<int, int> casos[] = {{127, 64},   {128, 64},   {200, 64},   {190, 127},  {256, 128},
                                         {257, 129},  {600, 300},  {1000, 900}, {1100, 100}, {2048, 1024},
                                         {2050, 1025}, {4000, 2000}, {2200, 1100}, {1500, 1024}, {3000, 1000}};
    for (auto [df, dg] : casos) {
        ok = ok && divide_igual(P::random(K, df), P::random(K, dg));
        ok = ok && divide_igual(todo_menos_uno(K, df), todo_menos_uno(K, dg));
        ok = ok && divide_igual(disperso(K, df), disperso(K, dg));
        // Resto nulo: f = g·q.
        const P g = P::random(K, dg);
        const P f = g * P::random(K, df - dg);
        ok = ok && divide_igual(f, g) && (f % g).is_zero();
        // Resto de grado bajo: f = g·q + c.
        ok = ok && divide_igual(f + P::constant(K, 5), g);
    }
    CHECK_MSG(ok, (std::string("divmod ") + nombre).c_str());
}

// ---------------------------------------------------------------- 3. Reducer
void test_reducer(const Fp64& K, const char* nombre) {
    bool ok = true;
    for (int n : {64, 65, 127, 128, 129, 255, 256, 257, 300, 511, 512, 513, 1024, 1025}) {
        const P m = P::random(K, n);
        const P::Reducer red(m);
        ok = ok && red.reduce(P::zero(K)).is_zero();
        for (int da : {n - 1, n, n + 1, n + 63, n + 64, n + 128, n + 200, 2 * n - 2, 2 * n - 1}) {
            const P a = P::random(K, da);
            ok = ok && red.reduce(a) == a.divmod_school(m).second;
        }
        const P sq = P::random(K, n - 1) * P::random(K, n - 1);
        ok = ok && red.reduce(sq) == sq.divmod_school(m).second;
        const P mm = todo_menos_uno(K, 2 * n - 2);
        ok = ok && red.reduce(mm) == mm.divmod_school(m).second;
        // Copia del Reducer (comparte las transformadas) y múltiplo exacto.
        const P::Reducer copia = red;
        const P mult = m * P::random(K, n - 1);
        ok = ok && copia.reduce(mult).is_zero() && copia.reduce(sq) == red.reduce(sq);
        // Módulo con coeficientes p - 1 y disperso.
        for (const P& m2 : {todo_menos_uno(K, n), disperso(K, n)}) {
            const P::Reducer r2(m2);
            const P b = P::random(K, 2 * n - 2);
            ok = ok && r2.reduce(b) == b.divmod_school(m2).second;
        }
    }
    CHECK_MSG(ok, (std::string("Reducer ") + nombre).c_str());
}

// Reducer creado con el cuerpo frío (inverso por la ruta anterior) que sigue
// reduciendo cuando el cuerpo ya tiene tablas (ruta nueva con ese inverso).
void test_reducer_frio_caliente() {
    Fp64 K(998244353);
    const P m = P::random(K, 700);
    const P::Reducer red(m);
    const P a = P::random(K, 1398), b = P::random(K, 1390);
    bool ok = red.reduce(a) == a.divmod_school(m).second && K.fft_cache_.empty();
    calienta(K);
    ok = ok && !K.fft_cache_.empty();
    ok = ok && red.reduce(b) == b.divmod_school(m).second && red.reduce(a) == a.divmod_school(m).second;
    CHECK_MSG(ok, "Reducer frío -> caliente");
}

// ------------------------------------------------------ 4. pow_mod y Composer
void test_pow_compose(const Fp64& K, const char* nombre) {
    bool ok = true;
    for (int n : {64, 129, 300, 513}) {
        const P m = P::random(K, n);
        const P::Reducer red(m);
        const P a = P::random(K, n + 10);
        ok = ok && P::pow_mod(a, K.cardinality(), red) == pow_mod_ref(a, K.cardinality(), m);
        ok = ok && P::pow_mod(a, BigInt(1000003), m) == pow_mod_ref(a, BigInt(1000003), m);
        const P h = P::random(K, n - 1), h2 = P::random(K, 2 * n);  // h2: grado >= 2n - 1
        const P g0 = P::random(K, n - 1);
        const P ref0 = compose_ref(g0, h, m);
        for (std::size_t k : {std::size_t(1), std::size_t(7), std::size_t(18), static_cast<std::size_t>(n)}) {
            const P::Composer comp(h, red, k);
            ok = ok && comp.compose(g0) == ref0;
            ok = ok && comp.compose(P::zero(K)).is_zero();
            ok = ok && comp.compose(P::constant(K, 3)) == P::constant(K, 3);
        }
        const P::Composer c2(h2, red, 12);
        const P g = P::random(K, n / 2);
        ok = ok && c2.compose(g) == compose_ref(g, h2, m);
        // h constante: las potencias y el paso gigante son constantes.
        const P::Composer c3(P::constant(K, 2), red, 9);
        ok = ok && c3.compose(g) == compose_ref(g, P::constant(K, 2), m);
    }
    CHECK_MSG(ok, (std::string("pow_mod y Composer ") + nombre).c_str());
}

// ---------------------------------------------------------------- 6. RNG
void test_rng() {
    // divmod, inverse_series, Reducer y Composer en un cuerpo recién creado: ni
    // consumen el generador ni construyen tablas (como antes: solo mul_no_rng).
    for (const char* ps : {"998244353", "4179340454199820289", "65537"}) {
        set_seed(777);
        Fp64 K1{BigInt(ps)};
        const P f1 = P::random(K1, 4000), g1 = P::random(K1, 2000), m1 = P::random(K1, 1000);
        const auto [q1, r1] = f1.divmod(g1);
        const P::Reducer red1(m1);
        const P x1 = red1.reduce(P::random(K1, 1998));
        const P::Composer c1(P::random(K1, 999), red1, 32);
        const P y1 = c1.compose(P::random(K1, 999));
        const BigInt s1 = rand_below(BigInt(1000000007));

        set_seed(777);
        Fp64 K2{BigInt(ps)};
        const P f2 = P::random(K2, 4000), g2 = P::random(K2, 2000), m2 = P::random(K2, 1000);
        const P a2 = P::random(K2, 1998), h2 = P::random(K2, 999), gg2 = P::random(K2, 999);
        const BigInt s2 = rand_below(BigInt(1000000007));
        const bool vacias = K1.fft_cache_.empty() && K2.fft_cache_.empty();
        const bool mismos = q1 == f2.divmod_school(g2).first && r1 == f2.divmod_school(g2).second &&
                            x1 == a2.divmod_school(m2).second && y1 == compose_ref(gg2, h2, m2);
        CHECK_MSG(s1 == s2 && vacias && mismos,
                  (std::string("RNG divmod/Reducer/Composer en frío p = ") + ps).c_str());
    }
    // pow_mod en frío: sus productos piden la raíz en el mismo momento que antes.
    for (const char* ps : {"998244353", "7681"}) {
        set_seed(31337);
        Fp64 K1{BigInt(ps)};
        const P m1 = P::random(K1, 400), a1 = P::random(K1, 399);
        const P e1 = P::pow_mod(a1, K1.cardinality(), m1);
        const BigInt s1 = rand_below(BigInt(1000000007));
        set_seed(31337);
        Fp64 K2{BigInt(ps)};
        const P m2 = P::random(K2, 400), a2 = P::random(K2, 399);
        const P e2 = pow_mod_ref(a2, K2.cardinality(), m2);
        const BigInt s2 = rand_below(BigInt(1000000007));
        CHECK_MSG(s1 == s2 && e1.coeffs() == e2.coeffs() && K1.fft_cache_.empty() == K2.fft_cache_.empty(),
                  (std::string("RNG pow_mod en frío p = ") + ps).c_str());
    }
}

}  // namespace

int main() {
    set_seed(20261108);
    Fp64 w998(998244353), w62(BigInt("4179340454199820289")), w65537(65537), w7681(7681), w101(101),
        wm61((BigInt(1) << 61) - 1);
    for (const Fp64* K : {&w998, &w62, &w65537, &w7681, &w101, &wm61}) calienta(*K);
    const std::pair<const Fp64*, const char*> cuerpos[] = {
        {&w998, "F_998244353"}, {&w62, "F_(29·2^57+1)"}, {&w65537, "F_65537"},
        {&w7681, "F_7681"},     {&w101, "F_101"},        {&wm61, "F_(2^61-1)"}};
    for (auto [K, nombre] : cuerpos) {
        test_inversion(*K, nombre);
        test_division(*K, nombre);
        test_reducer(*K, nombre);
        test_pow_compose(*K, nombre);
    }
    test_reducer_frio_caliente();
    test_rng();
    return check::summary("test_fase6_newton");
}
