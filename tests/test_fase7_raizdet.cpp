// test_fase7_raizdet.cpp — Fase 7, propuesta 3d: tablas NTT con raíz
// determinista para las rutas que no consumen el generador (Fp64).
//
//   1. fft::raiz_det tiene orden exactamente 2^capacidad y no depende del
//      cuerpo concreto (dos cuerpos con el mismo p dan la misma raíz).
//   2. En frío (caché FFT vacía) divmod, inverse_series, gcd y Composer dan lo mismo que las versiones escolares/clásicas, NO tocan el
//      generador (el siguiente rand_below coincide con el de un cuerpo en el
//      que no se hace nada) y fft_cache_ sigue vacía; las tablas usadas están
//      en fft_cache_det_.  pow_mod se comprueba aparte (consume el generador
//      como el original).
//   3. Mezcla de raíces: un Reducer y un Composer creados en frío (transformadas
//      guardadas con la raíz determinista) siguen dando el resultado correcto
//      después de que la caché se llene con la raíz del generador.
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <string>
#include <vector>

using namespace alcp;
using P = Poly<Fp64>;

namespace {

const int kGrados[] = {1500, 700, 1500, 1490, 300, 299, 299};
const char* const kPrimos[] = {"998244353", "4179340454199820289", "65537", "7681", "97"};

// g(h) mod m por Horner con el % de siempre.
P compone_ref(const P& g, const P& h, const P& m) {
    const Fp64& K = m.field();
    P acc = P::zero(K);
    const auto c = g.coeffs();
    for (std::size_t i = c.size(); i-- > 0;) acc = acc.mul_school(h).divmod_school(m).second + P::constant(K, c[i]);
    return acc.divmod_school(m).second;
}

// ------------------------------------------------------------ 1. raíz
void test_raiz() {
    for (const char* ps : kPrimos) {
        const Fp64 K{BigInt(ps)}, K2{BigInt(ps)};
        const unsigned long cap = K.fft_capacity();
        const auto r = fft::raiz_det(K);
        // Orden exactamente 2^cap: r^(2^cap) = 1 y r^(2^(cap-1)) = -1.
        const bool ok = K.pow(r, BigInt(1) << cap) == K.one() &&
                        K.pow(r, BigInt(1) << (cap - 1)) == K.neg(K.one()) && fft::raiz_det(K2) == r;
        CHECK_MSG(ok, (std::string("raiz_det de orden 2^cap, p = ") + ps).c_str());
    }
}

// ------------------------------------------------------------ 2. en frío
void test_frio(const char* ps) {
    const std::string tag = std::string(" p = ") + ps;
    // Referencia del generador: misma semilla, mismos sorteos, nada más.
    set_seed(777);
    {
        const Fp64 R{BigInt(ps)};
        for (int d : kGrados) (void)P::random(R, d);
    }
    const BigInt sig_ref = rand_below(BigInt(1000000007));

    set_seed(777);
    const Fp64 K{BigInt(ps)};
    // Mismos sorteos que la referencia (kGrados).
    const P f = P::random(K, 1500), g = P::random(K, 700), a = P::random(K, 1500), b = P::random(K, 1490);
    const P m = P::random(K, 300), h = P::random(K, 299), gc = P::random(K, 299);
    bool ok = true;
    // divmod (Newton aislado) frente a la escolar.
    const auto [q, r] = f.divmod(g);
    const auto [qs, rs] = f.divmod_school(g);
    ok = ok && q == qs && r == rs;
    CHECK_MSG(ok, ("divmod en frío = escolar" + tag).c_str());
    // inverse_series: a·x ≡ 1 mod x^n.
    auto serie = a.coeffs();
    serie[0] = K.one();  // a[0] != 0 (con p = 97 el sorteo puede dar 0)
    const auto inv = P::inverse_series(K, serie, 1000);
    const P prod = P(K, serie).mul_school(P(K, inv));
    bool ok_inv = prod.coeffs().size() >= 1000 && prod.coeffs()[0] == K.one();
    for (std::size_t i = 1; ok_inv && i < 1000; ++i) ok_inv = prod.coeffs()[i] == K.zero();
    CHECK_MSG(ok_inv, ("inverse_series en frío" + tag).c_str());
    // gcd (half-GCD) frente al clásico.
    CHECK_MSG(P::gcd(a, b) == P::gcd_classic(a, b), ("gcd en frío = clásico" + tag).c_str());
    // Composer frente a la referencia escolar.
    const P::Reducer red(m);
    const P::Composer comp(h, red, 18);
    CHECK_MSG(comp.compose(gc) == compone_ref(gc, h, m), ("Composer en frío" + tag).c_str());

    // Nada de lo anterior consume el generador ni llena fft_cache_ (las
    // referencias usan solo mul_school y divmod_school).
    const bool vacia = K.fft_cache_.empty();
    const BigInt sig = rand_below(BigInt(1000000007));
    CHECK_MSG(vacia, ("fft_cache_ sigue vacía" + tag).c_str());
    CHECK_MSG(sig == sig_ref, ("el generador no se consume en frío" + tag).c_str());
    // pow_mod va al final: multiplica con mul (regla Auto), que construye las
    // tablas con el generador igual que el original (ver test_fase6_newton).
    const P pm = P::pow_mod(h, BigInt(1000003), m);
    P pr = P::constant(K, K.one()), base = h.divmod_school(m).second;
    for (BigInt e(1000003); e > 0; e >>= 1) {
        if (mpz_odd_p(e.get_mpz_t())) pr = pr.mul_school(base).divmod_school(m).second;
        base = base.mul_school(base).divmod_school(m).second;
    }
    CHECK_MSG(pm == pr, ("pow_mod en frío" + tag).c_str());
}

// ------------------------------------------------------- 3. mezcla de raíces
void test_mezcla(const char* ps) {
    const std::string tag = std::string(" p = ") + ps;
    set_seed(99);
    const Fp64 K{BigInt(ps)};
    // Hace falta NTT de tamaño 2^12 (productos de grado ~1400); con p = 7681
    // o 97 no la hay y no hay transformadas que mezclar.
    if (K.fft_capacity() < 12) return;
    const P m = P::random(K, 700), h = P::random(K, 699);
    const P a = P::random(K, 1398), b = P::random(K, 1350), c = P::random(K, 1000);
    const P::Reducer red(m);
    bool ok = red.reduce(a) == a.divmod_school(m).second;
    const P::Composer comp(h, red, 30);
    const P g1 = P::random(K, 699);
    ok = ok && comp.compose(g1) == compone_ref(g1, h, m);
    ok = ok && K.fft_cache_.empty() && !K.fft_cache_det_.empty();
    // La caché se llena con la raíz del generador.
    (void)fft::tables(K, 12);
    ok = ok && !K.fft_cache_.empty();
    // Las transformadas guardadas (raíz determinista) y las nuevas (raíz del
    // generador) conviven: cada una se combina con sus tablas.
    ok = ok && red.reduce(b) == b.divmod_school(m).second && red.reduce(a) == a.divmod_school(m).second;
    ok = ok && red.reduce(c) == c.divmod_school(m).second;
    const P g2 = P::random(K, 699);
    ok = ok && comp.compose(g2) == compone_ref(g2, h, m) && comp.compose(g1) == compone_ref(g1, h, m);
    const P::Composer comp2(h, red, 30);  // nuevo, ya con la caché llena
    ok = ok && comp2.compose(g2) == compone_ref(g2, h, m);
    CHECK_MSG(ok, ("Reducer y Composer frío -> caliente" + tag).c_str());
}

}  // namespace

int main() {
    test_raiz();
    for (const char* ps : kPrimos) {
        test_frio(ps);
        test_mezcla(ps);
    }
    return check::summary("test_fase7_raizdet");
}
