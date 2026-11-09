// test_fase8_cabos.cpp — Fase 8, propuesta 6: cabos sueltos de la revisión de
// la fase 7.
//
//   1. fft::convolve (Fp64) rechaza coeficientes no canónicos (>= p) con
//      std::invalid_argument, ANTES de pedir las tablas: el generador no se
//      consume y fft_cache_ sigue vacía.  Con entradas canónicas, el mismo
//      resultado que el producto escolar.
//   2. Las tablas deterministas (fft_cache_det_) se liberan en cuanto la caché
//      normal tiene tablas y se vuelve a pedir una tabla sin generador; un
//      Reducer y un Composer creados en frío (que guardan sus tablas en
//      medio::Fija::tab) siguen dando el resultado correcto después.
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <stdexcept>
#include <string>
#include <vector>

using namespace alcp;
using P = Poly<Fp64>;
using V = std::vector<std::uint64_t>;

namespace {

// g(h) mod m por Horner, solo con mul_school y divmod_school.
P compone_ref(const P& g, const P& h, const P& m) {
    const Fp64& K = m.field();
    P acc = P::zero(K);
    const auto c = g.coeffs();
    for (std::size_t i = c.size(); i-- > 0;) acc = acc.mul_school(h).divmod_school(m).second + P::constant(K, c[i]);
    return acc.divmod_school(m).second;
}

void test_convolve_canonico() {
    for (const char* ps : {"998244353", "65537", "4179340454199820289"}) {
        const std::string tag = std::string(" p = ") + ps;
        set_seed(11);
        const Fp64 R{BigInt(ps)};
        const BigInt sig_ref = rand_below(BigInt(1000000007));

        set_seed(11);
        const Fp64 K{BigInt(ps)};
        const std::uint64_t p = K.modulus64();
        bool lanza = false;
        try {
            (void)fft::convolve(K, V{1, p, 2}, V{3, 4});
        } catch (const std::invalid_argument&) {
            lanza = true;
        }
        bool lanza_b = false;
        try {
            (void)fft::convolve(K, V{1, 2}, V{p + 5});
        } catch (const std::invalid_argument&) {
            lanza_b = true;
        }
        const bool vacia = K.fft_cache_.empty();
        CHECK_MSG(lanza && lanza_b, ("convolve rechaza coeficientes >= p" + tag).c_str());
        CHECK_MSG(vacia && rand_below(BigInt(1000000007)) == sig_ref,
                  ("una entrada inválida no consume el generador" + tag).c_str());
        // Entradas canónicas (incluido p - 1): igual que el escolar.
        V a(300), b(257);
        for (std::size_t i = 0; i < a.size(); ++i) a[i] = (i % 7 == 0) ? p - 1 : (i * 2654435761u) % p;
        for (std::size_t i = 0; i < b.size(); ++i) b[i] = (i % 5 == 0) ? p - 1 : (i * 40503u + 7) % p;
        const P ref = P(K, a).mul_school(P(K, b));
        CHECK_MSG(P(K, fft::convolve(K, a, b)) == ref, ("convolve canónico = escolar" + tag).c_str());
    }
}

void test_liberacion() {
    for (const char* ps : {"998244353", "4179340454199820289", "65537"}) {
        const std::string tag = std::string(" p = ") + ps;
        set_seed(21);
        const Fp64 K{BigInt(ps)};
        const P m = P::random(K, 700), h = P::random(K, 699);
        const P a = P::random(K, 1398), b = P::random(K, 1350), g1 = P::random(K, 699);
        const P::Reducer red(m);
        const P::Composer comp(h, red, 30);
        bool ok = red.reduce(a) == a.divmod_school(m).second && comp.compose(g1) == compone_ref(g1, h, m);
        ok = ok && K.fft_cache_.empty() && !K.fft_cache_det_.empty();
        CHECK_MSG(ok, ("en frío se usan tablas deterministas" + tag).c_str());
        // Se llena la caché normal y se pide una tabla sin generador: las
        // deterministas se liberan.
        (void)fft::tables(K, 12);
        (void)fft::tablas_sin_rng(K, 11);
        CHECK_MSG(!K.fft_cache_.empty() && K.fft_cache_det_.empty(),
                  ("caché determinista liberada al llenarse la normal" + tag).c_str());
        // Las transformadas guardadas en frío siguen siendo válidas.
        ok = red.reduce(b) == b.divmod_school(m).second && red.reduce(a) == a.divmod_school(m).second;
        const P g2 = P::random(K, 699);
        ok = ok && comp.compose(g2) == compone_ref(g2, h, m) && comp.compose(g1) == compone_ref(g1, h, m);
        CHECK_MSG(ok, ("Reducer y Composer de frío tras liberar" + tag).c_str());
    }
}

}  // namespace

int main() {
    test_convolve_canonico();
    test_liberacion();
    return check::summary("test_fase8_cabos");
}
