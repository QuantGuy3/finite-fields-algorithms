// test_fase8_reducer.cpp — Fase 8, propuesta 1: reducción por tabla de
// x^(n+i) mod m en el Reducer de Poly<Fp64> (p <= 2^16, sin NTT para el
// tamaño, con AVX2; namespace tabla de polynomial.hpp).
//
//   1. tabla::reduce directamente frente a divmod_school: grados de m límite
//      (2^j y 2^j ± 1, kMin y kMax), deg a = n..2n-1, módulos no mónicos,
//      coeficientes p - 1 (sumas máximas) y la frontera entre las sumas de 32
//      y de 64 bits (p = 4093, n = 256 / 257; p = 4099, n = 255 / 256).
//   1b. Revisión: suma máxima sintética en esa frontera (tabla y a con todas
//      las entradas p - 1) frente al valor exacto, y las versiones escalares
//      (crea_escalar, combina_escalar: las que se usan sin AVX2) frente a las
//      vectoriales y a divmod_school.
//   2. Reducer::reduce frente a divmod_school con el mismo Reducer muchas
//      veces (la tabla se construye a partir de la llamada 2 + n/128) y con
//      sus copias (comparten la tabla): deg a = n-1..2n-1, 2n, 2n+5 y cero.
//      Incluye cuerpos sin tabla (F_65537, F_(2^61-1)) y con NTT pequeña
//      (F_257, F_193: unas reducciones por la NTT y otras por la tabla).
//   2b. Revisión: copias de un Reducer usadas a la vez desde varios hilos.
//   3. pow_mod y Composer::compose frente a referencias escolares.
//   4. En frío: Reducer, pow_mod (sin FFT en el cuerpo) y Composer no
//      consumen el generador y dejan fft_cache_ vacía.
//   5. Diferencial frente a 04c1909: una huella (FNV-1a) de factor(),
//      is_irreducible(), pow_mod, Composer, el siguiente rand_below y
//      fft_cache_.empty() en F_2, F_3, F_7, F_101, F_65521, F_65537 (en frío),
//      F_257 y F_(2^61-1).  Las constantes se obtuvieron compilando ESTE
//      fichero con -DALCP_HUELLAS (solo esta sección, que imprime las huellas)
//      contra las cabeceras de 04c1909.
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <cstdint>
#include <cstdio>
#include <string>
#include <thread>
#include <vector>

using namespace alcp;
using P = Poly<Fp64>;

namespace {

// a con coeficientes aleatorios de grado exactamente d, o todos p - 1.
P aleatorio(const Fp64& K, int d, bool maximo) {
    if (!maximo) return P::random(K, d);
    return P(K, std::vector<std::uint64_t>(static_cast<std::size_t>(d) + 1, K.modulus64() - 1));
}

// g(h) mod m por Horner con el producto y la división escolares.
P compone_ref(const P& g, const P& h, const P& m) {
    const Fp64& K = m.field();
    P acc = P::zero(K);
    const auto c = g.coeffs();
    for (std::size_t i = c.size(); i-- > 0;) acc = acc.mul_school(h).divmod_school(m).second + P::constant(K, c[i]);
    return acc.divmod_school(m).second;
}

// a^k mod m por cuadrados con el producto y la división escolares.
P pow_ref(P a, BigInt k, const P& m) {
    const Fp64& K = m.field();
    P r = P::one(K).divmod_school(m).second;
    a = a.divmod_school(m).second;
    while (k > 0) {
        if (mpz_odd_p(k.get_mpz_t())) r = r.mul_school(a).divmod_school(m).second;
        k >>= 1;
        if (k > 0) a = a.mul_school(a).divmod_school(m).second;
    }
    return r;
}

#ifndef ALCP_HUELLAS
// ------------------------------------------------------------ 1. tabla directa
void test_tabla_directa() {
    struct Caso { const char* p; std::vector<int> ns; };
    const Caso casos[] = {
        {"2", {8, 9, 15, 16, 17, 64, 100, 384}},
        {"101", {8, 9, 31, 32, 33, 63, 64, 65, 100, 127, 128, 129, 255, 256, 257, 383, 384}},
        {"4093", {255, 256, 257}},  // frontera de las sumas de 32 bits
        {"4099", {255, 256}},       // (4099: v32 con 255, v64 con 256)
        {"65521", {8, 33, 100, 128, 129, 384}},
    };
    for (const Caso& c : casos) {
        set_seed(81);
        const Fp64 K{BigInt(c.p)};
        for (int n : c.ns) {
            const std::size_t nn = static_cast<std::size_t>(n);
            if (!tabla::aplicable(K, nn)) {
                // Sin AVX2 (o a -O0) la tabla no se usa: nada que comprobar.
                CHECK_MSG(!fft::word::avx2_disponible(), "tabla aplicable con AVX2");
                continue;
            }
            for (int maximo = 0; maximo < 2; ++maximo) {
                // Módulo no mónico (lc = 2 si p > 2) y, con maximo, coeficientes p - 1.
                P m = aleatorio(K, n, maximo != 0);
                if (K.modulus64() > 2) m = m.monic().scalar_mul(K.add(K.one(), K.one()));
                const tabla::Tabla t = tabla::crea(K, m.coeffs());
                bool ok = !t.vacia() && t.n == nn;
                for (int da : {n, n + 1, n + n / 2, 2 * n - 2, 2 * n - 1}) {
                    const P a = aleatorio(K, da, maximo != 0);
                    ok = ok && P(K, tabla::reduce(K, t, a.coeffs())) == a.divmod_school(m).second;
                }
                CHECK_MSG(ok, ("tabla::reduce = escolar, p = " + std::string(c.p) + ", n = " + std::to_string(n) +
                               (maximo ? " (p-1)" : ""))
                                  .c_str());
            }
        }
    }
    // Fuera de rango: grados < kMin o > kMax y p > 2^16.
    const Fp64 K101{BigInt(101)}, K65537{BigInt(65537)};
    CHECK(!tabla::aplicable(K101, tabla::kMin - 1));
    CHECK(!tabla::aplicable(K101, tabla::kMax + 1));
    CHECK(!tabla::aplicable(K65537, 100));
}

// ------------------------------------------------------------ 1b. frontera y escalares
void test_frontera() {
    struct Caso { const char* p; int n; bool v32; };
    // n·(p-1)² + (p-1) < 2^32: 4093 hasta n = 256, 4099 hasta n = 255.
    const Caso casos[] = {{"4093", 255, true}, {"4093", 256, true}, {"4093", 257, false},
                          {"4099", 255, true}, {"4099", 256, false}};
    for (const Caso& c : casos) {
        const std::string tag = std::string(" p = ") + c.p + ", n = " + std::to_string(c.n);
        set_seed(17);
        const Fp64 K{BigInt(c.p)};
        const std::uint64_t p = K.modulus64();
        const std::size_t n = static_cast<std::size_t>(c.n);
        // Tabla real: v32 según la cota, y la escalar idéntica a la vectorial.
        const P m = P::random(K, c.n);
        const tabla::Tabla t = tabla::crea(K, m.coeffs()), te = tabla::crea_escalar(K, m.coeffs());
        CHECK_MSG(t.v32 == c.v32 && te.v32 == c.v32, ("cota de 32 bits" + tag).c_str());
        CHECK_MSG(t.buf == te.buf, ("crea = crea_escalar" + tag).c_str());
        // Sintética: todas las entradas de la tabla y de a valen p - 1, la
        // suma máxima n·(p-1)² + (p-1) en cada columna.  Valor exacto con u128.
        tabla::Tabla s = t;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) s.buf[8 + i * s.ld + j] = static_cast<std::uint32_t>(p - 1);
        const std::vector<std::uint64_t> a(2 * n, p - 1);
        const unsigned __int128 exacto = static_cast<unsigned __int128>(n) * (p - 1) * (p - 1) + (p - 1);
        const std::uint64_t esperado = static_cast<std::uint64_t>(exacto % p);
        const auto r = tabla::reduce(K, s, a);
        const auto re = tabla::combina_escalar(K, s, a.data(), a.data() + n, n);
        bool ok = r.size() == n && r == re;
        for (auto x : r) ok = ok && x == esperado;
        CHECK_MSG(ok, ("suma máxima sintética" + tag).c_str());
        // Tabla real con a de coeficientes p - 1 y grados n..2n-1, frente a la
        // escolar, por la vectorial y por la escalar.
        bool ok2 = true;
        for (int da : {c.n, c.n + 1, 2 * c.n - 1}) {
            const P A = aleatorio(K, da, true);
            const P esp = A.divmod_school(m).second;
            const auto& ac = A.coeffs();
            ok2 = ok2 && P(K, tabla::reduce(K, t, ac)) == esp &&
                  P(K, tabla::combina_escalar(K, te, ac.data(), ac.data() + n, ac.size() - n)) == esp;
        }
        CHECK_MSG(ok2, ("tabla real, a con p - 1, = escolar" + tag).c_str());
    }
}

// ------------------------------------------------------------ 2. Reducer
void test_reducer(const char* ps, const std::vector<int>& ns) {
    set_seed(5);
    const Fp64 K{BigInt(ps)};
    for (int n : ns) {
        const P m = P::random(K, n);
        const P::Reducer red(m);
        const P::Reducer copia = red;  // comparte el estado de la tabla
        std::vector<P> as;
        for (int da = n - 1; da <= 2 * n - 1; da += std::max(1, n / 7)) as.push_back(P::random(K, da));
        as.push_back(P::random(K, 2 * n - 1));
        as.push_back(aleatorio(K, 2 * n - 1, true));
        as.push_back(P::random(K, 2 * n));       // d1 > n: división escolar
        as.push_back(P::random(K, 2 * n + 5));
        as.push_back(P::zero(K));
        bool ok = true;
        // Tres vueltas: la primera sin tabla, después con ella (y por la copia).
        for (int vuelta = 0; vuelta < 3; ++vuelta)
            for (const P& a : as) {
                const P esperado = a.divmod_school(m).second;
                ok = ok && red.reduce(a) == esperado && copia.reduce(a) == esperado;
            }
        CHECK_MSG(ok, ("Reducer::reduce = escolar, p = " + std::string(ps) + ", n = " + std::to_string(n)).c_str());
    }
}

// ------------------------------------------------------------ 2b. hilos
// Revisión: copias de un Reducer recién creado (comparten el estado de la
// tabla) usadas a la vez desde varios hilos; todas deben construir/usar la
// misma tabla sin carreras (contador atómico y std::call_once).
void test_hilos() {
    set_seed(31);
    const Fp64 K{BigInt(101)};
    for (int n : {16, 100, 300}) {
        const P m = P::random(K, n);
        std::vector<P> as, esperados;
        for (int i = 0; i < 24; ++i) {
            as.push_back(P::random(K, n + (i * n) / 24));
            esperados.push_back(as.back().divmod_school(m).second);
        }
        const P::Reducer red(m);
        std::vector<P::Reducer> copias(6, red);
        std::vector<int> fallos(copias.size(), 0);
        std::vector<std::thread> hilos;
        for (std::size_t h = 0; h < copias.size(); ++h)
            hilos.emplace_back([&, h] {
                for (int v = 0; v < 4; ++v)
                    for (std::size_t i = 0; i < as.size(); ++i)
                        if (!(copias[h].reduce(as[i]) == esperados[i])) ++fallos[h];
            });
        for (auto& t : hilos) t.join();
        int total = 0;
        for (int x : fallos) total += x;
        CHECK_MSG(total == 0, ("Reducer compartido entre hilos, n = " + std::to_string(n)).c_str());
    }
}

// ------------------------------------------------------------ 3. pow_mod y Composer
void test_pow_compose(const char* ps, const std::vector<int>& ns) {
    set_seed(9);
    const Fp64 K{BigInt(ps)};
    for (int n : ns) {
        const P m = P::random(K, n);
        const P::Reducer red(m);
        const P a = P::random(K, n - 1), h = P::random(K, n - 1), g = P::random(K, n - 1);
        bool ok = P::pow_mod(a, K.cardinality(), red) == pow_ref(a, K.cardinality(), m);
        ok = ok && P::pow_mod(a, BigInt(1000003), m) == pow_ref(a, BigInt(1000003), m);
        const P::Composer comp(h, red, 7);
        ok = ok && comp.compose(g) == compone_ref(g, h, m) && comp.compose(a) == compone_ref(a, h, m);
        CHECK_MSG(ok, ("pow_mod y Composer, p = " + std::string(ps) + ", n = " + std::to_string(n)).c_str());
    }
}

// ------------------------------------------------------------ 4. generador en frío
void test_frio(const char* ps, bool con_pow) {
    const std::string tag = std::string(" p = ") + ps;
    set_seed(4242);
    {
        const Fp64 R{BigInt(ps)};
        for (int d : {100, 99, 99, 99, 200}) (void)P::random(R, d);
    }
    const BigInt sig_ref = rand_below(BigInt(1000000007));

    set_seed(4242);
    const Fp64 K{BigInt(ps)};
    const P m = P::random(K, 100), a = P::random(K, 99), h = P::random(K, 99), g = P::random(K, 99);
    const P b = P::random(K, 200);
    const P::Reducer red(m);
    bool ok = true;
    for (int i = 0; i < 4; ++i) ok = ok && red.reduce(b) == b.divmod_school(m).second;
    const P::Composer comp(h, red, 10);
    ok = ok && comp.compose(g) == compone_ref(g, h, m);
    if (con_pow) ok = ok && P::pow_mod(a, K.cardinality(), red) == pow_ref(a, K.cardinality(), m);
    CHECK_MSG(ok, ("Reducer/Composer en frío" + tag).c_str());
    const bool vacia = K.fft_cache_.empty();
    CHECK_MSG(vacia, ("fft_cache_ sigue vacía" + tag).c_str());
    CHECK_MSG(rand_below(BigInt(1000000007)) == sig_ref, ("generador intacto" + tag).c_str());
}

#endif  // ALCP_HUELLAS

// ------------------------------------------------------------ 5. huellas
struct Huella {
    std::uint64_t h = 1469598103934665603ull;
    void mete(std::uint64_t x) {
        for (int i = 0; i < 8; ++i) { h ^= (x >> (8 * i)) & 0xff; h *= 1099511628211ull; }
    }
    void mete(const P& f) {
        mete(static_cast<std::uint64_t>(f.coeffs().size()));
        for (auto c : f.coeffs()) mete(c);
    }
    void mete(const BigInt& x) { mete(static_cast<std::uint64_t>(mpz_getlimbn(x.get_mpz_t(), 0))); }
};

std::uint64_t huella(const char* ps, const std::vector<int>& grados) {
    set_seed(2026);
    const Fp64 K{BigInt(ps)};  // en frío
    Huella H;
    for (int n : grados) {
        P f = P::random(K, n);
        // Un factor repetido de vez en cuando (multiplicidades y radical).
        if (n % 3 == 0) f = f * P::random(K, 2) * P::random(K, 2);
        const auto F = f.factor();
        H.mete(F.leading);
        H.mete(static_cast<std::uint64_t>(F.factors.size()));
        for (const auto& [u, e] : F.factors) { H.mete(u); H.mete(static_cast<std::uint64_t>(e)); }
        H.mete(static_cast<std::uint64_t>(f.is_irreducible()));
        H.mete(static_cast<std::uint64_t>(f.monic().is_irreducible()));
        const P m = P::random(K, n), a = P::random(K, n - 1), h = P::random(K, n - 1);
        const P::Reducer red(m);
        H.mete(P::pow_mod(a, K.cardinality(), red));
        H.mete(P::Composer(h, red, 6).compose(a));
        H.mete(rand_below(BigInt(1000000007)));
        H.mete(static_cast<std::uint64_t>(K.fft_cache_.empty()));
    }
    return H.h;
}

void test_huellas() {
    struct Caso { const char* p; std::vector<int> grados; std::uint64_t esperada; };
    // Constantes de 04c1909 (ver la cabecera del fichero).
    const Caso casos[] = {
        {"2", {8, 16, 33, 64, 100, 129, 200}, 0xa8f691b935e7e9b9ull},
        {"3", {8, 17, 32, 65, 100, 150}, 0xd50375f5b991a9e0ull},
        {"7", {9, 16, 40, 64, 100, 128, 129, 300}, 0x06d382429673df06ull},
        {"101", {8, 15, 16, 32, 40, 63, 64, 100, 128, 129, 200, 384}, 0xbc158a51da45ed9aull},
        {"65521", {8, 33, 64, 100}, 0x9a938f2f15cb5d5dull},
        {"65537", {8, 33, 64, 100}, 0x33348958ef13e4f6ull},
        {"257", {8, 64, 100, 130}, 0x753b09585fc00965ull},
        {"2305843009213693951", {8, 33, 64}, 0xe31004ae43d54296ull},
    };
    for (const Caso& c : casos) {
        const std::uint64_t h = huella(c.p, c.grados);
#ifdef ALCP_HUELLAS
        std::printf("huella p = %s: 0x%016llxull\n", c.p, static_cast<unsigned long long>(h));
#endif
        CHECK_MSG(h == c.esperada, ("huella igual a 04c1909, p = " + std::string(c.p)).c_str());
    }
}

}  // namespace

int main() {
#ifndef ALCP_HUELLAS
    test_tabla_directa();
    test_frontera();
    test_reducer("2", {8, 9, 16, 17, 64, 100});
    test_reducer("3", {8, 31, 32, 33, 100, 129});
    test_reducer("7", {7, 8, 9, 15, 16, 17, 63, 64, 65, 100, 255, 256, 257, 383, 384, 385});
    test_reducer("101", {1, 2, 7, 8, 9, 31, 32, 33, 63, 64, 65, 79, 80, 81, 100, 127, 128, 129, 255, 256, 257,
                         383, 384, 385, 511, 512, 513});
    test_reducer("4093", {255, 256, 257});
    test_reducer("65521", {8, 64, 100, 129, 384});
    test_reducer("257", {32, 64, 65, 100, 128, 129, 200});
    test_reducer("193", {32, 33, 64, 100});
    test_reducer("65537", {32, 100, 129});
    test_reducer("2305843009213693951", {32, 100});
    test_hilos();
    test_pow_compose("101", {8, 33, 64, 100, 129, 200});
    test_pow_compose("7", {16, 100, 256});
    test_pow_compose("65521", {33, 100});
    test_pow_compose("2", {64, 100});
    test_frio("101", true);
    test_frio("7", true);
    test_frio("2305843009213693951", true);
    test_frio("65537", false);
    test_frio("257", false);
    test_frio("193", false);
    test_frio("998244353", false);
#endif
    test_huellas();
    return check::summary("test_fase8_reducer");
}
