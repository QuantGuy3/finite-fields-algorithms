// test_fase8_fq.cpp — Fase 8: productos de FqSwar (NTT módulo 998244353,
// huecos de 16 bits alineados, ancho justo, reducción módulo g con tabla y
// umbral de Kronecker 8) y constantes de FqTab.  Todo debe dar EXACTAMENTE lo
// mismo que antes:
//
//   1. Productos sueltos (regla Auto, Kronecker, escolar, Karatsuba y
//      cuadrados) de muchas formas (equilibrados, desequilibrados, en los
//      bordes de cada camino de kronecker.hpp: m16_, kMinHuecos, kMinM) frente
//      al Kronecker de Fq64, cuyo código no cambia en la fase 8.
//   2. Batería de polinomios en frío (cuerpo recién creado, misma semilla):
//      productos, divmod (escolar y Newton), gcd, gcd_ext, pow_mod, compose,
//      is_irreducible y factor, con el MISMO estado del generador después
//      (rand_below) y la misma caché FFT (vacía o no) que Fq64, en cuerpos con
//      FFT útil (F_(3^16) y F_(3^32) en FqSwar, F_(97^2) en FqTab) y sin ella.
//   3. Huella: un resumen de todos los resultados que se imprime para poder
//      compilar este mismo programa contra las cabeceras de 04c1909 y
//      comparar (debe salir la misma huella).
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <array>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

using namespace alcp;
using V = Fq64::Elem;

namespace {

std::uint64_t huella = 1469598103934665603ULL;  // FNV-1a de todos los resultados
void mezcla(std::uint64_t x) {
    for (int i = 0; i < 8; ++i) {
        huella ^= (x >> (8 * i)) & 0xFF;
        huella *= 1099511628211ULL;
    }
}

V a_vec(const FqTab& K, FqTab::Elem e) { return K.decode(e); }
V a_vec(const FqSwar&, const V& e) { return e; }
V a_vec(const Fq64&, const V& e) { return e; }

template <class K>
std::vector<V> pvec(const K& F, const Poly<K>& p) {
    std::vector<V> r;
    for (const auto& c : p.coeffs()) r.push_back(a_vec(F, c));
    for (const auto& v : r)
        for (auto x : v) mezcla(x);
    mezcla(r.size());
    return r;
}
BigInt testigo() { return rand_below(BigInt(1) << 64); }

std::vector<BigInt> big(const std::vector<long>& v) {
    std::vector<BigInt> r;
    for (long x : v) r.emplace_back(x);
    return r;
}

// Polinomio aleatorio con na coeficientes y el principal no nulo, en Fq64.
std::vector<V> aleatorio(const Fq64& R, std::size_t na) {
    std::vector<V> c;
    for (std::size_t i = 0; i < na; ++i) c.push_back(R.random());
    if (R.is_zero(c.back())) c.back() = R.one();
    return c;
}
template <class K>
Poly<K> en(const K& F, const std::vector<V>& c) {
    std::vector<typename K::Elem> r;
    for (const auto& x : c) {
        if constexpr (std::is_same_v<K, FqTab>) r.push_back(F.encode(x));
        else r.push_back(x);
    }
    return Poly<K>(F, std::move(r));
}

// ------------------------------------------------------- 1. productos sueltos
template <class K>
void productos(const char* nombre, const Fq64& R, const K& F, const std::vector<std::pair<std::size_t, std::size_t>>& formas,
               bool escolar = true) {
    bool ok = true;
    set_seed(41);
    for (auto [na, nb] : formas) {
        const auto ca = aleatorio(R, na), cb = aleatorio(R, nb);
        const Poly<Fq64> A(R, ca), B(R, cb);
        const auto ref = A.mul(B, MulAlgo::Kronecker).coeffs();
        const auto ref2 = A.mul(A, MulAlgo::Kronecker).coeffs();
        const Poly<K> a = en(F, ca), b = en(F, cb);
        bool bien = pvec(F, a * b) == ref && pvec(F, b * a) == ref && pvec(F, a.mul(b, MulAlgo::Kronecker)) == ref &&
                    pvec(F, a * a) == ref2 && pvec(F, a.mul(a, MulAlgo::Kronecker)) == ref2;
        if (escolar && na * nb <= 70000)
            bien = bien && pvec(F, a.mul(b, MulAlgo::Schoolbook)) == ref && pvec(F, a.mul(b, MulAlgo::Karatsuba)) == ref;
        if (!bien) std::printf("    %s: falla %zu x %zu\n", nombre, na, nb);
        ok = ok && bien;
    }
    CHECK_MSG(ok, (std::string(nombre) + ": productos").c_str());
}

// --------------------------------------------------------- 2. batería en frío
template <class F>
struct Calc {
    using P = Poly<F>;
    static std::vector<V> num(long x) { return {V(1, static_cast<std::uint64_t>(x))}; }
    static void factor(const F& K, const P& f, std::vector<std::vector<V>>& out) {
        auto r = f.factor();
        out.push_back({a_vec(K, r.leading)});
        for (auto& [u, e] : r.factors) {
            out.push_back(pvec(K, u));
            out.push_back(num(e));
        }
    }
    static std::vector<std::vector<V>> bateria(const F& K, int d1, int d2, int df) {
        std::vector<std::vector<V>> out;
        const P a = P::random(K, d1), b = P::random(K, d2);
        out.push_back(pvec(K, a * b));
        out.push_back(pvec(K, a * a));
        out.push_back(pvec(K, a.mul(b, MulAlgo::Kronecker)));
        auto [qq, rr] = a.divmod(b);
        out.push_back(pvec(K, qq));
        out.push_back(pvec(K, rr));
        out.push_back(pvec(K, P::gcd(a, b)));
        auto e = P::gcd_ext(a, b);
        out.push_back(pvec(K, e.g));
        out.push_back(pvec(K, e.x));
        out.push_back(pvec(K, e.y));
        out.push_back(pvec(K, P::pow_mod(a, K.cardinality() + 3, b)));
        out.push_back(pvec(K, a.compose(P::random(K, 3)) % b));
        out.push_back(num(b.is_irreducible() ? 1 : 0));
        out.push_back(num(P::random(K, df).monic().is_irreducible() ? 1 : 0));
        factor(K, P::random(K, df), out);
        const P g = P::random(K, 2) * P::random(K, 2) * P::random(K, 1).pow(BigInt(3)) * P::random(K, df / 2 + 1);
        factor(K, g, out);
        out.push_back(pvec(K, g.radical()));
        return out;
    }
};

template <class K>
void frio(const char* nombre, const std::function<Fq64()>& crear, const std::function<K(Fq64)>& envolver,
          const std::vector<std::array<int, 3>>& casos) {
    unsigned long semilla = 500;
    for (auto [d1, d2, df] : casos) {
        const std::string t = std::string(nombre) + " " + std::to_string(d1) + "/" + std::to_string(d2) + "/" +
                              std::to_string(df);
        const Fq64 R = crear();
        const K F = envolver(crear());
        set_seed(++semilla);
        const auto ra = Calc<Fq64>::bateria(R, d1, d2, df);
        const BigInt ta = testigo();
        set_seed(semilla);
        const auto rb = Calc<K>::bateria(F, d1, d2, df);
        const BigInt tb = testigo();
        CHECK_MSG(ra == rb, (t + ": resultado").c_str());
        CHECK_MSG(ta == tb, (t + ": generador después").c_str());
        CHECK_MSG(R.fft_cache_.empty() == F.fft_cache_.empty(), (t + ": caché FFT").c_str());
    }
}

}  // namespace

int main() {
    using FT = std::function<FqTab(Fq64)>;
    using FS = std::function<FqSwar(Fq64)>;
    const FT tab = [](Fq64 r) { return FqTab(std::move(r)); };
    const FS swar = [](Fq64 r) { return FqSwar(std::move(r)); };
    auto mk = [](long p, std::vector<long> g) {
        return std::function<Fq64()>([p, g] { return Fq64(Fp64(BigInt(p)), big(g)); });
    };
    auto mkq = [](long p, unsigned long n, unsigned long semilla) {
        set_seed(semilla);
        const Fq64 R(BigInt(p), n);
        return std::function<Fq64()>([R] { return Fq64(R.base_field(), R.modulus()); });
    };
    const auto f3e30 = mk(3, {2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1});
    const auto f7e4 = mk(7, {3, 4, 1, 2, 1});
    const auto aes = mk(2, {1, 1, 0, 1, 1, 0, 0, 0, 1});
    const auto f3 = mk(3, {1, 1});              // n = 1
    const auto f97e2 = mk(97, {5, 0, 1});       // capacidad FFT 6
    const auto f5e10 = mkq(5, 10, 1);
    const auto f2e17 = mkq(2, 17, 2);
    const auto f17e3 = mkq(17, 3, 3);   // m16_ = 254: huecos de 16 bits hasta 84 coeficientes
    const auto f127e2 = mkq(127, 2, 4);  // m16_ = 3: casi todo por el ancho justo
    const auto f23e64 = mkq(23, 64, 5);  // n máximo y el margen de 16 bits más justo
    const auto f2e64 = mkq(2, 64, 6);
    const auto f3e16 = mkq(3, 16, 9);   // capacidad FFT 6
    const auto f3e32 = mkq(3, 32, 9);   // capacidad FFT 7
    const auto f2e12 = mkq(2, 12, 7);
    const auto f3e7 = mkq(3, 7, 8);
    const auto f3e5 = mkq(3, 5, 10);

    // Formas: escolar (< 8), huecos de 16 bits, ancho justo, NTT (huecos >=
    // 3072 y m >= 768) y sus bordes, desequilibradas y cuadradas.
    const std::vector<std::pair<std::size_t, std::size_t>> comunes = {
        {1, 1}, {1, 7}, {2, 2}, {3, 9}, {7, 7}, {7, 8}, {8, 8}, {8, 9}, {9, 30}, {12, 12}, {16, 200},
        {20, 20}, {25, 26}, {26, 26}, {31, 31}, {32, 32}, {40, 300}, {77, 77}, {78, 78}, {101, 101}};
    {
        const Fq64 R = f3e30();
        const FqSwar F(f3e30());
        auto formas = comunes;
        for (auto f : std::vector<std::pair<std::size_t, std::size_t>>{{1, 1000}, {4, 1000}, {26, 27}, {201, 201},
                                                                       {544, 545}, {545, 545}, {546, 546}, {1001, 1001}})
            formas.push_back(f);
        productos("F_(3^30) swar", R, F, formas);
    }
    productos("F_(5^10) swar", f5e10(), FqSwar(f5e10()), comunes);
    productos("F_(2^17) swar", f2e17(), FqSwar(f2e17()), comunes);
    {
        auto formas = comunes;
        for (auto f : std::vector<std::pair<std::size_t, std::size_t>>{{84, 84}, {85, 85}, {84, 600}, {85, 600}, {255, 255},
                                                                       {256, 256}, {300, 700}})
            formas.push_back(f);
        productos("F_(17^3) swar", f17e3(), FqSwar(f17e3()), formas);
    }
    {
        auto formas = comunes;
        for (auto f : std::vector<std::pair<std::size_t, std::size_t>>{{383, 383}, {384, 384}, {400, 900}})
            formas.push_back(f);
        productos("F_(127^2) swar", f127e2(), FqSwar(f127e2()), formas);
    }
    productos("F_(23^64) swar", f23e64(), FqSwar(f23e64()), {{1, 1}, {1, 5}, {2, 2}, {2, 3}, {8, 8}, {11, 40}, {12, 12}, {12, 40}, {30, 30}, {60, 61}});
    productos("F_(2^64) swar", f2e64(), FqSwar(f2e64()), {{1, 1}, {2, 9}, {8, 8}, {12, 12}, {13, 13}, {30, 30}, {60, 60}});
    productos("F_3 swar", f3(), FqSwar(f3()), comunes);
    productos("F_(3^16) swar", f3e16(), FqSwar(f3e16()), comunes);
    // FqTab: su producto no cambia en la fase 8, pero se comprueba igual.
    // Grandes: Karatsuba en paralelo desde 512 coeficientes con p impar y
    // desde 1024 con p = 2 (par_karatsuba_min).
    {
        auto formas = comunes;
        for (auto f : std::vector<std::pair<std::size_t, std::size_t>>{{511, 511}, {512, 512}, {600, 1100}, {1024, 1030}})
            formas.push_back(f);
        productos("F_256 tab", aes(), FqTab(aes()), formas);
        productos("F_(7^4) tab", f7e4(), FqTab(f7e4()), formas);
        productos("F_(3^5) tab", f3e5(), FqTab(f3e5()), formas);
    }

    // Batería en frío (d1, d2, grado de factor).
    frio<FqSwar>("F_(3^30) swar", f3e30, swar, {{9, 8, 5}, {40, 33, 5}, {300, 140, 8}});
    frio<FqSwar>("F_(5^10) swar", f5e10, swar, {{12, 9, 6}, {200, 120, 10}});
    frio<FqSwar>("F_(17^3) swar", f17e3, swar, {{100, 90, 12}, {300, 260, 7}});
    frio<FqSwar>("F_(127^2) swar", f127e2, swar, {{20, 9, 8}});
    frio<FqSwar>("F_(23^64) swar", f23e64, swar, {{14, 12, 3}});
    frio<FqSwar>("F_(3^16) swar", f3e16, swar, {{20, 9, 6}, {70, 40, 9}, {100, 60, 5}});
    frio<FqSwar>("F_(3^32) swar", f3e32, swar, {{40, 33, 4}});
    frio<FqTab>("F_256 tab", aes, tab, {{40, 33, 30}, {300, 140, 20}});
    frio<FqTab>("F_(2^12) tab", f2e12, tab, {{60, 40, 25}});
    frio<FqTab>("F_(3^7) tab", f3e7, tab, {{60, 40, 20}});
    frio<FqTab>("F_(7^4) tab", f7e4, tab, {{20, 9, 8}, {120, 70, 20}});
    frio<FqTab>("F_(97^2) tab", f97e2, tab, {{70, 40, 12}, {200, 100, 6}});

    // Ganchos nuevos de FqSwar: límites (con «requires» para que este mismo
    // programa compile también contra las cabeceras de 04c1909, sin ellos).
    auto limites = [](const auto& F, const auto& G) {  // genérica: el if constexpr descarta de verdad
        if constexpr (requires { F.kronecker16(std::size_t(0)); }) {
            CHECK(F.kronecker16(16354) && !F.kronecker16(16355) && F.palabras_bloque16() == 15);
            CHECK(F.bits_justos(1) == 3 && F.bits_justos(6030) == 15 && F.bits_justos(std::size_t(1) << 48) == 63);
            CHECK(G.kronecker16(3) && !G.kronecker16(4));
        }
    };
    limites(FqSwar(f3e30()), FqSwar(f127e2()));
    std::printf("huella: %016llx\n", static_cast<unsigned long long>(huella));
    return check::summary("test_fase8_fq");
}
