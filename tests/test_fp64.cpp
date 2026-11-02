// test_fp64.cpp — Backend de 64 bits (Fp64, Fq64 = FqT<Fp64>) frente a GMP.
//
// El backend de 64 bits debe ser INDISTINGUIBLE del de GMP salvo en velocidad:
// mismos valores, mismas excepciones, mismo consumo del generador aleatorio
// (random, fft_base_root, factor, Fq(q) aleatorio) y mismo orden canónico.
// Cada prueba hace la misma operación con los dos backends y compara los
// resultados convertidos a BigInt.
//
// Cubre:
//   1. Elementos de Fp64: add/sub/neg/mul/inv/div/pow (exponentes negativos y
//      grandes), from_int (negativos, enormes), to_bigint/to_string, con p
//      pequeños, medianos y cercanos a 2^62 (reducción de Barrett en el límite).
//   2. Acumulación perezosa con desbordamiento de 128 bits (pliegue por acarreo).
//   3. Polinomios: mul (school/karatsuba/fft/auto), divmod, gcd, gcd_ext,
//      pow_mod, is_irreducible, factor (con la misma semilla), derivative,
//      evaluate, compose.
//   4. FFT: misma raíz base y mismo consumo del generador.
//   5. Fq64 frente a Fq: aritmética, inverso, from_vector/from_index/to_index,
//      to_string, construcción aleatoria Fq64(q) == Fq(q) con la misma semilla,
//      factorización.
//   6. Excepciones y validación: p no primo, p demasiado grande, invertir 0.
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <utility>
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

// Coeficientes como BigInt, para comparar backends.
std::vector<BigInt> big(const Poly<Fp>& p) { return p.coeffs(); }
std::vector<BigInt> big(const Poly<Fp64>& p) {
    std::vector<BigInt> r;
    for (auto c : p.coeffs()) r.push_back(p.field().to_bigint(c));
    return r;
}
std::vector<BigInt> big(const Poly<Fq>& p) {
    std::vector<BigInt> r;
    for (const auto& c : p.coeffs()) r.push_back(p.field().to_index(c));
    return r;
}
std::vector<BigInt> big(const Poly<Fq64>& p) {
    std::vector<BigInt> r;
    for (const auto& c : p.coeffs()) r.push_back(p.field().to_index(c));
    return r;
}

// Polinomio aleatorio (generador propio, independiente del de la librería).
struct Rng {
    std::uint64_t s;
    std::uint64_t next() {
        std::uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    BigInt big_below(const BigInt& n) {
        BigInt r = (BigInt(static_cast<unsigned long>(next() >> 33)) << 62) + BigInt(static_cast<unsigned long>(next() >> 33)) * 7919 +
                   BigInt(static_cast<unsigned long>(next() & 0xFFFFFFFFu));
        r = (r << 64) + BigInt(static_cast<unsigned long>(next() & 0xFFFFFFFFu));
        return r % n;
    }
};

std::vector<BigInt> rand_ints(Rng& R, const BigInt& n, std::size_t len) {
    std::vector<BigInt> v;
    for (std::size_t i = 0; i < len; ++i) v.push_back(R.big_below(n));
    return v;
}

// ---------------------------------------------------------------------------
// 1-2. Elementos de Fp64
// ---------------------------------------------------------------------------
void test_elementos(const BigInt& p) {
    Fp G(p);
    Fp64 F(p);
    Rng R{std::uint64_t(mpz_get_ui(p.get_mpz_t())) * 31 + 7};
    bool ok = true;
    for (int it = 0; it < 2000 && ok; ++it) {
        const BigInt x = R.big_below(p), y = R.big_below(p);
        const auto a = F.from_int(x), b = F.from_int(y);
        ok = ok && F.to_bigint(a) == x && F.to_bigint(b) == y;
        ok = ok && F.to_bigint(F.add(a, b)) == G.add(x, y);
        ok = ok && F.to_bigint(F.sub(a, b)) == G.sub(x, y);
        ok = ok && F.to_bigint(F.neg(a)) == G.neg(x);
        ok = ok && F.to_bigint(F.mul(a, b)) == G.mul(x, y);
        if (y != 0) {
            ok = ok && F.to_bigint(F.inv(b)) == G.inv(y);
            ok = ok && F.to_bigint(F.div(a, b)) == G.div(x, y);
        }
        const BigInt k = R.big_below(p * p * 1000) - p * p * 500;  // exponentes con signo y grandes
        if (x != 0 || k >= 0) ok = ok && F.to_bigint(F.pow(a, k)) == G.pow(x, k);
        ok = ok && F.to_string(a) == G.to_string(x);
        // from_int de enteros negativos y enormes
        const BigInt z = (R.big_below(p) - p * 3) * (BigInt(1) << 200) - 5;
        ok = ok && F.to_bigint(F.from_int(z)) == G.from_int(z);
        ok = ok && F.to_bigint(F.from_int(-z)) == G.from_int(-z);
    }
    CHECK_MSG(ok, ("elementos de Fp64 = Fp con p = " + p.get_str()).c_str());
    CHECK(F.pow(F.zero(), BigInt(0)) == 1);
    CHECK(F.pow(F.zero(), BigInt(5)) == 0);
    CHECK(F.is_one(F.one()) && F.is_zero(F.zero()) && F.eq(F.one(), F.one()));
    CHECK(F.cardinality() == p && F.characteristic() == p);
    CHECK(F.fft_capacity() == G.fft_capacity());
    CHECK(F.pth_root(F.from_int(BigInt(3))) == F.from_int(BigInt(3)));

    // Acumulación perezosa: muchos productos máximos ((p-1)^2) seguidos, que
    // desbordan 128 bits cuando p es grande.
    Fp64::Acc s;
    F.acc_clear(s);
    BigInt exact = 0;
    const auto m = F.neg(F.one());  // p - 1
    for (int i = 0; i < 100; ++i) {
        F.acc_fma(s, m, m);
        exact += (p - 1) * (p - 1);
    }
    CHECK_MSG(F.to_bigint(F.acc_reduce(s)) == exact % p, "acumulación con desbordamiento");
    F.acc_set(s, F.from_int(BigInt(5)));
    CHECK(F.acc_reduce(s) == F.from_int(BigInt(5)));
}

// Producto (reducción de Barrett) y potencia en los límites de cada tamaño de
// p: el menor primo > 2^j y el mayor primo < 2^(j+1), j = 1..61 (todos los
// valores de k = bits de p), con los operandos extremos 1, 2, p-2, p-1 y
// pseudoaleatorios, frente a GMP.
void test_barrett_limites() {
    bool ok = true;
    int n = 0;
    Rng R{12345};
    for (unsigned j = 1; j <= 61 && ok; ++j) {
        BigInt lo = BigInt(1) << j, hi = BigInt(1) << (j + 1);
        mpz_nextprime(lo.get_mpz_t(), lo.get_mpz_t());
        mpz_prevprime(hi.get_mpz_t(), hi.get_mpz_t());
        for (const BigInt& p : {lo, hi}) {
            Fp64 F(p);
            Fp G(p);
            std::vector<BigInt> xs = {BigInt(1), BigInt(2) % p, p - 2, p - 1};
            for (int i = 0; i < 40; ++i) xs.push_back(R.big_below(p));
            for (const BigInt& x : xs)
                for (const BigInt& y : xs) ok = ok && F.to_bigint(F.mul(F.from_int(x), F.from_int(y))) == G.mul(x, y);
            for (const BigInt& x : xs) ok = ok && F.to_bigint(F.pow(F.from_int(x), p - 2)) == G.pow(x, p - 2);
            ++n;
        }
    }
    CHECK_MSG(ok && n == 122, "Barrett en los límites de cada tamaño de p");
}

// ---------------------------------------------------------------------------
// 3-4. Polinomios sobre Fp64 frente a Fp
// ---------------------------------------------------------------------------
void test_polinomios(const BigInt& p, unsigned long semilla) {
    Fp G(p);
    Fp64 F(p);
    Rng R{semilla};
    using PG = Poly<Fp>;
    using PF = Poly<Fp64>;
    const std::string tag = " (p = " + p.get_str() + ")";
    for (int d : {0, 3, 31, 33, 70, 150}) {
        const auto ai = rand_ints(R, p, static_cast<std::size_t>(d) + 1);
        const auto bi = rand_ints(R, p, static_cast<std::size_t>(d / 2) + 2);
        PG a = PG::from_ints(G, ai), b = PG::from_ints(G, bi);
        PF fa = PF::from_ints(F, ai), fb = PF::from_ints(F, bi);
        CHECK_MSG(big(fa) == big(a), ("from_ints" + tag).c_str());
        for (MulAlgo al : {MulAlgo::Auto, MulAlgo::Schoolbook, MulAlgo::Karatsuba, MulAlgo::FFT})
            CHECK_MSG(big(fa.mul(fb, al)) == big(a.mul(b, al)), ("mul" + tag).c_str());
        CHECK_MSG(big(fa + fb) == big(a + b) && big(fa - fb) == big(a - b) && big(fb - fa) == big(b - a),
                  ("suma/resta" + tag).c_str());
        CHECK_MSG(big(-fa) == big(-a), ("opuesto" + tag).c_str());
        if (!b.is_zero()) {
            auto [q, r] = a.divmod(b);
            auto [fq, fr] = fa.divmod(fb);
            CHECK_MSG(big(fq) == big(q) && big(fr) == big(r), ("divmod" + tag).c_str());
        }
        CHECK_MSG(big(PF::gcd(fa, fb)) == big(PG::gcd(a, b)), ("gcd" + tag).c_str());
        auto e = PG::gcd_ext(a, b);
        auto fe = PF::gcd_ext(fa, fb);
        CHECK_MSG(big(fe.g) == big(e.g) && big(fe.x) == big(e.x) && big(fe.y) == big(e.y), ("gcd_ext" + tag).c_str());
        CHECK_MSG(big(fa.derivative()) == big(a.derivative()), ("derivative" + tag).c_str());
        CHECK_MSG(big(fa.monic()) == big(a.monic()), ("monic" + tag).c_str());
        const BigInt x0 = R.big_below(p);
        CHECK_MSG(F.to_bigint(fa.evaluate(F.from_int(x0))) == a.evaluate(G.from_int(x0)), ("evaluate" + tag).c_str());
        if (d <= 33) CHECK_MSG(big(fb.compose(fb)) == big(b.compose(b)), ("compose" + tag).c_str());
        if (b.degree() >= 1) {
            const BigInt k = p * p + 12345;
            CHECK_MSG(big(PF::pow_mod(fa, k, fb)) == big(PG::pow_mod(a, k, b)), ("pow_mod" + tag).c_str());
        }
    }
    // Factorización e irreducibilidad, con la misma semilla en los dos backends.
    for (int d : {1, 6, 13}) {
        auto ci = rand_ints(R, p, static_cast<std::size_t>(d) + 1);
        ci.back() = 1;
        PG f = PG::from_ints(G, ci);
        f = f * f * PG::from_ints(G, rand_ints(R, p, 4));
        if (f.is_zero()) continue;
        PF ff = PF::from_ints(F, big(f));
        set_seed(semilla + 1);
        auto fac = f.factor();
        set_seed(semilla + 1);
        auto ffac = ff.factor();
        bool ok = fac.factors.size() == ffac.factors.size() && F.to_bigint(ffac.leading) == fac.leading;
        for (std::size_t i = 0; ok && i < fac.factors.size(); ++i)
            ok = big(ffac.factors[i].first) == big(fac.factors[i].first) &&
                 ffac.factors[i].second == fac.factors[i].second;
        CHECK_MSG(ok, ("factor (mismo orden canónico)" + tag).c_str());
        CHECK_MSG(ff.is_irreducible() == f.is_irreducible() && ff.is_squarefree() == f.is_squarefree(),
                  ("is_irreducible/is_squarefree" + tag).c_str());
        // El generador queda en el mismo estado: el siguiente aleatorio coincide.
        set_seed(semilla + 2);
        (void)f.factor();
        const BigInt sig = rand_below(BigInt(1) << 64);
        set_seed(semilla + 2);
        (void)ff.factor();
        CHECK_MSG(rand_below(BigInt(1) << 64) == sig, ("mismo consumo del generador en factor" + tag).c_str());
    }
    // random y fft_base_root consumen el generador igual que Fp.
    // (Cuerpos nuevos: G y F ya tienen la raíz en caché por las FFT de arriba.)
    set_seed(semilla + 3);
    Fp G2(p);
    const BigInt r1 = G2.random();
    const BigInt raiz = G2.fft_base_root();
    const BigInt sig = rand_below(BigInt(1) << 64);
    set_seed(semilla + 3);
    Fp64 F2(p);
    CHECK_MSG(F2.to_bigint(F2.random()) == r1, ("random" + tag).c_str());
    CHECK_MSG(F2.to_bigint(F2.fft_base_root()) == raiz, ("fft_base_root" + tag).c_str());
    CHECK_MSG(rand_below(BigInt(1) << 64) == sig, ("consumo del generador en fft_base_root" + tag).c_str());
}

// ---------------------------------------------------------------------------
// 5. Fq64 frente a Fq
// ---------------------------------------------------------------------------
void test_extension(const BigInt& p, const std::vector<BigInt>& g, unsigned long semilla) {
    Fq G(Fp(p), g);
    Fq64 F(Fp64(p), g);
    const std::string tag = " (F_" + p.get_str() + "^" + std::to_string(g.size() - 1) + ")";
    CHECK_MSG(F.modulus() == G.modulus() && F.cardinality() == G.cardinality() &&
                  F.extension_degree() == G.extension_degree() && F.characteristic() == G.characteristic(),
              ("datos del cuerpo" + tag).c_str());
    CHECK(F.fft_capacity() == G.fft_capacity());
    Rng R{semilla};
    const BigInt q = G.cardinality();
    bool ok = true;
    for (int it = 0; it < 60 && ok; ++it) {
        const BigInt x = R.big_below(q), y = R.big_below(q);
        const auto a = F.from_index(x), b = F.from_index(y);
        const auto ga = G.from_index(x), gb = G.from_index(y);
        ok = ok && F.to_index(a) == x;
        ok = ok && F.to_index(F.add(a, b)) == G.to_index(G.add(ga, gb));
        ok = ok && F.to_index(F.sub(a, b)) == G.to_index(G.sub(ga, gb));
        ok = ok && F.to_index(F.neg(a)) == G.to_index(G.neg(ga));
        ok = ok && F.to_index(F.mul(a, b)) == G.to_index(G.mul(ga, gb));
        ok = ok && F.to_string(a) == G.to_string(ga);
        if (y != 0) ok = ok && F.to_index(F.inv(b)) == G.to_index(G.inv(gb));
        ok = ok && F.to_index(F.pth_root(a)) == G.to_index(G.pth_root(ga));
        ok = ok && F.to_index(F.from_int(x)) == G.to_index(G.from_int(x));
        std::vector<BigInt> v = rand_ints(R, p * 5, g.size() + 3);
        v[0] -= p * 2;
        ok = ok && F.to_index(F.from_vector(v)) == G.to_index(G.from_vector(v));
    }
    CHECK_MSG(ok, ("elementos de Fq64 = Fq" + tag).c_str());
    CHECK_THROWS(division_by_zero, F.inv(F.zero()));
    CHECK_THROWS(division_by_zero, F.pow(F.zero(), BigInt(-1)));

    // Polinomios y factorización con la misma semilla.
    using PG = Poly<Fq>;
    using PF = Poly<Fq64>;
    std::vector<BigInt> ai = rand_ints(R, q, 40), bi = rand_ints(R, q, 23);
    auto conv = [](const auto& K, const std::vector<BigInt>& v) {
        std::vector<typename std::decay_t<decltype(K)>::Elem> e;
        for (const auto& x : v) e.push_back(K.from_index(x));
        return e;
    };
    PG a(G, conv(G, ai)), b(G, conv(G, bi));
    PF fa(F, conv(F, ai)), fb(F, conv(F, bi));
    for (MulAlgo al : {MulAlgo::Schoolbook, MulAlgo::Karatsuba, MulAlgo::FFT})
        CHECK_MSG(big(fa.mul(fb, al)) == big(a.mul(b, al)), ("mul" + tag).c_str());
    auto [q1, r1] = a.divmod(b);
    auto [q2, r2] = fa.divmod(fb);
    CHECK_MSG(big(q2) == big(q1) && big(r2) == big(r1), ("divmod" + tag).c_str());
    PG f = PG(G, conv(G, rand_ints(R, q, 4))) * PG(G, conv(G, rand_ints(R, q, 3)));
    f = f * f;
    if (f.is_zero()) f = PG::variable(G);
    PF ff(F, conv(F, big(f)));
    set_seed(semilla);
    auto fac = f.factor();
    set_seed(semilla);
    auto ffac = ff.factor();
    ok = fac.factors.size() == ffac.factors.size() && F.to_index(ffac.leading) == G.to_index(fac.leading);
    for (std::size_t i = 0; ok && i < fac.factors.size(); ++i)
        ok = big(ffac.factors[i].first) == big(fac.factors[i].first) && ffac.factors[i].second == fac.factors[i].second;
    CHECK_MSG(ok, ("factor" + tag).c_str());
}

// Construcción aleatoria: Fq64(q) y Fq64(p, n) eligen el MISMO módulo que Fq.
void test_construccion_aleatoria() {
    for (auto [p, n] : std::vector<std::pair<long, unsigned long>>{{2, 8}, {3, 5}, {7, 4}, {101, 3}}) {
        set_seed(4242 + static_cast<unsigned long>(p));
        Fq G(BigInt(p), n);
        set_seed(4242 + static_cast<unsigned long>(p));
        Fq64 F(BigInt(p), n);
        CHECK_MSG(F.modulus() == G.modulus(), ("Fq64(p, n) = Fq(p, n), p = " + std::to_string(p)).c_str());
        set_seed(99);
        Fq Gq(ipow(BigInt(p), n));
        set_seed(99);
        Fq64 Fq_(ipow(BigInt(p), n));
        CHECK_MSG(Fq_.modulus() == Gq.modulus(), ("Fq64(q) = Fq(q), p = " + std::to_string(p)).c_str());
        set_seed(7);
        const BigInt r = G.to_index(G.fft_base_root());
        set_seed(7);
        CHECK(F.to_index(F.fft_base_root()) == r);
    }
}

// ---------------------------------------------------------------------------
// 6. Validación
// ---------------------------------------------------------------------------
void test_validacion() {
    CHECK(Fp64::fits(BigInt(2)) && Fp64::fits((BigInt(1) << 62) - 1));
    CHECK(!Fp64::fits(BigInt(1) << 62) && !Fp64::fits(BigInt(0)) && !Fp64::fits(BigInt(-7)));
    CHECK_THROWS(std::invalid_argument, Fp64(BigInt(15)));
    CHECK_THROWS(std::invalid_argument, Fp64(BigInt(1)));
    // 2^89 - 1 es primo pero no cabe.
    CHECK_THROWS(std::invalid_argument, Fp64((BigInt(1) << 89) - 1));
    // q = p con p primo de 89 bits: la característica no cabe en Fp64.
    CHECK_THROWS(std::invalid_argument, Fq64((BigInt(1) << 89) - 1));
    Fp64 F(BigInt(101));
    CHECK_THROWS(division_by_zero, F.inv(0));
    CHECK_THROWS(division_by_zero, F.pow(0, BigInt(-3)));
    CHECK_THROWS(division_by_zero, F.div(1, 0));
    CHECK(same_field(F, Fp64(BigInt(101))) && !same_field(F, Fp64(BigInt(103))));
    const Fp64 F103(BigInt(103));
    CHECK_THROWS(std::invalid_argument, Poly<Fp64>::variable(F) + Poly<Fp64>::variable(F103));
    // Coeficiente principal nulo y módulo reducible (x^2 + 1 = (x + 2)(x + 3) en F_5).
    CHECK_THROWS(std::invalid_argument, Fq64(Fp64(BigInt(7)), std::vector<BigInt>{1, 0, 1, 0}));
    CHECK_THROWS(std::invalid_argument, Fq64(Fp64(BigInt(5)), std::vector<BigInt>{1, 0, 1}));
    CHECK_THROWS(std::invalid_argument, Fq64(BigInt(7), 0));
    Fq64 K(Fp64(BigInt(7)), std::vector<BigInt>{3, 4, 1, 2, 1});
    Fq64 K2(Fp64(BigInt(7)), std::vector<BigInt>{3, 4, 1, 2, 1});
    CHECK(same_field(K, K2) && K == K2);
}

}  // namespace

// Milisegundos desde t0 (informativo).
double ms_desde(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

int main() {
    set_seed(1);
    auto t0 = std::chrono::steady_clock::now();
    const BigInt m61 = (BigInt(1) << 61) - 1;
    // Mayor primo < 2^62: 2^62 - 57.
    const BigInt top = (BigInt(1) << 62) - 57;
    CHECK(is_prime(top) && Fp64::fits(top));
    for (const BigInt& p : {BigInt(2), BigInt(3), BigInt(7), BigInt(101), BigInt(998244353), BigInt("4294967311"),
                            m61, top})
        test_elementos(p);
    test_barrett_limites();
    std::printf("  elementos: %.0f ms\n", ms_desde(t0));
    unsigned long s = 10;
    for (const BigInt& p : {BigInt(2), BigInt(3), BigInt(7), BigInt(101), BigInt(998244353), m61, top})
        test_polinomios(p, s++);
    std::printf("  polinomios: %.0f ms\n", ms_desde(t0));
    test_extension(BigInt(2), {1, 1, 0, 1, 1, 0, 0, 0, 1}, 21);
    test_extension(BigInt(7), {3, 4, 1, 2, 1}, 22);
    test_extension(BigInt(3), {2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1},
                   23);
    test_extension(m61, {1, 0, 1}, 24);
    std::printf("  extensiones: %.0f ms\n", ms_desde(t0));
    test_construccion_aleatoria();
    test_validacion();
    std::printf("test_fp64: %.0f ms\n", ms_desde(t0));
    return check::summary("fp64");
}
