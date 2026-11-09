// test_fase7_fpn.cpp — Fase 7, propuesta 3: cuerpo primo de dos palabras FpN2
// (include/alcp/fpn.hpp) frente a Fp (GMP).
//
// FpN2 debe ser INDISTINGUIBLE de Fp salvo en velocidad: mismos valores,
// mismas excepciones y mensajes, mismo consumo del generador aleatorio y
// mismo orden canónico.  Cada prueba hace la misma operación en los dos
// cuerpos y compara los resultados convertidos a BigInt.
//
//   1. Elementos: add/sub/neg/mul/inv/div/pow (exponentes negativos y grandes),
//      from_int (negativos, enormes, de 1 a 4 palabras), to_bigint/to_string,
//      con p de 63 y 64 bits (inv con la palabra alta nula), 65 bits, 100,
//      127 y 128 bits (hasta 2^128 - 159, el mayor primo de 128 bits) y primos
//      pequeños con NTT (el constructor admite cualquier primo impar < 2^128).
//   2. Límites: el menor primo > 2^j y el mayor < 2^(j+1) para j = 62..127, con
//      operandos extremos (todas las ramas de acarreo de la doble REDC).
//   3. Acumulador perezoso: contenidos arbitrarios de las 5 palabras (el
//      pliegue de la quinta da una o dos vueltas), cadenas de (p-1)^2 y
//      reduce_limbs (el desempaquetado de Kronecker) con 1 a 7 palabras.
//   4. Polinomios Poly<FpN2> frente a Poly<Fp>: mul (todos los algoritmos),
//      divmod (escolar y Newton), gcd/gcd_ext (clásico y half-GCD), pow_mod,
//      compose, inv_mod, factor/is_irreducible/distinct_degree/equal_degree/
//      radical con la misma semilla, y P::random / K.random.
//   5. FFT: p = 81·2^121 + 1 (128 bits, capacidad 121), 2^64 - 2^32 + 1
//      (64 bits, capacidad 32) y 998244353, 193, 257, 577: fft_base_root,
//      fft::tables, mul con FFT y el generador en frío (mismo rand_below
//      posterior y mismo estado de fft_cache_ que Fp).
//   6. Validación: p no primo (mismo mensaje que Fp), p = 2, p >= 2^128,
//      invertir 0, 0^k con k < 0, fits().
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <cstdint>
#include <cstdio>
#include <exception>
#include <string>
#include <typeinfo>
#include <vector>

using namespace alcp;

namespace {

using PG = Poly<Fp>;
using PF = Poly<FpN2>;

template <class Ex, class Fn>
bool throws_exactly(Fn&& f, std::string* msg = nullptr) {
    try {
        f();
    } catch (const std::exception& e) {
        if (msg) *msg = e.what();
        return typeid(e) == typeid(Ex);
    }
    return false;
}
#define CHECK_THROWS(Ex, expr) CHECK_MSG(throws_exactly<Ex>([&] { (void)(expr); }), "lanza " #Ex ": " #expr)

std::vector<BigInt> big(const PG& p) { return p.coeffs(); }
std::vector<BigInt> big(const PF& p) {
    std::vector<BigInt> r;
    for (const auto& c : p.coeffs()) r.push_back(p.field().to_bigint(c));
    return r;
}

// Generador propio (independiente del de la librería).
struct Rng {
    std::uint64_t s;
    std::uint64_t next() {
        std::uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    BigInt bits(unsigned n) {  // entero aleatorio de n bits como mucho
        BigInt r = 0;
        for (unsigned i = 0; i < n; i += 64) {
            const std::uint64_t x = next();
            BigInt w;
            mpz_import(w.get_mpz_t(), 1, -1, 8, 0, 0, &x);
            r = (r << 64) + w;
        }
        return r >> ((n + 63) / 64 * 64 - n);
    }
    BigInt below(const BigInt& n) { return bits(static_cast<unsigned>(bit_length(n)) + 64) % n; }
};

std::vector<BigInt> rand_ints(Rng& R, const BigInt& n, std::size_t len) {
    std::vector<BigInt> v;
    for (std::size_t i = 0; i < len; ++i) v.push_back(R.below(n));
    return v;
}

BigInt pot2(unsigned k) { return BigInt(1) << k; }

// Primos de prueba (todos los rangos y formas especiales).
std::vector<BigInt> primos() {
    return {
        pot2(62) + 135,                    // el menor primo de 63 bits
        pot2(63) - 25,                     // el mayor primo de 63 bits
        pot2(64) - pot2(32) + 1,           // 64 bits, capacidad FFT 32
        pot2(64) - 59,                     // el mayor primo de 64 bits
        pot2(64) + 13,                     // 65 bits
        BigInt(43) * pot2(94) + 1,         // 100 bits, capacidad 94
        pot2(127) - 1,                     // primo de Mersenne
        pot2(127) + 29,                    // el menor primo de 128 bits
        BigInt(81) * pot2(121) + 1,        // 128 bits, capacidad 121
        pot2(128) - 159,                   // el mayor primo de 128 bits
        BigInt(3), BigInt(193), BigInt(257), BigInt(577), BigInt(998244353),
    };
}

// ---------------------------------------------------------------------------
// 1. Elementos
// ---------------------------------------------------------------------------
void test_elementos(const BigInt& p) {
    Fp G(p);
    FpN2 F(p);
    Rng R{static_cast<std::uint64_t>(mpz_get_ui(p.get_mpz_t())) * 31 + 7};
    const std::string tag = " (p = " + p.get_str() + ")";
    // Operandos extremos y aleatorios.
    std::vector<BigInt> xs = {BigInt(0), BigInt(1), BigInt(2) % p, p - 2, p - 1, (p - 1) / 2, (p + 1) / 2};
    for (unsigned k : {32u, 63u, 64u, 65u, 96u, 127u})
        if (bit_length(p) > k) xs.push_back(pot2(k) % p), xs.push_back((pot2(k) - 1) % p);
    for (int i = 0; i < 300; ++i) xs.push_back(R.below(p));
    bool ok = true;
    for (std::size_t i = 0; i < xs.size() && ok; ++i) {
        const BigInt& x = xs[i];
        const BigInt& y = xs[(i * 7 + 3) % xs.size()];
        const auto a = F.from_int(x), b = F.from_int(y);
        ok = ok && F.to_bigint(a) == x && F.to_string(a) == G.to_string(x);
        ok = ok && F.to_bigint(F.add(a, b)) == G.add(x, y);
        ok = ok && F.to_bigint(F.sub(a, b)) == G.sub(x, y) && F.to_bigint(F.sub(b, a)) == G.sub(y, x);
        ok = ok && F.to_bigint(F.neg(a)) == G.neg(x);
        ok = ok && F.to_bigint(F.mul(a, b)) == G.mul(x, y) && F.to_bigint(F.mul(a, a)) == G.mul(x, x);
        ok = ok && F.eq(a, b) == G.eq(x, y) && (a < b) == (x < y);
        ok = ok && F.is_zero(a) == G.is_zero(x) && F.is_one(a) == G.is_one(x);
        if (y != 0) {
            ok = ok && F.to_bigint(F.inv(b)) == G.inv(y);
            ok = ok && F.to_bigint(F.div(a, b)) == G.div(x, y);
        }
        // exponentes con signo, de 1 y de varias palabras
        for (const BigInt& k : std::vector<BigInt>{BigInt(0), BigInt(1), BigInt(2), BigInt(-1), BigInt(-3), p - 2, p - 1, p, p * p + 12345,
                                R.bits(200) - pot2(199), -(p * p * 3)}) {
            if (x != 0 || k >= 0) ok = ok && F.to_bigint(F.pow(a, k)) == G.pow(x, k);
        }
        // from_int de enteros negativos, enormes y de 1 a 4 palabras
        const BigInt z = (R.below(p) - p * 3) * pot2(200) - 5;
        for (const BigInt& w : std::vector<BigInt>{z, BigInt(-z), x + p, x + p * 5, BigInt(-x), BigInt(-x - p * 7), R.bits(64), R.bits(128),
                                R.bits(129), R.bits(192), R.bits(256), BigInt(-R.bits(256)), R.bits(257),
                                pot2(128) - 1, pot2(256) - 1})
            ok = ok && F.to_bigint(F.from_int(w)) == G.from_int(w);
    }
    CHECK_MSG(ok, ("elementos de FpN2 = Fp" + tag).c_str());
    // from_i64 (conversión directa del binding) = from_int, también en los bordes.
    bool oki = true;
    for (long long v : {0LL, 1LL, -1LL, 7LL, -7LL, (1LL << 62), -(1LL << 62), (1LL << 62) + 135, -((1LL << 62) + 135),
                        9223372036854775807LL, -9223372036854775807LL - 1, -9223372036854775807LL,
                        static_cast<long long>(R.next()), static_cast<long long>(R.next())}) {
        BigInt b;
        mpz_set_si(b.get_mpz_t(), 0);
        const std::uint64_t mag = v >= 0 ? static_cast<std::uint64_t>(v) : 0ULL - static_cast<std::uint64_t>(v);
        mpz_import(b.get_mpz_t(), 1, -1, 8, 0, 0, &mag);
        if (v < 0) b = -b;
        oki = oki && F.from_i64(v) == F.from_int(b) && F.to_bigint(F.from_i64(v)) == G.from_int(b);
    }
    CHECK_MSG(oki, ("from_i64" + tag).c_str());
    CHECK(F.to_bigint(F.pow(F.zero(), BigInt(0))) == 1);
    CHECK(F.is_zero(F.pow(F.zero(), BigInt(5))));
    CHECK(F.is_one(F.one()) && F.is_zero(F.zero()) && F.eq(F.one(), F.one()));
    CHECK(F.cardinality() == p && F.characteristic() == p);
    CHECK_MSG(F.fft_capacity() == G.fft_capacity(), ("fft_capacity" + tag).c_str());
    const auto tres = F.from_int(BigInt(3));
    CHECK(F.pth_root(tres) == tres);
}

// ---------------------------------------------------------------------------
// 2. Límites de cada tamaño de p (ramas de acarreo de la doble REDC)
// ---------------------------------------------------------------------------
void test_limites() {
    bool ok = true;
    int n = 0;
    Rng R{777};
    for (unsigned j = 62; j <= 127 && ok; ++j) {
        BigInt lo = pot2(j), hi = pot2(j + 1);
        mpz_nextprime(lo.get_mpz_t(), lo.get_mpz_t());
        mpz_prevprime(hi.get_mpz_t(), hi.get_mpz_t());
        for (const BigInt& p : {lo, hi}) {
            FpN2 F(p);
            Fp G(p);
            std::vector<BigInt> xs = {BigInt(1), BigInt(2), p - 2, p - 1, (p + 1) / 2};
            for (int i = 0; i < 24; ++i) xs.push_back(R.below(p));
            for (const BigInt& x : xs) {
                for (const BigInt& y : xs) {
                    ok = ok && F.to_bigint(F.mul(F.from_int(x), F.from_int(y))) == G.mul(x, y);
                    ok = ok && F.to_bigint(F.add(F.from_int(x), F.from_int(y))) == G.add(x, y);
                }
                ok = ok && F.to_bigint(F.pow(F.from_int(x), p - 2)) == G.pow(x, p - 2);
                ok = ok && F.to_bigint(F.inv(F.from_int(x))) == G.inv(x);
            }
            ++n;
        }
    }
    CHECK_MSG(ok && n == 132, "doble REDC en los límites de cada tamaño de p (63-128 bits)");
}

// ---------------------------------------------------------------------------
// 3. Acumulador perezoso y reduce_limbs
// ---------------------------------------------------------------------------
BigInt valor_acc(const FpN2::Acc& s) {
    BigInt v = 0;
    for (int i = 4; i >= 0; --i) {
        BigInt w;
        mpz_import(w.get_mpz_t(), 1, -1, 8, 0, 0, &s.w[i]);
        v = (v << 64) + w;
    }
    return v;
}

void test_acumulador(const BigInt& p) {
    FpN2 F(p);
    Rng R{99 + static_cast<std::uint64_t>(mpz_get_ui(p.get_mpz_t()))};
    const std::string tag = " (p = " + p.get_str() + ")";
    bool ok = true;
    const std::uint64_t M = ~std::uint64_t(0);
    // Contenidos arbitrarios: incluido todo unos (dos vueltas de pliegue).
    std::vector<FpN2::Acc> casos = {FpN2::Acc{{M, M, M, M, M}}, FpN2::Acc{{M, M, M, M, 1}}, FpN2::Acc{{0, 0, 0, 0, M}},
                                    FpN2::Acc{{M, M, M, M, 0}}, FpN2::Acc{{0, 0, 0, M, M}}, FpN2::Acc{{1, 0, 0, 0, 0}},
                                    FpN2::Acc{}};
    for (int i = 0; i < 300; ++i) casos.push_back(FpN2::Acc{{R.next(), R.next(), R.next(), R.next(), R.next() >> (i % 64)}});
    for (auto s : casos) {
        const BigInt v = valor_acc(s);
        ok = ok && F.to_bigint(F.acc_reduce(s)) == v % p;
    }
    // Cadenas de productos máximos (p-1)^2 y aleatorios con acarreos.
    for (int len : {1, 2, 3, 17, 100, 1000}) {
        FpN2::Acc s;
        F.acc_clear(s);
        BigInt exact = 0;
        for (int i = 0; i < len; ++i) {
            const BigInt x = (i % 3 == 0) ? p - 1 : R.below(p), y = (i % 2 == 0) ? p - 1 : R.below(p);
            F.acc_fma(s, F.from_int(x), F.from_int(y));
            exact += x * y;
        }
        ok = ok && valor_acc(s) == exact && F.to_bigint(F.acc_reduce(s)) == exact % p;
    }
    FpN2::Acc s;
    F.acc_set(s, F.from_int(p - 1));
    F.acc_fma(s, F.from_int(p - 1), F.from_int(p - 1));
    ok = ok && F.to_bigint(F.acc_reduce(s)) == (p - 1 + (p - 1) * (p - 1)) % p;
    // reduce_limbs con 1 a 7 palabras (más de 5: camino de GMP).
    for (std::size_t n = 1; n <= 7; ++n)
        for (int i = 0; i < 20; ++i) {
            std::vector<mp_limb_t> t(n);
            for (auto& w : t) w = (i == 0) ? M : R.next();
            if (i == 1) t.back() = 0;
            BigInt v;
            mpz_import(v.get_mpz_t(), n, -1, 8, 0, 0, t.data());
            ok = ok && F.to_bigint(F.reduce_limbs(t.data(), n)) == v % p;
        }
    CHECK_MSG(ok, ("acumulador perezoso y reduce_limbs" + tag).c_str());
}

// ---------------------------------------------------------------------------
// 4. Polinomios
// ---------------------------------------------------------------------------
template <class A, class B>
bool mismos_factores(const A& fa, const B& fb, const FpN2& F) {
    bool ok = fa.factors.size() == fb.factors.size() && F.to_bigint(fb.leading) == fa.leading;
    for (std::size_t i = 0; ok && i < fa.factors.size(); ++i)
        ok = big(fb.factors[i].first) == big(fa.factors[i].first) && fb.factors[i].second == fa.factors[i].second;
    return ok;
}

void test_polinomios(const BigInt& p, unsigned long semilla, bool grande) {
    Fp G(p);
    FpN2 F(p);
    Rng R{semilla};
    const std::string tag = " (p = " + p.get_str() + ")";
    std::vector<int> grados = {0, 1, 3, 20, 31, 33, 70, 140};
    if (grande) grados.push_back(600);
    for (int d : grados) {
        const auto ai = rand_ints(R, p, static_cast<std::size_t>(d) + 1);
        const auto bi = rand_ints(R, p, static_cast<std::size_t>(d / 2) + 2);
        PG a = PG::from_ints(G, ai), b = PG::from_ints(G, bi);
        PF fa = PF::from_ints(F, ai), fb = PF::from_ints(F, bi);
        CHECK_MSG(big(fa) == big(a), ("from_ints" + tag).c_str());
        for (MulAlgo al : {MulAlgo::Auto, MulAlgo::Schoolbook, MulAlgo::Karatsuba, MulAlgo::FFT, MulAlgo::Kronecker})
            CHECK_MSG(big(fa.mul(fb, al)) == big(a.mul(b, al)), ("mul" + tag).c_str());
        CHECK_MSG(big(fa * fa) == big(a * a), ("cuadrado" + tag).c_str());
        CHECK_MSG(big(fa + fb) == big(a + b) && big(fa - fb) == big(a - b) && big(fb - fa) == big(b - a),
                  ("suma/resta" + tag).c_str());
        CHECK_MSG(big(-fa) == big(-a), ("opuesto" + tag).c_str());
        if (!b.is_zero()) {
            auto [q, r] = a.divmod(b);
            auto [fq, fr] = fa.divmod(fb);
            auto [nq, nr] = fa.divmod_newton(fb);
            CHECK_MSG(big(fq) == big(q) && big(fr) == big(r) && big(nq) == big(q) && big(nr) == big(r),
                      ("divmod" + tag).c_str());
        }
        CHECK_MSG(big(PF::gcd(fa, fb)) == big(PG::gcd(a, b)), ("gcd" + tag).c_str());
        auto e = PG::gcd_ext(a, b);
        auto fe = PF::gcd_ext(fa, fb);
        CHECK_MSG(big(fe.g) == big(e.g) && big(fe.x) == big(e.x) && big(fe.y) == big(e.y), ("gcd_ext" + tag).c_str());
        CHECK_MSG(big(fa.derivative()) == big(a.derivative()) && big(fa.monic()) == big(a.monic()),
                  ("derivative/monic" + tag).c_str());
        const BigInt x0 = R.below(p);
        CHECK_MSG(F.to_bigint(fa.evaluate(F.from_int(x0))) == a.evaluate(G.from_int(x0)), ("evaluate" + tag).c_str());
        if (d <= 33) CHECK_MSG(big(fb.compose(fb)) == big(b.compose(b)), ("compose" + tag).c_str());
        if (b.degree() >= 1 && d <= 140) {
            const BigInt k = p * p + 12345;
            CHECK_MSG(big(PF::pow_mod(fa, k, fb)) == big(PG::pow_mod(a, k, b)), ("pow_mod" + tag).c_str());
            CHECK_MSG(big(fb.pow(BigInt(3))) == big(b.pow(BigInt(3))), ("pow" + tag).c_str());
            const PG g = PG::gcd(a, b);
            if (g.degree() == 0) CHECK_MSG(big(PF::inv_mod(fa, fb)) == big(PG::inv_mod(a, b)), ("inv_mod" + tag).c_str());
        }
        CHECK_MSG(fa.to_string() == a.to_string(), ("to_string" + tag).c_str());
    }
    // Half-GCD y Newton (grados por encima de los umbrales de FpN2).
    if (grande) {
        const auto ai = rand_ints(R, p, 1300), bi = rand_ints(R, p, 1100);
        PG a = PG::from_ints(G, ai), b = PG::from_ints(G, bi);
        PF fa = PF::from_ints(F, ai), fb = PF::from_ints(F, bi);
        auto e = PG::gcd_ext_classic(a, b);
        auto fe = PF::gcd_ext(fa, fb);
        CHECK_MSG(big(fe.g) == big(e.g) && big(fe.x) == big(e.x) && big(fe.y) == big(e.y),
                  ("gcd_ext con half-GCD" + tag).c_str());
        const PG c = a * b;
        const PF fc = fa * fb;
        CHECK_MSG(big(fc) == big(c), ("producto de Kronecker 1300x1100" + tag).c_str());
        CHECK_MSG(big(PF::gcd(fc, fa * fa)) == big(PG::gcd(c, a * a)), ("gcd de grado ~2600" + tag).c_str());
        auto [q, r] = (c * c).divmod(a);
        auto [fq, fr] = (fc * fc).divmod(fa);
        CHECK_MSG(big(fq) == big(q) && big(fr) == big(r), ("divmod de Newton" + tag).c_str());
    }
    // Factorización e irreducibilidad, con la misma semilla en los dos cuerpos.
    for (int d : {1, 6, 13, 20}) {
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
        CHECK_MSG(mismos_factores(fac, ffac, F), ("factor (mismo orden canónico)" + tag).c_str());
        CHECK_MSG(ff.is_irreducible() == f.is_irreducible() && ff.is_squarefree() == f.is_squarefree(),
                  ("is_irreducible/is_squarefree" + tag).c_str());
        CHECK_MSG(big(ff.radical()) == big(f.radical()), ("radical" + tag).c_str());
        const PG rg = f.radical().monic();
        const PF rf = ff.radical().monic();
        const auto dg = rg.distinct_degree();
        const auto df = rf.distinct_degree();
        bool okd = dg.size() == df.size();
        for (std::size_t i = 0; okd && i < dg.size(); ++i) okd = dg[i].first == df[i].first && big(dg[i].second) == big(df[i].second);
        CHECK_MSG(okd, ("distinct_degree" + tag).c_str());
        for (std::size_t i = 0; okd && i < dg.size(); ++i) {
            set_seed(semilla + 5);
            const auto eg = dg[i].second.equal_degree(dg[i].first);
            set_seed(semilla + 5);
            const auto ef = df[i].second.equal_degree(df[i].first);
            bool oke = eg.size() == ef.size();
            for (std::size_t j = 0; oke && j < eg.size(); ++j) oke = big(eg[j]) == big(ef[j]);
            CHECK_MSG(oke, ("equal_degree" + tag).c_str());
        }
        for (const auto& [u, m] : fac.factors) {
            const PF fu = PF::from_ints(F, big(u));
            CHECK_MSG(ff.multiplicity(fu) == f.multiplicity(u) && fu.is_irreducible() == u.is_irreducible(),
                      ("multiplicity/is_irreducible de los factores" + tag).c_str());
        }
        // El generador queda en el mismo estado: el siguiente aleatorio coincide.
        set_seed(semilla + 2);
        (void)f.factor();
        const BigInt sig = rand_below(pot2(64));
        set_seed(semilla + 2);
        (void)ff.factor();
        CHECK_MSG(rand_below(pot2(64)) == sig, ("mismo consumo del generador en factor" + tag).c_str());
    }
    // P::random y K.random con la misma semilla.
    set_seed(semilla + 3);
    const PG rg = PG::random(G, 50);
    const BigInt rgx = G.random();
    const BigInt sig = rand_below(pot2(64));
    set_seed(semilla + 3);
    const PF rf = PF::random(F, 50);
    const BigInt rfx = F.to_bigint(F.random());
    CHECK_MSG(big(rf) == big(rg) && rfx == rgx && rand_below(pot2(64)) == sig, ("P::random y K.random" + tag).c_str());
}

// ---------------------------------------------------------------------------
// 5. FFT y generador en frío
// ---------------------------------------------------------------------------
void test_fft(const BigInt& p, unsigned long semilla) {
    const std::string tag = " (p = " + p.get_str() + ")";
    Rng R{semilla};
    // 5a. fft_base_root y tablas en cuerpos nuevos con la misma semilla.
    {
        set_seed(semilla);
        Fp G(p);
        const BigInt raiz = G.fft_base_root();
        const BigInt sig = rand_below(pot2(64));
        set_seed(semilla);
        FpN2 F(p);
        CHECK_MSG(F.to_bigint(F.fft_base_root()) == raiz, ("fft_base_root" + tag).c_str());
        CHECK_MSG(rand_below(pot2(64)) == sig, ("consumo del generador en fft_base_root" + tag).c_str());
        const unsigned long kmax = std::min<unsigned long>(G.fft_capacity(), 8);
        bool ok = true;
        for (unsigned long k = 0; k <= kmax; ++k) {
            const auto& tg = fft::tables(G, k);
            const auto& tf = fft::tables(F, k);
            ok = ok && tg.fwd.size() == tf.fwd.size() && G.to_bigint(tg.inv_n) == F.to_bigint(tf.inv_n);
            for (std::size_t i = 0; ok && i < tg.fwd.size(); ++i)
                ok = tg.fwd[i] == F.to_bigint(tf.fwd[i]) && tg.inv[i] == F.to_bigint(tf.inv[i]);
        }
        CHECK_MSG(ok, ("fft::tables" + tag).c_str());
        // dft / idft (fast.hpp) con la misma raíz.
        if (kmax >= 3) {
            const auto ai = rand_ints(R, p, 8);
            std::vector<FpN2::Elem> af;
            for (const auto& x : ai) af.push_back(F.from_int(x));
            const BigInt g = G.pow(raiz, pot2(static_cast<unsigned>(G.fft_capacity() - 3)));
            const auto dg = fast::dft(G, g, 3, ai);
            const auto df = fast::dft(F, F.from_int(g), 3, af);
            bool okd = dg.size() == df.size();
            for (std::size_t i = 0; okd && i < dg.size(); ++i) okd = dg[i] == F.to_bigint(df[i]);
            const auto vf = fast::idft(F, F.from_int(g), 3, df);
            for (std::size_t i = 0; okd && i < vf.size(); ++i) okd = F.to_bigint(vf[i]) == ai[i];
            CHECK_MSG(okd, ("fast::dft/idft" + tag).c_str());
        }
    }
    // 5b. Producto en frío con la regla Auto: L >= 64 pide la raíz (aunque
    // FpN2 no use la FFT), igual que Fp; L < 64 y el escolar no la piden.
    // Tamaños que en Fp cubren todas las ramas: sin raíz (20x20, L < 64;
    // 40x40 y 100x100, Kronecker), FFT usada (20x60 con p < 2^62; 20x200 con
    // la NTT multi-limb) y raíz pedida sin usarla (5x200 con la multi-limb).
    int construidas = 0;
    for (auto [na, nb] : std::vector<std::pair<std::size_t, std::size_t>>{{20, 20}, {40, 40}, {100, 100}, {20, 60}, {20, 200}, {5, 200}, {10, 150}}) {
        const auto ai = rand_ints(R, p, na), bi = rand_ints(R, p, nb);
        set_seed(semilla + 10);
        Fp G(p);
        const PG c = PG::from_ints(G, ai).mul(PG::from_ints(G, bi), MulAlgo::Schoolbook);
        const bool vacia_g0 = G.fft_cache_.empty();
        const PG d = PG::from_ints(G, ai) * PG::from_ints(G, bi);
        const bool vacia_g = G.fft_cache_.empty();
        const BigInt sig = rand_below(pot2(64));
        set_seed(semilla + 10);
        FpN2 F(p);
        const PF fc = PF::from_ints(F, ai).mul(PF::from_ints(F, bi), MulAlgo::Schoolbook);
        const bool vacia_f0 = F.fft_cache_.empty();
        const PF fd = PF::from_ints(F, ai) * PF::from_ints(F, bi);
        CHECK_MSG(big(fc) == big(c) && big(fd) == big(d), ("producto en frío" + tag).c_str());
        CHECK_MSG(vacia_g0 && vacia_f0 && vacia_g == F.fft_cache_.empty(), ("estado de fft_cache_ en frío" + tag).c_str());
        CHECK_MSG(rand_below(pot2(64)) == sig, ("consumo del generador en mul en frío" + tag).c_str());
        construidas += vacia_g ? 0 : 1;
        // Con la caché ya llena, mul con FFT explícita (si cabe) da lo mismo.
        CHECK_MSG(big(PF::from_ints(F, ai).mul(PF::from_ints(F, bi), MulAlgo::FFT)) ==
                      big(PG::from_ints(G, ai).mul(PG::from_ints(G, bi), MulAlgo::FFT)),
                  ("mul con FFT" + tag).c_str());
    }
    if (Fp(p).fft_capacity() >= 8) CHECK_MSG(construidas > 0, ("algún producto en frío construye las tablas" + tag).c_str());
    // 5c. Rutas «sin generador» en frío: divmod, gcd, gcd_ext, pow_mod no
    // construyen tablas que Fp no construya.
    {
        const auto ai = rand_ints(R, p, 300), bi = rand_ints(R, p, 140);
        set_seed(semilla + 20);
        Fp G(p);
        const PG a = PG::from_ints(G, ai), b = PG::from_ints(G, bi);
        auto [q, r] = a.divmod(b);
        const PG g = PG::gcd(a, b);
        auto e = PG::gcd_ext(a, b);
        const bool vacia_g = G.fft_cache_.empty();
        const BigInt sig = rand_below(pot2(64));
        set_seed(semilla + 20);
        FpN2 F(p);
        const PF fa = PF::from_ints(F, ai), fb = PF::from_ints(F, bi);
        auto [fq, fr] = fa.divmod(fb);
        const PF fg = PF::gcd(fa, fb);
        auto fe = PF::gcd_ext(fa, fb);
        const bool vacia_f = F.fft_cache_.empty();
        CHECK_MSG(big(fq) == big(q) && big(fr) == big(r) && big(fg) == big(g) && big(fe.g) == big(e.g) &&
                      big(fe.x) == big(e.x) && big(fe.y) == big(e.y),
                  ("divmod/gcd/gcd_ext en frío" + tag).c_str());
        CHECK_MSG(vacia_g == vacia_f, ("fft_cache_ tras divmod/gcd en frío" + tag).c_str());
        CHECK_MSG(rand_below(pot2(64)) == sig, ("consumo del generador en divmod/gcd en frío" + tag).c_str());
    }
}

// pow_mod en frío (consume el generador si algún producto tiene L >= 64).
void test_pow_mod_frio(const BigInt& p, unsigned long semilla) {
    const std::string tag = " (p = " + p.get_str() + ")";
    Rng R{semilla};
    for (std::size_t n : {10u, 40u, 70u}) {
        const auto ai = rand_ints(R, p, n), mi = rand_ints(R, p, n + 1);
        set_seed(semilla);
        Fp G(p);
        const PG r = PG::pow_mod(PG::from_ints(G, ai), p * 5 + 1, PG::from_ints(G, mi));
        const bool vg = G.fft_cache_.empty();
        const BigInt sig = rand_below(pot2(64));
        set_seed(semilla);
        FpN2 F(p);
        const PF fr = PF::pow_mod(PF::from_ints(F, ai), p * 5 + 1, PF::from_ints(F, mi));
        CHECK_MSG(big(fr) == big(r) && vg == F.fft_cache_.empty() && rand_below(pot2(64)) == sig,
                  ("pow_mod en frío (valor, caché y generador)" + tag).c_str());
    }
}

// ---------------------------------------------------------------------------
// 6. Validación
// ---------------------------------------------------------------------------
void test_validacion() {
    std::string mg, mf;
    const BigInt np = pot2(100) + 1;  // no primo (divisible por 17)
    CHECK(throws_exactly<std::invalid_argument>([&] { Fp G(np); }, &mg));
    CHECK(throws_exactly<std::invalid_argument>([&] { FpN2 F(np); }, &mf));
    CHECK_MSG(mg == mf, "mismo mensaje que Fp con p no primo");
    CHECK(throws_exactly<std::invalid_argument>([&] { FpN2 F(BigInt(1)); }, &mf));
    CHECK(throws_exactly<std::invalid_argument>([&] { Fp G(BigInt(1)); }, &mg) && mg == mf);
    CHECK(throws_exactly<std::invalid_argument>([&] { FpN2 F(BigInt(-7)); }, &mf));
    CHECK(throws_exactly<std::invalid_argument>([&] { Fp G(BigInt(-7)); }, &mg) && mg == mf);
    CHECK_THROWS(std::invalid_argument, FpN2(BigInt(2)));
    BigInt q = pot2(128);
    mpz_nextprime(q.get_mpz_t(), q.get_mpz_t());
    CHECK_THROWS(std::invalid_argument, FpN2(q));
    CHECK(!FpN2::fits(BigInt(0)) && !FpN2::fits(BigInt(-5)) && !FpN2::fits(pot2(62) - 1) && FpN2::fits(pot2(62)) &&
          FpN2::fits(pot2(128) - 1) && !FpN2::fits(pot2(128)));
    CHECK(!Fp64::fits(pot2(62)) && Fp64::fits(pot2(62) - 1));  // los rangos se tocan sin solaparse
    FpN2 F(pot2(127) - 1);
    Fp G(pot2(127) - 1);
    CHECK(throws_exactly<division_by_zero>([&] { (void)F.inv(F.zero()); }, &mf));
    CHECK(throws_exactly<division_by_zero>([&] { (void)G.inv(G.zero()); }, &mg) && mg == mf);
    CHECK_THROWS(division_by_zero, F.div(F.one(), F.zero()));
    CHECK(throws_exactly<division_by_zero>([&] { (void)F.pow(F.zero(), BigInt(-1)); }, &mf));
    CHECK(throws_exactly<division_by_zero>([&] { (void)G.pow(G.zero(), BigInt(-1)); }, &mg) && mg == mf);
    // Igualdad de cuerpos y polinomios de cuerpos distintos.
    FpN2 F2(pot2(127) - 1), F3(pot2(127) + 29);
    CHECK(F == F2 && !(F == F3) && same_field(F, F2) && !same_field(F, F3));
    const PF a = PF::variable(F), b = PF::variable(F2), c = PF::variable(F3);
    CHECK(a == b && !(a == c));
    CHECK_THROWS(std::invalid_argument, a + c);
    // Copias del cuerpo (el binding lo guarda en un std::variant).
    FpN2 copia = F;
    CHECK(copia == F && F.to_bigint(copia.mul(copia.from_int(BigInt(3)), copia.from_int(BigInt(5)))) == 15);
}

}  // namespace

int main() {
    for (const BigInt& p : primos()) {
        test_elementos(p);
        test_acumulador(p);
    }
    test_limites();
    unsigned long s = 1;
    for (const BigInt& p : primos()) {
        const bool grande = p == pot2(127) - 1 || p == pot2(128) - 159;
        test_polinomios(p, 100 + s, grande);
        ++s;
    }
    for (const BigInt& p : std::vector<BigInt>{BigInt(81) * pot2(121) + 1, pot2(64) - pot2(32) + 1, BigInt(43) * pot2(94) + 1,
                            BigInt(998244353), BigInt(193), BigInt(257), BigInt(577), pot2(127) - 1, pot2(128) - 159}) {
        test_fft(p, 500 + s);
        test_pow_mod_frio(p, 600 + s);
        ++s;
    }
    test_validacion();
    return check::summary("test_fase7_fpn");
}
