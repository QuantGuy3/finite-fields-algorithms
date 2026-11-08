// test_fase6_kronecker.cpp — Fase 6, propuesta 5 (q pequeño): mismas salidas
// que antes con
//
//   1. ks::mul_primo en Fp64 sin «u128 % p»: huecos de <= 64 bits leídos de una
//      palabra y reducidos con Barrett de 64 bits, y huecos de dos palabras
//      reducidos por trozos.  Se compara con el producto exacto (convolución en
//      BigInt) y con la ruta de Fp (GMP), que no cambia, en p = 2, 3, 101,
//      primos justo por debajo / por encima de 2^16, 2^30, 2^32, 2^40, 2^61 y
//      2^62, con coeficientes p - 1 (peor caso de la cota) y longitudes que
//      dejan el hueco justo en 64 / 65 y 128 / 129 bits.
//   2. ks::mul_extension en Fq64 con desempaquetado directo al acumulador:
//      frente a la ruta de Fq (GMP) en F_9, F_256, F_(3^30) y una extensión
//      cuadrática con p ~ 2^30 en la frontera de 64 bits.
//   3. Acumulador de 64 bits de Poly<Fp64> (mul_school, divmod_school y
//      Composer) y ganchos de Fp64 (acc_fma sin vigilar el acarreo con
//      p < 2^32, acc_reduce con Barrett de 64 bits): frente a Poly<Fp> y a la
//      reducción exacta en BigInt, en las fronteras de la cota
//      terminos·(p-1)^2 < 2^64.
//   4. Generador aleatorio en frío: cuerpos recién creados, misma semilla;
//      factor() y productos dan lo mismo en Fp64 que en Fp (y en Fq64 que en
//      Fq), el siguiente rand_below coincide y fft_cache_ queda igual de vacía.
//
// Los apartados 1, 2 y 4 solo usan la API pública anterior a la fase 6: el
// mismo fichero compilado con las cabeceras de antes también pasa.
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <cstdio>
#include <string>
#include <vector>

using namespace alcp;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

namespace {

BigInt potencia2(unsigned k) { return BigInt(1) << k; }

// Mayor primo < x y menor primo > x.
BigInt primo_bajo(BigInt x) {
    do x -= 1; while (!is_prime(x));
    return x;
}
BigInt primo_sobre(BigInt x) {
    do x += 1; while (!is_prime(x));
    return x;
}

std::size_t bitlen(const BigInt& x) { return sgn(x) == 0 ? 0 : mpz_sizeinbase(x.get_mpz_t(), 2); }

// Bits de hueco que usa ks::mul_primo (y mul_extension con m = min·e).
std::size_t bits_hueco(const BigInt& p, std::size_t m) {
    return 2 * bitlen(p - 1) + bitlen(BigInt(static_cast<unsigned long>(m))) + 1;
}

// Producto exacto, coeficiente a coeficiente, en BigInt y reducido módulo p.
std::vector<BigInt> conv_exacta(const std::vector<BigInt>& a, const std::vector<BigInt>& b, const BigInt& p) {
    std::vector<BigInt> r(a.size() + b.size() - 1, BigInt(0));
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j) r[i + j] += a[i] * b[j];
    for (auto& x : r) x %= p;
    return r;
}

std::vector<BigInt> a_big(const std::vector<u64>& v) {
    std::vector<BigInt> r;
    r.reserve(v.size());
    for (u64 x : v) r.push_back(detail::u64_to_big(x));
    return r;
}

// Vector de n coeficientes: todos p - 1 (peor caso) o aleatorios.
std::vector<u64> coefs(const Fp64& K, std::size_t n, bool maximos) {
    std::vector<u64> v(n);
    for (auto& x : v) x = maximos ? K.modulus64() - 1 : K.random();
    return v;
}

// ---------------------------------------------------------------- 1. mul_primo
bool mul_primo_ok(const Fp64& K, const Fp& KG, const std::vector<u64>& a, const std::vector<u64>& b) {
    const auto r = ks::mul_primo(K, a.data(), a.size(), b.data(), b.size());
    const auto ab = a_big(a), bb = a_big(b);
    const auto ref = conv_exacta(ab, bb, K.characteristic());
    const auto rg = ks::mul_primo(KG, ab.data(), ab.size(), bb.data(), bb.size());
    return a_big(r) == ref && rg == ref;
}

void test_mul_primo(const BigInt& p) {
    const Fp64 K(p);
    const Fp KG(p);
    const std::string nombre = "mul_primo Fp64 p = " + p.get_str();
    bool ok = true;
    const std::vector<std::pair<std::size_t, std::size_t>> tam = {
        {1, 1}, {1, 7}, {2, 2}, {5, 9}, {31, 33}, {64, 64}, {100, 37}, {1, 300}, {300, 300}};
    for (auto [na, nb] : tam)
        for (bool maximos : {true, false}) {
            const auto a = coefs(K, na, maximos), b = coefs(K, nb, maximos);
            ok = ok && mul_primo_ok(K, KG, a, b);
            // Cuadrado: el mismo puntero (rama mpn_sqr).
            const auto r = ks::mul_primo(K, a.data(), a.size(), a.data(), a.size());
            const auto ab = a_big(a);
            ok = ok && a_big(r) == conv_exacta(ab, ab, p);
        }
    // Longitudes que dejan el hueco justo en la frontera de 64 y de 128 bits:
    // m = 2^j - 1 (bitlen j) y m = 2^j (bitlen j + 1).
    for (std::size_t limite : {64u, 128u}) {
        const std::size_t base = 2 * bitlen(p - 1) + 1;
        if (limite <= base) continue;
        const std::size_t j = limite - base;  // bitlen(m) que da exactamente `limite` bits
        if (j >= 12) continue;                 // longitudes demasiado grandes para el test
        for (std::size_t m : {(std::size_t(1) << j) - 1, std::size_t(1) << j}) {
            if (m == 0) continue;
            const auto a = coefs(K, m, true), b = coefs(K, m + 3, true);
            ok = ok && mul_primo_ok(K, KG, a, b);
            ok = ok && (bits_hueco(p, m) == limite) == (m == (std::size_t(1) << j) - 1);
        }
    }
    // Longitud grande (Poly con Kronecker forzado) frente a Fp.
    {
        Poly<Fp64> A(K, coefs(K, 1000, true)), B(K, coefs(K, 999, false));
        Poly<Fp> AG(KG, a_big(A.coeffs())), BG(KG, a_big(B.coeffs()));
        ok = ok && a_big(A.mul(B, MulAlgo::Kronecker).coeffs()) == AG.mul(BG, MulAlgo::Kronecker).coeffs();
    }
    CHECK_MSG(ok, nombre.c_str());
}

// ---------------------------------------------------------------- 2. mul_extension
template <class KGt, class P64>
Poly<KGt> a_fq(const KGt& KG, const P64& A) {
    std::vector<typename KGt::Elem> v;
    for (const auto& x : A.coeffs()) v.push_back(KG.from_index(A.field().to_index(x)));
    return Poly<KGt>(KG, std::move(v));
}

void test_mul_extension(const BigInt& p, const std::vector<BigInt>& g, const char* nombre,
                        const std::vector<std::pair<int, int>>& tam) {
    const Fq64 K(Fp64(p), g);
    const Fq KG(Fp(p), g);
    using P = Poly<Fq64>;
    using PG = Poly<Fq>;
    bool ok = true;
    for (auto [na, nb] : tam)
        for (bool maximos : {true, false}) {
            std::vector<Fq64::Elem> va(na), vb(nb);
            const u64 pm1 = K.base_field().modulus64() - 1;
            for (auto& x : va) x = maximos ? Fq64::Elem(K.extension_degree(), pm1) : K.random();
            for (auto& x : vb) x = maximos ? Fq64::Elem(K.extension_degree(), pm1) : K.random();
            const auto r = ks::mul_extension(K, va.data(), va.size(), vb.data(), vb.size());
            const auto r2 = ks::mul_extension(K, va.data(), va.size(), va.data(), va.size());  // cuadrado
            P A(K, va), B(K, vb);
            const PG AG = a_fq(KG, A), BG = a_fq(KG, B);
            const auto& ca = AG.coeffs();
            const auto& cb = BG.coeffs();
            const auto rg = ks::mul_extension(KG, ca.data(), ca.size(), cb.data(), cb.size());
            const auto rg2 = ks::mul_extension(KG, ca.data(), ca.size(), ca.data(), ca.size());
            ok = ok && r.size() == rg.size() && r2.size() == rg2.size();
            for (std::size_t i = 0; ok && i < r.size(); ++i) ok = K.to_index(r[i]) == KG.to_index(rg[i]);
            for (std::size_t i = 0; ok && i < r2.size(); ++i) ok = K.to_index(r2[i]) == KG.to_index(rg2[i]);
            // Y el producto de Poly<Fq64> (algoritmo automático) frente al escolar de Fq.
            ok = ok && a_fq(KG, A * B) == AG.mul_school(BG);
        }
    CHECK_MSG(ok, nombre);
}

// Polinomio x^2 - c irreducible sobre F_p (c no residuo cuadrático).
std::vector<BigInt> cuadratica(const BigInt& p) {
    for (BigInt c = 2;; c += 1) {
        BigInt e = (p - 1) / 2, r;
        mpz_powm(r.get_mpz_t(), c.get_mpz_t(), e.get_mpz_t(), p.get_mpz_t());
        if (r == p - 1) return {p - c, BigInt(0), BigInt(1)};
    }
}

// ---------------------------------------------------------------- 3. acumulador de 64 bits
void test_acumulador(const BigInt& p, const std::vector<std::size_t>& longitudes) {
    const Fp64 K(p);
    const Fp KG(p);
    using P = Poly<Fp64>;
    using PG = Poly<Fp>;
    const std::string nombre = "acumulador Fp64 p = " + p.get_str();
    bool ok = true, okm = true, okd = true, okc = true;
    for (std::size_t n : longitudes)
        for (bool maximos : {true, false}) {
            // Producto escolar: el término central recibe min(na, nb) sumandos.
            P A(K, coefs(K, n, maximos)), B(K, coefs(K, n + 2, maximos)), C(K, coefs(K, n, maximos));
            const PG AG(KG, a_big(A.coeffs())), BG(KG, a_big(B.coeffs())), CG(KG, a_big(C.coeffs()));
            okm = okm && a_big(A.mul_school(B).coeffs()) == AG.mul_school(BG).coeffs();
            okm = okm && a_big(A.mul_school(C).coeffs()) == conv_exacta(AG.coeffs(), CG.coeffs(), p);
            // División escolar: nr - db + 1 sumandos por acumulador.
            // (coeficiente principal p - 1: grados exactos aunque sean aleatorios)
            auto vD = coefs(K, 2 * n + 1, maximos), vd = coefs(K, n / 2 + 2, maximos);
            vD.back() = vd.back() = K.modulus64() - 1;
            P D(K, vD), d(K, vd);
            const PG DG(KG, a_big(D.coeffs())), dG(KG, a_big(d.coeffs()));
            auto [q, r] = D.divmod_school(d);
            auto [qg, rg] = DG.divmod_school(dG);
            okd = okd && a_big(q.coeffs()) == qg.coeffs() && a_big(r.coeffs()) == rg.coeffs();
            okd = okd && q * d + r == D;
            // Divisor de grado 0 (sin resto) y dividendo más corto.
            auto [q0, r0] = D.divmod_school(P(K, std::vector<u64>{K.modulus64() - 1}));
            okd = okd && r0.is_zero() && q0 * P(K, std::vector<u64>{K.modulus64() - 1}) == D;
            auto [q1, r1] = d.divmod_school(D);
            okd = okd && q1.is_zero() && r1 == d;
        }
    // Composición de Brent-Kung (combina): cnt = k sumandos por coeficiente.
    for (std::size_t k : {std::size_t(1), std::size_t(5), std::size_t(17), std::size_t(40)}) {
        const int n = 40;
        P f = P::random(K, n);
        typename P::Reducer red(f);
        P h(K, coefs(K, n, true));
        typename P::Composer comp(h, red, k);
        const PG fG(KG, a_big(f.coeffs())), hG(KG, a_big(h.coeffs()));
        typename PG::Reducer redG(fG);
        typename PG::Composer compG(hG, redG, k);
        for (int dg : {0, 1, 17, 39, 120}) {
            P g(K, coefs(K, static_cast<std::size_t>(dg) + 1, dg % 2 == 1));
            const PG gG(KG, a_big(g.coeffs()));
            okc = okc && a_big(comp.compose(g).coeffs()) == compG.compose(gG).coeffs();
        }
    }
    // Ganchos de Fp64: (a) suma de muchos productos máximos, (b) valores de
    // 128 bits arbitrarios reducidos con acc_reduce, frente al resto exacto.
    {
        const u64 pm1 = K.modulus64() - 1;
        Fp64::Acc s;
        K.acc_clear(s);
        BigInt exacto(0);
        const BigInt pm1b = detail::u64_to_big(pm1);
        for (int i = 0; i < 300; ++i) {
            K.acc_fma(s, pm1, pm1);
            exacto += pm1b * pm1b;
            if (i % 37 == 0) {
                Fp64::Acc t = s;
                ok = ok && detail::u64_to_big(K.acc_reduce(t)) == exacto % p;
            }
        }
        ok = ok && detail::u64_to_big(K.acc_reduce(s)) == exacto % p;
        u64 x = 0x9e3779b97f4a7c15ULL;
        for (int i = 0; i < 2000; ++i) {
            x ^= x << 13, x ^= x >> 7, x ^= x << 17;
            const u64 hi = (i % 3 == 0) ? 0 : (i % 3 == 1 ? (x >> 40) : x * 31);
            const u64 lo = x * 0x2545F4914F6CDD1DULL;
            Fp64::Acc v = (static_cast<u128>(hi) << 64) | lo;
            const BigInt vb = (detail::u64_to_big(hi) << 64) + detail::u64_to_big(lo);
            ok = ok && detail::u64_to_big(K.acc_reduce(v)) == vb % p;
        }
        // Frontera exacta: 2^64 - 1, 2^64, y múltiplos de p alrededor.
        for (u128 v : {u128(~u64(0)), u128(1) << 64, (u128(1) << 64) - K.modulus64(),
                       u128(K.modulus64()) * K.modulus64(), u128(K.modulus64()) * K.modulus64() - 1}) {
            Fp64::Acc w = v;
            const BigInt vb = (detail::u64_to_big(static_cast<u64>(v >> 64)) << 64) + detail::u64_to_big(static_cast<u64>(v));
            ok = ok && detail::u64_to_big(K.acc_reduce(w)) == vb % p;
        }
    }
    CHECK_MSG(okm, (nombre + ": producto escolar").c_str());
    CHECK_MSG(okd, (nombre + ": división escolar").c_str());
    CHECK_MSG(okc, (nombre + ": composición").c_str());
    CHECK_MSG(ok, (nombre + ": ganchos acc_fma / acc_reduce").c_str());
}

// ---------------------------------------------------------------- 4. generador en frío
template <class PA, class PB>
bool mismo_factor(const PA& a, const PB& b) {
    auto ra = a.factor();
    auto rb = b.factor();
    if (ra.factors.size() != rb.factors.size()) return false;
    if (a.field().to_bigint(ra.leading) != b.field().to_bigint(rb.leading)) return false;
    for (std::size_t j = 0; j < ra.factors.size(); ++j) {
        const auto& fa = ra.factors[j].first.coeffs();
        const auto& fb = rb.factors[j].first.coeffs();
        if (fa.size() != fb.size() || ra.factors[j].second != rb.factors[j].second) return false;
        for (std::size_t i = 0; i < fa.size(); ++i)
            if (a.field().to_bigint(fa[i]) != b.field().to_bigint(fb[i])) return false;
    }
    return true;
}

// Cuerpos recién creados (caché FFT vacía), misma semilla: polinomios
// aleatorios, productos de varios tamaños (rutas sin generador y con él) y
// factor(); después, el siguiente número del generador y el estado de la caché.
void test_rng_primo(const BigInt& p, int grado) {
    bool ok = true;
    BigInt sig64, sigG;
    bool vacia64 = false, vaciaG = false;
    // Recorrido completo, cada cuerpo en frío con la misma semilla.
    auto recorrido = [&](auto& K, BigInt& sig, bool& vacia) {
        using F = std::remove_cvref_t<decltype(K)>;
        using P = Poly<F>;
        set_seed(4242);
        P f = P::random(K, grado), g = P::random(K, 2 * grado + 50), h = P::random(K, 300);
        std::vector<BigInt> huella;
        auto anota = [&](const P& r) {
            for (const auto& c : r.coeffs()) huella.push_back(K.to_bigint(c));
            huella.push_back(BigInt(-1));
        };
        anota(g % h);         // división y mcd: productos «sin generador» (mul_no_rng)
        anota(P::gcd(g, h));  //   antes de que la caché FFT tenga tablas
        anota(f * g);
        anota(g * h);
        anota(g % h);
        auto Fz = (f * f * g).factor();
        for (const auto& [fac, mult] : Fz.factors) {
            anota(fac);
            huella.push_back(BigInt(mult));
        }
        sig = rand_below(BigInt(1) << 64);
        vacia = K.fft_cache_.empty();
        return huella;
    };
    const Fp64 K(p);
    const Fp KG(p);
    const auto h64 = recorrido(K, sig64, vacia64);
    const auto hG = recorrido(KG, sigG, vaciaG);
    ok = ok && h64 == hG && sig64 == sigG && vacia64 == vaciaG;
    const std::string nombre = "generador en frío Fp64 = Fp, p = " + p.get_str();
    CHECK_MSG(ok, nombre.c_str());
}

void test_rng_extension(const BigInt& p, const std::vector<BigInt>& g, int grado, const char* nombre) {
    auto recorrido = [&](auto& K, BigInt& sig, bool& vacia) {
        using F = std::remove_cvref_t<decltype(K)>;
        using P = Poly<F>;
        set_seed(9191);
        P f = P::random(K, grado), h = P::random(K, 40);
        std::vector<BigInt> huella;
        auto anota = [&](const P& r) {
            for (const auto& c : r.coeffs()) huella.push_back(K.to_index(c));
            huella.push_back(BigInt(-1));
        };
        anota(h * h);
        anota(h % f);
        anota(P::gcd(h, f));
        auto Fz = f.factor();
        for (const auto& [fac, mult] : Fz.factors) {
            anota(fac);
            huella.push_back(BigInt(mult));
        }
        sig = rand_below(BigInt(1) << 64);
        vacia = K.fft_cache_.empty();
        return huella;
    };
    BigInt s1, s2;
    bool v1 = false, v2 = false;
    const Fq64 K(Fp64(p), g);
    const Fq KG(Fp(p), g);
    const auto h1 = recorrido(K, s1, v1);
    const auto h2 = recorrido(KG, s2, v2);
    CHECK_MSG(h1 == h2 && s1 == s2 && v1 == v2, nombre);
}

}  // namespace

int main() {
    std::printf("test_fase6_kronecker: mul_primo y mul_extension sin división, acumulador de 64 bits\n");
    set_seed(2026);

    // 1. mul_primo en Fp64.
    std::vector<BigInt> primos = {BigInt(2), BigInt(3), BigInt(101)};
    for (unsigned k : {16u, 28u, 30u, 31u, 32u, 40u, 60u, 61u, 62u}) primos.push_back(primo_bajo(potencia2(k)));
    for (unsigned k : {16u, 32u}) primos.push_back(primo_sobre(potencia2(k)));
    primos.push_back(BigInt(998244353));
    for (const auto& p : primos) test_mul_primo(p);

    // 2. mul_extension en Fq64.
    const std::vector<std::pair<int, int>> tam = {{1, 1}, {1, 5}, {2, 3}, {17, 17}, {40, 33}, {201, 201}};
    test_mul_extension(BigInt(3), {1, 0, 1}, "mul_extension F_9", tam);
    test_mul_extension(BigInt(2), {1, 1, 0, 1, 1, 0, 0, 0, 1}, "mul_extension F_256", tam);
    const std::vector<int> mod330 = {2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0,
                                     2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1};
    std::vector<BigInt> g330;
    for (int x : mod330) g330.emplace_back(x);
    test_mul_extension(BigInt(3), g330, "mul_extension F_(3^30)", {{1, 1}, {3, 2}, {20, 21}, {201, 201}});
    {
        // p ~ 2^30, e = 2: bits = 61 + bitlen(2·m): 64 bits hasta m = 3 y 65 desde m = 4.
        const BigInt p30 = primo_bajo(potencia2(30));
        test_mul_extension(p30, cuadratica(p30), "mul_extension F_(p^2), p ~ 2^30 (frontera de 64 bits)",
                           {{1, 1}, {2, 2}, {3, 3}, {3, 9}, {4, 4}, {5, 8}, {40, 40}});
        const BigInt p16 = primo_bajo(potencia2(16));
        test_mul_extension(p16, cuadratica(p16), "mul_extension F_(p^2), p ~ 2^16", tam);
    }

    // 3. Acumulador de 64 bits: fronteras de terminos·(p-1)^2 < 2^64.
    test_acumulador(BigInt(2), {1, 2, 40, 300});
    test_acumulador(BigInt(3), {1, 2, 40, 300});
    test_acumulador(BigInt(101), {1, 2, 40, 300});
    test_acumulador(primo_bajo(potencia2(16)), {1, 40, 300});
    test_acumulador(primo_bajo(potencia2(28)), {1, 200, 255, 256, 257, 300});  // ~256 términos
    test_acumulador(primo_bajo(potencia2(30)), {1, 14, 15, 16, 17, 18, 40});    // ~16 términos
    test_acumulador(primo_bajo(potencia2(32)), {1, 2, 3, 40});                  // 1 término
    test_acumulador(primo_sobre(potencia2(32)), {1, 2, 40});                    // nunca; sin «pequeño»
    test_acumulador(primo_bajo(potencia2(62)), {1, 2, 40});                     // acarreo de 128 bits

    // 4. Generador aleatorio en frío.
    test_rng_primo(BigInt(101), 40);
    test_rng_primo(BigInt(3), 30);
    test_rng_primo(BigInt(998244353), 30);  // cuerpo con NTT: la raíz sale del generador
    test_rng_primo(primo_bajo(potencia2(30)), 20);
    test_rng_extension(BigInt(3), {1, 0, 1}, 12, "generador en frío F_9: Fq64 = Fq");
    test_rng_extension(BigInt(2), {1, 1, 0, 1, 1, 0, 0, 0, 1}, 10, "generador en frío F_256: Fq64 = Fq");
    test_rng_extension(BigInt(3), g330, 4, "generador en frío F_(3^30): Fq64 = Fq");

    return check::summary("test_fase6_kronecker");
}
