// test_fase7_fq.cpp — Fase 7, propuesta 7: los cuerpos de extensión rápidos
// FqTab (q <= 2^16, tablas de logaritmos) y FqSwar (p pequeño, huecos de 16
// bits) dan EXACTAMENTE lo mismo que Fq64 con el mismo módulo:
//
//   1. Elementos: random (misma secuencia y mismo estado del generador
//      después), from_int, from_index/to_index, from_vector, to_string, add,
//      sub, neg, mul, inv, div, pow (k negativo, cero, grande), pth_root,
//      is_one, el orden (operator< del Elem frente al lexicográfico de Fq64) y
//      la acumulación perezosa.  En los cuerpos pequeños, todos los pares.
//   2. Polinomios: productos (Auto, escolar, Karatsuba, Kronecker y FFT
//      cuando el cuerpo la admite), divmod, gcd, gcd_ext, pow_mod, compose,
//      radical, is_irreducible, distinct_degree, equal_degree y factor
//      (factores, multiplicidades, orden y coeficiente principal), con el
//      MISMO estado del generador después de cada cálculo.
//   3. Generador en frío: cada prueba crea los dos cuerpos desde cero con la
//      misma semilla y compara rand_below posterior y fft_cache_.empty(), en
//      cuerpos con FFT útil (F_(97^2): capacidad 6; F_(193^2): 7; F_(3^16) y
//      F_(3^32) en FqSwar: 6 y 7) y sin ella (F_256, F_(3^30)).
//   4. Construcción con módulo aleatorio (FqTab(q), FqSwar(p, n)): el mismo
//      módulo y el mismo generador que Fq64; mismas excepciones y mensajes.
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <algorithm>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

using namespace alcp;
using V = Fq64::Elem;

namespace {

// Elemento de K en la representación de Fq64.
V a_vec(const FqTab& K, FqTab::Elem e) { return K.decode(e); }
V a_vec(const FqSwar&, const V& e) { return e; }
FqTab::Elem de_vec(const FqTab& K, const V& v) { return K.encode(v); }
V de_vec(const FqSwar&, const V& v) { return v; }

template <class K>
std::vector<V> pvec(const K& F, const Poly<K>& p) {
    std::vector<V> r;
    for (const auto& c : p.coeffs()) r.push_back(a_vec(F, c));
    return r;
}
template <class K>
Poly<K> desde(const K& F, const Poly<Fq64>& p) {
    std::vector<typename K::Elem> c;
    for (const auto& x : p.coeffs()) c.push_back(de_vec(F, x));
    return Poly<K>(F, std::move(c));
}
template <class K>
bool igual(const K& F, const Poly<K>& a, const Poly<Fq64>& b) {
    return pvec(F, a) == b.coeffs();
}
BigInt testigo() { return rand_below(BigInt(1) << 64); }

std::vector<BigInt> big(const std::vector<long>& v) {
    std::vector<BigInt> r;
    for (long x : v) r.emplace_back(x);
    return r;
}

template <class Fn>
std::string mensaje(Fn&& fn) {
    try {
        fn();
    } catch (const std::exception& e) {
        return std::string(typeid(e).name()) + ": " + e.what();
    }
    return "(sin excepción)";
}

// ----------------------------------------------------------------- 1. elementos
template <class K>
void elementos(const char* nombre, const Fq64& R, const K& F, int pares) {
    bool ok = true;
    // random: misma secuencia y mismo generador después.
    set_seed(11);
    std::vector<V> ra;
    for (int i = 0; i < 64; ++i) ra.push_back(R.random());
    const BigInt ta = testigo();
    set_seed(11);
    std::vector<V> rb;
    for (int i = 0; i < 64; ++i) rb.push_back(a_vec(F, F.random()));
    CHECK_MSG(ra == rb && ta == testigo(), (std::string(nombre) + ": random").c_str());

    const BigInt q = R.cardinality();
    const bool todos = q * q <= BigInt(pares);
    std::vector<V> xs;
    if (todos) {
        for (BigInt i = 0; i < q; ++i) xs.push_back(R.from_index(i));
    } else {
        set_seed(12);
        xs.push_back(R.zero());
        xs.push_back(R.one());
        while (static_cast<int>(xs.size() * xs.size()) < pares) xs.push_back(R.random());
    }
    for (const V& x : xs) {
        const auto X = de_vec(F, x);
        ok = ok && a_vec(F, X) == x;
        ok = ok && F.to_index(X) == R.to_index(x) && a_vec(F, F.from_index(R.to_index(x))) == x;
        ok = ok && F.to_string(X) == R.to_string(x) && F.is_one(X) == R.is_one(x) && F.is_zero(X) == R.is_zero(x);
        ok = ok && a_vec(F, F.neg(X)) == R.neg(x) && a_vec(F, F.pth_root(X)) == R.pth_root(x);
        for (long k : {0L, 1L, 2L, 5L, -1L, -3L, 1000003L}) {
            if (R.is_zero(x) && k < 0) {
                ok = ok && mensaje([&] { (void)F.pow(X, BigInt(k)); }) == mensaje([&] { (void)R.pow(x, BigInt(k)); });
                continue;
            }
            ok = ok && a_vec(F, F.pow(X, BigInt(k))) == R.pow(x, BigInt(k));
        }
        ok = ok && a_vec(F, F.pow(X, q * q + 7)) == R.pow(x, q * q + 7);
        if (R.is_zero(x))
            ok = ok && mensaje([&] { (void)F.inv(X); }) == mensaje([&] { (void)R.inv(x); });
        else
            ok = ok && a_vec(F, F.inv(X)) == R.inv(x);
        for (const V& y : xs) {
            const auto Y = de_vec(F, y);
            ok = ok && a_vec(F, F.add(X, Y)) == R.add(x, y) && a_vec(F, F.sub(X, Y)) == R.sub(x, y);
            ok = ok && a_vec(F, F.mul(X, Y)) == R.mul(x, y) && (X < Y) == (x < y) && F.eq(X, Y) == R.eq(x, y);
            if (!R.is_zero(y)) ok = ok && a_vec(F, F.div(X, Y)) == R.div(x, y);
        }
        if (!ok) break;
    }
    CHECK_MSG(ok, (std::string(nombre) + ": operaciones de elementos").c_str());

    // from_int, from_vector, acumulación perezosa (con muchos términos para
    // cruzar la normalización de FqSwar).
    bool ok2 = true;
    for (long a : {0L, 1L, 2L, 3L, -1L, -7L, 123456789L})
        ok2 = ok2 && a_vec(F, F.from_int(BigInt(a))) == R.from_int(BigInt(a));
    for (const auto& v : {big({}), big({1}), big({-1, 5, 7}), big({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                                                    21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35})})
        ok2 = ok2 && a_vec(F, F.from_vector(v)) == R.from_vector(v);
    for (int rep = 0; rep < 8; ++rep) {
        typename K::Acc s;
        Fq64::Acc t;
        const V x0 = xs[static_cast<std::size_t>(rep) % xs.size()];
        F.acc_set(s, de_vec(F, x0));
        R.acc_set(t, x0);
        for (int i = 0; i < 300; ++i) {
            const V& a = xs[static_cast<std::size_t>(i * 7 + rep) % xs.size()];
            const V& b = xs[static_cast<std::size_t>(i * 13 + 3) % xs.size()];
            F.acc_fma(s, de_vec(F, a), de_vec(F, b));
            R.acc_fma(t, a, b);
        }
        ok2 = ok2 && a_vec(F, F.acc_reduce(s)) == R.acc_reduce(t);
    }
    CHECK_MSG(ok2, (std::string(nombre) + ": from_int, from_vector y acumulación").c_str());
    CHECK_MSG(F.fft_capacity() == R.fft_capacity() && F.cardinality() == R.cardinality() &&
                  F.characteristic() == R.characteristic() && F.extension_degree() == R.extension_degree() &&
                  F.modulus() == R.modulus(),
              (std::string(nombre) + ": datos del cuerpo").c_str());
}

// ------------------------------------------------------- 2-3. polinomios en frío
// Cada caso crea los dos cuerpos desde cero (raíz FFT sin calcular, caché
// vacía), fija la semilla, calcula y compara resultado, generador y caché.
template <class K>
void caso(const char* nombre, const std::function<Fq64()>& crear, const std::function<K(Fq64)>& envolver,
          const std::function<std::vector<std::vector<V>>(const Fq64&)>& calc_r,
          const std::function<std::vector<std::vector<V>>(const K&)>& calc_k, unsigned long semilla) {
    const Fq64 R = crear();
    const K F = envolver(crear());
    set_seed(semilla);
    const auto ra = calc_r(R);
    const BigInt ta = testigo();
    set_seed(semilla);
    const auto rb = calc_k(F);
    const BigInt tb = testigo();
    CHECK_MSG(ra == rb, (std::string(nombre) + ": resultado").c_str());
    if (ra != rb)
        for (std::size_t i = 0; i < std::min(ra.size(), rb.size()); ++i)
            if (ra[i] != rb[i]) {
                std::printf("    (primera diferencia en la posición %zu de %zu/%zu)\n", i, ra.size(), rb.size());
                break;
            }
    CHECK_MSG(ta == tb, (std::string(nombre) + ": generador después").c_str());
    CHECK_MSG(R.fft_cache_.empty() == F.fft_cache_.empty(), (std::string(nombre) + ": caché FFT").c_str());
}

// Resultados de un cálculo como lista de polinomios (en la representación Fq64).
template <class F>
struct Calc {
    using P = Poly<F>;
    static std::vector<V> vec(const F& K, const P& p) {
        std::vector<V> r;
        for (const auto& c : p.coeffs()) {
            if constexpr (std::is_same_v<F, Fq64>) r.push_back(c);
            else r.push_back(a_vec(K, c));
        }
        return r;
    }
    static std::vector<V> elem(const F& K, const typename F::Elem& e) {
        if constexpr (std::is_same_v<F, Fq64>) return {e};
        else return {a_vec(K, e)};
    }
    static std::vector<V> num(long x) { return {V(1, static_cast<std::uint64_t>(x))}; }
    static std::vector<std::vector<V>> factor(const F& K, const P& f) {
        auto r = f.factor();
        std::vector<std::vector<V>> out{elem(K, r.leading)};
        for (auto& [u, e] : r.factors) {
            out.push_back(vec(K, u));
            out.push_back(num(e));
        }
        return out;
    }
    // Batería sobre polinomios aleatorios de grados d1 >= d2 (todo sale del
    // generador, en el mismo orden en los dos cuerpos).
    static std::vector<std::vector<V>> bateria(const F& K, int d1, int d2) {
        std::vector<std::vector<V>> out;
        const P a = P::random(K, d1), b = P::random(K, d2);
        out.push_back(vec(K, a * b));
        out.push_back(vec(K, a.mul(b, MulAlgo::Schoolbook)));
        out.push_back(vec(K, a.mul(b, MulAlgo::Karatsuba)));
        out.push_back(vec(K, a.mul(b, MulAlgo::Kronecker)));
        auto [qq, rr] = a.divmod(b);
        out.push_back(vec(K, qq));
        out.push_back(vec(K, rr));
        out.push_back(vec(K, P::gcd(a, b)));
        auto e = P::gcd_ext(a, b);
        out.push_back(vec(K, e.g));
        out.push_back(vec(K, e.x));
        out.push_back(vec(K, e.y));
        out.push_back(vec(K, P::pow_mod(a, K.cardinality() + 3, b)));
        out.push_back(vec(K, a.compose(b) % P::random(K, d2 + 1)));
        out.push_back(num(b.is_irreducible() ? 1 : 0));
        return out;
    }
    // factor de un aleatorio, de un producto con multiplicidades y DDF/EDF.
    static std::vector<std::vector<V>> factores(const F& K, int d) {
        std::vector<std::vector<V>> out;
        const P f = P::random(K, d);
        for (auto& v : factor(K, f)) out.push_back(v);
        const P g = P::random(K, 2) * P::random(K, 2) * P::random(K, 1).pow(BigInt(3)) * P::random(K, 3);
        for (auto& v : factor(K, g)) out.push_back(v);
        out.push_back(vec(K, g.radical()));
        for (auto& [dd, prod] : f.radical().monic().distinct_degree()) {
            out.push_back(num(dd));
            for (auto& u : prod.equal_degree(dd)) out.push_back(vec(K, u));
        }
        return out;
    }
};

template <class K>
void polinomios(const char* nombre, const std::function<Fq64()>& crear, const std::function<K(Fq64)>& envolver,
                const std::vector<std::pair<int, int>>& grados, const std::vector<int>& grados_factor) {
    unsigned long s = 100;
    for (auto [d1, d2] : grados) {
        const std::string t = std::string(nombre) + " batería " + std::to_string(d1) + "/" + std::to_string(d2);
        caso<K>(t.c_str(), crear, envolver, [=](const Fq64& R) { return Calc<Fq64>::bateria(R, d1, d2); },
                [=](const K& F) { return Calc<K>::bateria(F, d1, d2); }, ++s);
    }
    for (int d : grados_factor) {
        const std::string t = std::string(nombre) + " factor " + std::to_string(d);
        caso<K>(t.c_str(), crear, envolver, [=](const Fq64& R) { return Calc<Fq64>::factores(R, d); },
                [=](const K& F) { return Calc<K>::factores(F, d); }, ++s);
    }
}

// FFT explícita (MulAlgo::FFT) y regla Auto en frío en un cuerpo con capacidad
// >= 6: la primera transformada calcula la raíz con el generador.
template <class K>
void prueba_fft(const char* nombre, const std::function<Fq64()>& crear, const std::function<K(Fq64)>& envolver) {
    for (int d : {20, 40, 70}) {
        const std::string t = std::string(nombre) + " FFT/Auto " + std::to_string(d);
        auto calc = [d](const auto& F) {
            using F_t = std::decay_t<decltype(F)>;
            using P = Poly<F_t>;
            std::vector<std::vector<V>> out;
            const P a = P::random(F, d), b = P::random(F, d);
            out.push_back(Calc<F_t>::vec(F, a * b));  // Auto: pide las tablas si L >= 64
            out.push_back(Calc<F_t>::vec(F, a.mul(b, MulAlgo::FFT)));
            out.push_back(Calc<F_t>::vec(F, a.divmod(P::random(F, d / 2 + 1)).second));
            out.push_back(Calc<F_t>::elem(F, F.fft_base_root()));
            return out;
        };
        caso<K>(t.c_str(), crear, envolver, calc, calc, 300 + static_cast<unsigned long>(d));
    }
}

// -------------------------------------------------------- 4. construcción
template <class K>
void construccion(const char* nombre, const BigInt& p, unsigned long n) {
    set_seed(77);
    const Fq64 R(p, n);
    const BigInt ta = testigo();
    set_seed(77);
    const K F(p, n);
    const BigInt tb = testigo();
    CHECK_MSG(R.modulus() == F.modulus() && ta == tb, (std::string(nombre) + ": módulo aleatorio (p, n)").c_str());
    const BigInt q = ipow(p, n);
    set_seed(78);
    const Fq64 R2(q);
    const BigInt tc = testigo();
    set_seed(78);
    const K F2(q);
    CHECK_MSG(R2.modulus() == F2.modulus() && tc == testigo(), (std::string(nombre) + ": módulo aleatorio (q)").c_str());
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

    const auto aes = mk(2, {1, 1, 0, 1, 1, 0, 0, 0, 1});
    const auto f4 = mk(2, {1, 1, 1});
    const auto f2 = mk(2, {1, 1});
    const auto f3 = mk(3, {1, 1});            // n = 1
    const auto f9 = mk(3, {1, 0, 1});
    const auto f7e4 = mk(7, {3, 4, 1, 2, 1});  // el de bench_suite
    const auto f97e2 = mk(97, {5, 0, 1});      // x^2 + 5 (-5 no es cuadrado mod 97): capacidad FFT 6
    const auto f193e2 = mk(193, {5, 0, 1});    // capacidad FFT 7
    const auto f3e30 = mk(3, {2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1});

    // ---------------------------------------------------------------- FqTab
    elementos("F_256 tab", aes(), FqTab(aes()), 1 << 17);
    elementos("F_4 tab", f4(), FqTab(f4()), 1 << 17);
    elementos("F_2 tab", f2(), FqTab(f2()), 1 << 17);
    elementos("F_3 tab", f3(), FqTab(f3()), 1 << 17);
    elementos("F_9 tab", f9(), FqTab(f9()), 1 << 17);
    elementos("F_(7^4) tab", f7e4(), FqTab(f7e4()), 1 << 14);
    elementos("F_(97^2) tab", f97e2(), FqTab(f97e2()), 1 << 14);
    elementos("F_(193^2) tab", f193e2(), FqTab(f193e2()), 1 << 14);
    {
        set_seed(5);
        const Fq64 R(BigInt(2), 16);
        elementos("F_(2^16) tab", R, FqTab(R), 1 << 12);
        set_seed(5);
        const Fq64 R3(BigInt(3), 10);
        elementos("F_(3^10) tab", R3, FqTab(R3), 1 << 12);
    }
    polinomios<FqTab>("F_256 tab", aes, tab, {{12, 5}, {40, 33}, {90, 70}, {300, 140}}, {5, 30, 45});
    polinomios<FqTab>("F_4 tab", f4, tab, {{20, 9}, {80, 40}}, {10, 40, 70});
    polinomios<FqTab>("F_2 tab", f2, tab, {{20, 9}, {100, 60}}, {20, 60});
    polinomios<FqTab>("F_9 tab", f9, tab, {{20, 9}, {80, 40}}, {12, 40});
    polinomios<FqTab>("F_(7^4) tab", f7e4, tab, {{20, 9}, {120, 70}}, {8, 20, 36});
    polinomios<FqTab>("F_(97^2) tab", f97e2, tab, {{20, 9}, {70, 40}, {200, 100}, {300, 200}}, {6, 20, 40, 70});
    polinomios<FqTab>("F_(193^2) tab", f193e2, tab, {{40, 33}, {150, 90}, {300, 200}}, {6, 34, 64});
    prueba_fft<FqTab>("F_(97^2) tab", f97e2, tab);
    prueba_fft<FqTab>("F_(193^2) tab", f193e2, tab);
    construccion<FqTab>("FqTab F_256", BigInt(2), 8);
    construccion<FqTab>("FqTab F_(5^3)", BigInt(5), 3);
    construccion<FqTab>("FqTab F_(2^16)", BigInt(2), 16);

    // --------------------------------------------------------------- FqSwar
    elementos("F_(3^30) swar", f3e30(), FqSwar(f3e30()), 1 << 14);
    elementos("F_(7^4) swar", f7e4(), FqSwar(f7e4()), 1 << 14);
    elementos("F_9 swar", f9(), FqSwar(f9()), 1 << 17);
    elementos("F_3 swar", f3(), FqSwar(f3()), 1 << 17);
    elementos("F_4 swar", f4(), FqSwar(f4()), 1 << 17);
    {
        set_seed(6);
        const Fq64 R2(BigInt(2), 64);
        elementos("F_(2^64) swar", R2, FqSwar(R2), 1 << 12);
        set_seed(6);
        const Fq64 R5(BigInt(5), 37);
        elementos("F_(5^37) swar", R5, FqSwar(R5), 1 << 12);
        set_seed(6);
        const Fq64 R17(BigInt(17), 3);  // 2·3·256 + 17 < 2^16
        elementos("F_(17^3) swar", R17, FqSwar(R17), 1 << 12);
    }
    polinomios<FqSwar>("F_(3^30) swar", f3e30, swar, {{12, 5}, {40, 33}, {200, 120}}, {3, 5, 8});
    polinomios<FqSwar>("F_(7^4) swar", f7e4, swar, {{20, 9}, {120, 70}}, {8, 20});
    {
        set_seed(9);
        const Fq64 R16(BigInt(3), 16);  // capacidad FFT 6
        set_seed(9);
        const Fq64 R32(BigInt(3), 32);  // capacidad FFT 7
        const auto c16 = [R16] { return Fq64(R16.base_field(), R16.modulus()); };
        const auto c32 = [R32] { return Fq64(R32.base_field(), R32.modulus()); };
        CHECK(R16.fft_capacity() == 6 && R32.fft_capacity() == 7);
        polinomios<FqSwar>("F_(3^16) swar", c16, swar, {{40, 33}, {100, 60}}, {4, 9});
        polinomios<FqSwar>("F_(3^32) swar", c32, swar, {{40, 33}}, {3});
        prueba_fft<FqSwar>("F_(3^16) swar", c16, swar);
        prueba_fft<FqSwar>("F_(3^32) swar", c32, swar);
    }
    construccion<FqSwar>("FqSwar F_(3^30)", BigInt(3), 30);
    construccion<FqSwar>("FqSwar F_(2^40)", BigInt(2), 40);

    // ------------------------------------------------- excepciones y límites
    {
        const auto m64 = mensaje([] { Fq64 k(Fp64(BigInt(2)), big({1, 0, 1})); });  // x^2 + 1 = (x + 1)^2
        CHECK(mensaje([] { FqTab k(Fp64(BigInt(2)), big({1, 0, 1})); }) == m64);
        CHECK(mensaje([] { FqSwar k(Fp64(BigInt(2)), big({1, 0, 1})); }) == m64);
        CHECK(mensaje([] { FqTab k(BigInt(12)); }) == mensaje([] { Fq64 k(BigInt(12)); }));
        CHECK(mensaje([] { FqTab k(BigInt(4), 2); }) == mensaje([] { Fq64 k(BigInt(4), 2); }));
        CHECK(mensaje([] { FqSwar k(BigInt(3), 0); }) == mensaje([] { Fq64 k(BigInt(3), 0); }));
        CHECK(FqTab::fits(BigInt(2), 16) && !FqTab::fits(BigInt(2), 17) && FqTab::fits(BigInt(251), 2) &&
              !FqTab::fits(BigInt(257), 2) && !FqTab::fits(BigInt(2), 0));
        CHECK(FqSwar::fits(BigInt(3), 30) && FqSwar::fits(BigInt(2), 64) && !FqSwar::fits(BigInt(2), 65) &&
              FqSwar::fits(BigInt(101), 2) && !FqSwar::fits(BigInt(101), 4) && FqSwar::fits(BigInt(23), 64) && !FqSwar::fits(BigInt(29), 64));
        set_seed(3);
        CHECK(mensaje([] { FqTab k(BigInt(2), 17); }).find("FqTab") != std::string::npos);
        CHECK(mensaje([] { FqSwar k(BigInt(101), 4); }).find("FqSwar") != std::string::npos);
        const FqTab T(aes());
        CHECK(T == FqTab(aes()) && same_field(T, FqTab(aes())) && !(T == FqTab(f4())));
        CHECK(T.to_coeffs(T.from_index(BigInt(0x53))) == aes().from_index(BigInt(0x53)));
    }
    return check::summary("test_fase7_fq");
}
