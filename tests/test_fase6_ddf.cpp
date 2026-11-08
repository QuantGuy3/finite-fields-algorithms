// test_fase6_ddf.cpp — Fase 6, propuesta 6: la DDF de la fase 4 compone de
// inmediato tras cada factor encontrado (regla de alquilar o comprar con el
// gasto acumulado) en lugar de rehacer un pow_mod con el nuevo módulo.
//
//   1. Misma salida: distinct_degree_classic y distinct_degree frente a una
//      copia de la DDF clásica con el Frobenius anterior (el de la fase 5, que
//      pone la cuenta a cero en cada cambio de módulo), en F_2, F_3, F_9,
//      F_101, F_193, F_257, F_998244353 (Fp64 y Fp GMP) y F_(2^127-1), con
//      polinomios aleatorios, con muchos factores pequeños (muchos cambios de
//      módulo), con productos de lineales y con entradas dispersas.
//   2. Generador aleatorio en frío: sobre un cuerpo recién creado y con la
//      misma semilla, la versión nueva y la copia consumen el generador igual
//      (mismo rand_below posterior, misma caché FFT vacía o no), en cuerpos
//      con FFT (998244353, 193, 257) y sin ella (101, 2^127-1, F_9).
//   3. factor() sigue reconstruyendo la entrada con factores irreducibles.
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <optional>
#include <string>
#include <type_traits>
#include <vector>

using namespace alcp;

namespace {

// ------------------------------------------------- copia del algoritmo anterior
// Frobenius de la fase 5 (antes de la propuesta 6): cambiar_modulo pone
// pasos_ = 0 y hace falta al menos un pow_mod antes de volver a componer.
template <class Field>
class FrobeniusAnterior {
    using P = Poly<Field>;
    static constexpr std::size_t kMinProductos = 2;
    static constexpr std::size_t kPresupuesto =
        std::is_arithmetic_v<typename Field::Elem> ? (std::size_t{1} << 22) : (std::size_t{1} << 18);

public:
    explicit FrobeniusAnterior(const P& f) : red_(f), q_(f.field().cardinality()) {
        productos_q_ = mpz_sizeinbase(q_.get_mpz_t(), 2) + mpz_popcount(q_.get_mpz_t()) - 1;
        const std::size_t n = static_cast<std::size_t>(std::max(1, f.degree()));
        k_ = std::max<std::size_t>(1, std::min(n, kPresupuesto / n));
    }
    P siguiente(const P& H) {
        if (comp_) return comp_->compose(H);
        P r = P::pow_mod(H, q_, red_);
        ++pasos_;
        if (!h1_) h1_ = r;
        if (h1_ && quiere_componer()) comp_.emplace(*h1_, red_, k_);
        return r;
    }
    void cambiar_modulo(const P& g) {
        red_ = typename P::Reducer(g);
        const std::size_t n = static_cast<std::size_t>(std::max(1, g.degree()));
        k_ = std::max<std::size_t>(1, std::min(n, kPresupuesto / n));
        if (h1_) h1_ = *h1_ % g;
        comp_.reset();
        pasos_ = 0;
    }

private:
    typename P::Reducer red_;
    BigInt q_;
    std::size_t productos_q_ = 1, k_ = 1, pasos_ = 0;
    std::optional<P> h1_;
    std::optional<typename P::Composer> comp_;

    bool quiere_componer() const {
        if (productos_q_ < kMinProductos) return false;
        const std::size_t n = static_cast<std::size_t>(std::max(1, red_.modulus().degree()));
        const std::size_t gigantes = (n + k_ - 1) / k_;
        if (gigantes > productos_q_) return false;
        return pasos_ * productos_q_ >= k_ + 1;
    }
};

// distinct_degree_classic de la fase 5 (distinct_degree_impl(false)).
template <class Field>
std::vector<std::pair<int, Poly<Field>>> ddf_anterior(const Poly<Field>& f) {
    using P = Poly<Field>;
    const Field& K = f.field();
    std::vector<std::pair<int, P>> out;
    P g = f.monic();
    P x = P::variable(K);
    P H = x % g;
    FrobeniusAnterior<Field> frob(g);
    int d = 1;
    while (g.degree() >= 2 * d) {
        H = frob.siguiente(H);
        P hd = P::gcd(H - x, g);
        if (hd.degree() > 0) {
            out.emplace_back(d, hd.monic());
            g = g / hd;
            H = H % g;
            frob.cambiar_modulo(g);
        }
        ++d;
    }
    if (g.degree() > 0) out.emplace_back(g.degree(), g.monic());
    return out;
}

template <class Field>
bool iguales(const std::vector<std::pair<int, Poly<Field>>>& a, const std::vector<std::pair<int, Poly<Field>>>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i].first != b[i].first || !(a[i].second == b[i].second)) return false;
    return true;
}

// ------------------------------------------------------------- entradas
// Producto de polinomios aleatorios mónicos de grados 1, 2, ..., hasta grado
// total >= n, libre de cuadrados (muchos factores: muchos cambios de módulo).
template <class Field>
Poly<Field> muchos_factores(const Field& K, int n, int dmax) {
    using P = Poly<Field>;
    P f = P::one(K);
    int d = 1;
    while (f.degree() < n) {
        f = f * P::random(K, d).monic();
        d = d % dmax + 1;
    }
    return f.radical();
}

// x^m + c por un producto de lineales: H_i dispersos (productos pequeños en
// los primeros pow_mod) y cambios de módulo con módulos aún grandes.
template <class Field>
Poly<Field> dispersa(const Field& K, int m, int lineales) {
    using P = Poly<Field>;
    std::vector<typename Field::Elem> c(static_cast<std::size_t>(m) + 1, K.zero());
    c[0] = K.random();
    if (K.is_zero(c[0])) c[0] = K.one();
    c.back() = K.one();
    P f(K, std::move(c));
    for (int i = 0; i < lineales; ++i) f = f * P::random(K, 1).monic();
    return f.radical();
}

// --------------------------------------------------------- 1. misma salida
template <class Field>
void test_salida(const Field& K, const std::vector<int>& grados, int veces, const char* nombre) {
    using P = Poly<Field>;
    bool ok = true;
    auto prueba = [&](const P& f) {
        if (f.degree() < 1) return;
        const auto ref = ddf_anterior(f);
        ok = ok && iguales(f.distinct_degree_classic(), ref);
        ok = ok && iguales(f.distinct_degree(), ref);
    };
    // Casos límite: grado 1, 2 y 3.
    prueba(P::random(K, 1).monic());
    prueba(P::random(K, 2).monic().radical());
    prueba(P::random(K, 3).monic().radical());
    for (int n : grados)
        for (int t = 0; t < veces; ++t) {
            prueba(P::random(K, n).monic().radical());
            prueba(muchos_factores(K, n, 4));
            prueba(muchos_factores(K, n, 12));
        }
    // Producto de lineales distintos (un cambio de módulo en d = 1).
    prueba(muchos_factores(K, 20, 1));
    prueba(dispersa(K, 36, 3));
    prueba(dispersa(K, 20, 2));
    CHECK_MSG(ok, nombre);
}

// ------------------------------------------------ 2. generador en frío
// Crea el cuerpo con `crea` tras set_seed(semilla), construye la entrada con
// `entrada` y aplica la DDF nueva o la copia; devuelve el siguiente
// rand_below, si la caché FFT quedó vacía y la salida.
template <class Field, class Crea, class Entrada>
bool rng_igual(Crea crea, Entrada entrada, unsigned semilla) {
    set_seed(semilla);
    Field K1 = crea();
    Poly<Field> f1 = entrada(K1);
    const auto d1 = f1.distinct_degree_classic();
    const BigInt t1 = rand_below(BigInt(1000000007));
    const bool v1 = K1.fft_cache_.empty();

    set_seed(semilla);
    Field K2 = crea();
    Poly<Field> f2 = entrada(K2);
    const auto d2 = ddf_anterior(f2);
    const BigInt t2 = rand_below(BigInt(1000000007));
    const bool v2 = K2.fft_cache_.empty();

    // distinct_degree (que por debajo de los umbrales es el mismo camino).
    set_seed(semilla);
    Field K3 = crea();
    Poly<Field> f3 = entrada(K3);
    const auto d3 = f3.distinct_degree();
    const BigInt t3 = rand_below(BigInt(1000000007));
    const bool v3 = K3.fft_cache_.empty();
    const bool ok3 = f3.degree() >= 40 || (t3 == t2 && v3 == v2);  // >= 40: KS (no cambia)

    return t1 == t2 && v1 == v2 && ok3 && d1.size() == d2.size() && d3.size() == d2.size();
}

template <class Field, class Crea>
void test_rng(Crea crea, const char* nombre) {
    using P = Poly<Field>;
    bool ok = true;
    unsigned s = 7;
    for (int n : {10, 20, 32, 33, 34, 39, 45, 64}) {
        ok = ok && rng_igual<Field>(crea, [n](const Field& K) { return P::random(K, n).monic().radical(); }, ++s);
        ok = ok && rng_igual<Field>(crea, [n](const Field& K) { return muchos_factores(K, n, 3); }, ++s);
        ok = ok && rng_igual<Field>(crea, [n](const Field& K) { return dispersa(K, n, 2); }, ++s);
        ok = ok && rng_igual<Field>(crea, [n](const Field& K) { return dispersa(K, n - 4, 4); }, ++s);
    }
    CHECK_MSG(ok, nombre);
}

// ------------------------------------------------------------- 3. factor
template <class Field>
void test_factor(const Field& K, const char* nombre) {
    using P = Poly<Field>;
    bool ok = true;
    for (int n : {5, 10, 20, 39}) {
        P r = muchos_factores(K, n, 5);
        auto F = r.factor();
        P prod = P::one(K);
        for (auto& [u, e] : F.factors) {
            ok = ok && u.is_irreducible();
            for (int i = 0; i < e; ++i) prod = prod * u;
        }
        ok = ok && prod.scalar_mul(F.leading) == r;
    }
    CHECK_MSG(ok, nombre);
}

}  // namespace

int main() {
    set_seed(20261108);
    Fp64 w2(2), w3(3), w101(101), w193(193), w257(257), wntt(998244353);
    Fp g101(BigInt(101)), gntt(BigInt(998244353)), g127((BigInt(1) << 127) - 1);
    Fq64 w9(Fp64(BigInt(3)), std::vector<BigInt>{1, 0, 1});

    test_salida(w2, {10, 20, 39, 60}, 3, "DDF F_2");
    test_salida(w3, {10, 20, 39, 60}, 3, "DDF F_3");
    test_salida(w9, {10, 20, 39}, 3, "DDF F_9");
    test_salida(w101, {10, 20, 39, 45}, 3, "DDF F_101 (Fp64)");
    test_salida(g101, {10, 20, 39}, 2, "DDF F_101 (Fp GMP)");
    test_salida(w193, {10, 33, 39, 50}, 3, "DDF F_193");
    test_salida(w257, {10, 33, 39}, 3, "DDF F_257");
    test_salida(wntt, {10, 20, 39, 50}, 3, "DDF F_998244353 (Fp64)");
    test_salida(gntt, {10, 20, 39}, 2, "DDF F_998244353 (Fp GMP)");
    test_salida(g127, {10, 20, 39}, 2, "DDF F_(2^127-1)");

    test_rng<Fp64>([] { return Fp64(998244353); }, "RNG en frío F_998244353 (Fp64)");
    test_rng<Fp64>([] { return Fp64(193); }, "RNG en frío F_193 (Fp64)");
    test_rng<Fp64>([] { return Fp64(257); }, "RNG en frío F_257 (Fp64)");
    test_rng<Fp>([] { return Fp(BigInt(998244353)); }, "RNG en frío F_998244353 (Fp GMP)");
    test_rng<Fp>([] { return Fp(BigInt(193)); }, "RNG en frío F_193 (Fp GMP)");
    test_rng<Fp64>([] { return Fp64(101); }, "RNG en frío F_101 (sin FFT)");
    test_rng<Fp>([] { return Fp((BigInt(1) << 127) - 1); }, "RNG en frío F_(2^127-1) (sin FFT)");
    test_rng<Fq64>([] { return Fq64(Fp64(BigInt(3)), std::vector<BigInt>{1, 0, 1}); }, "RNG en frío F_9 (sin FFT)");

    test_factor(w2, "factor F_2");
    test_factor(w101, "factor F_101");
    test_factor(w9, "factor F_9");
    test_factor(g127, "factor F_(2^127-1)");
    return check::summary("fase6_ddf");
}
