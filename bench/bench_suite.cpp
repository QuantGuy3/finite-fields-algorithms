// bench/bench_suite.cpp — Suite de benchmarks reproducible de alcp.
//
// Este archivo compila SIN CAMBIOS contra las cabeceras originales (etiqueta git
// `baseline-original`) y contra las cabeceras mejoradas del repositorio: solo usa
// API pública que ya existía en la versión original (Fp, Fq, Poly<F>, fast::*).
//
//   g++ -std=c++23 -O3 -Wall -Wextra -I<include> bench/bench_suite.cpp -lgmpxx -lgmp -o bench_suite
//   bench_suite --list
//   bench_suite --filter "fp\.ntt" --reps 7 --min-time 0.3 --json salida.json
//
// Metodología (ver bench/README.md):
//   * Entradas DETERMINISTAS e idénticas entre versiones: se generan con un PRNG
//     propio (splitmix64), independiente del RNG de la librería, de modo que un
//     cambio en el consumo de números aleatorios de alcp no altera las entradas.
//     Además se llama a set_seed(semilla fija del escenario) antes de preparar el
//     escenario y antes de cada repetición, para que la aleatoriedad interna de la
//     librería (Cantor–Zassenhaus, raíz de la FFT, Fq(p, n)) siga siempre el mismo
//     camino dentro de una misma versión.
//   * Los cuerpos F_q usan módulos irreducibles FIJOS escritos en el código.
//   * Las entradas se generan antes de cronometrar; solo se mide la operación.
//   * 1 calentamiento (configurable) + repeticiones hasta alcanzar a la vez un
//     mínimo de repeticiones (--reps) y un tiempo acumulado mínimo (--min-time).
//     Se guardan TODAS las muestras.
//   * Cada escenario calcula un checksum (FNV-1a de 64 bits) del resultado para
//     verificar que todas las variantes producen exactamente el mismo resultado.
//     En las factorizaciones el checksum es independiente del orden de los
//     factores (se ordenan los hashes de cada par (factor, multiplicidad)) y,
//     además, se reconstruye el producto: si no coincide con la entrada, el
//     checksum lleva el prefijo "INVALIDO-".
#include "alcp/alcp.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <functional>
#include <memory>
#include <optional>
#include <regex>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace alcp;

namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

// ===========================================================================
// Hash FNV-1a de 64 bits para los checksums.  Los enteros se vuelcan byte a
// byte en orden explícito, así que el resultado no depende de la representación
// interna de GMP ni de la versión de alcp.
// ===========================================================================
class Hasher {
public:
    void byte(unsigned char b) {
        h_ ^= b;
        h_ *= 1099511628211ULL;
    }
    void u64(std::uint64_t v) {
        for (int i = 0; i < 8; ++i) byte(static_cast<unsigned char>(v >> (8 * i)));
    }
    void i64(std::int64_t v) { u64(static_cast<std::uint64_t>(v)); }
    void big(const BigInt& b) {
        byte(static_cast<unsigned char>(sgn(b) + 1));
        const std::size_t nbytes = (mpz_sizeinbase(b.get_mpz_t(), 2) + 7) / 8;
        std::vector<unsigned char> buf(nbytes + 1);
        std::size_t count = 0;
        // Magnitud en big-endian, palabras de 1 byte.
        mpz_export(buf.data(), &count, 1, 1, 1, 0, b.get_mpz_t());
        u64(count);
        for (std::size_t i = 0; i < count; ++i) byte(buf[i]);
    }
    std::uint64_t value() const { return h_; }
    std::string hex() const {
        char s[24];
        std::snprintf(s, sizeof s, "%016llx", static_cast<unsigned long long>(h_));
        return s;
    }

private:
    std::uint64_t h_ = 1469598103934665603ULL;
};

std::uint64_t fnv_str(const std::string& s) {
    Hasher h;
    for (unsigned char c : s) h.byte(c);
    return h.value();
}

// ===========================================================================
// PRNG propio (splitmix64), independiente del RNG de alcp.  Genera las
// entradas de todos los escenarios.
// ===========================================================================
class Rng {
public:
    explicit Rng(std::uint64_t seed) : s_(seed) {}
    std::uint64_t next() {
        std::uint64_t z = (s_ += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    // Entero en [0, n), n > 0: bits(n) + 64 bits aleatorios reducidos módulo n
    // (sesgo despreciable; solo importa que sea determinista).
    BigInt below(const BigInt& n) {
        const std::size_t bits = mpz_sizeinbase(n.get_mpz_t(), 2) + 64;
        const std::size_t words = (bits + 63) / 64;
        std::vector<std::uint64_t> w(words);
        for (auto& x : w) x = next();
        BigInt r;
        mpz_import(r.get_mpz_t(), words, -1, sizeof(std::uint64_t), 0, 0, w.data());
        mpz_mod(r.get_mpz_t(), r.get_mpz_t(), n.get_mpz_t());
        return r;
    }

private:
    std::uint64_t s_;
};

// ===========================================================================
// Conversión elemento -> entero para hash y generación de elementos.  Se usan
// solo funciones de la API original: Fp::to_bigint/from_int, Fq::to_index/from_index.
// ===========================================================================
void hash_elem(Hasher& h, const Fp& K, const Fp::Elem& e) { h.big(K.to_bigint(e)); }
void hash_elem(Hasher& h, const Fq& K, const Fq::Elem& e) { h.big(K.to_index(e)); }

Fp::Elem rand_elem(const Fp& K, Rng& R) { return K.from_int(R.below(K.characteristic())); }
Fq::Elem rand_elem(const Fq& K, Rng& R) { return K.from_index(R.below(K.cardinality())); }

template <class F>
typename F::Elem rand_nonzero(const F& K, Rng& R) {
    while (true) {
        auto e = rand_elem(K, R);
        if (!K.is_zero(e)) return e;
    }
}

template <class F>
void hash_poly(Hasher& h, const Poly<F>& p) {
    h.i64(p.degree());
    for (const auto& c : p.coeffs()) hash_elem(h, p.field(), c);
}

template <class F>
void hash_vec(Hasher& h, const F& K, const std::vector<typename F::Elem>& v) {
    h.u64(v.size());
    for (const auto& c : v) hash_elem(h, K, c);
}

// Polinomio aleatorio de grado exacto d (mónico si se pide).
template <class F>
Poly<F> rand_poly(const F& K, Rng& R, int d, bool monic = false) {
    std::vector<typename F::Elem> v;
    v.reserve(static_cast<std::size_t>(d) + 1);
    for (int i = 0; i < d; ++i) v.push_back(rand_elem(K, R));
    v.push_back(monic ? K.one() : rand_nonzero(K, R));
    return Poly<F>(K, std::move(v));
}

// Polinomio mónico irreducible de grado d, por búsqueda (fuera del cronómetro).
// La búsqueda depende solo del PRNG propio y de que is_irreducible sea correcto.
template <class F>
Poly<F> rand_irreducible(const F& K, Rng& R, int d) {
    while (true) {
        Poly<F> f = rand_poly(K, R, d, true);
        if (f.is_irreducible()) return f;
    }
}

// ===========================================================================
// Checksums de resultados.
// ===========================================================================
template <class F>
std::string sum_poly(const Poly<F>& p) {
    Hasher h;
    hash_poly(h, p);
    return h.hex();
}

template <class F>
std::string sum_vec(const F& K, const std::vector<typename F::Elem>& v) {
    Hasher h;
    hash_vec(h, K, v);
    return h.hex();
}

template <class F>
std::string sum_elem(const F& K, const typename F::Elem& e) {
    Hasher h;
    hash_elem(h, K, e);
    return h.hex();
}

// Checksum de una factorización independiente del orden de los factores, con
// verificación de que leading * prod(u^e) reconstruye f.
template <class F, class Fac>
std::string sum_factor(const Poly<F>& f, const Fac& fac) {
    const F& K = f.field();
    std::vector<std::uint64_t> hs;
    Poly<F> prod = Poly<F>::constant(K, fac.leading);
    for (const auto& [u, e] : fac.factors) {
        Hasher h;
        hash_poly(h, u);
        h.i64(static_cast<std::int64_t>(e));
        hs.push_back(h.value());
        prod = prod * u.pow(BigInt(static_cast<long>(e)));
    }
    std::sort(hs.begin(), hs.end());
    Hasher t;
    t.u64(hs.size());
    for (auto x : hs) t.u64(x);
    hash_elem(t, K, fac.leading);
    return (prod == f ? std::string() : std::string("INVALIDO-")) + t.hex();
}

// ===========================================================================
// Registro de escenarios.
// ===========================================================================
struct Runner {
    std::function<void()> run;                 // operación cronometrada
    std::function<void()> reset;               // fuera del cronómetro: descarta el resultado previo
    std::function<std::string()> checksum;     // fuera del cronómetro, al final
};

using Maker = std::function<Runner(std::uint64_t seed)>;

struct Scenario {
    std::string name;   // identificador único (sin espacios), p. ej. "fp.ntt.mul_fft.d4096"
    std::string group;  // agrupación por cuerpo, p. ej. "fp.ntt"
    std::string op;     // agrupación por operación, p. ej. "mul_fft"
    std::string desc;   // descripción legible
    long inner;         // operaciones por repetición (informativo)
    Maker make;         // preparación: genera las entradas y devuelve el Runner
};

std::vector<Scenario>& registry() {
    static std::vector<Scenario> r;
    return r;
}

void add(std::string name, std::string group, std::string op, std::string desc, long inner, Maker make) {
    registry().push_back({std::move(name), std::move(group), std::move(op), std::move(desc), inner,
                          std::move(make)});
}

// Construye un Runner a partir de la operación `op` (devuelve un resultado R) y
// de la función de checksum `sum(const R&)`.  El resultado de cada repetición
// se guarda en `slot`; su destrucción ocurre en reset(), fuera del cronómetro.
template <class Op, class Sum>
Runner runner(Op op, Sum sum) {
    using R = std::decay_t<std::invoke_result_t<Op&>>;
    auto slot = std::make_shared<std::optional<R>>();
    auto opp = std::make_shared<Op>(std::move(op));
    Runner r;
    r.run = [slot, opp] { slot->emplace((*opp)()); };
    r.reset = [slot] { slot->reset(); };
    r.checksum = [slot, sum] {
        return slot->has_value() ? sum(**slot) : std::string("sin-resultado");
    };
    return r;
}

// Familia de escenarios sobre un cuerpo: prefijo del nombre, grupo, descripción
// del cuerpo y fábrica (cada escenario crea su propio cuerpo, con su propia
// caché de FFT, para que los escenarios sean independientes entre sí).
template <class F>
struct Fam {
    std::string prefix;
    std::string fdesc;
    std::function<std::shared_ptr<const F>()> mk;
};

std::string S(long v) { return std::to_string(v); }

// ---------------------------------------------------------------------------
// Micro-benchmarks de elemento.  Además del valor final se guarda una "foto"
// del acumulador cada kSnap iteraciones (coste despreciable), para que el
// checksum distinga resultados incluso en cuerpos con muy pocos elementos.
// ---------------------------------------------------------------------------
constexpr long kSnap = 1024;

template <class F>
void elem_muladd(const Fam<F>& fam, long n) {
    auto mk = fam.mk;
    add(fam.prefix + ".elem_muladd.n" + S(n), fam.prefix, "elem_muladd",
        S(n) + " iteraciones encadenadas x <- x*a + b en " + fam.fdesc, n,
        [mk, n](std::uint64_t seed) {
            auto K = mk();
            Rng R(seed);
            auto a = rand_nonzero(*K, R), b = rand_nonzero(*K, R), x0 = rand_nonzero(*K, R);
            return runner(
                [K, a, b, x0, n] {
                    std::vector<typename F::Elem> snap;
                    snap.reserve(static_cast<std::size_t>(n / kSnap) + 2);
                    auto x = x0;
                    for (long i = 0; i < n; ++i) {
                        x = K->add(K->mul(x, a), b);
                        if (i % kSnap == 0) snap.push_back(x);
                    }
                    snap.push_back(x);
                    return snap;
                },
                [K](const std::vector<typename F::Elem>& v) { return sum_vec(*K, v); });
        });
}

template <class F>
void elem_inv(const Fam<F>& fam, long n) {
    auto mk = fam.mk;
    add(fam.prefix + ".elem_inv.n" + S(n), fam.prefix, "elem_inv",
        S(n) + " inversos (acumulados con add) en " + fam.fdesc, n,
        [mk, n](std::uint64_t seed) {
            auto K = mk();
            Rng R(seed);
            std::vector<typename F::Elem> v;
            for (int i = 0; i < 256; ++i) v.push_back(rand_nonzero(*K, R));
            return runner(
                [K, v, n] {
                    std::vector<typename F::Elem> snap;
                    snap.reserve(static_cast<std::size_t>(n / kSnap) + 2);
                    auto acc = K->zero();
                    for (long i = 0; i < n; ++i) {
                        acc = K->add(acc, K->inv(v[static_cast<std::size_t>(i) & 255]));
                        if (i % kSnap == 0) snap.push_back(acc);
                    }
                    snap.push_back(acc);
                    return snap;
                },
                [K](const std::vector<typename F::Elem>& r) { return sum_vec(*K, r); });
        });
}

template <class F>
void elem_pow(const Fam<F>& fam, long n) {
    auto mk = fam.mk;
    add(fam.prefix + ".elem_pow.n" + S(n), fam.prefix, "elem_pow",
        S(n) + " potencias a^(q-2) (acumuladas con add) en " + fam.fdesc, n,
        [mk, n](std::uint64_t seed) {
            auto K = mk();
            Rng R(seed);
            std::vector<typename F::Elem> v;
            for (int i = 0; i < 256; ++i) v.push_back(rand_nonzero(*K, R));
            BigInt e = K->cardinality() - 2;
            return runner(
                [K, v, e, n] {
                    std::vector<typename F::Elem> snap;
                    snap.reserve(static_cast<std::size_t>(n / kSnap) + 2);
                    auto acc = K->zero();
                    for (long i = 0; i < n; ++i) {
                        acc = K->add(acc, K->pow(v[static_cast<std::size_t>(i) & 255], e));
                        if (i % kSnap == 0) snap.push_back(acc);
                    }
                    snap.push_back(acc);
                    return snap;
                },
                [K](const std::vector<typename F::Elem>& r) { return sum_vec(*K, r); });
        });
}

// ---------------------------------------------------------------------------
// Polinomios.
// ---------------------------------------------------------------------------
const char* algo_name(MulAlgo a) {
    switch (a) {
        case MulAlgo::Schoolbook: return "school";
        case MulAlgo::Karatsuba: return "karatsuba";
        case MulAlgo::FFT: return "fft";
        default: return "auto";
    }
}

// Producto de dos polinomios aleatorios de grado `deg`, `inner` veces por repetición.
// Con MulAlgo::Auto se usa operator* (despacho automático).
template <class F>
void poly_mul(const Fam<F>& fam, MulAlgo algo, int deg, long inner = 1) {
    auto mk = fam.mk;
    const std::string an = algo_name(algo);
    add(fam.prefix + ".mul_" + an + ".d" + S(deg), fam.prefix, "mul_" + an,
        "producto (" + an + ") de dos polinomios de grado " + S(deg) + " en " + fam.fdesc +
            (inner > 1 ? " x" + S(inner) : ""),
        inner, [mk, algo, deg, inner](std::uint64_t seed) {
            auto K = mk();
            Rng R(seed);
            Poly<F> a = rand_poly(*K, R, deg), b = rand_poly(*K, R, deg);
            return runner(
                [K, a, b, algo, inner] {
                    Poly<F> r = (algo == MulAlgo::Auto) ? a * b : a.mul(b, algo);
                    for (long i = 1; i < inner; ++i) r = (algo == MulAlgo::Auto) ? a * b : a.mul(b, algo);
                    return r;
                },
                [](const Poly<F>& r) { return sum_poly(r); });
        });
}

template <class F>
void poly_divmod(const Fam<F>& fam, int dnum, int dden) {
    auto mk = fam.mk;
    add(fam.prefix + ".divmod.d" + S(dnum) + "_" + S(dden), fam.prefix, "divmod",
        "divmod grado " + S(dnum) + " / grado " + S(dden) + " en " + fam.fdesc, 1,
        [mk, dnum, dden](std::uint64_t seed) {
            auto K = mk();
            Rng R(seed);
            Poly<F> a = rand_poly(*K, R, dnum), b = rand_poly(*K, R, dden);
            return runner([K, a, b] { return a.divmod(b); },
                          [](const std::pair<Poly<F>, Poly<F>>& qr) {
                              Hasher h;
                              hash_poly(h, qr.first);
                              hash_poly(h, qr.second);
                              return h.hex();
                          });
        });
}

// Pareja (a, b) de grado ~deg con un factor común de grado deg/8, para que el
// mcd no sea trivial.
template <class F>
std::pair<Poly<F>, Poly<F>> gcd_inputs(const F& K, Rng& R, int deg) {
    const int dc = std::max(1, deg / 8);
    Poly<F> c = rand_poly(K, R, dc, true);
    Poly<F> a = c * rand_poly(K, R, deg - dc);
    Poly<F> b = c * rand_poly(K, R, deg - dc - 1);
    return {a, b};
}

template <class F>
void poly_gcd(const Fam<F>& fam, int deg) {
    auto mk = fam.mk;
    add(fam.prefix + ".gcd.d" + S(deg), fam.prefix, "gcd",
        "mcd de dos polinomios de grado ~" + S(deg) + " con factor común en " + fam.fdesc, 1,
        [mk, deg](std::uint64_t seed) {
            auto K = mk();
            Rng R(seed);
            auto [a, b] = gcd_inputs(*K, R, deg);
            return runner([K, a, b] { return Poly<F>::gcd(a, b); },
                          [](const Poly<F>& g) { return sum_poly(g); });
        });
}

template <class F>
void poly_gcd_ext(const Fam<F>& fam, int deg) {
    auto mk = fam.mk;
    add(fam.prefix + ".gcd_ext.d" + S(deg), fam.prefix, "gcd_ext",
        "mcd extendido de dos polinomios de grado ~" + S(deg) + " en " + fam.fdesc, 1,
        [mk, deg](std::uint64_t seed) {
            auto K = mk();
            Rng R(seed);
            auto [a, b] = gcd_inputs(*K, R, deg);
            return runner([K, a, b] { return Poly<F>::gcd_ext(a, b); },
                          [](const auto& e) {
                              Hasher h;
                              hash_poly(h, e.g);
                              hash_poly(h, e.x);
                              hash_poly(h, e.y);
                              return h.hex();
                          });
        });
}

// x^(q^k) mod f, con f mónico aleatorio de grado deg.
template <class F>
void poly_pow_mod(const Fam<F>& fam, int deg, unsigned long k) {
    auto mk = fam.mk;
    add(fam.prefix + ".pow_mod.d" + S(deg) + "_k" + S(static_cast<long>(k)), fam.prefix, "pow_mod",
        "x^(q^" + S(static_cast<long>(k)) + ") mod f, f de grado " + S(deg) + " en " + fam.fdesc, 1,
        [mk, deg, k](std::uint64_t seed) {
            auto K = mk();
            Rng R(seed);
            Poly<F> f = rand_poly(*K, R, deg, true);
            Poly<F> x = Poly<F>::variable(*K);
            BigInt e = ipow(K->cardinality(), k);
            return runner([K, f, x, e] { return Poly<F>::pow_mod(x, e, f); },
                          [](const Poly<F>& r) { return sum_poly(r); });
        });
}

// Test de irreducibilidad sobre un polinomio irreducible (peor caso: test completo).
template <class F>
void poly_is_irreducible(const Fam<F>& fam, int deg) {
    auto mk = fam.mk;
    add(fam.prefix + ".is_irreducible.d" + S(deg), fam.prefix, "is_irreducible",
        "is_irreducible de un irreducible de grado " + S(deg) + " en " + fam.fdesc, 1,
        [mk, deg](std::uint64_t seed) {
            auto K = mk();
            Rng R(seed);
            Poly<F> f = rand_irreducible(*K, R, deg);  // búsqueda fuera del cronómetro
            // El checksum incluye la entrada: si una versión con un is_irreducible
            // defectuoso eligiera otro polinomio en la búsqueda, se detecta.
            return runner([K, f] { return f.is_irreducible(); },
                          [f](const bool& ok) {
                              Hasher h;
                              hash_poly(h, f);
                              h.byte(ok ? 1 : 0);
                              return h.hex();
                          });
        });
}

// Factorización de f = c * prod(u_i^{m_i}) con u_i mónicos aleatorios (no
// necesariamente irreducibles) de grados dados.  `shape` = {(grado, multiplicidad)}.
template <class F>
void poly_factor(const Fam<F>& fam, std::vector<std::pair<int, int>> shape, const std::string& extra = "") {
    auto mk = fam.mk;
    int deg = 0;
    for (auto [d, m] : shape) deg += d * m;
    add(fam.prefix + ".factor.d" + S(deg), fam.prefix, "factor",
        "factorización completa de un polinomio de grado " + S(deg) + " en " + fam.fdesc + extra, 1,
        [mk, shape](std::uint64_t seed) {
            auto K = mk();
            Rng R(seed);
            Poly<F> f = Poly<F>::constant(*K, rand_nonzero(*K, R));
            for (auto [d, m] : shape) {
                Poly<F> u = rand_poly(*K, R, d, true);
                for (int i = 0; i < m; ++i) f = f * u;
            }
            return runner([K, f] { return f.factor(); },
                          [f](const auto& fac) { return sum_factor(f, fac); });
        });
}

// ---------------------------------------------------------------------------
// Algoritmos rápidos (fast.hpp) sobre F_p con raíz primitiva conocida.
// ---------------------------------------------------------------------------
void fast_dft(const Fam<Fp>& fam, unsigned long gen, unsigned long k, bool inverse, long inner) {
    auto mk = fam.mk;
    const std::string op = inverse ? "idft" : "dft";
    add(fam.prefix + "." + op + ".k" + S(static_cast<long>(k)), "fast", op,
        op + " recursiva de tamaño 2^" + S(static_cast<long>(k)) + " en " + fam.fdesc +
            (inner > 1 ? " x" + S(inner) : ""),
        inner, [mk, gen, k, inverse, inner](std::uint64_t seed) {
            auto K = mk();
            Rng R(seed);
            // Raíz de orden exacto 2^k, FIJA (no aleatoria): gen^((p-1)/2^k).
            const BigInt cof = (K->characteristic() - 1) >> k;
            auto g = K->pow(K->from_int(BigInt(gen)), cof);
            std::vector<Fp::Elem> a;
            for (std::size_t i = 0; i < (std::size_t(1) << k); ++i) a.push_back(rand_elem(*K, R));
            return runner(
                [K, g, k, a, inverse, inner] {
                    auto r = inverse ? fast::idft(*K, g, k, a) : fast::dft(*K, g, k, a);
                    for (long i = 1; i < inner; ++i) r = inverse ? fast::idft(*K, g, k, a) : fast::dft(*K, g, k, a);
                    return r;
                },
                [K](const std::vector<Fp::Elem>& v) { return sum_vec(*K, v); });
        });
}

template <class F>
void fast_toeplitz_inverse(const Fam<F>& fam, std::size_t n) {
    auto mk = fam.mk;
    add(fam.prefix + ".toeplitz_lower_inverse.n" + S(static_cast<long>(n)), "fast", "toeplitz_lower_inverse",
        "inversa de Toeplitz triangular inferior n=" + S(static_cast<long>(n)) + " en " + fam.fdesc, 1,
        [mk, n](std::uint64_t seed) {
            auto K = mk();
            Rng R(seed);
            std::vector<typename F::Elem> a;
            a.push_back(rand_nonzero(*K, R));
            for (std::size_t i = 1; i < n; ++i) a.push_back(rand_elem(*K, R));
            return runner([K, n, a] { return fast::toeplitz_lower_inverse(*K, n, a); },
                          [K](const std::vector<typename F::Elem>& v) { return sum_vec(*K, v); });
        });
}

template <class F>
void fast_toeplitz_vec(const Fam<F>& fam, std::size_t n) {
    auto mk = fam.mk;
    add(fam.prefix + ".toeplitz_lower_vec.n" + S(static_cast<long>(n)), "fast", "toeplitz_lower_vec",
        "Toeplitz triangular inferior por vector n=" + S(static_cast<long>(n)) + " en " + fam.fdesc, 1,
        [mk, n](std::uint64_t seed) {
            auto K = mk();
            Rng R(seed);
            std::vector<typename F::Elem> a, b;
            for (std::size_t i = 0; i < n; ++i) a.push_back(rand_elem(*K, R));
            for (std::size_t i = 0; i < n; ++i) b.push_back(rand_elem(*K, R));
            return runner([K, n, a, b] { return fast::toeplitz_lower_vec(*K, n, a, b); },
                          [K](const std::vector<typename F::Elem>& v) { return sum_vec(*K, v); });
        });
}

template <class F>
void fast_toeplitz_divmod(const Fam<F>& fam, int dnum, int dden) {
    auto mk = fam.mk;
    add(fam.prefix + ".toeplitz_divmod.d" + S(dnum) + "_" + S(dden), "fast", "toeplitz_divmod",
        "división vía sistema de Toeplitz, grado " + S(dnum) + " / " + S(dden) + " en " + fam.fdesc, 1,
        [mk, dnum, dden](std::uint64_t seed) {
            auto K = mk();
            Rng R(seed);
            Poly<F> a = rand_poly(*K, R, dnum), b = rand_poly(*K, R, dden);
            return runner([K, a, b] { return fast::toeplitz_divmod(*K, a, b); },
                          [](const std::pair<Poly<F>, Poly<F>>& qr) {
                              Hasher h;
                              hash_poly(h, qr.first);
                              hash_poly(h, qr.second);
                              return h.hex();
                          });
        });
}

// ---------------------------------------------------------------------------
// Construcción de F_q con módulo aleatorio: Fq(p, n).  El checksum NO depende
// del módulo concreto elegido (otra versión puede consumir el RNG de otra forma
// y elegir otro irreducible igualmente válido): comprueba p, n y que el módulo
// sea mónico, de grado n e irreducible.
// ---------------------------------------------------------------------------
void fq_build(unsigned long p, unsigned long n) {
    const std::string tag = "q" + S(static_cast<long>(p)) + "e" + S(static_cast<long>(n));
    add("fq.build." + tag, "fq.build", "build",
        "construcción de F_" + S(static_cast<long>(p)) + "^" + S(static_cast<long>(n)) +
            " con módulo irreducible aleatorio (Fq(p, n))",
        1, [p, n](std::uint64_t) {
            return runner([p, n] { return std::make_shared<Fq>(BigInt(p), n); },
                          [](const std::shared_ptr<Fq>& F) {
                              const auto& g = F->modulus();
                              Poly<Fp> gp = Poly<Fp>::from_ints(F->base_field(), g);
                              const bool ok = gp.degree() == static_cast<int>(F->extension_degree()) &&
                                              gp.field().is_one(gp.lead()) && gp.is_irreducible();
                              Hasher h;
                              h.big(F->characteristic());
                              h.u64(F->extension_degree());
                              h.byte(ok ? 1 : 0);
                              return h.hex();
                          });
        });
}

// ===========================================================================
// Catálogo de escenarios.  Los tamaños están elegidos para que cada repetición
// tarde entre ~5 ms y ~1 s con la versión ORIGINAL compilada con -O3.
// ===========================================================================

// Primos de prueba (decimal).
constexpr const char* P_M61 = "2305843009213693951";                       // 2^61 - 1 (Mersenne, sin NTT)
constexpr const char* P_M127 = "170141183460469231731687303715884105727";  // 2^127 - 1 (Mersenne)
// Orden del subgrupo de BLS12-381 (255 bits, 2^32 | r - 1: admite NTT grandes).
constexpr const char* P_256 =
    "52435875175126190479447740508185965837690552500527637822603658699938581184513";

Fam<Fp> fp_family(const std::string& tag, const char* p, const std::string& fdesc) {
    const std::string ps = p;
    return {"fp." + tag, fdesc, [ps] { return std::make_shared<const Fp>(BigInt(ps)); }};
}

Fam<Fq> fq_family(const std::string& tag, const char* p, std::vector<long> modulus, const std::string& fdesc) {
    const std::string ps = p;
    return {"fq." + tag, fdesc, [ps, modulus] {
                std::vector<BigInt> g;
                for (long c : modulus) g.emplace_back(c);
                return std::make_shared<const Fq>(Fp(BigInt(ps)), std::move(g));
            }};
}

void register_all() {
    // ---------------- F_p, p = 7 ----------------
    {
        auto F = fp_family("p7", "7", "F_7");
        elem_muladd(F, 400000);
        elem_inv(F, 100000);
        elem_pow(F, 100000);
        poly_mul(F, MulAlgo::Schoolbook, 64, 20);
        poly_mul(F, MulAlgo::Schoolbook, 512);
        poly_mul(F, MulAlgo::Karatsuba, 512);
        poly_mul(F, MulAlgo::Karatsuba, 4096);
        poly_mul(F, MulAlgo::Auto, 2048);
        poly_divmod(F, 2000, 1000);
        poly_gcd(F, 500);
        poly_gcd_ext(F, 400);
        poly_pow_mod(F, 128, 48);   // 7^48 ~ 2^135
        poly_is_irreducible(F, 64);
        poly_factor(F, {{1, 2}, {2, 1}, {3, 2}, {5, 1}, {8, 1}, {12, 1}, {13, 1}, {14, 1}});  // grado 62
    }
    // ---------------- F_p, p = 101 ----------------
    {
        auto F = fp_family("p101", "101", "F_101");
        elem_muladd(F, 400000);
        elem_inv(F, 100000);
        elem_pow(F, 100000);
        poly_mul(F, MulAlgo::Schoolbook, 64, 20);
        poly_mul(F, MulAlgo::Schoolbook, 512);
        poly_mul(F, MulAlgo::Karatsuba, 512);
        poly_mul(F, MulAlgo::Karatsuba, 4096);
        poly_mul(F, MulAlgo::Auto, 2048);
        poly_divmod(F, 2000, 1000);
        poly_gcd(F, 500);
        poly_gcd_ext(F, 400);
        poly_pow_mod(F, 128, 20);   // 101^20 ~ 2^133
        poly_is_irreducible(F, 48);
        poly_factor(F, {{1, 3}, {2, 2}, {4, 1}, {6, 1}, {9, 1}, {14, 1}, {20, 1}, {40, 1}});  // grado 100
    }
    // ---------------- F_p, p = 998244353 (NTT) ----------------
    {
        auto F = fp_family("ntt", "998244353", "F_998244353");
        elem_muladd(F, 400000);
        elem_inv(F, 100000);
        elem_pow(F, 50000);
        poly_mul(F, MulAlgo::Schoolbook, 64, 20);
        poly_mul(F, MulAlgo::Schoolbook, 512);
        poly_mul(F, MulAlgo::Karatsuba, 512);
        poly_mul(F, MulAlgo::Karatsuba, 4096);
        poly_mul(F, MulAlgo::FFT, 4096);
        poly_mul(F, MulAlgo::FFT, 32768);
        poly_mul(F, MulAlgo::Auto, 2048);
        poly_divmod(F, 2000, 1000);
        poly_gcd(F, 500);
        poly_gcd_ext(F, 400);
        poly_pow_mod(F, 128, 5);    // p^5 ~ 2^150
        poly_is_irreducible(F, 32);
        poly_factor(F, {{1, 2}, {2, 1}, {3, 1}, {5, 2}, {8, 1}, {15, 1}});  // grado 40
    }
    // ---------------- F_p, p = 2^61 - 1 (sin NTT) ----------------
    {
        auto F = fp_family("m61", P_M61, "F_(2^61-1)");
        elem_muladd(F, 400000);
        elem_inv(F, 100000);
        elem_pow(F, 20000);
        poly_mul(F, MulAlgo::Schoolbook, 64, 20);
        poly_mul(F, MulAlgo::Schoolbook, 512);
        poly_mul(F, MulAlgo::Karatsuba, 512);
        poly_mul(F, MulAlgo::Karatsuba, 4096);
        poly_mul(F, MulAlgo::Auto, 2048);
        poly_divmod(F, 2000, 1000);
        poly_gcd(F, 500);
        poly_gcd_ext(F, 400);
        poly_pow_mod(F, 64, 2);     // p^2 ~ 2^122
        poly_is_irreducible(F, 24);
        poly_factor(F, {{1, 2}, {2, 1}, {4, 1}, {6, 1}, {12, 1}});  // grado 26
    }
    // ---------------- F_p, p = 2^127 - 1 ----------------
    {
        auto F = fp_family("m127", P_M127, "F_(2^127-1)");
        elem_muladd(F, 200000);
        elem_inv(F, 50000);
        elem_pow(F, 10000);
        poly_mul(F, MulAlgo::Schoolbook, 64, 20);
        poly_mul(F, MulAlgo::Schoolbook, 512);
        poly_mul(F, MulAlgo::Karatsuba, 512);
        poly_mul(F, MulAlgo::Karatsuba, 2048);
        poly_mul(F, MulAlgo::Auto, 2048);
        poly_divmod(F, 2000, 1000);
        poly_gcd(F, 400);
        poly_gcd_ext(F, 300);
        poly_pow_mod(F, 64, 1);     // p ~ 2^127
        poly_is_irreducible(F, 16);
        poly_factor(F, {{1, 2}, {2, 1}, {3, 1}, {5, 1}, {7, 1}});  // grado 20
    }
    // ---------------- F_p, p de 255 bits (BLS12-381 r, NTT) ----------------
    {
        auto F = fp_family("p256", P_256, "F_r (r de 255 bits, BLS12-381)");
        elem_muladd(F, 200000);
        elem_inv(F, 50000);
        elem_pow(F, 2000);
        poly_mul(F, MulAlgo::Schoolbook, 64, 20);
        poly_mul(F, MulAlgo::Schoolbook, 512);
        poly_mul(F, MulAlgo::Karatsuba, 512);
        poly_mul(F, MulAlgo::Karatsuba, 2048);
        poly_mul(F, MulAlgo::FFT, 4096);
        poly_mul(F, MulAlgo::FFT, 8192);
        poly_mul(F, MulAlgo::Auto, 2048);
        poly_divmod(F, 1000, 500);
        poly_gcd(F, 300);
        poly_gcd_ext(F, 200);
        poly_pow_mod(F, 64, 1);     // p ~ 2^255
        poly_is_irreducible(F, 12);
        poly_factor(F, {{1, 2}, {2, 1}, {3, 1}, {5, 1}});  // grado 12
    }

    // ---------------- F_q, q = 3^30 (escenario histórico) ----------------
    {
        // Módulo irreducible fijo de grado 30 sobre F_3 (calculado con la versión original).
        auto F = fq_family("q3e30", "3",
                           {2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1},
                           "F_(3^30)");
        elem_muladd(F, 4000);
        elem_inv(F, 100);
        poly_mul(F, MulAlgo::Schoolbook, 32);
        poly_mul(F, MulAlgo::Karatsuba, 128);
        poly_divmod(F, 64, 32);
        poly_factor(F, {{5, 1}}, " (histórico: aleatorio de grado 5)");
    }
    // ---------------- F_q, q = 2^8 (módulo AES) ----------------
    {
        auto F = fq_family("q2e8", "2", {1, 1, 0, 1, 1, 0, 0, 0, 1}, "F_(2^8) (x^8+x^4+x^3+x+1)");
        elem_muladd(F, 50000);
        elem_inv(F, 5000);
        poly_mul(F, MulAlgo::Schoolbook, 128);
        poly_mul(F, MulAlgo::Karatsuba, 512);
        poly_is_irreducible(F, 16);
        poly_factor(F, {{1, 2}, {2, 1}, {3, 1}, {4, 2}, {6, 1}, {7, 1}});  // grado 30
    }
    // ---------------- F_q, q = 7^4 ----------------
    {
        auto F = fq_family("q7e4", "7", {3, 4, 1, 2, 1}, "F_(7^4)");
        elem_muladd(F, 50000);
        elem_inv(F, 5000);
        poly_mul(F, MulAlgo::Schoolbook, 128);
        poly_mul(F, MulAlgo::Karatsuba, 512);
        poly_gcd(F, 200);
        poly_factor(F, {{1, 1}, {2, 2}, {3, 1}, {5, 1}, {7, 1}});  // grado 20
    }
    // ---------------- F_q, q = p^2, p = 2^61 - 1 (x^2 + 1; admite FFT) ----------------
    {
        auto F = fq_family("m61e2", P_M61, {1, 0, 1}, "F_((2^61-1)^2) (x^2+1)");
        elem_muladd(F, 100000);
        elem_inv(F, 1000);
        poly_mul(F, MulAlgo::Schoolbook, 128);
        poly_mul(F, MulAlgo::Karatsuba, 512);
        poly_mul(F, MulAlgo::FFT, 1024);
        poly_mul(F, MulAlgo::FFT, 4096);
        poly_mul(F, MulAlgo::Auto, 2048);
        poly_divmod(F, 500, 250);
        poly_gcd(F, 150);
        poly_factor(F, {{1, 2}, {2, 1}, {3, 1}});  // grado 7
    }
    // ---------------- Construcción de F_q con módulo aleatorio ----------------
    fq_build(3, 30);
    fq_build(7, 20);

    // ---------------- Algoritmos rápidos (fast.hpp) ----------------
    {
        auto F = fp_family("ntt", "998244353", "F_998244353");
        F.prefix = "fast.ntt";
        // 3 es raíz primitiva de 998244353.
        fast_dft(F, 3, 10, false, 8);
        fast_dft(F, 3, 11, false, 4);
        fast_dft(F, 3, 12, false, 2);
        fast_dft(F, 3, 12, true, 2);
        fast_toeplitz_inverse(F, 1024);
        fast_toeplitz_inverse(F, 8192);
        fast_toeplitz_vec(F, 4096);
        fast_toeplitz_divmod(F, 2000, 1000);
    }
    {
        auto F = fp_family("p7", "7", "F_7");
        F.prefix = "fast.p7";
        fast_toeplitz_inverse(F, 1024);
        fast_toeplitz_divmod(F, 2000, 1000);
    }
}

// ===========================================================================
// Ejecución, estadística y salida JSON.
// ===========================================================================
struct Options {
    bool list = false;
    std::string filter;
    long reps = 7;
    double min_time = 0.3;
    long warmup = 1;
    long max_reps = 10000;
    std::string json;
};

struct Result {
    const Scenario* sc = nullptr;
    double setup_s = 0;
    std::vector<double> warmup_s, samples_s;
    std::string checksum, error;
};

double median(std::vector<double> v) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

std::string json_str(const std::string& s) {
    std::string o = "\"";
    for (unsigned char c : s) {
        switch (c) {
            case '"': o += "\\\""; break;
            case '\\': o += "\\\\"; break;
            case '\n': o += "\\n"; break;
            case '\r': o += "\\r"; break;
            case '\t': o += "\\t"; break;
            default:
                if (c < 0x20) {
                    char b[8];
                    std::snprintf(b, sizeof b, "\\u%04x", c);
                    o += b;
                } else {
                    o += static_cast<char>(c);
                }
        }
    }
    return o + "\"";
}

std::string json_num(double x) {
    char b[40];
    std::snprintf(b, sizeof b, "%.9g", x);
    return b;
}

std::string json_arr(const std::vector<double>& v) {
    std::string o = "[";
    for (std::size_t i = 0; i < v.size(); ++i) o += (i ? ", " : "") + json_num(v[i]);
    return o + "]";
}

std::string compiler_string() {
#if defined(__clang__)
    return std::string("clang ") + __clang_version__;
#elif defined(__GNUC__)
    return std::string("g++ ") + __VERSION__;
#elif defined(_MSC_VER)
    return "msvc " + std::to_string(_MSC_VER);
#else
    return "desconocido";
#endif
}

bool write_json(const std::string& path, const Options& o, const std::vector<Result>& rs, double total_s) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    std::string s = "{\n";
    s += "  \"suite\": \"alcp-bench-suite\",\n";
    s += "  \"suite_version\": 1,\n";
    s += "  \"compiler\": " + json_str(compiler_string()) + ",\n";
    s += "  \"cplusplus\": " + std::to_string(static_cast<long>(__cplusplus)) + ",\n";
    s += "  \"gmp_version\": " + json_str(gmp_version) + ",\n";
#ifdef NDEBUG
    s += "  \"ndebug\": true,\n";
#else
    s += "  \"ndebug\": false,\n";
#endif
    s += "  \"reps\": " + std::to_string(o.reps) + ",\n";
    s += "  \"min_time_s\": " + json_num(o.min_time) + ",\n";
    s += "  \"warmup\": " + std::to_string(o.warmup) + ",\n";
    s += "  \"filter\": " + json_str(o.filter) + ",\n";
    s += "  \"total_s\": " + json_num(total_s) + ",\n";
    s += "  \"scenarios\": [\n";
    for (std::size_t i = 0; i < rs.size(); ++i) {
        const auto& r = rs[i];
        s += "    {\"name\": " + json_str(r.sc->name) + ", \"group\": " + json_str(r.sc->group) +
             ", \"op\": " + json_str(r.sc->op) + ", \"desc\": " + json_str(r.sc->desc) +
             ", \"inner\": " + std::to_string(r.sc->inner) + ",\n";
        s += "     \"setup_s\": " + json_num(r.setup_s) + ", \"warmup_s\": " + json_arr(r.warmup_s) +
             ", \"median_s\": " + json_num(median(r.samples_s)) + ", \"min_s\": " +
             json_num(r.samples_s.empty() ? 0 : *std::min_element(r.samples_s.begin(), r.samples_s.end())) +
             ",\n";
        s += "     \"checksum\": " + json_str(r.checksum) + ", \"error\": " + json_str(r.error) + ",\n";
        s += "     \"samples_s\": " + json_arr(r.samples_s) + "}" + (i + 1 < rs.size() ? "," : "") + "\n";
    }
    s += "  ]\n}\n";
    const bool ok = std::fwrite(s.data(), 1, s.size(), f) == s.size();
    return std::fclose(f) == 0 && ok;
}

void usage() {
    std::puts(
        "uso: bench_suite [--list] [--filter REGEX] [--reps N] [--min-time S]\n"
        "                 [--warmup N] [--max-reps N] [--json RUTA]\n"
        "  --list         lista los escenarios (nombre, grupo, operación, descripción)\n"
        "  --filter RE    solo escenarios cuyo nombre contiene una coincidencia de RE (ECMAScript)\n"
        "  --reps N       repeticiones mínimas cronometradas (por defecto 7)\n"
        "  --min-time S   tiempo acumulado mínimo por escenario en segundos (por defecto 0.3)\n"
        "  --warmup N     ejecuciones de calentamiento no medidas (por defecto 1)\n"
        "  --max-reps N   tope de repeticiones por escenario (por defecto 10000)\n"
        "  --json RUTA    guarda todas las muestras y checksums en JSON");
}

bool parse(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto need = [&](const char* what) -> const char* {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "falta el valor de %s\n", what);
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "--list") o.list = true;
        else if (a == "--filter") o.filter = need("--filter");
        else if (a == "--reps") o.reps = std::max(1L, std::atol(need("--reps")));
        else if (a == "--min-time") o.min_time = std::max(0.0, std::atof(need("--min-time")));
        else if (a == "--warmup") o.warmup = std::max(0L, std::atol(need("--warmup")));
        else if (a == "--max-reps") o.max_reps = std::max(1L, std::atol(need("--max-reps")));
        else if (a == "--json") o.json = need("--json");
        else if (a == "-h" || a == "--help") { usage(); std::exit(0); }
        else {
            std::fprintf(stderr, "argumento desconocido: %s\n", a.c_str());
            usage();
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    Options opt;
    if (!parse(argc, argv, opt)) return 2;
    register_all();

    std::regex re;
    try {
        re = std::regex(opt.filter.empty() ? std::string(".*") : opt.filter, std::regex::ECMAScript);
    } catch (const std::regex_error& e) {
        std::fprintf(stderr, "expresión regular inválida '%s': %s\n", opt.filter.c_str(), e.what());
        return 2;
    }
    std::vector<const Scenario*> sel;
    for (const auto& sc : registry())
        if (std::regex_search(sc.name, re)) sel.push_back(&sc);

    if (opt.list) {
        for (const auto* sc : sel)
            std::printf("%s\t%s\t%s\t%s\n", sc->name.c_str(), sc->group.c_str(), sc->op.c_str(), sc->desc.c_str());
        return 0;
    }
    if (sel.empty()) {
        std::fprintf(stderr, "ningún escenario coincide con el filtro '%s'\n", opt.filter.c_str());
        return 2;
    }

    std::printf("# alcp bench_suite | %s | GMP %s | %zu escenarios | reps>=%ld, min-time %.2f s\n",
                compiler_string().c_str(), gmp_version, sel.size(), opt.reps, opt.min_time);
    std::fflush(stdout);

    const auto t_all = Clock::now();
    std::vector<Result> results;
    bool any_error = false;
    for (const auto* sc : sel) {
        Result r;
        r.sc = sc;
        const std::uint64_t seed = fnv_str(sc->name);
        const auto lib_seed = static_cast<unsigned long>(seed ^ (seed >> 32));
        try {
            set_seed(lib_seed);
            auto t0 = Clock::now();
            Runner run = sc->make(seed);
            r.setup_s = seconds_since(t0);
            for (long w = 0; w < opt.warmup; ++w) {
                run.reset();
                set_seed(lib_seed);
                t0 = Clock::now();
                run.run();
                r.warmup_s.push_back(seconds_since(t0));
            }
            double total = 0;
            while ((static_cast<long>(r.samples_s.size()) < opt.reps || total < opt.min_time) &&
                   static_cast<long>(r.samples_s.size()) < opt.max_reps) {
                run.reset();
                set_seed(lib_seed);
                t0 = Clock::now();
                run.run();
                const double dt = seconds_since(t0);
                r.samples_s.push_back(dt);
                total += dt;
            }
            r.checksum = run.checksum();
        } catch (const std::exception& e) {
            r.error = e.what();
        }
        if (!r.error.empty() || r.checksum.rfind("INVALIDO", 0) == 0) any_error = true;

        if (r.error.empty()) {
            const double mn = *std::min_element(r.samples_s.begin(), r.samples_s.end());
            std::printf("%-46s %11.3f ms  min %11.3f ms  n=%-5zu %s\n", sc->name.c_str(), 1e3 * median(r.samples_s),
                        1e3 * mn, r.samples_s.size(), r.checksum.c_str());
        } else {
            std::printf("%-46s ERROR: %s\n", sc->name.c_str(), r.error.c_str());
        }
        std::fflush(stdout);
        results.push_back(std::move(r));
    }
    const double total_s = seconds_since(t_all);
    std::printf("# total %.2f s\n", total_s);

    if (!opt.json.empty()) {
        if (!write_json(opt.json, opt, results, total_s)) {
            std::fprintf(stderr, "no se pudo escribir %s\n", opt.json.c_str());
            return 3;
        }
    }
    return any_error ? 1 : 0;
}
