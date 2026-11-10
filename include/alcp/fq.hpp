// alcp/fq.hpp — The extension field F_q = F_p[x]/(g), q = p^n.
//
// Port of the Python `cuerpo_fq` kernel.  An element is a length-n vector of
// base-field coefficients in [0, p), i.e. a residue modulo the monic
// irreducible g.  Fq models the SAME Field interface as Fp, so Poly<Fq> and
// every algorithm in polynomial.hpp / fft.hpp apply to F_q[x] with no extra code.
//
// Implementación genérica FqT<Base> sobre el cuerpo primo Base:
//   * Fq   = FqT<Fp>    coeficientes BigInt (GMP), cualquier p: el tipo de siempre,
//                       con la misma API (Elem = std::vector<BigInt>);
//   * Fq64 = FqT<Fp64>  coeficientes de 64 bits, p < 2^62 (fp64.hpp).
// Las dos versiones hacen exactamente las mismas operaciones: la aritmética de
// los coeficientes la aporta Base (incluidos sus ganchos de acumulación
// perezosa, field.hpp), y la API pública habla en BigInt (modulus(),
// from_vector, from_index/to_index, characteristic, cardinality) en ambas.
//
// The base field, prime power and a monic irreducible modulus g are fixed at
// construction; if g is not supplied, one is found by randomized search.
//
// Validación (antes de calcular nada): el grado de extensión n debe cumplir
// 1 <= n <= Fq::kMaxExtensionDegree; en caso contrario se lanza
// std::invalid_argument (antes, un módulo vacío provocaba un desbordamiento de
// g.size() - 1 y Fq(p, 0) un bucle infinito en build_random).
//
// Elementos: como en Fp, se asume (precondición, no se comprueba por
// rendimiento) que los Elem que se pasan a las operaciones tienen longitud n y
// coeficientes en [0, p).  Para convertir un vector arbitrario de enteros en un
// elemento válido está from_vector().
#pragma once

#include "bigint.hpp"
#include "errors.hpp"
#include "fft_tables.hpp"
#include "field.hpp"
#include "fp.hpp"
#include "fp64.hpp"
#include "polynomial.hpp"

#include <atomic>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace alcp {

template <class Base>
    requires FiniteField<Base> && LazyAccumulation<Base>
class FqT {
public:
    using BaseField = Base;
    using BaseElem = typename Base::Elem;
    using Elem = std::vector<BaseElem>;  // length n, each coeff in [0, p)

    // Cota de cordura para el grado de extensión n.  No es un límite de
    // viabilidad (construir F_{p^n} con n de varios miles ya es lentísimo),
    // sino un tope que rechaza de inmediato valores absurdos (p. ej. n = 2^32-1)
    // que agotarían la memoria al reservar el módulo o al calcular q = p^n.
    static constexpr unsigned long kMaxExtensionDegree = 1UL << 16;

    // q = p^n with p prime; a random monic irreducible g of degree n is chosen.
    explicit FqT(const BigInt& q) : FqT(decompose(q)) {}
    // explicit p and n.
    FqT(const BigInt& p, unsigned long n) : FqT(build_random(Base(p), n)) {}
    // explicit base field and modulus (coeffs length n+1; made monic; must be
    // irreducible).  This is the core constructor.
    FqT(Base base, std::vector<BigInt> g)
        : base_(std::move(base)),
          p_(base_.characteristic()),
          n_(checked_degree(g)),  // valida ANTES de calcular q = p^n
          q_(ipow(p_, n_)),
          g_(std::move(g)) {
        normalize_and_check();
    }

    // --- field data ---------------------------------------------------------
    const Base& base_field() const { return base_; }
    const BigInt& characteristic() const { return p_; }
    BigInt cardinality() const { return q_; }
    unsigned long extension_degree() const { return n_; }
    const std::vector<BigInt>& modulus() const { return g_; }  // length n+1, monic

    // Igualdad de cuerpos: dos F_q son el mismo cuerpo (la misma presentación
    // F_p[x]/(g)) si tienen la misma característica y el mismo módulo mónico.
    // Dos módulos irreducibles distintos dan cuerpos isomorfos pero con
    // elementos codificados de forma distinta, así que NO se consideran iguales.
    bool operator==(const FqT& o) const { return p_ == o.p_ && g_ == o.g_; }

    // --- constants & coercion ----------------------------------------------
    Elem zero() const { return Elem(n_, base_.zero()); }
    Elem one() const { Elem e(n_, base_.zero()); e[0] = base_.one(); return e; }

    // Scalar embedding of an integer: (a mod p) * 1.  (NOT the base-p index
    // encoding — that is from_index/to_index below.)
    Elem from_int(const BigInt& a) const {
        Elem e = zero();
        e[0] = base_.from_int(a);
        return e;
    }
    // n llamadas a Base::random (rand_below(p) cada una), como siempre.
    Elem random() const {
        Elem e(n_);
        for (auto& x : e) x = base_.random();
        return e;
    }

    // Conversión generosa de un vector de enteros (coeficientes en a, orden
    // creciente: v[i] es el coeficiente de a^i) a un elemento válido:
    //   * cada coeficiente se reduce módulo p (admite negativos y valores >= p);
    //   * si hay más de n coeficientes, el polinomio se reduce módulo g;
    //   * si hay menos, se rellena con ceros.
    // El vector vacío representa el 0.
    Elem from_vector(std::vector<BigInt> v) const {
        for (auto& c : v) {
            c %= p_;
            if (sgn(c) < 0) c += p_;
        }
        // Reducción módulo el g mónico, de mayor a menor grado:
        // c·a^k = c·a^{k-n}·a^n y a^n = -(g_0 + ... + g_{n-1} a^{n-1}).
        for (std::size_t k = v.size(); k-- > n_; ) {
            const BigInt c = v[k];
            if (c == 0) continue;
            const std::size_t base = k - n_;
            for (unsigned long t = 0; t < n_; ++t) {
                BigInt& d = v[base + t];
                d -= c * g_[t];
                d %= p_;
                if (sgn(d) < 0) d += p_;
            }
        }
        v.resize(n_, BigInt(0));
        Elem e;
        e.reserve(n_);
        for (const auto& c : v) e.push_back(base_.from_int(c));
        return e;
    }

    // base-p digit encoding, used to enumerate the field (tables, tests, pybind).
    Elem from_index(const BigInt& idx) const {
        Elem e(n_);
        BigInt t = idx % q_;
        if (sgn(t) < 0) t += q_;
        BigInt d;
        for (unsigned long i = 0; i < n_; ++i) {
            mpz_fdiv_qr(t.get_mpz_t(), d.get_mpz_t(), t.get_mpz_t(), p_.get_mpz_t());
            e[i] = base_.from_int(d);
        }
        return e;
    }
    BigInt to_index(const Elem& a) const {
        // Horner en base p: c_0 + p·(c_1 + p·(c_2 + ...)).
        BigInt acc = 0;
        for (unsigned long i = n_; i-- > 0; ) {
            acc *= p_;
            acc += base_.to_bigint(a[i]);
        }
        return acc;
    }

    // --- arithmetic ---------------------------------------------------------
    Elem add(const Elem& a, const Elem& b) const {
        Elem r;
        r.reserve(n_);
        for (unsigned long i = 0; i < n_; ++i) r.push_back(base_.add(a[i], b[i]));
        return r;
    }
    Elem sub(const Elem& a, const Elem& b) const {
        Elem r;
        r.reserve(n_);
        for (unsigned long i = 0; i < n_; ++i) r.push_back(base_.sub(a[i], b[i]));
        return r;
    }
    Elem neg(const Elem& a) const {
        Elem r;
        r.reserve(n_);
        for (unsigned long i = 0; i < n_; ++i) r.push_back(base_.neg(a[i]));
        return r;
    }
    // Producto: convolución de los coeficientes y reducción módulo el g mónico,
    // las mismas operaciones que antes.  Detalles de implementación:
    //   * la convolución vive en un búfer que se reutiliza entre llamadas
    //     (con GMP sus mpz conservan la memoria reservada: no hay malloc por
    //     coeficiente en el régimen estacionario); si otro hilo lo está
    //     usando (región paralela), se usa un acumulador local;
    //   * a[i]·b[j] se acumula sin reducir con el gancho acc_fma de Base
    //     (mpz_addmul en Fp; producto de 128 bits en Fp64);
    //   * cada coeficiente se reduce módulo p una sola vez.
    Elem mul(const Elem& a, const Elem& b) const {
        Scratch& s = conv_scratch();
        if (!s.busy.exchange(true, std::memory_order_acquire)) {
            // Búfer compartido conseguido: se libera también si algo lanza.
            struct Liberar {
                std::atomic<bool>& f;
                ~Liberar() { f.store(false, std::memory_order_release); }
            } liberar{s.busy};
            acc_clear(s.buf);
            acc_fma(s.buf, a, b);
            return acc_reduce(s.buf);
        }
        Acc conv;  // búfer ocupado por otro hilo: acumulador propio
        acc_clear(conv);
        acc_fma(conv, a, b);
        return acc_reduce(conv);
    }

    Elem pow(const Elem& a, const BigInt& k) const {
        Elem base = a;
        BigInt e = k;
        if (e < 0) { base = inv(base); e = -e; }
        Elem result = one();
        while (e > 0) {
            if (mpz_odd_p(e.get_mpz_t())) result = mul(result, base);
            e >>= 1;
            if (e > 0) base = mul(base, base);
        }
        return result;
    }
    // Invertir 0 (también pow(0, k) con k < 0) lanza alcp::division_by_zero.
    // Fase 4: algoritmo de Euclides extendido en F_p[t] sobre (g, a).  Como g
    // es irreducible y a != 0 (mod g), mcd(g, a) es una constante c y la
    // identidad de Bézout s·a + u·g = c da a^{-1} = s·c^{-1} (mod g).  Cuesta
    // O(n^2) operaciones de F_p, frente a las ~2·log2(q) multiplicaciones de
    // F_q (cada una O(n^2)) de Fermat.  El inverso es único: el resultado es
    // el mismo que el de inv_fermat.
    Elem inv(const Elem& a) const {
        if (is_zero(a)) throw division_by_zero("0 is not invertible in F_q");
        // Restos r0 = g, r1 = a y coeficientes de Bézout de a: s0 = 0, s1 = 1
        // (vectores de coeficientes en F_p, orden creciente, sin ceros finales).
        std::vector<BaseElem> r0 = g_base_, r1 = recortado(a);
        std::vector<BaseElem> s0, s1{base_.one()};
        while (r1.size() > 1) {
            // División escolar de r0 entre r1 en F_p[t], acumulando el cociente
            // directamente en s0 -= cociente·s1.
            const BaseElem inv_lead = base_.inv(r1.back());
            const std::size_t nb = r1.size();
            while (r0.size() >= nb) {
                const BaseElem f = base_.mul(r0.back(), inv_lead);
                const std::size_t pos = r0.size() - nb;
                for (std::size_t j = 0; j < nb; ++j)
                    r0[pos + j] = base_.sub(r0[pos + j], base_.mul(f, r1[j]));
                // s0 -= f·t^pos·s1
                if (s0.size() < pos + s1.size()) s0.resize(pos + s1.size(), base_.zero());
                for (std::size_t j = 0; j < s1.size(); ++j)
                    s0[pos + j] = base_.sub(s0[pos + j], base_.mul(f, s1[j]));
                recortar(r0);
            }
            recortar(s0);
            std::swap(r0, r1);  // (r0, r1) <- (r1, r0 mod r1)
            std::swap(s0, s1);
        }
        // r1 es la constante c != 0; a^{-1} = s1·c^{-1}, con deg s1 < n.
        const BaseElem ci = base_.inv(r1[0]);
        Elem r(n_, base_.zero());
        for (std::size_t i = 0; i < s1.size(); ++i) r[i] = base_.mul(s1[i], ci);
        return r;
    }
    // Inverso por el pequeño teorema de Fermat, a^{q-2} (el algoritmo
    // original; se conserva para los tests de equivalencia y las medidas).
    Elem inv_fermat(const Elem& a) const {
        if (is_zero(a)) throw division_by_zero("0 is not invertible in F_q");
        return pow(a, q_ - 2);
    }
    Elem div(const Elem& a, const Elem& b) const { return mul(a, inv(b)); }

    // Unique p-th root: a^{q/p} (inverse Frobenius).  Used by square-free decomp.
    Elem pth_root(const Elem& a) const { return pow(a, q_ / p_); }

    // --- acumulación perezosa (ganchos opcionales, ver field.hpp) -----------
    // El acumulador es la convolución SIN reducir (2n - 1 acumuladores de
    // Base): sumar productos a·b de F_q es sumar sus convoluciones, y la
    // reducción módulo g y módulo p (lo caro) se hace una sola vez, en
    // acc_reduce.
    using Acc = std::vector<typename Base::Acc>;
    void acc_clear(Acc& s) const {
        const std::size_t L = conv_len();
        if (s.size() < L) s.resize(L);
        for (std::size_t k = 0; k < L; ++k) base_.acc_clear(s[k]);
    }
    void acc_set(Acc& s, const Elem& a) const {
        acc_clear(s);
        for (unsigned long i = 0; i < n_; ++i) base_.acc_set(s[i], a[i]);
    }
    // Fase 6: con Base = Fp64 se usa acc_fma_rama (producto de 64 bits si
    // p < 2^32, comprobado por término; ver fp64.hpp).
    void base_fma(typename Base::Acc& s, const BaseElem& a, const BaseElem& b) const {
        if constexpr (requires { base_.acc_fma_rama(s, a, b); }) base_.acc_fma_rama(s, a, b);
        else base_.acc_fma(s, a, b);
    }
    void acc_fma(Acc& s, const Elem& a, const Elem& b) const {
        for (unsigned long i = 0; i < n_; ++i) {
            const BaseElem& ai = a[i];
            if (base_.is_zero(ai)) continue;
            for (unsigned long j = 0; j < n_; ++j) base_fma(s[i + j], ai, b[j]);
        }
    }
    // Reduce la convolución acumulada (la consume: s queda como búfer de
    // trabajo).  Módulo el g mónico: x^n = -(g_0 + ... + g_{n-1} x^{n-1}); el
    // término c·x^k (k >= n) se sustituye sumando c·(-g_t) a x^{k-n+t}.
    Elem acc_reduce(Acc& s) const {
        const std::size_t L = conv_len();
        for (std::size_t k = L; k-- > n_; ) {
            const BaseElem c = base_.acc_reduce(s[k]);
            if (base_.is_zero(c)) continue;
            const std::size_t base = k - n_;
            for (unsigned long t = 0; t < n_; ++t) base_fma(s[base + t], c, neg_g_[t]);
        }
        Elem r;
        r.reserve(n_);
        for (unsigned long i = 0; i < n_; ++i) r.push_back(base_.acc_reduce(s[i]));
        return r;
    }

    // --- predicates ---------------------------------------------------------
    bool eq(const Elem& a, const Elem& b) const { return a == b; }
    bool is_zero(const Elem& a) const {
        for (const auto& x : a) if (!base_.is_zero(x)) return false;
        return true;
    }
    bool is_one(const Elem& a) const {
        if (!base_.is_one(a[0])) return false;
        for (unsigned long i = 1; i < n_; ++i) if (!base_.is_zero(a[i])) return false;
        return true;
    }

    // --- conversion ---------------------------------------------------------
    std::string to_string(const Elem& a) const {
        std::string s;
        bool first = true;
        for (unsigned long i = 0; i < n_; ++i) {
            if (base_.is_zero(a[i])) continue;
            if (!first) s += "+";
            first = false;
            std::string cs = base_.to_string(a[i]);
            if (i == 0) s += cs;
            else if (i == 1) s += (cs == "1" ? "" : cs) + "a";
            else s += (cs == "1" ? "" : cs) + "a^" + std::to_string(i);
        }
        return first ? "0" : s;
    }

    // --- FFT support --------------------------------------------------------
    unsigned long fft_capacity() const {
        if (q_ == 2) return 0;
        return v2(q_ - 1);
    }
    const Elem& fft_base_root() const {
        if (!base_root_ready_) {
            unsigned long cap = fft_capacity();
            BigInt order = BigInt(1) << cap;
            BigInt cofactor = (q_ - 1) / order;
            BigInt half = (cap == 0) ? BigInt(0) : (BigInt(1) << (cap - 1));
            while (true) {
                Elem x = random();
                if (is_zero(x)) continue;
                Elem cand = pow(x, cofactor);
                if (cap == 0) { base_root_ = cand; break; }
                if (!is_one(pow(cand, half))) { base_root_ = cand; break; }
            }
            base_root_ready_ = true;
        }
        return base_root_;
    }
    mutable FFTCache<Elem> fft_cache_;

private:
    Base base_;
    BigInt p_;
    unsigned long n_;
    BigInt q_;
    std::vector<BigInt> g_;           // length n+1, monic (forma canónica, API pública)
    std::vector<BaseElem> neg_g_;     // -g_0, ..., -g_{n-1} en Base (para la reducción)
    std::vector<BaseElem> g_base_;    // g_0, ..., g_n en Base (para el inverso por Euclides)
    mutable bool base_root_ready_ = false;
    mutable Elem base_root_{};

    std::size_t conv_len() const { return 2 * static_cast<std::size_t>(n_) - 1; }

    // Quita los coeficientes nulos finales (forma recortada de un polinomio de F_p[t]).
    void recortar(std::vector<BaseElem>& v) const {
        while (!v.empty() && base_.is_zero(v.back())) v.pop_back();
    }
    std::vector<BaseElem> recortado(const Elem& a) const {
        std::vector<BaseElem> v(a.begin(), a.end());
        recortar(v);
        return v;
    }

    // Búfer de la convolución de mul(), compartido por todos los cuerpos
    // FqT<Base> del proceso y protegido por un indicador atómico (mul no es
    // reentrante: no llama a nada que lo use).  Crece hasta el mayor tamaño
    // pedido y se libera al terminar el programa.
    // No es thread_local a propósito: con MinGW-w64 (g++ 15.2, enlace
    // estático), un thread_local con destructor no trivial corrompe el montón
    // cuando terminan los hilos de std::async de las regiones paralelas
    // (0xC0000374).  Así, el hilo que consigue el búfer lo reutiliza y los
    // demás usan un acumulador local; el resultado es el mismo.
    struct Scratch {
        std::atomic<bool> busy{false};
        Acc buf;
    };
    static Scratch& conv_scratch() {
        static Scratch s;
        return s;
    }

    // Delegating constructor from a (base field, modulus) pair.
    FqT(std::pair<Base, std::vector<BigInt>> pr) : FqT(std::move(pr.first), std::move(pr.second)) {}

    // Comprueba el grado de extensión n antes de usarlo (bucles, reservas de
    // memoria, ipow).  Lanza std::invalid_argument si n == 0 o si n es absurdo.
    static unsigned long check_extension_degree(unsigned long n) {
        if (n == 0)
            throw std::invalid_argument("F_q: el grado de extensión n debe ser >= 1");
        if (n > kMaxExtensionDegree)
            throw std::invalid_argument("F_q: grado de extensión n = " + std::to_string(n) +
                                        " demasiado grande (máximo " +
                                        std::to_string(kMaxExtensionDegree) + ")");
        return n;
    }
    // Grado del módulo g (g.size() - 1), validado.  Un módulo con menos de dos
    // coeficientes (vacío o constante) no define una extensión.
    static unsigned long checked_degree(const std::vector<BigInt>& g) {
        if (g.size() < 2)
            throw std::invalid_argument("F_q: modulus degree must be >= 1");
        if (g.size() - 1 > kMaxExtensionDegree)
            throw std::invalid_argument("F_q: grado del módulo " + std::to_string(g.size() - 1) +
                                        " demasiado grande (máximo " +
                                        std::to_string(kMaxExtensionDegree) + ")");
        return static_cast<unsigned long>(g.size() - 1);
    }

    void normalize_and_check() {
        // make monic
        BigInt lead = g_.back() % p_;
        if (sgn(lead) < 0) lead += p_;
        // Un coeficiente principal múltiplo de p (no solo el 0 literal) también
        // anula el grado; antes acababa en un domain_error de invmod.
        if (lead == 0) throw std::invalid_argument("F_q: modulus has zero leading coeff");
        BigInt inv_lead = invmod(lead, p_);
        for (auto& c : g_) { c = (c % p_) * inv_lead % p_; if (sgn(c) < 0) c += p_; }
        // verify irreducibility via F_p[x]
        Poly<Base> gp = Poly<Base>::from_ints(base_, g_);
        if (gp.degree() != static_cast<int>(n_) || !gp.is_irreducible())
            throw std::invalid_argument("F_q: modulus is not irreducible of degree n");
        neg_g_.reserve(n_);
        for (unsigned long t = 0; t < n_; ++t) neg_g_.push_back(base_.neg(base_.from_int(g_[t])));
        g_base_.reserve(n_ + 1);
        for (unsigned long t = 0; t <= n_; ++t) g_base_.push_back(base_.from_int(g_[t]));
    }

    static std::pair<Base, std::vector<BigInt>> decompose(const BigInt& q) {
        auto [p, n] = prime_power_factor(q);
        return build_random(Base(p), to_ulong(n));
    }
    static std::pair<Base, std::vector<BigInt>> build_random(Base base, unsigned long n) {
        // Con n == 0 ningún polinomio de grado 0 es irreducible: el bucle de
        // abajo no terminaría nunca.  Se valida antes de buscar.
        check_extension_degree(n);
        const BigInt& p = base.characteristic();
        while (true) {
            std::vector<BigInt> coeffs(n + 1);
            for (unsigned long i = 0; i < n; ++i) coeffs[i] = rand_below(p);
            coeffs[n] = 1;  // monic
            Poly<Base> g = Poly<Base>::from_ints(base, coeffs);
            if (g.degree() == static_cast<int>(n) && g.is_irreducible())
                return {std::move(base), std::move(coeffs)};
        }
    }
};

// El cuerpo F_q de siempre (coeficientes GMP) y su variante de 64 bits.
using Fq = FqT<Fp>;
using Fq64 = FqT<Fp64>;

// ¿Son a y b el mismo cuerpo F_q?  Camino rápido: el mismo objeto; si no,
// mismos parámetros (misma p y mismo módulo mónico).
template <class Base>
bool same_field(const FqT<Base>& a, const FqT<Base>& b) { return &a == &b || a == b; }

// Fq y Fq64 cumplen la interfaz de cuerpo finito (concepto de field.hpp).
static_assert(FiniteField<Fq>);
static_assert(LazyAccumulation<Fq>);
static_assert(FiniteField<Fq64>);
static_assert(LazyAccumulation<Fq64>);

}  // namespace alcp
