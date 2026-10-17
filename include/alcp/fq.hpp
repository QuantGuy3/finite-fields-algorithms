// alcp/fq.hpp — The extension field F_q = F_p[x]/(g), q = p^n.
//
// Port of the Python `cuerpo_fq` kernel.  An element is a length-n vector of
// BigInt coefficients in [0, p), i.e. a residue modulo the monic irreducible g.
// Fq models the SAME Field interface as Fp, so Poly<Fq> and every algorithm in
// polynomial.hpp / fft.hpp apply to F_q[x] with no extra code.
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
#include "polynomial.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace alcp {

class Fq {
public:
    using Elem = std::vector<BigInt>;  // length n, each coeff in [0, p)

    // Cota de cordura para el grado de extensión n.  No es un límite de
    // viabilidad (construir F_{p^n} con n de varios miles ya es lentísimo),
    // sino un tope que rechaza de inmediato valores absurdos (p. ej. n = 2^32-1)
    // que agotarían la memoria al reservar el módulo o al calcular q = p^n.
    static constexpr unsigned long kMaxExtensionDegree = 1UL << 16;

    // q = p^n with p prime; a random monic irreducible g of degree n is chosen.
    explicit Fq(const BigInt& q) : Fq(decompose(q)) {}
    // explicit p and n.
    Fq(const BigInt& p, unsigned long n) : Fq(build_random(Fp(p), n)) {}
    // explicit base field and modulus (coeffs length n+1; made monic; must be
    // irreducible).  This is the core constructor.
    Fq(Fp base, std::vector<BigInt> g)
        : base_(std::move(base)),
          p_(base_.characteristic()),
          n_(checked_degree(g)),  // valida ANTES de calcular q = p^n
          q_(ipow(p_, n_)),
          g_(std::move(g)) {
        normalize_and_check();
    }

    // --- field data ---------------------------------------------------------
    const Fp& base_field() const { return base_; }
    const BigInt& characteristic() const { return p_; }
    BigInt cardinality() const { return q_; }
    unsigned long extension_degree() const { return n_; }
    const std::vector<BigInt>& modulus() const { return g_; }  // length n+1, monic

    // Igualdad de cuerpos: dos F_q son el mismo cuerpo (la misma presentación
    // F_p[x]/(g)) si tienen la misma característica y el mismo módulo mónico.
    // Dos módulos irreducibles distintos dan cuerpos isomorfos pero con
    // elementos codificados de forma distinta, así que NO se consideran iguales.
    bool operator==(const Fq& o) const { return p_ == o.p_ && g_ == o.g_; }

    // --- constants & coercion ----------------------------------------------
    Elem zero() const { return Elem(n_, BigInt(0)); }
    Elem one() const { Elem e(n_, BigInt(0)); e[0] = 1; return e; }

    // Scalar embedding of an integer: (a mod p) * 1.  (NOT the base-p index
    // encoding — that is from_index/to_index below.)
    Elem from_int(const BigInt& a) const {
        Elem e = zero();
        BigInt r = a % p_;
        if (sgn(r) < 0) r += p_;
        e[0] = r;
        return e;
    }
    Elem random() const {
        Elem e(n_);
        for (auto& x : e) x = rand_below(p_);
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
        return v;
    }

    // base-p digit encoding, used to enumerate the field (tables, tests, pybind).
    Elem from_index(const BigInt& idx) const {
        Elem e(n_);
        BigInt t = idx % q_;
        if (sgn(t) < 0) t += q_;
        for (unsigned long i = 0; i < n_; ++i) { e[i] = t % p_; t /= p_; }
        return e;
    }
    BigInt to_index(const Elem& a) const {
        BigInt acc = 0, base = 1;
        for (unsigned long i = 0; i < n_; ++i) { acc += a[i] * base; base *= p_; }
        return acc;
    }

    // --- arithmetic ---------------------------------------------------------
    Elem add(const Elem& a, const Elem& b) const {
        Elem r(n_);
        for (unsigned long i = 0; i < n_; ++i) { r[i] = a[i] + b[i]; if (r[i] >= p_) r[i] -= p_; }
        return r;
    }
    Elem sub(const Elem& a, const Elem& b) const {
        Elem r(n_);
        for (unsigned long i = 0; i < n_; ++i) { r[i] = a[i] - b[i]; if (sgn(r[i]) < 0) r[i] += p_; }
        return r;
    }
    Elem neg(const Elem& a) const {
        Elem r(n_);
        for (unsigned long i = 0; i < n_; ++i) r[i] = (a[i] == 0) ? BigInt(0) : (p_ - a[i]);
        return r;
    }
    Elem mul(const Elem& a, const Elem& b) const {
        std::vector<BigInt> conv(2 * n_ - 1, BigInt(0));
        for (unsigned long i = 0; i < n_; ++i) {
            if (a[i] == 0) continue;
            for (unsigned long j = 0; j < n_; ++j)
                conv[i + j] += a[i] * b[j];
        }
        // reduce modulo the monic g: x^n = -(g_0 + ... + g_{n-1} x^{n-1}).
        for (std::size_t k = conv.size(); k-- > n_; ) {
            BigInt c = conv[k] % p_;
            if (c == 0) continue;
            std::size_t base = k - n_;
            for (unsigned long t = 0; t < n_; ++t)
                conv[base + t] -= c * g_[t];
        }
        Elem r(n_);
        for (unsigned long i = 0; i < n_; ++i) {
            BigInt v = conv[i] % p_;
            if (sgn(v) < 0) v += p_;
            r[i] = v;
        }
        return r;
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
    Elem inv(const Elem& a) const {
        if (is_zero(a)) throw division_by_zero("0 is not invertible in F_q");
        return pow(a, q_ - 2);  // Fermat: a^{q-2} = a^{-1}
    }
    Elem div(const Elem& a, const Elem& b) const { return mul(a, inv(b)); }

    // Unique p-th root: a^{q/p} (inverse Frobenius).  Used by square-free decomp.
    Elem pth_root(const Elem& a) const { return pow(a, q_ / p_); }

    // --- predicates ---------------------------------------------------------
    bool eq(const Elem& a, const Elem& b) const { return a == b; }
    bool is_zero(const Elem& a) const {
        for (const auto& x : a) if (x != 0) return false;
        return true;
    }
    bool is_one(const Elem& a) const {
        if (a[0] != 1) return false;
        for (unsigned long i = 1; i < n_; ++i) if (a[i] != 0) return false;
        return true;
    }

    // --- conversion ---------------------------------------------------------
    std::string to_string(const Elem& a) const {
        std::string s;
        bool first = true;
        for (unsigned long i = 0; i < n_; ++i) {
            if (a[i] == 0) continue;
            if (!first) s += "+";
            first = false;
            std::string cs = a[i].get_str();
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
    Fp base_;
    BigInt p_;
    unsigned long n_;
    BigInt q_;
    std::vector<BigInt> g_;  // length n+1, monic
    mutable bool base_root_ready_ = false;
    mutable Elem base_root_{};

    // Delegating constructor from a (base field, modulus) pair.
    Fq(std::pair<Fp, std::vector<BigInt>> pr) : Fq(std::move(pr.first), std::move(pr.second)) {}

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
        Poly<Fp> gp = Poly<Fp>::from_ints(base_, g_);
        if (gp.degree() != static_cast<int>(n_) || !gp.is_irreducible())
            throw std::invalid_argument("F_q: modulus is not irreducible of degree n");
    }

    static std::pair<Fp, std::vector<BigInt>> decompose(const BigInt& q) {
        auto [p, n] = prime_power_factor(q);
        return build_random(Fp(p), to_ulong(n));
    }
    static std::pair<Fp, std::vector<BigInt>> build_random(Fp base, unsigned long n) {
        // Con n == 0 ningún polinomio de grado 0 es irreducible: el bucle de
        // abajo no terminaría nunca.  Se valida antes de buscar.
        check_extension_degree(n);
        const BigInt& p = base.characteristic();
        while (true) {
            std::vector<BigInt> coeffs(n + 1);
            for (unsigned long i = 0; i < n; ++i) coeffs[i] = rand_below(p);
            coeffs[n] = 1;  // monic
            Poly<Fp> g = Poly<Fp>::from_ints(base, coeffs);
            if (g.degree() == static_cast<int>(n) && g.is_irreducible())
                return {std::move(base), std::move(coeffs)};
        }
    }
};

// ¿Son a y b el mismo cuerpo F_q?  Camino rápido: el mismo objeto; si no,
// mismos parámetros (misma p y mismo módulo mónico).
inline bool same_field(const Fq& a, const Fq& b) { return &a == &b || a == b; }

// Fq cumple la interfaz de cuerpo finito (concepto de field.hpp).
static_assert(FiniteField<Fq>);

}  // namespace alcp
