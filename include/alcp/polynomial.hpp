// alcp/polynomial.hpp — Generic univariate polynomials over a finite field.
//
// A single template Poly<Field> carries a pointer to its field K and a
// little-endian, trimmed coefficient vector (c[i] is the coefficient of x^i;
// the zero polynomial has an empty vector, degree -1).  Because every
// operation goes through the Field interface, this one file provides the whole
// algebra of F_p[x] AND F_q[x] (and any other field implementing the
// interface), fulfilling the "algoritmos generales" requirement:
//
//   * ring arithmetic (+, -, *, divmod, /, %, scalar mul, monic, derivative)
//   * multiplication that AUTO-DISPATCHES schoolbook / Karatsuba / FFT by size
//   * gcd, extended gcd, modular inverse, modular exponentiation
//   * Rabin irreducibility test
//   * Cantor-Zassenhaus factorization
//        radical (square-free part) -> distinct-degree -> equal-degree,
//        with per-factor multiplicities  (port of factorizacion.py, unified)
#pragma once

#include "bigint.hpp"
#include "fft.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace alcp {

// Which multiplication algorithm to use.  Auto picks by operand size; the
// explicit modes exist so tests can cross-check them against each other.
enum class MulAlgo { Auto, Schoolbook, Karatsuba, FFT };

template <class Field>
class Poly {
public:
    using Elem = typename Field::Elem;

    // Below this many coefficients in the smaller operand, schoolbook wins.
    static constexpr std::size_t kKaratsubaCutoff = 32;

    // ---- construction ------------------------------------------------------
    explicit Poly(const Field& field) : K_(&field) {}  // zero polynomial
    Poly(const Field& field, std::vector<Elem> coeffs) : K_(&field), c_(std::move(coeffs)) { trim(); }

    static Poly zero(const Field& F) { return Poly(F); }
    static Poly one(const Field& F) { return Poly(F, {F.one()}); }
    static Poly constant(const Field& F, const Elem& v) { return Poly(F, {v}); }
    static Poly variable(const Field& F) { return Poly(F, {F.zero(), F.one()}); }  // x

    // c * x^deg
    static Poly monomial(const Field& F, std::size_t deg, const Elem& coef) {
        std::vector<Elem> v(deg + 1, F.zero());
        v[deg] = coef;
        return Poly(F, std::move(v));
    }
    // Build from a list of integers (handy for F_p[x]); each is coerced.
    static Poly from_ints(const Field& F, const std::vector<BigInt>& ints) {
        std::vector<Elem> v;
        v.reserve(ints.size());
        for (const auto& a : ints) v.push_back(F.from_int(a));
        return Poly(F, std::move(v));
    }
    // Uniformly random polynomial of exact degree d (leading coeff nonzero).
    static Poly random(const Field& F, int d) {
        if (d < 0) return zero(F);
        std::vector<Elem> v(d + 1);
        for (int i = 0; i < d; ++i) v[i] = F.random();
        do { v[d] = F.random(); } while (F.is_zero(v[d]));
        return Poly(F, std::move(v));
    }

    // ---- basics ------------------------------------------------------------
    const Field& field() const { return *K_; }
    const std::vector<Elem>& coeffs() const { return c_; }
    int degree() const { return static_cast<int>(c_.size()) - 1; }
    bool is_zero() const { return c_.empty(); }
    bool is_one() const { return c_.size() == 1 && K_->is_one(c_[0]); }
    const Elem& lead() const { return c_.back(); }  // precondition: !is_zero()

    bool operator==(const Poly& o) const {
        if (c_.size() != o.c_.size()) return false;
        for (std::size_t i = 0; i < c_.size(); ++i)
            if (!K_->eq(c_[i], o.c_[i])) return false;
        return true;
    }
    bool operator!=(const Poly& o) const { return !(*this == o); }

    // ---- additive structure ------------------------------------------------
    Poly operator+(const Poly& o) const {
        std::size_t n = std::max(c_.size(), o.c_.size());
        std::vector<Elem> r(n, K_->zero());
        for (std::size_t i = 0; i < c_.size(); ++i) r[i] = c_[i];
        for (std::size_t i = 0; i < o.c_.size(); ++i) r[i] = K_->add(r[i], o.c_[i]);
        return Poly(*K_, std::move(r));
    }
    Poly operator-() const {
        std::vector<Elem> r(c_.size());
        for (std::size_t i = 0; i < c_.size(); ++i) r[i] = K_->neg(c_[i]);
        return Poly(*K_, std::move(r));
    }
    Poly operator-(const Poly& o) const {
        std::size_t n = std::max(c_.size(), o.c_.size());
        std::vector<Elem> r(n, K_->zero());
        for (std::size_t i = 0; i < c_.size(); ++i) r[i] = c_[i];
        for (std::size_t i = 0; i < o.c_.size(); ++i) r[i] = K_->sub(r[i], o.c_[i]);
        return Poly(*K_, std::move(r));
    }

    // ---- multiplication ----------------------------------------------------
    Poly operator*(const Poly& o) const { return mul(o, MulAlgo::Auto); }

    Poly mul(const Poly& o, MulAlgo algo) const {
        if (is_zero() || o.is_zero()) return zero(*K_);
        const std::size_t na = c_.size(), nb = o.c_.size();
        if (algo == MulAlgo::Auto) {
            const std::size_t L = na + nb - 1;
            const std::size_t small = std::min(na, nb);
            if (fft::fft_exponent(*K_, L) >= 0 && L >= 64)
                algo = MulAlgo::FFT;
            else if (small >= kKaratsubaCutoff)
                algo = MulAlgo::Karatsuba;
            else
                algo = MulAlgo::Schoolbook;
        }
        switch (algo) {
            case MulAlgo::FFT: {
                if (fft::fft_exponent(*K_, na + nb - 1) < 0)
                    return mul_karatsuba(o);  // field cannot host this size
                return Poly(*K_, fft::convolve(*K_, c_, o.c_));
            }
            case MulAlgo::Karatsuba: return mul_karatsuba(o);
            case MulAlgo::Schoolbook:
            default: return mul_school(o);
        }
    }

    Poly mul_school(const Poly& o) const {
        if (is_zero() || o.is_zero()) return zero(*K_);
        std::vector<Elem> r(c_.size() + o.c_.size() - 1, K_->zero());
        for (std::size_t i = 0; i < c_.size(); ++i) {
            if (K_->is_zero(c_[i])) continue;
            for (std::size_t j = 0; j < o.c_.size(); ++j)
                r[i + j] = K_->add(r[i + j], K_->mul(c_[i], o.c_[j]));
        }
        return Poly(*K_, std::move(r));
    }

    // Recursive Karatsuba (port of fp_x_mult_karatsuba / fq_x_mult_karatsuba).
    Poly mul_karatsuba(const Poly& o) const {
        if (is_zero() || o.is_zero()) return zero(*K_);
        if (c_.size() <= kKaratsubaCutoff || o.c_.size() <= kKaratsubaCutoff)
            return mul_school(o);

        std::size_t m = std::max(c_.size(), o.c_.size()) / 2;
        Poly aL = slice(0, m), aH = slice(m, c_.size());
        Poly bL = o.slice(0, m), bH = o.slice(m, o.c_.size());

        Poly z0 = aL.mul_karatsuba(bL);
        Poly z2 = aH.mul_karatsuba(bH);
        Poly z1 = (aL + aH).mul_karatsuba(bL + bH) - z0 - z2;

        return z0 + z1.shift(m) + z2.shift(2 * m);
    }

    // ---- Euclidean division ------------------------------------------------
    // Returns (q, r) with *this = q*o + r and (r == 0 or deg r < deg o).
    std::pair<Poly, Poly> divmod(const Poly& o) const {
        if (o.is_zero()) throw std::domain_error("division by zero polynomial");
        if (degree() < o.degree()) return {zero(*K_), *this};

        const int db = o.degree();
        Elem inv_lead = K_->inv(o.lead());
        std::vector<Elem> r = c_;  // working remainder
        std::vector<Elem> q(degree() - db + 1, K_->zero());

        for (int dr = static_cast<int>(r.size()) - 1; dr >= db; ) {
            // strip any leading zeros that appeared
            while (dr >= db && K_->is_zero(r[dr])) --dr;
            if (dr < db) break;
            int pos = dr - db;
            Elem factor = K_->mul(r[dr], inv_lead);
            q[pos] = factor;
            for (int j = 0; j <= db; ++j)
                r[pos + j] = K_->sub(r[pos + j], K_->mul(factor, o.c_[j]));
            --dr;
        }
        return {Poly(*K_, std::move(q)), Poly(*K_, std::move(r))};
    }
    Poly operator/(const Poly& o) const { return divmod(o).first; }
    Poly operator%(const Poly& o) const { return divmod(o).second; }

    // ---- scalar / shape ----------------------------------------------------
    Poly scalar_mul(const Elem& e) const {
        std::vector<Elem> r(c_.size());
        for (std::size_t i = 0; i < c_.size(); ++i) r[i] = K_->mul(c_[i], e);
        return Poly(*K_, std::move(r));
    }
    Poly monic() const {
        if (is_zero()) return *this;
        return scalar_mul(K_->inv(lead()));
    }
    Poly derivative() const {
        if (c_.size() <= 1) return zero(*K_);
        std::vector<Elem> r(c_.size() - 1);
        for (std::size_t i = 1; i < c_.size(); ++i)
            r[i - 1] = K_->mul(c_[i], K_->from_int(BigInt(static_cast<unsigned long>(i))));
        return Poly(*K_, std::move(r));
    }

    // Plain (non-modular) power.
    Poly pow(const BigInt& k) const {
        if (k < 0) throw std::domain_error("negative power of polynomial");
        Poly result = one(*K_), base = *this;
        BigInt e = k;
        while (e > 0) {
            if (mpz_odd_p(e.get_mpz_t())) result = result * base;
            e >>= 1;
            if (e > 0) base = base * base;
        }
        return result;
    }

    // ---- gcd family --------------------------------------------------------
    static Poly gcd(Poly a, Poly b) {
        while (!b.is_zero()) {
            Poly r = a % b;
            a = std::move(b);
            b = std::move(r);
        }
        return a.is_zero() ? a : a.monic();
    }

    struct ExtGcd { Poly g, x, y; };  // g = x*a + y*b, g monic
    static ExtGcd gcd_ext(const Poly& a, const Poly& b) {
        const Field& K = a.field();
        Poly x0 = one(K), x1 = zero(K), y0 = zero(K), y1 = one(K), r0 = a, r1 = b;
        while (!r1.is_zero()) {
            auto [q, r] = r0.divmod(r1);
            r0 = std::move(r1); r1 = std::move(r);
            Poly tx = x0 - q * x1; x0 = std::move(x1); x1 = std::move(tx);
            Poly ty = y0 - q * y1; y0 = std::move(y1); y1 = std::move(ty);
        }
        if (r0.is_zero()) return {zero(K), zero(K), zero(K)};
        Elem inv = K.inv(r0.lead());
        return {r0.scalar_mul(inv), x0.scalar_mul(inv), y0.scalar_mul(inv)};
    }

    // a^{-1} mod m  (requires gcd(a, m) = 1).
    static Poly inv_mod(const Poly& a, const Poly& m) {
        auto e = gcd_ext(a, m);
        if (e.g.degree() != 0 || !a.field().is_one(e.g.c_[0]))
            throw std::domain_error("polynomial not invertible modulo m");
        return e.x % m;
    }

    // a^k mod m.
    static Poly pow_mod(Poly a, BigInt k, const Poly& m) {
        const Field& K = a.field();
        if (k < 0) { a = inv_mod(a, m); k = -k; }
        Poly result = one(K) % m;
        a = a % m;
        while (k > 0) {
            if (mpz_odd_p(k.get_mpz_t())) result = (result * a) % m;
            k >>= 1;
            if (k > 0) a = (a * a) % m;
        }
        return result;
    }

    // ---- Rabin irreducibility test -----------------------------------------
    // f is irreducible over F_q iff x^{q^n} == x (mod f) and, for every prime
    // divisor's worth of degrees up to n/2, gcd(f, x^{q^i} - x) = 1.
    bool is_irreducible() const {
        const int n = degree();
        if (n <= 0) return false;
        const BigInt q = K_->cardinality();
        Poly f = *this;
        Poly x = variable(*K_);
        Poly H = x % f;                       // H = x^{q^0}
        for (int i = 1; i <= n / 2; ++i) {
            H = pow_mod(H, q, f);             // H = x^{q^i} (mod f)
            if (gcd(f, H - x).degree() != 0) return false;
        }
        Poly Hn = pow_mod(x, ipow(q, static_cast<unsigned long>(n)), f);
        return (Hn - x) % f == zero(*K_);
    }

    // ---- factorization (Cantor-Zassenhaus) ---------------------------------
    struct Factorization {
        std::vector<std::pair<Poly, int>> factors;  // (monic irreducible, multiplicity)
        Elem leading;                                // unit: leading coeff of original f
    };

    // Square-free part: the product of the distinct monic irreducible factors
    // of *this (a.k.a. radical).  Port of sqfree_fact_*.
    Poly radical() const {
        const Field& K = *K_;
        const BigInt p = K.characteristic();
        const unsigned long pu = fits_ulong(p) ? to_ulong(p) : 0;  // step for x^p case

        if (is_zero()) return zero(K);

        Poly df = derivative();
        if (df.is_zero()) {
            // f(x) = g(x^p); recover g via the p-th root of every p-th coeff.
            std::vector<Elem> g;
            for (std::size_t i = 0; i < c_.size(); i += pu) g.push_back(K.pth_root(c_[i]));
            return Poly(K, std::move(g)).radical();
        }

        Poly c = gcd(*this, df);
        Poly g = (*this / c).monic();

        // Strip every factor of g out of f; what remains is a perfect p-th power.
        Poly residue = *this;
        while (true) {
            Poly y = gcd(residue, g);
            if (y.degree() <= 0) break;
            residue = residue / y;
        }
        if (residue.degree() > 0) {
            std::vector<Elem> h;
            for (std::size_t i = 0; i < residue.c_.size(); i += pu) h.push_back(K.pth_root(residue.c_[i]));
            Poly g2 = Poly(K, std::move(h)).radical();
            return (g * g2) / gcd(g, g2);  // lcm(g, g2)
        }
        return g;
    }

    // Distinct-degree factorization of a square-free monic g.  Returns pairs
    // (d, product-of-all-degree-d-irreducible-factors), only for d that occur.
    std::vector<std::pair<int, Poly>> distinct_degree() const {
        const Field& K = *K_;
        const BigInt q = K.cardinality();
        std::vector<std::pair<int, Poly>> out;
        Poly g = monic();
        Poly x = variable(K);
        Poly H = x % g;  // x^{q^0}
        int d = 1;
        while (g.degree() >= 2 * d) {
            H = pow_mod(H, q, g);            // H = x^{q^d} (mod g)
            Poly hd = gcd(H - x, g);
            if (hd.degree() > 0) {
                out.emplace_back(d, hd.monic());
                g = g / hd;
                H = H % g;
            }
            ++d;
        }
        if (g.degree() > 0) out.emplace_back(g.degree(), g.monic());
        return out;
    }

    // Equal-degree factorization: h is a product of distinct monic irreducible
    // factors all of degree r; returns that list of irreducibles.  Works for
    // odd q (the (q^r-1)/2 split) and for q even (the trace split).
    std::vector<Poly> equal_degree(int r) const {
        std::vector<Poly> factors;
        std::vector<Poly> stack;
        stack.push_back(monic());
        while (!stack.empty()) {
            Poly g = std::move(stack.back());
            stack.pop_back();
            if (g.degree() <= 0) continue;
            if (g.degree() == r) { factors.push_back(g.monic()); continue; }
            Poly d = equal_degree_split(g, r);
            stack.push_back(d);
            stack.push_back(g / d);
        }
        return factors;
    }

    // Largest e with u^e | *this.  Port of multiplicidad_*.
    int multiplicity(const Poly& u) const {
        int e = 0;
        Poly f = *this;
        while (!f.is_zero()) {
            auto [q, r] = f.divmod(u);
            if (!r.is_zero()) break;
            ++e;
            f = std::move(q);
        }
        return e;
    }

    // Full Cantor-Zassenhaus factorization into monic irreducibles with
    // multiplicities, plus the leading coefficient (the field unit).
    Factorization factor() const {
        if (is_zero()) throw std::domain_error("cannot factor the zero polynomial");
        Factorization F;
        F.leading = lead();
        if (degree() == 0) return F;

        Poly g = radical();
        std::vector<Poly> irreducibles;
        for (auto& [d, prod] : g.distinct_degree()) {
            auto part = prod.equal_degree(d);
            irreducibles.insert(irreducibles.end(), part.begin(), part.end());
        }
        for (auto& u : irreducibles)
            F.factors.emplace_back(u, multiplicity(u));
        return F;
    }

    // ---- formatting --------------------------------------------------------
    std::string to_string(const std::string& var = "x") const {
        if (is_zero()) return "0";
        std::string s;
        bool first = true;
        for (std::size_t i = 0; i < c_.size(); ++i) {
            if (K_->is_zero(c_[i])) continue;
            if (!first) s += " + ";
            first = false;
            std::string cs = "(" + K_->to_string(c_[i]) + ")";
            if (i == 0) s += cs;
            else if (i == 1) s += cs + "*" + var;
            else s += cs + "*" + var + "^" + std::to_string(i);
        }
        return s;
    }

private:
    const Field* K_;
    std::vector<Elem> c_;

    void trim() {
        while (!c_.empty() && K_->is_zero(c_.back())) c_.pop_back();
    }
    // coefficients [lo, hi) as a polynomial (drops the x^lo..x^{hi-1} window down to x^0).
    Poly slice(std::size_t lo, std::size_t hi) const {
        if (lo >= c_.size()) return zero(*K_);
        hi = std::min(hi, c_.size());
        return Poly(*K_, std::vector<Elem>(c_.begin() + lo, c_.begin() + hi));
    }
    // multiply by x^s
    Poly shift(std::size_t s) const {
        if (is_zero()) return *this;
        std::vector<Elem> r(c_.size() + s, K_->zero());
        for (std::size_t i = 0; i < c_.size(); ++i) r[i + s] = c_[i];
        return Poly(*K_, std::move(r));
    }

    // One Cantor-Zassenhaus splitting step: given g (product of distinct
    // degree-r irreducibles, deg g > r), return a nontrivial factor.
    static Poly equal_degree_split(const Poly& g, int r) {
        const Field& K = g.field();
        const BigInt q = K.cardinality();
        const bool even_char = (K.characteristic() == 2);
        // For char 2 we need the number of squarings in the trace map:
        // s = log2(q), trace over F_{q^r}/F_2 uses s*r squarings.
        const unsigned long s = even_char ? (bit_length(q) - 1) : 0;

        while (true) {
            Poly a = Poly::random(K, g.degree() - 1) % g;
            if (a.degree() <= 0) continue;
            Poly b(K);
            if (!even_char) {
                BigInt e = (ipow(q, static_cast<unsigned long>(r)) - 1) / 2;
                b = pow_mod(a, e, g) - one(K);
            } else {
                Poly trace = zero(K), cur = a;
                for (unsigned long i = 0; i < s * static_cast<unsigned long>(r); ++i) {
                    trace = (trace + cur) % g;
                    cur = (cur * cur) % g;
                }
                b = trace;
            }
            Poly d = gcd(b, g);
            if (d.degree() > 0 && d.degree() < g.degree()) return d;
        }
    }
};

}  // namespace alcp
