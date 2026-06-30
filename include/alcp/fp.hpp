// alcp/fp.hpp — The prime field F_p.
//
// Port of the Python `cuerpo_fp` kernel.  An element is just a BigInt kept
// normalized in [0, p).  Fp models the generic Field interface consumed by
// Poly<Field>, PolyRing<Field> and the FFT engine, so every polynomial
// algorithm is written once and works over F_p, F_q, ... unchanged.
//
// Field interface (the contract every field type fulfils):
//   using Elem;
//   Elem  zero() / one() / from_int(BigInt) / random();
//   Elem  add/sub/neg/mul/inv/div(...) ;  Elem pow(Elem, BigInt);
//   bool  eq/is_zero/is_one(...);
//   BigInt cardinality() / characteristic();
//   std::string to_string(Elem);
//   unsigned long fft_capacity();  Elem fft_base_root();
//   FFTCache<Elem> fft_cache_;     // mutable, owned by the field
#pragma once

#include "bigint.hpp"
#include "fft_tables.hpp"

#include <stdexcept>
#include <string>

namespace alcp {

class Fp {
public:
    using Elem = BigInt;

    explicit Fp(BigInt p) : p_(std::move(p)) {
        if (!is_prime(p_))
            throw std::invalid_argument(p_.get_str() + " is not prime; cannot build F_p");
    }

    // --- field data ---------------------------------------------------------
    const BigInt& characteristic() const { return p_; }
    BigInt cardinality() const { return p_; }

    // --- constants & coercion ----------------------------------------------
    Elem zero() const { return Elem(0); }
    Elem one() const { return Elem(1); }

    Elem from_int(const BigInt& a) const {
        Elem r = a % p_;
        if (sgn(r) < 0) r += p_;
        return r;
    }

    Elem random() const { return rand_below(p_); }

    // --- arithmetic ---------------------------------------------------------
    Elem add(const Elem& a, const Elem& b) const {
        Elem r = a + b;
        if (r >= p_) r -= p_;
        return r;
    }
    Elem sub(const Elem& a, const Elem& b) const {
        Elem r = a - b;
        if (sgn(r) < 0) r += p_;
        return r;
    }
    Elem neg(const Elem& a) const { return (a == 0) ? Elem(0) : Elem(p_ - a); }
    Elem mul(const Elem& a, const Elem& b) const { return (a * b) % p_; }

    Elem inv(const Elem& a) const {
        if (a == 0) throw std::domain_error("0 is not invertible in F_p");
        return invmod(a, p_);
    }
    Elem div(const Elem& a, const Elem& b) const { return mul(a, inv(b)); }

    // a^k for any integer k (negative -> inverse first).
    Elem pow(const Elem& a, const BigInt& k) const {
        if (k < 0) {
            if (a == 0) throw std::domain_error("0 is not invertible in F_p");
            return powmod(inv(a), -k, p_);
        }
        return powmod(a, k, p_);
    }

    // --- predicates ---------------------------------------------------------
    bool eq(const Elem& a, const Elem& b) const { return a == b; }
    bool is_zero(const Elem& a) const { return a == 0; }
    bool is_one(const Elem& a) const { return a == 1; }

    // --- conversion ---------------------------------------------------------
    BigInt to_bigint(const Elem& a) const { return a; }
    std::string to_string(const Elem& a) const { return a.get_str(); }

    // Unique p-th root.  In F_p the Frobenius x -> x^p is the identity, so
    // every element is its own p-th root.  (General formula a^{q/p} reduces to
    // a^1 = a here.)  Used by the square-free decomposition.
    Elem pth_root(const Elem& a) const { return a; }

    // --- FFT support --------------------------------------------------------
    // Largest k with 2^k | (p - 1): the biggest power-of-two transform F_p can
    // support.  For p = 2 this is 0 (no nontrivial FFT; callers fall back).
    unsigned long fft_capacity() const {
        if (p_ == 2) return 0;
        return v2(p_ - 1);
    }

    // A fixed primitive 2^capacity-th root of unity (computed once, cached).
    const Elem& fft_base_root() const {
        if (!base_root_ready_) {
            unsigned long cap = fft_capacity();
            BigInt order = BigInt(1) << cap;          // 2^cap
            BigInt cofactor = (p_ - 1) / order;        // (p-1) / 2^cap
            // candidate = x^cofactor has order dividing 2^cap; it is primitive
            // iff candidate^(2^{cap-1}) != 1.  Random x succeeds with prob 1/2.
            BigInt half = (cap == 0) ? BigInt(0) : (BigInt(1) << (cap - 1));
            while (true) {
                Elem x = rand_below(p_ - 1) + 1;       // x in [1, p-1]
                Elem cand = powmod(x, cofactor, p_);
                if (cap == 0) { base_root_ = cand; break; }
                if (powmod(cand, half, p_) != 1) { base_root_ = cand; break; }
            }
            base_root_ready_ = true;
        }
        return base_root_;
    }

    mutable FFTCache<Elem> fft_cache_;

private:
    BigInt p_;
    mutable bool base_root_ready_ = false;
    mutable Elem base_root_{};
};

}  // namespace alcp
