// alcp/bigint.hpp — Arbitrary-precision integer backend and number-theoretic helpers.
//
// The BigInt type is abstracted behind a single typedef so the backend can be
// swapped.  The default backend is GMP's mpz_class (gold-standard, already
// installed on this machine), exposed as `alcp::BigInt`.
//
// This header also provides the integer-level algorithms that the original
// Python `cuerpos_finitos.py` implemented by hand:
//   * es_primo        -> is_prime
//   * factor_potencia -> prime_power_factor
// plus the small toolbox (modular exponentiation, 2-adic valuation, randomness)
// that the field and polynomial layers rely on.
#pragma once

#include "errors.hpp"

#include <gmpxx.h>

#include <climits>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>

namespace alcp {

// ---------------------------------------------------------------------------
// BigInt backend.  Swap this single line to change the arbitrary-precision
// integer implementation; everything downstream uses `alcp::BigInt`.
// ---------------------------------------------------------------------------
using BigInt = mpz_class;

// ---------------------------------------------------------------------------
// Randomness.  A single process-wide GMP random state, seedable for
// reproducible tests via set_seed().
// ---------------------------------------------------------------------------
namespace detail {
inline gmp_randclass& gmp_rng() {
    static gmp_randclass state(gmp_randinit_mt);
    static bool seeded = [] {
        state.seed(0x9E3779B9UL);  // deterministic default seed
        return true;
    }();
    (void)seeded;
    return state;
}
}  // namespace detail

inline void set_seed(unsigned long s) { detail::gmp_rng().seed(s); }

// Uniform random integer in [0, n).  Requires n > 0.
inline BigInt rand_below(const BigInt& n) {
    if (n <= 0) throw std::invalid_argument("rand_below: n must be positive");
    return detail::gmp_rng().get_z_range(n);
}

// ---------------------------------------------------------------------------
// Small integer helpers built directly on the GMP C API for speed.
// ---------------------------------------------------------------------------

// base^exp for a non-negative machine-word exponent.
inline BigInt ipow(const BigInt& base, unsigned long exp) {
    BigInt r;
    mpz_pow_ui(r.get_mpz_t(), base.get_mpz_t(), exp);
    return r;
}

// base^exp mod m  (m > 0, exp >= 0).
inline BigInt powmod(const BigInt& base, const BigInt& exp, const BigInt& m) {
    if (exp < 0) throw std::invalid_argument("powmod: negative exponent");
    // mpz_powm con módulo 0 divide por cero dentro de GMP (aborta el proceso).
    if (sgn(m) <= 0) throw std::invalid_argument("powmod: el módulo debe ser positivo");
    BigInt r;
    mpz_powm(r.get_mpz_t(), base.get_mpz_t(), exp.get_mpz_t(), m.get_mpz_t());
    return r;
}

// 2-adic valuation v2(m): the largest k with 2^k | m.  v2(0) is undefined here.
inline unsigned long v2(const BigInt& m) {
    if (m == 0) throw std::invalid_argument("v2: undefined for 0");
    return mpz_scan1(m.get_mpz_t(), 0);
}

// Number of bits in |m| (0 -> 0).
inline std::size_t bit_length(const BigInt& m) {
    if (m == 0) return 0;
    return mpz_sizeinbase(m.get_mpz_t(), 2);
}

// Does m fit in an unsigned long?
inline bool fits_ulong(const BigInt& m) {
    return mpz_fits_ulong_p(m.get_mpz_t()) != 0;
}

inline unsigned long to_ulong(const BigInt& m) {
    if (!fits_ulong(m)) throw std::overflow_error("to_ulong: value too large");
    return mpz_get_ui(m.get_mpz_t());
}

// ---------------------------------------------------------------------------
// Primality.  GMP's Baillie-PSW + Miller-Rabin test (`mpz_probab_prime_p`)
// with 40 rounds: 2 means "proven prime" (small numbers), 1 "probably prime"
// with error probability < 4^-40, 0 "composite".  This replaces the original
// O(sqrt(p)) trial division, which was hopeless for the large primes a BigInt
// field is meant to support.
// ---------------------------------------------------------------------------
inline bool is_prime(const BigInt& p) {
    if (p < 2) return false;
    return mpz_probab_prime_p(p.get_mpz_t(), 40) != 0;
}

// ---------------------------------------------------------------------------
// Prime-power factorization.  Given q >= 2, return (p, n) such that q = p^n
// with p prime, or throw if q is not a prime power.  Mirrors the Python
// `factor_potencia`, but the heavy lifting (root extraction, primality) is
// done with GMP so it stays correct and fast for large q.
//
// Detección de potencia perfecta (sustituye a la antigua división por
// tentativa, que no terminaba en la práctica para bases primas grandes como
// (2^61-1)^3): si q = p^n con p >= 2, entonces 2^n <= q < 2^bits, luego
// n <= bits - 1.  Se recorre n de mayor a menor, se extrae la raíz n-ésima
// entera con mpz_root y, si es exacta y prima, se devuelve (r, n).  Si q = p^n
// con p primo, la única raíz exacta con base prima es la de exponente n, así
// que el primer acierto (el mayor n) es la respuesta.  Coste: O(log q) raíces
// enteras, y ninguna si mpz_perfect_power_p descarta q de entrada.
// ---------------------------------------------------------------------------
inline std::pair<BigInt, BigInt> prime_power_factor(const BigInt& q) {
    if (q < 2) throw std::invalid_argument("prime_power_factor: q must be >= 2");
    if (is_prime(q)) return {q, BigInt(1)};

    if (mpz_perfect_power_p(q.get_mpz_t()) != 0) {
        const std::size_t bits = bit_length(q);
        // mpz_root recibe el exponente como unsigned long (32 bits en Windows);
        // un q con más de ULONG_MAX bits no cabe en memoria, pero se acota igual.
        const unsigned long top = (bits - 1 > static_cast<std::size_t>(ULONG_MAX))
                                      ? ULONG_MAX
                                      : static_cast<unsigned long>(bits - 1);
        BigInt r;
        for (unsigned long n = top; n >= 2; --n) {
            if (mpz_root(r.get_mpz_t(), q.get_mpz_t(), n) != 0 && is_prime(r))
                return {r, BigInt(n)};
        }
    }
    throw std::invalid_argument(q.get_str() + " is not a prime power");
}

// Modular inverse: a^{-1} mod m, throws if gcd(a, m) != 1.
inline BigInt invmod(const BigInt& a, const BigInt& m) {
    // Con m == 0 el comportamiento de mpz_invert no está definido.
    if (sgn(m) <= 0) throw std::invalid_argument("invmod: el módulo debe ser positivo");
    BigInt r;
    if (mpz_invert(r.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t()) == 0)
        throw std::domain_error("invmod: value not invertible");
    return r;
}

}  // namespace alcp
