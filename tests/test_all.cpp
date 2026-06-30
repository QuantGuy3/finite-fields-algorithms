// test_all.cpp — the verification suite.
//
// "Formal-ish" verification strategy:
//   * EXHAUSTIVE axiom checking over small fields (every element / every pair):
//     the field and ring axioms are proven by brute force, not sampled.
//   * RANDOMIZED property tests over large fields: Frobenius, Bezout identity,
//     divmod identity, multiplication-algorithm equivalence, FFT round-trips.
//   * A THEOREM cross-check: the number of monic irreducibles of each degree
//     over F_p must equal Gauss's necklace formula (1/n) sum_{d|n} mu(d) p^{n/d}.
//   * Factorization round-trips: re-multiplying the returned factors must
//     recover the input, and every returned factor must be irreducible.
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <iostream>
#include <vector>

using namespace alcp;

// ----------------------------------------------------------------------------
// number theory helpers used only by the tests
// ----------------------------------------------------------------------------
static int mobius(int n) {
    if (n == 1) return 1;
    int prod = 1;
    for (int p = 2; p * p <= n; ++p) {
        if (n % p == 0) {
            n /= p;
            if (n % p == 0) return 0;  // squared factor
            prod = -prod;
        }
    }
    if (n > 1) prod = -prod;  // remaining prime factor
    return prod;
}

// Gauss's formula for the count of monic irreducibles of degree deg over F_q.
static BigInt num_monic_irreducible(const BigInt& q, int deg) {
    BigInt acc = 0;
    for (int d = 1; d <= deg; ++d) {
        if (deg % d == 0) acc += mobius(d) * ipow(q, deg / d);
    }
    return acc / deg;
}

// ----------------------------------------------------------------------------
// 1. integer helpers
// ----------------------------------------------------------------------------
static void test_integers() {
    for (int n : {2, 3, 5, 7, 11, 13, 9973})
        CHECK(is_prime(BigInt(n)));
    for (int n : {0, 1, 4, 6, 9, 15, 21, 100})
        CHECK(!is_prime(BigInt(n)));
    CHECK(is_prime(BigInt("1000000000000000003")));

    auto pp1 = prime_power_factor(BigInt(27));
    CHECK(pp1.first == 3 && pp1.second == 3);
    auto pp2 = prime_power_factor(BigInt(1024));
    CHECK(pp2.first == 2 && pp2.second == 10);
    auto pp3 = prime_power_factor(BigInt(13));
    CHECK(pp3.first == 13 && pp3.second == 1);
    bool threw = false;
    try { prime_power_factor(BigInt(12)); } catch (...) { threw = true; }
    CHECK(threw);  // 12 is not a prime power
}

// ----------------------------------------------------------------------------
// 2. F_p axioms, exhaustive over all elements / pairs / triples
// ----------------------------------------------------------------------------
static void test_fp_axioms() {
    for (int pp : {2, 3, 5, 7, 13}) {
        Fp F(pp);
        BigInt p = F.characteristic();
        for (long i = 0; i < pp; ++i) {
            auto a = F.from_int(BigInt(i));
            CHECK(F.add(a, F.zero()) == a);
            CHECK(F.mul(a, F.one()) == a);
            CHECK(F.add(a, F.neg(a)) == F.zero());
            if (!F.is_zero(a)) {
                CHECK(F.mul(a, F.inv(a)) == F.one());
                CHECK(F.pow(a, p - 1) == F.one());           // Fermat
                CHECK(F.pow(a, p) == a);                     // x^p = x
            }
            for (long j = 0; j < pp; ++j) {
                auto b = F.from_int(BigInt(j));
                CHECK(F.add(a, b) == F.add(b, a));           // commutative +
                CHECK(F.mul(a, b) == F.mul(b, a));           // commutative *
                CHECK(F.sub(F.add(a, b), b) == a);           // (a+b)-b = a
                for (long k = 0; k < pp; ++k) {
                    auto c = F.from_int(BigInt(k));
                    CHECK(F.add(F.add(a, b), c) == F.add(a, F.add(b, c)));  // assoc +
                    CHECK(F.mul(F.mul(a, b), c) == F.mul(a, F.mul(b, c)));  // assoc *
                    CHECK(F.mul(a, F.add(b, c)) ==
                          F.add(F.mul(a, b), F.mul(a, c)));                 // distributive
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------
// 3. F_q axioms, exhaustive over all elements / pairs (triples for tiny q)
// ----------------------------------------------------------------------------
static void test_fq_axioms() {
    struct Spec { int p; unsigned long n; };
    for (Spec s : {Spec{2, 2}, Spec{2, 3}, Spec{3, 2}, Spec{5, 2}, Spec{3, 3}, Spec{7, 2}}) {
        Fq F(BigInt(s.p), s.n);
        long q = static_cast<long>(F.cardinality().get_si());
        BigInt p = F.characteristic();
        for (long i = 0; i < q; ++i) {
            auto a = F.from_index(BigInt(i));
            CHECK(F.add(a, F.zero()) == a);
            CHECK(F.mul(a, F.one()) == a);
            CHECK(F.add(a, F.neg(a)) == F.zero());
            CHECK(F.to_index(a) == i);                       // index <-> elem bijection
            if (!F.is_zero(a)) {
                CHECK(F.mul(a, F.inv(a)) == F.one());
                CHECK(F.pow(a, F.cardinality() - 1) == F.one());
                CHECK(F.pth_root(F.pow(a, p)) == a);         // p-th root undoes Frobenius
            }
            for (long j = 0; j < q; ++j) {
                auto b = F.from_index(BigInt(j));
                CHECK(F.mul(a, b) == F.mul(b, a));
                // Frobenius automorphism: (a+b)^p = a^p + b^p
                CHECK(F.pow(F.add(a, b), p) ==
                      F.add(F.pow(a, p), F.pow(b, p)));
                if (q <= 27) {
                    for (long k = 0; k < q; ++k) {
                        auto c = F.from_index(BigInt(k));
                        CHECK(F.mul(F.mul(a, b), c) == F.mul(a, F.mul(b, c)));
                        CHECK(F.mul(a, F.add(b, c)) ==
                              F.add(F.mul(a, b), F.mul(a, c)));
                    }
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------
// 4. polynomial ring laws + divmod identity (randomized)
// ----------------------------------------------------------------------------
template <class Field>
static void test_poly_ring(const Field& F, int rounds) {
    using P = Poly<Field>;
    for (int t = 0; t < rounds; ++t) {
        P a = P::random(F, rand() % 8);
        P b = P::random(F, rand() % 8);
        P c = P::random(F, rand() % 8);
        CHECK(a + b == b + a);
        CHECK(a * b == b * a);
        CHECK(a * (b + c) == a * b + a * c);
        CHECK((a + b) - b == a);
        if (!b.is_zero()) {
            auto [q, r] = a.divmod(b);
            CHECK(a == q * b + r);
            CHECK(r.is_zero() || r.degree() < b.degree());
        }
        if (!a.is_zero() && !b.is_zero())
            CHECK((a * b).degree() == a.degree() + b.degree());
    }
}

// ----------------------------------------------------------------------------
// 5. gcd, extended gcd (Bezout), modular inverse
// ----------------------------------------------------------------------------
template <class Field>
static void test_gcd(const Field& F, int rounds) {
    using P = Poly<Field>;
    for (int t = 0; t < rounds; ++t) {
        P a = P::random(F, 1 + rand() % 6);
        P b = P::random(F, 1 + rand() % 6);
        auto e = P::gcd_ext(a, b);
        CHECK(e.g == e.x * a + e.y * b);                 // Bezout identity
        if (!e.g.is_zero()) {
            CHECK((a % e.g).is_zero() && (b % e.g).is_zero());  // g divides a, b
            CHECK(F.is_one(e.g.lead()));                        // g monic
        }
        // modular inverse against an irreducible modulus
        P m = P::random(F, 3);
        if (m.is_irreducible() && !(a % m).is_zero()) {
            P ai = P::inv_mod(a, m);
            CHECK(((a * ai) % m).is_one());
        }
    }
}

// ----------------------------------------------------------------------------
// 6. irreducible counts over F_p match Gauss's formula (exhaustive enumeration)
// ----------------------------------------------------------------------------
static void test_irreducible_counts() {
    struct Spec { int p, maxdeg; };
    for (Spec s : {Spec{2, 7}, Spec{3, 4}, Spec{5, 3}, Spec{7, 3}}) {
        Fp F(s.p);
        for (int deg = 1; deg <= s.maxdeg; ++deg) {
            BigInt total = ipow(BigInt(s.p), deg);  // number of monic degree-deg polys
            long count = 0;
            for (BigInt idx = 0; idx < total; ++idx) {
                std::vector<BigInt> coeffs(deg + 1);
                BigInt t = idx;
                for (int i = 0; i < deg; ++i) { coeffs[i] = t % s.p; t /= s.p; }
                coeffs[deg] = 1;  // monic
                if (Poly<Fp>::from_ints(F, coeffs).is_irreducible()) ++count;
            }
            CHECK(BigInt(count) == num_monic_irreducible(BigInt(s.p), deg));
        }
    }
}

// ----------------------------------------------------------------------------
// 7. factorization round-trips
// ----------------------------------------------------------------------------
template <class Field>
static void test_factor(const Field& F, int rounds) {
    using P = Poly<Field>;
    for (int t = 0; t < rounds; ++t) {
        // assemble a known product of small random factors with multiplicities
        P comp = P::constant(F, F.random());
        while (comp.is_zero()) comp = P::constant(F, F.random());
        for (int i = 0; i < 3; ++i) {
            int d = 1 + rand() % 3;
            int e = 1 + rand() % 3;
            comp = comp * P::random(F, d).monic().pow(BigInt(e));
        }
        if (comp.degree() <= 0) continue;
        auto fac = comp.factor();
        // every factor irreducible
        for (auto& [u, e] : fac.factors) {
            CHECK(u.is_irreducible());
            CHECK(F.is_one(u.lead()));  // monic
        }
        // distinct factors
        for (std::size_t i = 0; i < fac.factors.size(); ++i)
            for (std::size_t j = i + 1; j < fac.factors.size(); ++j)
                CHECK(fac.factors[i].first != fac.factors[j].first);
        // re-multiply to recover comp
        P recon = P::constant(F, fac.leading);
        for (auto& [u, e] : fac.factors) recon = recon * u.pow(BigInt(e));
        CHECK(recon == comp);
    }
}

// ----------------------------------------------------------------------------
// 8. multiplication algorithms agree (schoolbook == Karatsuba == FFT == auto)
// ----------------------------------------------------------------------------
template <class Field>
static void test_mul_algos(const Field& F, int da, int db) {
    using P = Poly<Field>;
    P a = P::random(F, da), b = P::random(F, db);
    P school = a.mul(b, MulAlgo::Schoolbook);
    P kara = a.mul(b, MulAlgo::Karatsuba);
    P autom = a * b;
    CHECK(school == kara);
    CHECK(school == autom);
    if (fft::fft_exponent(F, a.coeffs().size() + b.coeffs().size() - 1) >= 0) {
        P viafft = a.mul(b, MulAlgo::FFT);
        CHECK(school == viafft);
    }
}

// ----------------------------------------------------------------------------
// 9. Toeplitz matrix-vector products and triangular inverse
// ----------------------------------------------------------------------------
template <class Field>
static void test_toeplitz(const Field& F, std::size_t n) {
    using Elem = typename Field::Elem;
    auto rvec = [&](std::size_t len) {
        std::vector<Elem> v(len);
        for (auto& x : v) x = F.random();
        return v;
    };

    // lower-triangular: T[i][j] = a[i-j] for i>=j
    {
        auto a = rvec(n), b = rvec(n);
        auto got = fast::toeplitz_lower_vec(F, n, a, b);
        for (std::size_t i = 0; i < n; ++i) {
            Elem acc = F.zero();
            for (std::size_t j = 0; j <= i; ++j) acc = F.add(acc, F.mul(a[i - j], b[j]));
            CHECK(got[i] == acc);
        }
    }
    // upper-triangular: T[i][j] = a[j-i] for j>=i
    {
        auto a = rvec(n), b = rvec(n);
        auto got = fast::toeplitz_upper_vec(F, n, a, b);
        for (std::size_t i = 0; i < n; ++i) {
            Elem acc = F.zero();
            for (std::size_t j = i; j < n; ++j) acc = F.add(acc, F.mul(a[j - i], b[j]));
            CHECK(got[i] == acc);
        }
    }
    // full Toeplitz with the 2n-1 layout (first row, then first column minus corner)
    {
        auto a = rvec(2 * n - 1), b = rvec(n);
        // diagonal values t_{i-j}: t_0=a[0]; t_{-s}=a[s]; t_{+s}=a[n-1+s]
        auto tval = [&](long diff) -> Elem {
            if (diff <= 0) return a[-diff];
            return a[n - 1 + diff];
        };
        auto got = fast::toeplitz_vec(F, n, a, b);
        for (std::size_t i = 0; i < n; ++i) {
            Elem acc = F.zero();
            for (std::size_t j = 0; j < n; ++j)
                acc = F.add(acc, F.mul(tval((long)i - (long)j), b[j]));
            CHECK(got[i] == acc);
        }
    }
    // lower-triangular inverse: convolving a with its inverse column gives e_0
    {
        std::vector<Elem> a = rvec(n);
        while (F.is_zero(a[0])) a[0] = F.random();
        auto x = fast::toeplitz_lower_inverse(F, n, a);
        auto prod = fast::toeplitz_lower_vec(F, n, a, x);
        CHECK(prod[0] == F.one());
        for (std::size_t i = 1; i < n; ++i) CHECK(F.is_zero(prod[i]));
    }
}

// ----------------------------------------------------------------------------
// 10. explicit DFT / IDFT round-trip and correctness
// ----------------------------------------------------------------------------
template <class Field>
static void test_dft(const Field& F, unsigned long k) {
    using Elem = typename Field::Elem;
    unsigned long cap = F.fft_capacity();
    CHECK(k <= cap);
    Elem g = F.pow(F.fft_base_root(), BigInt(1) << (cap - k));  // root of order 2^k
    std::size_t n = std::size_t(1) << k;

    std::vector<Elem> a(n);
    for (auto& x : a) x = F.random();

    auto A = fast::dft(F, g, k, a);
    // naive DFT: A[i] = sum_j a[j] g^{ij}
    for (std::size_t i = 0; i < n; ++i) {
        Elem acc = F.zero();
        Elem gi = F.pow(g, BigInt(static_cast<unsigned long>(i)));
        Elem w = F.one();
        for (std::size_t j = 0; j < n; ++j) { acc = F.add(acc, F.mul(a[j], w)); w = F.mul(w, gi); }
        CHECK(A[i] == acc);
    }
    auto back = fast::idft(F, g, k, A);
    CHECK(back == a);  // inverse recovers the input
}

int main() {
    set_seed(20260629);
    srand(12345);

    std::cout << "running alcp verification suite...\n";

    test_integers();
    test_fp_axioms();
    test_fq_axioms();

    Fp F7(7), F101(101), Fbig(998244353);  // last one is FFT-friendly (2^23 | p-1)
    test_poly_ring(F7, 200);
    test_poly_ring(F101, 200);
    test_gcd(F7, 200);
    test_gcd(F101, 100);
    test_factor(F7, 60);
    test_factor(F101, 40);

    Fq F9(BigInt(3), 2), F27(BigInt(3), 3), F25(BigInt(5), 2);
    test_poly_ring(F9, 100);
    test_gcd(F9, 80);
    test_factor(F9, 30);
    test_factor(F27, 20);
    test_factor(F25, 20);

    test_irreducible_counts();

    // multiplication algorithm equivalence, incl. an FFT-capable field
    test_mul_algos(F7, 40, 35);
    test_mul_algos(F9, 40, 30);
    test_mul_algos(Fbig, 600, 500);   // exercises the FFT path
    test_mul_algos(Fbig, 50, 3);

    // Toeplitz over both kinds of field
    test_toeplitz(F7, 9);
    test_toeplitz(F101, 12);
    test_toeplitz(F9, 7);

    // DFT/IDFT: F_17 (2^4 | 16) and F_9 (2^3 | 8)
    Fp F17(17);
    test_dft(F17, 4);
    test_dft(Fbig, 10);
    test_dft(F9, 3);

    return check::summary("alcp");
}
