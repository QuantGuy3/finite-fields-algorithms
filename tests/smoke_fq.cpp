// Smoke test for F_q and F_q[x]: arithmetic, inverse, Frobenius, irreducibility
// and factorization, all via the same generic Poly<Field> code.
#include "alcp/bigint.hpp"
#include "alcp/fq.hpp"
#include "alcp/polynomial.hpp"

#include <cassert>
#include <iostream>

using namespace alcp;

int main() {
    set_seed(2024);

    // F_9 = F_3[x]/(x^2 + 1)  (x^2+1 is irreducible over F_3)
    Fp F3(3);
    Fq F9(F3, std::vector<BigInt>{1, 0, 1});
    assert(F9.cardinality() == 9);
    assert(F9.characteristic() == 3);

    // field axioms by exhaustion over all 9 elements
    for (long i = 0; i < 9; ++i) {
        auto a = F9.from_index(i);
        assert(F9.add(a, F9.neg(a)) == F9.zero());
        assert(F9.add(a, F9.zero()) == a);
        assert(F9.mul(a, F9.one()) == a);
        if (!F9.is_zero(a)) {
            assert(F9.mul(a, F9.inv(a)) == F9.one());
            assert(F9.pow(a, F9.cardinality() - 1) == F9.one());  // a^(q-1)=1
        }
        // Frobenius is a field automorphism: (a+b)^p = a^p + b^p
        for (long j = 0; j < 9; ++j) {
            auto b = F9.from_index(j);
            auto lhs = F9.pow(F9.add(a, b), F9.characteristic());
            auto rhs = F9.add(F9.pow(a, F9.characteristic()), F9.pow(b, F9.characteristic()));
            assert(lhs == rhs);
        }
        // p-th root undoes Frobenius
        assert(F9.pth_root(F9.pow(a, 3)) == a);
    }

    // F_q[x] over F_9
    using PQ = Poly<Fq>;
    PQ x = PQ::variable(F9);
    PQ f = (x * x) + PQ::constant(F9, F9.from_index(2));  // x^2 + 2

    // divmod consistency on random polynomials
    for (int t = 0; t < 50; ++t) {
        PQ A = PQ::random(F9, 6);
        PQ B = PQ::random(F9, 3);
        auto [q, r] = A.divmod(B);
        assert(A == q * B + r);
        assert(r.is_zero() || r.degree() < B.degree());
    }

    // factorization round-trip over F_9
    PQ p1 = PQ::random(F9, 2);
    PQ p2 = PQ::random(F9, 3);
    PQ comp = p1 * p1 * p2;
    auto fac = comp.factor();
    PQ recon = PQ::constant(F9, fac.leading);
    for (auto& [u, e] : fac.factors) {
        assert(u.is_irreducible());
        recon = recon * u.pow(BigInt(e));
    }
    assert(recon == comp);

    // also test auto-built F_q from a prime power
    Fq F27(BigInt(27));
    assert(F27.cardinality() == 27);
    auto z = F27.random();
    if (!F27.is_zero(z)) assert(F27.mul(z, F27.inv(z)) == F27.one());

    std::cout << "smoke_fq OK: F_9 modulus g, factored deg-" << comp.degree()
              << " poly into " << fac.factors.size() << " irreducibles\n";
    return 0;
}
