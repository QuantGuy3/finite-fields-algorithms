// Minimal smoke test: F_p and F_p[x] basics, to catch compile/logic errors
// before the full suite.
//
// Usa el arnés CHECK de check.hpp (no assert): las comprobaciones se ejecutan
// también en Release (-DNDEBUG) y el código de salida refleja los fallos.
#include "alcp/bigint.hpp"
#include "alcp/fp.hpp"
#include "alcp/polynomial.hpp"
#include "check.hpp"

#include <iostream>

using namespace alcp;

int main() {
    set_seed(12345);

    // ---- integer helpers ----
    CHECK(is_prime(BigInt(7)));
    CHECK(!is_prime(BigInt(9)));
    auto [p, n] = prime_power_factor(BigInt(27));
    CHECK(p == 3 && n == 3);

    // ---- F_p ----
    Fp F(7);
    CHECK(F.add(F.from_int(5), F.from_int(4)) == 2);          // 9 mod 7
    CHECK(F.mul(F.from_int(3), F.from_int(5)) == 1);          // 15 mod 7
    CHECK(F.mul(F.from_int(3), F.inv(F.from_int(3))) == 1);   // a * a^-1
    CHECK(F.pow(F.from_int(3), BigInt(6)) == 1);              // Fermat

    using P = Poly<Fp>;
    // f = x^2 + 1, g = x + 2  over F_7
    P f = P::from_ints(F, {1, 0, 1});
    P g = P::from_ints(F, {2, 1});
    P prod = f * g;
    auto [q, r] = prod.divmod(g);
    CHECK(q == f && r.is_zero());                              // (f*g)/g == f

    // gcd(x^2-1, x-1) = x-1 (monic)
    P a = P::from_ints(F, {6, 0, 1});  // x^2 - 1
    P b = P::from_ints(F, {6, 1});     // x - 1
    P gg = P::gcd(a, b);
    CHECK(gg == b.monic());

    // irreducibility: x^2+1 over F_7 (7 = 3 mod 4 -> -1 is a non-residue -> irreducible)
    CHECK(f.is_irreducible());
    // x^2 - 1 = (x-1)(x+1) reducible
    CHECK(!a.is_irreducible());

    // factor (x-1)^2 (x^2+1) over F_7
    P composite = b * b * f;
    auto fac = composite.factor();
    BigInt total_deg = 0;
    for (auto& [u, e] : fac.factors) {
        CHECK(u.is_irreducible());
        total_deg += u.degree() * e;
    }
    CHECK(total_deg == composite.degree());

    // reconstruct: leading * prod(factors^mult) == composite
    P recon = P::constant(F, fac.leading);
    for (auto& [u, e] : fac.factors) recon = recon * u.pow(BigInt(e));
    CHECK(recon == composite);

    std::cout << "smoke: composite = " << composite.to_string()
              << "  ->  " << fac.factors.size() << " irreducible factor(s)\n";
    for (auto& [u, e] : fac.factors)
        std::cout << "   (" << u.to_string() << ")^" << e << "\n";
    return check::summary("smoke");
}
