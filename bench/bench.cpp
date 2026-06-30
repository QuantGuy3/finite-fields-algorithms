// bench.cpp — timing for the scenarios called out in the original C++ notes
// ("17s to factor a random degree-5 polynomial over F_q with q = 3^30").
#include "alcp/alcp.hpp"

#include <chrono>
#include <cstdio>

using namespace alcp;
using clk = std::chrono::steady_clock;

template <class F>
static double ms(F&& f) {
    auto t0 = clk::now();
    f();
    return std::chrono::duration<double, std::milli>(clk::now() - t0).count();
}

int main() {
    set_seed(1);

    // The headline case: F_q with q = 3^30, factor a random degree-5 poly.
    printf("building F_q, q = 3^30 (random irreducible modulus of degree 30)...\n");
    double t_build = ms([] { volatile Fq f(BigInt(3), 30); (void)f; });
    Fq Fq330(BigInt(3), 30);
    printf("  field built in %.1f ms  (q = %s)\n", t_build, Fq330.cardinality().get_str().c_str());

    using PQ = Poly<Fq>;
    PQ f = PQ::random(Fq330, 5);
    double t_fac = ms([&] {
        auto fac = f.factor();
        long parts = 0;
        for (auto& pr : fac.factors) parts += pr.second * pr.first.degree();
        if (parts != f.degree()) printf("  !! degree mismatch\n");
    });
    printf("  factor random deg-5 poly over F_(3^30): %.1f ms\n\n", t_fac);

    // FFT vs schoolbook on a big multiplication over an NTT-friendly prime.
    Fp P(998244353);
    using PF = Poly<Fp>;
    PF a = PF::random(P, 4000), b = PF::random(P, 4000);
    double t_school = ms([&] { volatile auto r = a.mul(b, MulAlgo::Schoolbook); (void)r; });
    double t_kara = ms([&] { volatile auto r = a.mul(b, MulAlgo::Karatsuba); (void)r; });
    double t_fft = ms([&] { volatile auto r = a.mul(b, MulAlgo::FFT); (void)r; });
    printf("multiply two degree-4000 polys over F_998244353:\n");
    printf("  schoolbook %.1f ms | karatsuba %.1f ms | FFT %.1f ms  (auto picks FFT)\n",
           t_school, t_kara, t_fft);

    // Factor a larger polynomial over F_7.
    Fp F7(7);
    PF big = PF::one(F7);
    for (int i = 0; i < 6; ++i) big = big * PF::random(F7, 8).monic();
    double t_big = ms([&] { volatile auto r = big.factor(); (void)r; });
    printf("\nfactor a degree-%d poly over F_7: %.1f ms\n", big.degree(), t_big);
    return 0;
}
