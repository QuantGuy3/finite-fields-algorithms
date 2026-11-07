// calibra_hgcd.cpp — gcd / gcd_ext con half-GCD frente al Euclides clásico (fase 4, nivel B).
//
// Pares aleatorios de grados n y n - 1 (sucesión de restos normal), rotando
// entre 8 pares distintos.  Las tablas FFT se construyen antes.
//
//   g++ -std=c++23 -O3 -Iinclude bench/calibra_hgcd.cpp -lgmpxx -lgmp -static -o build/calibra_hgcd.exe
#include "alcp/alcp.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <vector>

using namespace alcp;
using reloj = std::chrono::steady_clock;

template <class Fn>
double medir_us(Fn&& fn) {
    double mejor = 1e300;
    for (int tanda = 0; tanda < 3; ++tanda) {
        int veces = 0;
        auto t0 = reloj::now();
        double dt = 0;
        do {
            fn(veces % 8);
            ++veces;
            dt = std::chrono::duration<double, std::micro>(reloj::now() - t0).count();
        } while (dt < 30000);
        mejor = std::min(mejor, dt / veces);
    }
    return mejor;
}

template <class Field>
void calibrar(const Field& K, const char* nombre, int nmax) {
    using P = Poly<Field>;
    std::printf("\n%s\n%7s %12s %12s %8s %14s %14s %8s\n", nombre, "n", "gcd clás", "gcd hgcd", "x", "gcd_ext clás",
                "gcd_ext hgcd", "x");
    (void)(P::random(K, 300) * P::random(K, 300));
    for (int n : {128, 192, 256, 384, 512, 1024, 2048, 4096}) {
        if (n > nmax) break;
        std::vector<P> as, bs;
        for (int i = 0; i < 8; ++i) {
            as.push_back(P::random(K, n));
            bs.push_back(P::random(K, n - 1));
        }
        const double c = medir_us([&](int i) { (void)P::gcd_classic(as[i], bs[i]); });
        const double h = medir_us([&](int i) { (void)P::gcd(as[i], bs[i]); });
        const double ce = medir_us([&](int i) { (void)P::gcd_ext_classic(as[i], bs[i]); });
        const double he = medir_us([&](int i) { (void)P::gcd_ext(as[i], bs[i]); });
        std::printf("%7d %12.1f %12.1f %8.2f %14.1f %14.1f %8.2f\n", n, c, h, c / h, ce, he, ce / he);
    }
}

int main() {
    set_seed(1);
    const std::vector<BigInt> mod_3_30{2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0,
                                       2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1};
    Fp64 w101(101), wntt(998244353), wm61((BigInt(1) << 61) - 1);
    Fp g127((BigInt(1) << 127) - 1);
    Fq64 w9(Fp64(BigInt(3)), std::vector<BigInt>{1, 0, 1});
    Fq64 w256(Fp64(BigInt(2)), std::vector<BigInt>{1, 1, 0, 1, 1, 0, 0, 0, 1});
    Fq64 w3e30(Fp64(BigInt(3)), mod_3_30);
    calibrar(wntt, "Fp64 p = 998244353", 4096);
    calibrar(w101, "Fp64 p = 101", 4096);
    calibrar(wm61, "Fp64 p = 2^61 - 1", 4096);
    calibrar(g127, "Fp GMP p = 2^127 - 1", 2048);
    calibrar(w9, "Fq64 F_9", 2048);
    calibrar(w256, "Fq64 F_256", 1024);
    calibrar(w3e30, "Fq64 F_(3^30)", 512);
    return 0;
}
