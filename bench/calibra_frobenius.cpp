// calibra_frobenius.cpp — Frobenius por potencias frente a composición modular (fase 4, nivel B).
//
// is_irreducible de un polinomio IRREDUCIBLE de grado n (el caso caro: el
// bucle de Ben-Or llega hasta n/2) con FrobeniusAlgo::Potencia (pow_mod),
// Composicion (Brent-Kung desde el primer paso) y Auto (alquilar o comprar).
//
//   g++ -std=c++23 -O3 -Iinclude bench/calibra_frobenius.cpp -lgmpxx -lgmp -static -o build/calibra_frobenius.exe
#include "alcp/alcp.hpp"

#include <chrono>
#include <cstdio>
#include <vector>

using namespace alcp;
using reloj = std::chrono::steady_clock;

template <class Fn>
double medir_ms(Fn&& fn) {
    double mejor = 1e300;
    for (int t = 0; t < 2; ++t) {
        auto t0 = reloj::now();
        fn();
        mejor = std::min(mejor, std::chrono::duration<double, std::milli>(reloj::now() - t0).count());
    }
    return mejor;
}

template <class Field>
void calibrar(const Field& K, const char* nombre, const std::vector<int>& grados) {
    using P = Poly<Field>;
    std::printf("\n%s\n%6s %12s %12s %12s %8s %8s\n", nombre, "n", "Potencia ms", "Compos. ms", "Auto ms", "P/C",
                "P/Auto");
    (void)(P::random(K, 300) * P::random(K, 300));
    for (int n : grados) {
        P f = P::random(K, n).monic();
        while (!f.is_irreducible()) f = P::random(K, n).monic();
        const double tp = medir_ms([&] { (void)f.is_irreducible(FrobeniusAlgo::Potencia); });
        const double tc = medir_ms([&] { (void)f.is_irreducible(FrobeniusAlgo::Composicion); });
        const double ta = medir_ms([&] { (void)f.is_irreducible(); });
        std::printf("%6d %12.2f %12.2f %12.2f %8.2f %8.2f\n", n, tp, tc, ta, tp / tc, tp / ta);
    }
}

int main() {
    set_seed(1);
    const std::vector<BigInt> mod_3_30{2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0,
                                       2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1};
    Fp64 w2(2), w101(101), wntt(998244353), wm61((BigInt(1) << 61) - 1);
    Fp g127((BigInt(1) << 127) - 1);
    Fq64 w9(Fp64(BigInt(3)), std::vector<BigInt>{1, 0, 1});
    Fq64 w256(Fp64(BigInt(2)), std::vector<BigInt>{1, 1, 0, 1, 1, 0, 0, 0, 1});
    Fq64 w3e30(Fp64(BigInt(3)), mod_3_30);
    calibrar(w2, "F_2 (q = 2)", {64, 128, 256});
    calibrar(w101, "F_101", {32, 64, 128, 200, 400});
    calibrar(wntt, "F_998244353", {32, 64, 128, 256, 400});
    calibrar(wm61, "F_(2^61-1)", {16, 32, 64, 128});
    calibrar(g127, "F_(2^127-1) (GMP)", {8, 16, 32, 64});
    calibrar(w9, "F_9", {32, 64, 128});
    calibrar(w256, "F_256", {16, 32, 64});
    calibrar(w3e30, "F_(3^30)", {6, 10, 16, 24});
    return 0;
}
