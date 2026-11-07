// Coste aislado de la NTT en palabras (Fp64, p = 998244353).
//   g++ -std=c++23 -O3 -march=native -Iinclude bench/experimentos/perfil_ntt.cpp -lgmpxx -lgmp -static -o build/perfil_ntt.exe
//   (investigación de la fase 5: docs/optimizacion/investigacion-fase5.md)
#include "alcp/alcp.hpp"
#include <chrono>
#include <cstdio>
#include <vector>
using namespace alcp;
using reloj = std::chrono::steady_clock;
int main() {
    Fp64 K(998244353); using P = Poly<Fp64>;
    (void)(P::random(K, 300) * P::random(K, 300));
    for (unsigned k : {10u, 11u, 12u, 13u, 16u}) {
        const auto& t = fft::tables(K, k);
        const std::size_t n = std::size_t(1) << k;
        std::vector<std::uint64_t> a(n);
        for (auto& x : a) x = K.random();
        int v = 0; auto t0 = reloj::now(); double dt = 0;
        do { fft::word::transform(t, a.data(), n, false); ++v; dt = std::chrono::duration<double, std::micro>(reloj::now() - t0).count(); } while (dt < 100000);
        const double us = dt / v;
        std::printf("n = 2^%u: %8.1f us por transformada, %.2f ns por mariposa\n", k, us, us * 1e3 / (n / 2 * k));
    }
}
