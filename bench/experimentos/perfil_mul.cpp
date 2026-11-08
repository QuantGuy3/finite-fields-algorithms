// perfil_mul.cpp — Producto por algoritmo en Fp64 (auto, NTT, Kronecker, Karatsuba) y resto a % b.
//   g++ -std=c++23 -O3 -march=native -Iinclude bench/experimentos/perfil_mul.cpp -lgmpxx -lgmp -static -o build/perfil_mul.exe
//   (investigación de la fase 5: docs/optimizacion/investigacion-fase5.md)
#include "alcp/alcp.hpp"
#include <chrono>
#include <cstdio>
#include <vector>
#include <string>
using namespace alcp;
using reloj = std::chrono::steady_clock;
template <class Fn> double medir_us(Fn&& fn) {
    double mejor = 1e300;
    for (int t = 0; t < 3; ++t) { int v = 0; auto t0 = reloj::now(); double dt = 0;
        do { fn(v % 8); ++v; dt = std::chrono::duration<double, std::micro>(reloj::now() - t0).count(); } while (dt < 30000);
        mejor = std::min(mejor, dt / v); }
    return mejor;
}
int main() {
    for (std::uint64_t p : {998244353ull, 101ull, (1ull << 61) - 1}) {
        Fp64 K(BigInt(std::to_string(p))); using P = Poly<Fp64>;
        (void)(P::random(K, 300) * P::random(K, 300));
        std::printf("p=%llu\n%6s %9s %9s %9s %9s %9s\n", (unsigned long long)p, "n", "auto", "fft", "kron", "karat", "gcd%");
        for (int n : {100, 300, 1000, 2000, 8000}) {
            std::vector<P> a, b; for (int i = 0; i < 8; ++i) { a.push_back(P::random(K, n)); b.push_back(P::random(K, n)); }
            double ta = medir_us([&](int i){ (void)a[i].mul(b[i], MulAlgo::Auto); });
            double tf = fft::fft_exponent(K, 2*n+2) >= 0 ? medir_us([&](int i){ (void)a[i].mul(b[i], MulAlgo::FFT); }) : 0;
            double tk = medir_us([&](int i){ (void)a[i].mul(b[i], MulAlgo::Kronecker); });
            double tr = n <= 2000 ? medir_us([&](int i){ (void)a[i].mul(b[i], MulAlgo::Karatsuba); }) : 0;
            double tm = n <= 1000 ? medir_us([&](int i){ (void)(a[i] % b[(i+1)%8].monic()); }) : 0;
            std::printf("%6d %9.1f %9.1f %9.1f %9.1f %9.1f\n", n, ta, tf, tk, tr, tm);
        }
    }
}
