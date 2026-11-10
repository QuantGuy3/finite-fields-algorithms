// calibra_ntt_kronecker.cpp — Fase 5: frontera NTT en palabras / Kronecker /
// Karatsuba en Fp64 con NTT, tras las mariposas perezosas.  Entradas variadas
// (16 pares).  Ejecutar también con ALCP_THREADS=1.
//
//   g++ -std=c++23 -O3 -Iinclude bench/calibra_ntt_kronecker.cpp -lgmpxx -lgmp -static -o build/calibra_ntt_kronecker.exe
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
    for (int tanda = 0; tanda < 5; ++tanda) {
        int veces = 0;
        auto t0 = reloj::now();
        double dt = 0;
        do {
            fn(veces % 16);
            ++veces;
            dt = std::chrono::duration<double, std::micro>(reloj::now() - t0).count();
        } while (dt < 20000);
        mejor = std::min(mejor, dt / veces);
    }
    return mejor;
}

int main() {
    set_seed(1);
    for (const char* ps : {"998244353", "2305843009213693951"}) {
        Fp64 K{BigInt(ps)};
        using P = Poly<Fp64>;
        (void)(P::random(K, 300) * P::random(K, 300));
        std::printf("\np = %s\n%6s %10s %10s %10s %10s\n", ps, "n", "escolar", "Karatsuba", "NTT", "Kronecker");
        for (int n : {24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 4096}) {
            std::vector<P> a, b;
            for (int i = 0; i < 16; ++i) { a.push_back(P::random(K, n - 1)); b.push_back(P::random(K, n - 1)); }
            const bool fft = fft::fft_exponent(K, 2 * n) >= 0;
            const double ts = n <= 256 ? medir_us([&](int i) { (void)a[i].mul(b[i], MulAlgo::Schoolbook); }) : 0;
            const double tk = medir_us([&](int i) { (void)a[i].mul(b[i], MulAlgo::Karatsuba); });
            const double tf = fft ? medir_us([&](int i) { (void)a[i].mul(b[i], MulAlgo::FFT); }) : 0;
            const double tr = medir_us([&](int i) { (void)a[i].mul(b[i], MulAlgo::Kronecker); });
            std::printf("%6d %10.1f %10.1f %10.1f %10.1f\n", n, ts, tk, tf, tr);
        }
    }
}
