// Experimento: mariposa actual frente a mariposa sin saltos y mariposa perezosa
// de Harvey (valores en [0, 4p)).  Se comprueba que las tres dan lo mismo.
//   g++ -std=c++23 -O3 -march=native -Iinclude bench/experimentos/ntt_sin_saltos.cpp -lgmpxx -lgmp -static -o build/ntt_sin_saltos.exe
//   (investigación de la fase 5: docs/optimizacion/investigacion-fase5.md)
#include "alcp/alcp.hpp"
#include <chrono>
#include <cstdio>
#include <vector>
using namespace alcp;
using u64 = std::uint64_t; using u128 = unsigned __int128;
using reloj = std::chrono::steady_clock;

// Variante 1: sin saltos (máscaras en lugar de comparaciones).
template <class T> void sin_saltos(const T& t, u64* a, std::size_t n) {
    const auto& W = t.word; const u64 p = W.p;
    for (const auto& [i, j] : t.swaps) std::swap(a[i], a[j]);
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len >> 1;
        const u64* tw = W.fwd.data() + (half - 1); const u64* ts = W.fwd_sh.data() + (half - 1);
        for (std::size_t i = 0; i < n; i += len) {
            u64* x = a + i; u64* y = a + i + half;
            for (std::size_t j = 0; j < half; ++j) {
                const u64 q = static_cast<u64>((static_cast<u128>(ts[j]) * y[j]) >> 64);
                u64 v = tw[j] * y[j] - q * p;
                v -= p & (0 - static_cast<u64>(v >= p));
                const u64 u = x[j];
                u64 s = u + v - p; s += p & (0 - (s >> 63));
                u64 d = u - v;     d += p & (0 - (d >> 63));
                x[j] = s; y[j] = d;
            }
        }
    }
}
// Variante 2: perezosa (Harvey 2014): entradas y salidas en [0, 4p), una sola
// corrección por mariposa; normalización al final.  Requiere 4p < 2^64.
template <class T> void perezosa(const T& t, u64* a, std::size_t n) {
    const auto& W = t.word; const u64 p = W.p, p2 = 2 * p;
    for (const auto& [i, j] : t.swaps) std::swap(a[i], a[j]);
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len >> 1;
        const u64* tw = W.fwd.data() + (half - 1); const u64* ts = W.fwd_sh.data() + (half - 1);
        for (std::size_t i = 0; i < n; i += len) {
            u64* x = a + i; u64* y = a + i + half;
            for (std::size_t j = 0; j < half; ++j) {
                u64 u = x[j];
                u -= p2 & (0 - static_cast<u64>(u >= p2));               // u en [0, 2p)
                const u64 q = static_cast<u64>((static_cast<u128>(ts[j]) * y[j]) >> 64);
                const u64 v = tw[j] * y[j] - q * p;                       // v en [0, 2p)
                x[j] = u + v; y[j] = u - v + p2;                          // en [0, 4p)
            }
        }
    }
    for (std::size_t i = 0; i < n; ++i) { u64 z = a[i]; z -= p2 & (0 - (u64)(z >= p2)); z -= p & (0 - (u64)(z >= p)); a[i] = z; }
}
int main() {
    Fp64 K(998244353); using P = Poly<Fp64>;
    (void)(P::random(K, 300) * P::random(K, 300));
    for (unsigned k : {11u, 13u, 16u}) {
        const auto& t = fft::tables(K, k);
        const std::size_t n = std::size_t(1) << k;
        std::vector<u64> a(n); for (auto& x : a) x = K.random();
        std::vector<u64> r0 = a, r1 = a, r2 = a;
        fft::word::transform(t, r0.data(), n, false); sin_saltos(t, r1.data(), n); perezosa(t, r2.data(), n);
        const bool ok = r0 == r1 && r0 == r2;
        auto mide = [&](auto&& f) { std::vector<u64> b = a; int v = 0; auto t0 = reloj::now(); double dt = 0;
            do { f(b.data()); ++v; dt = std::chrono::duration<double, std::micro>(reloj::now() - t0).count(); } while (dt < 100000);
            return dt / v * 1e3 / (n / 2 * k); };
        // Nota: se transforma repetidamente el mismo vector (datos pseudoaleatorios en cada pasada).
        const double c0 = mide([&](u64* b){ fft::word::transform(t, b, n, false); });
        const double c1 = mide([&](u64* b){ sin_saltos(t, b, n); });
        const double c2 = mide([&](u64* b){ perezosa(t, b, n); });
        std::printf("n = 2^%u  iguales=%d  actual %.2f ns  sin saltos %.2f ns  perezosa %.2f ns por mariposa\n", k, ok, c0, c1, c2);
    }
}
