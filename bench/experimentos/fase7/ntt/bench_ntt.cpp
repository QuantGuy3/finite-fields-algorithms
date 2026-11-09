// Fase 7, propuesta 4 («ntt»): medición del producto por NTT antes y después.
//
// Compilación (misma fuente contra las cabeceras nuevas o contra las de 80b983d
// extraídas a otro directorio, para el A/B):
//   g++ -std=c++23 -O2 -Iinclude bench/experimentos/fase7/ntt/bench_ntt.cpp -lgmpxx -lgmp -static -o build/bench_ntt.exe
// Ejecución: ALCP_THREADS=1 build/bench_ntt.exe
//
// Mide fft::convolve (producto y cuadrado) con transformadas de tamaño
// n = 256, 2048 y 32768 (grado 100, 1000 y 10000) en F_998244353 y en el primo
// de 62 bits 4179340454199820289, con Fp64 y con Fp (GMP).  8 entradas
// rotando, mejor de 5 tandas de >= 30 ms.  Imprime además una suma de control
// de los resultados para comprobar que ambas versiones calculan lo mismo.
#include "alcp/alcp.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <vector>

using namespace alcp;
using reloj = std::chrono::steady_clock;
static volatile unsigned long long sumidero = 0;

template <class F>
static double medir(F&& f) {
    double mejor = 1e300;
    for (int tanda = 0; tanda < 5; ++tanda) {
        long reps = 0;
        unsigned long long acc = 0;
        double dt = 0;
        const auto t0 = reloj::now();
        do {
            acc += f(static_cast<int>(reps % 8));
            ++reps;
            dt = std::chrono::duration<double, std::micro>(reloj::now() - t0).count();
        } while (dt < 30e3);
        sumidero = sumidero + acc;
        mejor = std::min(mejor, dt / reps);
    }
    return mejor;
}

template <class Field>
static void caso(const char* nombre, const Field& K) {
    set_seed(12345);
    for (std::size_t m : {std::size_t(100), std::size_t(1000), std::size_t(10000)}) {
        using E = typename Field::Elem;
        std::vector<std::vector<E>> A(8), B(8);
        for (int e = 0; e < 8; ++e) {
            for (std::size_t i = 0; i <= m; ++i) A[e].push_back(K.random());
            for (std::size_t i = 0; i <= m; ++i) B[e].push_back(K.random());
        }
        auto h = [&](const std::vector<E>& r) {
            unsigned long long s = 0;
            for (std::size_t i = 0; i < r.size(); i += 7)
                s = s * 1000003u + static_cast<unsigned long long>(fft::word::to_word(r[i]));
            return s;
        };
        unsigned long long chk = 0;
        for (int e = 0; e < 8; ++e) chk ^= h(fft::convolve(K, A[e], B[e])) ^ (h(fft::convolve(K, A[e], A[e])) << 1);
        const double tp = medir([&](int e) { return h(fft::convolve(K, A[e], B[e])); });
        const double tc = medir([&](int e) { return h(fft::convolve(K, A[e], A[e])); });
        std::printf("%-26s n = %5zu  producto %9.1f us  cuadrado %9.1f us  control %016llx\n", nombre,
                    std::size_t(1) << fft::ceil_log2(2 * m + 1), tp, tc, chk);
    }
}

int main() {
    caso("Fp64 F_998244353", Fp64(998244353));
    caso("Fp64 F_4179340454199820289", Fp64(BigInt("4179340454199820289")));
    caso("Fp   F_998244353", Fp(BigInt(998244353)));
    caso("Fp   F_4179340454199820289", Fp(BigInt("4179340454199820289")));
    std::printf("# sumidero %llu\n", static_cast<unsigned long long>(sumidero));
}
