// calibra_ddf.cpp — Fase 5: DDF de Kaltofen–Shoup (distinct_degree) frente a
// la DDF paso a paso de la fase 4 (distinct_degree_classic).
//
// Entradas: partes libres de cuadrados de polinomios mónicos aleatorios de
// grado n (rotando 8).  Comprueba que las dos salidas son idénticas y mide los
// tiempos para fijar kDdfKsMin.  Para medir el camino nuevo por debajo del
// umbral, compilar con -DALCP_DDF_KS_MIN=33 (el mínimo permitido).
//
//   g++ -std=c++23 -O3 -DALCP_DDF_KS_MIN=33 -Iinclude bench/calibra_ddf.cpp -lgmpxx -lgmp -static -o build/calibra_ddf.exe
#include "alcp/alcp.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <vector>

using namespace alcp;
using reloj = std::chrono::steady_clock;

template <class Fn>
double medir_ms(Fn&& fn) {
    double mejor = 1e300;
    for (int tanda = 0; tanda < 3; ++tanda) {
        int veces = 0;
        auto t0 = reloj::now();
        double dt = 0;
        do {
            fn(veces % 8);
            ++veces;
            dt = std::chrono::duration<double, std::milli>(reloj::now() - t0).count();
        } while (dt < 60 && veces < 1000);
        mejor = std::min(mejor, dt / veces);
    }
    return mejor;
}

template <class Field>
void calibrar(const Field& K, const char* nombre, std::vector<int> grados) {
    using P = Poly<Field>;
    std::printf("\n%s\n%6s %12s %12s %8s %8s\n", nombre, "n", "fase 4 ms", "KS ms", "x", "iguales");
    (void)(P::random(K, 300) * P::random(K, 300));  // tablas FFT, si las hay
    for (int n : grados) {
        std::vector<P> f;
        for (int i = 0; i < 8; ++i) f.push_back(P::random(K, n).monic().radical());
        bool iguales = true;
        for (int i = 0; i < 8; ++i) {
            const auto a = f[i].distinct_degree(), b = f[i].distinct_degree_classic();
            iguales = iguales && a.size() == b.size();
            for (std::size_t t = 0; iguales && t < a.size(); ++t)
                iguales = a[t].first == b[t].first && a[t].second == b[t].second;
        }
        const double tc = medir_ms([&](int i) { (void)f[i].distinct_degree_classic(); });
        const double tk = medir_ms([&](int i) { (void)f[i].distinct_degree(); });
        std::printf("%6d %12.3f %12.3f %8.2f %8s\n", n, tc, tk, tc / tk, iguales ? "sí" : "NO");
    }
}

int main() {
    set_seed(3);
    Fp64 w2(2), w101(101), wntt(998244353), wm61((BigInt(1) << 61) - 1);
    Fp g127((BigInt(1) << 127) - 1);
    Fq64 w9(Fp64(BigInt(3)), std::vector<BigInt>{1, 0, 1});
    const std::vector<int> g{24, 32, 40, 48, 64, 96, 150, 300};
    calibrar(w2, "F_2", {24, 32, 48, 64, 96, 150, 300, 600});
    calibrar(w101, "F_101", g);
    calibrar(wntt, "F_998244353", g);
    calibrar(wm61, "F_(2^61-1)", {24, 32, 48, 64, 96, 150});
    calibrar(w9, "F_9 (Fq64)", {24, 32, 48, 64, 96, 150});
    calibrar(g127, "F_(2^127-1) (GMP)", {24, 32, 48, 64, 96});
    return 0;
}
