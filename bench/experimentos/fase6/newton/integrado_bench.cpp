// integrado_bench.cpp — Fase 6, tema «newton»: medición de la integración en la
// librería (producto medio en Newton y Reducer/Composer con transformadas
// precalculadas) en F_998244353 con Fp64.  Se compila igual contra la librería
// anterior y la nueva para comparar antes/después (mismo checksum = mismos
// resultados).  Cuerpo «caliente» (tablas FFT ya construidas): es el único caso
// en que cambia algo; el frío sigue por Kronecker como antes.
//
//   g++ -std=c++23 -O2 -Iinclude bench/experimentos/fase6/newton/integrado_bench.cpp -lgmpxx -lgmp -static -o build/integrado_bench.exe
//   ALCP_THREADS=1 build/integrado_bench.exe
#include "alcp/alcp.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace alcp;
using P = Poly<Fp64>;
using reloj = std::chrono::steady_clock;

// g_chk: checksum de los resultados (una vez por entrada, fuera de los bucles
// cronometrados, que repiten un número variable de veces); g_sumidero solo
// impide que el compilador elimine el trabajo cronometrado.
static std::uint64_t g_chk = 0, g_sumidero = 0;
static void usa(const P& v) { g_sumidero += v.coeffs().size() + (v.is_zero() ? 0 : v.coeffs()[0]); }
static void chk(const P& v) {
    for (auto x : v.coeffs()) g_chk = g_chk * 1000003u + x;
    g_chk = g_chk * 1000003u + v.coeffs().size();
}

// Mejor de 5 tandas de >= 40 ms, rotando 8 entradas.
template <class Fn>
double medir_us(Fn&& fn) {
    double mejor = 1e300;
    for (int tanda = 0; tanda < 5; ++tanda) {
        int veces = 0;
        auto t0 = reloj::now();
        double dt = 0;
        do {
            fn(veces % 8);
            ++veces;
            dt = std::chrono::duration<double, std::micro>(reloj::now() - t0).count();
        } while (dt < 40000);
        mejor = std::min(mejor, dt / veces);
    }
    return mejor;
}

int main() {
    Fp64 K(BigInt("998244353"));
    set_seed(12345);
    (void)(P::random(K, 300) * P::random(K, 300));  // tablas FFT (cuerpo caliente)

    // 4000/2000: caso favorable (todos los niveles de Newton justo por debajo de
    // una potencia de 2); 2200/1100 y 4100/2050: desfavorables (d1 algo por
    // encima de una potencia de 2).
    for (auto [df, dg] : {std::pair{4000, 2000}, std::pair{2200, 1100}, std::pair{4100, 2050}}) {
        std::vector<P> f, g;
        for (int i = 0; i < 8; ++i) { f.push_back(P::random(K, df)); g.push_back(P::random(K, dg)); }
        const double t = medir_us([&](int i) { auto [q, r] = f[i].divmod(g[i]); usa(q); usa(r); });
        for (int i = 0; i < 8; ++i) { auto [q, r] = f[i].divmod(g[i]); chk(q); chk(r); }
        std::printf("divmod %d/%d                : %9.1f us\n", df, dg, t);
    }
    // n = 513 = 2^9 + 1: tamaño desfavorable para el Reducer.
    for (int n : {300, 513, 1000}) {
        P m = P::random(K, n);
        P::Reducer red(m);
        std::vector<P> a;
        for (int i = 0; i < 8; ++i) a.push_back(P::random(K, n - 1) * P::random(K, n - 1));
        usa(red.reduce(a[0]));  // inverso (y transformadas) ya calculados
        const double tr = medir_us([&](int i) { usa(red.reduce(a[i])); });
        std::vector<P> b;
        for (int i = 0; i < 8; ++i) b.push_back(P::random(K, n - 1));
        const BigInt q = K.cardinality();
        const double tp = medir_us([&](int i) { usa(P::pow_mod(b[i], q, red)); });
        const std::size_t k = static_cast<std::size_t>(std::ceil(std::sqrt(double(n))));
        P::Composer comp(b[0], red, k);
        const double tc = medir_us([&](int i) { usa(comp.compose(b[(i + 1) % 8])); });
        for (int i = 0; i < 8; ++i) { chk(red.reduce(a[i])); chk(P::pow_mod(b[i], q, red)); chk(comp.compose(b[(i + 1) % 8])); }
        std::printf("n = %4d: reduce %8.1f us | pow_mod(q) %9.1f us | compose (k=%zu) %9.1f us\n", n, tr, tp, k,
                    tc);
    }
    std::printf("checksum %016llx  (sumidero %llu)\n", static_cast<unsigned long long>(g_chk),
                static_cast<unsigned long long>(g_sumidero % 1000));
}
