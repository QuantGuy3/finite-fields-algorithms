// perfil_factor.cpp — Desglose del tiempo de las operaciones que siguen
// lejos de FLINT, para priorizar la investigación de la fase 5.
//
// Mide, con entradas variadas (rotando 8), en F_998244353 (Fp64, con NTT) y
// F_101:
//   * coste de las primitivas: producto, cuadrado, divmod, gcd, pow_mod;
//   * factor() desglosado en radical, distinct_degree (Frobenius + MCD) y
//     equal_degree (Cantor–Zassenhaus).
//
//   g++ -std=c++23 -O3 -march=native -Iinclude bench/experimentos/perfil_factor.cpp -lgmpxx -lgmp -static -o build/perfil_factor.exe
#include "alcp/alcp.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <vector>

using namespace alcp;
using reloj = std::chrono::steady_clock;

// Mejor de 3 tandas de al menos 30 ms; fn recibe el índice de la entrada (0..7).
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
void primitivas(const Field& K, const char* nombre) {
    using P = Poly<Field>;
    std::printf("\n== %s: primitivas (us)\n%6s %10s %10s %12s %10s %12s\n", nombre, "n", "mul", "sqr",
                "divmod 2n/n", "gcd n", "pow_mod q");
    (void)(P::random(K, 300) * P::random(K, 300));  // tablas FFT
    const BigInt q = K.cardinality();
    for (int n : {100, 300, 1000, 2000}) {
        std::vector<P> a, b, c, m;
        for (int i = 0; i < 8; ++i) {
            a.push_back(P::random(K, n));
            b.push_back(P::random(K, n));
            c.push_back(P::random(K, 2 * n));
            m.push_back(P::random(K, n).monic());
        }
        const double tm = medir_us([&](int i) { (void)(a[i] * b[i]); });
        const double ts = medir_us([&](int i) { (void)(a[i] * a[i]); });
        const double td = medir_us([&](int i) { (void)c[i].divmod(b[i]); });
        const double tg = medir_us([&](int i) { (void)P::gcd(a[i], b[i]); });
        typename P::Reducer red(m[0]);
        const double tp = n <= 1000 ? medir_us([&](int i) { (void)P::pow_mod(a[i] % m[0], q, red); }) : 0.0;
        std::printf("%6d %10.1f %10.1f %12.1f %10.1f %12.1f\n", n, tm, ts, td, tg, tp);
    }
}

template <class Field>
void desglose_factor(const Field& K, int n, const char* nombre) {
    using P = Poly<Field>;
    std::vector<P> f;
    for (int i = 0; i < 8; ++i) f.push_back(P::random(K, n).monic());
    double t_rad = 0, t_ddf = 0, t_edf = 0, t_mult = 0, t_frob = 0, t_gcd = 0;
    int veces = 0, n_gcd = 0;
    const auto ahora = [] { return reloj::now(); };
    const auto us = [](auto a, auto b) { return std::chrono::duration<double, std::micro>(b - a).count(); };
    for (int rep = 0; rep < 3; ++rep)
        for (int i = 0; i < 8; ++i) {
            auto t0 = ahora();
            P g = f[i].radical();
            auto t1 = ahora();
            // distinct_degree reproducido con el mismo código para separar Frobenius y MCD
            std::vector<std::pair<int, P>> dd;
            {
                P gg = g.monic(), x = P::variable(K), H = x % gg;
                typename P::Frobenius frob(gg);
                int d = 1;
                while (gg.degree() >= 2 * d) {
                    auto a = ahora();
                    H = frob.siguiente(H);
                    auto b = ahora();
                    P hd = P::gcd(H - x, gg);
                    auto c = ahora();
                    t_frob += us(a, b);
                    t_gcd += us(b, c);
                    ++n_gcd;
                    if (hd.degree() > 0) {
                        dd.emplace_back(d, hd.monic());
                        gg = gg / hd;
                        H = H % gg;
                        frob.cambiar_modulo(gg);
                    }
                    ++d;
                }
                if (gg.degree() > 0) dd.emplace_back(gg.degree(), gg.monic());
            }
            auto t2 = ahora();
            std::vector<P> irr;
            for (auto& [d, prod] : dd) {
                auto part = prod.equal_degree(d);
                irr.insert(irr.end(), part.begin(), part.end());
            }
            auto t3 = ahora();
            for (auto& u : irr) (void)f[i].multiplicity(u);
            auto t4 = ahora();
            t_rad += us(t0, t1);
            t_ddf += us(t1, t2);
            t_edf += us(t2, t3);
            t_mult += us(t3, t4);
            ++veces;
        }
    const double tot = t_rad + t_ddf + t_edf + t_mult;
    std::printf("\n== %s: factor grado %d (media de %d, ms)\n", nombre, n, veces);
    std::printf("  total %.2f | radical %.2f (%.0f%%) | DDF %.2f (%.0f%%: Frobenius %.2f, MCD %.2f en %d MCD) | "
                "EDF %.2f (%.0f%%) | multiplicidad %.2f (%.0f%%)\n",
                tot / veces / 1e3, t_rad / veces / 1e3, 100 * t_rad / tot, t_ddf / veces / 1e3, 100 * t_ddf / tot,
                t_frob / veces / 1e3, t_gcd / veces / 1e3, n_gcd / veces, t_edf / veces / 1e3, 100 * t_edf / tot,
                t_mult / veces / 1e3, 100 * t_mult / tot);
}

int main() {
    set_seed(5);
    Fp64 wntt(998244353), w101(101);
    primitivas(wntt, "F_998244353");
    primitivas(w101, "F_101");
    desglose_factor(wntt, 300, "F_998244353");
    desglose_factor(w101, 100, "F_101");
    desglose_factor(w101, 300, "F_101");
    return 0;
}
