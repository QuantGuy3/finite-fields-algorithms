// calibra_newton.cpp — Calibración de los umbrales de la división de Newton (fase 4).
//
// Para cada tipo de cuerpo y cada tamaño n (divisor de grado n, cociente de n
// coeficientes, es decir, dividendo de grado 2n - 1, el caso de pow_mod) mide:
//   * escolar  : divmod_school
//   * newton   : divmod_newton (calcula el inverso de la serie en cada llamada)
//   * reductor : divmod_with_inverse con el inverso ya calculado (lo que hace
//                Reducer dentro de pow_mod)
// Las tablas FFT del cuerpo se construyen antes (como ocurre dentro de pow_mod,
// cuyos productos ya las piden).  Imprime el mínimo de 5 tandas de >= 30 ms.
//
// ENTRADAS VARIADAS: cada medida rota entre 16 dividendos distintos (y, en las
// divisiones aisladas, 16 divisores distintos).  Con una sola entrada repetida
// el predictor de saltos del procesador aprende las restas condicionales de la
// NTT y la FFT parece hasta 3 veces más rápida de lo que es con datos reales
// (medido en Zen 3); la primera calibración cayó en ese sesgo.
// La columna esc/red compara, con el MISMO divisor fijo, la escolar con la
// reducción con el inverso precalculado.
//
//   g++ -std=c++23 -O3 -Iinclude bench/calibra_newton.cpp -lgmpxx -lgmp -static -o build/calibra_newton.exe
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
        } while (dt < 30000);
        mejor = std::min(mejor, dt / veces);
    }
    return mejor;
}

template <class Field>
void calibrar(const Field& K, const char* nombre) {
    using P = Poly<Field>;
    std::printf("\n%s\n%8s %12s %12s %12s %9s %9s\n", nombre, "n", "escolar us", "newton us", "reductor us",
                "esc/new", "esc/red");
    // Tablas FFT construidas (si el cuerpo admite FFT).
    (void)(P::random(K, 300) * P::random(K, 300));
    for (int n : {32, 64, 96, 128, 192, 256, 384, 512, 1024, 2048}) {
        std::vector<P> gs, fs;
        for (int i = 0; i < 16; ++i) {
            gs.push_back(P::random(K, n));
            fs.push_back(P::random(K, 2 * n - 1));
        }
        // Inverso de rev(g0) calculado una vez (lo que hace Reducer en pow_mod).
        const P& g0 = gs[0];
        std::vector<typename Field::Elem> rg;
        for (int j = 0; j < n; ++j) rg.push_back(g0.coeffs()[static_cast<std::size_t>(n - j)]);
        const auto ginv = P::inverse_series(K, rg, static_cast<std::size_t>(n));
        const double te = medir_us([&](int i) { (void)fs[i].divmod_school(gs[i]); });
        const double tn = medir_us([&](int i) { (void)fs[i].divmod_newton(gs[i]); });
        const double tr0 = medir_us([&](int i) { (void)fs[i].divmod_school(g0); });
        const double tr = medir_us([&](int i) { (void)fs[i].divmod_with_inverse(g0, ginv); });
        std::printf("%8d %12.1f %12.1f %12.1f %9.2f %9.2f\n", n, te, tn, tr, te / tn, tr0 / tr);
    }
}

int main() {
    set_seed(1);
    const std::vector<BigInt> mod_3_30{2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0,
                                       2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1};
    Fp64 w101(101), wntt(998244353), wm61((BigInt(1) << 61) - 1);
    Fp g127((BigInt(1) << 127) - 1);
    Fp g256(BigInt("52435875175126190479447740508185965837690552500527637822603658699938581184513"));
    Fq64 w9(Fp64(BigInt(3)), std::vector<BigInt>{1, 0, 1});
    Fq64 w256(Fp64(BigInt(2)), std::vector<BigInt>{1, 1, 0, 1, 1, 0, 0, 0, 1});
    Fq64 w3e30(Fp64(BigInt(3)), mod_3_30);
    Fq64 wm61e2(Fp64((BigInt(1) << 61) - 1), std::vector<BigInt>{BigInt(1), BigInt(0), BigInt(1)});  // x^2 + 1 (p = 3 mod 4)
    calibrar(wntt, "Fp64 p = 998244353 (NTT en palabras)");
    calibrar(w101, "Fp64 p = 101 (sin NTT: Karatsuba)");
    calibrar(wm61, "Fp64 p = 2^61 - 1 (NTT pequeña)");
    calibrar(g127, "Fp GMP p = 2^127 - 1");
    calibrar(g256, "Fp GMP p de 255 bits (NTT multi-limb)");
    calibrar(w9, "Fq64 F_9");
    calibrar(w256, "Fq64 F_256 (sin FFT)");
    calibrar(w3e30, "Fq64 F_(3^30)");
    calibrar(wm61e2, "Fq64 F_((2^61-1)^2) (con FFT)");
    return 0;
}
