// calibra_kronecker.cpp — Calibración del producto por sustitución de Kronecker (fase 4, nivel B).
//
// Para cada cuerpo y tamaño n (dos factores de n coeficientes) compara el
// despacho Auto actual (escolar / Karatsuba / FFT, con sus hilos) con el
// producto de Kronecker sobre Z con GMP.  Las entradas rotan entre 16 pares
// distintos (con una sola entrada repetida el predictor de saltos sesga las
// medidas de la NTT; ver docs/optimizacion/fase4.md).  Comprueba además que los
// dos productos coinciden.  Ejecutar también con ALCP_THREADS=1.
//
//   g++ -std=c++23 -O3 -Iinclude bench/calibra_kronecker.cpp -lgmpxx -lgmp -static -o build/calibra_kronecker.exe
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

template <class Field>
void calibrar(const Field& K, const char* nombre, int nmax) {
    using P = Poly<Field>;
    std::printf("\n%s\n%8s %12s %12s %9s\n", nombre, "n", "Auto us", "Kronecker us", "Auto/KS");
    (void)(P::random(K, 300) * P::random(K, 300));  // tablas FFT, si las hay
    for (int n = 16; n <= nmax; n *= 2) {
        std::vector<P> as, bs;
        for (int i = 0; i < 16; ++i) {
            as.push_back(P::random(K, n - 1));
            bs.push_back(P::random(K, n - 1));
        }
        if (!(as[0] * bs[0] == as[0].mul(bs[0], MulAlgo::Kronecker)) ||
            !(as[1] * as[1] == as[1].mul(as[1], MulAlgo::Kronecker)))
            std::printf("  ¡DISTINTOS en n = %d!\n", n);
        const double ta = medir_us([&](int i) { (void)(as[i] * bs[i]); });
        const double tk = medir_us([&](int i) { (void)as[i].mul(bs[i], MulAlgo::Kronecker); });
        std::printf("%8d %12.1f %12.1f %9.2f\n", n, ta, tk, ta / tk);
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
    Fq64 wm61e2(Fp64((BigInt(1) << 61) - 1), std::vector<BigInt>{1, 0, 1});
    calibrar(wntt, "Fp64 p = 998244353 (NTT)", 8192);
    calibrar(w101, "Fp64 p = 101", 8192);
    calibrar(wm61, "Fp64 p = 2^61 - 1", 8192);
    calibrar(g127, "Fp GMP p = 2^127 - 1", 4096);
    calibrar(g256, "Fp GMP p de 255 bits (NTT multi-limb)", 4096);
    calibrar(w9, "Fq64 F_9", 4096);
    calibrar(w256, "Fq64 F_256", 2048);
    calibrar(w3e30, "Fq64 F_(3^30)", 1024);
    calibrar(wm61e2, "Fq64 F_((2^61-1)^2) (FFT genérica)", 2048);
    return 0;
}
