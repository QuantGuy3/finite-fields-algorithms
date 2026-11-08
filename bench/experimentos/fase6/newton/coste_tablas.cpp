// coste_tablas.cpp — Fase 6, tema «newton»: ¿cuánto cuesta construir las tablas
// FFT de un cuerpo Fp64 recién creado?  Sirve para valorar la propuesta (d):
// que mul_no_rng (productos de Newton) use tablas PROPIAS con una raíz
// determinista (sin el RNG) cuando la caché principal aún está vacía, como pasa
// en divmod desde Python (cuerpo nuevo que nunca ha hecho un producto Auto).
// Aquí se mide fft::tables (que sí pide la raíz al RNG: en un prototipo da igual;
// la búsqueda determinista de la raíz cuesta lo mismo, unas pocas potencias).
// También comprueba que la raíz NO influye en el resultado: convoluciona con las
// tablas de dos cuerpos con raíces distintas y compara.
//
//   g++ -std=c++23 -O3 -Iinclude bench/experimentos/fase6/newton/coste_tablas.cpp -lgmpxx -lgmp -static -o build/fase6/newton/coste_tablas.exe
//   ALCP_THREADS=1 build/fase6/newton/coste_tablas.exe
#include "alcp/alcp.hpp"

#include <chrono>
#include <cstdio>
#include <memory>
#include <vector>

using namespace alcp;
using reloj = std::chrono::steady_clock;

int main() {
    const BigInt p("998244353");
    unsigned long long chk = 0;
    // Coste de construir las tablas 2^1..2^K en un cuerpo nuevo (mejor de 3 tandas >= 30 ms).
    for (unsigned K : {11u, 12u, 13u}) {
        double mejor = 1e300;
        for (int tanda = 0; tanda < 3; ++tanda) {
            int veces = 0;
            auto t0 = reloj::now();
            double dt = 0;
            do {
                auto F = std::make_unique<Fp64>(p);
                for (unsigned k = 1; k <= K; ++k) chk += fft::tables(*F, k).word.fwd[0];
                ++veces;
                dt = std::chrono::duration<double, std::micro>(reloj::now() - t0).count();
            } while (dt < 30000);
            mejor = std::min(mejor, dt / veces);
        }
        std::printf("tablas 2^1..2^%u en un cuerpo nuevo: %8.1f us\n", K, mejor);
    }
    // La raíz no influye en el producto: dos cuerpos (dos raíces aleatorias distintas).
    Fp64 A(p), B(p);
    (void)fft::tables(A, 12);
    (void)fft::tables(B, 12);
    int dif = 0;
    for (int i = 0; i < 8; ++i) {
        auto f = Poly<Fp64>::random(A, 1500).coeffs(), g = Poly<Fp64>::random(A, 1200).coeffs();
        if (fft::convolve(A, f, g) != fft::convolve(B, f, g)) ++dif;
    }
    std::printf("raíces distintas: %s;  convoluciones distintas: %d de 8   (checksum %llx)\n",
                A.fft_base_root() != B.fft_base_root() ? "sí" : "no (coincidencia)", dif, chk);
    return dif;
}
