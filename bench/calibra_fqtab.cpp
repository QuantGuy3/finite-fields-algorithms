// calibra_fqtab.cpp — Calibración de las constantes de Poly<FqTab> (fase 8,
// revisión de la fase 7): umbral de Karatsuba en paralelo (par_karatsuba_min)
// y presupuesto de las potencias de Frobenius (kFrobPresupuesto), que FqTab
// heredaba de Fp64 sin medir.
//
// Las constantes son static constexpr de Poly, así que cada valor candidato se
// mide compilando ESTE programa contra una copia de include/ con la constante
// cambiada (solo para FqTab, Elem = uint32_t) y ejecutando los binarios
// intercalados varias tandas (mejor tanda de cada uno).  Modos:
//   mul  productos de grado 512..4096 (Auto = Karatsuba, con hilos);
//   irr  is_irreducible de un irreducible de grado 800 (p = 2) o 600 (p impar);
//   fac  factor() de un aleatorio de grado 200 y 400.
// Cuerpos: F_256 (módulo de AES), F_(2^12), F_(3^7) y F_(7^4).
//
//   g++ -std=c++23 -O3 -Iinclude bench/calibra_fqtab.cpp -lgmpxx -lgmp -static -o build/calibra_fqtab.exe
//   build/calibra_fqtab.exe mul
#include "alcp/alcp.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <vector>

using namespace alcp;
using reloj = std::chrono::steady_clock;

template <class F>
double mejor_ms(F&& f, int reps) {
    double mejor = 1e300;
    for (int r = 0; r < reps; ++r) {
        auto t0 = reloj::now();
        f();
        mejor = std::min(mejor, std::chrono::duration<double, std::milli>(reloj::now() - t0).count());
    }
    return mejor;
}

std::vector<BigInt> big(const std::vector<long>& v) {
    std::vector<BigInt> r;
    for (long x : v) r.emplace_back(x);
    return r;
}

int main(int argc, char** argv) {
    const char* modo = argc > 1 ? argv[1] : "mul";
    set_seed(1);
    FqTab f256(Fp64(BigInt(2)), big({1, 1, 0, 1, 1, 0, 0, 0, 1}));
    set_seed(2);
    FqTab f2e12(BigInt(2), 12);
    set_seed(3);
    FqTab f3e7(BigInt(3), 7);
    set_seed(4);
    FqTab f7e4(Fp64(BigInt(7)), big({3, 4, 1, 2, 1}));
    const FqTab* K[] = {&f256, &f2e12, &f3e7, &f7e4};
    const char* nom[] = {"F_256", "F_2^12", "F_3^7", "F_7^4"};
    for (int i = 0; i < 4; ++i) {
        using P = Poly<FqTab>;
        set_seed(10 + static_cast<unsigned long>(i));
        if (!std::strcmp(modo, "mul")) {
            for (int d : {512, 1023, 1024, 1500, 2048, 4096}) {
                P a = P::random(*K[i], d), b = P::random(*K[i], d);
                const double t = mejor_ms([&] { volatile auto n = (a * b).coeffs().size(); (void)n; }, 7);
                std::printf("%s mul d=%d %.3f\n", nom[i], d, t);
            }
        } else if (!std::strcmp(modo, "irr")) {
            for (int n : {i < 2 ? 800 : 600}) {  // con p impar, 800 tarda demasiado
                P f = P::random(*K[i], n).monic();
                while (!f.is_irreducible()) f = P::random(*K[i], n).monic();
                const double t = mejor_ms([&] { volatile bool b = f.is_irreducible(); (void)b; }, 2);
                std::printf("%s irr n=%d %.1f\n", nom[i], n, t);
            }
        } else {
            for (int n : {200, 400}) {
                P f = P::random(*K[i], n).monic();
                const double t = mejor_ms([&] { volatile auto s = f.factor().factors.size(); (void)s; }, 2);
                std::printf("%s fac n=%d %.1f\n", nom[i], n, t);
            }
        }
        std::fflush(stdout);
    }
}
