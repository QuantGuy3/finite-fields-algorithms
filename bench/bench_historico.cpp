// bench_historico.cpp — Casos históricos de bench.cpp con repeticiones y los dos backends.
//
// bench.cpp (que se conserva tal cual) mide una sola vez:
//   * factor() de un polinomio aleatorio de grado 5 sobre F_(3^30);
//   * el producto de dos polinomios de grado 4000 en F_998244353 con los
//     algoritmos escolar, Karatsuba y FFT.
// Este programa mide exactamente esos casos (mismas llamadas, misma semilla), pero
// repite cada medida y da la mediana y el mínimo, y, si la cabecera fp64.hpp
// existe, también con el backend de 64 bits (Fp64/Fq64).  Compila sin cambios
// contra las cabeceras originales (baseline-original), que no tienen fp64.hpp.
//
// Uso: bench_historico.exe [repeticiones]   (por defecto 5)
// Salida: una línea por caso, «nombre<TAB>mediana_ms<TAB>min_ms<TAB>n».
#include "alcp/alcp.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

#if __has_include("alcp/fp64.hpp")
#define ALCP_HISTORICO_FP64 1
#endif

using namespace alcp;
using clk = std::chrono::steady_clock;

namespace {

int g_reps = 5;

// Mide f() g_reps veces (tras un calentamiento) e imprime mediana y mínimo.
template <class F>
void medir(const char* nombre, F&& f) {
    f();  // calentamiento (tablas FFT, cachés)
    std::vector<double> t;
    for (int i = 0; i < g_reps; ++i) {
        const auto t0 = clk::now();
        f();
        t.push_back(std::chrono::duration<double, std::milli>(clk::now() - t0).count());
    }
    std::sort(t.begin(), t.end());
    const double med = (t.size() % 2) ? t[t.size() / 2] : (t[t.size() / 2 - 1] + t[t.size() / 2]) / 2;
    std::printf("%s\t%.3f\t%.3f\t%d\n", nombre, med, t.front(), g_reps);
    std::fflush(stdout);
}

// factor() del caso histórico: F_(3^30) con módulo aleatorio (semilla 1, como
// bench.cpp) y un polinomio aleatorio de grado 5.  Antes de cada factorización
// se fija la semilla, para que todas las repeticiones sigan el mismo camino.
template <class FQ>
void caso_factor(const char* nombre) {
    set_seed(1);
    FQ K(BigInt(3), 30);
    using PQ = Poly<FQ>;
    const PQ f = PQ::random(K, 5);
    medir(nombre, [&] {
        set_seed(12345);
        auto fac = f.factor();
        long grado = 0;
        for (auto& pr : fac.factors) grado += pr.second * pr.first.degree();
        if (grado != f.degree()) std::printf("  !! grado distinto\n");
    });
}

// Producto de dos polinomios aleatorios de grado 4000 en F_998244353.
template <class FP>
void caso_mul(const char* sufijo) {
    set_seed(2);
    FP P(BigInt(998244353));
    using PF = Poly<FP>;
    const PF a = PF::random(P, 4000), b = PF::random(P, 4000);
    char n[96];
    std::snprintf(n, sizeof n, "mul4000.school.%s", sufijo);
    medir(n, [&] { volatile auto r = a.mul(b, MulAlgo::Schoolbook).degree(); (void)r; });
    std::snprintf(n, sizeof n, "mul4000.karatsuba.%s", sufijo);
    medir(n, [&] { volatile auto r = a.mul(b, MulAlgo::Karatsuba).degree(); (void)r; });
    std::snprintf(n, sizeof n, "mul4000.fft.%s", sufijo);
    medir(n, [&] { volatile auto r = a.mul(b, MulAlgo::FFT).degree(); (void)r; });
    std::snprintf(n, sizeof n, "mul4000.auto.%s", sufijo);
    medir(n, [&] { volatile auto r = (a * b).degree(); (void)r; });
}

}  // namespace

int main(int argc, char** argv) {
    if (argc > 1) g_reps = std::max(1, std::atoi(argv[1]));
    caso_factor<Fq>("factor.q3e30.d5.gmp");
    caso_mul<Fp>("gmp");
#ifdef ALCP_HISTORICO_FP64
    caso_factor<Fq64>("factor.q3e30.d5.u64");
    caso_mul<Fp64>("u64");
#endif
    return 0;
}
