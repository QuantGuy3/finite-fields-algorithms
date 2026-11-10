// bench_cpp.cpp — Fase 6, tema «python»: referencia en C++ puro del coste de las
// mismas operaciones que se miden desde Python (mide_python.py).
//
// Mide, con Poly<Fp64> sobre F_998244353 y las MISMAS entradas que genera
// mide_python.py (splitmix64, 16 pares rotando):
//   * producto f*g de grado 10, 100 y 1000 (operator*, MulAlgo::Auto, como __mul__);
//   * construcción de un Poly desde un std::vector<u64> ya reducido (copia + trim);
//   * operaciones de elementos (mul, add) de Fp64.
// Imprime ns por operación (mejor de 3 tandas de >= 30 ms) y un hash de los
// productos que mide_python.py compara con el de alcp_native.
//
// Compilación (desde la raíz del repositorio; portable, como el módulo):
//   g++ -std=c++23 -O3 -Iinclude bench/experimentos/fase6/python/bench_cpp.cpp
//       -lgmpxx -lgmp -static -o build/fase6/python/bench_cpp.exe
// Ejecución:
//   ALCP_THREADS=1 build/fase6/python/bench_cpp.exe
#include "alcp/alcp.hpp"

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

using alcp::Fp64;
using Poly = alcp::Poly<Fp64>;
using u64 = std::uint64_t;

// Generador splitmix64 (idéntico al de mide_python.py).
struct SplitMix {
    u64 s;
    u64 operator()() {
        u64 z = (s += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
};

// Coeficientes de un polinomio de grado d (coeficiente principal no nulo).
static std::vector<u64> coefs(SplitMix& r, u64 p, int d) {
    std::vector<u64> c(d + 1);
    for (int i = 0; i < d; ++i) c[i] = r() % p;
    c[d] = r() % (p - 1) + 1;
    return c;
}

// Mejor de 3 tandas de >= 30 ms; fn(i) hace una operación con la entrada i % 16
// y devuelve un valor que se acumula (evita que el compilador elimine el trabajo).
template <class Fn>
static double medir(Fn&& fn, u64& chk) {
    using reloj = std::chrono::steady_clock;
    double mejor = 1e300;
    for (int tanda = 0; tanda < 3; ++tanda) {
        long n = 0;
        const auto t0 = reloj::now();
        double dt = 0;
        do {
            for (int k = 0; k < 64; ++k, ++n) chk += fn(static_cast<unsigned>(n));
            dt = std::chrono::duration<double>(reloj::now() - t0).count();
        } while (dt < 0.030);
        mejor = std::min(mejor, dt / n);
    }
    return mejor * 1e9;
}

int main() {
    const u64 p = 998244353;
    Fp64 K{alcp::BigInt(std::to_string(p))};
    u64 chk = 0;
    std::printf("C++ puro, Poly<Fp64> en F_998244353 (ns por operación, mejor de 3 tandas >= 30 ms, 16 entradas rotando)\n");
    for (int d : {10, 100, 1000}) {
        SplitMix r{static_cast<u64>(d)};
        std::vector<Poly> F, G;
        std::vector<std::vector<u64>> vf;
        for (int i = 0; i < 16; ++i) {
            vf.push_back(coefs(r, p, d));
            F.emplace_back(K, vf.back());
            G.emplace_back(K, coefs(r, p, d));
        }
        // Hash de los 16 productos (comparable con el de mide_python.py).
        u64 h = 0;
        for (int i = 0; i < 16; ++i) {
            const Poly P = F[i] * G[i];
            const auto& c = P.coeffs();
            for (std::size_t j = 0; j < c.size(); ++j) h += c[j] * (j + 1);
        }
        const double t_mul = medir([&](unsigned i) { return (F[i & 15] * G[i & 15]).coeffs()[d]; }, chk);
        const double t_add = medir([&](unsigned i) { return (F[i & 15] + G[i & 15]).coeffs()[d / 2]; }, chk);
        const double t_ctor = medir([&](unsigned i) { return Poly(K, vf[i & 15]).coeffs()[d / 2]; }, chk);
        std::printf("d=%-5d mul %10.1f ns   add %8.1f ns   Poly(K, vector) %8.1f ns   hash=%llu\n", d, t_mul, t_add,
                    t_ctor, static_cast<unsigned long long>(h));
    }
    // Elementos: 16 pares rotando.
    SplitMix r{7};
    u64 a[16], b[16];
    for (int i = 0; i < 16; ++i) a[i] = r() % p, b[i] = r() % p;
    const double t_emul = medir([&](unsigned i) { return K.mul(a[i & 15], b[(i + 3) & 15]); }, chk);
    const double t_eadd = medir([&](unsigned i) { return K.add(a[i & 15], b[(i + 3) & 15]); }, chk);
    std::printf("elementos: mul %.2f ns   add %.2f ns\n", t_emul, t_eadd);
    std::printf("(checksum %llu)\n", static_cast<unsigned long long>(chk));
}
