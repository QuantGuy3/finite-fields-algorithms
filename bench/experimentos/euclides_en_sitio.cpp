// Experimento: Euclides clásico de la librería frente a un Euclides en sitio
// sobre vectores u64 con multiplicador de Shoup y sin saltos (mismo algoritmo,
// mismos cocientes; se comprueba que el mcd es idéntico).
//   g++ -std=c++23 -O3 -march=native -Iinclude bench/experimentos/euclides_en_sitio.cpp -lgmpxx -lgmp -static -o build/euclides_en_sitio.exe
//   (investigación de la fase 5: docs/optimizacion/investigacion-fase5.md)
#include "alcp/alcp.hpp"
#include <chrono>
#include <cstdio>
#include <vector>
using namespace alcp;
using u64 = std::uint64_t; using u128 = unsigned __int128;
using reloj = std::chrono::steady_clock;

static u64 inv_mod(u64 a, u64 p) {  // Euclides extendido en 64 bits
    std::int64_t t = 0, nt = 1; u64 r = p, nr = a;
    while (nr) { u64 q = r / nr; std::int64_t tt = t - (std::int64_t)q * nt; t = nt; nt = tt; u64 rr = r - q * nr; r = nr; nr = rr; }
    return t < 0 ? (u64)(t + (std::int64_t)p) : (u64)t;
}

// a <- a mod b en sitio; devuelve el nuevo grado de a (-1 si es cero).
static int resto_en_sitio(std::vector<u64>& a, int da, const std::vector<u64>& b, int db, u64 p) {
    const u64 il = inv_mod(b[db], p);
    for (int dr = da; dr >= db; --dr) {
        const u64 lr = a[dr];
        if (lr == 0) continue;
        const u64 c = p - static_cast<u64>((static_cast<u128>(lr) * il) % p);   // -factor
        const u64 cs = static_cast<u64>((static_cast<u128>(c) << 64) / p);       // cociente de Shoup
        u64* r = a.data() + (dr - db);
        for (int j = 0; j < db; ++j) {
            const u64 q = static_cast<u64>((static_cast<u128>(cs) * b[j]) >> 64);
            u64 v = c * b[j] - q * p;  v -= p & (0 - (u64)(v >= p));
            u64 s = r[j] + v - p;      s += p & (0 - (s >> 63));
            r[j] = s;
        }
        a[dr] = 0;
    }
    int d = db - 1;
    while (d >= 0 && a[d] == 0) --d;
    return d;
}
static std::vector<u64> gcd_en_sitio(std::vector<u64> a, std::vector<u64> b, u64 p) {
    int da = (int)a.size() - 1, db = (int)b.size() - 1;
    while (db >= 0) { int dr = resto_en_sitio(a, da, b, db, p); std::swap(a, b); da = db; db = dr; }
    a.resize(da + 1);
    const u64 il = inv_mod(a[da], p);
    for (auto& x : a) x = (u64)((u128)x * il % p);
    return a;
}
int main() {
    for (u64 p : {998244353ull, 101ull, (1ull << 61) - 1}) {
        Fp64 K(BigInt(std::to_string(p))); using P = Poly<Fp64>;
        for (int n : {300, 1000}) {
            std::vector<P> A, B; for (int i = 0; i < 8; ++i) { A.push_back(P::random(K, n)); B.push_back(P::random(K, n - 1)); }
            bool ok = true;
            for (int i = 0; i < 8; ++i) ok = ok && gcd_en_sitio(A[i].coeffs(), B[i].coeffs(), p) == P::gcd_classic(A[i], B[i]).coeffs();
            auto mide = [&](auto&& f) { double mejor = 1e300; for (int t = 0; t < 3; ++t) { int v = 0; auto t0 = reloj::now(); double dt = 0;
                do { f(v % 8); ++v; dt = std::chrono::duration<double, std::micro>(reloj::now() - t0).count(); } while (dt < 50000);
                mejor = std::min(mejor, dt / v); } return mejor; };
            const double t0 = mide([&](int i){ (void)P::gcd_classic(A[i], B[i]); });
            const double t1 = mide([&](int i){ (void)gcd_en_sitio(A[i].coeffs(), B[i].coeffs(), p); });
            std::printf("p=%llu n=%d iguales=%d  librería %.0f us  en sitio %.0f us  (%.1fx, %.2f ns por operación)\n",
                        (unsigned long long)p, n, ok, t0, t1, t0 / t1, t1 * 1e3 / (double(n) * n));
        }
    }
}
