// medir.cpp — Fase 7, propuesta 7: Fq64 frente a FqTab y FqSwar integrados en
// la librería, con las mismas 16 entradas de la fase 6
// (bench/results/fase6/factor_entradas.txt, las que mide python-flint).
//
// Medidor corregido (reserva del revisor de la fase 6): mejor de 3 tandas de
// >= 30 ms, y cada tanda recorre un número ENTERO de vueltas a las 16
// entradas (antes paraba a los 30 ms y en los casos lentos solo medía las
// primeras).  Las dos variantes se miden intercaladas (A, B, A, B...) y antes
// se comprueba que los resultados y el generador posterior son idénticos.
//
// Compilación (desde la raíz del repo):
//   g++ -std=c++23 -O3 -Iinclude bench/experimentos/fase7/fq/medir.cpp -lgmpxx -lgmp -static -o build/fase7/medir.exe
// Ejecución:  ALCP_THREADS=1 build/fase7/medir.exe
#include "alcp/alcp.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace alcp;
using reloj = std::chrono::steady_clock;

static volatile unsigned long long sumidero = 0;

template <class Fn>
double una_tanda(Fn&& fn, int nent) {
    int veces = 0;
    const auto t0 = reloj::now();
    double dt = 0;
    do {
        fn(veces % nent);
        ++veces;
        dt = std::chrono::duration<double, std::micro>(reloj::now() - t0).count();
    } while (dt < 30000 || veces % nent != 0);
    return dt / veces;
}

// Mejor de 3 tandas de cada variante, intercaladas.
template <class FA, class FB>
std::pair<double, double> medir_ab(FA&& a, FB&& b, int nent = 16) {
    double ma = 1e300, mb = 1e300;
    for (int t = 0; t < 3; ++t) {
        ma = std::min(ma, una_tanda(a, nent));
        mb = std::min(mb, una_tanda(b, nent));
    }
    return {ma, mb};
}

static std::map<std::string, std::vector<std::vector<BigInt>>> leer_entradas(const char* ruta) {
    std::map<std::string, std::vector<std::vector<BigInt>>> m;
    std::ifstream in(ruta);
    std::string linea;
    while (std::getline(in, linea)) {
        std::istringstream ss(linea);
        std::string caso, idx, c;
        if (!(ss >> caso >> idx)) continue;
        std::vector<BigInt> v;
        while (ss >> c) v.emplace_back(c);
        m[caso].push_back(std::move(v));
    }
    return m;
}

template <class K>
Poly<K> poly_de(const K& F, const std::vector<BigInt>& idx) {
    std::vector<typename K::Elem> v;
    for (const auto& i : idx) v.push_back(F.from_index(i));
    return Poly<K>(F, std::move(v));
}

// Coeficientes en la representación de Fq64.
std::vector<Fq64::Elem> vec(const Poly<Fq64>& p) { return p.coeffs(); }
std::vector<Fq64::Elem> vec(const Poly<FqSwar>& p) { return p.coeffs(); }
std::vector<Fq64::Elem> vec(const Poly<FqTab>& p) {
    std::vector<Fq64::Elem> r;
    for (auto c : p.coeffs()) r.push_back(p.field().decode(c));
    return r;
}
Fq64::Elem el(const Fq64&, const Fq64::Elem& e) { return e; }
Fq64::Elem el(const FqTab& K, FqTab::Elem e) { return K.decode(e); }

template <class PA, class PB>
bool mismo_factor(const PA& a, const PB& b) {
    set_seed(1234);
    auto ra = a.factor();
    const BigInt ta = rand_below(BigInt(1) << 64);
    set_seed(1234);
    auto rb = b.factor();
    const BigInt tb = rand_below(BigInt(1) << 64);
    if (ta != tb || el(a.field(), ra.leading) != el(b.field(), rb.leading) || ra.factors.size() != rb.factors.size())
        return false;
    for (std::size_t j = 0; j < ra.factors.size(); ++j)
        if (vec(ra.factors[j].first) != vec(rb.factors[j].first) || ra.factors[j].second != rb.factors[j].second)
            return false;
    return true;
}

static std::vector<BigInt> a_big(const std::vector<int>& v) {
    std::vector<BigInt> r;
    for (int x : v) r.emplace_back(x);
    return r;
}

template <class KB>
void bench_factor(const char* titulo, const std::vector<std::vector<BigInt>>& ent, const Fq64& ka, const KB& kb,
                  const char* nombre) {
    std::vector<Poly<Fq64>> A;
    std::vector<Poly<KB>> B;
    bool ok = true;
    for (const auto& c : ent) {
        A.push_back(poly_de(ka, c));
        B.push_back(poly_de(kb, c));
        ok = ok && mismo_factor(A.back(), B.back());
    }
    auto [ta, tb] = medir_ab([&](int i) { sumidero += A[i].factor().factors.size(); },
                             [&](int i) { sumidero += B[i].factor().factors.size(); }, static_cast<int>(ent.size()));
    std::printf("\n== factor %s (%zu entradas)\n  Fq64        %10.1f us\n  %-10s  %10.1f us  (%.2fx)  idéntico (factores, orden y generador): %s\n",
                titulo, ent.size(), ta, nombre, tb, ta / tb, ok ? "sí" : "NO");
}

template <class KB>
void bench_mul(const char* titulo, const std::vector<std::vector<BigInt>>& ent, const Fq64& ka, const KB& kb,
               const char* nombre) {
    std::vector<Poly<Fq64>> A1, A2;
    std::vector<Poly<KB>> B1, B2;
    bool ok = true;
    for (std::size_t i = 0; i + 1 < ent.size(); i += 2) {
        A1.push_back(poly_de(ka, ent[i]));
        A2.push_back(poly_de(ka, ent[i + 1]));
        B1.push_back(poly_de(kb, ent[i]));
        B2.push_back(poly_de(kb, ent[i + 1]));
        ok = ok && vec(A1.back() * A2.back()) == vec(B1.back() * B2.back());
    }
    const int n = static_cast<int>(A1.size());
    auto [ta, tb] = medir_ab([&](int i) { sumidero += (A1[i] * A2[i]).degree(); },
                             [&](int i) { sumidero += (B1[i] * B2[i]).degree(); }, n);
    std::printf("\n== producto %s (%d pares)\n  Fq64        %10.1f us\n  %-10s  %10.1f us  (%.2fx)  idéntico: %s\n", titulo, n,
                ta, nombre, tb, ta / tb, ok ? "sí" : "NO");
}

int main() {
    auto ent = leer_entradas("bench/results/fase6/factor_entradas.txt");
    if (ent.empty()) {
        std::printf("No se encuentra bench/results/fase6/factor_entradas.txt\n");
        return 1;
    }
    const std::vector<int> mod28 = {1, 1, 0, 1, 1, 0, 0, 0, 1};
    const std::vector<int> mod330 = {2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0,
                                     2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1};
    const Fq64 k256(Fp64(BigInt(2)), a_big(mod28));
    const FqTab t256(Fp64(BigInt(2)), a_big(mod28));
    const Fq64 k3(Fp64(BigInt(3)), a_big(mod330));
    const FqSwar s3(Fp64(BigInt(3)), a_big(mod330));
    std::printf("Fase 7, propuesta 7: Fq64 frente a FqTab / FqSwar (un hilo, mejor de 3 tandas de >= 30 ms\n"
                "con vueltas completas a las entradas, A/B intercalado)\n");
    bench_factor("F_256 grado 30", ent["f2e8_d30"], k256, t256, "FqTab");
    bench_mul("F_256 grado 1000", ent["f2e8_mul_d1000"], k256, t256, "FqTab");
    bench_factor("F_(3^30) grado 5", ent["f3e30_d5"], k3, s3, "FqSwar");
    bench_mul("F_(3^30) grado 200", ent["f3e30_mul_d200"], k3, s3, "FqSwar");
    std::printf("\n(sumidero %llu)\n", static_cast<unsigned long long>(sumidero));
    return 0;
}
