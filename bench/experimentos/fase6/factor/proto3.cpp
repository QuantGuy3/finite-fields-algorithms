// proto3.cpp — Fase 6, tema «factor»: F_p con p pequeño (F_101).
//
//   G) Fp64S: Fp64 con acumulador perezoso de 64 bits (Acc = uint64_t, s += a·b
//      sin acarreo ni comprobación) y reducción final con Barrett de 64 bits.
//      Válido si (p-1)^2 · (número de términos) < 2^64; con p < 2^16 hacen
//      falta más de 2^32 términos para desbordar (polinomios de > 2^32
//      coeficientes).  Lo usan mul_school, divmod_school (Reducer) y la
//      composición de Brent-Kung de la librería.
//   H) Fp64RS: G + la mul_primo rápida de proto.cpp (desempaquetado directo).
//
// Se comparan con la LIBRERÍA (Poly<Fp64>): factores, multiplicidades y estado
// del generador tras factor().
//
// Compilación (desde la raíz del repo, portable como el módulo):
//   g++ -std=c++23 -O3 -Iinclude bench/experimentos/fase6/factor/proto3.cpp -lgmpxx -lgmp -static -o build/fase6/factor/proto3.exe
// Ejecución:  ALCP_THREADS=1 build/fase6/factor/proto3.exe
#define PROTO_SIN_MAIN
#include "proto.cpp"

struct Fp64S : Fp64 {
    using Acc = u64;
    explicit Fp64S(BigInt p) : Fp64(p), red_(modulus64()) {
        if (modulus64() >= (u64(1) << 16)) throw std::invalid_argument("Fp64S: prototipo solo para p < 2^16");
    }
    void acc_clear(Acc& s) const { s = 0; }
    void acc_set(Acc& s, const Elem& a) const { s = a; }
    void acc_fma(Acc& s, const Elem& a, const Elem& b) const { s += a * b; }
    Elem acc_reduce(Acc& s) const { return red_(s); }

private:
    Barrett64 red_;
};
struct Fp64RS : Fp64S {
    using Fp64S::Fp64S;
};
namespace alcp::ks {
template <>
std::vector<u64> mul_primo<Fp64RS>(const Fp64RS& K, const u64* a, std::size_t na, const u64* b, std::size_t nb) {
    return proto::mul_primo_rapida(K, a, na, b, nb);
}
}  // namespace alcp::ks

int main() {
    auto ent = leer_entradas("bench/results/fase6/factor_entradas.txt");
    if (ent.empty()) {
        std::printf("No se encuentra bench/results/fase6/factor_entradas.txt (ejecute flint_ref.py)\n");
        return 1;
    }
    set_seed(6);
    Fp64 f101(101);
    Fp64R f101r(101);
    Fp64S f101s(101);
    Fp64RS f101rs(101);
    std::printf("Prototipos fase 6 (3), tema factor: F_101 (un hilo, mejor de 3 tandas de >= 30 ms, 16 entradas rotando)\n");
    bench_factor("F_101 grado 40 (producto)", ent["fp101_d40"], f101, std::pair{"mul_primo rápida", &f101r},
                 std::pair{"Acc de 64 bits", &f101s}, std::pair{"las dos", &f101rs});
    bench_factor("F_101 grado 100", ent["fp101_d100"], f101, std::pair{"mul_primo rápida", &f101r},
                 std::pair{"Acc de 64 bits", &f101s}, std::pair{"las dos", &f101rs});
    // Productos sueltos de F_101 por tamaños (Kronecker desde 32 coeficientes).
    std::printf("\n== producto en F_101 (16 pares aleatorios por tamaño)\n");
    for (int n : {16, 40, 100, 300, 1000}) {
        std::vector<Poly<Fp64>> A, B;
        std::vector<Poly<Fp64RS>> C, D;
        for (int i = 0; i < 16; ++i) {
            A.push_back(Poly<Fp64>::random(f101, n));
            B.push_back(Poly<Fp64>::random(f101, n));
            C.emplace_back(f101rs, A.back().coeffs());
            D.emplace_back(f101rs, B.back().coeffs());
        }
        bool ok = true;
        for (int i = 0; i < 16; ++i) ok = ok && (A[i] * B[i]).coeffs() == (C[i] * D[i]).coeffs();
        const double t1 = medir_us([&](int i) { sumidero += (A[i] * B[i]).degree(); });
        const double t2 = medir_us([&](int i) { sumidero += (C[i] * D[i]).degree(); });
        std::printf("  grado %5d: librería %8.2f us | prototipo %8.2f us (%.2fx) | idéntico: %s\n", n, t1, t2, t1 / t2,
                    ok ? "sí" : "NO");
    }
    std::printf("\n(sumidero %llu)\n", static_cast<unsigned long long>(sumidero));
    return 0;
}
