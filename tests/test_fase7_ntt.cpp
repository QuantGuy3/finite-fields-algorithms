// test_fase7_ntt.cpp — Fase 7, propuesta 4 («ntt»): la convolución en palabras
// con NTT radix 4 sin permutación (y el producto punto a punto de Montgomery
// con n^{-1} fusionado) da exactamente los mismos coeficientes que antes.
//
//   1. fft::convolve (producto y cuadrado) frente a una copia de la
//      convolución anterior (word::transform + mul_mod + inversa) y frente a
//      mul_school, en longitudes de resultado 2^j - 1, 2^j y 2^j + 1 hasta
//      2^16 (hasta la capacidad FFT del cuerpo), con operandos equilibrados,
//      desequilibrados y de un solo coeficiente, y con coeficientes todos p - 1.
//   2. Cuerpos Fp64: F_998244353, F_7681, F_65537, F_193, F_257, F_577, F_17,
//      F_5, F_3, F_2, F_1004535809 (< 2^30, cerca de la cota del camino AVX2),
//      F_2013265921 (> 2^30: camino escalar) y el primo de 62 bits 4179340454199820289; y Fp (GMP)
//      con 998244353 y el primo de 62 bits (comparten convolve_word).
//   3. fft::transform (API pública, orden natural) sigue siendo la DFT: se
//      compara con la definición en tamaños pequeños e ida y vuelta.
//   4. Camino AVX2 (p < 2^30, n >= 16, si la CPU lo tiene) frente al camino
//      escalar radix 4 forzado y frente a la convolución anterior; las tablas
//      de 32 bits solo existen si el camino se usa.
//   5. Generador aleatorio en frío: el producto por FFT en un cuerpo recién
//      creado consume el generador igual que antes (lo único que lo consume es
//      la primera llamada a fft::tables) y la división (mul_no_rng) en frío
//      no construye tablas.
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <cstdio>
#include <string>
#include <vector>

using namespace alcp;
using u64 = std::uint64_t;
using V = std::vector<u64>;

namespace {

// --------------------------------------------- copia de la convolución anterior
// convolve_word de 80b983d: transformada directa natural (permutación + DIT),
// producto punto a punto con mul_mod y transformada inversa con n^{-1}.
template <class Field>
V conv_anterior(const Field& K, const V& A, const V& B) {
    const std::size_t L = A.size() + B.size() - 1;
    const unsigned long k = fft::ceil_log2(L);
    const std::size_t n = std::size_t(1) << k;
    const auto& t = fft::tables(K, k);
    const u64 p = t.word.p;
    V wa(n, 0), wb(n, 0);
    std::copy(A.begin(), A.end(), wa.begin());
    std::copy(B.begin(), B.end(), wb.begin());
    fft::word::transform(t, wa.data(), n, false);
    fft::word::transform(t, wb.data(), n, false);
    for (std::size_t i = 0; i < n; ++i) wa[i] = fft::word::mul_mod(wa[i], wb[i], p);
    fft::word::transform(t, wa.data(), n, true);
    wa.resize(L);
    return wa;
}

// Camino escalar radix 4 forzado: copia de las tablas sin la parte de 32 bits
// (en una CPU con AVX2 fft::convolve usa el camino vectorial cuando p < 2^30 y
// n >= 16).  B puede ser el mismo objeto que A (cuadrado).
V conv_escalar(const Fp64& K, const V& A, const V& B) {
    const std::size_t L = A.size() + B.size() - 1;
    const unsigned long k = fft::ceil_log2(L);
    FFTTables<u64> t = fft::tables(K, k);
    t.word.w32 = Word32NTTTables{};
    return fft::detail::convolve_word(t, A, B, std::size_t(1) << k, L);
}

bool usa_avx2(const Fp64& K, std::size_t L) {
    return fft::tables(K, fft::ceil_log2(L)).word.w32.p != 0;
}

// Producto escolar directo (referencia independiente de la NTT).
V school(u64 p, const V& A, const V& B) {
    V r(A.size() + B.size() - 1, 0);
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < B.size(); ++j)
            r[i + j] = static_cast<u64>((r[i + j] + static_cast<unsigned __int128>(A[i]) * B[j]) % p);
    return r;
}

V aleatorio(const Fp64& K, std::size_t len, bool todo_pm1) {
    V v(len);
    for (auto& x : v) x = todo_pm1 ? K.modulus64() - 1 : K.random();
    return v;
}

// Longitudes de resultado alrededor de potencias de 2 hasta 2^min(cap, 16).
std::vector<std::size_t> longitudes(const Fp64& K) {
    std::vector<std::size_t> Ls;
    const unsigned long cap = std::min<unsigned long>(K.fft_capacity(), 16);
    for (unsigned long j = 0; j <= cap; ++j) {
        const std::size_t n = std::size_t(1) << j;
        for (std::size_t L : {n - 1, n, n + 1})
            if (L >= 1 && fft::fft_exponent(K, L) >= 0 && (Ls.empty() || Ls.back() < L)) Ls.push_back(L);
    }
    return Ls;
}

// ----------------------------------------------------- 1-2. Fp64: producto y cuadrado
void test_fp64(const Fp64& K, const char* nombre) {
    const u64 p = K.modulus64();
    int pruebas = 0, fallos = 0, vectoriales = 0;
    auto cmp = [&](const V& r, const V& ref, const char* que, std::size_t L) {
        ++pruebas;
        if (r != ref) {
            ++fallos;
            if (fallos <= 5) std::printf("  distinto: %s %s L=%zu\n", nombre, que, L);
        }
    };
    for (std::size_t L : longitudes(K)) {
        for (int modo = 0; modo < 2; ++modo) {  // 0: aleatorio, 1: todo p - 1
            // Repartos de L + 1 = |A| + |B|: equilibrado, 1/3 y un solo coeficiente.
            const std::size_t repartos[] = {(L + 1) / 2, std::max<std::size_t>(1, (L + 1) / 3), 1};
            for (std::size_t la : repartos) {
                const std::size_t lb = L + 1 - la;
                if (lb == 0) continue;
                const V A = aleatorio(K, la, modo), B = aleatorio(K, lb, modo);
                const V r = fft::convolve(K, A, B);
                cmp(r, conv_anterior(K, A, B), "producto/anterior", L);
                cmp(conv_escalar(K, A, B), r, "producto/escalar", L);
                vectoriales += usa_avx2(K, L);
                if (L <= 2048) cmp(r, school(p, A, B), "producto/escolar", L);
            }
            // Cuadrado (mismo objeto: una sola transformada directa).
            if (L % 2 == 1) {
                const V A = aleatorio(K, (L + 1) / 2, modo);
                const V r = fft::convolve(K, A, A);
                const V copia = A;
                cmp(r, conv_anterior(K, A, copia), "cuadrado/anterior", L);
                cmp(r, fft::convolve(K, A, copia), "cuadrado/producto", L);
                cmp(conv_escalar(K, A, A), r, "cuadrado/escalar", L);
                if (L <= 2048) cmp(r, school(p, A, A), "cuadrado/escolar", L);
            }
        }
    }
    // Poly: producto y cuadrado (a*a) frente al escolar de la librería.
    for (int d : {40, 63, 64, 100, 255, 256, 511, 700}) {
        if (fft::fft_exponent(K, 2 * static_cast<std::size_t>(d) + 1) < 0) continue;
        const Poly<Fp64> a = Poly<Fp64>::random(K, d), b = Poly<Fp64>::random(K, d / 2 + 1);
        ++pruebas;
        if (a * b != a.mul_school(b) || a * a != a.mul_school(a)) ++fallos;
    }
    std::printf("  %s: %d comparaciones, %d distintas (%d productos por AVX2)\n", nombre, pruebas, fallos,
                vectoriales);
    CHECK_MSG(fallos == 0, (std::string("convolve Fp64 ") + nombre).c_str());
}

// ----------------------------------------------------- 2. Fp (GMP), mismo convolve_word
void test_fp(const char* primo) {
    const Fp K{BigInt(primo)};
    const Fp64 K64{BigInt(primo)};
    int fallos = 0, pruebas = 0;
    for (std::size_t L : {std::size_t(63), std::size_t(64), std::size_t(65), std::size_t(1023),
                          std::size_t(2049), std::size_t(32768)}) {
        for (int modo = 0; modo < 2; ++modo) {
            const std::size_t la = (L + 1) / 2, lb = L + 1 - la;
            std::vector<BigInt> A(la), B(lb);
            V a(la), b(lb);
            for (std::size_t i = 0; i < la; ++i) { A[i] = modo ? K.characteristic() - 1 : K.random(); a[i] = fft::word::to_word(A[i]); }
            for (std::size_t i = 0; i < lb; ++i) { B[i] = modo ? K.characteristic() - 1 : K.random(); b[i] = fft::word::to_word(B[i]); }
            const V ref = conv_anterior(K64, a, b);
            const V ref2 = conv_anterior(K64, a, a);
            auto words = [](const std::vector<BigInt>& r) {
                V w(r.size());
                for (std::size_t i = 0; i < r.size(); ++i) w[i] = fft::word::to_word(r[i]);
                return w;
            };
            pruebas += 2;
            if (words(fft::convolve(K, A, B)) != ref) ++fallos;
            if (words(fft::convolve(K, A, A)) != ref2) ++fallos;
        }
    }
    std::printf("  Fp %s: %d comparaciones, %d distintas\n", primo, pruebas, fallos);
    CHECK_MSG(fallos == 0, (std::string("convolve Fp ") + primo).c_str());
}

// ----------------------------------------------------- 3. fft::transform intacta
void test_transform(const Fp64& K, const char* nombre) {
    bool ok = true;
    for (unsigned k = 0; k <= std::min<unsigned long>(K.fft_capacity(), 6); ++k) {
        const std::size_t n = std::size_t(1) << k;
        const V a = aleatorio(K, n, false);
        // DFT por definición con la raíz de las tablas: g = fwd[1] (n >= 4) o -1.
        const auto& t = fft::tables(K, k);
        const u64 g = n >= 4 ? t.fwd[1] : (n == 2 ? K.modulus64() - 1 : 1);
        V ref(n, 0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                ref[i] = K.add(ref[i], K.mul(a[j], K.pow(g, BigInt(static_cast<unsigned long>(i * j % n)))));
        V b = a;
        fft::transform(K, b, false);
        ok = ok && b == ref;
        fft::transform(K, b, true);
        ok = ok && b == a;
    }
    CHECK_MSG(ok, (std::string("fft::transform natural ") + nombre).c_str());
}

// ----------------------------------------------------- 4. tablas de 32 bits
void test_tablas32() {
    const bool avx2 = fft::word::avx2_disponible();
    std::printf("  AVX2 disponible y activo: %s\n", avx2 ? "sí" : "no");
    bool ok = true;
    for (const char* q : {"998244353", "7681", "65537", "4179340454199820289", "1004535809", "2013265921"}) {
        const Fp64 K{BigInt(q)};
        const bool pequeno = K.modulus64() < (u64(1) << 30);
        for (unsigned k = 0; k <= std::min<unsigned long>(K.fft_capacity(), 6); ++k) {
            const auto& W = fft::tables(K, k).word;
            const bool espera = avx2 && pequeno && k >= 4;
            ok = ok && ((W.w32.p != 0) == espera);
            if (W.w32.p != 0) {  // mismos twiddles que las tablas de 64 bits
                for (std::size_t i = 0; i < W.fwd.size(); ++i)
                    ok = ok && W.w32.fwd[i] == W.fwd[i] && W.w32.inv[i] == W.inv[i];
            }
        }
    }
    CHECK_MSG(ok, "tablas de 32 bits: solo con AVX2, p < 2^30 y n >= 16; mismos twiddles");
}

// ----------------------------------------------------- 5. generador en frío
void test_generador() {
    for (unsigned long p : {998244353UL, 193UL, 257UL, 577UL, 101UL}) {
        // a) Producto por FFT en un cuerpo recién creado frente a pedir las
        //    tablas a mano en otro recién creado: mismo consumo del generador.
        set_seed(4242);
        const Fp64 K1(p);
        const Poly<Fp64> a(K1, aleatorio(K1, 150, false)), b(K1, aleatorio(K1, 150, false));
        const Poly<Fp64> c = a * b;
        const BigInt r1 = rand_below(BigInt(1000000007));

        set_seed(4242);
        const Fp64 K2(p);
        const Poly<Fp64> a2(K2, aleatorio(K2, 150, false)), b2(K2, aleatorio(K2, 150, false));
        if (fft::fft_exponent(K2, 299) >= 0) (void)fft::tables(K2, fft::ceil_log2(299));
        const Poly<Fp64> c2 = a2.mul_school(b2);
        const BigInt r2 = rand_below(BigInt(1000000007));
        CHECK_MSG(r1 == r2 && c.coeffs() == c2.coeffs(),
                  (std::string("generador en frío, producto, p = ") + std::to_string(p)).c_str());

        // b) Las rutas «sin generador» (división, que usa mul_no_rng) en frío
        //    no construyen tablas ni consumen el generador.
        set_seed(99);
        const Fp64 K3(p);
        const Poly<Fp64> x(K3, aleatorio(K3, 600, false)), y(K3, aleatorio(K3, 300, false));
        const BigInt s0 = rand_below(BigInt(1000000007));
        set_seed(99);
        const Fp64 K4(p);
        const Poly<Fp64> x4(K4, aleatorio(K4, 600, false)), y4(K4, aleatorio(K4, 300, false));
        const auto qr = x4.divmod(y4);
        const BigInt s1 = rand_below(BigInt(1000000007));
        CHECK_MSG(s0 == s1 && K4.fft_cache_.empty() && qr == x4.divmod_school(y4),
                  (std::string("generador en frío, divmod, p = ") + std::to_string(p)).c_str());
    }
}

}  // namespace

int main() {
    set_seed(20261108);
    const char* primos[] = {"998244353", "7681", "65537", "193", "257", "577", "17", "5", "3", "2", "1004535809", "2013265921",
                            "4179340454199820289"};
    for (const char* q : primos) {
        const Fp64 K{BigInt(q)};
        test_fp64(K, q);
        test_transform(K, q);
    }
    test_fp("998244353");
    test_fp("4179340454199820289");
    test_tablas32();
    test_generador();
    return check::summary("fase7_ntt");
}
