// Fase 6 (investigación), tema «ntt»: coste por partes de la convolución NTT en
// palabras de alcp y prototipos de variantes más rápidas, FUERA de la librería.
//
// Compilación (portable, como el módulo; el AVX2 se elige en tiempo de ejecución):
//   g++ -std=c++23 -O3 -Iinclude bench/experimentos/fase6/ntt/ntt_fase6.cpp -lgmpxx -lgmp -static -o build/fase6/ntt/ntt_fase6.exe
// Ejecución (desde la raíz del repositorio):
//   ALCP_THREADS=1 build/fase6/ntt/ntt_fase6_c.exe   (el .exe se reenlazó con -DRELINK_3 por Smart App Control)
//
// Qué hace:
//   1) Perfil por partes de fft::convolve (camino word::transform) para
//      longitudes de resultado ~2000 (n = 2^11) y ~20000 (n = 2^15) en
//      F_998244353: reservas y conversiones, permutación (t.swaps), etapas
//      half = 1, half = 2, resto de etapas, normalizaciones, producto punto a
//      punto (mul_mod con %) y escalado por n^{-1}.
//   2) Variantes, todas comprobadas contra fft::convolve de la LIBRERÍA
//      (resultado idéntico coeficiente a coeficiente):
//        A  : directa DIF sin permutación + inversa DIT sin permutación
//             (el producto punto a punto no depende del orden).
//        A+M: A con producto punto a punto de Montgomery fusionado con n^{-1}
//             (desaparece la pasada final de escalado).
//        B4 : radix 4 (dos etapas fusionadas en registros) + M.
//        B4e: B4 con la última/primera pareja de etapas (len = 4) especializada
//             sin multiplicar por los twiddles iguales a 1.
//        AVX: p < 2^30, 8 carriles de 32 bits, Shoup en las mariposas y
//             Montgomery en el producto, con __attribute__((target("avx2")))
//             y detección con __builtin_cpu_supports("avx2").
//   3) Tiempo del producto Poly<Fp64> de grado 1000 y 10000 de la librería
//      (para estimar la ganancia de extremo a extremo).
//
// Metodología: 16 entradas distintas rotando, mejor de 3 tandas de >= 30 ms,
// checksum de los resultados en una variable volátil.
#include "alcp/alcp.hpp"

#include <immintrin.h>

#include <algorithm>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <vector>

using namespace alcp;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using u128 = unsigned __int128;
using reloj = std::chrono::steady_clock;

static volatile u64 sumidero = 0;
static constexpr int kEntradas = 16;

// Mejor de 3 tandas de >= 30 ms; f(idx) con idx rotando en [0, 16); devuelve
// nanosegundos por llamada.
template <class F>
static double medir(F&& f) {
    double mejor = 1e300;
    for (int tanda = 0; tanda < 3; ++tanda) {
        long reps = 0;
        u64 acc = 0;
        double dt = 0;
        const auto t0 = reloj::now();
        do {
            acc += f(static_cast<int>(reps % kEntradas));
            ++reps;
            if ((reps & 3) == 0)
                dt = std::chrono::duration<double, std::nano>(reloj::now() - t0).count();
        } while (dt < 30e6);
        dt = std::chrono::duration<double, std::nano>(reloj::now() - t0).count();
        sumidero = sumidero + acc;
        mejor = std::min(mejor, dt / reps);
    }
    return mejor;
}

// ---------------------------------------------------------------------------
// Utilidades escalares (mismas fórmulas que fft.hpp)
// ---------------------------------------------------------------------------
static inline u64 mulhi(u64 a, u64 b) { return static_cast<u64>((static_cast<u128>(a) * b) >> 64); }
// w·y con Shoup SIN corrección final: resultado en [0, 2p) para cualquier y < 2^64.
static inline u64 shoup_lazy(u64 w, u64 ws, u64 y, u64 p) { return w * y - mulhi(ws, y) * p; }
// x en [0, 2m) -> [0, m) sin saltos.
static inline u64 red(u64 x, u64 m) { return x - (m & (0 - static_cast<u64>(x >= m))); }

// ---------------------------------------------------------------------------
// Copia instrumentable de word::transform (por trozos, para el perfil)
// ---------------------------------------------------------------------------
// Etapas DIT con half en [hlo, hhi] (copia literal del bucle de la librería).
static void lib_etapas(const WordNTTTables& W, u64* a, std::size_t n, bool inv, std::size_t hlo,
                       std::size_t hhi) {
    const u64 p = W.p, p2 = 2 * p;
    const u64* tw_all = inv ? W.inv.data() : W.fwd.data();
    const u64* ts_all = inv ? W.inv_sh.data() : W.fwd_sh.data();
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len >> 1;
        if (half < hlo || half > hhi) continue;
        const u64* tw = tw_all + (half - 1);
        const u64* ts = ts_all + (half - 1);
        for (std::size_t i = 0; i < n; i += len) {
            u64* x = a + i;
            u64* y = a + i + half;
            for (std::size_t j = 0; j < half; ++j) {
                u64 u = x[j];
                u -= p2 & (0 - static_cast<u64>(u >= p2));
                const u64 q = mulhi(ts[j], y[j]);
                const u64 v = tw[j] * y[j] - q * p;
                x[j] = u + v;
                y[j] = u - v + p2;
            }
        }
    }
}
template <class T>
static void lib_swaps(const T& t, u64* a) {
    for (const auto& [i, j] : t.swaps) std::swap(a[i], a[j]);
}
static void lib_norm_fwd(const WordNTTTables& W, u64* a, std::size_t n) {
    const u64 p = W.p, p2 = 2 * p;
    for (std::size_t i = 0; i < n; ++i) a[i] = red(red(a[i], p2), p);
}
static void lib_escala_inv(const WordNTTTables& W, u64* a, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) a[i] = fft::word::mul_shoup(W.inv_n, W.inv_n_sh, a[i], W.p);
}
static void lib_punto(const WordNTTTables& W, u64* a, const u64* b, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) a[i] = fft::word::mul_mod(a[i], b[i], W.p);
}

// Copia completa de convolve_word (para comprobar que la copia cuesta lo mismo).
template <class T>
static std::vector<u64> conv_copia(const T& t, const std::vector<u64>& A, const std::vector<u64>& B,
                                   std::size_t n, std::size_t L) {
    std::vector<u64> wa(n, 0), wb(n, 0);
    std::copy(A.begin(), A.end(), wa.begin());
    std::copy(B.begin(), B.end(), wb.begin());
    fft::word::transform(t, wa.data(), n, false);
    fft::word::transform(t, wb.data(), n, false);
    lib_punto(t.word, wa.data(), wb.data(), n);
    fft::word::transform(t, wa.data(), n, true);
    std::vector<u64> out(L);
    std::copy(wa.begin(), wa.begin() + L, out.begin());
    return out;
}

// ---------------------------------------------------------------------------
// Variante A: DIF radix 2 sin permutación (entrada natural, salida en orden de
// bits invertido; valores en [0, 2p) a la salida, Harvey DIF perezosa).
// ---------------------------------------------------------------------------
static void dif_r2(const WordNTTTables& W, u64* a, std::size_t n, std::size_t len_min = 2) {
    const u64 p = W.p, p2 = 2 * p;
    for (std::size_t len = n; len >= len_min; len >>= 1) {
        const std::size_t half = len >> 1;
        const u64* tw = W.fwd.data() + (half - 1);
        const u64* ts = W.fwd_sh.data() + (half - 1);
        for (std::size_t i = 0; i < n; i += len) {
            u64* x = a + i;
            u64* y = a + i + half;
            for (std::size_t j = 0; j < half; ++j) {
                const u64 u = x[j], v = y[j];
                x[j] = red(u + v, p2);                        // [0, 2p)
                y[j] = shoup_lazy(tw[j], ts[j], u - v + p2, p);  // [0, 2p)
            }
        }
    }
}
// Inversa DIT sin permutación (entrada en orden de bits invertido): las etapas
// de la librería tal cual; valores en [0, 4p).
static void dit_r2(const WordNTTTables& W, u64* a, std::size_t n) {
    lib_etapas(W, a, n, true, 1, n);
}

// ---------------------------------------------------------------------------
// Producto punto a punto de Montgomery fusionado con n^{-1}
// ---------------------------------------------------------------------------
// REDC(a·b) = a·b·2^{-64} en [0, 2p); luego ·(n^{-1}·2^64 mod p) con Shoup
// perezoso -> a·b·n^{-1} en [0, 2p).  Así desaparece la pasada de escalado
// tras la inversa (solo queda normalizar [0, 4p) -> [0, p)).  Vale para
// p < 2^62 impar con a, b < 2p.
struct Mont {
    u64 p, pinv;  // pinv = -p^{-1} mod 2^64
    u64 c, cs;    // c = n^{-1}·2^64 mod p y su cociente de Shoup
};
static Mont hacer_mont(const WordNTTTables& W) {
    Mont m{};
    m.p = W.p;
    u64 x = W.p;
    for (int i = 0; i < 6; ++i) x *= 2 - W.p * x;  // p^{-1} mod 2^64
    m.pinv = 0 - x;
    const u64 r = static_cast<u64>((~u128(0)) % W.p + 1) % W.p;  // 2^128 mod p
    // c = inv_n · 2^64 mod p
    m.c = static_cast<u64>((static_cast<u128>(W.inv_n) << 64) % W.p);
    (void)r;
    m.cs = fft::word::shoup(m.c, W.p);
    return m;
}
static void punto_mont(const Mont& M, u64* a, const u64* b, std::size_t n) {
    const u64 p = M.p;
    for (std::size_t i = 0; i < n; ++i) {
        const u128 t = static_cast<u128>(a[i]) * b[i];
        const u64 m = static_cast<u64>(t) * M.pinv;
        const u64 u = static_cast<u64>((t + static_cast<u128>(m) * p) >> 64);  // [0, 2p)
        a[i] = shoup_lazy(M.c, M.cs, u, p);                                  // [0, 2p)
    }
}
static void norm_4p(u64* a, std::size_t n, u64 p) {
    for (std::size_t i = 0; i < n; ++i) a[i] = red(red(a[i], 2 * p), p);
}

// ---------------------------------------------------------------------------
// Variante B: radix 4 (dos etapas por pasada, en registros)
// ---------------------------------------------------------------------------
template <bool ESP4>
static void dif_r4(const WordNTTTables& W, u64* a, std::size_t n) {
    const u64 p = W.p, p2 = 2 * p;
    const u64* TW = W.fwd.data();
    const u64* TS = W.fwd_sh.data();
    std::size_t len = n;
    if (std::countr_zero(n) & 1) {  // k impar: una etapa radix 2 al principio
        dif_r2(W, a, n, n);
        len >>= 1;
    }
    for (; len >= 4; len >>= 2) {
        const std::size_t q = len >> 2;
        const u64 *w1 = TW + (2 * q - 1), *s1 = TS + (2 * q - 1);  // etapa len (half 2q)
        const u64 *w2 = TW + (q - 1), *s2 = TS + (q - 1);          // etapa len/2 (half q)
        if (ESP4 && len == 4) {
            // w_4^0 = w_2^0 = 1: una sola multiplicación por cada 4 valores.
            const u64 w4 = w1[1], w4s = s1[1];
            for (std::size_t i = 0; i < n; i += 4) {
                const u64 a0 = a[i], a1 = a[i + 1], a2 = a[i + 2], a3 = a[i + 3];
                const u64 b0 = red(a0 + a2, p2), b2 = red(a0 - a2 + p2, p2);
                const u64 b1 = red(a1 + a3, p2), b3 = shoup_lazy(w4, w4s, a1 - a3 + p2, p);
                a[i] = red(b0 + b1, p2);
                a[i + 1] = red(b0 - b1 + p2, p2);
                a[i + 2] = red(b2 + b3, p2);
                a[i + 3] = red(b2 - b3 + p2, p2);
            }
            continue;
        }
        for (std::size_t i = 0; i < n; i += len) {
            u64* x = a + i;
            for (std::size_t j = 0; j < q; ++j) {
                const u64 a0 = x[j], a1 = x[j + q], a2 = x[j + 2 * q], a3 = x[j + 3 * q];
                const u64 b0 = red(a0 + a2, p2), b2 = shoup_lazy(w1[j], s1[j], a0 - a2 + p2, p);
                const u64 b1 = red(a1 + a3, p2), b3 = shoup_lazy(w1[j + q], s1[j + q], a1 - a3 + p2, p);
                x[j] = red(b0 + b1, p2);
                x[j + q] = shoup_lazy(w2[j], s2[j], b0 - b1 + p2, p);
                x[j + 2 * q] = red(b2 + b3, p2);
                x[j + 3 * q] = shoup_lazy(w2[j], s2[j], b2 - b3 + p2, p);
            }
        }
    }
}
template <bool ESP4>
static void dit_r4(const WordNTTTables& W, u64* a, std::size_t n) {
    const u64 p = W.p, p2 = 2 * p;
    const u64* TW = W.inv.data();
    const u64* TS = W.inv_sh.data();
    std::size_t len = 4;
    for (; len <= n; len <<= 2) {
        const std::size_t q = len >> 2;
        const u64 *w1 = TW + (2 * q - 1), *s1 = TS + (2 * q - 1);  // half 2q
        const u64 *w2 = TW + (q - 1), *s2 = TS + (q - 1);          // half q
        if (ESP4 && len == 4) {
            const u64 w4 = w1[1], w4s = s1[1];
            for (std::size_t i = 0; i < n; i += 4) {
                const u64 a0 = red(a[i], p2), a1 = red(a[i + 1], p2);
                const u64 a2 = red(a[i + 2], p2), a3 = red(a[i + 3], p2);
                const u64 b0 = a0 + a1, b1 = a0 - a1 + p2;  // [0, 4p)
                const u64 b2 = a2 + a3, b3 = a2 - a3 + p2;
                const u64 u0 = red(b0, p2), v2 = red(b2, p2);
                a[i] = u0 + v2;
                a[i + 2] = u0 - v2 + p2;
                const u64 u1 = red(b1, p2), v3 = shoup_lazy(w4, w4s, b3, p);
                a[i + 1] = u1 + v3;
                a[i + 3] = u1 - v3 + p2;
            }
            continue;
        }
        for (std::size_t i = 0; i < n; i += len) {
            u64* x = a + i;
            for (std::size_t j = 0; j < q; ++j) {
                u64 a0 = red(x[j], p2), a2 = red(x[j + 2 * q], p2);
                u64 v = shoup_lazy(w2[j], s2[j], x[j + q], p);
                const u64 b0 = a0 + v, b1 = a0 - v + p2;
                v = shoup_lazy(w2[j], s2[j], x[j + 3 * q], p);
                const u64 b2 = a2 + v, b3 = a2 - v + p2;
                const u64 u0 = red(b0, p2), u1 = red(b1, p2);
                v = shoup_lazy(w1[j], s1[j], b2, p);
                x[j] = u0 + v;
                x[j + 2 * q] = u0 - v + p2;
                v = shoup_lazy(w1[j + q], s1[j + q], b3, p);
                x[j + q] = u1 + v;
                x[j + 3 * q] = u1 - v + p2;
            }
        }
    }
    if (std::countr_zero(n) & 1) lib_etapas(W, a, n, true, n / 2, n / 2);  // última etapa radix 2
}

// Convolución con variantes escalares.  modo: 0 = A (DIF+DIT, punto a punto
// con mul_mod tras normalizar, escalado final), 1 = A+M, 2 = B4+M, 3 = B4e+M.
template <class T>
static std::vector<u64> conv_var(const T& t, const Mont& M, int modo, const std::vector<u64>& A,
                                 const std::vector<u64>& B, std::size_t n, std::size_t L) {
    const WordNTTTables& W = t.word;
    const u64 p = W.p;
    std::vector<u64> wa(n, 0), wb(n, 0);
    std::copy(A.begin(), A.end(), wa.begin());
    std::copy(B.begin(), B.end(), wb.begin());
    if (modo <= 1) {
        dif_r2(W, wa.data(), n);
        dif_r2(W, wb.data(), n);
    } else if (modo == 2) {
        dif_r4<false>(W, wa.data(), n);
        dif_r4<false>(W, wb.data(), n);
    } else {
        dif_r4<true>(W, wa.data(), n);
        dif_r4<true>(W, wb.data(), n);
    }
    if (modo == 0) {
        for (std::size_t i = 0; i < n; ++i)
            wa[i] = fft::word::mul_mod(red(wa[i], p), red(wb[i], p), p);
        dit_r2(W, wa.data(), n);
        lib_escala_inv(W, wa.data(), n);
    } else {
        punto_mont(M, wa.data(), wb.data(), n);
        if (modo == 1) dit_r2(W, wa.data(), n);
        else if (modo == 2) dit_r4<false>(W, wa.data(), n);
        else dit_r4<true>(W, wa.data(), n);
        norm_4p(wa.data(), n, p);
    }
    std::vector<u64> out(L);
    std::copy(wa.begin(), wa.begin() + L, out.begin());
    return out;
}

// ---------------------------------------------------------------------------
// Variante AVX2 (p < 2^30): 8 carriles de 32 bits
// ---------------------------------------------------------------------------
// Con p < 2^30 los valores perezosos [0, 4p) caben en 32 bits.  Mariposas de
// Shoup con cociente de 32 bits: ws = floor(w·2^32 / p).
struct T32 {
    u32 p = 0, pinv = 0, c = 0, cs = 0;  // pinv = -p^{-1} mod 2^32; c = n^{-1}·2^32 mod p
    std::vector<u32> fw, fws, iw, iws;   // mismos twiddles que la librería (n - 1)
    alignas(32) u32 f4[8], f4s[8], f2[8], f2s[8], i4[8], i4s[8], i2[8], i2s[8];
};
static T32 hacer_t32(const WordNTTTables& W) {
    T32 T;
    T.p = static_cast<u32>(W.p);
    auto sh = [&](u64 w) { return static_cast<u32>((w << 32) / W.p); };
    for (std::size_t i = 0; i < W.fwd.size(); ++i) {
        T.fw.push_back(static_cast<u32>(W.fwd[i]));
        T.fws.push_back(sh(W.fwd[i]));
        T.iw.push_back(static_cast<u32>(W.inv[i]));
        T.iws.push_back(sh(W.inv[i]));
    }
    // Vectores de twiddles para las etapas half = 4 ([w0..w3] dos veces) y
    // half = 2 ([w0, w1] cuatro veces); la etapa half = 4 empieza en 3 y la
    // half = 2 en 1.
    for (int l = 0; l < 8; ++l) {
        T.f4[l] = T.fw[3 + (l & 3)]; T.f4s[l] = T.fws[3 + (l & 3)];
        T.f2[l] = T.fw[1 + (l & 1)]; T.f2s[l] = T.fws[1 + (l & 1)];
        T.i4[l] = T.iw[3 + (l & 3)]; T.i4s[l] = T.iws[3 + (l & 3)];
        T.i2[l] = T.iw[1 + (l & 1)]; T.i2s[l] = T.iws[1 + (l & 1)];
    }
    u32 x = T.p;
    for (int i = 0; i < 5; ++i) x *= 2 - T.p * x;
    T.pinv = 0u - x;
    T.c = static_cast<u32>((W.inv_n << 32) % W.p);
    T.cs = sh(T.c);
    return T;
}

#define AVX2 __attribute__((target("avx2")))
using V = __m256i;
AVX2 static inline V ld(const u32* q) { return _mm256_loadu_si256(reinterpret_cast<const V*>(q)); }
AVX2 static inline void st(u32* q, V v) { _mm256_storeu_si256(reinterpret_cast<V*>(q), v); }
// x en [0, 2m) -> [0, m): min(x, x - m) sin signo.
AVX2 static inline V vred(V x, V m) { return _mm256_min_epu32(x, _mm256_sub_epi32(x, m)); }
// w·y con Shoup de 32 bits, en [0, 2p).
AVX2 static inline V vshoup(V w, V ws, V y, V P) {
    const V pe = _mm256_mul_epu32(ws, y);
    const V po = _mm256_mul_epu32(_mm256_srli_epi64(ws, 32), _mm256_srli_epi64(y, 32));
    const V q = _mm256_blend_epi32(_mm256_srli_epi64(pe, 32), po, 0xAA);
    return _mm256_sub_epi32(_mm256_mullo_epi32(w, y), _mm256_mullo_epi32(q, P));
}
// Mariposa DIF perezosa: entradas y salidas en [0, 2p).
AVX2 static inline void vdif(V& x, V& y, V w, V ws, V P, V P2) {
    const V s = vred(_mm256_add_epi32(x, y), P2);
    const V d = _mm256_add_epi32(_mm256_sub_epi32(x, y), P2);
    x = s;
    y = vshoup(w, ws, d, P);
}
// DIF con twiddle 1 (sin multiplicar).
AVX2 static inline void vdif1(V& x, V& y, V P2) {
    const V s = vred(_mm256_add_epi32(x, y), P2);
    const V d = vred(_mm256_add_epi32(_mm256_sub_epi32(x, y), P2), P2);
    x = s;
    y = d;
}
// Mariposa DIT perezosa: entradas y salidas en [0, 4p).
AVX2 static inline void vdit(V& x, V& y, V w, V ws, V P, V P2) {
    const V u = vred(x, P2);
    const V v = vshoup(w, ws, y, P);
    x = _mm256_add_epi32(u, v);
    y = _mm256_add_epi32(_mm256_sub_epi32(u, v), P2);
}
AVX2 static inline void vdit1(V& x, V& y, V P2) {
    const V u = vred(x, P2);
    const V v = vred(y, P2);
    x = _mm256_add_epi32(u, v);
    y = _mm256_add_epi32(_mm256_sub_epi32(u, v), P2);
}
// Reparto de 16 valores (a0 = 0..7, a1 = 8..15) en pares para half = 4, 2, 1 y
// recomposición (ver comentarios de la investigación).
#define REP4(a0, a1, X, Y) X = _mm256_permute2x128_si256(a0, a1, 0x20); Y = _mm256_permute2x128_si256(a0, a1, 0x31)
#define REP2(a0, a1, X, Y) X = _mm256_unpacklo_epi64(a0, a1); Y = _mm256_unpackhi_epi64(a0, a1)
#define REP1(a0, a1, X, Y)                                                                      \
    X = _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(a0), _mm256_castsi256_ps(a1), 0x88)); \
    Y = _mm256_castps_si256(_mm256_shuffle_ps(_mm256_castsi256_ps(a0), _mm256_castsi256_ps(a1), 0xDD))
#define JUN1(X, Y, a0, a1) a0 = _mm256_unpacklo_epi32(X, Y); a1 = _mm256_unpackhi_epi32(X, Y)

AVX2 static void dif_avx2(const T32& T, u32* a, std::size_t n) {
    const V P = _mm256_set1_epi32(static_cast<int>(T.p)), P2 = _mm256_set1_epi32(static_cast<int>(2 * T.p));
    for (std::size_t len = n; len >= 32; len >>= 1) {
        const std::size_t half = len >> 1;
        const u32* tw = T.fw.data() + (half - 1);
        const u32* ts = T.fws.data() + (half - 1);
        for (std::size_t i = 0; i < n; i += len)
            for (std::size_t j = 0; j < half; j += 8) {
                V x = ld(a + i + j), y = ld(a + i + j + half);
                vdif(x, y, ld(tw + j), ld(ts + j), P, P2);
                st(a + i + j, x);
                st(a + i + j + half, y);
            }
    }
    // Etapas half = 8, 4, 2, 1 fusionadas sobre bloques de 16.
    const V W8 = ld(T.fw.data() + 7), W8s = ld(T.fws.data() + 7);
    const V W4 = ld(T.f4), W4s = ld(T.f4s), W2 = ld(T.f2), W2s = ld(T.f2s);
    for (std::size_t i = 0; i < n; i += 16) {
        V a0 = ld(a + i), a1 = ld(a + i + 8), X, Y;
        vdif(a0, a1, W8, W8s, P, P2);
        REP4(a0, a1, X, Y); vdif(X, Y, W4, W4s, P, P2); REP4(X, Y, a0, a1);
        REP2(a0, a1, X, Y); vdif(X, Y, W2, W2s, P, P2); REP2(X, Y, a0, a1);
        REP1(a0, a1, X, Y); vdif1(X, Y, P2); JUN1(X, Y, a0, a1);
        st(a + i, a0);
        st(a + i + 8, a1);
    }
}
AVX2 static void dit_avx2(const T32& T, u32* a, std::size_t n) {
    const V P = _mm256_set1_epi32(static_cast<int>(T.p)), P2 = _mm256_set1_epi32(static_cast<int>(2 * T.p));
    const V W8 = ld(T.iw.data() + 7), W8s = ld(T.iws.data() + 7);
    const V W4 = ld(T.i4), W4s = ld(T.i4s), W2 = ld(T.i2), W2s = ld(T.i2s);
    for (std::size_t i = 0; i < n; i += 16) {
        V a0 = ld(a + i), a1 = ld(a + i + 8), X, Y;
        REP1(a0, a1, X, Y); vdit1(X, Y, P2); JUN1(X, Y, a0, a1);
        REP2(a0, a1, X, Y); vdit(X, Y, W2, W2s, P, P2); REP2(X, Y, a0, a1);
        REP4(a0, a1, X, Y); vdit(X, Y, W4, W4s, P, P2); REP4(X, Y, a0, a1);
        vdit(a0, a1, W8, W8s, P, P2);
        st(a + i, a0);
        st(a + i + 8, a1);
    }
    for (std::size_t len = 32; len <= n; len <<= 1) {
        const std::size_t half = len >> 1;
        const u32* tw = T.iw.data() + (half - 1);
        const u32* ts = T.iws.data() + (half - 1);
        for (std::size_t i = 0; i < n; i += len)
            for (std::size_t j = 0; j < half; j += 8) {
                V x = ld(a + i + j), y = ld(a + i + j + half);
                vdit(x, y, ld(tw + j), ld(ts + j), P, P2);
                st(a + i + j, x);
                st(a + i + j + half, y);
            }
    }
}
// a = REDC(a·b)·c (Montgomery de 32 bits y Shoup), en [0, 2p); a, b < 2p.
AVX2 static void punto_avx2(const T32& T, u32* a, const u32* b, std::size_t n) {
    const V P = _mm256_set1_epi32(static_cast<int>(T.p)), PI = _mm256_set1_epi32(static_cast<int>(T.pinv));
    const V C = _mm256_set1_epi32(static_cast<int>(T.c)), CS = _mm256_set1_epi32(static_cast<int>(T.cs));
    for (std::size_t i = 0; i < n; i += 8) {
        const V x = ld(a + i), y = ld(b + i);
        const V te = _mm256_mul_epu32(x, y);
        const V me = _mm256_mul_epu32(te, PI);
        const V ue = _mm256_srli_epi64(_mm256_add_epi64(te, _mm256_mul_epu32(me, P)), 32);
        const V to = _mm256_mul_epu32(_mm256_srli_epi64(x, 32), _mm256_srli_epi64(y, 32));
        const V mo = _mm256_mul_epu32(to, PI);
        const V uo = _mm256_add_epi64(to, _mm256_mul_epu32(mo, P));
        const V u = _mm256_blend_epi32(ue, uo, 0xAA);
        st(a + i, vshoup(C, CS, u, P));
    }
}
// Conversión de entrada (u64 -> u32 con ceros), salida normalizada [0, 4p) -> [0, p).
AVX2 static std::vector<u64> conv_avx2(const T32& T, const std::vector<u64>& A, const std::vector<u64>& B,
                                       std::size_t n, std::size_t L) {
    std::vector<u32> wa(n, 0), wb(n, 0);
    for (std::size_t i = 0; i < A.size(); ++i) wa[i] = static_cast<u32>(A[i]);
    for (std::size_t i = 0; i < B.size(); ++i) wb[i] = static_cast<u32>(B[i]);
    dif_avx2(T, wa.data(), n);
    dif_avx2(T, wb.data(), n);
    punto_avx2(T, wa.data(), wb.data(), n);
    dit_avx2(T, wa.data(), n);
    const V P = _mm256_set1_epi32(static_cast<int>(T.p)), P2 = _mm256_set1_epi32(static_cast<int>(2 * T.p));
    for (std::size_t i = 0; i < n; i += 8) st(wa.data() + i, vred(vred(ld(wa.data() + i), P2), P));
    std::vector<u64> out(L);
    for (std::size_t i = 0; i < L; ++i) out[i] = wa[i];
    return out;
}
#undef AVX2

// ---------------------------------------------------------------------------
// Programa
// ---------------------------------------------------------------------------
static u64 suma(const std::vector<u64>& v) { return v.empty() ? 0 : v[0] ^ v[v.size() / 2] ^ v.back(); }

int main() {
    const bool avx2 = __builtin_cpu_supports("avx2");
    std::printf("# fase 6, tema ntt (prototipos fuera de la librería)\n");
    std::printf("# hilos alcp: %u (debe ser 1: ALCP_THREADS=1)   AVX2 en esta CPU: %s\n",
                par::max_threads(), avx2 ? "sí" : "no");
    std::printf("# 16 entradas rotando, mejor de 3 tandas de >= 30 ms\n\n");

    Fp64 K(998244353);
    using P = Poly<Fp64>;

    // ---- Comprobación de igualdad con la librería en muchos tamaños
    // (y también con un primo de 62 bits para las variantes escalares).
    {
        // Primo p = c·2^40 + 1 < 2^62 (para probar Montgomery de 64 bits).
        u64 p62 = 0;
        for (u64 c = (u64(1) << 22) - 1; c > 0; --c) {
            const u64 cand = (c << 40) | 1;
            if (is_prime(BigInt(std::to_string(cand)))) { p62 = cand; break; }
        }
        Fp64 K62(BigInt(std::to_string(p62)));
        int fallos = 0, pruebas = 0;
        for (int campo = 0; campo < 2; ++campo) {
            const Fp64& F = campo ? K62 : K;
            for (unsigned k = 4; k <= 16; ++k) {
                const auto& t = fft::tables(F, k);
                const Mont M = hacer_mont(t.word);
                T32 T;
                const bool usa_avx = avx2 && campo == 0;
                if (usa_avx) T = hacer_t32(t.word);
                const std::size_t n = std::size_t(1) << k;
                for (int rep = 0; rep < 4; ++rep) {
                    const std::size_t la = n / 4 + 1 + (rep * 37) % (n / 4);
                    const std::size_t lb = n / 2 - (rep * 11) % (n / 4);
                    std::vector<u64> A(la), B(lb);
                    for (auto& x : A) x = F.random();
                    for (auto& x : B) x = F.random();
                    if (rep == 3) { std::fill(A.begin(), A.end(), F.modulus64() - 1); }  // extremo
                    const std::size_t L = la + lb - 1;
                    const auto ref = fft::convolve(F, A, B);
                    auto chk = [&](const std::vector<u64>& r, const char* nom) {
                        ++pruebas;
                        if (r != ref) { ++fallos; std::printf("DISTINTO: %s k=%u campo=%d\n", nom, k, campo); }
                    };
                    chk(conv_copia(t, A, B, n, L), "copia");
                    for (int m = 0; m < 4; ++m) chk(conv_var(t, M, m, A, B, n, L), "variante");
                    if (usa_avx && n >= 16) chk(conv_avx2(T, A, B, n, L), "avx2");
                }
            }
        }
        std::printf("# comprobación contra fft::convolve (p = 998244353 y p = %llu): %d pruebas, %d distintas\n\n",
                    static_cast<unsigned long long>(p62), pruebas, fallos);
        if (fallos) return 1;
    }

    // ---- Perfil por partes y variantes, para L ~ 2000 y ~ 20000.
    for (std::size_t m : {std::size_t(100), std::size_t(1000), std::size_t(10000)}) {
        // Producto de dos polinomios de grado m: L = 2m + 1 coeficientes.
        const std::size_t L = 2 * m + 1;
        const unsigned k = static_cast<unsigned>(std::bit_width(L - 1));
        const std::size_t n = std::size_t(1) << k;
        const auto& t = fft::tables(K, k);
        const WordNTTTables& W = t.word;
        const Mont M = hacer_mont(W);
        T32 T;
        if (avx2) T = hacer_t32(W);

        std::vector<std::vector<u64>> As, Bs, buf(kEntradas), bufb(kEntradas);
        std::vector<P> PA, PB;
        for (int e = 0; e < kEntradas; ++e) {
            P a = P::random(K, static_cast<int>(m)), b = P::random(K, static_cast<int>(m));
            As.push_back(a.coeffs());
            Bs.push_back(b.coeffs());
            PA.push_back(a);
            PB.push_back(b);
            buf[e].assign(n, 0);
            bufb[e].assign(n, 0);
            for (std::size_t i = 0; i < n; ++i) { buf[e][i] = K.random(); bufb[e][i] = K.random(); }
        }
        const double mariposas = double(n / 2) * k;
        std::printf("== Producto de grado %zu (L = %zu, n = 2^%u = %zu) ==\n", m, L, k, n);

        const double t_poly = medir([&](int e) { return suma((PA[e] * PB[e]).coeffs()); });
        const double t_lib = medir([&](int e) { return suma(fft::convolve(K, As[e], Bs[e])); });
        const double t_copia = medir([&](int e) { return suma(conv_copia(t, As[e], Bs[e], n, L)); });
        std::printf("Poly<Fp64> a*b (librería)           %10.1f us\n", t_poly / 1e3);
        std::printf("fft::convolve (librería)            %10.1f us\n", t_lib / 1e3);
        std::printf("copia de convolve_word              %10.1f us\n", t_copia / 1e3);

        // --- 1) Perfil por partes (por convolución: 2 directas + 1 inversa).
        const double p_res = medir([&](int e) {
            std::vector<u64> wa(n, 0), wb(n, 0);
            std::copy(As[e].begin(), As[e].end(), wa.begin());
            std::copy(Bs[e].begin(), Bs[e].end(), wb.begin());
            return wa[e] + wb[e + 1];
        });
        const double p_sal = medir([&](int e) {
            std::vector<u64> out(L);
            std::copy(buf[e].begin(), buf[e].begin() + L, out.begin());
            return out[e];
        });
        const double p_swap = medir([&](int e) { lib_swaps(t, buf[e].data()); return buf[e][1]; });
        const double p_h1 = medir([&](int e) { lib_etapas(W, buf[e].data(), n, false, 1, 1); return buf[e][1]; });
        const double p_h2 = medir([&](int e) { lib_etapas(W, buf[e].data(), n, false, 2, 2); return buf[e][1]; });
        const double p_resto = medir([&](int e) { lib_etapas(W, buf[e].data(), n, false, 4, n); return buf[e][1]; });
        const double p_tr = medir([&](int e) { fft::word::transform(t, buf[e].data(), n, false); return buf[e][1]; });
        const double p_norm = medir([&](int e) { lib_norm_fwd(W, buf[e].data(), n); return buf[e][1]; });
        for (auto& v : buf) lib_norm_fwd(W, v.data(), n);  // [0, p) para mul_mod
        for (auto& v : bufb) lib_norm_fwd(W, v.data(), n);
        const double p_punto = medir([&](int e) {
            std::vector<u64>& a = buf[e];
            lib_punto(W, a.data(), bufb[e].data(), n);
            return a[1];
        });
        const double p_esc = medir([&](int e) { lib_escala_inv(W, buf[e].data(), n); return buf[e][1]; });
        const double p_mont = medir([&](int e) {
            punto_mont(M, buf[e].data(), bufb[e].data(), n);
            lib_norm_fwd(W, buf[e].data(), n);  // vuelve a [0, p) para la siguiente vuelta (se resta abajo)
            return buf[e][1];
        }) - p_norm;
        const double suma_partes = p_res + 3 * (p_swap + p_h1 + p_h2 + p_resto) + 2 * p_norm + p_punto + p_esc + p_sal;
        std::printf("-- 1) perfil por partes (us por convolución; x3 = dos directas + una inversa)\n");
        std::printf("  reservas + conversión de entrada  %8.1f us\n", p_res / 1e3);
        std::printf("  permutación t.swaps         x3    %8.1f us  (%.1f us cada una)\n", 3 * p_swap / 1e3, p_swap / 1e3);
        std::printf("  etapa half = 1              x3    %8.1f us  (%.2f ns/mariposa)\n", 3 * p_h1 / 1e3, p_h1 / (n / 2));
        std::printf("  etapa half = 2              x3    %8.1f us  (%.2f ns/mariposa)\n", 3 * p_h2 / 1e3, p_h2 / (n / 2));
        std::printf("  resto de etapas (%2u)        x3    %8.1f us  (%.2f ns/mariposa)\n", k - 2, 3 * p_resto / 1e3,
                    p_resto / (double(n / 2) * (k - 2)));
        std::printf("  normalización final directa x2    %8.1f us\n", 2 * p_norm / 1e3);
        std::printf("  producto punto a punto (mul_mod)  %8.1f us  (%.2f ns/elemento)\n", p_punto / 1e3, p_punto / n);
        std::printf("  escalado por n^-1 (inversa)       %8.1f us\n", p_esc / 1e3);
        std::printf("  vector de salida                  %8.1f us\n", p_sal / 1e3);
        std::printf("  suma de las partes                %8.1f us  (convolve medido: %.1f us)\n", suma_partes / 1e3, t_lib / 1e3);
        std::printf("  word::transform entera            %8.1f us  (%.2f ns/mariposa)\n", p_tr / 1e3, p_tr / mariposas);
        std::printf("  [alternativa] punto a punto Montgomery + n^-1 fusionado: %.1f us (%.2f ns/elemento)\n", p_mont / 1e3, p_mont / n);

        // --- 2) Transformadas directas sueltas.
        const double d_r2 = medir([&](int e) { dif_r2(W, buf[e].data(), n); return buf[e][1]; });
        const double d_r4 = medir([&](int e) { dif_r4<false>(W, buf[e].data(), n); return buf[e][1]; });
        const double d_r4e = medir([&](int e) { dif_r4<true>(W, buf[e].data(), n); return buf[e][1]; });
        const double i_r4e = medir([&](int e) { dit_r4<true>(W, buf[e].data(), n); norm_4p(buf[e].data(), n, W.p); return buf[e][1]; }) - p_norm;
        std::printf("-- 2) transformada directa suelta (us; ns/mariposa)\n");
        std::printf("  librería (swaps + DIT + norm.)    %8.1f us  %.2f\n", p_tr / 1e3, p_tr / mariposas);
        std::printf("  A: DIF radix 2 sin permutación    %8.1f us  %.2f\n", d_r2 / 1e3, d_r2 / mariposas);
        std::printf("  B4: DIF radix 4                   %8.1f us  %.2f\n", d_r4 / 1e3, d_r4 / mariposas);
        std::printf("  B4e: DIF radix 4, len=4 especial  %8.1f us  %.2f\n", d_r4e / 1e3, d_r4e / mariposas);
        std::printf("  B4e inversa DIT (sin norm.)       %8.1f us  %.2f\n", i_r4e / 1e3, i_r4e / mariposas);
        if (avx2) {
            std::vector<std::vector<u32>> b32(kEntradas, std::vector<u32>(n));
            for (int e = 0; e < kEntradas; ++e)
                for (std::size_t i = 0; i < n; ++i) b32[e][i] = static_cast<u32>(buf[e][i] % W.p);
            const double d_av = medir([&](int e) { dif_avx2(T, b32[e].data(), n); return u64(b32[e][1]); });
            const double i_av = medir([&](int e) { dit_avx2(T, b32[e].data(), n); return u64(b32[e][1]); });
            std::printf("  AVX2: DIF 8x32 bits               %8.1f us  %.2f\n", d_av / 1e3, d_av / mariposas);
            std::printf("  AVX2: DIT 8x32 bits (sin norm.)   %8.1f us  %.2f\n", i_av / 1e3, i_av / mariposas);
        }

        // --- Convoluciones completas.
        std::printf("-- 2) convolución completa (us; aceleración frente a fft::convolve)\n");
        const char* nombres[4] = {"A   DIF+DIT sin permutación      ", "A+M + Montgomery fusionado       ",
                                  "B4  radix 4 + M                  ", "B4e radix 4 + len=4 especial + M "};
        for (int modo = 0; modo < 4; ++modo) {
            const double tv = medir([&](int e) { return suma(conv_var(t, M, modo, As[e], Bs[e], n, L)); });
            std::printf("  %s %8.1f us  %.2fx\n", nombres[modo], tv / 1e3, t_lib / tv);
        }
        if (avx2) {
            const double tv = medir([&](int e) { return suma(conv_avx2(T, As[e], Bs[e], n, L)); });
            std::printf("  AVX AVX2 8x32 (DIF+DIT+Montgomery) %8.1f us  %.2fx\n", tv / 1e3, t_lib / tv);
        }
        std::printf("\n");
    }

    // ---- e) Coste del salto de potencia de 2: L = 2048 frente a L = 2049.
    {
        std::printf("== Salto de potencia de 2 (librería): L = 2048 frente a 2049, 4096 frente a 4097 ==\n");
        for (std::size_t L : {std::size_t(2048), std::size_t(2049), std::size_t(4096), std::size_t(4097)}) {
            const std::size_t la = L / 2 + 1, lb = L + 1 - la;
            std::vector<std::vector<u64>> As(kEntradas), Bs(kEntradas);
            for (int e = 0; e < kEntradas; ++e) {
                As[e].resize(la); Bs[e].resize(lb);
                for (auto& x : As[e]) x = K.random();
                for (auto& x : Bs[e]) x = K.random();
            }
            (void)fft::tables(K, std::bit_width(L - 1));
            const double tl = medir([&](int e) { return suma(fft::convolve(K, As[e], Bs[e])); });
            std::printf("  L = %5zu  %8.1f us\n", L, tl / 1e3);
        }
    }
    std::printf("\n# sumidero %llu\n", static_cast<unsigned long long>(sumidero));
}
