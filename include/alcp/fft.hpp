// alcp/fft.hpp — Generic Cooley-Tukey FFT over any finite field that owns a
// power-of-two root of unity, with PRECOMPUTED, CACHED twiddle factors.
//
// Port and generalisation of fp_fft / fp_ifft / fq_fft / fq_ifft from
// algoritmos_rapidos.py.  The original recomputed the root powers on every
// call; here the "ints necessary for the (I)FFT" — the forward and inverse
// twiddle tables and n^{-1} — are built once per size and stored in the field's
// FFTCache, then reused by every multiplication.
//
// The transform is written purely against the Field interface (add/sub/mul/
// pow/inv/...), so the *same* code runs over F_p and over F_q.
#pragma once

#include "bigint.hpp"
#include "fft_tables.hpp"
#include "fp.hpp"

#include <algorithm>
#include <bit>
#include <concepts>
#include <cstdint>
#include <cstdlib>
#include <future>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#if defined(__GNUC__) && defined(__x86_64__)
#include <immintrin.h>  // fase 7: NTT AVX2 con despacho en tiempo de ejecución
#endif

// ---------------------------------------------------------------------------
// Paralelismo (fase 3)
// ---------------------------------------------------------------------------
// Ayudantes mínimos sobre std::async (sin grupo de hilos persistente: así no
// quedan hilos vivos al descargar el módulo de Python).  Solo se usan en
// regiones sin estado compartido mutable: las operaciones de cuerpo son puras,
// las tablas FFT se construyen ANTES de entrar en paralelo y el RNG no se usa
// en la multiplicación.  El resultado no depende del número de hilos.
// Número de hilos: std::thread::hardware_concurrency(), limitable con la
// variable de entorno ALCP_THREADS (ALCP_THREADS=1 desactiva el paralelismo).
namespace alcp::par {

inline unsigned max_threads() {
    static const unsigned n = [] {
        unsigned hw = std::thread::hardware_concurrency();
        if (hw == 0) hw = 1;
        if (const char* s = std::getenv("ALCP_THREADS")) {
            const long v = std::strtol(s, nullptr, 10);
            if (v >= 1) return static_cast<unsigned>(v) < hw ? static_cast<unsigned>(v) : hw;
        }
        return hw;
    }();
    return n;
}

// f(lo, hi) sobre [0, total) troceado en bloques de al menos min_chunk
// elementos; el último bloque lo ejecuta el hilo llamante.  Las excepciones
// se propagan (future::get).
template <class F>
void parallel_for(std::size_t total, std::size_t min_chunk, F&& f) {
    std::size_t T = max_threads();
    if (min_chunk == 0) min_chunk = 1;
    if (total / min_chunk < T) T = total / min_chunk;
    if (T <= 1) { f(std::size_t(0), total); return; }
    std::vector<std::future<void>> fut;
    fut.reserve(T - 1);
    const std::size_t chunk = (total + T - 1) / T;
    std::size_t lo = 0;
    for (std::size_t t = 0; t + 1 < T && lo < total; ++t, lo += chunk) {
        const std::size_t hi = std::min(total, lo + chunk);
        fut.push_back(std::async(std::launch::async, [&f, lo, hi] { f(lo, hi); }));
    }
    if (lo < total) f(lo, total);
    for (auto& x : fut) x.get();
}

}  // namespace alcp::par

namespace alcp::fft {

// ---------------------------------------------------------------------------
// NTT en palabras de máquina (F_p con p < 2^62)
// ---------------------------------------------------------------------------
// Cambio de REPRESENTACIÓN, no de algoritmo: la transformada es la misma NTT
// radix-2 DIT, con la misma permutación, los mismos twiddles (copiados de las
// tablas BigInt) y las mismas operaciones del cuerpo (v = w·y; y = x - v;
// x = x + v; escalado final por n^{-1}).  Solo cambia cómo se hace cada
// operación: sobre uint64_t en lugar de mpz (sin reservas de memoria), con el
// producto por un twiddle fijo reducido mediante el cociente precalculado de
// Shoup.  Desde la fase 5 los valores intermedios de la transformada viven en
// [0, 4p) (mariposas perezosas, ver word::transform) y se normalizan al final;
// el resultado es idéntico bit a bit.
namespace word {

// Disponible si los limbs de GMP son de 64 bits sin «nails».
inline constexpr bool kLimbs64 = (GMP_NUMB_BITS == 64 && GMP_NAIL_BITS == 0);
// Cota para p: con p < 2^62 la reducción de Shoup y las sumas caben de sobra.
inline constexpr unsigned kMaxBits = 62;

using u64 = std::uint64_t;
using u128 = unsigned __int128;

// Valor de un BigInt en [0, 2^64) (precondición: 0 <= x < p < 2^62).
inline u64 to_word(const BigInt& x) {
    return mpz_size(x.get_mpz_t()) ? static_cast<u64>(mpz_getlimbn(x.get_mpz_t(), 0)) : 0;
}
// r = v, reutilizando la memoria de r si ya tiene un limb.
inline void set_word(BigInt& r, u64 v) {
    if (v == 0) { mpz_set_ui(r.get_mpz_t(), 0); return; }
    mp_limb_t* d = mpz_limbs_write(r.get_mpz_t(), 1);
    d[0] = static_cast<mp_limb_t>(v);
    mpz_limbs_finish(r.get_mpz_t(), 1);
}

// Mismas operaciones para cuerpos cuyo elemento ya es una palabra canónica
// en [0, p) (Fp64): la conversión es la identidad.
inline u64 to_word(u64 x) { return x; }
inline void set_word(u64& r, u64 v) { r = v; }

// Cociente de Shoup de w: floor(w·2^64 / p), con 0 <= w < p.
inline u64 shoup(u64 w, u64 p) { return static_cast<u64>((static_cast<u128>(w) << 64) / p); }
// w·x mod p con el cociente de Shoup ws de w (x < 2^64, p < 2^63): el resto
// aproximado está en [0, 2p) y una resta condicional lo deja en [0, p).
inline u64 mul_shoup(u64 w, u64 ws, u64 x, u64 p) {
    const u64 q = static_cast<u64>((static_cast<u128>(ws) * x) >> 64);
    u64 r = w * x - q * p;
    return (r >= p) ? r - p : r;
}
// a·b mod p para a, b en [0, p).
inline u64 mul_mod(u64 a, u64 b, u64 p) {
    if (p <= 0xFFFFFFFFu) return (a * b) % p;  // el producto cabe en 64 bits
    return static_cast<u64>((static_cast<u128>(a) * b) % p);
}

// ---------------------------------------------------------------------------
// Despacho AVX2 (fase 7)
// ---------------------------------------------------------------------------
// La librería se compila sin -mavx2 (el módulo de Python es «portable»): las
// rutinas vectoriales llevan __attribute__((target("avx2"))) y solo se usan si
// __builtin_cpu_supports("avx2"), que se consulta una sola vez.  Sin
// optimización (-O0: compilaciones de depuración y de cobertura) el camino se
// desactiva: GCC guarda entonces los __m256i locales en la pila con
// movimientos alineados a 32 bytes, y en MinGW la pila solo está alineada a 16
// (GCC, bug 54412).
#if defined(__GNUC__) && defined(__x86_64__) && defined(__OPTIMIZE__)
#define ALCP_NTT_AVX2 1
#else
#define ALCP_NTT_AVX2 0
#endif

inline bool avx2_disponible() {
#if ALCP_NTT_AVX2
    static const bool si = [] {
        __builtin_cpu_init();
        return __builtin_cpu_supports("avx2") != 0;
    }();
    return si;
#else
    return false;
#endif
}

// Copia de 32 bits de los mismos twiddles (W.fwd / W.inv ya rellenas) para la
// convolución AVX2; requiere p < 2^30 impar y n >= 16.
inline void build_tables32(WordNTTTables& W, u64 p) {
    using u32 = std::uint32_t;
    Word32NTTTables& T = W.w32;
    auto sh = [p](u64 x) { return static_cast<u32>((x << 32) / p); };  // floor(x·2^32 / p)
    const std::size_t m = W.fwd.size();
    T.fwd.resize(m);
    T.fwd_sh.resize(m);
    T.inv.resize(m);
    T.inv_sh.resize(m);
    for (std::size_t i = 0; i < m; ++i) {
        T.fwd[i] = static_cast<u32>(W.fwd[i]);
        T.fwd_sh[i] = sh(W.fwd[i]);
        T.inv[i] = static_cast<u32>(W.inv[i]);
        T.inv_sh[i] = sh(W.inv[i]);
    }
    // Etapa half = 4 (desde la posición 3): [w0..w3] dos veces; etapa half = 2
    // (desde la 1): [w0, w1] cuatro veces.
    for (int l = 0; l < 8; ++l) {
        T.f4[l] = T.fwd[3 + (l & 3)];
        T.f4_sh[l] = T.fwd_sh[3 + (l & 3)];
        T.f2[l] = T.fwd[1 + (l & 1)];
        T.f2_sh[l] = T.fwd_sh[1 + (l & 1)];
        T.i4[l] = T.inv[3 + (l & 3)];
        T.i4_sh[l] = T.inv_sh[3 + (l & 3)];
        T.i2[l] = T.inv[1 + (l & 1)];
        T.i2_sh[l] = T.inv_sh[1 + (l & 1)];
    }
    u32 x = static_cast<u32>(p);  // p^{-1} mod 2^32 por Newton
    for (int i = 0; i < 5; ++i) x *= 2 - static_cast<u32>(p) * x;
    T.pinv = 0u - x;
    T.c = static_cast<u32>((W.inv_n << 32) % p);  // n^{-1}·2^32 mod p
    T.c_sh = sh(T.c);
    T.p = static_cast<u32>(p);
}

// Rellena las tablas en palabras a partir de las tablas BigInt ya calculadas.
template <class Elem>
void build_tables(FFTTables<Elem>& t, u64 p) {
    const std::size_t half = t.fwd.size();
    WordNTTTables& w = t.word;
    w.fwd.reserve(2 * half);
    w.inv.reserve(2 * half);
    // Etapas len = 2..n/2 (de fwd_st/inv_st) y la última, len = n (fwd/inv).
    for (const auto& x : t.fwd_st) w.fwd.push_back(to_word(x));
    for (const auto& x : t.fwd) w.fwd.push_back(to_word(x));
    for (const auto& x : t.inv_st) w.inv.push_back(to_word(x));
    for (const auto& x : t.inv) w.inv.push_back(to_word(x));
    w.fwd_sh.resize(w.fwd.size());
    w.inv_sh.resize(w.inv.size());
    for (std::size_t i = 0; i < w.fwd.size(); ++i) w.fwd_sh[i] = shoup(w.fwd[i], p);
    for (std::size_t i = 0; i < w.inv.size(); ++i) w.inv_sh[i] = shoup(w.inv[i], p);
    w.inv_n = to_word(t.inv_n);
    w.inv_n_sh = shoup(w.inv_n, p);
    // Fase 7: constantes de Montgomery para la convolución radix 4 (solo p
    // impar; p = 2 solo admite n = 1 y sigue el camino anterior).
    if (p & 1) {
        u64 x = p;  // p·p ≡ 1 mod 8; cada paso de Newton duplica los bits
        for (int i = 0; i < 6; ++i) x *= 2 - p * x;
        w.pinv = 0 - x;
        w.mont_c = static_cast<u64>((static_cast<u128>(w.inv_n) << 64) % p);
        w.mont_c_sh = shoup(w.mont_c, p);
        // Fase 7: copia de 32 bits para la convolución AVX2 (n >= 16).
        if (p < (u64(1) << 30) && w.fwd.size() >= 15 && avx2_disponible()) build_tables32(w, p);
    }
    w.p = p;
}

// Transformada en sitio sobre a[0..n) (n = 2^k, valores en [0, p)).
template <class Elem>
void transform(const FFTTables<Elem>& t, u64* a, std::size_t n, bool inverse) {
    const WordNTTTables& W = t.word;
    const u64 p = W.p;
    const u64* tw_all = inverse ? W.inv.data() : W.fwd.data();
    const u64* ts_all = inverse ? W.inv_sh.data() : W.fwd_sh.data();

    for (const auto& [i, j] : t.swaps) std::swap(a[i], a[j]);

    // Fase 5: mariposas perezosas de Harvey ("Faster arithmetic for
    // number-theoretic transforms", 2014) y sin saltos.  Antes cada mariposa
    // hacía tres comparaciones con datos impredecibles y el predictor de saltos
    // fallaba la mitad de las veces (14 ns por mariposa; ver
    // docs/optimizacion/investigacion-fase5.md).  Ahora:
    //   * los valores viven en [0, 4p) entre etapas (4p < 2^64 porque p < 2^62);
    //   * u se lleva a [0, 2p) con una única corrección por máscara;
    //   * v = w·y con Shoup SIN la corrección final está en [0, 2p) para
    //     cualquier y < 2^64;
    //   * x' = u + v e y' = u - v + 2p quedan en [0, 4p).
    // La aritmética es la misma módulo p: solo cambia el representante
    // intermedio, y al final se normaliza a [0, p).  El resultado es idéntico
    // bit a bit al de la versión anterior.
    const u64 p2 = 2 * p;
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len >> 1;
        const u64* tw = tw_all + (half - 1);
        const u64* ts = ts_all + (half - 1);
        for (std::size_t i = 0; i < n; i += len) {
            u64* x = a + i;
            u64* y = a + i + half;
            for (std::size_t j = 0; j < half; ++j) {
                u64 u = x[j];
                u -= p2 & (0 - static_cast<u64>(u >= p2));  // u en [0, 2p)
                const u64 q = static_cast<u64>((static_cast<u128>(ts[j]) * y[j]) >> 64);
                const u64 v = tw[j] * y[j] - q * p;          // v ≡ w·y, en [0, 2p)
                x[j] = u + v;
                y[j] = u - v + p2;
            }
        }
    }

    if (inverse) {
        // mul_shoup admite cualquier entrada < 2^64 y devuelve [0, p).
        for (std::size_t i = 0; i < n; ++i) a[i] = mul_shoup(W.inv_n, W.inv_n_sh, a[i], p);
    } else {
        for (std::size_t i = 0; i < n; ++i) {
            u64 z = a[i];
            z -= p2 & (0 - static_cast<u64>(z >= p2));
            z -= p & (0 - static_cast<u64>(z >= p));
            a[i] = z;
        }
    }
}

// ---------------------------------------------------------------------------
// Convolución radix 4 sin permutación (fase 7)
// ---------------------------------------------------------------------------
// Rutinas PROPIAS de la convolución (transform, que es API pública vía
// fft::transform y devuelve el orden natural, no cambia).  En una convolución
// el orden de los valores transformados da igual mientras el producto punto a
// punto empareje los mismos índices, así que:
//   * la directa es DIF (Gentleman-Sande): entrada natural, salida en orden de
//     bits invertido, sin la permutación t.swaps;
//   * la inversa es DIT sobre esa entrada en orden de bits invertido y deja la
//     salida en orden natural, también sin permutación;
//   * las etapas se agrupan de dos en dos (radix 4: cuatro valores en registros
//     por cada pasada de memoria) y la última pareja (len = 4), cuyos twiddles
//     son 1 salvo uno, se especializa;
//   * el producto punto a punto es de Montgomery y lleva fusionado el factor
//     n^{-1} (desaparece la pasada de escalado de la inversa).
// Mismos twiddles que transform (tablas fwd/inv por etapas).  La aritmética es
// exacta módulo p y la salida se normaliza a [0, p): el resultado es idéntico
// al del camino anterior.  Medido en la fase 6 (bench/results/fase6/ntt.txt,
// variante «B4e»): 1,5-1,8x en la convolución completa.
//
// Representantes perezosos (Harvey): la directa trabaja en [0, 2p), el producto
// devuelve [0, 2p) y la inversa trabaja en [0, 4p); 4p < 2^64 porque p < 2^62.

// x en [0, 2m) -> [0, m) sin saltos.
inline u64 red(u64 x, u64 m) { return x - (m & (0 - static_cast<u64>(x >= m))); }
// w·y mod p con Shoup SIN la corrección final: en [0, 2p) para cualquier y < 2^64.
inline u64 shoup_lazy(u64 w, u64 ws, u64 y, u64 p) {
    return w * y - static_cast<u64>((static_cast<u128>(ws) * y) >> 64) * p;
}

// Directa DIF radix 4 en sitio: a[0..n) en [0, 2p) (o [0, p)) en orden
// natural -> transformada en orden de bits invertido, en [0, 2p).
inline void dif_r4(const WordNTTTables& W, u64* a, std::size_t n) {
    const u64 p = W.p, p2 = 2 * p;
    const u64* TW = W.fwd.data();
    const u64* TS = W.fwd_sh.data();
    std::size_t len = n;
    if (n >= 2 && (std::countr_zero(n) & 1)) {
        // k impar: una etapa radix 2 (len = n) y el resto, de dos en dos.
        const std::size_t half = n >> 1;
        const u64* tw = TW + (half - 1);
        const u64* ts = TS + (half - 1);
        for (std::size_t j = 0; j < half; ++j) {
            const u64 u = a[j], v = a[j + half];
            a[j] = red(u + v, p2);
            a[j + half] = shoup_lazy(tw[j], ts[j], u - v + p2, p);
        }
        len = half;
    }
    for (; len > 4; len >>= 2) {
        // Etapas len (twiddles w1, half = 2q) y len/2 (w2, half = q) fusionadas.
        const std::size_t q = len >> 2;
        const u64 *w1 = TW + (2 * q - 1), *s1 = TS + (2 * q - 1);
        const u64 *w2 = TW + (q - 1), *s2 = TS + (q - 1);
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
    if (len == 4) {
        // Última pareja de etapas: w_4^0 = w_2^0 = 1, una sola multiplicación
        // (por w_4^1) por cada cuatro valores.
        const u64 w4 = TW[2], w4s = TS[2];
        for (std::size_t i = 0; i < n; i += 4) {
            const u64 a0 = a[i], a1 = a[i + 1], a2 = a[i + 2], a3 = a[i + 3];
            const u64 b0 = red(a0 + a2, p2), b2 = red(a0 - a2 + p2, p2);
            const u64 b1 = red(a1 + a3, p2), b3 = shoup_lazy(w4, w4s, a1 - a3 + p2, p);
            a[i] = red(b0 + b1, p2);
            a[i + 1] = red(b0 - b1 + p2, p2);
            a[i + 2] = red(b2 + b3, p2);
            a[i + 3] = red(b2 - b3 + p2, p2);
        }
    }
}

// Inversa DIT radix 4 en sitio SIN el factor n^{-1} (va en el producto punto a
// punto): a[0..n) en [0, 4p) en orden de bits invertido -> orden natural, en
// [0, 4p).  Es la exacta inversa de dif_r4 salvo el factor n.
inline void dit_r4(const WordNTTTables& W, u64* a, std::size_t n) {
    const u64 p = W.p, p2 = 2 * p;
    const u64* TW = W.inv.data();
    const u64* TS = W.inv_sh.data();
    std::size_t len = 4;
    if (n >= 4) {
        const u64 w4 = TW[2], w4s = TS[2];
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
        len = 16;
    }
    for (; len <= n; len <<= 2) {
        const std::size_t q = len >> 2;
        const u64 *w1 = TW + (2 * q - 1), *s1 = TS + (2 * q - 1);  // half 2q
        const u64 *w2 = TW + (q - 1), *s2 = TS + (q - 1);          // half q
        for (std::size_t i = 0; i < n; i += len) {
            u64* x = a + i;
            for (std::size_t j = 0; j < q; ++j) {
                const u64 a0 = red(x[j], p2), a2 = red(x[j + 2 * q], p2);
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
    if (n >= 2 && (std::countr_zero(n) & 1)) {
        // k impar: falta la última etapa radix 2 (len = n).
        const std::size_t half = n >> 1;
        const u64* tw = TW + (half - 1);
        const u64* ts = TS + (half - 1);
        for (std::size_t j = 0; j < half; ++j) {
            const u64 u = red(a[j], p2);
            const u64 v = shoup_lazy(tw[j], ts[j], a[j + half], p);
            a[j] = u + v;
            a[j + half] = u - v + p2;
        }
    }
}

// a[i] <- a[i]·b[i]·n^{-1} mod p, en [0, 2p), con a[i], b[i] en [0, 2p).
// REDC(a·b) = a·b·2^{-64} está en [0, 2p) porque a·b < 4p^2 y 4p < 2^64; el
// producto por mont_c = n^{-1}·2^64 (Shoup perezoso) quita el 2^{-64} y pone
// el n^{-1}.  Requiere p impar (W.pinv != 0).
inline u64 mul_mont_n(const WordNTTTables& W, u64 a, u64 b) {
    const u128 t = static_cast<u128>(a) * b;
    const u64 m = static_cast<u64>(t) * W.pinv;
    const u64 u = static_cast<u64>((t + static_cast<u128>(m) * W.p) >> 64);
    return shoup_lazy(W.mont_c, W.mont_c_sh, u, W.p);
}

// ---------------------------------------------------------------------------
// Convolución AVX2 (fase 7): 8 carriles de 32 bits, p < 2^30
// ---------------------------------------------------------------------------
// Con p < 2^30 los representantes perezosos [0, 4p) caben en 32 bits y cada
// registro lleva 8 valores.  Mismo esquema que la convolución radix 4 (DIF sin
// permutación, producto de Montgomery con n^{-1}, DIT sin permutación), con
// mariposas radix 2 vectoriales y las cuatro etapas pequeñas (half = 8, 4, 2,
// 1) fusionadas sobre bloques de 16 valores en dos registros.  Mismos
// twiddles (Word32NTTTables).  Medido en la fase 6 (variante «AVX» de
// bench/experimentos/fase6/ntt/ntt_fase6.cpp): 7,8-9,9x en la convolución.
// Todas las cargas y almacenamientos son no alineados (loadu/storeu: nada
// supone alineación de 32 bytes) y las funciones con __m256i en su interfaz
// se expanden siempre en línea, así que ningún __m256i cruza una llamada.
#if ALCP_NTT_AVX2
namespace avx2 {

using u32 = std::uint32_t;
using V = __m256i;

#define ALCP_AVX2_FN __attribute__((target("avx2")))
#define ALCP_AVX2_IN __attribute__((target("avx2"), always_inline)) inline

ALCP_AVX2_IN V ld(const u32* q) { return _mm256_loadu_si256(reinterpret_cast<const V*>(q)); }
ALCP_AVX2_IN void st(u32* q, V v) { _mm256_storeu_si256(reinterpret_cast<V*>(q), v); }
ALCP_AVX2_IN V set1(u32 x) { return _mm256_set1_epi32(static_cast<int>(x)); }
// x en [0, 2m) -> [0, m): min(x, x - m) sin signo (si x < m, x - m da la vuelta).
ALCP_AVX2_IN V vred(V x, V m) { return _mm256_min_epu32(x, _mm256_sub_epi32(x, m)); }
// w·y con Shoup de 32 bits (ws = floor(w·2^32 / p)), en [0, 2p) para y < 2^32.
ALCP_AVX2_IN V vshoup(V w, V ws, V y, V P) {
    const V pe = _mm256_mul_epu32(ws, y);                                             // carriles pares
    const V po = _mm256_mul_epu32(_mm256_srli_epi64(ws, 32), _mm256_srli_epi64(y, 32));  // impares
    const V q = _mm256_blend_epi32(_mm256_srli_epi64(pe, 32), po, 0xAA);              // parte alta
    return _mm256_sub_epi32(_mm256_mullo_epi32(w, y), _mm256_mullo_epi32(q, P));
}
// Mariposa DIF perezosa: entradas y salidas en [0, 2p).
ALCP_AVX2_IN void vdif(V& x, V& y, V w, V ws, V P, V P2) {
    const V s = vred(_mm256_add_epi32(x, y), P2);
    const V d = _mm256_add_epi32(_mm256_sub_epi32(x, y), P2);
    x = s;
    y = vshoup(w, ws, d, P);
}
// Ídem con twiddle 1.
ALCP_AVX2_IN void vdif1(V& x, V& y, V P2) {
    const V s = vred(_mm256_add_epi32(x, y), P2);
    const V d = vred(_mm256_add_epi32(_mm256_sub_epi32(x, y), P2), P2);
    x = s;
    y = d;
}
// Mariposa DIT perezosa: entradas y salidas en [0, 4p).
ALCP_AVX2_IN void vdit(V& x, V& y, V w, V ws, V P, V P2) {
    const V u = vred(x, P2);
    const V v = vshoup(w, ws, y, P);
    x = _mm256_add_epi32(u, v);
    y = _mm256_add_epi32(_mm256_sub_epi32(u, v), P2);
}
ALCP_AVX2_IN void vdit1(V& x, V& y, V P2) {
    const V u = vred(x, P2);
    const V v = vred(y, P2);
    x = _mm256_add_epi32(u, v);
    y = _mm256_add_epi32(_mm256_sub_epi32(u, v), P2);
}
// Reparto de un bloque de 16 valores (a0 = 0..7, a1 = 8..15) en los pares de
// las etapas half = 4, 2 y 1 (X = primeros elementos, Y = segundos) y vuelta.
// rep4 y rep2 son involuciones (aplicadas dos veces dejan el bloque igual);
// para half = 1 la vuelta es jun1.
ALCP_AVX2_IN void rep4(V a0, V a1, V& X, V& Y) {
    X = _mm256_permute2x128_si256(a0, a1, 0x20);
    Y = _mm256_permute2x128_si256(a0, a1, 0x31);
}
ALCP_AVX2_IN void rep2(V a0, V a1, V& X, V& Y) {
    X = _mm256_unpacklo_epi64(a0, a1);
    Y = _mm256_unpackhi_epi64(a0, a1);
}
ALCP_AVX2_IN void rep1(V a0, V a1, V& X, V& Y) {
    const __m256 f0 = _mm256_castsi256_ps(a0), f1 = _mm256_castsi256_ps(a1);
    X = _mm256_castps_si256(_mm256_shuffle_ps(f0, f1, 0x88));
    Y = _mm256_castps_si256(_mm256_shuffle_ps(f0, f1, 0xDD));
}
ALCP_AVX2_IN void jun1(V X, V Y, V& a0, V& a1) {
    a0 = _mm256_unpacklo_epi32(X, Y);
    a1 = _mm256_unpackhi_epi32(X, Y);
}

// Directa DIF en sitio: a[0..n) en [0, p), n >= 16 -> orden de bits invertido,
// en [0, 2p).
ALCP_AVX2_FN inline void dif(const Word32NTTTables& T, u32* a, std::size_t n) {
    const V P = set1(T.p), P2 = set1(2 * T.p);
    for (std::size_t len = n; len >= 32; len >>= 1) {
        const std::size_t half = len >> 1;
        const u32* tw = T.fwd.data() + (half - 1);
        const u32* ts = T.fwd_sh.data() + (half - 1);
        for (std::size_t i = 0; i < n; i += len)
            for (std::size_t j = 0; j < half; j += 8) {
                V x = ld(a + i + j), y = ld(a + i + j + half);
                vdif(x, y, ld(tw + j), ld(ts + j), P, P2);
                st(a + i + j, x);
                st(a + i + j + half, y);
            }
    }
    const V W8 = ld(T.fwd.data() + 7), W8s = ld(T.fwd_sh.data() + 7);
    const V W4 = ld(T.f4), W4s = ld(T.f4_sh), W2 = ld(T.f2), W2s = ld(T.f2_sh);
    for (std::size_t i = 0; i < n; i += 16) {
        V a0 = ld(a + i), a1 = ld(a + i + 8), X, Y;
        vdif(a0, a1, W8, W8s, P, P2);
        rep4(a0, a1, X, Y); vdif(X, Y, W4, W4s, P, P2); rep4(X, Y, a0, a1);
        rep2(a0, a1, X, Y); vdif(X, Y, W2, W2s, P, P2); rep2(X, Y, a0, a1);
        rep1(a0, a1, X, Y); vdif1(X, Y, P2); jun1(X, Y, a0, a1);
        st(a + i, a0);
        st(a + i + 8, a1);
    }
}

// Inversa DIT en sitio sin n^{-1}: a[0..n) en [0, 4p) en orden de bits
// invertido, n >= 16 -> orden natural, en [0, 4p).
ALCP_AVX2_FN inline void dit(const Word32NTTTables& T, u32* a, std::size_t n) {
    const V P = set1(T.p), P2 = set1(2 * T.p);
    {
        const V W8 = ld(T.inv.data() + 7), W8s = ld(T.inv_sh.data() + 7);
        const V W4 = ld(T.i4), W4s = ld(T.i4_sh), W2 = ld(T.i2), W2s = ld(T.i2_sh);
        for (std::size_t i = 0; i < n; i += 16) {
            V a0 = ld(a + i), a1 = ld(a + i + 8), X, Y;
            rep1(a0, a1, X, Y); vdit1(X, Y, P2); jun1(X, Y, a0, a1);
            rep2(a0, a1, X, Y); vdit(X, Y, W2, W2s, P, P2); rep2(X, Y, a0, a1);
            rep4(a0, a1, X, Y); vdit(X, Y, W4, W4s, P, P2); rep4(X, Y, a0, a1);
            vdit(a0, a1, W8, W8s, P, P2);
            st(a + i, a0);
            st(a + i + 8, a1);
        }
    }
    for (std::size_t len = 32; len <= n; len <<= 1) {
        const std::size_t half = len >> 1;
        const u32* tw = T.inv.data() + (half - 1);
        const u32* ts = T.inv_sh.data() + (half - 1);
        for (std::size_t i = 0; i < n; i += len)
            for (std::size_t j = 0; j < half; j += 8) {
                V x = ld(a + i + j), y = ld(a + i + j + half);
                vdit(x, y, ld(tw + j), ld(ts + j), P, P2);
                st(a + i + j, x);
                st(a + i + j + half, y);
            }
    }
}

// a[i] <- a[i]·b[i]·n^{-1} en [0, 2p), con a[i], b[i] en [0, 2p) (n múltiplo
// de 8): REDC de 32 bits (a·b < 4p^2 < 2^62, resultado en [0, 2p)) y Shoup por
// c = n^{-1}·2^32.  b puede ser a (cuadrado).
ALCP_AVX2_FN inline void punto(const Word32NTTTables& T, u32* a, const u32* b, std::size_t n) {
    const V P = set1(T.p), PI = set1(T.pinv), C = set1(T.c), CS = set1(T.c_sh);
    for (std::size_t i = 0; i < n; i += 8) {
        const V x = ld(a + i), y = ld(b + i);
        const V te = _mm256_mul_epu32(x, y);  // carriles pares, 64 bits
        const V me = _mm256_mul_epu32(te, PI);
        const V ue = _mm256_srli_epi64(_mm256_add_epi64(te, _mm256_mul_epu32(me, P)), 32);
        const V to = _mm256_mul_epu32(_mm256_srli_epi64(x, 32), _mm256_srli_epi64(y, 32));  // impares
        const V mo = _mm256_mul_epu32(to, PI);
        const V uo = _mm256_add_epi64(to, _mm256_mul_epu32(mo, P));  // la parte alta ya está en su carril
        st(a + i, vshoup(C, CS, _mm256_blend_epi32(ue, uo, 0xAA), P));
    }
}

#undef ALCP_AVX2_FN
#undef ALCP_AVX2_IN

}  // namespace avx2
#endif  // ALCP_NTT_AVX2

}  // namespace word

// ---------------------------------------------------------------------------
// NTT multi-limb (F_p con p >= 2^62)
// ---------------------------------------------------------------------------
// El mismo cambio de representación para primos grandes: cada coeficiente
// ocupa N limbs contiguos de un único bloque y se opera con las funciones mpn
// de GMP sobre memoria propia (sin reservas por operación).  El producto por
// un twiddle se reduce con Montgomery: el twiddle se guarda como w·R mod p y
// REDC(y · w·R) = y·w mod p, ya en [0, p).  El producto punto a punto se
// reduce con una división mpn_tdiv_qr.  Todos los valores son canónicos: el
// resultado es idéntico bit a bit al del camino BigInt.
namespace limb {

using u64 = std::uint64_t;

// Copia el valor de x (0 <= x < p) en d[0..N).
inline void to_limbs(const BigInt& x, u64* d, std::size_t N) {
    const std::size_t s = mpz_size(x.get_mpz_t());
    for (std::size_t i = 0; i < N; ++i) d[i] = (i < s) ? mpz_getlimbn(x.get_mpz_t(), i) : 0;
}
// r = valor de d[0..N), reutilizando la memoria de r.
inline void from_limbs(BigInt& r, const u64* d, std::size_t N) {
    mp_limb_t* w = mpz_limbs_write(r.get_mpz_t(), static_cast<mp_size_t>(N));
    for (std::size_t i = 0; i < N; ++i) w[i] = d[i];
    mpz_limbs_finish(r.get_mpz_t(), static_cast<mp_size_t>(N));  // normaliza
}

// Contexto de operación: primo, constante de Montgomery y espacio de trabajo.
struct Ctx {
    std::size_t N;
    const u64* p;
    u64 pinv;
    std::vector<u64> t;   // 2N + 1 limbs
    std::vector<u64> q;   // N + 1 limbs (cociente de tdiv_qr)
    std::vector<u64> v;   // N limbs
    explicit Ctx(const LimbNTTTables& L) : N(L.N), p(L.p.data()), pinv(L.pinv), t(2 * L.N + 1), q(L.N + 1), v(L.N) {}

    const mp_limb_t* mp(const u64* x) const { return reinterpret_cast<const mp_limb_t*>(x); }
    mp_limb_t* mp(u64* x) const { return reinterpret_cast<mp_limb_t*>(x); }
    mp_size_t sz() const { return static_cast<mp_size_t>(N); }

    // out = REDC(t[0..2N)) = t·R^{-1} mod p, con t < p·R; out en [0, p).
    void redc(u64* out) {
        u64 hi = 0;
        for (std::size_t i = 0; i < N; ++i) {
            const u64 m = t[i] * pinv;
            const u64 c = mpn_addmul_1(mp(t.data() + i), mp(p), sz(), m);
            hi += mpn_add_1(mp(t.data() + i + N), mp(t.data() + i + N),
                            static_cast<mp_size_t>(N - i), c);
        }
        if (hi || mpn_cmp(mp(t.data() + N), mp(p), sz()) >= 0)
            mpn_sub_n(mp(out), mp(t.data() + N), mp(p), sz());
        else
            for (std::size_t i = 0; i < N; ++i) out[i] = t[N + i];
    }
    // out = y·w mod p, con wm = w·R mod p (twiddle en forma de Montgomery).
    void mul_mont(u64* out, const u64* y, const u64* wm) {
        mpn_mul_n(mp(t.data()), mp(y), mp(wm), sz());
        redc(out);
    }
    // out = a·b mod p (división).
    void mul_mod(u64* out, const u64* a, const u64* b) {
        mpn_mul_n(mp(t.data()), mp(a), mp(b), sz());
        mpn_tdiv_qr(mp(q.data()), mp(out), 0, mp(t.data()), 2 * sz(), mp(p), sz());
    }
};

// Rellena las tablas multi-limb a partir de las tablas BigInt.
template <class Elem>
void build_tables(FFTTables<Elem>& t, const BigInt& P) {
    LimbNTTTables& L = t.limb;
    const std::size_t N = mpz_size(P.get_mpz_t());
    L.p.resize(N);
    to_limbs(P, L.p.data(), N);
    // -p^{-1} mod 2^64 por Newton (p impar): cada paso duplica los bits.
    u64 x = L.p[0];  // p·p ≡ 1 mod 8
    for (int i = 0; i < 6; ++i) x *= 2 - L.p[0] * x;
    L.pinv = ~x + 1;
    const unsigned long rbits = 64UL * N;
    auto mont = [&](const BigInt& w, std::vector<u64>& dst) {
        BigInt m = (w << rbits) % P;  // w·R mod p
        const std::size_t off = dst.size();
        dst.resize(off + N);
        to_limbs(m, dst.data() + off, N);
    };
    for (const auto& w : t.fwd_st) mont(w, L.fwd);
    for (const auto& w : t.fwd) mont(w, L.fwd);
    for (const auto& w : t.inv_st) mont(w, L.inv);
    for (const auto& w : t.inv) mont(w, L.inv);
    mont(t.inv_n, L.inv_n);
    L.N = N;
}

// Transformada en sitio sobre a[0..n·N) (valores en [0, p)).
template <class Elem>
void transform(const FFTTables<Elem>& t, Ctx& C, u64* a, std::size_t n, bool inverse) {
    const std::size_t N = C.N;
    const u64* tw_all = inverse ? t.limb.inv.data() : t.limb.fwd.data();
    const mp_size_t sz = C.sz();

    for (const auto& [i, j] : t.swaps)
        std::swap_ranges(a + i * N, a + (i + 1) * N, a + j * N);

    u64* v = C.v.data();
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len >> 1;
        const u64* tw = tw_all + (half - 1) * N;
        for (std::size_t i = 0; i < n; i += len) {
            for (std::size_t j = 0; j < half; ++j) {
                u64* x = a + (i + j) * N;
                u64* y = a + (i + j + half) * N;
                C.mul_mont(v, y, tw + j * N);                     // v = w·y
                if (mpn_sub_n(C.mp(y), C.mp(x), C.mp(v), sz))     // y = x - v
                    mpn_add_n(C.mp(y), C.mp(y), C.mp(C.p), sz);
                if (mpn_add_n(C.mp(x), C.mp(x), C.mp(v), sz) ||   // x = x + v
                    mpn_cmp(C.mp(x), C.mp(C.p), sz) >= 0)
                    mpn_sub_n(C.mp(x), C.mp(x), C.mp(C.p), sz);
            }
        }
    }

    if (inverse)
        for (std::size_t i = 0; i < n; ++i) {
            C.mul_mont(v, a + i * N, t.limb.inv_n.data());
            for (std::size_t l = 0; l < N; ++l) a[i * N + l] = v[l];
        }
}

}  // namespace limb

// Cuerpo primo cuyo elemento ya es una palabra de máquina canónica en [0, p)
// con p < 2^62 (Fp64, fp64.hpp).  Para él la NTT en palabras no necesita
// conversiones: es el mismo camino que para Fp con p < 2^62.
template <class Field>
concept NativeWordField = std::is_same_v<typename Field::Elem, std::uint64_t> &&
                          requires(const Field& K) {
                              { K.modulus64() } -> std::same_as<std::uint64_t>;
                          };

// ¿Admite el cuerpo la NTT en palabras?  F_p con limbs de 64 bits y p < 2^62
// (Fp con p pequeño, o Fp64, que siempre cumple p < 2^62).
template <class Field>
bool word_ntt_eligible(const Field& K) {
    if constexpr (NativeWordField<Field>) {
        (void)K;
        return true;
    } else if constexpr (std::is_same_v<Field, Fp> && word::kLimbs64) {
        return mpz_sizeinbase(K.characteristic().get_mpz_t(), 2) <= word::kMaxBits;
    } else {
        (void)K;
        return false;
    }
}

// ¿Admite el cuerpo la NTT multi-limb?  F_p con limbs de 64 bits y p >= 2^62
// (impar, así que existe p^{-1} mod 2^64).
template <class Field>
bool limb_ntt_eligible(const Field& K) {
    if constexpr (std::is_same_v<Field, Fp> && word::kLimbs64) {
        return mpz_sizeinbase(K.characteristic().get_mpz_t(), 2) > word::kMaxBits &&
               mpz_odd_p(K.characteristic().get_mpz_t());
    } else {
        (void)K;
        return false;
    }
}

// Implementación de la NTT que usará el cuerpo (para los umbrales de despacho
// de Poly::mul, que dependen mucho de ella).
enum class Backend { Word, Limb, Generic };
template <class Field>
Backend backend(const Field& K) {
    if (word_ntt_eligible(K)) return Backend::Word;
    if (limb_ntt_eligible(K)) return Backend::Limb;
    return Backend::Generic;
}

// Smallest k with 2^k >= n   (ceil log2).  ceil_log2(0)=ceil_log2(1)=0.
inline unsigned long ceil_log2(std::size_t n) {
    unsigned long k = 0;
    std::size_t pw = 1;
    while (pw < n) { pw <<= 1; ++k; }
    return k;
}

// Build (or fetch) the precomputed twiddle tables for a transform of size
// n = 2^k in field K.  Requires k <= K.fft_capacity().  The result lives in
// K.fft_cache_ and is reused on every later transform of the same size.
// Precondición comprobada (std::invalid_argument): k <= K.fft_capacity() y
// 2^k representable en std::size_t.  Solo se comprueba al construir una tabla
// nueva: las que ya están en la caché se validaron al insertarlas.
namespace detail {

// Comprueba la precondición de tables()/tablas_sin_rng() (ver arriba).
template <class Field>
void valida_k(const Field& K, unsigned long k) {
    if (k >= static_cast<unsigned long>(std::numeric_limits<std::size_t>::digits) - 1)
        throw std::invalid_argument("fft::tables/tablas_sin_rng: k = " + std::to_string(k) + " demasiado grande");
    if (k > K.fft_capacity())
        throw std::invalid_argument("fft::tables/tablas_sin_rng: el cuerpo no tiene raíces de la unidad de orden 2^" +
                                    std::to_string(k) + " (capacidad FFT " +
                                    std::to_string(K.fft_capacity()) + ")");
}

// Construye las tablas de tamaño n = 2^k a partir de una raíz base de orden
// 2^capacidad (base).  Es el cuerpo de tables() de siempre, con la raíz como
// parámetro para poder usar también la raíz determinista (tablas_sin_rng).
template <class Field>
std::shared_ptr<const FFTTables<typename Field::Elem>> construye_tablas(const Field& K, unsigned long k,
                                                                         const typename Field::Elem& base) {
    using Elem = typename Field::Elem;
    auto t = std::make_shared<FFTTables<Elem>>();
    t->k = k;
    const std::size_t n = std::size_t(1) << k;
    const std::size_t half = n / 2;

    // Primitive n-th root of unity: base_root has order 2^capacity, so raising
    // it to 2^(capacity - k) yields an element of order exactly n.
    unsigned long cap = K.fft_capacity();
    Elem g = K.pow(base, BigInt(1) << (cap - k));
    Elem g_inv = K.inv(g);

    t->fwd.resize(half);
    t->inv.resize(half);
    if (half > 0) {
        t->fwd[0] = K.one();
        t->inv[0] = K.one();
        for (std::size_t i = 1; i < half; ++i) {
            t->fwd[i] = K.mul(t->fwd[i - 1], g);
            t->inv[i] = K.mul(t->inv[i - 1], g_inv);
        }
    }
    t->inv_n = K.inv(K.from_int(BigInt(static_cast<unsigned long>(n))));

    // Twiddles por etapa, contiguos (ver fft_tables.hpp): copia de los mismos
    // valores reordenados; no se calcula ninguna potencia nueva.
    if (half > 1) {
        t->fwd_st.reserve(half - 1);
        t->inv_st.reserve(half - 1);
        for (std::size_t len = 2; len < n; len <<= 1) {
            const std::size_t step = n / len;
            for (std::size_t j = 0; j < len / 2; ++j) {
                t->fwd_st.push_back(t->fwd[j * step]);
                t->inv_st.push_back(t->inv[j * step]);
            }
        }
    }

    // Permutación de inversión de bits, precalculada como lista de intercambios.
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) t->swaps.emplace_back(i, j);
    }

    // Copia en palabras de máquina de las mismas tablas (solo F_p, p < 2^62).
    // Copia multi-limb para el resto de F_p (p >= 2^62, impar).
    if constexpr (std::is_same_v<Field, Fp>) {
        if (word_ntt_eligible(K))
            word::build_tables(*t, word::to_word(K.characteristic()));
        else if (limb_ntt_eligible(K))
            limb::build_tables(*t, K.characteristic());
    } else if constexpr (NativeWordField<Field>) {
        word::build_tables(*t, K.modulus64());
    }
    return t;
}

}  // namespace detail

template <class Field>
const FFTTables<typename Field::Elem>& tables(const Field& K, unsigned long k) {
    // Acceso directo por índice: la caché es un vector indexado por k.
    if (k < K.fft_cache_.size() && K.fft_cache_[k]) return *K.fft_cache_[k];
    detail::valida_k(K, k);
    // fft_base_root() saca la raíz del generador aleatorio la primera vez.
    auto t = detail::construye_tablas(K, k, K.fft_base_root());
    if (K.fft_cache_.size() <= k) K.fft_cache_.resize(k + 1);
    K.fft_cache_[k] = t;
    return *t;
}

// ---------------------------------------------------------------------------
// Tablas sin el generador aleatorio (fase 7, propuesta 3d; solo Fp64)
// ---------------------------------------------------------------------------
// Las rutas que no deben consumir el generador (mul_no_rng, Newton, Reducer,
// Composer, half-GCD) solo podían usar la NTT si la caché ya tenía tablas, es
// decir, si la raíz base ya se había sacado del generador.  Desde Python, cada
// cuerpo nace con la caché vacía y esas rutas iban por Kronecker o Karatsuba.
//
// El producto calculado con la NTT no depende de la raíz primitiva elegida
// (es la convolución exacta módulo p), así que esas rutas pueden usar tablas
// construidas con una raíz DETERMINISTA, guardadas en una caché aparte
// (fft_cache_det_).  fft_base_root(), tables() y fft_cache_ no cambian: el
// generador se consume exactamente igual que antes.
//
// Raíz determinista: el primer a = 2, 3, ... tal que c = a^((p-1)/2^cap)
// cumple c^(2^(cap-1)) != 1.  Termina: c tiene orden 2^cap exactamente si y
// solo si a no es un resto cuadrático (a^((p-1)/2) = c^(2^(cap-1)) = -1), y
// hay no restos en [2, p-1] porque p es un primo impar.
//
// Cuidado: una transformada solo se puede combinar con otra hecha con las
// MISMAS tablas.  tablas_sin_rng() devuelve las de fft_cache_ si ya existen y
// las deterministas si no, y la caché puede llenarse entre dos llamadas (p. ej.
// un Reducer que vive durante un factor()).  Por eso medio::Fija guarda el
// puntero a las tablas con las que se transformó, y las rutas que no guardan
// transformadas piden las tablas una sola vez por operación.
template <class Field>
    requires NativeWordField<Field>
typename Field::Elem raiz_det(const Field& K) {
    const unsigned long cap = K.fft_capacity();
    if (cap == 0) throw std::invalid_argument("fft::raiz_det: el cuerpo no admite NTT (capacidad FFT 0)");
    const BigInt pm1 = K.characteristic() - 1;
    const BigInt cofactor = pm1 >> cap;
    const BigInt half = BigInt(1) << (cap - 1);  // cap >= 1 (comprobado arriba)
    for (unsigned long a = 2;; ++a) {
        const auto c = K.pow(K.from_int(BigInt(a)), cofactor);
        if (K.pow(c, half) != K.one()) return c;
    }
}

template <class Field>
    requires NativeWordField<Field>
std::shared_ptr<const FFTTables<typename Field::Elem>> tablas_sin_rng(const Field& K, unsigned long k) {
    if (!K.fft_cache_.empty()) {
        (void)tables(K, k);  // la caché no está vacía: la raíz ya existe, no consume el generador
        return K.fft_cache_[k];
    }
    if (k < K.fft_cache_det_.size() && K.fft_cache_det_[k]) return K.fft_cache_det_[k];
    detail::valida_k(K, k);
    auto t = detail::construye_tablas(K, k, raiz_det(K));
    if (K.fft_cache_det_.size() <= k) K.fft_cache_det_.resize(k + 1);
    K.fft_cache_det_[k] = t;
    return t;
}

namespace detail {

// Mariposas DIT de las etapas len = 2..n sobre `a` ya permutado.  Para la
// etapa len < n los twiddles se leen seguidos de `st` (desde len/2 - 1); para
// len = n, de `last` (la tabla completa, paso 1).  Las operaciones son las del
// algoritmo original: v = w·a[i+j+h]; a[i+j+h] = a[i+j] - v; a[i+j] = a[i+j] + v.
// Se calcula la resta antes que la suma para no tener que copiar a[i+j].
//
// Paralelismo: las n/2 mariposas de una etapa son independientes entre sí
// (cada una toca un par distinto); se reparten por bloques entre hilos y se
// espera a todos antes de la etapa siguiente.  Bloque mínimo: kParChunk
// mariposas (más grande si el elemento es un entero barato).
// Con elementos de F_q sobre Fp64 (vector de palabras) cada operación es
// barata y está dominada por la reserva de memoria del vector resultado, que
// escala mal con muchos hilos: se usa el mismo bloque que con BigInt.
template <class Elem>
inline constexpr std::size_t kParChunk =
    (std::is_same_v<Elem, BigInt> || std::is_same_v<Elem, std::vector<std::uint64_t>>) ? 1024 : 32;

template <class Field>
void butterflies(const Field& K, typename Field::Elem* a, std::size_t n,
                 const std::vector<typename Field::Elem>& st,
                 const std::vector<typename Field::Elem>& last) {
    using Elem = typename Field::Elem;
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len >> 1;
        const Elem* tw = (len == n) ? last.data() : st.data() + (half - 1);
        // Mariposas b en [lo, hi) de esta etapa: bloque i = (b / half)·len,
        // posición j = b % half.
        auto body = [&](std::size_t lo, std::size_t hi) {
            for (std::size_t b = lo; b < hi; ++b) {
                const std::size_t j = b % half;
                Elem* x = a + (b / half) * len + j;
                Elem* y = x + half;
                Elem v = K.mul(tw[j], *y);
                *y = K.sub(*x, v);
                *x = K.add(*x, v);
            }
        };
        par::parallel_for(n / 2, kParChunk<Elem>, body);
    }
}

}  // namespace detail

// In-place iterative radix-2 DIT transform of `a` (size must be 2^k).
// inverse=false: forward DFT.  inverse=true: inverse DFT (includes the n^{-1}
// scaling), exactly IDFT_{n,g} = n^{-1} * DFT_{n,g^{-1}}.
// Precondiciones comprobadas (std::invalid_argument): a.size() potencia de 2 y
// k <= K.fft_capacity() (en tables).  Se comprueban ANTES de tocar `a`, así que
// si se lanza la excepción el vector queda intacto.
template <class Field>
void transform(const Field& K, std::vector<typename Field::Elem>& a, bool inverse) {
    const std::size_t n = a.size();
    if (n <= 1) return;
    if ((n & (n - 1)) != 0)
        throw std::invalid_argument("fft::transform: la longitud " + std::to_string(n) +
                                    " no es potencia de 2");
    const unsigned long k = ceil_log2(n);

    // Las tablas no dependen de `a`: obtenerlas antes de la permutación no
    // cambia el resultado (ni el consumo del generador aleatorio).
    const auto& t = tables(K, k);

    if constexpr (NativeWordField<Field>) {
        // Fp64: los elementos ya son palabras canónicas; misma transformada
        // en palabras de máquina, en sitio y sin conversiones.
        word::transform(t, a.data(), n, inverse);
    } else {
        // F_p con GMP (palabras o multi-limb) y camino genérico.  Va en el
        // else para que no se instancie para Fp64.
        if constexpr (std::is_same_v<Field, Fp>) {
            if (t.word.p != 0) {  // misma transformada, en palabras de máquina
                std::vector<word::u64> w(n);
                for (std::size_t i = 0; i < n; ++i) w[i] = word::to_word(a[i]);
                word::transform(t, w.data(), n, inverse);
                for (std::size_t i = 0; i < n; ++i) word::set_word(a[i], w[i]);
                return;
            }
            if (t.limb.N != 0) {  // misma transformada, multi-limb
                const std::size_t N = t.limb.N;
                limb::Ctx C(t.limb);
                std::vector<limb::u64> w(n * N);
                for (std::size_t i = 0; i < n; ++i) limb::to_limbs(a[i], w.data() + i * N, N);
                limb::transform(t, C, w.data(), n, inverse);
                for (std::size_t i = 0; i < n; ++i) limb::from_limbs(a[i], w.data() + i * N, N);
                return;
            }
        }

        // Bit-reversal permutation (precalculada en las tablas).
        for (const auto& [i, j] : t.swaps) std::swap(a[i], a[j]);

        if (inverse) detail::butterflies(K, a.data(), n, t.inv_st, t.inv);
        else         detail::butterflies(K, a.data(), n, t.fwd_st, t.fwd);

        if (inverse)
            par::parallel_for(n, 2 * detail::kParChunk<typename Field::Elem>,
                              [&](std::size_t lo, std::size_t hi) {
                                  for (std::size_t i = lo; i < hi; ++i) a[i] = K.mul(a[i], t.inv_n);
                              });
    }
}

// Is an FFT-based convolution of total length `result_len` possible in K?
// Returns the transform exponent k (>=0) if 2^k >= result_len and the field
// supports a transform of that size; otherwise returns -1.
template <class Field>
long fft_exponent(const Field& K, std::size_t result_len) {
    if (result_len == 0) return -1;
    unsigned long k = ceil_log2(result_len);
    return (k <= K.fft_capacity()) ? static_cast<long>(k) : -1;
}

// Tamaño de transformada a partir del cual las dos transformadas directas de
// convolve van en hilos distintos (camino en palabras / multi-limb).
inline constexpr std::size_t kParWordMin = std::size_t(1) << 13;
inline constexpr std::size_t kParLimbMin = std::size_t(1) << 9;

namespace detail {

// Convolución por el camino en palabras de máquina (tablas t con t.word.p != 0;
// A y B no vacíos; n = 2^k >= L = |A| + |B| - 1).  Mismo esquema que el camino
// genérico: transformar A y B, producto punto a punto, transformada inversa.
// Se convierte una vez a la entrada y una vez a la salida (identidad si el
// elemento ya es una palabra, Fp64).
template <class Elem>
std::vector<Elem> convolve_word(const FFTTables<Elem>& t, const std::vector<Elem>& A,
                                const std::vector<Elem>& B, std::size_t n, std::size_t L) {
    using word::u64;
    const u64 p = t.word.p;
#if ALCP_NTT_AVX2
    if (t.word.w32.p != 0) {
        // Fase 7: el mismo esquema en 8 carriles de 32 bits (p < 2^30, n >= 16
        // y CPU con AVX2: si no, w32.p == 0).  Ver «Convolución AVX2».
        using u32 = std::uint32_t;
        const Word32NTTTables& T = t.word.w32;
        std::vector<u32> wa(n, 0);
        for (std::size_t i = 0; i < A.size(); ++i) wa[i] = static_cast<u32>(word::to_word(A[i]));
        if (&A == &B) {  // cuadrado: una sola transformada directa
            word::avx2::dif(T, wa.data(), n);
            word::avx2::punto(T, wa.data(), wa.data(), n);
        } else {
            std::vector<u32> wb(n, 0);
            for (std::size_t i = 0; i < B.size(); ++i) wb[i] = static_cast<u32>(word::to_word(B[i]));
            if (n >= kParWordMin && par::max_threads() > 1) {
                auto fb = std::async(std::launch::async, [&] { word::avx2::dif(T, wb.data(), n); });
                word::avx2::dif(T, wa.data(), n);
                fb.get();
            } else {
                word::avx2::dif(T, wa.data(), n);
                word::avx2::dif(T, wb.data(), n);
            }
            word::avx2::punto(T, wa.data(), wb.data(), n);
        }
        word::avx2::dit(T, wa.data(), n);
        // [0, 4p) -> [0, p) al copiar los L coeficientes del resultado.
        std::vector<Elem> out(L);
        for (std::size_t i = 0; i < L; ++i) word::set_word(out[i], word::red(word::red(wa[i], 2 * p), p));
        return out;
    }
#endif
    std::vector<u64> wa(n, 0);
    for (std::size_t i = 0; i < A.size(); ++i) wa[i] = word::to_word(A[i]);
    if (t.word.pinv != 0) {
        // Fase 7: directa DIF radix 4, producto de Montgomery con n^{-1} e
        // inversa DIT radix 4, sin permutaciones (ver «Convolución radix 4»).
        const WordNTTTables& W = t.word;
        if (&A == &B) {  // cuadrado: una sola transformada directa
            word::dif_r4(W, wa.data(), n);
            for (std::size_t i = 0; i < n; ++i) wa[i] = word::mul_mont_n(W, wa[i], wa[i]);
        } else {
            std::vector<u64> wb(n, 0);
            for (std::size_t i = 0; i < B.size(); ++i) wb[i] = word::to_word(B[i]);
            if (n >= kParWordMin && par::max_threads() > 1) {
                auto fb = std::async(std::launch::async, [&] { word::dif_r4(W, wb.data(), n); });
                word::dif_r4(W, wa.data(), n);
                fb.get();
            } else {
                word::dif_r4(W, wa.data(), n);
                word::dif_r4(W, wb.data(), n);
            }
            for (std::size_t i = 0; i < n; ++i) wa[i] = word::mul_mont_n(W, wa[i], wb[i]);
        }
        word::dit_r4(W, wa.data(), n);
        // Solo hacen falta los L primeros coeficientes: se normalizan al copiar.
        std::vector<Elem> out(L);
        for (std::size_t i = 0; i < L; ++i) word::set_word(out[i], word::red(word::red(wa[i], 2 * p), p));
        return out;
    }
    if (&A == &B) {  // cuadrado: una sola transformada directa
        word::transform(t, wa.data(), n, false);
        for (std::size_t i = 0; i < n; ++i) wa[i] = word::mul_mod(wa[i], wa[i], p);
    } else {
        std::vector<u64> wb(n, 0);
        for (std::size_t i = 0; i < B.size(); ++i) wb[i] = word::to_word(B[i]);
        // Las dos transformadas directas son independientes: con n
        // grande, la de B va en otro hilo (las tablas ya existen).
        if (n >= kParWordMin && par::max_threads() > 1) {
            auto fb = std::async(std::launch::async,
                                 [&] { word::transform(t, wb.data(), n, false); });
            word::transform(t, wa.data(), n, false);
            fb.get();
        } else {
            word::transform(t, wa.data(), n, false);
            word::transform(t, wb.data(), n, false);
        }
        for (std::size_t i = 0; i < n; ++i) wa[i] = word::mul_mod(wa[i], wb[i], p);
    }
    word::transform(t, wa.data(), n, true);
    std::vector<Elem> out(L);
    for (std::size_t i = 0; i < L; ++i) word::set_word(out[i], wa[i]);
    return out;
}

}  // namespace detail

// Linear convolution of coefficient vectors A and B via FFT.  Caller must have
// checked fft_exponent(K, A.size()+B.size()-1) >= 0.
// Si A y B son el MISMO vector (mismo objeto), el producto es un cuadrado y
// se transforma una sola vez.
template <class Field>
std::vector<typename Field::Elem> convolve(const Field& K,
                                           const std::vector<typename Field::Elem>& A,
                                           const std::vector<typename Field::Elem>& B) {
    using Elem = typename Field::Elem;
    if (A.empty() || B.empty()) return {};
    const std::size_t L = A.size() + B.size() - 1;
    const unsigned long k = ceil_log2(L);
    const std::size_t n = std::size_t(1) << k;

    if constexpr (NativeWordField<Field>) {
        // Fp64: siempre por el camino en palabras (p < 2^62).  tables() se
        // llama en el mismo punto que en el camino genérico (mismo consumo
        // del RNG).
        return detail::convolve_word(tables(K, k), A, B, n, L);
    } else {
        // F_p con GMP (palabras o multi-limb) y camino genérico.  Va en el
        // else para que no se instancie para Fp64.
        if constexpr (std::is_same_v<Field, Fp>) {
            // Mismo esquema (transformar A y B, producto punto a punto,
            // transformada inversa) en palabras de máquina: se convierte una vez
            // a la entrada y una vez a la salida.  tables() se llama en el mismo
            // punto que en el camino genérico (mismo consumo del RNG).
            const auto& t = tables(K, k);
            if (t.word.p != 0) return detail::convolve_word(t, A, B, n, L);
            if (t.limb.N != 0) {  // ídem, multi-limb
                using limb::u64;
                const std::size_t N = t.limb.N;
                limb::Ctx C(t.limb);
                std::vector<u64> wa(n * N, 0);
                for (std::size_t i = 0; i < A.size(); ++i) limb::to_limbs(A[i], wa.data() + i * N, N);
                if (&A == &B) {  // cuadrado: una sola transformada directa
                    limb::transform(t, C, wa.data(), n, false);
                    for (std::size_t i = 0; i < n; ++i)
                        C.mul_mod(wa.data() + i * N, wa.data() + i * N, wa.data() + i * N);
                } else {
                    std::vector<u64> wb(n * N, 0);
                    for (std::size_t i = 0; i < B.size(); ++i) limb::to_limbs(B[i], wb.data() + i * N, N);
                    // Transformadas de A y B en paralelo (cada hilo con su propio
                    // espacio de trabajo Ctx).
                    if (n >= kParLimbMin && par::max_threads() > 1) {
                        auto fb = std::async(std::launch::async, [&] {
                            limb::Ctx C2(t.limb);
                            limb::transform(t, C2, wb.data(), n, false);
                        });
                        limb::transform(t, C, wa.data(), n, false);
                        fb.get();
                    } else {
                        limb::transform(t, C, wa.data(), n, false);
                        limb::transform(t, C, wb.data(), n, false);
                    }
                    for (std::size_t i = 0; i < n; ++i)
                        C.mul_mod(wa.data() + i * N, wa.data() + i * N, wb.data() + i * N);
                }
                limb::transform(t, C, wa.data(), n, true);
                std::vector<Elem> out(L);
                for (std::size_t i = 0; i < L; ++i) limb::from_limbs(out[i], wa.data() + i * N, N);
                return out;
            }
        }

        // Copia directa de los coeficientes y relleno con ceros (sin construir
        // primero n ceros para sobrescribirlos después).
        std::vector<Elem> fa;
        fa.reserve(n);
        fa.assign(A.begin(), A.end());
        fa.resize(n, K.zero());
        transform(K, fa, false);

        const std::size_t chunk = 2 * detail::kParChunk<Elem>;
        if (&A == &B) {
            // Cuadrado (a*a, frecuente en pow_mod y Ben-Or): la transformada
            // directa de B sería idéntica a la de A, así que se reutiliza.  Es
            // evitar trabajo redundante idéntico: mismo algoritmo, mismo resultado.
            par::parallel_for(n, chunk, [&](std::size_t lo, std::size_t hi) {
                for (std::size_t i = lo; i < hi; ++i) fa[i] = K.mul(fa[i], fa[i]);
            });
        } else {
            std::vector<Elem> fb;
            fb.reserve(n);
            fb.assign(B.begin(), B.end());
            fb.resize(n, K.zero());
            transform(K, fb, false);
            par::parallel_for(n, chunk, [&](std::size_t lo, std::size_t hi) {
                for (std::size_t i = lo; i < hi; ++i) fa[i] = K.mul(fa[i], fb[i]);
            });
        }
        transform(K, fa, true);

        fa.resize(L);
        return fa;
    }
}

}  // namespace alcp::fft
