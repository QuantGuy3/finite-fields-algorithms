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

#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

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
// Shoup.  Todas las operaciones devuelven el representante canónico en
// [0, p), así que el resultado es idéntico bit a bit.
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

    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len >> 1;
        const u64* tw = tw_all + (half - 1);
        const u64* ts = ts_all + (half - 1);
        for (std::size_t i = 0; i < n; i += len) {
            u64* x = a + i;
            u64* y = a + i + half;
            for (std::size_t j = 0; j < half; ++j) {
                const u64 v = mul_shoup(tw[j], ts[j], y[j], p);
                const u64 u = x[j];
                const u64 s = u + v;
                x[j] = (s >= p) ? s - p : s;
                y[j] = (u >= v) ? u - v : u + p - v;
            }
        }
    }

    if (inverse)
        for (std::size_t i = 0; i < n; ++i) a[i] = mul_shoup(W.inv_n, W.inv_n_sh, a[i], p);
}

}  // namespace word

// ¿Admite el cuerpo la NTT en palabras?  F_p con limbs de 64 bits y p < 2^62.
template <class Field>
bool word_ntt_eligible(const Field& K) {
    if constexpr (std::is_same_v<Field, Fp> && word::kLimbs64) {
        return mpz_sizeinbase(K.characteristic().get_mpz_t(), 2) <= word::kMaxBits;
    } else {
        (void)K;
        return false;
    }
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
template <class Field>
const FFTTables<typename Field::Elem>& tables(const Field& K, unsigned long k) {
    using Elem = typename Field::Elem;
    // Acceso directo por índice: la caché es un vector indexado por k.
    if (k < K.fft_cache_.size() && K.fft_cache_[k]) return *K.fft_cache_[k];

    if (k >= static_cast<unsigned long>(std::numeric_limits<std::size_t>::digits) - 1)
        throw std::invalid_argument("fft::tables: k = " + std::to_string(k) + " demasiado grande");
    if (k > K.fft_capacity())
        throw std::invalid_argument("fft::tables: el cuerpo no tiene raíces de la unidad de orden 2^" +
                                    std::to_string(k) + " (capacidad FFT " +
                                    std::to_string(K.fft_capacity()) + ")");

    auto t = std::make_shared<FFTTables<Elem>>();
    t->k = k;
    const std::size_t n = std::size_t(1) << k;
    const std::size_t half = n / 2;

    // Primitive n-th root of unity: base_root has order 2^capacity, so raising
    // it to 2^(capacity - k) yields an element of order exactly n.
    unsigned long cap = K.fft_capacity();
    Elem g = K.pow(K.fft_base_root(), BigInt(1) << (cap - k));
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
    if constexpr (std::is_same_v<Field, Fp>)
        if (word_ntt_eligible(K)) word::build_tables(*t, word::to_word(K.characteristic()));

    if (K.fft_cache_.size() <= k) K.fft_cache_.resize(k + 1);
    K.fft_cache_[k] = t;
    return *t;
}

namespace detail {

// Mariposas DIT de las etapas len = 2..n sobre `a` ya permutado.  Para la
// etapa len < n los twiddles se leen seguidos de `st` (desde len/2 - 1); para
// len = n, de `last` (la tabla completa, paso 1).  Las operaciones son las del
// algoritmo original: v = w·a[i+j+h]; a[i+j+h] = a[i+j] - v; a[i+j] = a[i+j] + v.
// Se calcula la resta antes que la suma para no tener que copiar a[i+j].
template <class Field>
void butterflies(const Field& K, typename Field::Elem* a, std::size_t n,
                 const std::vector<typename Field::Elem>& st,
                 const std::vector<typename Field::Elem>& last) {
    using Elem = typename Field::Elem;
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len >> 1;
        const Elem* tw = (len == n) ? last.data() : st.data() + (half - 1);
        for (std::size_t i = 0; i < n; i += len) {
            Elem* x = a + i;
            Elem* y = a + i + half;
            for (std::size_t j = 0; j < half; ++j) {
                Elem v = K.mul(tw[j], y[j]);
                y[j] = K.sub(x[j], v);
                x[j] = K.add(x[j], v);
            }
        }
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

    if constexpr (std::is_same_v<Field, Fp>) {
        if (t.word.p != 0) {  // misma transformada, en palabras de máquina
            std::vector<word::u64> w(n);
            for (std::size_t i = 0; i < n; ++i) w[i] = word::to_word(a[i]);
            word::transform(t, w.data(), n, inverse);
            for (std::size_t i = 0; i < n; ++i) word::set_word(a[i], w[i]);
            return;
        }
    }

    // Bit-reversal permutation (precalculada en las tablas).
    for (const auto& [i, j] : t.swaps) std::swap(a[i], a[j]);

    if (inverse) detail::butterflies(K, a.data(), n, t.inv_st, t.inv);
    else         detail::butterflies(K, a.data(), n, t.fwd_st, t.fwd);

    if (inverse)
        for (auto& x : a) x = K.mul(x, t.inv_n);
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

    if constexpr (std::is_same_v<Field, Fp>) {
        // Mismo esquema (transformar A y B, producto punto a punto,
        // transformada inversa) en palabras de máquina: se convierte una vez
        // a la entrada y una vez a la salida.  tables() se llama en el mismo
        // punto que en el camino genérico (mismo consumo del RNG).
        const auto& t = tables(K, k);
        if (t.word.p != 0) {
            using word::u64;
            const u64 p = t.word.p;
            std::vector<u64> wa(n, 0);
            for (std::size_t i = 0; i < A.size(); ++i) wa[i] = word::to_word(A[i]);
            word::transform(t, wa.data(), n, false);
            if (&A == &B) {  // cuadrado: una sola transformada directa
                for (std::size_t i = 0; i < n; ++i) wa[i] = word::mul_mod(wa[i], wa[i], p);
            } else {
                std::vector<u64> wb(n, 0);
                for (std::size_t i = 0; i < B.size(); ++i) wb[i] = word::to_word(B[i]);
                word::transform(t, wb.data(), n, false);
                for (std::size_t i = 0; i < n; ++i) wa[i] = word::mul_mod(wa[i], wb[i], p);
            }
            word::transform(t, wa.data(), n, true);
            std::vector<Elem> out(L);
            for (std::size_t i = 0; i < L; ++i) word::set_word(out[i], wa[i]);
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

    if (&A == &B) {
        // Cuadrado (a*a, frecuente en pow_mod y Ben-Or): la transformada
        // directa de B sería idéntica a la de A, así que se reutiliza.  Es
        // evitar trabajo redundante idéntico: mismo algoritmo, mismo resultado.
        for (std::size_t i = 0; i < n; ++i) fa[i] = K.mul(fa[i], fa[i]);
    } else {
        std::vector<Elem> fb;
        fb.reserve(n);
        fb.assign(B.begin(), B.end());
        fb.resize(n, K.zero());
        transform(K, fb, false);
        for (std::size_t i = 0; i < n; ++i) fa[i] = K.mul(fa[i], fb[i]);
    }
    transform(K, fa, true);

    fa.resize(L);
    return fa;
}

}  // namespace alcp::fft
