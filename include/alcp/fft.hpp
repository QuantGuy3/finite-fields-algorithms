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
    // Copia multi-limb para el resto de F_p (p >= 2^62, impar).
    if constexpr (std::is_same_v<Field, Fp>) {
        if (word_ntt_eligible(K))
            word::build_tables(*t, word::to_word(K.characteristic()));
        else if (limb_ntt_eligible(K))
            limb::build_tables(*t, K.characteristic());
    }

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
        if (t.limb.N != 0) {  // ídem, multi-limb
            using limb::u64;
            const std::size_t N = t.limb.N;
            limb::Ctx C(t.limb);
            std::vector<u64> wa(n * N, 0);
            for (std::size_t i = 0; i < A.size(); ++i) limb::to_limbs(A[i], wa.data() + i * N, N);
            limb::transform(t, C, wa.data(), n, false);
            if (&A == &B) {  // cuadrado: una sola transformada directa
                for (std::size_t i = 0; i < n; ++i)
                    C.mul_mod(wa.data() + i * N, wa.data() + i * N, wa.data() + i * N);
            } else {
                std::vector<u64> wb(n * N, 0);
                for (std::size_t i = 0; i < B.size(); ++i) limb::to_limbs(B[i], wb.data() + i * N, N);
                limb::transform(t, C, wb.data(), n, false);
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
