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

#include <vector>

namespace alcp::fft {

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
template <class Field>
const FFTTables<typename Field::Elem>& tables(const Field& K, unsigned long k) {
    using Elem = typename Field::Elem;
    auto it = K.fft_cache_.find(k);
    if (it != K.fft_cache_.end()) return it->second;

    FFTTables<Elem> t;
    t.k = k;
    const std::size_t n = std::size_t(1) << k;
    const std::size_t half = n / 2;

    // Primitive n-th root of unity: base_root has order 2^capacity, so raising
    // it to 2^(capacity - k) yields an element of order exactly n.
    unsigned long cap = K.fft_capacity();
    Elem g = K.pow(K.fft_base_root(), BigInt(1) << (cap - k));
    Elem g_inv = K.inv(g);

    t.fwd.resize(half);
    t.inv.resize(half);
    if (half > 0) {
        t.fwd[0] = K.one();
        t.inv[0] = K.one();
        for (std::size_t i = 1; i < half; ++i) {
            t.fwd[i] = K.mul(t.fwd[i - 1], g);
            t.inv[i] = K.mul(t.inv[i - 1], g_inv);
        }
    }
    t.inv_n = K.inv(K.from_int(BigInt(static_cast<unsigned long>(n))));

    auto [pos, ok] = K.fft_cache_.emplace(k, std::move(t));
    (void)ok;
    return pos->second;
}

// In-place iterative radix-2 DIT transform of `a` (size must be 2^k).
// inverse=false: forward DFT.  inverse=true: inverse DFT (includes the n^{-1}
// scaling), exactly IDFT_{n,g} = n^{-1} * DFT_{n,g^{-1}}.
template <class Field>
void transform(const Field& K, std::vector<typename Field::Elem>& a, bool inverse) {
    const std::size_t n = a.size();
    if (n <= 1) return;
    const unsigned long k = ceil_log2(n);

    // Bit-reversal permutation.
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    const auto& t = tables(K, k);
    const auto& tw = inverse ? t.inv : t.fwd;  // length n/2, full-size twiddles

    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t step = n / len;       // stride into the twiddle table
        const std::size_t half = len >> 1;
        for (std::size_t i = 0; i < n; i += len) {
            for (std::size_t j = 0; j < half; ++j) {
                const auto& w = tw[j * step];
                typename Field::Elem u = a[i + j];
                typename Field::Elem v = K.mul(w, a[i + j + half]);
                a[i + j] = K.add(u, v);
                a[i + j + half] = K.sub(u, v);
            }
        }
    }

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
template <class Field>
std::vector<typename Field::Elem> convolve(const Field& K,
                                           const std::vector<typename Field::Elem>& A,
                                           const std::vector<typename Field::Elem>& B) {
    using Elem = typename Field::Elem;
    if (A.empty() || B.empty()) return {};
    const std::size_t L = A.size() + B.size() - 1;
    const unsigned long k = ceil_log2(L);
    const std::size_t n = std::size_t(1) << k;

    std::vector<Elem> fa(n, K.zero()), fb(n, K.zero());
    for (std::size_t i = 0; i < A.size(); ++i) fa[i] = A[i];
    for (std::size_t i = 0; i < B.size(); ++i) fb[i] = B[i];

    transform(K, fa, false);
    transform(K, fb, false);
    for (std::size_t i = 0; i < n; ++i) fa[i] = K.mul(fa[i], fb[i]);
    transform(K, fa, true);

    fa.resize(L);
    return fa;
}

}  // namespace alcp::fft
