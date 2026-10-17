// alcp/fast.hpp — Fast algorithms from algoritmos_rapidos.py (Práctica 3),
// written once and generic over the field.
//
//   * dft / idft               : Cooley-Tukey DFT with an explicitly given root
//                                of unity g of order n = 2^k  (the literal port
//                                of fp_fft/fp_ifft/fq_fft/fq_ifft).
//   * toeplitz_*_vec           : (lower / upper / full) Toeplitz matrix-vector
//                                products, computed via polynomial multiplication
//                                (which itself auto-dispatches schoolbook /
//                                Karatsuba / FFT).
//   * toeplitz_*_inverse       : first column/row of the inverse of a triangular
//                                Toeplitz matrix, by divide-and-conquer.
//   * toeplitz_divmod          : Euclidean division reformulated as a triangular
//                                Toeplitz system.
//
// Note: Karatsuba and the precomputed-twiddle FFT live in polynomial.hpp /
// fft.hpp and are reached automatically through Poly::operator*; these routines
// add the Toeplitz/DFT layer on top.
//
// Validación de argumentos (std::invalid_argument salvo que se indique otra
// cosa), siempre antes de calcular nada:
//   * dft / idft: k < 63 (2^k representable), a.size() == 2^k y g de orden
//     EXACTAMENTE 2^k (g^(2^k) == 1 y, si k >= 1, g^(2^(k-1)) != 1).
//   * toeplitz_lower_vec / toeplitz_upper_vec: n >= 1, a.size() >= n
//     (solo se usan a[0..n-1]) y b.size() == n.
//   * toeplitz_vec: n >= 1, a.size() == 2n-1 y b.size() == n.
//   * toeplitz_lower_inverse / toeplitz_upper_inverse: n >= 1 y a.size() >= n;
//     a[0] == 0 (matriz singular) lanza alcp::division_by_zero.
//   * toeplitz_divmod: f, g y K sobre el mismo cuerpo; g == 0 lanza
//     alcp::division_by_zero.
#pragma once

#include "bigint.hpp"
#include "errors.hpp"
#include "field.hpp"
#include "polynomial.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace alcp::fast {

namespace detail {

// Comprueba los argumentos de dft/idft (ver la cabecera del archivo).  Coste:
// k cuadrados en el cuerpo, despreciable frente a la propia transformada.
template <class Field>
void check_dft_args(const Field& K, const typename Field::Elem& g, unsigned long k,
                    std::size_t len, const char* who) {
    using Elem = typename Field::Elem;
    if (k >= static_cast<unsigned long>(std::numeric_limits<std::size_t>::digits) - 1)
        throw std::invalid_argument(std::string(who) + ": k = " + std::to_string(k) +
                                    " demasiado grande (se requiere k < 63)");
    const std::size_t n = std::size_t(1) << k;
    if (len != n) throw std::invalid_argument(std::string(who) + ": input length must be 2^k");
    // h = g^(2^(k-1)) por cuadrados sucesivos (h = g si k <= 1).
    Elem h = g;
    for (unsigned long i = 1; i < k; ++i) h = K.mul(h, h);
    if (k >= 1) {
        if (K.is_one(h))
            throw std::invalid_argument(std::string(who) + ": la raíz g no tiene orden exactamente 2^" +
                                        std::to_string(k) + " (g^(2^(k-1)) == 1)");
        h = K.mul(h, h);  // h = g^(2^k)
    }
    if (!K.is_one(h))
        throw std::invalid_argument(std::string(who) + ": la raíz g no tiene orden exactamente 2^" +
                                    std::to_string(k) + " (g^(2^k) != 1)");
}

// Recursión de la DFT sin validar (los argumentos se validan una sola vez en
// la llamada pública: validarlos en cada nivel costaría O(n·k) productos).
template <class Field>
std::vector<typename Field::Elem> dft_rec(const Field& K, const typename Field::Elem& g,
                                          unsigned long k, const std::vector<typename Field::Elem>& a) {
    using Elem = typename Field::Elem;
    const std::size_t n = std::size_t(1) << k;
    if (a.size() != n) throw std::invalid_argument("dft: input length must be 2^k");
    if (n == 1) return {a[0]};

    std::vector<Elem> even(n / 2), odd(n / 2);
    for (std::size_t i = 0; i < n / 2; ++i) { even[i] = a[2 * i]; odd[i] = a[2 * i + 1]; }

    Elem g2 = K.mul(g, g);
    std::vector<Elem> P = dft_rec(K, g2, k - 1, even);
    std::vector<Elem> I = dft_rec(K, g2, k - 1, odd);

    std::vector<Elem> A(n);
    Elem w = K.one();
    for (std::size_t j = 0; j < n / 2; ++j) {
        Elem t = K.mul(w, I[j]);
        A[j] = K.add(P[j], t);
        A[j + n / 2] = K.sub(P[j], t);
        w = K.mul(w, g);
    }
    return A;
}

// Validaciones comunes de las rutinas de Toeplitz.
inline void check_toeplitz_n(std::size_t n, const char* who) {
    if (n == 0) throw std::invalid_argument(std::string(who) + ": n debe ser >= 1");
}
inline void check_len_at_least(std::size_t len, std::size_t n, const char* who, const char* what) {
    if (len < n)
        throw std::invalid_argument(std::string(who) + ": " + what + " debe tener al menos n = " +
                                    std::to_string(n) + " elementos (tiene " + std::to_string(len) + ")");
}
inline void check_len_equal(std::size_t len, std::size_t n, const char* who, const char* what) {
    if (len != n)
        throw std::invalid_argument(std::string(who) + ": " + what + " debe tener exactamente " +
                                    std::to_string(n) + " elementos (tiene " + std::to_string(len) + ")");
}

}  // namespace detail

// DFT_{n,g}(a): discrete Fourier transform with the supplied root of unity g
// (which must have order n = 2^k).  Recursive Cooley-Tukey, exactly as in the
// Python reference.
template <class Field>
std::vector<typename Field::Elem> dft(const Field& K, const typename Field::Elem& g,
                                      unsigned long k, const std::vector<typename Field::Elem>& a) {
    detail::check_dft_args(K, g, k, a.size(), "dft");
    return detail::dft_rec(K, g, k, a);
}

// IDFT_{n,g}(a) = n^{-1} * DFT_{n,g^{-1}}(a).
template <class Field>
std::vector<typename Field::Elem> idft(const Field& K, const typename Field::Elem& g,
                                       unsigned long k, const std::vector<typename Field::Elem>& a) {
    using Elem = typename Field::Elem;
    detail::check_dft_args(K, g, k, a.size(), "idft");
    const std::size_t n = std::size_t(1) << k;
    Elem g_inv = K.inv(g);
    // g^{-1} tiene el mismo orden que g: no hace falta volver a validarlo.
    std::vector<Elem> A = detail::dft_rec(K, g_inv, k, a);
    Elem n_inv = K.inv(K.from_int(BigInt(static_cast<unsigned long>(n))));
    for (auto& x : A) x = K.mul(n_inv, x);
    return A;
}

namespace detail {
// First `n` coefficients of a polynomial as a length-n vector (zero-padded).
template <class Field>
std::vector<typename Field::Elem> head(const Field& K, const Poly<Field>& p, std::size_t n) {
    std::vector<typename Field::Elem> r(n, K.zero());
    const auto& c = p.coeffs();
    for (std::size_t i = 0; i < n && i < c.size(); ++i) r[i] = c[i];
    return r;
}
}  // namespace detail

// Lower-triangular Toeplitz (first column a, length n) times vector b.
template <class Field>
std::vector<typename Field::Elem> toeplitz_lower_vec(const Field& K, std::size_t n,
                                                     const std::vector<typename Field::Elem>& a,
                                                     const std::vector<typename Field::Elem>& b) {
    detail::check_toeplitz_n(n, "toeplitz_lower_vec");
    detail::check_len_at_least(a.size(), n, "toeplitz_lower_vec", "a");
    detail::check_len_equal(b.size(), n, "toeplitz_lower_vec", "b");
    Poly<Field> A(K, a), B(K, b);
    return detail::head(K, A * B, n);
}

// Upper-triangular Toeplitz (first row a) times vector b, via the reversal trick.
template <class Field>
std::vector<typename Field::Elem> toeplitz_upper_vec(const Field& K, std::size_t n,
                                                     const std::vector<typename Field::Elem>& a,
                                                     const std::vector<typename Field::Elem>& b) {
    detail::check_toeplitz_n(n, "toeplitz_upper_vec");
    detail::check_len_at_least(a.size(), n, "toeplitz_upper_vec", "a");
    detail::check_len_equal(b.size(), n, "toeplitz_upper_vec", "b");
    std::vector<typename Field::Elem> b_rev(b.rbegin(), b.rend());
    auto y_rev = toeplitz_lower_vec(K, n, a, b_rev);
    std::reverse(y_rev.begin(), y_rev.end());
    return y_rev;
}

// Full Toeplitz times vector.  a has length 2n-1: first row, then first column
// minus the corner (matching the Python layout).
template <class Field>
std::vector<typename Field::Elem> toeplitz_vec(const Field& K, std::size_t n,
                                               const std::vector<typename Field::Elem>& a,
                                               const std::vector<typename Field::Elem>& b) {
    using Elem = typename Field::Elem;
    detail::check_toeplitz_n(n, "toeplitz_vec");
    detail::check_len_equal(a.size(), 2 * n - 1, "toeplitz_vec", "a");
    detail::check_len_equal(b.size(), n, "toeplitz_vec", "b");
    std::vector<Elem> coeffs;
    coeffs.reserve(2 * n - 1);
    for (std::size_t i = n - 1; i >= 1; --i) coeffs.push_back(a[i]);  // reversed a[1..n-1]
    coeffs.push_back(a[0]);
    for (std::size_t i = n; i < 2 * n - 1; ++i) coeffs.push_back(a[i]);

    Poly<Field> A(K, coeffs), B(K, b);
    Poly<Field> C = A * B;
    std::vector<Elem> all = detail::head(K, C, 2 * n - 1);
    return std::vector<Elem>(all.begin() + (n - 1), all.begin() + (2 * n - 1));
}

// First column of the inverse of a lower-triangular Toeplitz matrix whose first
// column is a (length n, a[0] != 0).  Newton-style divide-and-conquer doubling.
// Con n == 0 la recursión no terminaba (k = (0+1)/2 = 0): ahora se rechaza.
template <class Field>
std::vector<typename Field::Elem> toeplitz_lower_inverse(const Field& K, std::size_t n,
                                                         const std::vector<typename Field::Elem>& a) {
    using Elem = typename Field::Elem;
    detail::check_toeplitz_n(n, "toeplitz_lower_inverse");
    detail::check_len_at_least(a.size(), n, "toeplitz_lower_inverse", "a");
    if (n == 1) {
        if (K.is_zero(a[0]))
            throw division_by_zero("toeplitz_lower_inverse: a[0] is zero (matriz singular)");
        return {K.inv(a[0])};
    }
    std::size_t k = (n + 1) / 2;
    std::vector<Elem> a_k(a.begin(), a.begin() + k);
    std::vector<Elem> x_k = toeplitz_lower_inverse(K, k, a_k);

    Poly<Field> Xk(K, x_k), An(K, a);
    Poly<Field> T = An * Xk;
    Poly<Field> two = Poly<Field>::constant(K, K.add(K.one(), K.one()));
    Poly<Field> R = (-T) + two;          // 2 - T
    Poly<Field> Xn = R * Xk;
    return detail::head(K, Xn, n);
}

// For a triangular Toeplitz matrix the inverse is triangular Toeplitz too, and
// the upper case reduces to the lower one (first row <-> first column).
template <class Field>
std::vector<typename Field::Elem> toeplitz_upper_inverse(const Field& K, std::size_t n,
                                                         const std::vector<typename Field::Elem>& a) {
    detail::check_toeplitz_n(n, "toeplitz_upper_inverse");
    detail::check_len_at_least(a.size(), n, "toeplitz_upper_inverse", "a");
    return toeplitz_lower_inverse(K, n, a);
}

// Euclidean division f = g*q + r reformulated as a lower-triangular Toeplitz
// system for the reversed coefficients (port of fp_x_divmod/fq_x_divmod).
template <class Field>
std::pair<Poly<Field>, Poly<Field>> toeplitz_divmod(const Field& K, const Poly<Field>& f,
                                                    const Poly<Field>& g) {
    using Elem = typename Field::Elem;
    if (!alcp::detail::fields_match(K, f.field()) || !alcp::detail::fields_match(K, g.field()))
        throw std::invalid_argument("polinomios sobre cuerpos distintos");
    if (g.is_zero()) throw division_by_zero("toeplitz_divmod: division by zero");
    const int m = f.degree();
    const int n = g.degree();
    if (m < n) return {Poly<Field>::zero(K), f};

    const auto& F = f.coeffs();
    const auto& G = g.coeffs();
    const int d = m - n;

    std::vector<Elem> f_rev(m + 1, K.zero()), g_rev(n + 1, K.zero());
    for (int i = 0; i <= m; ++i) f_rev[m - i] = F[i];
    for (int j = 0; j <= n; ++j) g_rev[n - j] = G[j];

    Elem g0_inv = K.inv(g_rev[0]);
    std::vector<Elem> q_hat(d + 1, K.zero());
    for (int i = 0; i <= d; ++i) {
        Elem s = K.zero();
        int jmax = std::min(i, n);
        for (int j = 1; j <= jmax; ++j) s = K.add(s, K.mul(g_rev[j], q_hat[i - j]));
        q_hat[i] = K.mul(K.sub(f_rev[i], s), g0_inv);
    }
    std::vector<Elem> q_coeffs(d + 1, K.zero());
    for (int k = 0; k <= d; ++k) q_coeffs[k] = q_hat[d - k];

    Poly<Field> q(K, q_coeffs);
    Poly<Field> r = f - g * q;
    return {q, r};
}

}  // namespace alcp::fast
