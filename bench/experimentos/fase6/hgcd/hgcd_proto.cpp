// hgcd_proto.cpp — Investigación de la fase 6, tema hgcd: half-GCD en Fp64.
//
// Prototipos FUERA de la librería (no se modifica include/).  Se usan los
// internos (privados, ver "#define private public" más abajo) de Poly<Fp64> (hgcd, hgcd_base_w, aplicar, mul_mat, paso,
// paso_w, div_xk, identidad, mul_no_rng) y fft::word::transform.
//
//   1) "perfil": copia de P::hgcd con cronómetros por categoría y nivel.
//   2) Variantes que dan los MISMOS cocientes (mcd y Bézout idénticos a
//      gcd_classic / gcd_ext_classic; se comprueba con las funciones de la
//      librería):
//        V1  hgcd que devuelve también los restos (A, B) = M·(a, b) y aplica la
//            matriz solo a la parte baja: (A, B) = (A1, B1)·x^m + R·(a mod x^m, b mod x^m);
//            gcd sin el aplicar final de tamaño completo.  Productos de la librería.
//        V2  V1 + matriz 2x2 con transformadas NTT reutilizadas (aplicar: 6 directas
//            + 2 inversas en lugar de 12; mul_mat: 8 directas + 4 inversas en lugar de 24).
//        V3  V2 + en la segunda mitad las transformadas de S se comparten entre
//            S·(partes bajas) y S·QR.
//        W   Strassen-Winograd (7 productos de la librería) en mul_mat.
//        V2e/V3e  V2/V3 + "envoltura": si la longitud L del producto pasa poco de 2^(k-1)
//            (p. ej. 513), convolución cíclica de tamaño 2^(k-1) y los coeficientes
//            altos se calculan directamente y se restan (mul_mat 2,3x más rápido).
//        gcd sin cofactores en el nivel superior (sin S·QR) en todas las V*.
//        gcd_ext: cola por debajo del umbral en sitio (como gcd_ext_classic_w).
//   3) Barrido de tamaños para el punto de corte frente al clásico.
//
// Uso:  hgcd_proto.exe [check|perfil|micro|gcd|ext|todo]   (ALCP_THREADS=1)
//
// Compilación (portable, como el módulo):
//   g++ -std=c++23 -O3 -Iinclude bench/experimentos/fase6/hgcd/hgcd_proto.cpp -lgmpxx -lgmp -static -o build/fase6/hgcd/hgcd_proto.exe
// (Smart App Control bloqueó versiones intermedias: el binario válido de la versión
//  final es build/fase6/hgcd/hgcd_proto3.exe, enlazado con -DRELINK_3 y el mismo código.)
// Medición:  ALCP_THREADS=1 ./build/fase6/hgcd/hgcd_proto3.exe todo   (salida en bench/results/fase6/hgcd_preliminar.txt)
// Los internos de hgcd (Mat2, hgcd_base_w, aplicar, mul_mat, paso, div_xk...)
// son privados de Poly: SOLO en este prototipo se hacen visibles con
// "#define private public" alrededor del #include (la librería no se toca).
// Las cabeceras estándar y de GMP se incluyen antes para no alterarlas.
#include <algorithm>
#include <array>
#include <atomic>
#include <climits>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <type_traits>
#include <gmp.h>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <functional>
#include <future>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <variant>
#include <vector>
#include <gmpxx.h>
#define private public
#include "alcp/alcp.hpp"
#undef private

using namespace alcp;
using P = Poly<Fp64>;
using Mat2 = P::Mat2;
using W = std::uint64_t;
using u128 = unsigned __int128;
using reloj = std::chrono::steady_clock;

// ===========================================================================
// Medición: mejor de 3 tandas de >= 30 ms, rotando entre las entradas.
// ===========================================================================
static W g_check = 0;  // checksum para que el compilador no elimine el trabajo
static void usa(const P& r) { g_check += r.coeffs().size() + (r.coeffs().empty() ? 0 : r.coeffs()[0]); }

template <class Fn>
double medir_us(int nentradas, Fn&& fn) {
    double mejor = 1e300;
    for (int tanda = 0; tanda < 3; ++tanda) {
        int veces = 0;
        auto t0 = reloj::now();
        double dt = 0;
        do {
            fn(veces % nentradas);
            ++veces;
            dt = std::chrono::duration<double, std::micro>(reloj::now() - t0).count();
        } while (dt < 30000);
        mejor = std::min(mejor, dt / veces);
    }
    return mejor;
}

static bool igual(const P& x, const P& y) { return x.coeffs() == y.coeffs(); }
static bool igual(const Mat2& x, const Mat2& y) {
    return igual(x.m00, y.m00) && igual(x.m01, y.m01) && igual(x.m10, y.m10) && igual(x.m11, y.m11);
}

// ===========================================================================
// 1) Perfil: copia literal de P::hgcd (con kHgcdMin como parámetro) y cronómetros.
// ===========================================================================
enum Cat { BASE, APLICAR, PASO, MULMAT, COPIAS, NCAT };
static const char* kNombreCat[NCAT] = {"caso base", "aplicar", "paso", "mul_mat", "copias/div_xk"};
struct Perfil {
    double t[8][NCAT] = {};  // [nivel][categoría] en µs
    long c[8][NCAT] = {};
};
static Perfil* g_pf = nullptr;
struct Crono {
    int niv, cat;
    reloj::time_point t0;
    Crono(int n, int c) : niv(n < 7 ? n : 7), cat(c), t0(reloj::now()) {}
    ~Crono() {
        if (!g_pf) return;
        g_pf->t[niv][cat] += std::chrono::duration<double, std::micro>(reloj::now() - t0).count();
        g_pf->c[niv][cat]++;
    }
};

static Mat2 hgcd_perfil(const P& a, const P& b, int minimo, int niv = 0) {
    const Fp64& K = a.field();
    const int n = a.degree();
    const int m = (n + 1) / 2;
    if (b.is_zero() || b.degree() < m) return P::identidad(K);
    if (n < minimo) {
        Crono c(niv, BASE);
        return P::hgcd_base_w(a, b, m);
    }
    P a0(K), b0(K);
    {
        Crono c(niv, COPIAS);
        a0 = P::div_xk(a, static_cast<std::size_t>(m));
        b0 = P::div_xk(b, static_cast<std::size_t>(m));
    }
    const Mat2 R = hgcd_perfil(a0, b0, minimo, niv + 1);
    P A(K), B(K);
    {
        Crono c(niv, COPIAS);
        A = a;
        B = b;
    }
    {
        Crono c(niv, APLICAR);
        P::aplicar(R, A, B);
    }
    if (B.is_zero() || B.degree() < m) return R;
    Mat2 QR = R;
    {
        Crono c(niv, PASO);
        P::paso(A, B, QR);
    }
    if (B.is_zero() || B.degree() < m) return QR;
    const int k = 2 * m - A.degree();
    if (k < 0) {
        while (!B.is_zero() && B.degree() >= m) P::paso(A, B, QR);
        return QR;
    }
    P c0(K), d0(K);
    {
        Crono c(niv, COPIAS);
        c0 = P::div_xk(A, static_cast<std::size_t>(k));
        d0 = P::div_xk(B, static_cast<std::size_t>(k));
    }
    const Mat2 S = hgcd_perfil(c0, d0, minimo, niv + 1);
    Crono c(niv, MULMAT);
    return P::mul_mat(S, QR);
}

// Copia del gcd de la librería con el umbral como parámetro (y hgcd_perfil).
static P gcd_lib(P a, P b, int umbral, int minimo) {
    while (!b.is_zero()) {
        if (a.degree() > b.degree() && a.degree() >= umbral) {
            const Mat2 R = hgcd_perfil(a, b, minimo);
            P::aplicar(R, a, b);
            if (b.is_zero()) break;
        }
        if (a.degree() < umbral) return P::gcd_classic(std::move(a), std::move(b));
        P r = a % b;
        a = std::move(b);
        b = std::move(r);
    }
    return a.is_zero() ? a : a.monic();
}
// Copia del gcd_ext de la librería con el umbral como parámetro (cache FFT ya llena).
static P::ExtGcd gcd_ext_lib(const P& a, const P& b, int umbral, int minimo) {
    const Fp64& K = a.field();
    if (std::max(a.degree(), b.degree()) < umbral) return P::gcd_ext_classic(a, b);
    P r0 = a, r1 = b;
    Mat2 M = P::identidad(K);
    while (!r1.is_zero()) {
        if (r0.degree() > r1.degree() && r0.degree() >= umbral) {
            const Mat2 R = hgcd_perfil(r0, r1, minimo);
            P::aplicar(R, r0, r1);
            M = P::mul_mat(R, M);
            if (r1.is_zero()) break;
        }
        P::paso(r0, r1, M);
    }
    if (r0.is_zero()) return {P::zero(K), P::zero(K), P::zero(K)};
    const W inv = K.inv(r0.lead());
    return {r0.scalar_mul(inv), M.m00.scalar_mul(inv), M.m01.scalar_mul(inv)};
}

// ===========================================================================
// 2) Variantes
// ===========================================================================
struct Op {
    int minimo = 256;      // grado del caso base de la recursión
    int modo = 0;          // 0 = productos de la librería, 1 = NTT reutilizada, 2 = NTT + compartir S
    bool winograd = false; // mul_mat con Strassen-Winograd (solo modo 0)
    int ntt_min = 64;      // NTT reutilizada solo si las entradas de la matriz tienen >= ntt_min coef.
    bool envolver = false; // si L se pasa poco de 2^(k-1): convolución cíclica de tamaño 2^(k-1) y
                           // corrección de los coeficientes altos calculados directamente
};

// ---- utilidades sobre polinomios ----
static P bajo(const P& p, std::size_t k) {  // p mod x^k
    const auto& c = p.coeffs();
    return P(p.field(), std::vector<W>(c.begin(), c.begin() + static_cast<std::ptrdiff_t>(std::min(k, c.size()))));
}
static P sube_suma(const P& hi, std::size_t k, const P& lo) {  // hi·x^k + lo
    const Fp64& K = hi.field();
    const auto& h = hi.coeffs();
    const auto& l = lo.coeffs();
    if (h.empty()) return lo;
    std::vector<W> r(std::max(h.size() + k, l.size()), 0);
    std::copy(l.begin(), l.end(), r.begin());
    for (std::size_t i = 0; i < h.size(); ++i) r[i + k] = K.add(r[i + k], h[i]);
    return P(K, std::move(r));
}

// ---- NTT en palabras: reducción de Barrett para el producto punto a punto ----
struct Red {
    W p, mu;  // mu = floor((2^64 - 1) / p)
    explicit Red(W p_) : p(p_), mu(~W(0) / p_) {}
    W red(W s) const {  // s mod p para cualquier s < 2^64
        W q = static_cast<W>((static_cast<u128>(s) * mu) >> 64);
        W r = s - q * p;
        r -= p & (0 - static_cast<W>(r >= p));
        r -= p & (0 - static_cast<W>(r >= p));
        return r;
    }
};

using Tablas = FFTTables<W>;
static std::size_t lon(const P& x) { return x.coeffs().size(); }
static std::size_t lprod(const P& x, const P& y) { return (lon(x) && lon(y)) ? lon(x) + lon(y) - 1 : 0; }
// Transformada directa de x mod (x^n - 1) (se pliega si x tiene más de n coeficientes).
static std::vector<W> directa(const Tablas& t, const P& x, std::size_t n) {
    std::vector<W> v(n, 0);
    const auto& c = x.coeffs();
    const W p = t.word.p;
    if (c.size() <= n) {
        std::copy(c.begin(), c.end(), v.begin());
    } else {
        for (std::size_t i = 0; i < c.size(); ++i) {
            W s = v[i & (n - 1)] + c[i];
            v[i & (n - 1)] = s >= p ? s - p : s;
        }
    }
    if (!c.empty()) fft::word::transform(t, v.data(), n, false);
    return v;
}

// Exponente k de la transformada para resultados de longitud <= L.  Con
// envolver, si L - 2^(k-1) <= 2^(k-1)/16 se usa 2^(k-1): la convolución es
// cíclica y los coeficientes j >= 2^(k-1) se calculan directamente (desenrolla).
static unsigned long elige_k(std::size_t L, bool envolver) {
    unsigned long k = fft::ceil_log2(L);
    if (envolver && k >= 6) {
        const std::size_t mitad = std::size_t(1) << (k - 1);
        if (L - mitad <= mitad / 16) --k;
    }
    return k;
}

// Coeficiente j de X·Y calculado directamente (acumulado en 128 bits).
static u128 coef(const P& X, const P& Y, std::size_t j) {
    const auto& x = X.coeffs();
    const auto& y = Y.coeffs();
    if (x.empty() || y.empty() || j > x.size() + y.size() - 2) return 0;
    const std::size_t i0 = j >= y.size() ? j - (y.size() - 1) : 0;
    const std::size_t i1 = std::min(j, x.size() - 1);
    u128 s = 0;
    for (std::size_t i = i0; i <= i1; ++i) s += static_cast<u128>(x[i]) * y[j - i];
    return s;
}

// Transformada inversa de v (tamaño n) y reconstrucción de X1·Y1 + X2·Y2 (Y2/X2 pueden ser nulos):
// si el resultado tiene más de n coeficientes, los altos se calculan directamente
// y se restan de los envueltos.
static P final_(const Fp64& K, const Tablas& t, std::vector<W>& v, std::size_t n, const P& X1, const P& Y1,
                const P& X2, const P& Y2) {
    fft::word::transform(t, v.data(), n, true);
    const std::size_t Lc = std::max(lprod(X1, Y1), lprod(X2, Y2));
    if (Lc <= n) {
        v.resize(Lc);
        return P(K, std::move(v));
    }
    const W p = K.modulus64();
    v.resize(Lc, 0);
    for (std::size_t j = n; j < Lc; ++j) {
        const W tj = static_cast<W>((coef(X1, Y1, j) % p + coef(X2, Y2, j) % p) % p);
        v[j] = tj;
        v[j - n] = v[j - n] >= tj ? v[j - n] - tj : v[j - n] + p - tj;
    }
    return P(K, std::move(v));
}
// out = x1·y1 + x2·y2 (punto a punto, transformadas en [0, p)); p < 2^31 en el camino rápido.
static void fma2(std::vector<W>& out, const std::vector<W>& x1, const std::vector<W>& y1, const std::vector<W>& x2,
                 const std::vector<W>& y2, std::size_t n, W p) {
    out.resize(n);
    if (p < (W(1) << 31)) {
        const Red R(p);
        for (std::size_t i = 0; i < n; ++i) out[i] = R.red(x1[i] * y1[i] + x2[i] * y2[i]);
    } else {
        for (std::size_t i = 0; i < n; ++i)
            out[i] = static_cast<W>((static_cast<u128>(x1[i]) * y1[i] % p + static_cast<u128>(x2[i]) * y2[i] % p) % p);
    }
}
static std::size_t max_ent(const Mat2& M) { return std::max({lon(M.m00), lon(M.m01), lon(M.m10), lon(M.m11)}); }

static bool puede_ntt(const Fp64& K, std::size_t L) { return !K.fft_cache_.empty() && fft::fft_exponent(K, L) >= 0; }

// (a, b) <- M·(a, b), con transformadas reutilizadas (mismo resultado que P::aplicar).
static void aplicar_v(const Op& op, const Mat2& M, P& a, P& b) {
    const Fp64& K = a.field();
    const std::size_t L = std::max({lprod(M.m00, a), lprod(M.m01, b), lprod(M.m10, a), lprod(M.m11, b)});
    if (op.modo == 0 || L == 0 || max_ent(M) < static_cast<std::size_t>(op.ntt_min) || !puede_ntt(K, L)) {
        P::aplicar(M, a, b);
        return;
    }
    const unsigned long k = elige_k(L, op.envolver);
    const std::size_t n = std::size_t(1) << k;
    const Tablas& t = fft::tables(K, k);  // la caché no está vacía: no consume el generador
    const auto fa = directa(t, a, n), fb = directa(t, b, n);
    const auto f00 = directa(t, M.m00, n), f01 = directa(t, M.m01, n), f10 = directa(t, M.m10, n),
               f11 = directa(t, M.m11, n);
    std::vector<W> na, nb;
    fma2(na, f00, fa, f01, fb, n, K.modulus64());
    fma2(nb, f10, fa, f11, fb, n, K.modulus64());
    P ra = final_(K, t, na, n, M.m00, a, M.m01, b);
    P rb = final_(K, t, nb, n, M.m10, a, M.m11, b);
    a = std::move(ra);
    b = std::move(rb);
}

// A·B con transformadas reutilizadas (8 directas + 4 inversas).
static Mat2 mul_mat_ntt(const Mat2& A, const Mat2& B, bool envolver = false) {
    const Fp64& K = A.m00.field();
    const std::size_t L = std::max({lprod(A.m00, B.m00), lprod(A.m01, B.m10), lprod(A.m00, B.m01), lprod(A.m01, B.m11),
                                    lprod(A.m10, B.m00), lprod(A.m11, B.m10), lprod(A.m10, B.m01), lprod(A.m11, B.m11)});
    const unsigned long k = elige_k(L, envolver);
    const std::size_t n = std::size_t(1) << k;
    const Tablas& t = fft::tables(K, k);
    const auto a00 = directa(t, A.m00, n), a01 = directa(t, A.m01, n), a10 = directa(t, A.m10, n),
               a11 = directa(t, A.m11, n);
    const auto b00 = directa(t, B.m00, n), b01 = directa(t, B.m01, n), b10 = directa(t, B.m10, n),
               b11 = directa(t, B.m11, n);
    const W p = K.modulus64();
    std::vector<W> c00, c01, c10, c11;
    fma2(c00, a00, b00, a01, b10, n, p);
    fma2(c01, a00, b01, a01, b11, n, p);
    fma2(c10, a10, b00, a11, b10, n, p);
    fma2(c11, a10, b01, a11, b11, n, p);
    return {final_(K, t, c00, n, A.m00, B.m00, A.m01, B.m10), final_(K, t, c01, n, A.m00, B.m01, A.m01, B.m11),
            final_(K, t, c10, n, A.m10, B.m00, A.m11, B.m10), final_(K, t, c11, n, A.m10, B.m01, A.m11, B.m11)};
}

// Strassen-Winograd (7 productos de la librería, 15 sumas).
static Mat2 mul_mat_winograd(const Mat2& A, const Mat2& B) {
    const P s1 = A.m10 + A.m11, s2 = s1 - A.m00, s3 = A.m00 - A.m10, s4 = A.m01 - s2;
    const P t1 = B.m01 - B.m00, t2 = B.m11 - t1, t3 = B.m11 - B.m01, t4 = t2 - B.m10;
    const P p1 = A.m00.mul_no_rng(B.m00), p2 = A.m01.mul_no_rng(B.m10), p3 = s4.mul_no_rng(B.m11),
            p4 = A.m11.mul_no_rng(t4), p5 = s1.mul_no_rng(t1), p6 = s2.mul_no_rng(t2), p7 = s3.mul_no_rng(t3);
    const P u2 = p1 + p6, u3 = u2 + p7, u4 = u2 + p5;
    return {p1 + p2, u4 + p3, u3 - p4, u3 + p5};
}

static Mat2 mul_mat_v(const Op& op, const Mat2& A, const Mat2& B) {
    const Fp64& K = A.m00.field();
    if (op.modo >= 1 && std::max(max_ent(A), max_ent(B)) >= static_cast<std::size_t>(op.ntt_min)) {
        const std::size_t L = max_ent(A) + max_ent(B);
        if (puede_ntt(K, L)) return mul_mat_ntt(A, B, op.envolver);
    }
    if (op.winograd) return mul_mat_winograd(A, B);
    return P::mul_mat(A, B);
}

// Segunda mitad compartiendo las transformadas de S: (a, b) <- S·(a, b) y SQ = S·Q.
static Mat2 segunda_compartida(const Op& op, const Mat2& S, P& a, P& b, const Mat2& Q) {
    const Fp64& K = a.field();
    const std::size_t L1 = std::max({lprod(S.m00, a), lprod(S.m01, b), lprod(S.m10, a), lprod(S.m11, b)});
    const std::size_t L2 = max_ent(S) + max_ent(Q);
    const std::size_t L = std::max(L1, L2);
    if (L1 == 0 || max_ent(S) < static_cast<std::size_t>(op.ntt_min) || !puede_ntt(K, L) ||
        elige_k(L1, op.envolver) != elige_k(L2, op.envolver)) {  // tamaños distintos: no compensa compartir
        aplicar_v(op, S, a, b);
        return mul_mat_v(op, S, Q);
    }
    const unsigned long k = elige_k(L, op.envolver);
    const std::size_t n = std::size_t(1) << k;
    const Tablas& t = fft::tables(K, k);
    const W p = K.modulus64();
    const auto s00 = directa(t, S.m00, n), s01 = directa(t, S.m01, n), s10 = directa(t, S.m10, n),
               s11 = directa(t, S.m11, n);
    const auto fa = directa(t, a, n), fb = directa(t, b, n);
    const auto q00 = directa(t, Q.m00, n), q01 = directa(t, Q.m01, n), q10 = directa(t, Q.m10, n),
               q11 = directa(t, Q.m11, n);
    std::vector<W> na, nb, c00, c01, c10, c11;
    fma2(na, s00, fa, s01, fb, n, p);
    fma2(nb, s10, fa, s11, fb, n, p);
    fma2(c00, s00, q00, s01, q10, n, p);
    fma2(c01, s00, q01, s01, q11, n, p);
    fma2(c10, s10, q00, s11, q10, n, p);
    fma2(c11, s10, q01, s11, q11, n, p);
    P ra = final_(K, t, na, n, S.m00, a, S.m01, b);
    P rb = final_(K, t, nb, n, S.m10, a, S.m11, b);
    a = std::move(ra);
    b = std::move(rb);
    return {final_(K, t, c00, n, S.m00, Q.m00, S.m01, Q.m10), final_(K, t, c01, n, S.m00, Q.m01, S.m01, Q.m11),
            final_(K, t, c10, n, S.m10, Q.m00, S.m11, Q.m10), final_(K, t, c11, n, S.m10, Q.m01, S.m11, Q.m11)};
}

// ---- hgcd con restos ----
struct Res {
    Mat2 M;
    P A, B;  // (A, B) = M·(a, b)
};

// Caso base: copia de P::hgcd_base_w que además devuelve los restos.
static Res base_rem(const P& a, const P& b, int m) {
    using euclides::W;
    const Fp64& K = a.field();
    const W p = K.modulus64();
    std::vector<W> A(a.coeffs().begin(), a.coeffs().end()), B(b.coeffs().begin(), b.coeffs().end()), Q;
    int da = a.degree(), db = b.degree();
    std::vector<W> m00{1}, m01, m10, m11{1};
    while (db >= m) {
        const int dr = P::paso_w(K, A, da, B, db, &Q);
        std::swap(A, B);
        da = db;
        db = dr;
        while (!Q.empty() && Q.back() == 0) Q.pop_back();
        std::vector<W> n10 = euclides::menos_producto(m00, Q, m10, p);
        std::vector<W> n11 = euclides::menos_producto(m01, Q, m11, p);
        m00 = std::move(m10);
        m01 = std::move(m11);
        m10 = std::move(n10);
        m11 = std::move(n11);
    }
    A.resize(static_cast<std::size_t>(da + 1));
    B.resize(static_cast<std::size_t>(db + 1));
    auto Pv = [&](std::vector<W>& v) { return P(K, std::move(v)); };
    return {{Pv(m00), Pv(m01), Pv(m10), Pv(m11)}, Pv(A), Pv(B)};
}

// Mismo resultado que (M = P::hgcd(a, b), (A, B) = M·(a, b)).  Con con_matriz =
// false (nivel superior de gcd) M no se calcula (se devuelve la identidad).
static Res hgcd_rem(const P& a, const P& b, const Op& op, bool con_matriz = true, int niv = 0) {
    // (los Crono solo acumulan si g_pf != nullptr: perfil de la variante)
    const Fp64& K = a.field();
    const int n = a.degree();
    const int m = (n + 1) / 2;
    if (b.is_zero() || b.degree() < m) return {P::identidad(K), a, b};
    if (n < op.minimo) {
        Crono c(niv, BASE);
        return base_rem(a, b, m);
    }
    const std::size_t um = static_cast<std::size_t>(m);
    P a0(K), b0(K);
    {
        Crono c(niv, COPIAS);
        a0 = P::div_xk(a, um);
        b0 = P::div_xk(b, um);
    }
    Res r1 = hgcd_rem(a0, b0, op, true, niv + 1);
    // (A, B) = (A1, B1)·x^m + R·(a mod x^m, b mod x^m)
    P al(K), bl(K);
    {
        Crono c(niv, COPIAS);
        al = bajo(a, um);
        bl = bajo(b, um);
    }
    {
        Crono c(niv, APLICAR);
        aplicar_v(op, r1.M, al, bl);
    }
    P A(K), B(K);
    {
        Crono c(niv, COPIAS);
        A = sube_suma(r1.A, um, al);
        B = sube_suma(r1.B, um, bl);
    }
    if (B.is_zero() || B.degree() < m) return {std::move(r1.M), std::move(A), std::move(B)};
    Mat2 QR = std::move(r1.M);
    {
        Crono c(niv, PASO);
        if (con_matriz) {
            P::paso(A, B, QR);
        } else {
            P r = A % B;
            A = std::move(B);
            B = std::move(r);
        }
    }
    if (B.is_zero() || B.degree() < m) return {std::move(QR), std::move(A), std::move(B)};
    const int k = 2 * m - A.degree();
    if (k < 0) {  // no debería ocurrir (Thull-Yap); por seguridad, pasos clásicos
        while (!B.is_zero() && B.degree() >= m) P::paso(A, B, QR);
        return {std::move(QR), std::move(A), std::move(B)};
    }
    const std::size_t uk = static_cast<std::size_t>(k);
    P c0(K), d0(K), cl(K), dl(K);
    {
        Crono c(niv, COPIAS);
        c0 = P::div_xk(A, uk);
        d0 = P::div_xk(B, uk);
        cl = bajo(A, uk);
        dl = bajo(B, uk);
    }
    Res r2 = hgcd_rem(c0, d0, op, true, niv + 1);
    Mat2 SQ = P::identidad(K);
    if (!con_matriz) {
        Crono c(niv, APLICAR);
        aplicar_v(op, r2.M, cl, dl);
    } else if (op.modo == 2) {
        Crono c(niv, MULMAT);  // aplicar + mul_mat compartidos: se cuentan en mul_mat
        SQ = segunda_compartida(op, r2.M, cl, dl, QR);
    } else {
        {
            Crono c(niv, APLICAR);
            aplicar_v(op, r2.M, cl, dl);
        }
        Crono c(niv, MULMAT);
        SQ = mul_mat_v(op, r2.M, QR);
    }
    Crono c(niv, COPIAS);
    return {std::move(SQ), sube_suma(r2.A, uk, cl), sube_suma(r2.B, uk, dl)};
}

// gcd: hgcd sin matriz en el nivel superior (restos directamente), un paso, y
// el Euclides clásico en sitio por debajo del umbral.  Mismos cocientes.
static P gcd_v(P a, P b, int umbral, const Op& op) {
    while (!b.is_zero()) {
        if (a.degree() > b.degree() && a.degree() >= umbral) {
            Res r = hgcd_rem(a, b, op, false);
            a = std::move(r.A);
            b = std::move(r.B);
            if (b.is_zero()) break;
        }
        if (a.degree() < umbral) return P::gcd_classic(std::move(a), std::move(b));
        P r = a % b;
        a = std::move(b);
        b = std::move(r);
    }
    return a.is_zero() ? a : a.monic();
}

// gcd_ext: hgcd con restos, mul_mat de la variante y, si cola_en_sitio, los
// pasos por debajo del umbral en sitio sobre palabras (copia del bucle de
// gcd_ext_classic_w partiendo de la matriz acumulada).  Requiere caché FFT
// no vacía (como el camino rápido de la librería).
static P::ExtGcd gcd_ext_v(const P& a, const P& b, int umbral, const Op& op, bool cola_en_sitio) {
    const Fp64& K = a.field();
    if (std::max(a.degree(), b.degree()) < umbral) return P::gcd_ext_classic(a, b);
    P r0 = a, r1 = b;
    Mat2 M = P::identidad(K);
    while (!r1.is_zero()) {
        if (r0.degree() > r1.degree() && r0.degree() >= umbral) {
            Res r = hgcd_rem(r0, r1, op, true);
            r0 = std::move(r.A);
            r1 = std::move(r.B);
            M = mul_mat_v(op, r.M, M);
            if (r1.is_zero()) break;
        }
        if (cola_en_sitio && r0.degree() < umbral) break;
        P::paso(r0, r1, M);
    }
    if (cola_en_sitio && !r1.is_zero()) {
        // Cola en sitio: (R0, R1) con cofactores (X0, X1) = (M00, M10), (Y0, Y1) = (M01, M11).
        const W p = K.modulus64();
        std::vector<W> R0(r0.coeffs()), R1(r1.coeffs()), Q;
        int d0 = r0.degree(), d1 = r1.degree();
        std::vector<W> X0(M.m00.coeffs()), X1(M.m10.coeffs()), Y0(M.m01.coeffs()), Y1(M.m11.coeffs());
        while (d1 >= 0) {
            const int dr = P::paso_w(K, R0, d0, R1, d1, &Q);
            std::swap(R0, R1);
            d0 = d1;
            d1 = dr;
            while (!Q.empty() && Q.back() == 0) Q.pop_back();
            auto actualiza = [&](std::vector<W>& x0, std::vector<W>& x1) {
                std::vector<W> t;
                if (Q.size() <= 64 || x1.empty()) {
                    t = euclides::menos_producto(x0, Q, x1, p);
                } else {
                    const P PT = P(K, x0) - P(K, Q).mul_no_rng(P(K, x1));
                    t = PT.coeffs();
                }
                x0 = std::move(x1);
                x1 = std::move(t);
            };
            actualiza(X0, X1);
            actualiza(Y0, Y1);
        }
        R0.resize(static_cast<std::size_t>(d0) + 1);
        r0 = P(K, std::move(R0));
        M.m00 = P(K, std::move(X0));
        M.m01 = P(K, std::move(Y0));
    }
    if (r0.is_zero()) return {P::zero(K), P::zero(K), P::zero(K)};
    const W inv = K.inv(r0.lead());
    return {r0.scalar_mul(inv), M.m00.scalar_mul(inv), M.m01.scalar_mul(inv)};
}

// ===========================================================================
// Programa
// ===========================================================================
struct Entradas {
    std::vector<P> as, bs;
};
static Entradas entradas(const Fp64& K, int n, int cuantas) {
    Entradas e;
    for (int i = 0; i < cuantas; ++i) {
        e.as.push_back(P::random(K, n));
        e.bs.push_back(P::random(K, n - 1));
    }
    return e;
}

static const Op kV1{256, 0, false, 64}, kV2{256, 1, false, 64}, kV3{256, 2, false, 64}, kW{256, 0, true, 64};
static const Op kV2e{256, 1, false, 64, true}, kV3e{256, 2, false, 64, true};

static void perfil(const Fp64& K) {
    std::printf("\n== Perfil de hgcd (copia de P::hgcd, kHgcdMin = 256), F_998244353, a grado n, b grado n-1 ==\n");
    for (int n : {1024, 2048, 4096}) {
        const int NE = 8;
        Entradas e = entradas(K, n, NE);
        Perfil pf;
        // tiempo total (sin cronómetros) para medir la sobrecarga
        const double total = medir_us(NE, [&](int i) { usa(hgcd_perfil(e.as[i], e.bs[i], 256).m00); });
        // perfil: 3 pasadas por cada entrada, se divide entre el número de llamadas
        g_pf = &pf;
        const int rep = 3;
        auto t0 = reloj::now();
        for (int r = 0; r < rep; ++r)
            for (int i = 0; i < NE; ++i) usa(hgcd_perfil(e.as[i], e.bs[i], 256).m00);
        const double con = std::chrono::duration<double, std::micro>(reloj::now() - t0).count() / (rep * NE);
        g_pf = nullptr;
        // aplicar final de gcd (fuera de hgcd): 4 productos de tamaño completo
        std::vector<Mat2> Rs;
        for (int i = 0; i < NE; ++i) Rs.push_back(P::hgcd(e.as[i], e.bs[i]));
        const double apl = medir_us(NE, [&](int i) {
            P A = e.as[i], B = e.bs[i];
            P::aplicar(Rs[i], A, B);
            usa(B);
        });
        std::printf("\nn = %d: hgcd %.1f us (con cronómetros %.1f us); aplicar final del gcd (fuera de hgcd) %.1f us\n", n,
                    total, con, apl);
        std::printf("  %-6s", "nivel");
        for (int c = 0; c < NCAT; ++c) std::printf(" %22s", kNombreCat[c]);
        std::printf("\n");
        double suma[NCAT] = {};
        for (int lv = 0; lv < 8; ++lv) {
            bool hay = false;
            for (int c = 0; c < NCAT; ++c) hay = hay || pf.c[lv][c];
            if (!hay) continue;
            std::printf("  %-6d", lv);
            for (int c = 0; c < NCAT; ++c) {
                const double t = pf.t[lv][c] / (rep * NE);
                suma[c] += t;
                std::printf(" %13.1f us (%4.1f)", t, double(pf.c[lv][c]) / (rep * NE));
            }
            std::printf("\n");
        }
        double s = 0;
        for (double x : suma) s += x;
        std::printf("  %-6s", "total");
        for (int c = 0; c < NCAT; ++c) std::printf(" %13.1f us %5.1f%%", suma[c], 100 * suma[c] / con);
        std::printf("\n  resto (control, Mat2 copias, returns): %.1f us (%.1f%%)   [(n) = llamadas por hgcd]\n", con - s,
                    100 * (con - s) / con);
        // Perfil de las variantes (hgcd con restos y matriz, como en gcd_ext)
        const std::pair<const char*, Op> vs[] = {{"V1", kV1}, {"V2", kV2}, {"V3", kV3}, {"V2 min128", Op{128, 1, false, 64}}, {"V2e", kV2e}, {"V3e", kV3e}, {"V3e min128", Op{128, 2, false, 64, true}}};
        for (const auto& [nom, op] : vs) {
            const double tv = medir_us(NE, [&](int i) { usa(hgcd_rem(e.as[i], e.bs[i], op, true).B); });
            Perfil q;
            g_pf = &q;
            for (int r = 0; r < rep; ++r)
                for (int i = 0; i < NE; ++i) usa(hgcd_rem(e.as[i], e.bs[i], op, true).B);
            g_pf = nullptr;
            double sq[NCAT] = {};
            for (int lv = 0; lv < 8; ++lv)
                for (int c = 0; c < NCAT; ++c) sq[c] += q.t[lv][c] / (rep * NE);
            std::printf("  %-10s hgcd+restos %8.1f us (%.2fx frente a hgcd de la librería):", nom, tv, total / tv);
            for (int c = 0; c < NCAT; ++c) std::printf(" %s %.0f;", kNombreCat[c], sq[c]);
            std::printf("\n");
        }
    }
}

// Micro: productos de matrices 2x2 aislados con matrices reales de hgcd.
// S = hgcd de un par de grado n/2 (entradas de grado ~n/4) y Q ídem de otro par;
// aplicar: S a (a mod x^(n/2), b mod x^(n/2)) de grado n/2 (como la parte baja).
static void micro(const Fp64& K) {
    std::printf("\n== Micro: producto de matrices 2x2 y aplicar (us), F_998244353, matrices reales de hgcd ==\n");
    std::printf("%6s %10s %10s %10s %10s %10s | %10s %10s %10s | %10s\n", "n", "mul_mat", "NTT reut.", "NTT+envol",
                "Winograd", "x env.", "aplicar", "NTT reut.", "NTT+envol", "mul n/4");
    for (int n : {512, 1024, 2048, 4096, 8192}) {
        const int NE = 8;
        std::vector<Mat2> S, Q;
        std::vector<P> al, bl, f, g;
        for (int i = 0; i < NE; ++i) {
            S.push_back(P::hgcd(P::random(K, n / 2), P::random(K, n / 2 - 1)));
            Q.push_back(P::hgcd(P::random(K, n / 2), P::random(K, n / 2 - 1)));
            al.push_back(P::random(K, n / 2 - 1));
            bl.push_back(P::random(K, n / 2 - 1));
            f.push_back(P::random(K, n / 4));
            g.push_back(P::random(K, n / 4));
        }
        for (int i = 0; i < NE; ++i) {  // comprobación
            if (!igual(mul_mat_ntt(S[i], Q[i]), P::mul_mat(S[i], Q[i])) ||
                !igual(mul_mat_ntt(S[i], Q[i], true), P::mul_mat(S[i], Q[i])) ||
                !igual(mul_mat_winograd(S[i], Q[i]), P::mul_mat(S[i], Q[i]))) {
                std::printf("FALLO mul_mat\n");
                return;
            }
        }
        const double t0 = medir_us(NE, [&](int i) { usa(P::mul_mat(S[i], Q[i]).m11); });
        const double t1 = medir_us(NE, [&](int i) { usa(mul_mat_ntt(S[i], Q[i]).m11); });
        const double t2 = medir_us(NE, [&](int i) { usa(mul_mat_winograd(S[i], Q[i]).m11); });
        const double t3 = medir_us(NE, [&](int i) { usa(mul_mat_ntt(S[i], Q[i], true).m11); });
        const double a0 = medir_us(NE, [&](int i) {
            P x = al[i], y = bl[i];
            P::aplicar(S[i], x, y);
            usa(y);
        });
        const double a1 = medir_us(NE, [&](int i) {
            P x = al[i], y = bl[i];
            aplicar_v(kV2, S[i], x, y);
            usa(y);
        });
        const double a2 = medir_us(NE, [&](int i) {
            P x = al[i], y = bl[i];
            aplicar_v(kV2e, S[i], x, y);
            usa(y);
        });
        const double m1 = medir_us(NE, [&](int i) { usa(f[i].mul_no_rng(g[i])); });
        std::printf("%6d %10.1f %10.1f %10.1f %10.1f %9.2fx | %10.1f %10.1f %10.1f | %10.1f\n", n, t0, t1, t3, t2,
                    t0 / t3, a0, a1, a2, m1);
    }
}

static bool comprobar(const Fp64& K, const char* nombre) {
    bool ok = true;
    const std::vector<Op> ops = {kV1, kV2, kV3, kW, {128, 1, false, 64}, {512, 1, false, 64}, {128, 2, false, 32},
                                 {64, 2, false, 32}, {256, 2, false, 128}, kV2e, kV3e,
                                 {128, 2, false, 32, true}, {64, 1, false, 16, true}};
    for (int n : {130, 300, 777, 1000, 1024, 2048, 3001}) {
        Entradas e = entradas(K, n, 4);
        // casos degenerados: mcd no trivial (factor común de grado n/3)
        P g = P::random(K, n / 3);
        e.as.push_back(e.as[0] * g);
        e.bs.push_back(e.bs[1] * g);
        e.as.push_back(e.as[2]);
        e.bs.push_back(P::random(K, n - 7));
        for (std::size_t i = 0; i < e.as.size(); ++i) {
            const P& a = e.as[i];
            const P& b = e.bs[i];
            const P gc = P::gcd_classic(a, b);
            const auto ec = P::gcd_ext_classic(a, b);
            const Mat2 H = P::hgcd(a, b);
            P HA = a, HB = b;
            P::aplicar(H, HA, HB);
            ok = ok && igual(gcd_lib(a, b, 0, 256), gc);
            for (const Op& op : ops) {
                const Res r = hgcd_rem(a, b, op, true);
                const Res r0 = hgcd_rem(a, b, op, false);
                ok = ok && igual(r.M, H) && igual(r.A, HA) && igual(r.B, HB) && igual(r0.A, HA) && igual(r0.B, HB);
                ok = ok && igual(gcd_v(a, b, 0, op), gc);
                for (bool cola : {false, true}) {
                    const auto ev = gcd_ext_v(a, b, 0, op, cola);
                    ok = ok && igual(ev.g, ec.g) && igual(ev.x, ec.x) && igual(ev.y, ec.y);
                }
                const auto ev = gcd_ext_v(a, b, 512, op, true);
                ok = ok && igual(ev.g, ec.g) && igual(ev.x, ec.x) && igual(ev.y, ec.y);
            }
            const auto el = gcd_ext_lib(a, b, 0, 256);
            ok = ok && igual(el.g, ec.g) && igual(el.x, ec.x) && igual(el.y, ec.y);
            // la librería tal cual (umbral de producción)
            ok = ok && igual(P::gcd(a, b), gc);
        }
    }
    std::printf("Comprobación %s (hgcd M y restos == P::hgcd + P::aplicar; gcd == gcd_classic; gcd_ext == "
                "gcd_ext_classic): %s\n",
                nombre, ok ? "OK" : "FALLO");
    return ok;
}

static void barrido_gcd(const Fp64& K) {
    std::printf("\n== gcd, F_998244353 (us; umbral 0 = half-GCD en todo el tamaño; 8 pares variados) ==\n");
    std::vector<std::array<double, 4>> extra;
    std::printf("%6s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n", "n", "clásico", "lib hgcd", "V1", "V2",
                "V3", "W", "V2e", "V3e", "V3e min128", "V3e min512", "lib prod");
    for (int n : {256, 384, 512, 768, 1000, 1024, 1536, 2048, 3072, 4096, 8192}) {
        const int NE = 8;
        Entradas e = entradas(K, n, NE);
        auto mide = [&](auto&& f) { return medir_us(NE, [&](int i) { usa(f(e.as[i], e.bs[i])); }); };
        const double c = mide([&](const P& a, const P& b) { return P::gcd_classic(a, b); });
        const double l = mide([&](const P& a, const P& b) { return gcd_lib(a, b, 0, 256); });
        const double v1 = mide([&](const P& a, const P& b) { return gcd_v(a, b, 0, kV1); });
        const double v2 = mide([&](const P& a, const P& b) { return gcd_v(a, b, 0, kV2); });
        const double v3 = mide([&](const P& a, const P& b) { return gcd_v(a, b, 0, kV3); });
        const double w = mide([&](const P& a, const P& b) { return gcd_v(a, b, 0, kW); });
        const double v2e = mide([&](const P& a, const P& b) { return gcd_v(a, b, 0, kV2e); });
        const double v3e = mide([&](const P& a, const P& b) { return gcd_v(a, b, 0, kV3e); });
        const double v3a = mide([&](const P& a, const P& b) { return gcd_v(a, b, 0, Op{128, 2, false, 64, true}); });
        const double v3b = mide([&](const P& a, const P& b) { return gcd_v(a, b, 0, Op{512, 2, false, 64, true}); });
        const double lp = mide([&](const P& a, const P& b) { return P::gcd(a, b); });
        double u[3];
        for (int j = 0; j < 3; ++j) {
            const int U = 512 << j;
            u[j] = mide([&](const P& a, const P& b) { return gcd_v(a, b, U, kV3e); });
        }
        extra.push_back({double(n), u[0], u[1], u[2]});
        std::printf("%6d %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n", n, c, l, v1,
                    v2, v3, w, v2e, v3e, v3a, v3b, lp);
        std::fflush(stdout);
    }
    std::printf("(lib prod = P::gcd con el umbral de producción, 4096 en F_998244353)\n");
    std::printf("V3e con umbral exterior U (clásico en sitio por debajo de U):\n%6s %10s %10s %10s\n", "n", "U=512", "U=1024", "U=2048");
    for (const auto& r : extra) std::printf("%6d %10.1f %10.1f %10.1f\n", static_cast<int>(r[0]), r[1], r[2], r[3]);
}

static void barrido_ext(const Fp64& K) {
    std::printf("\n== gcd_ext, F_998244353 (us; umbral 0 salvo indicación; 8 pares variados) ==\n");
    std::vector<double> extra;
    std::printf("%6s %10s %10s %10s %10s %10s %10s %10s %12s %12s %12s %10s\n", "n", "clásico", "lib hgcd", "V1", "V2",
                "V3", "V3+cola", "V3e+cola", "V3e+c min128", "V3e+c u512", "V3e+c u256", "lib prod");
    for (int n : {256, 384, 512, 768, 1000, 1024, 1536, 2048, 4096}) {
        const int NE = 8;
        Entradas e = entradas(K, n, NE);
        auto mide = [&](auto&& f) { return medir_us(NE, [&](int i) { usa(f(e.as[i], e.bs[i]).x); }); };
        const double c = mide([&](const P& a, const P& b) { return P::gcd_ext_classic(a, b); });
        const double l = mide([&](const P& a, const P& b) { return gcd_ext_lib(a, b, 0, 256); });
        const double v1 = mide([&](const P& a, const P& b) { return gcd_ext_v(a, b, 0, kV1, false); });
        const double v2 = mide([&](const P& a, const P& b) { return gcd_ext_v(a, b, 0, kV2, false); });
        const double v3 = mide([&](const P& a, const P& b) { return gcd_ext_v(a, b, 0, kV3, false); });
        const double v3c = mide([&](const P& a, const P& b) { return gcd_ext_v(a, b, 0, kV3, true); });
        const double v3u = mide([&](const P& a, const P& b) { return gcd_ext_v(a, b, 512, kV3e, true); });
        const double v3v = mide([&](const P& a, const P& b) { return gcd_ext_v(a, b, 256, kV3e, true); });
        const double v3e = mide([&](const P& a, const P& b) { return gcd_ext_v(a, b, 0, kV3e, true); });
        const double v3e128 = mide([&](const P& a, const P& b) { return gcd_ext_v(a, b, 0, Op{128, 2, false, 64, true}, true); });
        const double lp = mide([&](const P& a, const P& b) { return P::gcd_ext(a, b); });
        const double l128 = mide([&](const P& a, const P& b) { return gcd_ext_lib(a, b, 128, 256); });
        extra.push_back(l128);
        std::printf("%6d %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %12.1f %12.1f %12.1f %10.1f\n", n, c, l, v1,
                    v2, v3, v3c, v3e, v3e128, v3u, v3v, lp);
        std::fflush(stdout);
    }
    std::printf("(lib prod = P::gcd_ext con el umbral de producción, 2048 en F_998244353)\n");
    std::printf("copia de gcd_ext de la librería con umbral 128 (como calibra_hgcd -DALCP_HGCD_UMBRAL=128):");
    for (double x : extra) std::printf(" %.1f", x);
    std::printf("\n");
}

int main(int argc, char** argv) {
    const std::string modo = argc > 1 ? argv[1] : "todo";
    set_seed(1);
    Fp64 K(998244353), K101(101), K61((BigInt(1) << 61) - 1);
    // Tablas FFT construidas antes (como calibra_hgcd.cpp): el camino rápido de
    // gcd_ext y la NTT de mul_no_rng solo se usan con la caché no vacía.
    (void)(P::random(K, 300) * P::random(K, 300));
    (void)(P::random(K61, 300) * P::random(K61, 300));
    if (modo == "check" || modo == "todo") {
        bool ok = comprobar(K, "F_998244353");
        ok = comprobar(K101, "F_101 (sin NTT: caminos de reserva)") && ok;
        ok = comprobar(K61, "F_(2^61-1)") && ok;
        if (!ok) return 1;
    }
    if (modo == "perfil" || modo == "todo") perfil(K);
    if (modo == "micro" || modo == "todo") micro(K);
    if (modo == "gcd" || modo == "todo") barrido_gcd(K);
    if (modo == "ext" || modo == "todo") barrido_ext(K);
    std::printf("\nchecksum %llu\n", static_cast<unsigned long long>(g_check));
    return 0;
}
