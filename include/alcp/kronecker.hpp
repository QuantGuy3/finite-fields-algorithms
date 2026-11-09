// alcp/kronecker.hpp — Producto de polinomios por sustitución de Kronecker sobre Z (fase 4, nivel B).
//
// Idea: un polinomio a(x) = Σ a_i x^i con coeficientes enteros pequeños es el
// entero A = a(2^b) si cada coeficiente cabe en b bits.  Si b es lo bastante
// grande para que ningún coeficiente del producto se desborde, el producto de
// polinomios se lee en el producto de enteros A·B, que GMP calcula con sus
// algoritmos rápidos (Toom-Cook y FFT de Schönhage-Strassen).  Después cada
// trozo de b bits se reduce módulo p.
//
//   * F_p (Fp con GMP o Fp64): cada coeficiente ocupa un hueco.  El coeficiente
//     k del producto es Σ a_i b_{k-i} < m·(p-1)^2, con m = min(na, nb).
//   * F_q = F_p[t]/(g) de grado e: cada elemento es un polinomio en t de grado
//     < e; se coloca en un bloque de 2e-1 huecos (Kronecker en dos niveles), de
//     modo que el producto de dos elementos (grado <= 2e-2 en t) cae entero en
//     su bloque.  Cada hueco suma como mucho m·e productos.  Al desempaquetar,
//     cada bloque de 2e-1 coeficientes de F_p se reduce módulo g (los ganchos de
//     acumulación de FqT).
//
// No usa raíces de la unidad ni números aleatorios: no altera la secuencia del
// generador global.  El resultado es exactamente el producto (mismos valores
// que cualquier otro algoritmo de multiplicación).
#pragma once

#include "bigint.hpp"
#include "field.hpp"
#include "fp64.hpp"
#include "fpn.hpp"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

#include <gmp.h>

namespace alcp::ks {

using limb = mp_limb_t;
static_assert(sizeof(limb) == 8, "kronecker.hpp supone limbs de 64 bits");

// ¿Es Field una extensión FqT (tiene cuerpo base y grado de extensión)?
template <class Field>
concept Extension = requires(const Field& K) {
    K.base_field();
    K.extension_degree();
};

namespace detail {

// dst[bit .. bit + 64·nv) |= v[0..nv) (dst con ceros en esa zona).
inline void escribir(limb* dst, std::size_t bit, const limb* v, std::size_t nv) {
    const std::size_t w = bit >> 6, s = bit & 63;
    if (s == 0) {
        for (std::size_t i = 0; i < nv; ++i) dst[w + i] |= v[i];
        return;
    }
    for (std::size_t i = 0; i < nv; ++i) {
        dst[w + i] |= v[i] << s;
        dst[w + i + 1] |= v[i] >> (64 - s);
    }
}

// out[0..nout) = bits [bit, bit + b) de src[0..nsrc) (lo que no existe, a cero).
inline void leer(const limb* src, std::size_t nsrc, std::size_t bit, std::size_t b, limb* out,
                 std::size_t nout) {
    const std::size_t w = bit >> 6, s = bit & 63;
    auto at = [&](std::size_t i) -> limb { return i < nsrc ? src[i] : 0; };
    for (std::size_t i = 0; i < nout; ++i) {
        limb x = at(w + i) >> s;
        if (s) x |= at(w + i + 1) << (64 - s);
        out[i] = x;
    }
    // Quita los bits por encima de b.
    const std::size_t sobra = nout * 64 - b;
    if (sobra > 0 && sobra < 64) out[nout - 1] &= (~limb(0)) >> sobra;
}

// dst[bit .. bit + 64) |= v en dos palabras (dst con ceros en esa zona y con
// la palabra siguiente disponible).  Es escribir() con nv = 1 (fase 6).
inline void poner(limb* dst, std::size_t bit, limb v) {
    const std::size_t w = bit >> 6, s = bit & 63;
    dst[w] |= v << s;
    if (s) dst[w + 1] |= v >> (64 - s);
}

// Bits [bit, bit + b) de src con b <= 64, en una palabra (fase 6).  Lee
// src[bit/64] y, si el hueco cruza la frontera, src[bit/64 + 1]: el llamador
// garantiza que existen.  mask = 2^b - 1 (o ~0 con b = 64).
inline limb tomar(const limb* src, std::size_t bit, limb mask) {
    const std::size_t w = bit >> 6, s = bit & 63;
    limb x = src[w] >> s;
    if (s) x |= src[w + 1] << (64 - s);
    return x & mask;
}
inline limb mascara(std::size_t b) { return b >= 64 ? ~limb(0) : (limb(1) << b) - 1; }

// Primera palabra de x (x < 2^64).  No se usa mpz_get_ui: devuelve unsigned
// long, que en Windows tiene 32 bits y truncaría p.
inline std::uint64_t palabra(const BigInt& x) { return static_cast<std::uint64_t>(mpz_getlimbn(x.get_mpz_t(), 0)); }

inline std::size_t bitlen(const BigInt& x) { return sgn(x) == 0 ? 0 : mpz_sizeinbase(x.get_mpz_t(), 2); }

// Bits de cada hueco para un producto en que cada coeficiente suma a lo sumo m
// productos de elementos < p: 2·bitlen(p-1) + bitlen(m) + 1 (m·(p-1)^2 < 2^bits).
// Fase 6: con Fp64 se calcula con palabras, sin los dos mpz temporales (p - 1 y
// m) que costaban dos reservas de memoria por producto; mismo valor.
template <class Base>
std::size_t bits_hueco(const Base& K, std::size_t m) {
    if constexpr (std::is_same_v<typename Base::Elem, std::uint64_t>) {
        const std::uint64_t pm1 = static_cast<std::uint64_t>(mpz_getlimbn(K.characteristic().get_mpz_t(), 0)) - 1;
        return 2 * static_cast<std::size_t>(std::bit_width(pm1)) + static_cast<std::size_t>(std::bit_width(m)) + 1;
    } else {
        const BigInt& p = K.characteristic();
        return 2 * bitlen(p - 1) + bitlen(BigInt(static_cast<unsigned long>(m))) + 1;
    }
}

// Valor de un elemento de cuerpo primo como limbs (vista sin copia si se puede).
struct Limbs {
    const limb* p;
    std::size_t n;
    limb uno;  // almacenamiento para el caso de una palabra
};
inline Limbs como_limbs(const BigInt& x) {
    return {mpz_limbs_read(x.get_mpz_t()), static_cast<std::size_t>(mpz_size(x.get_mpz_t())), 0};
}
inline Limbs como_limbs(std::uint64_t x) { return {nullptr, x ? 1u : 0u, x}; }
// Fase 7: elemento de FpN2 (dos palabras contiguas, la baja primero).
inline Limbs como_limbs(const alcp::detail::Palabras2& x) {
    return {reinterpret_cast<const limb*>(x.w), x.w[1] ? 2u : (x.w[0] ? 1u : 0u), 0};
}

}  // namespace detail

// Producto de a[0..na) y b[0..nb) (na, nb >= 1) en el cuerpo PRIMO K, cuyos
// elementos son BigInt (Fp), palabras canónicas (Fp64) o pares de palabras
// (FpN2, fase 7).
template <class Field>
std::vector<typename Field::Elem> mul_primo(const Field& K, const typename Field::Elem* a, std::size_t na,
                                            const typename Field::Elem* b, std::size_t nb) {
    using Elem = typename Field::Elem;
    const BigInt& p = K.characteristic();
    const std::size_t m = std::min(na, nb);
    // Cota de cada coeficiente del producto: m·(p-1)^2 < 2^bits.
    const std::size_t bits = detail::bits_hueco(K, m);
    const std::size_t L = na + nb - 1;
    const bool cuadrado = (a == b && na == nb);
    const std::size_t la = (na * bits + 63) / 64 + 1, lb = cuadrado ? 0 : (nb * bits + 63) / 64 + 1;
    std::vector<limb> A(la, 0), B(lb, 0);
    // Fase 6: en Fp64 cada coeficiente es una palabra y se escribe directamente
    // (detail::poner); en Fp, sus limbs.
    auto empaquetar = [&](const Elem* x, std::size_t nx, limb* X) {
        for (std::size_t i = 0; i < nx; ++i) {
            if constexpr (std::is_same_v<Elem, std::uint64_t>) {
                detail::poner(X, i * bits, x[i]);
            } else {
                auto v = detail::como_limbs(x[i]);
                detail::escribir(X, i * bits, v.p ? v.p : &v.uno, v.n);
            }
        }
    };
    empaquetar(a, na, A.data());
    if (!cuadrado) empaquetar(b, nb, B.data());
    // Longitudes reales (sin limbs nulos al final; mpn_mul exige el último no nulo).
    std::size_t na_l = la, nb_l = cuadrado ? la : lb;
    while (na_l > 1 && A[na_l - 1] == 0) --na_l;
    const limb* Bp = cuadrado ? A.data() : B.data();
    while (nb_l > 1 && Bp[nb_l - 1] == 0) --nb_l;
    // Fase 6: dos limbs nulos de más al final, para que detail::tomar pueda leer
    // la palabra siguiente del último hueco sin comprobar el límite (el hueco
    // k empieza en el bit k·bits y L·bits puede pasar de 64·(na_l + nb_l)
    // cuando los coeficientes altos del producto son pequeños).
    std::vector<limb> C(std::max(na_l + nb_l, (L * bits + 63) / 64) + 2, 0);
    if (cuadrado)
        mpn_sqr(C.data(), A.data(), static_cast<mp_size_t>(na_l));
    else if (na_l >= nb_l)
        mpn_mul(C.data(), A.data(), static_cast<mp_size_t>(na_l), Bp, static_cast<mp_size_t>(nb_l));
    else
        mpn_mul(C.data(), Bp, static_cast<mp_size_t>(nb_l), A.data(), static_cast<mp_size_t>(na_l));

    std::vector<Elem> r(L);
    const std::size_t nw = (bits + 63) / 64;
    std::vector<limb> t(nw);
    if constexpr (std::is_same_v<Elem, std::uint64_t>) {
        const std::uint64_t pw = detail::palabra(p);
        // Fase 6 (F_101 y p pequeños en general): sin «u128 % p».
        //   * bits <= 64 (cada hueco cabe en una palabra; con p = 101 hasta
        //     m < 2^49 coeficientes, con p < 2^30 hasta m < 2^3): se lee el
        //     hueco directamente (tomar) y se reduce con Barrett de 64 bits.
        //   * 64 < bits <= 128: hueco de dos palabras, reducido por trozos
        //     (hi mod p con Barrett y después divq; Barrett64::dos_palabras).
        // El valor de cada hueco es el coeficiente exacto (< 2^bits, sin
        // acarreos entre huecos), así que el resto es el mismo que antes.
        const alcp::detail::Barrett64 red(pw);
        if (bits <= 64) {
            const limb mask = detail::mascara(bits);
            for (std::size_t k = 0; k < L; ++k) r[k] = red(detail::tomar(C.data(), k * bits, mask));
            return r;
        }
        if (nw == 2) {
            for (std::size_t k = 0; k < L; ++k) {
                detail::leer(C.data(), C.size(), k * bits, bits, t.data(), nw);
                r[k] = red.dos_palabras(t[1], t[0]);
            }
            return r;
        }
        limb pl = pw;
        std::vector<limb> q(nw);
        for (std::size_t k = 0; k < L; ++k) {
            detail::leer(C.data(), C.size(), k * bits, bits, t.data(), nw);
            std::size_t n = nw;
            while (n > 1 && t[n - 1] == 0) --n;
            r[k] = static_cast<std::uint64_t>(mpn_mod_1(t.data(), static_cast<mp_size_t>(n), pl));
        }
        return r;
    } else if constexpr (requires { K.reduce_limbs(t.data(), nw); }) {
        // Fase 7 (FpN2): cada hueco se reduce con el acumulador de 5 palabras
        // del cuerpo (doble REDC), sin pasar por un mpz.  Mismo resto exacto.
        for (std::size_t k = 0; k < L; ++k) {
            detail::leer(C.data(), C.size(), k * bits, bits, t.data(), nw);
            r[k] = K.reduce_limbs(t.data(), nw);
        }
        return r;
    } else {
        BigInt v;
        for (std::size_t k = 0; k < L; ++k) {
            detail::leer(C.data(), C.size(), k * bits, bits, t.data(), nw);
            mpz_import(v.get_mpz_t(), nw, -1, sizeof(limb), 0, 0, t.data());
            mpz_tdiv_r(v.get_mpz_t(), v.get_mpz_t(), p.get_mpz_t());
            r[k] = v;
        }
        return r;
    }
}

// Producto de a[0..na) y b[0..nb) en la EXTENSIÓN K = F_p[t]/(g) de grado e
// (FqT): Kronecker en dos niveles, bloques de 2e-1 huecos por coeficiente.
template <class Field>
    requires Extension<Field>
std::vector<typename Field::Elem> mul_extension(const Field& K, const typename Field::Elem* a, std::size_t na,
                                                const typename Field::Elem* b, std::size_t nb) {
    using Elem = typename Field::Elem;
    const auto& base = K.base_field();
    using BaseElem = typename std::decay_t<decltype(base)>::Elem;
    const std::size_t e = K.extension_degree();
    const std::size_t bloque = 2 * e - 1;
    const BigInt& p = K.characteristic();
    const std::size_t m = std::min(na, nb) * e;
    const std::size_t bits = detail::bits_hueco(base, m);
    const std::size_t L = na + nb - 1;
    auto empaquetar = [&](const Elem* x, std::size_t nx, std::vector<limb>& X) {
        X.assign((nx * bloque * bits + 63) / 64 + 1, 0);
        for (std::size_t i = 0; i < nx; ++i)
            for (std::size_t j = 0; j < e; ++j) {
                auto v = detail::como_limbs(x[i][j]);
                if (v.n) detail::escribir(X.data(), (i * bloque + j) * bits, v.p ? v.p : &v.uno, v.n);
            }
    };
    std::vector<limb> A, B;
    const bool cuadrado = (a == b && na == nb);
    empaquetar(a, na, A);
    if (!cuadrado) empaquetar(b, nb, B);
    std::size_t na_l = A.size();
    while (na_l > 1 && A[na_l - 1] == 0) --na_l;
    const std::vector<limb>& Bv = cuadrado ? A : B;
    std::size_t nb_l = Bv.size();
    while (nb_l > 1 && Bv[nb_l - 1] == 0) --nb_l;
    // Dos limbs nulos de más para detail::tomar (ver mul_primo; fase 6).
    std::vector<limb> C(std::max(na_l + nb_l, (L * bloque * bits + 63) / 64) + 2, 0);
    if (cuadrado)
        mpn_sqr(C.data(), A.data(), static_cast<mp_size_t>(na_l));
    else if (na_l >= nb_l)
        mpn_mul(C.data(), A.data(), static_cast<mp_size_t>(na_l), Bv.data(), static_cast<mp_size_t>(nb_l));
    else
        mpn_mul(C.data(), Bv.data(), static_cast<mp_size_t>(nb_l), A.data(), static_cast<mp_size_t>(na_l));

    typename Field::Acc acc;
    std::vector<Elem> r;
    r.reserve(L);
    if constexpr (std::is_same_v<std::decay_t<decltype(base)>, alcp::Fp64>) {
        // Fase 6, Fq64 con huecos de <= 64 bits (F_(3^30): 18 bits en grado
        // 200; F_256: 16-24 bits): cada hueco se lee directamente en su
        // acumulador de Fp64 (u128) SIN reducirlo módulo p.  El acumulador es
        // una suma sin reducir por definición y acc_reduce de FqT reduce módulo
        // p (Fp64::acc_reduce acepta cualquier valor de 128 bits) y módulo g,
        // así que el elemento resultante es el mismo; se ahorran la lectura
        // genérica (leer) y el mpn_mod_1 de cada hueco.
        if (bits <= 64) {
            const limb mask = detail::mascara(bits);
            for (std::size_t k = 0; k < L; ++k) {
                K.acc_clear(acc);
                const std::size_t b0 = k * bloque * bits;
                for (std::size_t j = 0; j < bloque; ++j) acc[j] = detail::tomar(C.data(), b0 + j * bits, mask);
                r.push_back(K.acc_reduce(acc));
            }
            return r;
        }
    }
    const std::size_t nw = (bits + 63) / 64;
    std::vector<limb> t(nw);
    BigInt v;
    const limb pw = static_cast<limb>(detail::palabra(p));  // solo se usa si p cabe en una palabra
    for (std::size_t k = 0; k < L; ++k) {
        K.acc_clear(acc);
        for (std::size_t j = 0; j < bloque; ++j) {
            detail::leer(C.data(), C.size(), (k * bloque + j) * bits, bits, t.data(), nw);
            BaseElem c;
            if constexpr (std::is_same_v<BaseElem, std::uint64_t>) {
                std::size_t n = nw;
                while (n > 1 && t[n - 1] == 0) --n;
                c = static_cast<std::uint64_t>(mpn_mod_1(t.data(), static_cast<mp_size_t>(n), pw));
            } else {
                mpz_import(v.get_mpz_t(), nw, -1, sizeof(limb), 0, 0, t.data());
                mpz_tdiv_r(v.get_mpz_t(), v.get_mpz_t(), p.get_mpz_t());
                c = v;
            }
            base.acc_set(acc[j], c);
        }
        r.push_back(K.acc_reduce(acc));
    }
    return r;
}

}  // namespace alcp::ks
