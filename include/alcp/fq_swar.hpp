// alcp/fq_swar.hpp — F_q = F_p[t]/(g) con p pequeño: aritmética «SWAR» con
// cuatro coeficientes de 16 bits por palabra (fase 7, propuesta 7).
//
// Motivo (docs/optimizacion/investigacion-fase6.md, tema 5 «factor»): en
// F_(3^30) el tiempo de factor() es casi todo producto, acumulación e inverso
// de elementos de Fq64, que trabajan coeficiente a coeficiente con
// acumuladores de 128 bits.  Con p pequeño (2·n·(p-1)^2 + p < 2^16) cada
// coeficiente de una convolución cabe en 16 bits, así que se pueden sumar
// cuatro a la vez en una palabra de 64 bits:
//   * producto y acumulación: el segundo factor se empaqueta en huecos de 16
//     bits (y sus cuatro desplazamientos de 0..3 huecos); cada a_i·B suma una
//     fila entera de la convolución con ceil(n/4) + 1 productos de palabra en
//     lugar de n (F_(3^30): 9 en vez de 30).  La reducción módulo g es igual
//     (fase 8: con la tabla empaquetada de t^k mod g, k = n..2n-2), y módulo p
//     se reduce una vez por hueco;
//   * Acc de tamaño fijo (std::array): sin memoria dinámica por coeficiente,
//     con un contador que normaliza los huecos módulo p antes de que puedan
//     desbordar (nunca con las cotas de factor() en F_(3^30));
//   * inverso: el mismo algoritmo de Euclides en F_p[t] de FqT::inv, en
//     arrays de la pila y con la reducción módulo p diferida (una vez por
//     ronda en lugar de una por operación).  El inverso es único: mismo valor.
//   * productos de Kronecker (kronecker.hpp, mul_extension): cada hueco del
//     producto de enteros se lee con una máscara y el bloque de 2n-1 huecos se
//     reduce con la reducción escalar de 64 bits (gancho reducir_huecos).
//     Fase 8: el producto de los huecos se hace con la NTT módulo 998244353
//     (mul_extension_ntt) o con huecos de 16 bits alineados que ya son el
//     acumulador de este cuerpo (mul_extension16), y cada bloque se reduce con
//     reducir_carriles / reducir_bloque16 (F_(3^30), grado 200: de 1,6 a
//     0,7 ms, con carga).
//
// MISMO CUERPO Y MISMOS ELEMENTOS QUE Fq64: FqSwar deriva de Fq64 (Elem =
// std::vector<uint64_t> con los coeficientes en [0, p)), hereda constructores,
// validación, random(), from_vector/from_index/to_index, to_string,
// fft_capacity y fft_base_root (mismo consumo del generador aleatorio), y solo
// sustituye las operaciones caras por versiones con el mismo valor exacto.
// Prototipo: bench/experimentos/fase6/factor/proto2.cpp; con el medidor
// corregido de la fase 7, factor() de grado 5 en F_(3^30) pasa de 3,4 a 1,7 ms.
#pragma once

#include "bigint.hpp"
#include "errors.hpp"
#include "field.hpp"
#include "fp64.hpp"
#include "fq.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace alcp {

class FqSwar : public Fq64 {
public:
    using u64 = std::uint64_t;
    // Grado de extensión máximo (búferes en la pila).
    static constexpr std::size_t kMaxE = 64;
    // Palabras de un elemento empaquetado (+1 por el desplazamiento de 0..3 huecos).
    static constexpr std::size_t kNB = kMaxE / 4 + 1;

    // Acumulador: la convolución sin reducir, en huecos de 16 bits, y cuántos
    // productos lleva desde la última normalización.
    struct Acc {
        std::array<u64, (2 * kMaxE) / 4 + 2> w{};
        u64 cnt = 0;
    };

    // ¿Cabe F_(p^n) en FqSwar?  Cada hueco de una convolución suma como mucho
    // n productos (p-1)^2; con margen para la reducción y la normalización, se
    // exige 2·n·(p-1)^2 + p < 2^16 (F_(3^30): 243; F_(7^4): 295).
    static bool fits(const BigInt& p, unsigned long n) {
        if (sgn(p) <= 0 || n == 0 || n > kMaxE || p >= BigInt(1UL << 16)) return false;
        const u64 pp = mpz_get_ui(p.get_mpz_t());
        return 2 * n * (pp - 1) * (pp - 1) + pp < (u64(1) << 16);
    }

    // Los mismos constructores que Fq64 (mismas excepciones y mismo consumo
    // del generador: se construye primero el Fq64).
    explicit FqSwar(const BigInt& q) : FqSwar(Fq64(q)) {}
    FqSwar(const BigInt& p, unsigned long n) : FqSwar(Fq64(p, n)) {}
    FqSwar(Fp64 base, std::vector<BigInt> g) : FqSwar(Fq64(std::move(base), std::move(g))) {}
    explicit FqSwar(Fq64 ref) : Fq64(std::move(ref)) { construir(); }

    // El mismo cuerpo visto como Fq64 (para el camino genérico de Kronecker).
    const Fq64& como_fq64() const { return *this; }

    // --- aritmética (mismo valor exacto que Fq64) ---------------------------
    Elem mul(const Elem& a, const Elem& b) const {
        std::array<u64, std::tuple_size_v<decltype(Acc::w)>> S{};
        conv16(S.data(), a.data(), b.data());
        return reducir16(S.data());
    }
    Elem div(const Elem& a, const Elem& b) const { return mul(a, inv(b)); }
    // El mismo algoritmo binario de FqT::pow con el producto rápido.
    Elem pow(const Elem& a, const BigInt& k) const {
        Elem base = a;
        BigInt e = k;
        if (e < 0) { base = inv(base); e = -e; }
        Elem result = one();
        while (e > 0) {
            if (mpz_odd_p(e.get_mpz_t())) result = mul(result, base);
            e >>= 1;
            if (e > 0) base = mul(base, base);
        }
        return result;
    }
    Elem pth_root(const Elem& a) const { return pow(a, cardinality() / characteristic()); }

    // Euclides extendido en F_p[t] sobre (g, a), como FqT::inv, con reducción
    // módulo p perezosa: cada entrada de los restos y de los coeficientes de
    // Bézout suma a lo sumo n productos < p^2 entre dos reducciones (cabe de
    // sobra en 64 bits con p < 2^16 y n <= 64).
    Elem inv(const Elem& a) const {
        if (is_zero(a)) throw division_by_zero("0 is not invertible in F_q");
        const int e = static_cast<int>(e_);
        std::array<u64, kMaxE + 1> A{}, B{}, S0{}, S1{};
        u64 *r0 = A.data(), *r1 = B.data(), *s0 = S0.data(), *s1 = S1.data();
        for (int i = 0; i <= e; ++i) r0[i] = g_[i];
        for (int i = 0; i < e; ++i) r1[i] = a[i];
        int d0 = e, d1 = e - 1;
        while (r1[d1] == 0) --d1;
        int l0 = 0, l1 = 1;  // longitudes de s0 (= 0) y s1 (= 1)
        s1[0] = 1;
        while (d1 > 0) {
            const u64 il = base_field().inv(r1[d1]);
            while (d0 >= d1) {
                const u64 lead = red_(r0[d0]);
                if (lead == 0) { r0[d0--] = 0; continue; }
                const u64 nf = p_ - red_(lead * il);  // -cociente (lead·il != 0 mod p)
                const int pos = d0 - d1;
                for (int j = 0; j < d1; ++j) r0[pos + j] += nf * r1[j];
                for (int j = 0; j < l1; ++j) s0[pos + j] += nf * s1[j];
                l0 = std::max(l0, pos + l1);
                r0[d0--] = 0;
            }
            for (int j = 0; j < d1; ++j) r0[j] = red_(r0[j]);
            for (int j = 0; j < l0; ++j) s0[j] = red_(s0[j]);
            d0 = d1 - 1;
            while (d0 >= 0 && r0[d0] == 0) --d0;
            while (l0 > 0 && s0[l0 - 1] == 0) --l0;
            std::swap(r0, r1);
            std::swap(d0, d1);
            std::swap(s0, s1);
            std::swap(l0, l1);
        }
        const u64 ci = base_field().inv(r1[0]);
        Elem r(e_, 0);
        for (int i = 0; i < l1; ++i) r[i] = red_(s1[i] * ci);
        return r;
    }

    // --- acumulación perezosa (ganchos de field.hpp) -------------------------
    void acc_clear(Acc& s) const {
        for (std::size_t i = 0; i < ns_; ++i) s.w[i] = 0;
        s.cnt = 0;
    }
    void acc_set(Acc& s, const Elem& a) const {
        acc_clear(s);
        for (std::size_t i = 0; i < e_; ++i) s.w[i / 4] |= a[i] << (16 * (i % 4));
    }
    void acc_fma(Acc& s, const Elem& a, const Elem& b) const {
        if (++s.cnt > lim16_) {  // normaliza cada hueco módulo p antes de poder desbordar
            for (std::size_t k = 0; k < 2 * e_ - 1; ++k) {
                const u64 v = red_(hueco(s.w.data(), k));
                s.w[k / 4] = (s.w[k / 4] & ~(u64(0xFFFF) << (16 * (k % 4)))) | (v << (16 * (k % 4)));
            }
            s.cnt = 1;
        }
        conv16(s.w.data(), a.data(), b.data());
    }
    Elem acc_reduce(Acc& s) const { return reducir16(s.w.data()); }

    // --- ganchos de kronecker.hpp, fase 8 -------------------------------------
    // Producto de Kronecker con m = min(na, nb)·n términos por hueco: ¿cabe
    // cada hueco en 16 bits, incluidos los (n-1)·(p-1)^2 que le suma la
    // reducción módulo g (reducir16)?  Entonces vale mul_extension16.
    // F_(3^30): hasta m = 16354 (545 coeficientes: 545·30 = 16350); F_(5^10): hasta 4086.
    bool kronecker16(std::size_t m) const { return m <= m16_; }
    // Palabras de un bloque de 2n-1 huecos de 16 bits.
    std::size_t palabras_bloque16() const { return ns_ - 1; }
    // Reduce módulo g y módulo p un bloque de 2n-1 huecos de 16 bits (lee las
    // palabras 0..palabras_bloque16()-1 y no las modifica).
    Elem reducir_bloque16(const u64* S) const { return reducir16(S); }
    // Reduce módulo p y módulo g un bloque de 2n-1 coeficientes de 64 bits
    // (mul_extension_ntt): cada uno se reduce módulo p y se empaqueta en un
    // hueco de 16 bits; después, reducir16 (cada hueco bajo recibe a lo sumo
    // (n-1)·(p-1)^2 más: cabe por fits()).
    Elem reducir_carriles(const u64* s) const {
        std::array<u64, std::tuple_size_v<decltype(Acc::w)>> S{};
        for (std::size_t k = 0; k < 2 * e_ - 1; ++k) S[k / 4] |= red_(s[k]) << (16 * (k % 4));
        return reducir16(S.data());
    }
    // Ancho justo de un hueco con m términos: m·(p-1)^2 < 2^bits (63 si no
    // cabe en una palabra, que el llamador trata como «demasiado ancho»).
    std::size_t bits_justos(std::size_t m) const {
        const u64 pm = (p_ - 1) * (p_ - 1);  // < 2^15 por fits()
        if (m >= (std::size_t(1) << 48)) return 63;
        return static_cast<std::size_t>(std::bit_width(static_cast<u64>(m) * pm));
    }

private:
    u64 p_ = 2;
    std::size_t e_ = 1, nb_ = 1, ns_ = 1;
    u64 lim16_ = 1;
    std::size_t m16_ = 0;  // máximo m de kronecker16
    detail::Barrett64 red_;
    std::array<u64, kMaxE + 1> ng_{}, g_{};          // -g_t y g_t en [0, p)
    // Fase 8: t^k mod g para k = n..2n-2, empaquetado (nb_ palabras por fila).
    std::vector<u64> rt_;

    void construir() {
        e_ = extension_degree();
        if (!fits(characteristic(), e_))
            throw std::invalid_argument("FqSwar: p o n demasiado grandes para huecos de 16 bits (use Fq64)");
        p_ = base_field().modulus64();
        red_ = detail::Barrett64(p_);
        for (std::size_t t = 0; t <= e_; ++t) {
            const u64 gt = mpz_get_ui(modulus()[t].get_mpz_t());  // ya en [0, p) y mónico
            g_[t] = gt;
            if (t < e_) ng_[t] = gt == 0 ? 0 : p_ - gt;
        }
        nb_ = (e_ + 3) / 4;
        ns_ = (2 * e_ - 1 + 3) / 4 + 1;
        const u64 pm = (p_ - 1) * (p_ - 1);
        // Cota de productos entre normalizaciones: un hueco normalizado (< p)
        // más lim16_ convoluciones (<= n·pm cada una) más la reducción módulo g
        // ((n-1)·pm) no llega a 2^16.  Con p = 2 y pm = 1 la cota es grande.
        lim16_ = pm == 0 ? ~u64(0) : ((u64(1) << 16) - 1 - e_ * pm - p_) / (e_ * pm);
        // Fase 8: m·pm + (n-1)·pm <= 2^16 - 1 (pm >= 1 porque p >= 2).
        m16_ = static_cast<std::size_t>(((u64(1) << 16) - 1 - (e_ - 1) * pm) / pm);
        // t^n mod g = -g (sin el término principal); t^(k+1) = t·t^k mod g.
        rt_.assign((e_ - 1) * nb_, 0);
        std::array<u64, kMaxE> cur{};
        for (std::size_t t = 0; t < e_; ++t) cur[t] = ng_[t];
        for (std::size_t k = 0; k + 1 < e_; ++k) {
            for (std::size_t t = 0; t < e_; ++t) rt_[k * nb_ + t / 4] |= cur[t] << (16 * (t % 4));
            const u64 top = cur[e_ - 1];
            for (std::size_t t = e_ - 1; t > 0; --t) cur[t] = red_(cur[t - 1] + top * ng_[t]);
            cur[0] = red_(top * ng_[0]);
        }
    }

    static u64 hueco(const u64* S, std::size_t k) { return (S[k / 4] >> (16 * (k % 4))) & 0xFFFF; }

    // dst[r] = src desplazado r huecos (16·r bits) hacia arriba, r = 0..3.
    void desplazar(const u64* src, std::array<std::array<u64, kNB + 1>, 4>& dst) const {
        for (int r = 0; r < 4; ++r)
            for (std::size_t w = 0; w <= nb_; ++w) {
                const u64 cur = w < nb_ ? src[w] : 0;
                const u64 prev = w > 0 ? src[w - 1] : 0;
                dst[r][w] = r == 0 ? cur : ((cur << (16 * r)) | (prev >> (64 - 16 * r)));
            }
    }
    // S (huecos de 16 bits) += a·b.
    void conv16(u64* __restrict S, const u64* a, const u64* b) const {
        std::array<u64, kNB + 1> bw{};
        for (std::size_t j = 0; j < e_; ++j) bw[j / 4] |= b[j] << (16 * (j % 4));
        std::array<std::array<u64, kNB + 1>, 4> bsh;
        desplazar(bw.data(), bsh);
        const std::size_t nb = nb_;
        for (std::size_t i = 0; i < e_; ++i) {
            const u64 ai = a[i];
            if (ai == 0) continue;
            const u64* __restrict Bp = bsh[i & 3].data();
            u64* __restrict D = S + (i >> 2);
            for (std::size_t w = 0; w <= nb; ++w) D[w] += ai * Bp[w];
        }
    }
    // Reduce la convolución S (2n-1 huecos) módulo g y módulo p; no modifica S.
    // Fase 8: con la tabla rt_ de t^k mod g (k = n..2n-2), el resultado es
    // S_bajo + Σ_k (S_k mod p)·(t^k mod g).  Antes se reducía en cascada
    // (cada hueco alto restaba c·g desplazado y el siguiente hueco dependía
    // del anterior): la cadena de dependencias Barrett -> productos -> hueco
    // siguiente limitaba la velocidad.  Ahora los n-1 coeficientes altos son
    // independientes.  Cada hueco bajo recibe a lo sumo (n-1)·(p-1)^2, la
    // misma cota que en cascada, así que valen lim16_ y m16_.  Mismo elemento
    // (el resto módulo g es único).
    Elem reducir16(const u64* S) const {
        const std::size_t e = e_, nb = nb_;
        std::array<u64, kNB> R;
        for (std::size_t w = 0; w < nb; ++w) R[w] = S[w];
        if (e % 4) R[nb - 1] &= (u64(1) << (16 * (e % 4))) - 1;  // fuera los huecos altos
        const u64* __restrict T = rt_.data();
        for (std::size_t k = e; k < 2 * e - 1; ++k, T += nb) {
            const u64 c = red_(hueco(S, k));
            for (std::size_t w = 0; w < nb; ++w) R[w] += c * T[w];
        }
        Elem r(e);
        for (std::size_t i = 0; i < e; ++i) r[i] = red_(hueco(R.data(), i));
        return r;
    }
};

static_assert(FiniteField<FqSwar>);
static_assert(LazyAccumulation<FqSwar>);

}  // namespace alcp
