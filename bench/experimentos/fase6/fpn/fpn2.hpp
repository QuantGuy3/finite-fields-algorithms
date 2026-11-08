// fpn2.hpp — PROTOTIPO (fase 6, tema fpn) de un cuerpo F_p de tamaño fijo con
// p de 2 palabras (2^64 < p < 2^128, p impar), sin reservas de memoria.
//
// No forma parte de la librería: incluye "alcp/alcp.hpp" y modela el concepto
// alcp::FiniteField (y LazyAccumulation), igual que Fp64 en fp64.hpp, para que
// Poly<FpN2> funcione sin tocar include/.
//
// Representación: el elemento es el representante CANÓNICO en [0, p) en dos
// palabras (lo, hi).  Así operator<, eq, is_one y el orden canónico de
// factor() son los mismos que con Fp (BigInt) y el resultado es idéntico bit a
// bit.  La reducción es Montgomery con R = 2^128 aplicada dos veces:
//     x mod p = REDC( REDC(x) · (R^2 mod p) )
// (la primera REDC da x·R^-1, la segunda multiplica por R^2·R^-1 = R).
// Cada REDC son dos pasos de una palabra (4 productos 64x64->128).
//
// Acumulación perezosa: Acc = 5 palabras (suma exacta de productos < 2^256;
// caben 2^64 términos).  acc_reduce pliega la quinta palabra con
// 2^256 ≡ R^2 (mod p) y reduce las cuatro restantes con la doble REDC
// (la REDC acepta aquí cualquier t < 2^256: el resultado es < R + p y se
// corrige con una resta antes de la segunda REDC, cuya entrada sí es < p·R).
//
// Generador aleatorio: random() y fft_base_root() hacen EXACTAMENTE las mismas
// llamadas a rand_below que Fp, de modo que factor() sigue el mismo camino.
#pragma once

#include <gmpxx.h>

#include <cstdint>
#include <stdexcept>
#include <string>

namespace fpn {

using u64 = std::uint64_t;
using u128 = unsigned __int128;

// Elemento: representante canónico en [0, p), palabra baja y alta (en este
// orden y sin relleno: se lee como 2 limbs de GMP).
struct E2 {
    u64 lo = 0, hi = 0;
    friend bool operator==(const E2&, const E2&) = default;
    friend bool operator<(const E2& a, const E2& b) { return a.hi != b.hi ? a.hi < b.hi : a.lo < b.lo; }
    // Solo para que compile el desempaquetado genérico de kronecker.hpp (rama
    // BigInt de mul_primo), que Poly instancia aunque Auto nunca lo elija
    // para este cuerpo (prefer_kronecker es false).  v ya está en [0, p).
    E2& operator=(const mpz_class& v) {
        const std::size_t n = mpz_size(v.get_mpz_t());
        lo = n > 0 ? static_cast<u64>(mpz_getlimbn(v.get_mpz_t(), 0)) : 0;
        hi = n > 1 ? static_cast<u64>(mpz_getlimbn(v.get_mpz_t(), 1)) : 0;
        return *this;
    }
};

}  // namespace fpn

// Declaración adelantada para kronecker.hpp (la llamada está calificada, así
// que la sobrecarga debe verse ANTES de la definición de la plantilla).
namespace alcp::ks::detail {
struct Limbs;
inline Limbs como_limbs(const fpn::E2& x);
}  // namespace alcp::ks::detail

#include "alcp/alcp.hpp"

namespace alcp::ks::detail {
inline Limbs como_limbs(const fpn::E2& x) {
    return {reinterpret_cast<const limb*>(&x.lo), x.hi ? 2u : (x.lo ? 1u : 0u), 0};
}
}  // namespace alcp::ks::detail

namespace fpn {

using alcp::BigInt;

inline u128 val(const E2& a) { return (static_cast<u128>(a.hi) << 64) | a.lo; }
inline E2 mk(u128 v) { return E2{static_cast<u64>(v), static_cast<u64>(v >> 64)}; }

// mpz (0 <= v < 2^128) -> E2
inline E2 from_mpz(mpz_srcptr z) {
    const std::size_t n = mpz_size(z);
    return E2{n > 0 ? static_cast<u64>(mpz_getlimbn(z, 0)) : 0, n > 1 ? static_cast<u64>(mpz_getlimbn(z, 1)) : 0};
}
inline BigInt to_big(const E2& a) {
    BigInt r;
    if (a.hi == 0 && a.lo == 0) return r;
    mp_limb_t* d = mpz_limbs_write(r.get_mpz_t(), 2);
    d[0] = a.lo;
    d[1] = a.hi;
    mpz_limbs_finish(r.get_mpz_t(), a.hi ? 2 : 1);
    return r;
}

class FpN2 {
public:
    using Elem = E2;
    struct Acc {
        u64 w[5] = {0, 0, 0, 0, 0};
    };

    static bool fits(const BigInt& p) {
        const std::size_t b = alcp::bit_length(p);
        return b > 64 && b <= 128 && mpz_odd_p(p.get_mpz_t());
    }

    explicit FpN2(BigInt p) : pbig_(std::move(p)) {
        if (!alcp::is_prime(pbig_))
            throw std::invalid_argument(pbig_.get_str() + " is not prime; cannot build F_p");
        if (!fits(pbig_)) throw std::invalid_argument("FpN2: p debe tener entre 65 y 128 bits");
        const E2 pe = from_mpz(pbig_.get_mpz_t());
        p0_ = pe.lo;
        p1_ = pe.hi;
        P_ = val(pe);
        u64 x = p0_;  // inverso de p módulo 2^64 por Newton (3 -> 96 bits correctos)
        for (int i = 0; i < 5; ++i) x *= 2 - p0_ * x;
        pinv_ = ~x + 1;  // -p^-1 mod 2^64
        BigInt r2 = (BigInt(1) << 256) % pbig_;
        R2_ = from_mpz(r2.get_mpz_t());
        mpz_init2(inv_tmp_, 256);
    }
    FpN2(const FpN2&) = delete;  // el prototipo no necesita copiarse
    ~FpN2() { mpz_clear(inv_tmp_); }

    // --- datos del cuerpo ---------------------------------------------------
    const BigInt& characteristic() const { return pbig_; }
    BigInt cardinality() const { return pbig_; }
    bool operator==(const FpN2& o) const { return P_ == o.P_; }

    // --- constantes y coerción ----------------------------------------------
    Elem zero() const { return {}; }
    Elem one() const { return E2{1, 0}; }
    Elem from_int(const BigInt& a) const {
        BigInt r;
        mpz_fdiv_r(r.get_mpz_t(), a.get_mpz_t(), pbig_.get_mpz_t());
        return from_mpz(r.get_mpz_t());
    }
    // Misma llamada al generador que Fp::random.
    Elem random() const { return from_mpz(alcp::rand_below(pbig_).get_mpz_t()); }

    // --- aritmética (sin saltos dependientes de los datos en add/sub) --------
    Elem add(const Elem& a, const Elem& b) const {
        const u128 x = val(a), s = x + val(b);
        const bool ge = (s < x) | (s >= P_);  // acarreo o s >= p
        return mk(s - (ge ? P_ : 0));
    }
    Elem sub(const Elem& a, const Elem& b) const {
        const u128 x = val(a), y = val(b);
        return mk(x - y + (x < y ? P_ : 0));
    }
    Elem neg(const Elem& a) const {
        const u128 x = val(a);
        return mk(x == 0 ? 0 : P_ - x);
    }
    Elem mul(const Elem& a, const Elem& b) const {
        u64 t[4];
        mul22(a.lo, a.hi, b.lo, b.hi, t);
        return reduce4(t);
    }

    // --- acumulación perezosa -------------------------------------------------
    void acc_clear(Acc& s) const { s = Acc{}; }
    void acc_set(Acc& s, const Elem& a) const { s = Acc{{a.lo, a.hi, 0, 0, 0}}; }
    void acc_fma(Acc& s, const Elem& a, const Elem& b) const {
        const u128 p00 = static_cast<u128>(a.lo) * b.lo, p01 = static_cast<u128>(a.lo) * b.hi;
        const u128 p10 = static_cast<u128>(a.hi) * b.lo, p11 = static_cast<u128>(a.hi) * b.hi;
        u128 c = static_cast<u128>(s.w[0]) + static_cast<u64>(p00);
        s.w[0] = static_cast<u64>(c);
        c = (c >> 64) + (p00 >> 64) + static_cast<u64>(p01) + static_cast<u64>(p10) + s.w[1];
        s.w[1] = static_cast<u64>(c);
        c = (c >> 64) + (p01 >> 64) + (p10 >> 64) + static_cast<u64>(p11) + s.w[2];
        s.w[2] = static_cast<u64>(c);
        c = (c >> 64) + (p11 >> 64) + s.w[3];
        s.w[3] = static_cast<u64>(c);
        s.w[4] += static_cast<u64>(c >> 64);
    }
    Elem acc_reduce(Acc& s) const {
        u64 t[4] = {s.w[0], s.w[1], s.w[2], s.w[3]};
        u64 w4 = s.w[4];
        // Pliegue: t + w4·2^256 ≡ t + w4·R2 (mod p); un acarreo vuelve a sumar R2.
        while (w4) {
            const u128 a = static_cast<u128>(w4) * R2_.lo, b = static_cast<u128>(w4) * R2_.hi;
            u128 c = static_cast<u128>(t[0]) + static_cast<u64>(a);
            t[0] = static_cast<u64>(c);
            c = (c >> 64) + (a >> 64) + static_cast<u64>(b) + t[1];
            t[1] = static_cast<u64>(c);
            c = (c >> 64) + (b >> 64) + t[2];
            t[2] = static_cast<u64>(c);
            c = (c >> 64) + t[3];
            t[3] = static_cast<u64>(c);
            w4 = static_cast<u64>(c >> 64);
        }
        return reduce4(t);
    }

    // Inverso con mpz_invert sobre vistas de solo lectura y un mpz de trabajo
    // ya reservado (sin reservas por llamada).
    Elem inv(const Elem& a) const {
        if (a.lo == 0 && a.hi == 0) throw alcp::division_by_zero("0 is not invertible in F_p");
        const mp_limb_t la[2] = {a.lo, a.hi}, lp[2] = {p0_, p1_};
        mpz_t za, zp;
        mpz_roinit_n(za, la, a.hi ? 2 : 1);
        mpz_roinit_n(zp, lp, 2);
        mpz_invert(inv_tmp_, za, zp);
        return from_mpz(inv_tmp_);
    }
    Elem div(const Elem& a, const Elem& b) const { return mul(a, inv(b)); }

    // a^k (k < 0: inverso primero), cuadrados y productos de izquierda a derecha.
    Elem pow(const Elem& a, const BigInt& k) const {
        Elem base = a;
        mpz_srcptr e = k.get_mpz_t();
        BigInt negk;
        if (mpz_sgn(e) < 0) {
            if (is_zero(a)) throw alcp::division_by_zero("0 is not invertible in F_p");
            base = inv(a);
            negk = -k;
            e = negk.get_mpz_t();
        }
        Elem r = one();
        for (std::size_t i = mpz_sgn(e) ? mpz_sizeinbase(e, 2) : 0; i-- > 0;) {
            r = mul(r, r);
            if (mpz_tstbit(e, i)) r = mul(r, base);
        }
        return r;
    }

    // --- predicados -------------------------------------------------------------
    bool eq(const Elem& a, const Elem& b) const { return a == b; }
    bool is_zero(const Elem& a) const { return (a.lo | a.hi) == 0; }
    bool is_one(const Elem& a) const { return a.lo == 1 && a.hi == 0; }

    // --- conversión ---------------------------------------------------------------
    BigInt to_bigint(const Elem& a) const { return to_big(a); }
    std::string to_string(const Elem& a) const { return to_big(a).get_str(); }
    Elem pth_root(const Elem& a) const { return a; }

    // --- FFT (mismo procedimiento y mismo consumo del generador que Fp) --------
    unsigned long fft_capacity() const { return alcp::v2(pbig_ - 1); }
    const Elem& fft_base_root() const {
        if (!base_root_ready_) {
            const unsigned long cap = fft_capacity();
            const BigInt cofactor = (pbig_ - 1) / (BigInt(1) << cap);
            const BigInt half = (cap == 0) ? BigInt(0) : (BigInt(1) << (cap - 1));
            while (true) {
                BigInt x = alcp::rand_below(pbig_ - 1) + 1;
                BigInt cand = alcp::powmod(x, cofactor, pbig_);
                if (cap == 0 || alcp::powmod(cand, half, pbig_) != 1) {
                    base_root_ = from_mpz(cand.get_mpz_t());
                    break;
                }
            }
            base_root_ready_ = true;
        }
        return base_root_;
    }
    mutable alcp::FFTCache<Elem> fft_cache_;

private:
    BigInt pbig_;
    u64 p0_ = 0, p1_ = 0, pinv_ = 0;
    u128 P_ = 0;
    E2 R2_{};  // 2^256 mod p
    mutable mpz_t inv_tmp_;
    mutable bool base_root_ready_ = false;
    mutable Elem base_root_{};

    // t[0..3] = (a1:a0)·(b1:b0)
    static void mul22(u64 a0, u64 a1, u64 b0, u64 b1, u64 t[4]) {
        const u128 p00 = static_cast<u128>(a0) * b0, p01 = static_cast<u128>(a0) * b1;
        const u128 p10 = static_cast<u128>(a1) * b0, p11 = static_cast<u128>(a1) * b1;
        t[0] = static_cast<u64>(p00);
        u128 c = (p00 >> 64) + static_cast<u64>(p01) + static_cast<u64>(p10);
        t[1] = static_cast<u64>(c);
        c = (c >> 64) + (p01 >> 64) + (p10 >> 64) + static_cast<u64>(p11);
        t[2] = static_cast<u64>(c);
        t[3] = static_cast<u64>((c >> 64) + (p11 >> 64));
    }

    // REDC de t (4 palabras, t < 2^256): devuelve u = (t + m·p)/R < R + p en
    // (u2:u1:u0) con u2 en {0, 1}.
    void redc(const u64 t[4], u64& u0, u64& u1, u64& u2) const {
        u64 t0 = t[0], t1 = t[1], t2 = t[2], t3 = t[3], t4 = 0;
        u64 m = t0 * pinv_;
        u128 c = static_cast<u128>(m) * p0_ + t0;  // palabra baja = 0
        c = (c >> 64) + static_cast<u128>(m) * p1_ + t1;
        t1 = static_cast<u64>(c);
        c = (c >> 64) + t2;
        t2 = static_cast<u64>(c);
        c = (c >> 64) + t3;
        t3 = static_cast<u64>(c);
        t4 = static_cast<u64>(c >> 64);
        m = t1 * pinv_;
        c = static_cast<u128>(m) * p0_ + t1;  // palabra baja = 0
        c = (c >> 64) + static_cast<u128>(m) * p1_ + t2;
        t2 = static_cast<u64>(c);
        c = (c >> 64) + t3;
        t3 = static_cast<u64>(c);
        t4 += static_cast<u64>(c >> 64);
        u0 = t2;
        u1 = t3;
        u2 = t4;
    }

    // x mod p para cualquier x = t[0..3] < 2^256.
    Elem reduce4(const u64 t[4]) const {
        u64 u0, u1, u2;
        redc(t, u0, u1, u2);  // u ≡ x·R^-1 (mod p), u < R + p
        u128 u = (static_cast<u128>(u1) << 64) | u0;
        if (u2) u -= P_;  // u >= R > p: u - p < R (resta módulo 2^128 exacta)
        // Ahora u < R, así que u·R2 < p·R y la segunda REDC da un valor < 2p.
        u64 v[4];
        mul22(static_cast<u64>(u), static_cast<u64>(u >> 64), R2_.lo, R2_.hi, v);
        redc(v, u0, u1, u2);
        const u128 r = (static_cast<u128>(u1) << 64) | u0;
        const bool ge = (u2 != 0) | (r >= P_);  // valor real = u2·2^128 + r
        return mk(r - (ge ? P_ : 0));
    }
};

inline bool same_field(const FpN2& a, const FpN2& b) { return &a == &b || a == b; }

static_assert(alcp::FiniteField<FpN2>);
static_assert(alcp::LazyAccumulation<FpN2>);

}  // namespace fpn
