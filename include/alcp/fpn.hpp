// alcp/fpn.hpp — F_p de tamaño fijo en N palabras de máquina (fase 7, propuesta 3).
//
// FpN<2> (= FpN2) es el cuerpo primo para p de DOS palabras: cubre el hueco entre
// Fp64 (p < 2^62, fp64.hpp) y Fp (GMP, cualquier p).  El binding lo elige para
// 2^62 <= p < 2^128 (fits()); el constructor admite además cualquier primo impar
// p < 2^128, lo que permite probarlo con primos pequeños con NTT (193, 257, 577,
// 998244353...) y comparar su consumo del generador con el de Fp.  p = 2 se
// rechaza (Montgomery exige p impar; Fp64 ya cubre p = 2).
//
// Igual que Fp64, modela el MISMO concepto FiniteField que Fp (field.hpp) con
// su misma interfaz pública (from_int(BigInt), to_bigint, cardinality() como
// BigInt, pow con exponente BigInt, fft_capacity/fft_base_root...), así que
// Poly<FpN2>, la FFT genérica y fast.hpp funcionan sin cambios.  Solo cambia la
// representación: un elemento son dos palabras (lo, hi) en lugar de un mpz, y
// ninguna operación reserva memoria (en Fp, cerca del 40 % del tiempo de
// factor() en F_(2^127-1) era malloc; docs/optimizacion/investigacion-fase6.md,
// tema 1).
//
// Comportamiento observable idéntico al de Fp:
//   * el elemento es el representante CANÓNICO en [0, p) (no la forma de
//     Montgomery), así que operator< (orden canónico de factor()), eq, is_one y
//     to_bigint son los mismos que con BigInt;
//   * las mismas excepciones y mensajes (p no primo, invertir 0, 0^k con k < 0);
//   * random() y fft_base_root() hacen EXACTAMENTE las mismas llamadas a
//     rand_below que Fp (una por elemento; rand_below(p - 1) por candidato a
//     raíz), de modo que factor(), equal_degree y la primera tabla FFT siguen
//     el mismo camino.  Los productos de polinomios deciden cuándo piden la
//     raíz con la regla «original» de Poly::mul (L >= 64), independiente del
//     backend de la NTT, así que el generador se consume en el mismo momento.
//
// Reducción: Montgomery con R = 2^128 aplicada DOS veces sobre el
// representante canónico (prototipo de la fase 6, bench/experimentos/fase6/fpn):
//     x mod p = REDC( REDC(x) · (R^2 mod p) )
// La primera REDC da x·R^-1 (mod p); la segunda multiplica por R^2·R^-1 = R.
// Cada REDC son dos pasos de una palabra.  Así se evita la división de 256 bits
// entre 128 y no hace falta convertir los elementos a la forma de Montgomery
// (que haría no triviales el orden, la igualdad y las conversiones).
// pow() sí trabaja en forma de Montgomery dentro de la cadena de cuadrados (una
// sola REDC por producto), como Fp64::pow: mismo resultado, otra representación
// solo dentro del bucle.
//
// Acumulación perezosa (ganchos de field.hpp): Acc = 5 palabras, suma exacta de
// productos < p^2 < 2^256 (caben 2^64 términos sin desbordar la quinta).
// acc_reduce pliega la quinta palabra con 2^256 ≡ R^2 (mod p) y reduce las
// cuatro restantes con la doble REDC, que acepta cualquier valor < 2^256.
//
// FFT: no hay NTT especializada.  fft::tables, fast::dft/idft y
// mul(MulAlgo::FFT) funcionan con la FFT genérica de fft.hpp (mariposas con
// K.mul) y dan lo mismo que Fp, pero la regla Auto de Poly::mul nunca la elige
// para FpN2: pierde con el escolar/Karatsuba y, desde 128 coeficientes, con
// Kronecker (polynomial.hpp, prefer_fft y prefer_kronecker).  Lo que sí se
// reproduce es CUÁNDO Poly<Fp> pediría la raíz (pide_raiz_como_fp), para
// consumir el generador igual.  F_(2^127-1) y la mayoría de primos de 128 bits
// tienen capacidad FFT 1; los de la forma k·2^m + 1 (81·2^121 + 1,
// 2^64 - 2^32 + 1, 43·2^94 + 1) sí la admiten y se prueban en
// tests/test_fase7_fpn.cpp.
//
// Requiere un compilador con unsigned __int128 (GCC/Clang de 64 bits) y una GMP
// con limbs de 64 bits.
#pragma once

#include "bigint.hpp"
#include "errors.hpp"
#include "fft_tables.hpp"
#include "field.hpp"

#include <compare>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace alcp {

static_assert(GMP_NUMB_BITS == 64, "FpN requiere una GMP con limbs de 64 bits");

namespace detail {

// Elemento de FpN<2>: representante canónico en [0, p) en dos palabras,
// w[0] la baja y w[1] la alta (en este orden y sin relleno, como dos limbs de
// GMP: kronecker.hpp lo lee como un número de dos limbs).  operator< compara
// primero la palabra alta: es el orden numérico, el mismo que el de BigInt.
struct Palabras2 {
    std::uint64_t w[2] = {0, 0};
    friend bool operator==(const Palabras2& a, const Palabras2& b) { return a.w[0] == b.w[0] && a.w[1] == b.w[1]; }
    friend bool operator<(const Palabras2& a, const Palabras2& b) {
        return a.w[1] != b.w[1] ? a.w[1] < b.w[1] : a.w[0] < b.w[0];
    }
};

}  // namespace detail

// F_p con p de N palabras.  Solo está implementado N = 2 (la propuesta 3 de la
// fase 6 midió la ganancia para p de 65-128 bits; N = 3, 4 quedan pendientes).
template <std::size_t N>
class FpN {
    static_assert(N == 2, "FpN: solo está implementado N = 2 (FpN2)");
};

template <>
class FpN<2> {
    using u64 = std::uint64_t;
    using u128 = unsigned __int128;

public:
    using Elem = detail::Palabras2;
    // Acumulador perezoso: suma exacta en 5 palabras (w[4] la más alta).
    struct Acc {
        u64 w[5] = {0, 0, 0, 0, 0};
    };

    // Rango que elige el binding (2^62 <= p < 2^128): por debajo, Fp64.
    static constexpr unsigned kMinBits = 63;
    static constexpr unsigned kMaxBits = 128;

    // ¿Es p del rango de FpN2 (63 a 128 bits)?  Solo el tamaño; la primalidad
    // la comprueba el constructor.
    static bool fits(const BigInt& p) {
        if (sgn(p) <= 0) return false;
        const std::size_t b = bit_length(p);
        return b >= kMinBits && b <= kMaxBits;
    }

    explicit FpN(BigInt p) : pbig_(std::move(p)) {
        // Mismo mensaje y misma excepción que Fp para un p no primo.
        if (!is_prime(pbig_))
            throw std::invalid_argument(pbig_.get_str() + " is not prime; cannot build F_p");
        if (bit_length(pbig_) > kMaxBits || !mpz_odd_p(pbig_.get_mpz_t()))
            throw std::invalid_argument("FpN2: p = " + pbig_.get_str() +
                                        " debe ser un primo impar de a lo sumo 128 bits (use Fp)");
        const Elem pe = from_mpz(pbig_.get_mpz_t());
        p0_ = pe.w[0];
        p1_ = pe.w[1];
        P_ = val(pe);
        // -p^(-1) mod 2^64 por iteración de Newton (p·p ≡ 1 mod 8: 3 bits
        // correctos; cada paso los duplica: 3 -> 96).
        u64 x = p0_;
        for (int i = 0; i < 5; ++i) x *= 2 - p0_ * x;
        pinv_ = ~x + 1;
        const BigInt r2 = (BigInt(1) << 256) % pbig_;  // R^2 mod p
        R2_ = from_mpz(r2.get_mpz_t());
        const BigInt r1 = (BigInt(1) << 128) % pbig_;  // R mod p (el 1 de Montgomery)
        R1_ = from_mpz(r1.get_mpz_t());
    }

    // --- field data ---------------------------------------------------------
    const BigInt& characteristic() const { return pbig_; }
    BigInt cardinality() const { return pbig_; }

    // Igualdad de cuerpos: el mismo p.
    bool operator==(const FpN& o) const { return P_ == o.P_; }

    // --- constants & coercion ----------------------------------------------
    Elem zero() const { return {}; }
    Elem one() const { return mk(1); }

    // (a mod p) en [0, p), también para a negativo o enorme (como mpz_fdiv_r).
    Elem from_int(const BigInt& a) const {
        mpz_srcptr z = a.get_mpz_t();
        const std::size_t n = mpz_size(z);
        if (mpz_sgn(z) >= 0 && n <= 4) {
            // Valor no negativo de hasta 4 palabras: la doble REDC acepta
            // cualquier x < 2^256 (y con x < p basta copiarlo).
            u64 t[4] = {0, 0, 0, 0};
            for (std::size_t i = 0; i < n; ++i) t[i] = static_cast<u64>(mpz_getlimbn(z, static_cast<mp_size_t>(i)));
            if (n <= 2) {
                const u128 v = (static_cast<u128>(t[1]) << 64) | t[0];
                if (v < P_) return mk(v);
            }
            return reduce4(t);
        }
        BigInt r;
        mpz_fdiv_r(r.get_mpz_t(), z, pbig_.get_mpz_t());  // en [0, p)
        return from_mpz(r.get_mpz_t());
    }

    // (v mod p) en [0, p) para un entero de máquina con signo, sin pasar por
    // un mpz (conversión directa del binding, fase 7).  Mismo valor que
    // from_int(BigInt(v)).
    Elem from_i64(long long v) const {
        // |v| en u64 (también para LLONG_MIN, cuyo opuesto no cabe en long long).
        const u64 mag = v >= 0 ? static_cast<u64>(v) : 0ULL - static_cast<u64>(v);
        Elem r = mk(mag);
        if (static_cast<u128>(mag) >= P_) {  // solo con p < 2^63
            const u64 t[4] = {mag, 0, 0, 0};
            r = reduce4(t);
        }
        return v >= 0 ? r : neg(r);
    }

    // Misma llamada al generador que Fp::random (misma secuencia).
    Elem random() const { return from_mpz(rand_below(pbig_).get_mpz_t()); }

    // --- arithmetic ---------------------------------------------------------
    // Sin saltos dependientes de los datos (el compilador usa cmov).
    Elem add(const Elem& a, const Elem& b) const {
        const u128 x = val(a), s = x + val(b);
        const bool ge = (s < x) | (s >= P_);  // acarreo de 2^128 o s >= p
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
        mul22(a.w[0], a.w[1], b.w[0], b.w[1], t);
        return reduce4(t);
    }

    // Inverso con mpz_invert sobre vistas de solo lectura de a y de p (las
    // longitudes se normalizan aquí: con p de 63-64 bits la palabra alta es 0
    // y GMP exige que el limb más alto de una vista no sea nulo).  El
    // resultado va a un mpz local (una reserva por inversión, despreciable
    // frente a los ~600 ns del propio mpz_invert); no hay estado compartido,
    // así que es seguro en las regiones paralelas de Karatsuba.
    Elem inv(const Elem& a) const {
        if (is_zero(a)) throw division_by_zero("0 is not invertible in F_p");
        const mp_limb_t la[2] = {a.w[0], a.w[1]}, lp[2] = {p0_, p1_};
        mpz_t za, zp;
        mpz_roinit_n(za, la, a.w[1] ? 2 : 1);
        mpz_roinit_n(zp, lp, p1_ ? 2 : 1);
        BigInt r;
        mpz_invert(r.get_mpz_t(), za, zp);  // en [0, p): el mismo valor que invmod en Fp
        return from_mpz(r.get_mpz_t());
    }
    Elem div(const Elem& a, const Elem& b) const { return mul(a, inv(b)); }

    // a^k para cualquier entero k (negativo -> inverso primero; 0^k con k < 0
    // lanza alcp::division_by_zero), por cuadrados y productos sucesivos.
    Elem pow(const Elem& a, const BigInt& k) const {
        if (sgn(k) < 0) {
            if (is_zero(a)) throw division_by_zero("0 is not invertible in F_p");
            const BigInt mk_ = -k;
            return pow_nonneg(inv(a), mk_.get_mpz_t());
        }
        return pow_nonneg(a, k.get_mpz_t());
    }

    // --- acumulación perezosa (ganchos opcionales, ver field.hpp) -----------
    void acc_clear(Acc& s) const { s = Acc{}; }
    void acc_set(Acc& s, const Elem& a) const { s = Acc{{a.w[0], a.w[1], 0, 0, 0}}; }
    // s += a·b con acarreo exacto en las 5 palabras.  La quinta solo crece en
    // los acarreos de 2^256: harían falta más de 2^64 términos para
    // desbordarla, y cada acumulador recibe a lo sumo uno por coeficiente.
    void acc_fma(Acc& s, const Elem& a, const Elem& b) const {
        const u128 p00 = static_cast<u128>(a.w[0]) * b.w[0], p01 = static_cast<u128>(a.w[0]) * b.w[1];
        const u128 p10 = static_cast<u128>(a.w[1]) * b.w[0], p11 = static_cast<u128>(a.w[1]) * b.w[1];
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
    // Valor de s módulo p, para CUALQUIER contenido de las 5 palabras.
    // Pliegue: t + w4·2^256 ≡ t + w4·R2 (mod p).  w4·R2 < 2^192, así que la
    // suma desborda 2^256 a lo sumo una vez; en ese caso el resto es < 2^192 y
    // la segunda vuelta (w4 = 1) ya no puede desbordar: a lo sumo dos vueltas.
    Elem acc_reduce(Acc& s) const {
        u64 t[4] = {s.w[0], s.w[1], s.w[2], s.w[3]};
        u64 w4 = s.w[4];
        while (w4) {
            const u128 a = static_cast<u128>(w4) * R2_.w[0], b = static_cast<u128>(w4) * R2_.w[1];
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

    // Valor de t[0..n) (n limbs, sin signo) módulo p.  Lo usa el
    // desempaquetado de Kronecker (kronecker.hpp, mul_primo): cada hueco es un
    // coeficiente exacto del producto, de hasta 2·128 + bitlen(m) + 1 bits, que
    // cabe en las 5 palabras del acumulador; si no (imposible en la práctica),
    // se reduce con GMP.
    Elem reduce_limbs(const mp_limb_t* t, std::size_t n) const {
        if (n <= 5) {
            Acc s;
            for (std::size_t i = 0; i < n; ++i) s.w[i] = static_cast<u64>(t[i]);
            return acc_reduce(s);
        }
        mpz_t z;
        while (n > 0 && t[n - 1] == 0) --n;
        mpz_roinit_n(z, t, static_cast<mp_size_t>(n));
        BigInt r;
        mpz_tdiv_r(r.get_mpz_t(), z, pbig_.get_mpz_t());
        return from_mpz(r.get_mpz_t());
    }

    // --- predicates ---------------------------------------------------------
    bool eq(const Elem& a, const Elem& b) const { return a == b; }
    bool is_zero(const Elem& a) const { return (a.w[0] | a.w[1]) == 0; }
    bool is_one(const Elem& a) const { return a.w[0] == 1 && a.w[1] == 0; }

    // --- conversion ---------------------------------------------------------
    BigInt to_bigint(const Elem& a) const {
        BigInt r;
        if (is_zero(a)) return r;
        mp_limb_t* d = mpz_limbs_write(r.get_mpz_t(), 2);
        d[0] = a.w[0];
        d[1] = a.w[1];
        mpz_limbs_finish(r.get_mpz_t(), a.w[1] ? 2 : 1);
        return r;
    }
    std::string to_string(const Elem& a) const { return to_bigint(a).get_str(); }

    // En F_p la raíz p-ésima de a es a (ver Fp::pth_root).
    Elem pth_root(const Elem& a) const { return a; }

    // --- FFT support --------------------------------------------------------
    // Mayor k con 2^k | (p - 1).  p es impar, así que p - 1 = (p1, p0 - 1) sin
    // préstamo.
    unsigned long fft_capacity() const {
        const u64 lo = p0_ - 1;
        if (lo != 0) return static_cast<unsigned long>(__builtin_ctzll(lo));
        return 64UL + static_cast<unsigned long>(__builtin_ctzll(p1_));  // p = p1·2^64 + 1, p1 != 0
    }

    // Mismo procedimiento y MISMO consumo del generador que Fp::fft_base_root:
    // x = rand_below(p - 1) + 1 y candidato x^((p-1)/2^cap) hasta que sea
    // primitivo (cand^(2^(cap-1)) != 1).
    const Elem& fft_base_root() const {
        if (!base_root_ready_) {
            const unsigned long cap = fft_capacity();
            const BigInt cofactor = (pbig_ - 1) >> cap;
            const BigInt half = (cap == 0) ? BigInt(0) : (BigInt(1) << (cap - 1));
            const BigInt pm1 = pbig_ - 1;
            while (true) {
                const BigInt xb = rand_below(pm1) + 1;  // x en [1, p-1]
                const Elem cand = pow(from_mpz(xb.get_mpz_t()), cofactor);
                if (cap == 0 || !is_one(pow(cand, half))) {
                    base_root_ = cand;
                    break;
                }
            }
            base_root_ready_ = true;
        }
        return base_root_;
    }

    mutable FFTCache<Elem> fft_cache_;

private:
    BigInt pbig_;
    u64 p0_ = 0, p1_ = 0;  // palabras de p
    u64 pinv_ = 0;         // -p^(-1) mod 2^64
    u128 P_ = 0;           // p como entero de 128 bits
    Elem R2_{};            // R^2 mod p = 2^256 mod p
    Elem R1_{};            // R mod p = 2^128 mod p (1 en forma de Montgomery)
    mutable bool base_root_ready_ = false;
    mutable Elem base_root_{};

    static u128 val(const Elem& a) { return (static_cast<u128>(a.w[1]) << 64) | a.w[0]; }
    static Elem mk(u128 v) { return Elem{{static_cast<u64>(v), static_cast<u64>(v >> 64)}}; }
    // mpz con 0 <= v < 2^128 (precondición) -> Elem.
    static Elem from_mpz(mpz_srcptr z) {
        const std::size_t n = mpz_size(z);
        return Elem{{n > 0 ? static_cast<u64>(mpz_getlimbn(z, 0)) : 0, n > 1 ? static_cast<u64>(mpz_getlimbn(z, 1)) : 0}};
    }

    // t[0..3] = (a1:a0)·(b1:b0), producto completo de 256 bits.
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

    // REDC de t (4 palabras, cualquier t < 2^256): u = (t + m·p) / R con
    // m = -t·p^(-1) mod R, de modo que u ≡ t·R^(-1) (mod p) y
    // u < (2^256 + R·p) / R = R + p.  Devuelve u en (u2:u1:u0), u2 en {0, 1}.
    // Si t < p·R, u < 2p.
    void redc(const u64 t[4], u64& u0, u64& u1, u64& u2) const {
        u64 t1 = t[1], t2 = t[2], t3 = t[3], t4 = 0;
        u64 m = t[0] * pinv_;
        u128 c = static_cast<u128>(m) * p0_ + t[0];  // palabra baja = 0
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

    // u2·2^128 + (u1:u0) menos p si es >= p (entrada < 2p): valor en [0, p).
    Elem corrige(u64 u0, u64 u1, u64 u2) const {
        const u128 r = (static_cast<u128>(u1) << 64) | u0;
        const bool ge = (u2 != 0) | (r >= P_);
        return mk(r - (ge ? P_ : 0));
    }

    // x mod p para cualquier x = t[0..3] < 2^256 (doble REDC).
    Elem reduce4(const u64 t[4]) const {
        u64 u0, u1, u2;
        redc(t, u0, u1, u2);  // u ≡ x·R^-1 (mod p), u < R + p
        u128 u = (static_cast<u128>(u1) << 64) | u0;
        if (u2) u -= P_;  // u >= R > p: u - p < R (resta módulo 2^128 exacta)
        // Ahora u < R, así que u·R2 < R·p y la segunda REDC da un valor < 2p.
        u64 v[4];
        mul22(static_cast<u64>(u), static_cast<u64>(u >> 64), R2_.w[0], R2_.w[1], v);
        redc(v, u0, u1, u2);
        return corrige(u0, u1, u2);
    }

    // Producto de Montgomery: a·b·R^(-1) mod p para a, b < p (a·b < p·R).
    Elem mont_mul(const Elem& a, const Elem& b) const {
        u64 t[4];
        mul22(a.w[0], a.w[1], b.w[0], b.w[1], t);
        u64 u0, u1, u2;
        redc(t, u0, u1, u2);
        return corrige(u0, u1, u2);
    }

    // a^e (e >= 0) por cuadrados y productos de izquierda a derecha, con r y a
    // en forma de Montgomery dentro del bucle (la misma cadena que Fp64).
    Elem pow_nonneg(const Elem& a, mpz_srcptr e) const {
        if (mpz_sgn(e) == 0) return one();
        const Elem am = mont_mul(a, R2_);  // a·R mod p
        Elem r = R1_;                      // 1·R mod p
        if (mpz_size(e) == 1) {
            // Exponente de una palabra: los bits se leen de un registro.
            const u64 x = static_cast<u64>(mpz_getlimbn(e, 0));
            for (int i = 63 - __builtin_clzll(x); i >= 0; --i) {
                r = mont_mul(r, r);
                if ((x >> i) & 1U) r = mont_mul(r, am);
            }
        } else {
            for (std::size_t i = mpz_sizeinbase(e, 2); i-- > 0;) {
                r = mont_mul(r, r);
                if (mpz_tstbit(e, i)) r = mont_mul(r, am);
            }
        }
        return mont_mul(r, mk(1));  // r·R^(-1) mod p
    }
};

using FpN2 = FpN<2>;

// ¿Son a y b el mismo cuerpo?  El mismo objeto o el mismo p.
inline bool same_field(const FpN2& a, const FpN2& b) { return &a == &b || a == b; }

static_assert(FiniteField<FpN2>);
static_assert(LazyAccumulation<FpN2>);

}  // namespace alcp
