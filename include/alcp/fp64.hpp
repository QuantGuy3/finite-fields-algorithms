// alcp/fp64.hpp — F_p con p < 2^62 en palabras de máquina (backend de 64 bits).
//
// Fp64 modela el MISMO concepto FiniteField que Fp (field.hpp) y tiene su misma
// interfaz pública (from_int(BigInt), to_bigint, cardinality() como BigInt,
// pow con exponente BigInt, fft_capacity/fft_base_root...), así que Poly<Fp64>,
// la FFT, fast.hpp y la extensión FqT<Fp64> (= Fq64, fq.hpp) funcionan sin
// cambios.  Solo cambia la representación: un elemento es un std::uint64_t en
// [0, p) en lugar de un mpz.
//
// Comportamiento observable idéntico al de Fp:
//   * los mismos valores (representante canónico en [0, p)), luego el mismo
//     orden (operator< de Elem), la misma igualdad y las mismas conversiones;
//   * las mismas excepciones (p no primo, invertir 0, 0^k con k < 0);
//   * random() y fft_base_root() consumen el generador global EXACTAMENTE igual
//     que Fp (una llamada a rand_below(p), resp. rand_below(p - 1), por
//     elemento), de modo que la factorización y Fq(p, n) aleatorios siguen el
//     mismo camino y dan el mismo resultado.
//
// Elección de la reducción: Barrett (no Montgomery) para los elementos.  Con Montgomery los
// elementos vivirían en la forma a·R mod p: habría que convertir en cada
// entrada/salida y, sobre todo, operator<, eq, is_one, to_bigint y el orden
// canónico de factor() dejarían de ser triviales.  Barrett trabaja sobre el
// representante canónico, cuesta dos productos de 64x64->128 bits y uno de
// 64x64->64 por multiplicación, sin divisiones, y su resultado es exacto tras a
// lo sumo dos restas condicionales.
//
// Excepción interna: pow() hace la MISMA cadena de cuadrados y productos, pero
// con los valores intermedios en forma de Montgomery (a·2^64 mod p), que
// convierte al entrar y al salir (p impar).  El resultado es el mismo; solo
// cambia la representación dentro del bucle, donde REDC tiene menos latencia
// que Barrett y la cadena de cuadrados es secuencial.
//
//   k = bit_length(p) <= 62,  m = floor(2^(2k) / p)  (m <= 2^(k+1) <= 2^63)
//   x = a·b < p^2 < 2^(2k)
//   q = floor( floor(x / 2^(k-1)) · m / 2^(k+1) )  ->  floor(x/p) - 2 <= q <= floor(x/p)
//   r = x - q·p  in [0, 3p)  (3p < 2^64: se calcula en aritmética de 64 bits)
//
// Acumulación perezosa (ganchos de field.hpp): Acc = unsigned __int128.  Cada
// producto es < 2^124, así que caben al menos 16 sin reducir; si una suma
// desborda 2^128 (solo posible con p > 2^60 y muchos términos) se detecta por
// el acarreo y se pliega: t + 2^128 ≡ (t mod p) + (2^128 mod p).
// Fase 6 (q pequeño): con p < 2^32 el producto es de 64 bits y no se vigila el
// acarreo (cota en acc_fma), y una suma < 2^64 se reduce con Barrett de 64 bits
// (detail::Barrett64) en lugar de dividir.
//
// Requiere un compilador con unsigned __int128 (GCC/Clang de 64 bits) y una GMP
// con limbs de 64 bits.
#pragma once

#include "bigint.hpp"
#include "errors.hpp"
#include "fft_tables.hpp"
#include "field.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>

namespace alcp {

static_assert(GMP_NUMB_BITS == 64, "Fp64 requiere una GMP con limbs de 64 bits");

namespace detail {

using u128 = unsigned __int128;

// Valor de un mpz con 0 <= v < 2^64 (precondición) como entero de 64 bits.
inline std::uint64_t mpz_to_u64(mpz_srcptr v) { return static_cast<std::uint64_t>(mpz_getlimbn(v, 0)); }

// Entero de 64 bits -> mpz (independiente del tamaño de unsigned long, que es
// de 32 bits en Windows).
inline void mpz_set_u64(mpz_ptr z, std::uint64_t x) {
    if (x == 0) {
        mpz_set_ui(z, 0);
        return;
    }
    mp_limb_t* d = mpz_limbs_write(z, 1);
    d[0] = static_cast<mp_limb_t>(x);
    mpz_limbs_finish(z, 1);
}

inline BigInt u64_to_big(std::uint64_t x) {
    BigInt r;
    mpz_set_u64(r.get_mpz_t(), x);
    return r;
}

// Reducción de Barrett para CUALQUIER x < 2^64 (fase 6), sin división.
//   m = floor((2^64 - 1) / p)  ->  m·p >= 2^64 - p, luego x·m / 2^64 >= x/p - x/2^64 > x/p - 1
//   q = floor(x·m / 2^64) > x/p - 2  ->  r = x - q·p en [0, 2p): basta una resta.
// La usan los acumuladores de Fp64 cuando la suma no cabe en 2^(2k) pero sí en
// una palabra, y el desempaquetado de Kronecker (kronecker.hpp), que antes
// dividían (u128 % p o divq).
struct Barrett64 {
    std::uint64_t p = 2, m = 0;
    Barrett64() = default;
    explicit Barrett64(std::uint64_t pp) : p(pp), m(~std::uint64_t(0) / pp) {}
    std::uint64_t operator()(std::uint64_t x) const {
        const std::uint64_t q = static_cast<std::uint64_t>((static_cast<u128>(x) * m) >> 64);
        const std::uint64_t r = x - q * p;
        return r >= p ? r - p : r;
    }
    // (hi·2^64 + lo) mod p para cualquier valor de dos palabras: primero hi mod
    // p (Barrett) y luego la división 128/64 del hardware, cuyo cociente cabe
    // porque hi < p.  Evita el «u128 % p» de software (__umodti3).
    std::uint64_t dos_palabras(std::uint64_t hi, std::uint64_t lo) const {
        if (hi == 0) return (*this)(lo);
        hi = (*this)(hi);
#if defined(__x86_64__)
        std::uint64_t q, r;
        __asm__("divq %4" : "=a"(q), "=d"(r) : "0"(lo), "1"(hi), "rm"(p) : "cc");
        (void)q;
        return r;
#else
        return static_cast<std::uint64_t>(((static_cast<u128>(hi) << 64) | lo) % p);
#endif
    }
};

}  // namespace detail

class Fp64 {
public:
    using Elem = std::uint64_t;
    using Acc = detail::u128;

    // Mayor número de bits de p que admite el backend (p < 2^62).
    static constexpr unsigned kMaxBits = 62;

    // ¿Se puede construir Fp64(p)?  (Solo el tamaño; la primalidad la comprueba
    // el constructor.)  p >= 2 y p < 2^62.
    static bool fits(const BigInt& p) { return sgn(p) > 0 && bit_length(p) <= kMaxBits; }

    explicit Fp64(BigInt p) : pbig_(std::move(p)) {
        // Mismo mensaje y misma excepción que Fp para un p no primo.
        if (!is_prime(pbig_))
            throw std::invalid_argument(pbig_.get_str() + " is not prime; cannot build F_p");
        if (!fits(pbig_))
            throw std::invalid_argument("Fp64: p = " + pbig_.get_str() + " no cabe en " +
                                        std::to_string(kMaxBits) + " bits (use Fp)");
        p_ = detail::mpz_to_u64(pbig_.get_mpz_t());
        k_ = static_cast<unsigned>(bit_length(pbig_));
        m_ = static_cast<std::uint64_t>((detail::u128(1) << (2 * k_)) / p_);
        lim_ = detail::u128(1) << (2 * k_);
        bar64_ = detail::Barrett64(p_);
        pequeno_ = k_ <= 32;
        {
            const detail::u128 c = static_cast<detail::u128>(p_ - 1) * (p_ - 1);  // >= 1
            max_terminos64_ = (c >> 64) ? 0 : static_cast<std::size_t>(~std::uint64_t(0) / static_cast<std::uint64_t>(c));
        }
        // 2^128 mod p = ((2^128 - 1) mod p + 1) mod p
        c128_ = static_cast<std::uint64_t>((~detail::u128(0)) % p_);
        c128_ = (c128_ + 1 == p_) ? 0 : c128_ + 1;
        // Montgomery (solo para pow, p impar): pinv = -p^(-1) mod 2^64 por
        // iteración de Newton (cada paso duplica los bits correctos: 3 -> 96).
        if (p_ & 1U) {
            std::uint64_t x = p_;  // p·p ≡ 1 (mod 8): correcto en 3 bits
            for (int i = 0; i < 5; ++i) x *= 2 - p_ * x;
            pinv_ = ~x + 1;
        }
    }

    // --- field data ---------------------------------------------------------
    const BigInt& characteristic() const { return pbig_; }
    BigInt cardinality() const { return pbig_; }
    std::uint64_t modulus64() const { return p_; }

    // Igualdad de cuerpos: el mismo p.
    bool operator==(const Fp64& o) const { return p_ == o.p_; }

    // --- constants & coercion ----------------------------------------------
    Elem zero() const { return 0; }
    Elem one() const { return 1; }

    Elem from_int(const BigInt& a) const {
        mpz_srcptr z = a.get_mpz_t();
        if (mpz_sgn(z) >= 0 && mpz_size(z) <= 1) return detail::mpz_to_u64(z) % p_;
        BigInt r;
        mpz_fdiv_r(r.get_mpz_t(), z, pbig_.get_mpz_t());  // en [0, p)
        return detail::mpz_to_u64(r.get_mpz_t());
    }

    // Misma llamada al generador que Fp::random (misma secuencia).
    Elem random() const { return detail::mpz_to_u64(rand_below(pbig_).get_mpz_t()); }

    // --- arithmetic ---------------------------------------------------------
    Elem add(const Elem& a, const Elem& b) const {
        const Elem s = a + b;  // < 2^63: sin desbordamiento
        return s >= p_ ? s - p_ : s;
    }
    Elem sub(const Elem& a, const Elem& b) const { return a >= b ? a - b : a + (p_ - b); }
    Elem neg(const Elem& a) const { return a == 0 ? 0 : p_ - a; }
    Elem mul(const Elem& a, const Elem& b) const { return barrett(static_cast<detail::u128>(a) * b); }

    // Inverso por el algoritmo de Euclides extendido (como invmod/mpz_invert en
    // Fp), en aritmética de 64 bits con signo: |t| <= p < 2^62.
    Elem inv(const Elem& a) const {
        if (a == 0) throw division_by_zero("0 is not invertible in F_p");
        std::int64_t t = 0, nt = 1;
        std::uint64_t r = p_, nr = a;
        while (nr != 0) {
            const std::uint64_t q = r / nr;
            const std::int64_t tt = t - static_cast<std::int64_t>(q) * nt;
            t = nt;
            nt = tt;
            const std::uint64_t rr = r - q * nr;
            r = nr;
            nr = rr;
        }
        return t < 0 ? static_cast<Elem>(t + static_cast<std::int64_t>(p_)) : static_cast<Elem>(t);
    }
    Elem div(const Elem& a, const Elem& b) const { return mul(a, inv(b)); }

    // a^k para cualquier entero k (negativo -> inverso primero; 0^k con k < 0
    // lanza alcp::division_by_zero), por cuadrados y productos sucesivos.
    Elem pow(const Elem& a, const BigInt& k) const {
        if (sgn(k) < 0) {
            if (a == 0) throw division_by_zero("0 is not invertible in F_p");
            const BigInt mk = -k;
            return pow_nonneg(inv(a), mk.get_mpz_t());
        }
        return pow_nonneg(a, k.get_mpz_t());
    }

    // --- acumulación perezosa (ganchos opcionales, ver field.hpp) -----------
    void acc_clear(Acc& s) const { s = 0; }
    void acc_set(Acc& s, const Elem& a) const { s = a; }
    void acc_fma(Acc& s, const Elem& a, const Elem& b) const {
        const Acc t = s + static_cast<Acc>(a) * b;
        // Acarreo (t < s): se ha perdido 2^128; se pliega a un valor < 2p.
        s = (t < s) ? static_cast<Acc>(reduce128(t)) + c128_ : t;
    }
    // Fase 6, p < 2^32 (F_101, cuerpos base de F_(3^30), F_256...): a·b <
    // p^2 < 2^64 se calcula con el producto de 64 bits y se suma sin vigilar
    // el acarreo.  Cota: para desbordar 2^128 harían falta más de 2^64
    // términos en un mismo acumulador, y cada acumulador recibe a lo sumo
    // un término por coeficiente de un polinomio en memoria (escolar,
    // división, composición) o por coeficiente de F_q (FqT: convolución de
    // 2e-1 huecos más e pliegues por hueco); todos son < 2^61.  Mismo valor
    // exacto que acc_fma, luego el mismo resultado tras reducir.  Solo se usa
    // si acc_pequeno(), y la comprobación se hace fuera del bucle
    // (detail::con_acc_fma): dentro de acc_fma frenaba los bucles de p grande.
    bool acc_pequeno() const { return pequeno_; }
    void acc_fma_peq(Acc& s, const Elem& a, const Elem& b) const { s += static_cast<std::uint64_t>(a * b); }
    // Con la comprobación por término: para FqT<Fp64>, cuyos bucles internos
    // son de longitud e (el grado de la extensión) y se llaman por cada
    // producto de F_q; ahí el despacho fuera del bucle medía un 3-8 % peor.
    void acc_fma_rama(Acc& s, const Elem& a, const Elem& b) const {
        if (pequeno_) acc_fma_peq(s, a, b);
        else acc_fma(s, a, b);
    }
    Elem acc_reduce(Acc& s) const { return reduce128(s); }

    // --- acumulador de 64 bits (fase 6, p pequeño) ---------------------------
    // Los bucles escolares de Poly<Fp64> (producto, división y combinación de
    // Brent-Kung, polynomial.hpp) acumulan en una palabra cuando la suma no
    // puede desbordarla: un acumulador que recibe como mucho `terminos`
    // sumandos, cada uno <= (p-1)^2, vale <= terminos·(p-1)^2, y eso es
    // < 2^64 si terminos <= floor((2^64 - 1) / (p-1)^2).  Con p = 101, hasta
    // 1,8·10^15 términos (cualquier polinomio); con p < 2^16, más de 2^32; con
    // p >= 2^32 + 1 nunca.  Un valor de partida < p (división) cuenta como un
    // término más, pues p - 1 <= (p-1)^2.
    bool acc64_cabe(std::size_t terminos) const { return terminos <= max_terminos64_; }
    Elem reduce64(std::uint64_t x) const { return bar64_(x); }

    // --- predicates ---------------------------------------------------------
    bool eq(const Elem& a, const Elem& b) const { return a == b; }
    bool is_zero(const Elem& a) const { return a == 0; }
    bool is_one(const Elem& a) const { return a == 1; }

    // --- conversion ---------------------------------------------------------
    BigInt to_bigint(const Elem& a) const { return detail::u64_to_big(a); }
    std::string to_string(const Elem& a) const { return std::to_string(a); }

    // En F_p la raíz p-ésima de a es a (ver Fp::pth_root).
    Elem pth_root(const Elem& a) const { return a; }

    // --- FFT support --------------------------------------------------------
    unsigned long fft_capacity() const {
        if (p_ == 2) return 0;
        return static_cast<unsigned long>(__builtin_ctzll(p_ - 1));
    }

    // Mismo procedimiento y MISMO consumo del generador que Fp::fft_base_root:
    // x = rand_below(p - 1) + 1 y candidato x^((p-1)/2^cap) hasta que sea
    // primitivo.
    const Elem& fft_base_root() const {
        if (!base_root_ready_) {
            const unsigned long cap = fft_capacity();
            const BigInt cofactor = detail::u64_to_big((p_ - 1) >> cap);
            const BigInt half = (cap == 0) ? BigInt(0) : (BigInt(1) << (cap - 1));
            const BigInt pm1 = pbig_ - 1;
            while (true) {
                const Elem x = detail::mpz_to_u64(rand_below(pm1).get_mpz_t()) + 1;  // x en [1, p-1]
                const Elem cand = pow(x, cofactor);
                if (cap == 0) { base_root_ = cand; break; }
                if (pow(cand, half) != 1) { base_root_ = cand; break; }
            }
            base_root_ready_ = true;
        }
        return base_root_;
    }

    mutable FFTCache<Elem> fft_cache_;

private:
    BigInt pbig_;
    std::uint64_t p_ = 0;
    unsigned k_ = 0;         // bits de p
    std::uint64_t m_ = 0;    // floor(2^(2k) / p)
    detail::u128 lim_ = 0;   // 2^(2k): cota de validez de barrett()
    std::uint64_t c128_ = 0; // 2^128 mod p (= R^2 mod p para Montgomery, R = 2^64)
    std::uint64_t pinv_ = 0; // -p^(-1) mod 2^64 (solo p impar)
    detail::Barrett64 bar64_;  // x mod p para x < 2^64 (fase 6)
    bool pequeno_ = false;     // p < 2^32: a·b cabe en una palabra (acc_fma, fase 6)
    std::size_t max_terminos64_ = 0;  // floor((2^64-1)/(p-1)^2): ver acc64_cabe (fase 6)
    mutable bool base_root_ready_ = false;
    mutable Elem base_root_ = 0;

    // x mod p para x < 2^(2k) (en particular, x = a·b con a, b < p).
    // Los desplazamientos de 128 bits se escriben con dos palabras de 64 bits:
    // con 2 <= k <= 62, todas las cuentas (k-1, 65-k, k+1, 63-k) están en
    // [1, 63], así que no hace falta el caso «desplazamiento >= 64» que g++
    // genera (con cmov) para un u128 >> n de n variable.  Mismo resultado.
    Elem barrett(detail::u128 x) const {
        const std::uint64_t lo = static_cast<std::uint64_t>(x);
        const std::uint64_t hi = static_cast<std::uint64_t>(x >> 64);
        const std::uint64_t x1 = (hi << (65 - k_)) | (lo >> (k_ - 1));  // x >> (k-1), < 2^(k+1)
        const detail::u128 t = static_cast<detail::u128>(x1) * m_;
        const std::uint64_t q = (static_cast<std::uint64_t>(t >> 64) << (63 - k_)) |
                                (static_cast<std::uint64_t>(t) >> (k_ + 1));  // t >> (k+1)
        std::uint64_t r = lo - q * p_;  // en [0, 3p), exacto módulo 2^64
        if (r >= p_) r -= p_;
        if (r >= p_) r -= p_;
        return r;
    }

    // x mod p para cualquier x de 128 bits.  Fase 6: con p pequeño la suma de
    // un acumulador casi nunca cabe en 2^(2k) (con p = 101, en 2^14) pero sí en
    // una palabra; ese caso se reduce con Barrett de 64 bits en lugar de las
    // dos divisiones de antes (hi % p y divq).  Mismo resto exacto.
    Elem reduce128(detail::u128 x) const {
        if (x < lim_) return barrett(x);
        // (Medido: volver aquí a la división de antes, hi % p y divq, hace
        // divmod de grado 2000 en F_(2^61-1) un 15 % más lento que esto.)
        return bar64_.dos_palabras(static_cast<std::uint64_t>(x >> 64), static_cast<std::uint64_t>(x));
    }

    // REDC: t·2^(-64) mod p para t < p·2^64 (aquí t = a·b con a, b < p).
    // t + m·p < 2^124 + 2^126 no desborda; el cociente es < 2p.
    Elem redc(detail::u128 t) const {
        const std::uint64_t m = static_cast<std::uint64_t>(t) * pinv_;
        const std::uint64_t r = static_cast<std::uint64_t>((t + static_cast<detail::u128>(m) * p_) >> 64);
        return r >= p_ ? r - p_ : r;
    }
    Elem mont_mul(Elem a, Elem b) const { return redc(static_cast<detail::u128>(a) * b); }

    // a^e (e >= 0) por cuadrados y productos de izquierda a derecha.  Con p
    // impar, r y a viven en forma de Montgomery dentro del bucle (mismas
    // operaciones, otra representación); con p = 2, Barrett.
    Elem pow_nonneg(Elem a, mpz_srcptr e) const {
        if (mpz_sgn(e) == 0) return 1;
        if ((p_ & 1U) == 0) return pow_bits(a, e, 1, [this](Elem x, Elem y) { return mul(x, y); });
        const Elem am = mont_mul(a, c128_);     // a·R mod p
        const Elem uno = mont_mul(1, c128_);    // R mod p
        const Elem r = pow_bits(am, e, uno, [this](Elem x, Elem y) { return mont_mul(x, y); });
        return redc(r);                         // r·R^(-1) mod p
    }

    // Cadena de cuadrados y productos con el producto mul_ (Barrett o REDC).
    template <class Mul>
    static Elem pow_bits(Elem a, mpz_srcptr e, Elem r, Mul mul_) {
        if (mpz_size(e) == 1) {
            // Exponente de una palabra: los bits se leen de un registro.
            const std::uint64_t x = detail::mpz_to_u64(e);
            for (int i = 63 - __builtin_clzll(x); i >= 0; --i) {
                r = mul_(r, r);
                if ((x >> i) & 1U) r = mul_(r, a);
            }
            return r;
        }
        for (std::size_t i = mpz_sizeinbase(e, 2); i-- > 0; ) {
            r = mul_(r, r);
            if (mpz_tstbit(e, i)) r = mul_(r, a);
        }
        return r;
    }
};

// ¿Son a y b el mismo cuerpo?  El mismo objeto o el mismo p.
inline bool same_field(const Fp64& a, const Fp64& b) { return &a == &b || a == b; }

static_assert(FiniteField<Fp64>);
static_assert(LazyAccumulation<Fp64>);

}  // namespace alcp
