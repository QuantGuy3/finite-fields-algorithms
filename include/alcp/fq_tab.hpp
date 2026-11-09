// alcp/fq_tab.hpp — F_q pequeño (q = p^n <= 2^16) con el elemento en UNA palabra
// y la aritmética por tablas de logaritmos (fase 7, propuesta 7).
//
// Motivo (docs/optimizacion/investigacion-fase6.md, tema 5 «factor»): en F_256
// y cuerpos parecidos casi todo el tiempo de factor() es aritmética de
// elementos de Fq64, donde cada elemento es un std::vector de n coeficientes
// (memoria dinámica en cada operación y una convolución con reducción módulo g
// en cada producto).  Con q <= 2^16 cabe todo el grupo multiplicativo en una
// tabla: el producto y el inverso son dos o tres lecturas de tabla.  El
// prototipo (bench/experimentos/fase6/factor/proto2.cpp), con el medidor
// corregido de la fase 7, pasa factor() de grado 30 en F_256 de 5,4 ms a
// 0,16 ms.
//
// MISMO CUERPO QUE Fq64: FqTab guarda un Fq64 de referencia (ref_) construido
// con los mismos argumentos.  Así la validación (mensajes de error), la
// búsqueda del módulo aleatorio de FqTab(q) / FqTab(p, n) (y su consumo del
// generador), modulus(), from_vector, from_index/to_index y to_string son
// literalmente los de Fq64.  Solo cambia la representación interna:
//
//   * Elem = código del elemento en una palabra: c_0·p^(n-1) + c_1·p^(n-2) +
//     ... + c_{n-1}, es decir, los coeficientes en base p con c_0 como dígito
//     MÁS significativo.  Con esta elección el orden numérico de los códigos
//     es el orden lexicográfico de los std::vector de Fq64 (c_0 primero), y
//     canonical_less (orden de los factores de factor()) ordena igual.  Con
//     p = 2 el código son los n bits de los coeficientes y la suma es XOR.
//   * exp_/log_: potencias de un generador fijo del grupo multiplicativo (el
//     primer elemento primitivo de una lista determinista: t + c y después
//     todos por índice; no usa el generador aleatorio).  mul, inv, pow y
//     pth_root son lecturas de tabla.
//   * p impar: la suma usa logaritmos de Zech, a + b = a·(1 + b/a):
//     zech_[d] = log(1 + g^d), con un centinela si 1 + g^d = 0; neg suma
//     log(-1) = (q-1)/2.  Sin tablas de q^2 entradas.
//
// Mismo consumo del generador aleatorio que Fq64: random() hace las mismas n
// llamadas a Base::random (rand_below(p)) en el mismo orden (c_0 primero), y
// fft_capacity / fft_base_root siguen el procedimiento de FqT con los mismos
// sorteos, así que la raíz base es el mismo elemento.  No es una extensión
// para kronecker.hpp (no ofrece base_field()): Poly<FqTab> multiplica por
// Karatsuba/escolar sobre elementos de una palabra, que con estas tablas ya
// es más rápido que el Kronecker de Fq64 (producto de grado 1000 en F_256:
// 803 -> 302 µs en el prototipo).
#pragma once

#include "bigint.hpp"
#include "errors.hpp"
#include "fft_tables.hpp"
#include "field.hpp"
#include "fp64.hpp"
#include "fq.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace alcp {

class FqTab {
public:
    using Elem = std::uint32_t;
    using Acc = std::uint32_t;  // acumulación: ya es exacta (ver acc_fma)
    using Ref = Fq64;           // el cuerpo de referencia, misma presentación

    // Cardinal máximo: los logaritmos (< q - 1) y los códigos (< q) caben en
    // tablas de 16 bits; con q = 2^16 las tres tablas ocupan ~384 KiB.
    static constexpr unsigned long kMaxCardinality = 1UL << 16;

    // ¿Cabe F_(p^n) en FqTab?  Solo el tamaño (la primalidad de p y la
    // irreducibilidad del módulo las comprueba Fq64 al construir).
    static bool fits(const BigInt& p, unsigned long n) {
        if (sgn(p) <= 0 || n == 0 || n > 16 || p > BigInt(kMaxCardinality)) return false;
        return ipow(p, n) <= BigInt(kMaxCardinality);
    }

    // Los mismos constructores que Fq64 (mismas excepciones y mismo consumo
    // del generador: se construye primero el Fq64).
    explicit FqTab(const BigInt& q) : FqTab(Fq64(q)) {}
    FqTab(const BigInt& p, unsigned long n) : FqTab(Fq64(p, n)) {}
    FqTab(Fp64 base, std::vector<BigInt> g) : FqTab(Fq64(std::move(base), std::move(g))) {}
    explicit FqTab(Fq64 ref) : ref_(std::move(ref)) { construir(); }

    // --- datos del cuerpo (los de Fq64) -------------------------------------
    const Fq64& reference() const { return ref_; }
    const BigInt& characteristic() const { return ref_.characteristic(); }
    BigInt cardinality() const { return ref_.cardinality(); }
    unsigned long extension_degree() const { return ref_.extension_degree(); }
    const std::vector<BigInt>& modulus() const { return ref_.modulus(); }
    bool operator==(const FqTab& o) const { return ref_ == o.ref_; }

    // --- conversión con la representación de Fq64 --------------------------
    Elem encode(const Fq64::Elem& v) const {
        Elem c = 0;
        for (unsigned long i = 0; i < n_; ++i) c = c * static_cast<Elem>(p_) + static_cast<Elem>(v[i]);
        return c;
    }
    Fq64::Elem decode(Elem c) const {
        Fq64::Elem v(n_);
        for (unsigned long i = n_; i-- > 0;) {
            v[i] = c % p_;
            c /= static_cast<Elem>(p_);
        }
        return v;
    }
    // Coeficientes c_0..c_{n-1} (en [0, p)), como los de Fq64.
    std::vector<std::uint64_t> to_coeffs(Elem c) const { return decode(c); }

    // --- constantes y coerción ----------------------------------------------
    Elem zero() const { return 0; }
    Elem one() const { return one_; }
    // (a mod p)·1: el dígito c_0 vale a mod p.
    Elem from_int(const BigInt& a) const {
        return static_cast<Elem>(ref_.base_field().from_int(a)) * one_;
    }
    // n llamadas a Base::random (rand_below(p)), c_0 primero, como FqT::random.
    Elem random() const {
        Elem c = 0;
        for (unsigned long i = 0; i < n_; ++i)
            c = c * static_cast<Elem>(p_) + static_cast<Elem>(ref_.base_field().random());
        return c;
    }
    Elem from_vector(std::vector<BigInt> v) const { return encode(ref_.from_vector(std::move(v))); }
    // Índice en base p (c_0 + c_1·p + ...), el de Fq64; tablas de q entradas.
    Elem from_index(const BigInt& idx) const {
        BigInt t;
        mpz_fdiv_r_ui(t.get_mpz_t(), idx.get_mpz_t(), q_);
        return idx2code_[mpz_get_ui(t.get_mpz_t())];
    }
    BigInt to_index(Elem a) const { return BigInt(static_cast<unsigned long>(code2idx_[a])); }

    // --- aritmética ---------------------------------------------------------
    Elem add(Elem a, Elem b) const {
        if (car2_) return a ^ b;
        if (a == 0) return b;
        if (b == 0) return a;
        // a + b = a·(1 + b/a) = g^(la + zech(lb - la)).
        const std::uint32_t la = log_[a], lb = log_[b];
        const std::uint32_t d = lb >= la ? lb - la : lb + nq_ - la;
        const std::uint32_t z = zech_[d];
        if (z == kSinLog) return 0;
        return exp_[la + z];
    }
    Elem neg(Elem a) const {
        if (car2_ || a == 0) return a;
        return exp_[log_[a] + mitad_];  // -1 = g^((q-1)/2)
    }
    Elem sub(Elem a, Elem b) const { return car2_ ? (a ^ b) : add(a, neg(b)); }
    Elem mul(Elem a, Elem b) const {
        if (a == 0 || b == 0) return 0;
        return exp_[static_cast<std::uint32_t>(log_[a]) + log_[b]];
    }
    // Mismo mensaje y misma excepción que FqT::inv.
    Elem inv(Elem a) const {
        if (a == 0) throw division_by_zero("0 is not invertible in F_q");
        return exp_[nq_ - log_[a]];
    }
    Elem div(Elem a, Elem b) const { return mul(a, inv(b)); }
    // Como FqT::pow: con k < 0 se invierte primero (0 lanza); 0^0 = 1.
    Elem pow(Elem a, const BigInt& k) const {
        if (a == 0) {
            if (sgn(k) < 0) throw division_by_zero("0 is not invertible in F_q");
            return sgn(k) == 0 ? one_ : 0;
        }
        // a^k = a^(k mod (q-1)) para a != 0 (mpz_fdiv_ui da el resto en [0, q-1)).
        const unsigned long r = mpz_fdiv_ui(k.get_mpz_t(), nq_);
        return exp_[static_cast<std::uint32_t>((static_cast<std::uint64_t>(log_[a]) * r) % nq_)];
    }
    // Raíz p-ésima: a^(q/p).
    Elem pth_root(Elem a) const {
        if (a == 0) return 0;
        return exp_[static_cast<std::uint32_t>((static_cast<std::uint64_t>(log_[a]) * qp_) % nq_)];
    }

    // --- acumulación (ganchos de field.hpp) ----------------------------------
    // La suma por tablas ya es exacta y barata: no hay nada que diferir.
    void acc_clear(Acc& s) const { s = 0; }
    void acc_set(Acc& s, Elem a) const { s = a; }
    void acc_fma(Acc& s, Elem a, Elem b) const { s = add(s, mul(a, b)); }
    Elem acc_reduce(Acc& s) const { return s; }

    // --- predicados ---------------------------------------------------------
    bool eq(Elem a, Elem b) const { return a == b; }
    bool is_zero(Elem a) const { return a == 0; }
    bool is_one(Elem a) const { return a == one_; }

    std::string to_string(Elem a) const { return ref_.to_string(decode(a)); }

    // --- FFT: el mismo procedimiento que FqT (mismos sorteos, misma raíz) ---
    unsigned long fft_capacity() const { return ref_.fft_capacity(); }
    const Elem& fft_base_root() const {
        if (!base_root_ready_) {
            const unsigned long cap = fft_capacity();
            const BigInt order = BigInt(1) << cap;
            const BigInt cofactor = (cardinality() - 1) / order;
            const BigInt half = (cap == 0) ? BigInt(0) : (BigInt(1) << (cap - 1));
            while (true) {
                const Elem x = random();
                if (is_zero(x)) continue;
                const Elem cand = pow(x, cofactor);
                if (cap == 0) { base_root_ = cand; break; }
                if (!is_one(pow(cand, half))) { base_root_ = cand; break; }
            }
            base_root_ready_ = true;
        }
        return base_root_;
    }
    mutable FFTCache<Elem> fft_cache_;

private:
    static constexpr std::uint16_t kSinLog = 0xFFFF;  // log(0) en zech_ (los logs son < q - 1 <= 65535)

    Fq64 ref_;
    std::uint64_t p_ = 2;
    unsigned long n_ = 1;
    std::uint32_t q_ = 2, nq_ = 1, mitad_ = 0, qp_ = 1;  // q, q - 1, (q-1)/2, (q/p) mod (q-1)
    Elem one_ = 1;
    bool car2_ = true;
    std::vector<std::uint16_t> exp_;   // exp_[k] = código de g^k, k en [0, 2(q-1))
    std::vector<std::uint16_t> log_;   // log_[c] = log_g del código c != 0
    std::vector<std::uint16_t> zech_;  // p impar: zech_[d] = log(1 + g^d) o kSinLog
    std::vector<std::uint16_t> idx2code_, code2idx_;
    mutable bool base_root_ready_ = false;
    mutable Elem base_root_ = 0;

    void construir() {
        const BigInt& p = ref_.characteristic();
        n_ = ref_.extension_degree();
        if (!fits(p, n_))
            throw std::invalid_argument("FqTab: q = " + ref_.cardinality().get_str() +
                                        " demasiado grande (máximo 2^16; use Fq64)");
        p_ = mpz_get_ui(p.get_mpz_t());
        q_ = static_cast<std::uint32_t>(mpz_get_ui(ref_.cardinality().get_mpz_t()));
        nq_ = q_ - 1;
        car2_ = p_ == 2;
        mitad_ = nq_ / 2;
        one_ = 1;
        for (unsigned long i = 1; i < n_; ++i) one_ *= static_cast<Elem>(p_);
        qp_ = (q_ / static_cast<std::uint32_t>(p_)) % nq_;

        // Índice en base p (c_0 menos significativo) <-> código (c_0 más significativo).
        idx2code_.assign(q_, 0);
        code2idx_.assign(q_, 0);
        for (std::uint32_t idx = 0; idx < q_; ++idx) {
            std::uint32_t t = idx, c = 0;
            for (unsigned long i = 0; i < n_; ++i) {
                c = c * static_cast<std::uint32_t>(p_) + t % static_cast<std::uint32_t>(p_);
                t /= static_cast<std::uint32_t>(p_);
            }
            idx2code_[idx] = static_cast<std::uint16_t>(c);
            code2idx_[c] = static_cast<std::uint16_t>(idx);
        }

        // -g_0..-g_{n-1} en [0, p): x^n = Σ -g_i x^i.
        std::vector<std::uint64_t> ng(n_);
        for (unsigned long i = 0; i < n_; ++i) {
            const std::uint64_t gi = mpz_get_ui(ref_.modulus()[i].get_mpz_t());
            ng[i] = gi == 0 ? 0 : p_ - gi;
        }
        const std::uint64_t P = p_;
        // y = x·w en F_p[t]/(g) con coeficientes pequeños (solo al construir).
        auto mul_dig = [&](const std::vector<std::uint64_t>& x, const std::vector<std::uint64_t>& w) {
            std::vector<std::uint64_t> s(2 * n_ - 1, 0);
            for (unsigned long i = 0; i < n_; ++i)
                for (unsigned long j = 0; j < n_; ++j) s[i + j] = (s[i + j] + x[i] * w[j]) % P;
            for (unsigned long k = 2 * n_ - 1; k-- > n_;) {
                const std::uint64_t c = s[k];
                if (c == 0) continue;
                for (unsigned long t = 0; t < n_; ++t) s[k - n_ + t] = (s[k - n_ + t] + c * ng[t]) % P;
            }
            s.resize(n_);
            return s;
        };
        // Divisores primos de q - 1 (q - 1 < 2^16: división por tentativa).
        std::vector<std::uint32_t> primos;
        for (std::uint32_t m = nq_, d = 2; m > 1; ++d) {
            if (d * d > m) { primos.push_back(m); break; }
            if (m % d == 0) {
                primos.push_back(d);
                while (m % d == 0) m /= d;
            }
        }
        auto potencia = [&](std::vector<std::uint64_t> b, std::uint32_t k) {
            std::vector<std::uint64_t> r(n_, 0);
            r[0] = 1;
            while (k) {
                if (k & 1U) r = mul_dig(r, b);
                k >>= 1;
                if (k) b = mul_dig(b, b);
            }
            return r;
        };
        auto es_uno = [&](const std::vector<std::uint64_t>& v) {
            if (v[0] != 1) return false;
            for (unsigned long i = 1; i < n_; ++i)
                if (v[i] != 0) return false;
            return true;
        };
        auto primitivo = [&](const std::vector<std::uint64_t>& w) {
            bool nulo = true;
            for (auto c : w) nulo = nulo && c == 0;
            if (nulo) return false;
            for (std::uint32_t r : primos)
                if (es_uno(potencia(w, nq_ / r))) return false;
            return true;
        };
        // Generador: el primer elemento primitivo de la lista t + c (c = 0..p-1)
        // y después todos por índice 1, 2, ...  Es determinista y no usa el
        // generador aleatorio; cualquier generador da los mismos resultados.
        // Con gen = t + c, multiplicar por gen cuesta O(n) (desplazar y
        // reducir una vez), lo que abarata la tabla de q = 2^16.
        auto digitos = [&](std::uint32_t idx) {
            std::vector<std::uint64_t> w(n_);
            for (unsigned long i = 0; i < n_; ++i) {
                w[i] = idx % static_cast<std::uint32_t>(P);
                idx /= static_cast<std::uint32_t>(P);
            }
            return w;
        };
        std::vector<std::uint64_t> gen;
        long cte = -1;  // c si gen = t + c
        if (n_ >= 2)
            for (std::uint32_t c = 0; c < P && gen.empty(); ++c)
                if (primitivo(digitos(static_cast<std::uint32_t>(P) + c))) {
                    gen = digitos(static_cast<std::uint32_t>(P) + c);
                    cte = static_cast<long>(c);
                }
        for (std::uint32_t idx = 1; idx < q_ && gen.empty(); ++idx)
            if (primitivo(digitos(idx))) gen = digitos(idx);
        auto por_gen = [&](const std::vector<std::uint64_t>& x) {
            if (cte < 0) return mul_dig(x, gen);
            // x·(t + c) = x·t + c·x, con x·t = desplazar y sumar top·(-g).
            std::vector<std::uint64_t> y(n_);
            const std::uint64_t top = x[n_ - 1];
            for (unsigned long i = 0; i < n_; ++i)
                y[i] = ((i ? x[i - 1] : 0) + top * ng[i] + static_cast<std::uint64_t>(cte) * x[i]) % P;
            return y;
        };
        auto codigo = [&](const std::vector<std::uint64_t>& v) {
            std::uint32_t c = 0;
            for (unsigned long i = 0; i < n_; ++i) c = c * static_cast<std::uint32_t>(P) + static_cast<std::uint32_t>(v[i]);
            return c;
        };
        exp_.assign(2 * static_cast<std::size_t>(nq_), 0);
        log_.assign(q_, 0);
        std::vector<std::uint64_t> x(n_, 0);
        x[0] = 1;
        for (std::uint32_t k = 0; k < nq_; ++k) {
            const std::uint32_t c = codigo(x);
            exp_[k] = exp_[k + nq_] = static_cast<std::uint16_t>(c);
            log_[c] = static_cast<std::uint16_t>(k);
            x = por_gen(x);
        }
        if (!car2_) {
            // zech_[d] = log(1 + g^d); 1 + g^d = 0 exactamente cuando g^d = -1.
            zech_.assign(nq_, kSinLog);
            for (std::uint32_t d = 0; d < nq_; ++d) {
                const std::uint32_t c = exp_[d];
                // Sumar 1 al dígito c_0 (el más significativo del código).
                const std::uint32_t c0 = c / one_;
                const std::uint32_t s = c - c0 * one_ + ((c0 + 1) % static_cast<std::uint32_t>(P)) * one_;
                if (s != 0) zech_[d] = log_[s];
            }
        }
    }
};

// ¿Son a y b el mismo cuerpo?  Mismos p y módulo mónico (como Fq64).
inline bool same_field(const FqTab& a, const FqTab& b) { return &a == &b || a == b; }

static_assert(FiniteField<FqTab>);
static_assert(LazyAccumulation<FqTab>);

}  // namespace alcp
