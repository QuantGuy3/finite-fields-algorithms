// alcp/polynomial.hpp — Generic univariate polynomials over a finite field.
//
// A single template Poly<Field> carries a pointer to its field K and a
// little-endian, trimmed coefficient vector (c[i] is the coefficient of x^i;
// the zero polynomial has an empty vector, degree -1).  Because every
// operation goes through the Field interface, this one file provides the whole
// algebra of F_p[x] AND F_q[x] (and any other field implementing the
// interface), fulfilling the "algoritmos generales" requirement:
//
//   * ring arithmetic (+, -, *, divmod, /, %, scalar mul, monic, derivative)
//   * multiplication that AUTO-DISPATCHES schoolbook / Karatsuba / FFT by size
//   * gcd, extended gcd, modular inverse, modular exponentiation
//   * test de irreducibilidad de Ben-Or (gcd(f, x^{q^i} - x) para i <= n/2)
//   * Cantor-Zassenhaus factorization
//        radical (square-free part) -> distinct-degree -> equal-degree,
//        with per-factor multiplicities  (port of factorizacion.py, unified)
//   * evaluación (Horner), composición f(g(x)), test libre de cuadrados
//
// Cuerpo común: todas las operaciones binarias (+, -, *, mul, divmod, /, %,
// gcd, gcd_ext, inv_mod, pow_mod, multiplicity, compose) exigen que ambos
// operandos estén sobre el mismo cuerpo: el mismo objeto (camino rápido, una
// comparación de punteros) o un cuerpo con los mismos parámetros (same_field).
// Si no, lanzan std::invalid_argument("polinomios sobre cuerpos distintos").
// operator== NO lanza: dos polinomios sobre cuerpos distintos son distintos.
//
// Vida del cuerpo: Poly guarda un puntero al cuerpo, que debe sobrevivir al
// polinomio.  Por eso los constructores y factorías que reciben un cuerpo
// temporal (rvalue) están borrados (= delete): Poly<Fp>(Fp(7)) no compila.
//
// Errores: dividir por el polinomio cero lanza alcp::division_by_zero (deriva
// de std::domain_error); los argumentos inválidos, std::invalid_argument.
#pragma once

#include "bigint.hpp"
#include "errors.hpp"
#include "fft.hpp"
#include "field.hpp"
#include "kronecker.hpp"

#include <algorithm>
#include <cstddef>
#include <future>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace alcp {

// Which multiplication algorithm to use.  Auto picks by operand size; the
// explicit modes exist so tests can cross-check them against each other.
enum class MulAlgo { Auto, Schoolbook, Karatsuba, FFT, Kronecker };

// Cómo calcular la sucesión de Frobenius H_i = x^{q^i} mod f en is_irreducible
// y distinct_degree (fase 4, nivel B): Potencia = H_i = H_{i-1}^q (el original);
// Composicion = H_i = H_{i-1}(H_1) por composición modular de Brent-Kung con
// las potencias de H_1 precalculadas; Auto elige (ver Frobenius).
enum class FrobeniusAlgo { Auto, Potencia, Composicion };

template <class Field>
class Poly {
public:
    using Elem = typename Field::Elem;

    // Below this many coefficients in the smaller operand, schoolbook wins.
    static constexpr std::size_t kKaratsubaCutoff = 32;
    // Karatsuba reparte sus tres subproductos entre hilos solo si ambos
    // operandos tienen al menos tantos coeficientes (coste de crear hilos):
    // 256 si el elemento es un entero GMP (F_p), 1024 si es una palabra de
    // máquina (Fp64: el producto cuesta unos nanosegundos y con 64-512
    // coeficientes crear los hilos cuesta más que el propio producto) y 64 si
    // es un vector (F_q: cada producto de elementos ya cuesta O(n^2)).
    static constexpr std::size_t kParKaratsubaMin =
        std::is_same_v<Elem, BigInt> ? 256 : (std::is_arithmetic_v<Elem> ? 1024 : 64);

    // Número máximo de intentos aleatorios de equal_degree_split antes de
    // concluir que la entrada no cumple la precondición (ver equal_degree).
    // Con una entrada válida cada intento separa factores con probabilidad del
    // orden de 1/2 (análisis estándar de Cantor-Zassenhaus), de modo que fallar
    // 256 veces seguidas es prácticamente imposible; con una entrada inválida
    // (p. ej. un irreducible de grado mayor que r) evita un bucle infinito.
    // Mientras no se agota, la secuencia de intentos es idéntica a la original.
    static constexpr int kMaxSplitAttempts = 256;

    // ---- construction ------------------------------------------------------
    explicit Poly(const Field& field) : K_(&field) {}  // zero polynomial
    Poly(const Field& field, std::vector<Elem> coeffs) : K_(&field), c_(std::move(coeffs)) { trim(); }

    // Un cuerpo temporal dejaría el puntero K_ colgando: prohibido.
    explicit Poly(const Field&&) = delete;
    Poly(const Field&&, std::vector<Elem>) = delete;

    static Poly zero(const Field& F) { return Poly(F); }
    static Poly one(const Field& F) { return Poly(F, {F.one()}); }
    static Poly constant(const Field& F, const Elem& v) { return Poly(F, {v}); }
    static Poly variable(const Field& F) { return Poly(F, {F.zero(), F.one()}); }  // x

    // c * x^deg
    static Poly monomial(const Field& F, std::size_t deg, const Elem& coef) {
        std::vector<Elem> v(deg + 1, F.zero());
        v[deg] = coef;
        return Poly(F, std::move(v));
    }
    // Build from a list of integers (handy for F_p[x]); each is coerced.
    static Poly from_ints(const Field& F, const std::vector<BigInt>& ints) {
        std::vector<Elem> v;
        v.reserve(ints.size());
        for (const auto& a : ints) v.push_back(F.from_int(a));
        return Poly(F, std::move(v));
    }
    // Uniformly random polynomial of exact degree d (leading coeff nonzero).
    static Poly random(const Field& F, int d) {
        if (d < 0) return zero(F);
        std::vector<Elem> v(d + 1);
        for (int i = 0; i < d; ++i) v[i] = F.random();
        do { v[d] = F.random(); } while (F.is_zero(v[d]));
        return Poly(F, std::move(v));
    }

    // Las mismas factorías con un cuerpo temporal: prohibidas (puntero colgante).
    static Poly zero(const Field&&) = delete;
    static Poly one(const Field&&) = delete;
    static Poly constant(const Field&&, const Elem&) = delete;
    static Poly variable(const Field&&) = delete;
    static Poly monomial(const Field&&, std::size_t, const Elem&) = delete;
    static Poly from_ints(const Field&&, const std::vector<BigInt>&) = delete;
    static Poly random(const Field&&, int) = delete;

    // ---- basics ------------------------------------------------------------
    const Field& field() const { return *K_; }
    const std::vector<Elem>& coeffs() const { return c_; }
    int degree() const { return static_cast<int>(c_.size()) - 1; }
    bool is_zero() const { return c_.empty(); }
    bool is_one() const { return c_.size() == 1 && K_->is_one(c_[0]); }
    // Coeficiente principal.  El polinomio cero no tiene: lanza
    // std::domain_error (antes era comportamiento indefinido).
    const Elem& lead() const {
        if (c_.empty()) throw std::domain_error("lead: el polinomio cero no tiene coeficiente principal");
        return c_.back();
    }

    // Coeficiente de x^i (el cero del cuerpo si i > grado).
    Elem coeff(std::size_t i) const { return i < c_.size() ? c_[i] : K_->zero(); }

    // Igualdad de polinomios.  Si los cuerpos son distintos (ni el mismo objeto
    // ni same_field) devuelve false en lugar de lanzar: son elementos de anillos
    // distintos y, por tanto, distintos.
    bool operator==(const Poly& o) const {
        if (K_ != o.K_ && !detail::fields_match(*K_, *o.K_)) return false;
        if (c_.size() != o.c_.size()) return false;
        for (std::size_t i = 0; i < c_.size(); ++i)
            if (!K_->eq(c_[i], o.c_[i])) return false;
        return true;
    }
    bool operator!=(const Poly& o) const { return !(*this == o); }

    // ---- additive structure ------------------------------------------------
    // Suma y resta: el resultado se construye una sola vez, coeficiente a
    // coeficiente (sin rellenar primero con ceros y copiar), con las mismas
    // operaciones del cuerpo: add/sub en la parte común; en la parte que solo
    // tiene un operando, su copia (o -o_i = 0 - o_i en la resta).
    Poly operator+(const Poly& o) const {
        require_same_field(o);
        const std::size_t m = std::min(c_.size(), o.c_.size());
        const std::vector<Elem>& larga = c_.size() >= o.c_.size() ? c_ : o.c_;
        std::vector<Elem> r;
        r.reserve(larga.size());
        for (std::size_t i = 0; i < m; ++i) r.push_back(K_->add(c_[i], o.c_[i]));
        r.insert(r.end(), larga.begin() + static_cast<std::ptrdiff_t>(m), larga.end());
        return Poly(*K_, std::move(r));
    }
    Poly operator-() const {
        std::vector<Elem> r(c_.size());
        for (std::size_t i = 0; i < c_.size(); ++i) r[i] = K_->neg(c_[i]);
        return Poly(*K_, std::move(r));
    }
    Poly operator-(const Poly& o) const {
        require_same_field(o);
        const std::size_t m = std::min(c_.size(), o.c_.size());
        std::vector<Elem> r;
        r.reserve(std::max(c_.size(), o.c_.size()));
        for (std::size_t i = 0; i < m; ++i) r.push_back(K_->sub(c_[i], o.c_[i]));
        for (std::size_t i = m; i < c_.size(); ++i) r.push_back(c_[i]);
        for (std::size_t i = m; i < o.c_.size(); ++i) r.push_back(K_->neg(o.c_[i]));
        return Poly(*K_, std::move(r));
    }

    // ---- multiplication ----------------------------------------------------
    Poly operator*(const Poly& o) const { return mul(o, MulAlgo::Auto); }

    Poly mul(const Poly& o, MulAlgo algo) const {
        require_same_field(o);
        if (is_zero() || o.is_zero()) return zero(*K_);
        const std::size_t na = c_.size(), nb = o.c_.size();
        if (algo == MulAlgo::Auto) {
            const std::size_t L = na + nb - 1;
            const std::size_t small = std::min(na, nb);
            const bool can_fft = fft::fft_exponent(*K_, L) >= 0;
            // Regla original: FFT si cabe y L >= 64.  Se conserva solo para
            // reproducir el consumo del generador aleatorio (ver abajo).
            const bool old_fft = can_fft && L >= 64;
            bool use_fft = can_fft && prefer_fft(*K_, small, std::max(na, nb), L);
            // La primera transformada de un cuerpo calcula su raíz base con el
            // RNG.  Para que la secuencia aleatoria (y por tanto factor(),
            // Fq(p, n), etc.) sea idéntica a la original, la raíz se calcula
            // exactamente cuando la regla original habría usado la FFT: no se
            // adelanta (FFT nueva solo si la caché ya tiene tablas, es decir,
            // si la raíz ya existe) ni se omite (se pide aunque no se use).
            if (use_fft && !old_fft && K_->fft_cache_.empty()) use_fft = false;
            // (Se construyen las tablas, que piden la raíz: así la caché deja
            // de estar vacía y los productos siguientes pueden usar la FFT.)
            if (old_fft && !use_fft) (void)fft::tables(*K_, fft::ceil_log2(L));
            // Fase 4 (nivel B): Kronecker sobre Z cuando es más rápido.  Va después
            // de la línea anterior, así que el consumo del generador no cambia.
            if (prefer_kronecker(*K_, small, can_fft))
                algo = MulAlgo::Kronecker;
            else if (use_fft)
                algo = MulAlgo::FFT;
            else if (small >= kKaratsubaCutoff)
                algo = MulAlgo::Karatsuba;
            else
                algo = MulAlgo::Schoolbook;
        }
        switch (algo) {
            case MulAlgo::FFT: {
                if (fft::fft_exponent(*K_, na + nb - 1) < 0)
                    return mul_karatsuba(o);  // field cannot host this size
                return Poly(*K_, fft::convolve(*K_, c_, o.c_));
            }
            case MulAlgo::Karatsuba: return mul_karatsuba(o);
            case MulAlgo::Kronecker: return mul_kronecker(o);
            case MulAlgo::Schoolbook:
            default: return mul_school(o);
        }
    }

    // Producto por sustitución de Kronecker sobre Z con GMP (kronecker.hpp,
    // fase 4, nivel B): los coeficientes se empaquetan en un entero, se
    // multiplica con mpn_mul y se desempaqueta reduciendo módulo p (y módulo g
    // en F_q).  No usa raíces de la unidad ni el generador aleatorio.
    Poly mul_kronecker(const Poly& o) const {
        require_same_field(o);
        if (is_zero() || o.is_zero()) return zero(*K_);
        if constexpr (ks::Extension<Field>)
            return Poly(*K_, ks::mul_extension(*K_, c_.data(), c_.size(), o.c_.data(), o.c_.size()));
        else
            return Poly(*K_, ks::mul_primo(*K_, c_.data(), c_.size(), o.c_.data(), o.c_.size()));
    }

    Poly mul_school(const Poly& o) const {
        require_same_field(o);
        if (is_zero() || o.is_zero()) return zero(*K_);
        if constexpr (LazyAccumulation<Field>) {
            // Mismo doble bucle, con reducción perezosa: cada coeficiente del
            // producto acumula sus productos c_i·o_j sin reducir y se reduce
            // una sola vez al final (ganchos opcionales del cuerpo, field.hpp).
            const std::size_t L = c_.size() + o.c_.size() - 1;
            std::vector<typename Field::Acc> acc(L);
            for (auto& s : acc) K_->acc_clear(s);
            for (std::size_t i = 0; i < c_.size(); ++i) {
                if (K_->is_zero(c_[i])) continue;
                for (std::size_t j = 0; j < o.c_.size(); ++j)
                    K_->acc_fma(acc[i + j], c_[i], o.c_[j]);
            }
            std::vector<Elem> r;
            r.reserve(L);
            for (auto& s : acc) r.push_back(K_->acc_reduce(s));
            return Poly(*K_, std::move(r));
        }
        std::vector<Elem> r(c_.size() + o.c_.size() - 1, K_->zero());
        for (std::size_t i = 0; i < c_.size(); ++i) {
            if (K_->is_zero(c_[i])) continue;
            for (std::size_t j = 0; j < o.c_.size(); ++j)
                r[i + j] = K_->add(r[i + j], K_->mul(c_[i], o.c_[j]));
        }
        return Poly(*K_, std::move(r));
    }

    // Recursive Karatsuba (port of fp_x_mult_karatsuba / fq_x_mult_karatsuba).
    // Implementación sobre vistas (puntero + longitud) en karatsuba_rec: la
    // recursión y los subproductos son exactamente los del original
    // (aL·bL, aH·bH, (aL+aH)·(bL+bH), con las mismas mitades recortadas), pero
    // sin copiar las mitades a polinomios nuevos (slice) ni desplazar copiando
    // (shift): z1 se corrige en sitio y la recombinación acumula en el vector
    // resultado.
    Poly mul_karatsuba(const Poly& o) const {
        require_same_field(o);
        if (is_zero() || o.is_zero()) return zero(*K_);
        if (c_.size() <= kKaratsubaCutoff || o.c_.size() <= kKaratsubaCutoff)
            return mul_school(o);
        // Profundidad paralela: 2 niveles (hasta 9 ramas) si hay >= 9 hilos,
        // 1 nivel (3 ramas) si hay >= 3, ninguno en otro caso.
        const unsigned T = par::max_threads();
        const int depth = (T >= 9) ? 2 : (T >= 3 ? 1 : 0);
        return Poly(*K_, karatsuba_rec(*K_, c_.data(), c_.size(), o.c_.data(), o.c_.size(), depth));
    }

    // ---- Euclidean division ------------------------------------------------
    // Returns (q, r) with *this = q*o + r and (r == 0 or deg r < deg o).
    // Dividir por el polinomio cero lanza alcp::division_by_zero.
    // Despacho (fase 4): con divisor y cociente grandes (prefer_newton) se usa
    // la división por inversión de Newton (divmod_newton); si no, la división
    // escolar de siempre (divmod_school).  El cociente y el resto son únicos,
    // así que el resultado es idéntico por los dos caminos.
    std::pair<Poly, Poly> divmod(const Poly& o) const {
        require_same_field(o);
        if (o.is_zero()) throw division_by_zero("division by zero polynomial");
        if (degree() < o.degree()) return {zero(*K_), *this};
        if (prefer_newton(*K_, static_cast<std::size_t>(o.degree()),
                          static_cast<std::size_t>(degree() - o.degree()) + 1))
            return divmod_newton(o);
        return divmod_school(o);
    }

    // División escolar (el algoritmo original, port de divmod de cuerpos_finitos.py).
    std::pair<Poly, Poly> divmod_school(const Poly& o) const {
        require_same_field(o);
        if (o.is_zero()) throw division_by_zero("division by zero polynomial");
        if (degree() < o.degree()) return {zero(*K_), *this};

        const int db = o.degree();
        Elem inv_lead = K_->inv(o.lead());
        if constexpr (LazyAccumulation<Field>) {
            // Misma división escolar con reducción perezosa: el resto de trabajo
            // se guarda en acumuladores sin reducir; cada paso resta
            // factor·o (sumando (-factor)·o_j, que es lo mismo módulo p) y solo
            // se reduce el coeficiente principal cuando hace falta leerlo.  El
            // término j = db no se acumula: el algoritmo lo anula exactamente
            // y ese coeficiente ya no se vuelve a leer.
            const std::size_t nr = c_.size();
            std::vector<typename Field::Acc> acc(nr);
            for (std::size_t k = 0; k < nr; ++k) K_->acc_set(acc[k], c_[k]);
            std::vector<Elem> q(degree() - db + 1, K_->zero());
            for (int dr = static_cast<int>(nr) - 1; dr >= db; --dr) {
                Elem lr = K_->acc_reduce(acc[dr]);
                if (K_->is_zero(lr)) continue;  // cero que ha aparecido: se salta
                const int pos = dr - db;
                Elem factor = K_->mul(lr, inv_lead);
                const Elem nfactor = K_->neg(factor);
                q[pos] = std::move(factor);
                for (int j = 0; j < db; ++j) K_->acc_fma(acc[pos + j], nfactor, o.c_[j]);
            }
            std::vector<Elem> rem;
            rem.reserve(static_cast<std::size_t>(db));
            for (int k = 0; k < db; ++k) rem.push_back(K_->acc_reduce(acc[k]));
            return {Poly(*K_, std::move(q)), Poly(*K_, std::move(rem))};
        }
        std::vector<Elem> r = c_;  // working remainder
        std::vector<Elem> q(degree() - db + 1, K_->zero());

        for (int dr = static_cast<int>(r.size()) - 1; dr >= db; ) {
            // strip any leading zeros that appeared
            while (dr >= db && K_->is_zero(r[dr])) --dr;
            if (dr < db) break;
            int pos = dr - db;
            Elem factor = K_->mul(r[dr], inv_lead);
            q[pos] = factor;
            for (int j = 0; j <= db; ++j)
                r[pos + j] = K_->sub(r[pos + j], K_->mul(factor, o.c_[j]));
            --dr;
        }
        return {Poly(*K_, std::move(q)), Poly(*K_, std::move(r))};
    }

    // División por inversión de Newton (fase 4), con las piezas de la
    // Práctica 3 (algoritmos_rapidos.py): con m = deg f, n = deg g y
    // d = m - n, los coeficientes invertidos cumplen
    //     rev(q) = rev(f) · rev(g)^{-1}   (mod x^{d+1}),
    // donde rev(g)^{-1} es el inverso de una serie de potencias, que se calcula
    // con la iteración de Newton de fp_toep_inf_inv (inverse_series).  Después,
    // r = f - g·q.  Coste O(M(d) + M(n, d)) en lugar de O(n·d).
    std::pair<Poly, Poly> divmod_newton(const Poly& o) const {
        require_same_field(o);
        if (o.is_zero()) throw division_by_zero("division by zero polynomial");
        if (degree() < o.degree()) return {zero(*K_), *this};
        const std::size_t d1 = static_cast<std::size_t>(degree() - o.degree()) + 1;
        return divmod_with_inverse(o, inverse_series(*K_, reversed_head(o, d1), d1));
    }

    // Inverso de la serie de potencias a módulo x^n (a[0] != 0): la iteración
    // de Newton de fp_toep_inf_inv / fq_toep_inf_inv (Práctica 3) con la misma
    // recursión, k = ceil(n/2) y x_n = x_k·(2 - a·x_k) mod x^n.  Solo se usan
    // los n primeros coeficientes de a (los demás no influyen módulo x^n).
    // Los productos internos no adelantan el cálculo de la raíz de la FFT
    // (mul_no_rng), para conservar la secuencia del generador aleatorio.
    static std::vector<Elem> inverse_series(const Field& K, const std::vector<Elem>& a, std::size_t n) {
        if (n == 0) return {};
        if (a.empty() || K.is_zero(a[0]))
            throw division_by_zero("inverse_series: a[0] es cero (serie no invertible)");
        if (n == 1) return {K.inv(a[0])};
        const std::size_t k = (n + 1) / 2;
        std::vector<Elem> xk = inverse_series(K, a, k);
        Poly X(K, xk);
        Poly A(K, std::vector<Elem>(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(std::min(n, a.size()))));
        Poly T(K, head(K, A.mul_no_rng(X), n));
        Poly R = constant(K, K.add(K.one(), K.one())) - T;   // 2 - T
        return head(K, R.mul_no_rng(X), n);
    }

    // División de *this entre o dado ginv = rev(o)^{-1} mod x^{d1}, con
    // d1 = deg(*this) - deg(o) + 1 (deg(*this) >= deg(o)): sirve para dividir
    // muchas veces entre el mismo o calculando el inverso una sola vez (lo
    // que hace Reducer).  Precondición: ginv.size() == d1 exactamente.
    std::pair<Poly, Poly> divmod_with_inverse(const Poly& o, const std::vector<Elem>& ginv) const {
        const Field& K = *K_;
        const std::size_t m1 = c_.size(), d1 = ginv.size();
        std::vector<Elem> frev;
        frev.reserve(d1);
        for (std::size_t i = 0; i < d1; ++i) frev.push_back(c_[m1 - 1 - i]);
        const Poly F(K, std::move(frev)), G(K, ginv);
        const Poly FG = F.mul_no_rng(G);
        std::vector<Elem> qrev = head(K, FG, d1);
        std::reverse(qrev.begin(), qrev.end());
        Poly q(K, std::move(qrev));
        Poly r = *this - o.mul_no_rng(q);
        return {std::move(q), std::move(r)};
    }
    Poly operator/(const Poly& o) const { return divmod(o).first; }
    Poly operator%(const Poly& o) const { return divmod(o).second; }

    // ---- scalar / shape ----------------------------------------------------
    Poly scalar_mul(const Elem& e) const {
        std::vector<Elem> r(c_.size());
        for (std::size_t i = 0; i < c_.size(); ++i) r[i] = K_->mul(c_[i], e);
        return Poly(*K_, std::move(r));
    }
    Poly monic() const {
        if (is_zero()) return *this;
        return scalar_mul(K_->inv(lead()));
    }
    Poly derivative() const {
        if (c_.size() <= 1) return zero(*K_);
        std::vector<Elem> r(c_.size() - 1);
        for (std::size_t i = 1; i < c_.size(); ++i)
            r[i - 1] = K_->mul(c_[i], K_->from_int(BigInt(static_cast<unsigned long>(i))));
        return Poly(*K_, std::move(r));
    }

    // Plain (non-modular) power.
    Poly pow(const BigInt& k) const {
        if (k < 0) throw std::domain_error("negative power of polynomial");
        Poly result = one(*K_), base = *this;
        BigInt e = k;
        while (e > 0) {
            if (mpz_odd_p(e.get_mpz_t())) result = result * base;
            e >>= 1;
            if (e > 0) base = base * base;
        }
        return result;
    }

    // ---- evaluación y composición ------------------------------------------
    // f(x0) por la regla de Horner.  El polinomio cero vale 0 en todo punto.
    // Precondición (como en el resto de la librería): x0 es un elemento
    // normalizado de este mismo cuerpo.
    Elem evaluate(const Elem& x0) const {
        Elem acc = K_->zero();
        for (std::size_t i = c_.size(); i-- > 0; )
            acc = K_->add(K_->mul(acc, x0), c_[i]);
        return acc;
    }

    // Composición f(g(x)) por Horner con polinomios:
    //   f(g) = (...((c_n·g + c_{n-1})·g + c_{n-2})...)·g + c_0.
    // g debe estar sobre el mismo cuerpo (si no, std::invalid_argument).
    Poly compose(const Poly& g) const {
        require_same_field(g);
        Poly acc = zero(*K_);
        for (std::size_t i = c_.size(); i-- > 0; )
            acc = acc * g + constant(*K_, c_[i]);
        return acc;
    }

    // ¿Es f libre de cuadrados (sin factores irreducibles repetidos)?
    //   * el polinomio cero NO lo es (lo divide el cuadrado de cualquier cosa);
    //   * una constante no nula SÍ (es una unidad);
    //   * si deg f >= 1, f es libre de cuadrados  <=>  gcd(f, f') = 1.  En
    //     particular, si f' = 0 (f = g(x^p) = h^p en un cuerpo finito), f no es
    //     libre de cuadrados: gcd(f, 0) = f mónico tiene grado >= 1.
    bool is_squarefree() const {
        if (is_zero()) return false;
        if (degree() == 0) return true;
        Poly df = derivative();
        if (df.is_zero()) return false;
        return gcd(*this, df).degree() == 0;
    }

    // ---- gcd family --------------------------------------------------------
    // Fase 4 (nivel B): con grados grandes, gcd y gcd_ext usan el algoritmo
    // half-GCD (Thull-Yap), que obtiene la misma sucesión de cocientes que el
    // algoritmo de Euclides clásico en O(M(n) log n) en lugar de O(n^2).  Como
    // los cocientes son los mismos, el mcd y los coeficientes de Bézout son
    // idénticos a los de gcd_classic / gcd_ext_classic (los originales).
    static Poly gcd(Poly a, Poly b) {
        a.require_same_field(b);
        const int umbral = static_cast<int>(umbral_hgcd(a.field(), false));
        while (!b.is_zero()) {
            if (a.degree() > b.degree() && a.degree() >= umbral) {
                const Mat2 R = hgcd(a, b);
                aplicar(R, a, b);
                if (b.is_zero()) break;
            }
            Poly r = a % b;
            a = std::move(b);
            b = std::move(r);
        }
        return a.is_zero() ? a : a.monic();
    }

    // Algoritmo de Euclides clásico (el original).
    static Poly gcd_classic(Poly a, Poly b) {
        a.require_same_field(b);
        while (!b.is_zero()) {
            Poly r = a % b;
            a = std::move(b);
            b = std::move(r);
        }
        return a.is_zero() ? a : a.monic();
    }

    struct ExtGcd { Poly g, x, y; };  // g = x*a + y*b, g monic
    static ExtGcd gcd_ext(const Poly& a, const Poly& b) {
        a.require_same_field(b);
        const Field& K = a.field();
        const int umbral = static_cast<int>(umbral_hgcd(K, true));
        if (std::max(a.degree(), b.degree()) < umbral) return gcd_ext_classic(a, b);
        // El gcd_ext original multiplica q·x1 y q·y1 con la regla Auto, que puede
        // calcular la raíz de la FFT (consumiendo el generador) a mitad del
        // algoritmo.  Para que la secuencia aleatoria no cambie, el camino rápido
        // solo se usa si esa raíz ya existe o si el cuerpo no admite una FFT de
        // tamaño 64 (entonces el original nunca la pediría).
        if (K.fft_cache_.empty() && fft::fft_exponent(K, 64) >= 0) return gcd_ext_classic(a, b);
        // M acumula la matriz de cocientes: (r0, r1) = M·(a, b) en todo momento.
        Poly r0 = a, r1 = b;
        Mat2 M = identidad(K);
        while (!r1.is_zero()) {
            if (r0.degree() > r1.degree() && r0.degree() >= umbral) {
                const Mat2 R = hgcd(r0, r1);
                aplicar(R, r0, r1);
                M = mul_mat(R, M);
                if (r1.is_zero()) break;
            }
            paso(r0, r1, M);
        }
        if (r0.is_zero()) return {zero(K), zero(K), zero(K)};
        Elem inv = K.inv(r0.lead());
        return {r0.scalar_mul(inv), M.m00.scalar_mul(inv), M.m01.scalar_mul(inv)};
    }

    // Algoritmo de Euclides extendido clásico (el original).
    static ExtGcd gcd_ext_classic(const Poly& a, const Poly& b) {
        a.require_same_field(b);
        const Field& K = a.field();
        Poly x0 = one(K), x1 = zero(K), y0 = zero(K), y1 = one(K), r0 = a, r1 = b;
        while (!r1.is_zero()) {
            auto [q, r] = r0.divmod(r1);
            r0 = std::move(r1); r1 = std::move(r);
            Poly tx = x0 - q * x1; x0 = std::move(x1); x1 = std::move(tx);
            Poly ty = y0 - q * y1; y0 = std::move(y1); y1 = std::move(ty);
        }
        if (r0.is_zero()) return {zero(K), zero(K), zero(K)};
        Elem inv = K.inv(r0.lead());
        return {r0.scalar_mul(inv), x0.scalar_mul(inv), y0.scalar_mul(inv)};
    }

    // a^{-1} mod m  (requires gcd(a, m) = 1).
    static Poly inv_mod(const Poly& a, const Poly& m) {
        auto e = gcd_ext(a, m);  // comprueba el cuerpo común
        if (e.g.degree() != 0 || !a.field().is_one(e.g.c_[0]))
            throw std::domain_error("polynomial not invertible modulo m");
        return e.x % m;
    }

    // a^k mod m.
    // Fase 4: las reducciones módulo m de los productos (grado <= 2·deg m - 2)
    // usan un Reducer, que precalcula UNA vez el inverso de rev(m) si conviene
    // la división de Newton; si no, es el % de siempre.  Los productos y la
    // cadena de cuadrados no cambian.
    static Poly pow_mod(Poly a, BigInt k, const Poly& m) {
        a.require_same_field(m);
        if (m.is_zero()) throw division_by_zero("division by zero polynomial");
        return pow_mod(std::move(a), std::move(k), Reducer(m));
    }

    // Reducción repetida módulo un m fijo (ver pow_mod).  La primera vez que
    // prefer_newton lo aconseja (con el inverso amortizado), calcula y guarda
    // rev(m)^{-1} mod x^{deg m}; desde entonces reduce(a) divide con ese
    // inverso (a de grado <= 2·deg m - 1, el caso de los productos de
    // pow_mod).  El cálculo es perezoso porque la decisión depende de que el
    // cuerpo ya tenga tablas FFT, que los productos de pow_mod construyen.
    class Reducer {
    public:
        explicit Reducer(const Poly& m) : m_(m) {}
        const Poly& modulus() const { return m_; }
        Poly reduce(const Poly& a) const {
            const int n = m_.degree();
            if (a.degree() < n) return a;
            const std::size_t nn = static_cast<std::size_t>(n);
            const std::size_t d1 = static_cast<std::size_t>(a.degree() - n) + 1;
            if (d1 <= nn && prefer_newton(m_.field(), nn, d1, true)) {
                if (inv_.empty()) inv_ = inverse_series(m_.field(), reversed_head(m_, nn), nn);
                return a.divmod_with_inverse(
                            m_, std::vector<Elem>(inv_.begin(), inv_.begin() + static_cast<std::ptrdiff_t>(d1)))
                    .second;
            }
            return a.divmod_school(m_).second;
        }
    private:
        Poly m_;
        mutable std::vector<Elem> inv_;  // rev(m)^{-1} mod x^{deg m}, o vacío si aún no hace falta
    };

    // Composición modular de Brent-Kung con h fijo: precalcula h^0..h^{k-1}
    // mod f (pasos pequeños) y G = h^k mod f (paso gigante).  Para g de grado
    // < deg f, g(h) mod f = Σ_j B_j·G^j con B_j = Σ_{i<k} g_{jk+i}·h^i, que se
    // evalúa por Horner en G.  Cada B_j es una combinación lineal de las
    // potencias precalculadas (acumulación perezosa, ~k·n operaciones de
    // cuerpo); los pasos gigantes cuestan ceil(deg g / k) productos modulares.
    class Composer {
    public:
        Composer(const Poly& h, const Reducer& red, std::size_t k) : red_(red) {
            const Field& K = h.field();
            const std::size_t n = static_cast<std::size_t>(red_.modulus().degree());
            k_ = std::max<std::size_t>(1, std::min(k, n));
            Poly hr = red_.reduce(h.degree() >= 2 * static_cast<int>(n) - 1 ? h % red_.modulus() : h);
            pot_.reserve(k_);
            pot_.push_back(red_.reduce(one(K)));
            for (std::size_t i = 1; i < k_; ++i) pot_.push_back(red_.reduce(pot_.back().mul_no_rng(hr)));
            gigante_ = red_.reduce(pot_.back().mul_no_rng(hr));
        }
        std::size_t k() const { return k_; }
        // g(h) mod f
        Poly compose(const Poly& g) const {
            const Field& K = red_.modulus().field();
            if (g.is_zero()) return zero(K);
            const std::size_t L = g.c_.size();
            const std::size_t bloques = (L + k_ - 1) / k_;
            Poly acc = zero(K);
            for (std::size_t j = bloques; j-- > 0;) {
                Poly B = combina(g, j * k_, std::min(k_, L - j * k_));
                acc = (j + 1 == bloques) ? std::move(B) : red_.reduce(acc.mul_no_rng(gigante_)) + B;
            }
            return acc;
        }

    private:
        Reducer red_;
        std::size_t k_ = 1;
        std::vector<Poly> pot_;  // h^0..h^{k-1} mod f
        Poly gigante_{red_.modulus().field()};  // h^k mod f

        // Σ_{i<cnt} g_{ini+i}·h^i
        Poly combina(const Poly& g, std::size_t ini, std::size_t cnt) const {
            const Field& K = red_.modulus().field();
            std::size_t n = 0;
            for (std::size_t i = 0; i < cnt; ++i) n = std::max(n, pot_[i].c_.size());
            if constexpr (LazyAccumulation<Field>) {
                std::vector<typename Field::Acc> s(n);
                for (auto& a : s) K.acc_clear(a);
                for (std::size_t i = 0; i < cnt; ++i) {
                    const Elem& c = g.c_[ini + i];
                    if (K.is_zero(c)) continue;
                    const auto& P = pot_[i].c_;
                    for (std::size_t t = 0; t < P.size(); ++t) K.acc_fma(s[t], c, P[t]);
                }
                std::vector<Elem> r;
                r.reserve(n);
                for (auto& a : s) r.push_back(K.acc_reduce(a));
                return Poly(K, std::move(r));
            } else {
                std::vector<Elem> r(n, K.zero());
                for (std::size_t i = 0; i < cnt; ++i) {
                    const Elem& c = g.c_[ini + i];
                    if (K.is_zero(c)) continue;
                    const auto& P = pot_[i].c_;
                    for (std::size_t t = 0; t < P.size(); ++t) r[t] = K.add(r[t], K.mul(c, P[t]));
                }
                return Poly(K, std::move(r));
            }
        }
    };

    // Sucesión de Frobenius H_i = x^{q^i} mod f, paso a paso (is_irreducible,
    // distinct_degree).  Empieza con H_i = H_{i-1}^q (pow_mod, el original) y
    // pasa a la composición modular (Composer con h = H_1) cuando lo gastado en
    // potencias iguala el coste estimado de precalcular las potencias de H_1
    // (regla de alquilar o comprar: nunca más del doble de lo óptimo, y nada si
    // el bucle termina pronto, como con casi todos los reducibles).  Las
    // potencias guardadas se limitan a kFrobPresupuesto coeficientes; con menos
    // potencias hacen falta pasos gigantes, y no se compone si serían más que
    // los productos de un pow_mod.
    class Frobenius {
    public:
        Frobenius(const Poly& f, FrobeniusAlgo algo = FrobeniusAlgo::Auto)
            : red_(f), q_(f.field().cardinality()), algo_(algo) {
            productos_q_ = mpz_sizeinbase(q_.get_mpz_t(), 2) + mpz_popcount(q_.get_mpz_t()) - 1;
            const std::size_t n = static_cast<std::size_t>(std::max(1, f.degree()));
            k_ = std::max<std::size_t>(1, std::min(n, kFrobPresupuesto / n));
        }
        const Reducer& reducer() const { return red_; }
        // Siguiente término: de H = H_{i-1} a H_i.
        Poly siguiente(const Poly& H) {
            if (comp_) return comp_->compose(H);
            Poly r = pow_mod(H, q_, red_);
            ++pasos_;
            if (!h1_) h1_ = r;  // el primer término es H_1 = x^q mod f (se empieza en H_0 = x)
            if (h1_ && quiere_componer()) comp_.emplace(*h1_, red_, k_);
            return r;
        }
        // f cambia por un divisor g (distinct_degree): H_1 mod g sigue siendo x^q mod g.
        void cambiar_modulo(const Poly& g) {
            red_ = Reducer(g);
            const std::size_t n = static_cast<std::size_t>(std::max(1, g.degree()));
            k_ = std::max<std::size_t>(1, std::min(n, kFrobPresupuesto / n));
            if (h1_) h1_ = *h1_ % g;
            comp_.reset();
            pasos_ = 0;
        }

    private:
        Reducer red_;
        BigInt q_;
        FrobeniusAlgo algo_;
        std::size_t productos_q_ = 1, k_ = 1, pasos_ = 0;
        std::optional<Poly> h1_;
        std::optional<Composer> comp_;

        bool quiere_componer() const {
            if (algo_ == FrobeniusAlgo::Potencia) return false;
            if (algo_ == FrobeniusAlgo::Composicion) return true;
            if (productos_q_ < kFrobMinProductos) return false;
            // Coste de precalcular: k productos modulares (+ el paso gigante).
            const std::size_t n = static_cast<std::size_t>(std::max(1, red_.modulus().degree()));
            const std::size_t gigantes = (n + k_ - 1) / k_;
            if (gigantes > productos_q_) return false;  // cada paso no ahorraría productos
            return pasos_ * productos_q_ >= k_ + 1;
        }
    };

    // pow_mod con un Reducer ya construido (is_irreducible y distinct_degree
    // reducen muchas veces módulo el mismo polinomio).
    static Poly pow_mod(Poly a, BigInt k, const Reducer& red) {
        const Poly& m = red.modulus();
        a.require_same_field(m);
        const Field& K = a.field();
        if (k < 0) { a = inv_mod(a, m); k = -k; }
        Poly result = red.reduce(one(K));
        a = a % m;
        while (k > 0) {
            if (mpz_odd_p(k.get_mpz_t())) result = red.reduce(result * a);
            k >>= 1;
            if (k > 0) a = red.reduce(a * a);
        }
        return result;
    }

    // ---- Test de irreducibilidad de Ben-Or -----------------------------------
    // f de grado n es irreducible sobre F_q  <=>  no tiene ningún factor
    // irreducible de grado i <= n/2  <=>  gcd(f, x^{q^i} - x) = 1 para todo
    // i = 1..floor(n/2) (x^{q^i} - x es el producto de todos los irreducibles
    // mónicos de grado divisor de i).  (No es el test de Rabin, que usa los
    // divisores primos de n.)
    // Fase 4: el original (y la práctica de Python) comprobaba además, tras el
    // bucle, que x^{q^n} == x (mod f) calculando esa potencia desde cero.  Con
    // el bucle completo esa condición es redundante (lo anterior ya prueba la
    // irreducibilidad) y costaba unos dos tercios del tiempo con entradas
    // irreducibles: se suprime.  Las reducciones usan un Reducer (pow_mod).
    // Fase 4 (nivel B): la sucesión H_i se calcula con Frobenius (pow_mod o
    // composición modular, según convenga); los H_i son los mismos.
    bool is_irreducible(FrobeniusAlgo algo = FrobeniusAlgo::Auto) const {
        const int n = degree();
        if (n <= 0) return false;
        Poly f = *this;
        Poly x = variable(*K_);
        Frobenius frob(f, algo);
        Poly H = x % f;                       // H = x^{q^0}
        for (int i = 1; i <= n / 2; ++i) {
            H = frob.siguiente(H);            // H = x^{q^i} (mod f)
            if (gcd(f, H - x).degree() != 0) return false;
        }
        return true;
    }

    // ---- factorization (Cantor-Zassenhaus) ---------------------------------
    struct Factorization {
        std::vector<std::pair<Poly, int>> factors;  // (monic irreducible, multiplicity)
        Elem leading;                                // unit: leading coeff of original f
    };

    // Square-free part: the product of the distinct monic irreducible factors
    // of *this (a.k.a. radical).  Port of sqfree_fact_*.
    // Casos triviales: rad(0) = 0 y, para una constante no nula, rad(c) = 1
    // (antes, con una constante, la rama f' = 0 se llamaba a sí misma sin fin).
    Poly radical() const {
        const Field& K = *K_;
        if (is_zero()) return zero(K);
        if (degree() == 0) return one(K);

        const BigInt p = K.characteristic();
        const unsigned long pu = fits_ulong(p) ? to_ulong(p) : 0;  // step for x^p case

        Poly df = derivative();
        if (df.is_zero()) {
            // f(x) = g(x^p); recover g via the p-th root of every p-th coeff.
            require_pth_step(pu);
            std::vector<Elem> g;
            for (std::size_t i = 0; i < c_.size(); i += pu) g.push_back(K.pth_root(c_[i]));
            return Poly(K, std::move(g)).radical();
        }

        Poly c = gcd(*this, df);
        Poly g = (*this / c).monic();

        // Strip every factor of g out of f; what remains is a perfect p-th power.
        Poly residue = *this;
        while (true) {
            Poly y = gcd(residue, g);
            if (y.degree() <= 0) break;
            residue = residue / y;
        }
        if (residue.degree() > 0) {
            require_pth_step(pu);
            std::vector<Elem> h;
            for (std::size_t i = 0; i < residue.c_.size(); i += pu) h.push_back(K.pth_root(residue.c_[i]));
            Poly g2 = Poly(K, std::move(h)).radical();
            return (g * g2) / gcd(g, g2);  // lcm(g, g2)
        }
        return g;
    }

    // Distinct-degree factorization of a square-free monic g.  Returns pairs
    // (d, product-of-all-degree-d-irreducible-factors), only for d that occur.
    // Precondiciones: *this no nulo, de grado >= 1 (se comprueban: si no,
    // std::invalid_argument) y libre de cuadrados (NO se comprueba por coste;
    // con factores repetidos el resultado no es una factorización válida, pero
    // el bucle termina igualmente porque d está acotado por el grado).
    std::vector<std::pair<int, Poly>> distinct_degree() const {
        if (is_zero())
            throw std::invalid_argument("distinct_degree: el polinomio cero no se puede factorizar");
        if (degree() < 1)
            throw std::invalid_argument("distinct_degree: se requiere un polinomio de grado >= 1");
        const Field& K = *K_;
        std::vector<std::pair<int, Poly>> out;
        Poly g = monic();
        Poly x = variable(K);
        Poly H = x % g;  // x^{q^0}
        Frobenius frob(g);  // fase 4: pow_mod o composición modular (ver Frobenius)
        int d = 1;
        while (g.degree() >= 2 * d) {
            H = frob.siguiente(H);           // H = x^{q^d} (mod g)
            Poly hd = gcd(H - x, g);
            if (hd.degree() > 0) {
                out.emplace_back(d, hd.monic());
                g = g / hd;
                H = H % g;
                frob.cambiar_modulo(g);
            }
            ++d;
        }
        if (g.degree() > 0) out.emplace_back(g.degree(), g.monic());
        return out;
    }

    // Equal-degree factorization: h is a product of distinct monic irreducible
    // factors all of degree r; returns that list of irreducibles.  Works for
    // odd q (the (q^r-1)/2 split) and for q even (the trace split).
    // Precondiciones:
    //   * se comprueban (std::invalid_argument): r >= 1, *this no nulo, de grado
    //     >= 1 y múltiplo de r; también que cada trozo intermedio tenga grado
    //     múltiplo de r (en una ejecución válida siempre lo tiene);
    //   * no se comprueban por coste: que *this sea libre de cuadrados y que
    //     todos sus factores irreducibles tengan grado exactamente r.  Si eso
    //     falla, la división aleatoria puede no encontrar nunca un factor
    //     propio; para que no haya un bucle infinito, tras kMaxSplitAttempts
    //     intentos fallidos se lanza std::invalid_argument.
    std::vector<Poly> equal_degree(int r) const {
        if (r < 1)
            throw std::invalid_argument("equal_degree: r debe ser >= 1");
        if (is_zero())
            throw std::invalid_argument("equal_degree: el polinomio cero no se puede factorizar");
        if (degree() < 1)
            throw std::invalid_argument("equal_degree: se requiere un polinomio de grado >= 1");
        if (degree() % r != 0)
            throw std::invalid_argument("equal_degree: el grado debe ser múltiplo de r");
        std::vector<Poly> factors;
        std::vector<Poly> stack;
        stack.push_back(monic());
        while (!stack.empty()) {
            Poly g = std::move(stack.back());
            stack.pop_back();
            if (g.degree() <= 0) continue;
            if (g.degree() == r) { factors.push_back(g.monic()); continue; }
            if (g.degree() % r != 0)
                throw std::invalid_argument(
                    "equal_degree: la entrada no es producto de irreducibles distintos de grado r");
            Poly d = equal_degree_split(g, r);
            stack.push_back(d);
            stack.push_back(g / d);
        }
        return factors;
    }

    // Largest e with u^e | *this.  Port of multiplicidad_*.
    // Precondiciones (std::invalid_argument si no se cumplen): u de grado >= 1
    // (con u constante no nula u^e divide siempre: el bucle no terminaba; u = 0
    // no divide) y *this no nulo (u^e divide a 0 para todo e: multiplicidad
    // infinita).
    int multiplicity(const Poly& u) const {
        require_same_field(u);
        if (u.degree() < 1)
            throw std::invalid_argument("multiplicity: el divisor debe tener grado >= 1");
        if (is_zero())
            throw std::invalid_argument("multiplicity: la multiplicidad en el polinomio cero es infinita");
        int e = 0;
        Poly f = *this;
        while (!f.is_zero()) {
            auto [q, r] = f.divmod(u);
            if (!r.is_zero()) break;
            ++e;
            f = std::move(q);
        }
        return e;
    }

    // Full Cantor-Zassenhaus factorization into monic irreducibles with
    // multiplicities, plus the leading coefficient (the field unit).
    // Los factores se devuelven en ORDEN CANÓNICO (ver canonical_less): por
    // grado ascendente y, a igual grado, por coeficientes de mayor a menor
    // grado.  El algoritmo no cambia; solo se ordena el resultado al final.
    Factorization factor() const {
        if (is_zero()) throw std::domain_error("cannot factor the zero polynomial");
        Factorization F;
        F.leading = lead();
        if (degree() == 0) return F;

        Poly g = radical();
        std::vector<Poly> irreducibles;
        for (auto& [d, prod] : g.distinct_degree()) {
            auto part = prod.equal_degree(d);
            irreducibles.insert(irreducibles.end(), part.begin(), part.end());
        }
        for (auto& u : irreducibles)
            F.factors.emplace_back(u, multiplicity(u));
        std::sort(F.factors.begin(), F.factors.end(),
                  [](const std::pair<Poly, int>& a, const std::pair<Poly, int>& b) {
                      return canonical_less(a.first, b.first);
                  });
        return F;
    }

    // Orden canónico entre polinomios: primero el grado (ascendente); a igual
    // grado, se comparan los coeficientes desde el de mayor grado hasta el
    // término independiente con operator< de Elem (BigInt en F_p; en F_q,
    // std::vector<BigInt>, comparación lexicográfica de sus coordenadas).
    // Es un orden estricto débil; no comprueba el cuerpo.
    static bool canonical_less(const Poly& a, const Poly& b) {
        if (a.degree() != b.degree()) return a.degree() < b.degree();
        for (std::size_t i = a.c_.size(); i-- > 0; ) {
            if (a.c_[i] < b.c_[i]) return true;
            if (b.c_[i] < a.c_[i]) return false;
        }
        return false;
    }

    // ---- formatting --------------------------------------------------------
    std::string to_string(const std::string& var = "x") const {
        if (is_zero()) return "0";
        std::string s;
        bool first = true;
        for (std::size_t i = 0; i < c_.size(); ++i) {
            if (K_->is_zero(c_[i])) continue;
            if (!first) s += " + ";
            first = false;
            std::string cs = "(" + K_->to_string(c_[i]) + ")";
            if (i == 0) s += cs;
            else if (i == 1) s += cs + "*" + var;
            else s += cs + "*" + var + "^" + std::to_string(i);
        }
        return s;
    }

private:
    const Field* K_;
    std::vector<Elem> c_;

    void trim() {
        while (!c_.empty() && K_->is_zero(c_.back())) c_.pop_back();
    }

    // Comprobación de cuerpo común.  Camino rápido: el mismo objeto (una
    // comparación de punteros); si no, mismos parámetros (same_field).
    void require_same_field(const Poly& o) const {
        if (K_ != o.K_ && !detail::fields_match(*K_, *o.K_)) [[unlikely]]
            throw std::invalid_argument("polinomios sobre cuerpos distintos");
    }

    // En radical(), el paso de x^p a x exige recorrer los coeficientes de p en
    // p.  Si p no cabe en unsigned long (32 bits en Windows) esa rama solo se
    // alcanzaría con grado >= p, imposible en memoria; se protege igualmente
    // para no entrar nunca en un bucle con paso 0.
    static void require_pth_step(unsigned long pu) {
        if (pu == 0)
            throw std::length_error("radical: característica demasiado grande para la rama f(x) = g(x^p)");
    }

    // coefficients [lo, hi) as a polynomial (drops the x^lo..x^{hi-1} window down to x^0).
    // ---- Umbrales de despacho Auto ------------------------------------------
    // ¿Conviene la FFT (supuesto que cabe) para un producto de operandos con
    // `small` y `large` coeficientes (L = small + large - 1)?  Depende de la
    // implementación de la NTT del cuerpo y, sobre todo, del operando MENOR:
    // con uno diminuto (p. ej. q·x1 en gcd_ext, con q de grado 1) el escolar
    // cuesta O(large) y la FFT O(L log L).  Calibrado con medidas en Zen 3
    // (g++ 15.2 -O3, GMP 6.3) de escolar/Karatsuba/FFT en una rejilla de
    // tamaños sobre F_998244353, el r de 255 bits de BLS12-381 y
    // F_((2^61-1)^2):
    //   * Word (F_p, p < 2^62): la FFT gana desde L = 2-4 en todos los casos.
    //   * Limb (F_p, p >= 2^62): FFT si small >= 16, o small >= 8 y L <= 128.
    //   * Generic (F_q y resto): FFT si small >= 64, o small >= 32 y
    //     large <= 2·small; con small <= 16 la FFT era hasta 10-20x más lenta.
    //   * Fp64 (elementos en palabras de máquina, NTT en palabras): el escolar
    //     con acumulación perezosa de 128 bits es tan barato que la FFT solo
    //     gana con small >= 128, o con small >= 48 y large <= 4·small
    //     (medido en F_998244353 con g++ 15.2 -O3 en Zen 3).
    static bool prefer_fft(const Field& K, std::size_t small, std::size_t large, std::size_t L) {
        if constexpr (fft::NativeWordField<Field>) {
            (void)K; (void)L;
            return small >= 128 || (small >= 48 && large <= 4 * small);
        }
        switch (fft::backend(K)) {
            case fft::Backend::Word: return L >= 4;
            case fft::Backend::Limb: return small >= 16 || (small >= 8 && L <= 128);
            case fft::Backend::Generic:
            default: return small >= 64 || (small >= 32 && large <= 2 * small);
        }
    }

    // ---- Frobenius (fase 4, nivel B) ----------------------------------------
    // Productos modulares por paso de pow_mod a partir de los cuales Auto
    // considera la composición modular, y número máximo de coeficientes
    // guardados en las potencias de H_1 (k·n); calibrado con
    // bench/calibra_frobenius.cpp.
    // La composición ganó en todos los cuerpos medidos con entradas irreducibles
    // (de 1,2x con q = 2 a 19x con q = 2^61 - 1), así que el mínimo es 2.
    static constexpr std::size_t kFrobMinProductos = 2;
    static constexpr std::size_t kFrobPresupuesto = std::is_arithmetic_v<Elem> ? (std::size_t{1} << 22)
                                                                             : (std::size_t{1} << 18);

    // ---- Half-GCD (fase 4, nivel B) ------------------------------------------
    // Por debajo de este grado, la recursión de hgcd usa los pasos clásicos.
    static constexpr std::size_t kHgcdMin = 128;
    // Grado a partir del cual gcd (ext = false) o gcd_ext (ext = true) usan
    // half-GCD.  Medido con bench/calibra_hgcd.cpp (Zen 3, entradas variadas):
    // gcd_ext gana antes porque el clásico también arrastra los cofactores.
    //   * Fp64: depende del tamaño de p (cuánto gana Kronecker en los productos):
    //     p < 2^16: 384 / 128;  p < 2^32: 1024 / 384;  si no: 2048 / 512.
    //   * Fp con GMP: 384 / 128.   * F_q: 384 / 192.
    static std::size_t umbral_hgcd(const Field& K, bool ext) {
        if constexpr (fft::NativeWordField<Field>) {
            const std::size_t bits = mpz_sizeinbase(K.characteristic().get_mpz_t(), 2);
            if (bits <= 16) return ext ? 128 : 384;
            if (bits <= 32) return ext ? 384 : 1024;
            return ext ? 512 : 2048;
        } else if constexpr (std::is_same_v<Elem, BigInt>) {
            (void)K;
            return ext ? 128 : 384;
        } else {
            (void)K;
            return ext ? 192 : 384;
        }
    }

    // Matriz 2x2 de polinomios [[m00, m01], [m10, m11]].
    struct Mat2 { Poly m00, m01, m10, m11; };
    static Mat2 identidad(const Field& K) { return {one(K), zero(K), zero(K), one(K)}; }
    // A·B.  Los productos no calculan la raíz de la FFT (mul_no_rng): el gcd
    // original no multiplicaba y no debe consumir el generador aleatorio.
    static Mat2 mul_mat(const Mat2& A, const Mat2& B) {
        return {A.m00.mul_no_rng(B.m00) + A.m01.mul_no_rng(B.m10), A.m00.mul_no_rng(B.m01) + A.m01.mul_no_rng(B.m11),
                A.m10.mul_no_rng(B.m00) + A.m11.mul_no_rng(B.m10), A.m10.mul_no_rng(B.m01) + A.m11.mul_no_rng(B.m11)};
    }
    // (a, b) <- M·(a, b)
    static void aplicar(const Mat2& M, Poly& a, Poly& b) {
        Poly na = M.m00.mul_no_rng(a) + M.m01.mul_no_rng(b);
        Poly nb = M.m10.mul_no_rng(a) + M.m11.mul_no_rng(b);
        a = std::move(na);
        b = std::move(nb);
    }
    // Paso de Euclides (a, b) <- (b, a mod b) acumulado en M: M <- [[0, 1], [1, -q]]·M,
    // es decir, M' = [[M10, M11], [M00 - q·M10, M01 - q·M11]] (2 productos por q,
    // que suele ser de grado 1, en lugar de un producto de matrices completo).
    static void paso(Poly& a, Poly& b, Mat2& M) {
        auto [q, r] = a.divmod(b);
        a = std::move(b);
        b = std::move(r);
        Poly n10 = M.m00 - q.mul_no_rng(M.m10);
        Poly n11 = M.m01 - q.mul_no_rng(M.m11);
        M.m00 = std::move(M.m10);
        M.m01 = std::move(M.m11);
        M.m10 = std::move(n10);
        M.m11 = std::move(n11);
    }
    // p div x^k (los coeficientes de grado >= k, desplazados).
    static Poly div_xk(const Poly& p, std::size_t k) {
        if (p.c_.size() <= k) return zero(p.field());
        return Poly(p.field(), std::vector<Elem>(p.c_.begin() + static_cast<std::ptrdiff_t>(k), p.c_.end()));
    }

    // HGCD (Thull-Yap): con deg a = n > deg b, devuelve la matriz M de los
    // cocientes de Euclides tales que (r_j, r_{j+1}) = M·(a, b) cumple
    // deg r_j >= m > deg r_{j+1}, con m = ceil(n/2).  Los cocientes de esa
    // parte de la sucesión dependen solo de los coeficientes altos, por eso se
    // obtienen recursivamente de a div x^m y b div x^m.
    static Mat2 hgcd(const Poly& a, const Poly& b) {
        const Field& K = a.field();
        const int n = a.degree();
        const int m = (n + 1) / 2;
        if (b.is_zero() || b.degree() < m) return identidad(K);
        if (n < static_cast<int>(kHgcdMin)) {
            // Caso base: pasos clásicos hasta bajar de grado m.
            Mat2 M = identidad(K);
            Poly A = a, B = b;
            while (!B.is_zero() && B.degree() >= m) paso(A, B, M);
            return M;
        }
        const Mat2 R = hgcd(div_xk(a, static_cast<std::size_t>(m)), div_xk(b, static_cast<std::size_t>(m)));
        Poly A = a, B = b;
        aplicar(R, A, B);
        if (B.is_zero() || B.degree() < m) return R;
        Mat2 QR = R;
        paso(A, B, QR);
        if (B.is_zero() || B.degree() < m) return QR;
        const int k = 2 * m - A.degree();
        if (k < 0) {
            // No debería ocurrir (Thull-Yap); por seguridad, pasos clásicos.
            while (!B.is_zero() && B.degree() >= m) paso(A, B, QR);
            return QR;
        }
        const Mat2 S = hgcd(div_xk(A, static_cast<std::size_t>(k)), div_xk(B, static_cast<std::size_t>(k)));
        return mul_mat(S, QR);
    }

    // ---- Kronecker (fase 4, nivel B) ----------------------------------------
    // ¿Conviene el producto de Kronecker sobre Z (kronecker.hpp) para un
    // operando menor de `small` coeficientes?  Calibrado con
    // bench/calibra_kronecker.cpp (entradas variadas, con y sin hilos, Zen 3):
    //   * F_q (FqT) y F_p con GMP: gana desde 16-32 coeficientes (1,5x a 25x;
    //     a Karatsuba sobre elementos caros y a la NTT multi-limb).
    //   * Fp64 sin NTT: desde 32 si p < 2^32 (hasta 20x con p = 101); con p
    //     mayor los huecos son de ~2·62 bits y gana desde 256.
    //   * Fp64 con NTT: solo entre 128 y 1024 (1,5x a 2,5x); por debajo gana
    //     el escolar y por encima la NTT en palabras (con hilos) empata o gana.
    static bool prefer_kronecker(const Field& K, std::size_t small, bool can_fft) {
        if (small < 32) return false;
        if constexpr (ks::Extension<Field>) {
            (void)K; (void)can_fft;
            return true;
        } else if constexpr (std::is_same_v<Elem, BigInt>) {
            (void)K; (void)can_fft;
            return true;
        } else if constexpr (fft::NativeWordField<Field>) {
            if (can_fft) return small >= 128 && small <= 1024;
            return mpz_sizeinbase(K.characteristic().get_mpz_t(), 2) <= 32 || small >= 256;
        } else {
            (void)K; (void)can_fft;
            return false;
        }
    }

    // ---- División de Newton (fase 4) ----------------------------------------
    // ¿Conviene la división de Newton para un divisor de grado n y un cociente
    // de d1 coeficientes?  La escolar (con acumulación perezosa) cuesta ~n·d1
    // productos de elementos; Newton, unos pocos productos de polinomios de
    // tamaño max(n, d1): 2 con el inverso ya calculado (`amortizado`, Reducer)
    // y unos 6 si cada división lo calcula (divmod).  Solo compensa con
    // productos subcuadráticos rápidos, y el punto de corte depende de si hay
    // FFT y de lo que cueste cada producto de elementos.  Umbrales medidos con
    // bench/calibra_newton.cpp en Zen 3 CON ENTRADAS VARIADAS (con una entrada
    // repetida el predictor de saltos acelera la NTT hasta 3 veces y Newton
    // parece mejor de lo que es; ver docs/optimizacion/fase4.md).
    // Nivel B: los productos usan Kronecker (prefer_kronecker), mucho más rápido
    // con p pequeño, y los umbrales bajan (bench/results/fase4/calibra_newton.txt;
    // la medida del nivel A, sin Kronecker, está en calibra_newton_nivel_a.txt):
    //   * Fp64: p < 2^16: 64 amortizado / 192 aislada; p < 2^32: 192 / 768;
    //     p mayor (huecos de Kronecker de ~2·62 bits): 768 / 1536.
    //   * Fp (GMP): 96 / 256.
    //   * F_q de grado e: con p < 2^32, max(32, 192/e) / max(64, 1024/e)
    //     (F_9: 96 / 512; F_256: 32 / 128; F_(3^30): 32 / 64); con p mayor
    //     (F_((2^61-1)^2)): 384 / 1024.
    static constexpr std::size_t kNunca = static_cast<std::size_t>(-1);
    static bool prefer_newton(const Field& K, std::size_t n, std::size_t d1, bool amortizado = false) {
        const std::size_t small = std::min(n, d1);
        if (small < 32) return false;  // por debajo, la escolar gana en todos los cuerpos medidos
        const std::size_t bits = mpz_sizeinbase(K.characteristic().get_mpz_t(), 2);
        std::size_t umbral = kNunca;
        if constexpr (fft::NativeWordField<Field>) {
            if (bits <= 16) umbral = amortizado ? 64 : 192;
            else if (bits <= 32) umbral = amortizado ? 192 : 768;
            else umbral = amortizado ? 768 : 1536;
        } else if constexpr (std::is_same_v<Elem, BigInt>) {
            umbral = amortizado ? 96 : 256;
        } else if constexpr (requires { K.extension_degree(); }) {
            const std::size_t e = std::max<std::size_t>(1, K.extension_degree());
            if (bits <= 32) umbral = amortizado ? std::max<std::size_t>(32, 192 / e) : std::max<std::size_t>(64, 1024 / e);
            else umbral = amortizado ? 384 : 1024;
        }
        return small >= umbral;
    }

    // Los d1 primeros coeficientes de rev(g) (de mayor a menor grado de g),
    // rellenando con ceros si d1 > deg g + 1.
    static std::vector<Elem> reversed_head(const Poly& g, std::size_t d1) {
        const Field& K = g.field();
        std::vector<Elem> r;
        r.reserve(d1);
        const std::size_t n1 = g.c_.size();
        for (std::size_t j = 0; j < d1; ++j) r.push_back(j < n1 ? g.c_[n1 - 1 - j] : K.zero());
        return r;
    }

    // Los n primeros coeficientes de P (rellenando con ceros).
    static std::vector<Elem> head(const Field& K, const Poly& P, std::size_t n) {
        std::vector<Elem> r(P.c_.begin(), P.c_.begin() + static_cast<std::ptrdiff_t>(std::min(n, P.c_.size())));
        r.resize(n, K.zero());
        return r;
    }

    // Producto con la regla Auto pero SIN adelantar el cálculo de la raíz de la
    // FFT: solo usa la FFT si las tablas del cuerpo ya existen.  Los productos
    // de la división de Newton no estaban en el original y no deben consumir el
    // generador aleatorio en otro momento (ver mul).
    Poly mul_no_rng(const Poly& o) const {
        if (is_zero() || o.is_zero()) return zero(*K_);
        const std::size_t na = c_.size(), nb = o.c_.size();
        const std::size_t L = na + nb - 1, small = std::min(na, nb);
        const bool can_fft = fft::fft_exponent(*K_, L) >= 0;
        if (prefer_kronecker(*K_, small, can_fft)) return mul_kronecker(o);
        const bool use_fft = !K_->fft_cache_.empty() && can_fft && prefer_fft(*K_, small, std::max(na, nb), L);
        if (use_fft) return Poly(*K_, fft::convolve(*K_, c_, o.c_));
        if (small >= kKaratsubaCutoff) return mul_karatsuba(o);
        return mul_school(o);
    }

    // ---- Karatsuba sobre vistas -------------------------------------------
    // Longitud de la vista a[0..n) sin ceros finales (lo que haría trim()).
    static std::size_t trimmed_len(const Field& K, const Elem* a, std::size_t n) {
        while (n > 0 && K.is_zero(a[n - 1])) --n;
        return n;
    }

    // x[0..nx) + y[0..ny), recortado (como operator+ sobre las mitades).
    static std::vector<Elem> add_views(const Field& K, const Elem* x, std::size_t nx,
                                       const Elem* y, std::size_t ny) {
        std::vector<Elem> r;
        r.reserve(std::max(nx, ny));
        const std::size_t nmin = std::min(nx, ny);
        for (std::size_t i = 0; i < nmin; ++i) r.push_back(K.add(x[i], y[i]));
        for (std::size_t i = nmin; i < nx; ++i) r.push_back(x[i]);
        for (std::size_t i = nmin; i < ny; ++i) r.push_back(y[i]);
        r.resize(trimmed_len(K, r.data(), r.size()));
        return r;
    }

    // Producto Karatsuba de las vistas a[0..na) y b[0..nb), ambas no vacías y
    // recortadas (último coeficiente no nulo).  Reproduce la recursión de la
    // versión con Poly: m = max(na, nb)/2; aL = a[0..m) y aH = a[m..na)
    // recortadas (vacías si m >= na), ídem para b; z0 = aL·bL, z2 = aH·bH,
    // z1 = (aL+aH)·(bL+bH) - z0 - z2; resultado z0 + z1·x^m + z2·x^{2m}.  El
    // resultado tiene longitud exacta na + nb - 1 (coeficiente principal
    // lead(a)·lead(b) != 0).
    // Si a y b son la misma vista (cuadrado), aL+aH y bL+bH coinciden: la suma
    // se calcula una vez (trabajo redundante idéntico).
    // Paralelismo: con par_depth > 0 y ambos operandos de al menos
    // kParKaratsubaMin coeficientes, z0 y z2 se calculan en otros dos hilos
    // mientras este calcula z1 (los tres subproductos son independientes y las
    // operaciones de cuerpo son puras); los hijos reciben par_depth - 1.
    static std::vector<Elem> karatsuba_rec(const Field& K, const Elem* a, std::size_t na,
                                           const Elem* b, std::size_t nb, int par_depth = 0) {
        const bool same = (a == b && na == nb);
        if (na <= kKaratsubaCutoff || nb <= kKaratsubaCutoff) {
            // Caso base: exactamente mul_school (se reutiliza para heredar
            // cualquier mejora de su bucle).  La copia de las vistas cuesta
            // O(na + nb), despreciable frente a los na·nb productos.
            Poly pa(K, std::vector<Elem>(a, a + na));
            if (same) return std::move(pa.mul_school(pa).c_);
            Poly pb(K, std::vector<Elem>(b, b + nb));
            return std::move(pa.mul_school(pb).c_);
        }

        const std::size_t m = std::max(na, nb) / 2;
        const std::size_t naL = trimmed_len(K, a, std::min(m, na));
        const std::size_t naH = (m < na) ? na - m : 0;  // ya recortada (a lo está)
        const std::size_t nbL = same ? naL : trimmed_len(K, b, std::min(m, nb));
        const std::size_t nbH = (m < nb) ? nb - m : 0;

        const bool par = par_depth > 0 && std::min(na, nb) >= kParKaratsubaMin;
        const int child = par ? par_depth - 1 : 0;
        auto calc_z0 = [&]() -> std::vector<Elem> {
            if (naL > 0 && nbL > 0) return karatsuba_rec(K, a, naL, b, nbL, child);
            return {};
        };
        auto calc_z2 = [&]() -> std::vector<Elem> {
            if (naH > 0 && nbH > 0) return karatsuba_rec(K, a + m, naH, b + m, nbH, child);
            return {};
        };
        auto calc_z1 = [&]() -> std::vector<Elem> {
            std::vector<Elem> sa = (naH > 0) ? add_views(K, a, naL, a + m, naH)
                                             : std::vector<Elem>(a, a + naL);
            if (same) {
                if (!sa.empty()) return karatsuba_rec(K, sa.data(), sa.size(), sa.data(), sa.size(), child);
                return {};
            }
            std::vector<Elem> sb = (nbH > 0) ? add_views(K, b, nbL, b + m, nbH)
                                             : std::vector<Elem>(b, b + nbL);
            if (!sa.empty() && !sb.empty())
                return karatsuba_rec(K, sa.data(), sa.size(), sb.data(), sb.size(), child);
            return {};
        };

        std::vector<Elem> z0, z2, z1;
        if (par) {
            auto f0 = std::async(std::launch::async, calc_z0);
            auto f2 = std::async(std::launch::async, calc_z2);
            z1 = calc_z1();
            z0 = f0.get();
            z2 = f2.get();
        } else {
            z0 = calc_z0();
            z2 = calc_z2();
            z1 = calc_z1();
        }
        // z1 -= z0; z1 -= z2 (en sitio) y recorte.
        const std::size_t n1 = std::max({z1.size(), z0.size(), z2.size()});
        z1.resize(n1, K.zero());
        for (std::size_t i = 0; i < z0.size(); ++i) z1[i] = K.sub(z1[i], z0[i]);
        for (std::size_t i = 0; i < z2.size(); ++i) z1[i] = K.sub(z1[i], z2[i]);
        z1.resize(trimmed_len(K, z1.data(), z1.size()));

        // Recombinación en el vector resultado (se reutiliza el de z0).
        std::vector<Elem> r = std::move(z0);
        r.resize(na + nb - 1, K.zero());
        for (std::size_t i = 0; i < z1.size(); ++i) r[m + i] = K.add(r[m + i], z1[i]);
        for (std::size_t i = 0; i < z2.size(); ++i) r[2 * m + i] = K.add(r[2 * m + i], z2[i]);
        return r;
    }

    // One Cantor-Zassenhaus splitting step: given g (product of distinct
    // degree-r irreducibles, deg g > r), return a nontrivial factor.
    // Si tras kMaxSplitAttempts intentos no hay división, la entrada no cumple
    // la precondición (ver equal_degree) y se lanza std::invalid_argument.
    static Poly equal_degree_split(const Poly& g, int r) {
        const Field& K = g.field();
        const BigInt q = K.cardinality();
        const bool even_char = (K.characteristic() == 2);
        // For char 2 we need the number of squarings in the trace map:
        // s = log2(q), trace over F_{q^r}/F_2 uses s*r squarings.
        const unsigned long s = even_char ? (bit_length(q) - 1) : 0;

        for (int attempt = 0; attempt < kMaxSplitAttempts; ++attempt) {
            Poly a = Poly::random(K, g.degree() - 1) % g;
            if (a.degree() <= 0) continue;
            Poly b(K);
            if (!even_char) {
                BigInt e = (ipow(q, static_cast<unsigned long>(r)) - 1) / 2;
                b = pow_mod(a, e, g) - one(K);
            } else {
                Poly trace = zero(K), cur = a;
                for (unsigned long i = 0; i < s * static_cast<unsigned long>(r); ++i) {
                    trace = (trace + cur) % g;
                    cur = (cur * cur) % g;
                }
                b = trace;
            }
            Poly d = gcd(b, g);
            if (d.degree() > 0 && d.degree() < g.degree()) return d;
        }
        throw std::invalid_argument(
            "equal_degree: no se encuentra un factor propio; la entrada no es producto de "
            "irreducibles distintos de grado r");
    }
};

}  // namespace alcp
