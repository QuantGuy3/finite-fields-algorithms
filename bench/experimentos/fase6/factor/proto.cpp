// proto.cpp — Fase 6, tema «factor»: prototipos FUERA de la librería para
// factor() con q pequeño y en extensiones, y para el producto en F_q.
//
// Prototipos (cada uno es un tipo de cuerpo DERIVADO de los de la librería, de
// modo que Poly<...> de la librería, sin tocarla, los usa en factor() y mul):
//
//   A) FqPeq: F_q = F_p[t]/(g) con p pequeño (2·e·(p-1)^2 < 2^62).  Mismos
//      elementos (std::vector<uint64_t> de longitud e, coeficientes en [0, p))
//      y misma secuencia aleatoria (random() es el de FqT), pero:
//        * mul: convolución en acumuladores de 64 bits (sin u128 ni
//          comprobación de desbordamiento), reducción módulo g en 64 bits y
//          una sola reducción módulo p (Barrett de 64 bits) por coeficiente;
//          búferes en la pila (sin memoria dinámica salvo el resultado);
//        * acumulación perezosa (Acc) en 64 bits con un contador que
//          normaliza antes de poder desbordar (nunca, con p = 2, 3);
//        * inv: el mismo Euclides en F_p[t], en sitio sobre arrays de la pila.
//   B) mul_extension rápida: especialización explícita de ks::mul_extension
//      para los tipos *KS: mismo empaquetado y mismo mpn_mul, pero el
//      desempaquetado lee cada hueco (<= 63 bits) con un desplazamiento y una
//      máscara y lo deja SIN reducir en el acumulador (que ya reduce módulo p
//      al final), en lugar de leer() + mpn_mod_1 por hueco.
//   C) mul_primo rápida para F_p (Fp64R): especialización de ks::mul_primo con
//      desempaquetado directo y Barrett de 64 bits en lugar de u128 % p.
//   D) pow_mod con ventana deslizante (k = 4) frente al pow_mod binario de la
//      librería (mismo resultado; menos productos modulares).
//
// Todas las variantes se comparan con la LIBRERÍA (Poly<Fq64>, Poly<Fp64>):
// mismos factores, mismas multiplicidades, mismo coeficiente principal y mismo
// estado del generador aleatorio después de factor().  Entradas: las 16 de
// bench/results/fase6/factor_entradas.txt (generadas por flint_ref.py, las
// mismas que mide python-flint).  Mejor de 3 tandas de >= 30 ms rotando las 16.
//
// Compilación (desde la raíz del repo, portable como el módulo):
//   g++ -std=c++23 -O3 -Iinclude bench/experimentos/fase6/factor/proto.cpp -lgmpxx -lgmp -static -o build/fase6/factor/proto.exe
// Ejecución:  ALCP_THREADS=1 build/fase6/factor/proto.exe
#include "alcp/alcp.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace alcp;
using reloj = std::chrono::steady_clock;
using u64 = std::uint64_t;
using u128 = unsigned __int128;

static volatile u64 sumidero = 0;

template <class Fn>
double medir_us(Fn&& fn, int nent = 16) {
    double mejor = 1e300;
    for (int tanda = 0; tanda < 3; ++tanda) {
        int veces = 0;
        auto t0 = reloj::now();
        double dt = 0;
        do {
            fn(veces % nent);
            ++veces;
            dt = std::chrono::duration<double, std::micro>(reloj::now() - t0).count();
        } while (dt < 30000);
        mejor = std::min(mejor, dt / veces);
    }
    return mejor;
}

// ---------------------------------------------------------------------------
// Barrett de 64 bits: x mod p para cualquier x < 2^64 (p >= 2).
// m = floor((2^64-1)/p) da un cociente con error <= 1: r < 2p, una corrección.
struct Barrett64 {
    u64 p = 2, m = 0;
    Barrett64() = default;
    explicit Barrett64(u64 pp) : p(pp), m(~u64(0) / pp) {}
    u64 operator()(u64 x) const {
        const u64 q = static_cast<u64>((static_cast<u128>(x) * m) >> 64);
        u64 r = x - q * p;
        r -= p & (0 - static_cast<u64>(r >= p));
        return r;
    }
};

// ---------------------------------------------------------------------------
// A) FqPeq
// Cuerpo base con acumulador de 64 bits (lo usa mul_extension: acc_set).
struct BasePeq : Fp64 {
    using Fp64::Fp64;
    explicit BasePeq(const Fp64& b) : Fp64(b) {}
    using Acc = u64;
    void acc_clear(Acc& s) const { s = 0; }
    void acc_set(Acc& s, const Elem& a) const { s = a; }
};

// Acumulador de FqPeq: la convolución sin reducir en 64 bits y cuántos
// productos lleva desde la última normalización.
struct AccPeq {
    std::vector<u64> v;
    u64 cnt = 0;
    u64& operator[](std::size_t i) { return v[i]; }
    const u64& operator[](std::size_t i) const { return v[i]; }
};

struct FqPeq : Fq64 {
    static constexpr std::size_t kMaxE = 64;  // búferes en la pila
    using Acc = AccPeq;

    FqPeq(Fp64 base, std::vector<BigInt> g) : Fq64(base, g), bp_(base) {
        p_ = base.modulus64();
        e_ = extension_degree();
        red_ = Barrett64(p_);
        for (unsigned long t = 0; t < e_; ++t) {
            const u64 gt = base.from_int(g[t]);  // g ya es mónico en este uso
            ng_[t] = gt == 0 ? 0 : p_ - gt;
            g_[t] = gt;
        }
        g_[e_] = 1;
        const u128 prod = static_cast<u128>(e_) * (p_ - 1) * (p_ - 1) + 1;
        if (e_ > kMaxE || 2 * prod >= (u128(1) << 62))
            throw std::invalid_argument("FqPeq: p o e demasiado grandes para el prototipo");
        lim_ = static_cast<u64>((u128(1) << 62) / prod) - 2;
    }

    const BasePeq& base_field() const { return bp_; }

    Elem mul(const Elem& a, const Elem& b) const {
        std::array<u64, 2 * kMaxE> s{};
        conv(s.data(), a.data(), b.data());
        return reducir(s.data());
    }
    Elem div(const Elem& a, const Elem& b) const { return mul(a, inv(b)); }

    Elem inv(const Elem& a) const {
        if (is_zero(a)) throw division_by_zero("0 is not invertible in F_q");
        const std::size_t e = e_;
        std::array<u64, kMaxE + 1> A{}, B{}, S0{}, S1{};
        u64 *r0 = A.data(), *r1 = B.data(), *s0 = S0.data(), *s1 = S1.data();
        for (std::size_t i = 0; i <= e; ++i) r0[i] = g_[i];
        for (std::size_t i = 0; i < e; ++i) r1[i] = a[i];
        int d0 = static_cast<int>(e), d1 = static_cast<int>(e) - 1;
        while (r1[d1] == 0) --d1;
        int l0 = 0, l1 = 1;  // longitudes de s0 (= 0) y s1 (= 1)
        s1[0] = 1;
        while (d1 > 0) {
            const u64 il = bp_.inv(r1[d1]);
            while (d0 >= d1) {
                const u64 f = red_(r0[d0] * il);
                const u64 nf = p_ - f;  // f != 0: r0[d0] != 0
                const int pos = d0 - d1;
                for (int j = 0; j <= d1; ++j) r0[pos + j] = red_(r0[pos + j] + nf * r1[j]);
                for (int j = 0; j < l1; ++j) s0[pos + j] = red_(s0[pos + j] + nf * s1[j]);
                l0 = std::max(l0, pos + l1);
                while (d0 >= 0 && r0[d0] == 0) --d0;
            }
            while (l0 > 0 && s0[l0 - 1] == 0) --l0;
            std::swap(r0, r1);
            std::swap(d0, d1);
            std::swap(s0, s1);
            std::swap(l0, l1);
        }
        const u64 ci = bp_.inv(r1[0]);
        Elem r(e, 0);
        for (int i = 0; i < l1; ++i) r[i] = red_(s1[i] * ci);
        return r;
    }

    // --- acumulación perezosa ---
    void acc_clear(Acc& s) const {
        s.v.assign(2 * e_ - 1, 0);
        s.cnt = 0;
    }
    void acc_set(Acc& s, const Elem& a) const {
        acc_clear(s);
        for (std::size_t i = 0; i < e_; ++i) s.v[i] = a[i];
    }
    void acc_fma(Acc& s, const Elem& a, const Elem& b) const {
        if (++s.cnt >= lim_) {  // nunca con p pequeño; protege la cota
            for (auto& x : s.v) x = red_(x);
            s.cnt = 1;
        }
        conv(s.v.data(), a.data(), b.data());
    }
    Elem acc_reduce(Acc& s) const { return reducir(s.v.data()); }

protected:  // (FqSwar, en proto2.cpp, reutiliza conv y reducir)
    BasePeq bp_;
    u64 p_ = 2, lim_ = 1;
    std::size_t e_ = 1;
    Barrett64 red_;
    std::array<u64, kMaxE + 1> ng_{}, g_{};

    // s[0..2e-1) += a·b (convolución sin reducir).
    void conv(u64* __restrict s, const u64* __restrict a, const u64* __restrict b) const {
        const std::size_t e = e_;
        for (std::size_t i = 0; i < e; ++i) {
            const u64 ai = a[i];
            u64* __restrict si = s + i;
            for (std::size_t j = 0; j < e; ++j) si[j] += ai * b[j];
        }
    }
    // Reduce s (2e-1 huecos) módulo g y módulo p; consume s.
    Elem reducir(u64* s) const {
        const std::size_t e = e_;
        for (std::size_t k = 2 * e - 1; k-- > e;) {
            const u64 c = red_(s[k]);
            u64* __restrict d = s + (k - e);
            for (std::size_t t = 0; t < e; ++t) d[t] += c * ng_[t];
        }
        Elem r(e);
        for (std::size_t i = 0; i < e; ++i) r[i] = red_(s[i]);
        return r;
    }
};

// Variantes con la mul_extension rápida (B), por especialización explícita.
struct Fq64KS : Fq64 {
    using Fq64::Fq64;
};
struct FqPeqKS : FqPeq {
    using FqPeq::FqPeq;
};
// C) F_p con la mul_primo rápida.
struct Fp64R : Fp64 {
    using Fp64::Fp64;
};

// ---------------------------------------------------------------------------
// B) mul_extension rápida
namespace proto {
using ks::limb;

inline void poner(limb* X, std::size_t bit, u64 v) {
    const std::size_t w = bit >> 6, s = bit & 63;
    X[w] |= v << s;
    if (s) X[w + 1] |= v >> (64 - s);  // v < 2^63: si s + bits <= 64 esto es 0
}
inline u64 tomar(const limb* C, std::size_t bit, u64 mask) {
    const std::size_t w = bit >> 6, s = bit & 63;
    u64 x = C[w] >> s;
    if (s) x |= C[w + 1] << (64 - s);
    return x & mask;
}

// Mismo esquema que ks::mul_extension (bloques de 2e-1 huecos de `bits`).
// Desempaquetado directo, sin reducir módulo p (el acumulador lo hace).
template <class Field>
std::vector<typename Field::Elem> mul_extension_rapida(const Field& K, const typename Field::Elem* a,
                                                        std::size_t na, const typename Field::Elem* b,
                                                        std::size_t nb) {
    using Elem = typename Field::Elem;
    const std::size_t e = K.extension_degree();
    const std::size_t bloque = 2 * e - 1;
    const BigInt& p = K.characteristic();
    const std::size_t m = std::min(na, nb) * e;
    const std::size_t bits = 2 * ks::detail::bitlen(p - 1) + ks::detail::bitlen(BigInt(static_cast<unsigned long>(m))) + 1;
    if (bits > 60) return ks::mul_extension<Fq64>(static_cast<const Fq64&>(K), a, na, b, nb);  // fuera del prototipo
    const std::size_t L = na + nb - 1;
    auto empaquetar = [&](const Elem* x, std::size_t nx, std::vector<limb>& X) {
        X.assign((nx * bloque * bits + 63) / 64 + 2, 0);
        for (std::size_t i = 0; i < nx; ++i)
            for (std::size_t j = 0; j < e; ++j) poner(X.data(), (i * bloque + j) * bits, x[i][j]);
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
    std::vector<limb> C(na_l + nb_l + 2, 0);
    if (cuadrado)
        mpn_sqr(C.data(), A.data(), static_cast<mp_size_t>(na_l));
    else if (na_l >= nb_l)
        mpn_mul(C.data(), A.data(), static_cast<mp_size_t>(na_l), Bv.data(), static_cast<mp_size_t>(nb_l));
    else
        mpn_mul(C.data(), Bv.data(), static_cast<mp_size_t>(nb_l), A.data(), static_cast<mp_size_t>(na_l));
    // Los huecos del último bloque pueden pasar del final de C: se amplía con ceros.
    const std::size_t necesita = (L * bloque * bits + 63) / 64 + 2;
    if (C.size() < necesita) C.resize(necesita, 0);

    const u64 mask = (bits >= 64) ? ~u64(0) : ((u64(1) << bits) - 1);
    typename Field::Acc acc;
    std::vector<Elem> r;
    r.reserve(L);
    for (std::size_t k = 0; k < L; ++k) {
        K.acc_clear(acc);
        const std::size_t b0 = k * bloque * bits;
        for (std::size_t j = 0; j < bloque; ++j) acc[j] = tomar(C.data(), b0 + j * bits, mask);
        r.push_back(K.acc_reduce(acc));
    }
    return r;
}

// C) mul_primo con desempaquetado directo y Barrett de 64 bits (huecos <= 63 bits).
template <class Field>
std::vector<u64> mul_primo_rapida(const Field& K, const u64* a, std::size_t na, const u64* b, std::size_t nb) {
    const BigInt& p = K.characteristic();
    const std::size_t m = std::min(na, nb);
    const std::size_t bits = 2 * ks::detail::bitlen(p - 1) + ks::detail::bitlen(BigInt(static_cast<unsigned long>(m))) + 1;
    if (bits > 63) return ks::mul_primo<Fp64>(static_cast<const Fp64&>(K), a, na, b, nb);
    const std::size_t L = na + nb - 1;
    const std::size_t la = (na * bits + 63) / 64 + 2, lb = (nb * bits + 63) / 64 + 2;
    std::vector<limb> A(la, 0), B(lb, 0);
    for (std::size_t i = 0; i < na; ++i) poner(A.data(), i * bits, a[i]);
    const bool cuadrado = (a == b && na == nb);
    if (!cuadrado)
        for (std::size_t i = 0; i < nb; ++i) poner(B.data(), i * bits, b[i]);
    std::size_t na_l = la, nb_l = cuadrado ? la : lb;
    while (na_l > 1 && A[na_l - 1] == 0) --na_l;
    const limb* Bp = cuadrado ? A.data() : B.data();
    while (nb_l > 1 && Bp[nb_l - 1] == 0) --nb_l;
    std::vector<limb> C(std::max(na_l + nb_l, (L * bits + 63) / 64) + 2, 0);
    if (cuadrado)
        mpn_sqr(C.data(), A.data(), static_cast<mp_size_t>(na_l));
    else if (na_l >= nb_l)
        mpn_mul(C.data(), A.data(), static_cast<mp_size_t>(na_l), Bp, static_cast<mp_size_t>(nb_l));
    else
        mpn_mul(C.data(), Bp, static_cast<mp_size_t>(nb_l), A.data(), static_cast<mp_size_t>(na_l));
    const Barrett64 red(K.modulus64());
    const u64 mask = (u64(1) << bits) - 1;
    std::vector<u64> r(L);
    for (std::size_t k = 0; k < L; ++k) r[k] = red(tomar(C.data(), k * bits, mask));
    return r;
}
}  // namespace proto

// Especializaciones explícitas: Poly<Fq64KS>, Poly<FqPeqKS> y Poly<Fp64R> de la
// librería llaman a estas versiones (antes de cualquier instanciación).
namespace alcp::ks {
template <>
std::vector<Fq64KS::Elem> mul_extension<Fq64KS>(const Fq64KS& K, const Fq64KS::Elem* a, std::size_t na,
                                                const Fq64KS::Elem* b, std::size_t nb) {
    return proto::mul_extension_rapida(K, a, na, b, nb);
}
template <>
std::vector<FqPeqKS::Elem> mul_extension<FqPeqKS>(const FqPeqKS& K, const FqPeqKS::Elem* a, std::size_t na,
                                                  const FqPeqKS::Elem* b, std::size_t nb) {
    return proto::mul_extension_rapida(K, a, na, b, nb);
}
template <>
std::vector<u64> mul_primo<Fp64R>(const Fp64R& K, const u64* a, std::size_t na, const u64* b, std::size_t nb) {
    return proto::mul_primo_rapida(K, a, na, b, nb);
}
}  // namespace alcp::ks

// ---------------------------------------------------------------------------
// D) pow_mod con ventana deslizante de k bits (mismo resultado que pow_mod).
template <class P>
P pow_mod_ventana(P a, const BigInt& k, const typename P::Reducer& red, int w = 4) {
    const auto& K = a.field();
    const P& m = red.modulus();
    a = a % m;
    if (k == 0) return red.reduce(P::one(K));
    // Potencias impares a^1, a^3, ..., a^(2^w - 1).
    std::vector<P> imp(std::size_t(1) << (w - 1), P(K));
    imp[0] = a;
    const P a2 = red.reduce(a * a);
    for (std::size_t i = 1; i < imp.size(); ++i) imp[i] = red.reduce(imp[i - 1] * a2);
    long i = static_cast<long>(mpz_sizeinbase(k.get_mpz_t(), 2)) - 1;
    P r(K);
    bool hay = false;
    while (i >= 0) {
        if (!mpz_tstbit(k.get_mpz_t(), i)) {
            if (hay) r = red.reduce(r * r);
            --i;
            continue;
        }
        long j = std::max(0L, i - w + 1);
        while (!mpz_tstbit(k.get_mpz_t(), j)) ++j;
        unsigned v = 0;
        for (long t = i; t >= j; --t) v = (v << 1) | static_cast<unsigned>(mpz_tstbit(k.get_mpz_t(), t));
        if (hay)
            for (long t = i; t >= j; --t) r = red.reduce(r * r);
        r = hay ? red.reduce(r * imp[v >> 1]) : imp[v >> 1];
        hay = true;
        i = j - 1;
    }
    return r;
}

// ---------------------------------------------------------------------------
// Entradas y comparación
static std::map<std::string, std::vector<std::vector<BigInt>>> leer_entradas(const char* ruta) {
    std::map<std::string, std::vector<std::vector<BigInt>>> m;
    std::ifstream in(ruta);
    std::string linea;
    while (std::getline(in, linea)) {
        std::istringstream ss(linea);
        std::string caso, idx, c;
        if (!(ss >> caso >> idx)) continue;
        std::vector<BigInt> v;
        while (ss >> c) v.emplace_back(c);
        m[caso].push_back(std::move(v));
    }
    return m;
}

template <class K>
Poly<K> poly_de(const K& F, const std::vector<BigInt>& idx) {
    std::vector<typename K::Elem> v;
    for (const auto& i : idx) {
        if constexpr (requires { F.from_index(i); })
            v.push_back(F.from_index(i));
        else
            v.push_back(F.from_int(i));
    }
    return Poly<K>(F, std::move(v));
}

// ¿Mismo resultado de factor() y mismo consumo del generador?
template <class PA, class PB>
bool mismo_factor(const PA& a, const PB& b) {
    set_seed(1234);
    auto ra = a.factor();
    const BigInt ta = rand_below(BigInt(1) << 64);
    set_seed(1234);
    auto rb = b.factor();
    const BigInt tb = rand_below(BigInt(1) << 64);
    if (ta != tb || !(ra.leading == rb.leading) || ra.factors.size() != rb.factors.size()) return false;
    for (std::size_t j = 0; j < ra.factors.size(); ++j)
        if (ra.factors[j].first.coeffs() != rb.factors[j].first.coeffs() || ra.factors[j].second != rb.factors[j].second)
            return false;
    return true;
}

static const std::vector<int> MOD_2_8 = {1, 1, 0, 1, 1, 0, 0, 0, 1};
static const std::vector<int> MOD_3_30 = {2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0,
                                          2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1};
static std::vector<BigInt> a_big(const std::vector<int>& v) {
    std::vector<BigInt> r;
    for (int x : v) r.emplace_back(x);
    return r;
}

// factor() de la librería frente a las variantes, mismas entradas.
template <class KL, class... KV>
void bench_factor(const char* titulo, const std::vector<std::vector<BigInt>>& ent, const KL& kl,
                  const std::pair<const char*, KV*>&... var) {
    std::vector<Poly<KL>> fl;
    for (const auto& c : ent) fl.push_back(poly_de(kl, c));
    const double tl = medir_us([&](int i) { sumidero += fl[i].factor().factors.size(); });
    std::printf("\n== factor %s (16 entradas)\n  librería              %10.1f us\n", titulo, tl);
    auto una = [&](const char* nombre, const auto* kv) {
        using KVt = std::remove_cv_t<std::remove_pointer_t<decltype(kv)>>;
        std::vector<Poly<KVt>> fv;
        for (const auto& c : ent) fv.push_back(poly_de(*kv, c));
        bool ok = true;
        for (std::size_t i = 0; i < fl.size(); ++i) ok = ok && mismo_factor(fl[i], fv[i]);
        const double tv = medir_us([&](int i) { sumidero += fv[i].factor().factors.size(); });
        std::printf("  %-20s %10.1f us  (%.2fx)  idéntico a la librería (factores y RNG): %s\n", nombre, tv, tl / tv,
                    ok ? "sí" : "NO");
    };
    (una(var.first, var.second), ...);
}

template <class KL, class... KV>
void bench_mul(const char* titulo, const std::vector<std::vector<BigInt>>& ent, const KL& kl,
               const std::pair<const char*, KV*>&... var) {
    std::vector<Poly<KL>> A, B;
    for (std::size_t i = 0; i + 1 < ent.size(); i += 2) {
        A.push_back(poly_de(kl, ent[i]));
        B.push_back(poly_de(kl, ent[i + 1]));
    }
    const double tl = medir_us([&](int i) { sumidero += (A[i] * B[i]).degree(); });
    std::printf("\n== producto %s (16 pares)\n  librería              %10.1f us\n", titulo, tl);
    auto una = [&](const char* nombre, const auto* kv) {
        using KVt = std::remove_cv_t<std::remove_pointer_t<decltype(kv)>>;
        std::vector<Poly<KVt>> Av, Bv;
        for (std::size_t i = 0; i + 1 < ent.size(); i += 2) {
            Av.push_back(poly_de(*kv, ent[i]));
            Bv.push_back(poly_de(*kv, ent[i + 1]));
        }
        bool ok = true;
        for (std::size_t i = 0; i < A.size(); ++i) ok = ok && (A[i] * B[i]).coeffs() == (Av[i] * Bv[i]).coeffs();
        const double tv = medir_us([&](int i) { sumidero += (Av[i] * Bv[i]).degree(); });
        std::printf("  %-20s %10.1f us  (%.2fx)  idéntico: %s\n", nombre, tv, tl / tv, ok ? "sí" : "NO");
    };
    (una(var.first, var.second), ...);
}

// Desglose de ks::mul_extension de la librería: empaquetar, mpn_mul,
// desempaquetar; y el núcleo de KS2 (dos mpn_mul de la mitad de bits) como cota.
template <class K>
void desglose_mul_extension(const char* titulo, const K& F, const std::vector<std::vector<BigInt>>& ent) {
    std::vector<Poly<K>> A, B;
    for (std::size_t i = 0; i + 1 < ent.size(); i += 2) {
        A.push_back(poly_de(F, ent[i]));
        B.push_back(poly_de(F, ent[i + 1]));
    }
    const std::size_t na = A[0].coeffs().size(), e = F.extension_degree(), bloque = 2 * e - 1;
    const std::size_t bits = 2 * ks::detail::bitlen(F.characteristic() - 1) +
                             ks::detail::bitlen(BigInt(static_cast<unsigned long>(na * e))) + 1;
    const std::size_t limbs = (na * bloque * bits + 63) / 64;
    std::vector<std::vector<mp_limb_t>> X(16), Y(16), H1(16), H2(16);
    for (int i = 0; i < 16; ++i) {
        X[i].resize(limbs);
        Y[i].resize(limbs);
        H1[i].resize(limbs / 2 + 1);
        H2[i].resize(limbs / 2 + 1);
        for (auto& w : X[i]) w = (static_cast<u64>(rand_below(BigInt(1) << 62).get_ui()) << 7) ^ (i * 0x9E3779B97F4A7C15ULL);
        for (auto& w : Y[i]) w = (static_cast<u64>(rand_below(BigInt(1) << 62).get_ui()) << 9) ^ (i * 0xC2B2AE3D27D4EB4FULL);
        for (std::size_t t = 0; t < H1[i].size(); ++t) { H1[i][t] = X[i][t] | 1; H2[i][t] = Y[i][t] | 1; }
    }
    std::vector<mp_limb_t> C(2 * limbs + 4);
    const double t_mpn = medir_us([&](int i) {
        mpn_mul(C.data(), X[i].data(), static_cast<mp_size_t>(limbs), Y[i].data(), static_cast<mp_size_t>(limbs));
        sumidero += C[limbs];
    });
    const std::size_t h = limbs / 2 + 1;
    const double t_ks2 = medir_us([&](int i) {
        mpn_mul(C.data(), H1[i].data(), static_cast<mp_size_t>(h), H2[i].data(), static_cast<mp_size_t>(h));
        sumidero += C[h];
        mpn_mul(C.data(), H2[i].data(), static_cast<mp_size_t>(h), H1[i].data(), static_cast<mp_size_t>(h));
        sumidero += C[h];
    });
    const double t_tot = medir_us([&](int i) {
        sumidero += ks::mul_extension(F, A[i].coeffs().data(), na, B[i].coeffs().data(), na).size();
    });
    std::printf("  desglose ks::mul_extension %s: total %.1f us; mpn_mul %zu x %zu limbs (huecos de %zu bits) %.1f us "
                "(%.0f%%); resto (empaquetar + desempaquetar + reducir) %.1f us\n"
                "  núcleo KS2 (2 mpn_mul de %zu limbs): %.1f us  -> ahorro máximo de KS2 en el producto entero %.1f us\n",
                titulo, t_tot, limbs, limbs, bits, t_mpn, 100 * t_mpn / t_tot, t_tot - t_mpn, h, t_ks2,
                t_mpn - t_ks2);
}

// D) pow_mod ventana frente a pow_mod binario.
template <class K>
void bench_pow(const char* titulo, const K& F, int n, const BigInt& k) {
    using P = Poly<K>;
    P m = P::random(F, n).monic();
    typename P::Reducer red(m);
    std::vector<P> a;
    for (int i = 0; i < 16; ++i) a.push_back(P::random(F, n - 1));
    bool ok = true;
    for (int i = 0; i < 16; ++i) ok = ok && P::pow_mod(a[i], k, red) == pow_mod_ventana(a[i], k, red);
    const double t1 = medir_us([&](int i) { sumidero += P::pow_mod(a[i], k, red).degree(); });
    const double t2 = medir_us([&](int i) { sumidero += pow_mod_ventana(a[i], k, red).degree(); });
    std::printf("  pow_mod %s (bits %zu, popcount %lu): binario %.1f us | ventana 4 %.1f us (%.2fx) | idéntico: %s\n",
                titulo, mpz_sizeinbase(k.get_mpz_t(), 2), mpz_popcount(k.get_mpz_t()), t1, t2, t1 / t2,
                ok ? "sí" : "NO");
}

template <class K>
void elem(const char* titulo, const K& F) {
    std::vector<typename K::Elem> a, b;
    for (int i = 0; i < 16; ++i) {
        a.push_back(F.random());
        auto y = F.random();
        if (F.is_zero(y)) y = F.one();
        b.push_back(y);
    }
    typename K::Acc acc;
    const double tm = medir_us([&](int i) { sumidero += F.mul(a[i], b[i])[0]; });
    const double ti = medir_us([&](int i) { sumidero += F.inv(b[i])[0]; });
    const double tf = medir_us([&](int i) {
        F.acc_clear(acc);
        for (int t = 0; t < 8; ++t) F.acc_fma(acc, a[(i + t) & 15], b[(i + 2 * t) & 15]);
        sumidero += F.acc_reduce(acc)[0];
    });
    std::printf("  %-24s mul %6.0f ns | inv %6.0f ns | acc (clear + 8 fma + reduce) %6.0f ns\n", titulo, tm * 1e3,
                ti * 1e3, tf * 1e3);
}

#ifndef PROTO_SIN_MAIN  // proto2.cpp incluye este fichero sin su main
int main() {
    auto ent = leer_entradas("bench/results/fase6/factor_entradas.txt");
    if (ent.empty()) {
        std::printf("No se encuentra bench/results/fase6/factor_entradas.txt (ejecute flint_ref.py)\n");
        return 1;
    }
    set_seed(6);
    Fq64 k3(Fp64(3), a_big(MOD_3_30)), k256(Fp64(2), a_big(MOD_2_8));
    Fq64KS k3ks(Fp64(3), a_big(MOD_3_30)), k256ks(Fp64(2), a_big(MOD_2_8));
    FqPeq p3(Fp64(3), a_big(MOD_3_30)), p256(Fp64(2), a_big(MOD_2_8));
    FqPeqKS p3ks(Fp64(3), a_big(MOD_3_30)), p256ks(Fp64(2), a_big(MOD_2_8));
    Fp64 f101(101);
    Fp64R f101r(101);

    std::printf("Prototipos fase 6, tema factor (un hilo, mejor de 3 tandas de >= 30 ms, 16 entradas rotando)\n");
    std::printf("\n== Elementos (16 pares rotando)\n");
    elem("F_(3^30) librería", k3);
    elem("F_(3^30) FqPeq", p3);
    elem("F_256 librería", k256);
    elem("F_256 FqPeq", p256);
    {  // comprobación exhaustiva de elementos aleatorios
        bool ok = true;
        for (int i = 0; i < 2000; ++i) {
            auto a = k3.random(), b = k3.random();
            ok = ok && k3.mul(a, b) == p3.mul(a, b);
            if (!k3.is_zero(b)) ok = ok && k3.inv(b) == p3.inv(b);
            auto c = k256.random(), d = k256.random();
            ok = ok && k256.mul(c, d) == p256.mul(c, d);
            if (!k256.is_zero(d)) ok = ok && k256.inv(d) == p256.inv(d);
        }
        std::printf("  FqPeq mul/inv == librería en 2000 pares por cuerpo: %s\n", ok ? "sí" : "NO");
    }

    bench_factor("F_(3^30) grado 5", ent["f3e30_d5"], k3, std::pair{"FqPeq", &p3}, std::pair{"FqPeq + KS rápida", &p3ks});
    bench_factor("F_256 grado 30", ent["f2e8_d30"], k256, std::pair{"FqPeq", &p256}, std::pair{"FqPeq + KS rápida", &p256ks});
    bench_factor("F_101 grado 40 (producto)", ent["fp101_d40"], f101, std::pair{"mul_primo rápida", &f101r});
    bench_factor("F_101 grado 100", ent["fp101_d100"], f101, std::pair{"mul_primo rápida", &f101r});

    bench_mul("F_(3^30) grado 200", ent["f3e30_mul_d200"], k3, std::pair{"Fq64 + KS rápida", &k3ks},
              std::pair{"FqPeq", &p3}, std::pair{"FqPeq + KS rápida", &p3ks});
    desglose_mul_extension("F_(3^30) d200", k3, ent["f3e30_mul_d200"]);
    bench_mul("F_256 grado 1000", ent["f2e8_mul_d1000"], k256, std::pair{"Fq64 + KS rápida", &k256ks},
              std::pair{"FqPeq", &p256}, std::pair{"FqPeq + KS rápida", &p256ks});
    desglose_mul_extension("F_256 d1000", k256, ent["f2e8_mul_d1000"]);

    std::printf("\n== pow_mod (ventana deslizante, k = 4)\n");
    const BigInt q3 = k3.cardinality();
    bench_pow("F_(3^30) mod grado 5, exponente q", p3, 5, q3);
    bench_pow("F_(3^30) mod grado 5, exponente (q-1)/2", p3, 5, (q3 - 1) / 2);
    bench_pow("F_101 mod grado 100, exponente 101", f101, 100, BigInt(101));
    bench_pow("F_101 mod grado 20, exponente (101^3-1)/2", f101, 20, (BigInt(101) * 101 * 101 - 1) / 2);
    std::printf("\n(sumidero %llu)\n", static_cast<unsigned long long>(sumidero));
    return 0;
}
#endif  // PROTO_SIN_MAIN
