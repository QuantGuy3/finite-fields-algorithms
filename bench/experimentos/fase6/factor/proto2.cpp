// proto2.cpp — Fase 6, tema «factor»: segunda tanda de prototipos, sobre los
// de proto.cpp (que se incluye sin su main).
//
//   E) FqSwar (F_q con p pequeño, p. ej. F_(3^30)): como FqPeq, pero
//        * producto y acumulación «SWAR»: 4 coeficientes de 16 bits por
//          palabra; la convolución hace e·(e/4+1) productos de palabra en lugar
//          de e^2, y la reducción módulo g igual (e-1 filas de e/4+1 palabras);
//        * inverso: el mismo Euclides en F_p[t] con reducción módulo p
//          perezosa (una vez por ronda en lugar de una por operación);
//        * Acc de tamaño fijo (std::array): un vector<Acc> no hace una reserva
//          de memoria por coeficiente;
//        * productos de Kronecker con la mul_extension rápida (huecos leídos
//          directamente y reducidos con la reducción escalar de FqPeq).
//   F) FqTab (F_(2^e), e <= 16, p. ej. F_256): elemento = UNA palabra (los e
//      bits de los coeficientes, c_0 en el bit más alto para que el orden
//      numérico sea el lexicográfico de std::vector y canonical_less ordene
//      igual), suma = XOR y producto/inverso con tablas de logaritmos
//      (exp/log, 2(q-1) + q entradas).  random() hace las mismas e llamadas a
//      rand_below(2) que FqT, en el mismo orden.  Poly<FqTab> es la
//      librería sin cambios (Elem aritmético: productos grandes por
//      Karatsuba, sin Kronecker).
//
// Todas se comparan con la LIBRERÍA (Poly<Fq64>): factores, multiplicidades,
// orden de salida y estado del generador aleatorio tras factor().
//
// Compilación (desde la raíz del repo, portable como el módulo):
//   g++ -std=c++23 -O3 -Iinclude bench/experimentos/fase6/factor/proto2.cpp -lgmpxx -lgmp -static -o build/fase6/factor/proto2.exe
// Ejecución:  ALCP_THREADS=1 build/fase6/factor/proto2.exe
#define PROTO_SIN_MAIN
#include "proto.cpp"

// ---------------------------------------------------------------------------
// E) FqSwar
struct AccSwar {
    static constexpr std::size_t kW = (2 * FqPeq::kMaxE) / 4 + 2;
    std::array<u64, kW> w{};
    u64 cnt = 0;
};

struct FqSwar : FqPeq {
    using Acc = AccSwar;
    static constexpr std::size_t kNB = FqPeq::kMaxE / 4 + 1;  // palabras de un elemento (+1 por el desplazamiento)

    FqSwar(Fp64 base, std::vector<BigInt> g) : FqPeq(base, g) {
        nb_ = (e_ + 3) / 4;
        ns_ = (2 * e_ - 1 + 3) / 4 + 1;
        const u64 pm = (p_ - 1) * (p_ - 1);
        if (2 * e_ * pm + p_ >= (u64(1) << 16))
            throw std::invalid_argument("FqSwar: p o e demasiado grandes para huecos de 16 bits");
        lim16_ = ((u64(1) << 16) - 1 - e_ * pm - p_) / (e_ * pm) ;
        // -g desplazado 0..3 huecos
        std::array<u64, kNB + 1> base_w{};
        for (std::size_t t = 0; t < e_; ++t) base_w[t / 4] |= ng_[t] << (16 * (t % 4));
        desplazar(base_w.data(), ngsh_);
    }

    Elem mul(const Elem& a, const Elem& b) const {
        std::array<u64, AccSwar::kW> S{};
        conv16(S.data(), a.data(), b.data());
        return reducir16(S.data());
    }
    Elem div(const Elem& a, const Elem& b) const { return mul(a, inv(b)); }

    // Euclides en F_p[t] con reducción módulo p perezosa.
    Elem inv(const Elem& a) const {
        if (is_zero(a)) throw division_by_zero("0 is not invertible in F_q");
        const int e = static_cast<int>(e_);
        std::array<u64, kMaxE + 1> A{}, B{}, S0{}, S1{};
        u64 *r0 = A.data(), *r1 = B.data(), *s0 = S0.data(), *s1 = S1.data();
        for (int i = 0; i <= e; ++i) r0[i] = g_[i];
        for (int i = 0; i < e; ++i) r1[i] = a[i];
        int d0 = e, d1 = e - 1;
        while (r1[d1] == 0) --d1;
        int l0 = 0, l1 = 1;
        s1[0] = 1;
        while (d1 > 0) {
            const u64 il = bp_.inv(r1[d1]);
            while (d0 >= d1) {
                const u64 lead = red_(r0[d0]);
                if (lead == 0) { r0[d0--] = 0; continue; }
                const u64 nf = p_ - red_(lead * il);
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
        const u64 ci = bp_.inv(r1[0]);
        Elem r(e_, 0);
        for (int i = 0; i < l1; ++i) r[i] = red_(s1[i] * ci);
        return r;
    }

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

    // Reducción escalar de 2e-1 huecos de 64 bits (para la mul_extension rápida).
    Elem reducir_huecos(u64* s) const { return reducir(s); }

private:
    std::size_t nb_ = 1, ns_ = 1;
    u64 lim16_ = 1;
    std::array<std::array<u64, kNB + 1>, 4> ngsh_{};

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
    // S (huecos de 16 bits) += a·b
    void conv16(u64* __restrict S, const u64* a, const u64* b) const {
        std::array<u64, kNB + 1> bw{};
        for (std::size_t j = 0; j < e_; ++j) bw[j / 4] |= b[j] << (16 * (j % 4));
        std::array<std::array<u64, kNB + 1>, 4> bsh;
        desplazar(bw.data(), bsh);
        const std::size_t nb = nb_;
        for (std::size_t i = 0; i < e_; ++i) {
            const u64 ai = a[i];
            const u64* __restrict B = bsh[i & 3].data();
            u64* __restrict D = S + (i >> 2);
            for (std::size_t w = 0; w <= nb; ++w) D[w] += ai * B[w];
        }
    }
    Elem reducir16(u64* S) const {
        const std::size_t e = e_, nb = nb_;
        for (std::size_t k = 2 * e - 1; k-- > e;) {
            const u64 c = red_(hueco(S, k));
            const std::size_t base = k - e;
            const u64* __restrict G = ngsh_[base & 3].data();
            u64* __restrict D = S + (base >> 2);
            for (std::size_t w = 0; w <= nb; ++w) D[w] += c * G[w];
        }
        Elem r(e);
        for (std::size_t i = 0; i < e; ++i) r[i] = red_(hueco(S, i));
        return r;
    }
};

// mul_extension rápida para FqSwar: huecos leídos directamente y reducción escalar.
namespace proto {
template <class Field>
std::vector<typename Field::Elem> mul_extension_escalar(const Field& K, const typename Field::Elem* a,
                                                         std::size_t na, const typename Field::Elem* b,
                                                         std::size_t nb) {
    using Elem = typename Field::Elem;
    const std::size_t e = K.extension_degree();
    const std::size_t bloque = 2 * e - 1;
    const std::size_t m = std::min(na, nb) * e;
    const std::size_t bits = 2 * ks::detail::bitlen(K.characteristic() - 1) +
                             ks::detail::bitlen(BigInt(static_cast<unsigned long>(m))) + 1;
    if (bits > 60) throw std::invalid_argument("prototipo: huecos de más de 60 bits");
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
    std::vector<limb> C(std::max(na_l + nb_l, (L * bloque * bits + 63) / 64) + 2, 0);
    if (cuadrado)
        mpn_sqr(C.data(), A.data(), static_cast<mp_size_t>(na_l));
    else if (na_l >= nb_l)
        mpn_mul(C.data(), A.data(), static_cast<mp_size_t>(na_l), Bv.data(), static_cast<mp_size_t>(nb_l));
    else
        mpn_mul(C.data(), Bv.data(), static_cast<mp_size_t>(nb_l), A.data(), static_cast<mp_size_t>(na_l));
    const u64 mask = (u64(1) << bits) - 1;
    std::array<u64, 2 * FqPeq::kMaxE> s{};
    std::vector<Elem> r;
    r.reserve(L);
    for (std::size_t k = 0; k < L; ++k) {
        const std::size_t b0 = k * bloque * bits;
        for (std::size_t j = 0; j < bloque; ++j) s[j] = tomar(C.data(), b0 + j * bits, mask);
        r.push_back(K.reducir_huecos(s.data()));
    }
    return r;
}
}  // namespace proto

namespace alcp::ks {
template <>
std::vector<FqSwar::Elem> mul_extension<FqSwar>(const FqSwar& K, const FqSwar::Elem* a, std::size_t na,
                                                const FqSwar::Elem* b, std::size_t nb) {
    return proto::mul_extension_escalar(K, a, na, b, nb);
}
}  // namespace alcp::ks

// ---------------------------------------------------------------------------
// F) FqTab: F_(2^e) con e <= 16, elemento = palabra, tablas de logaritmos.
class FqTab {
public:
    using Elem = u64;
    using Acc = u64;

    FqTab(const std::vector<int>& g) : base_(BigInt(2)), e_(static_cast<unsigned>(g.size() - 1)) {
        if (e_ < 1 || e_ > 16 || g.back() != 1) throw std::invalid_argument("FqTab: e en [1, 16] y g mónico");
        q_ = BigInt(1) << e_;
        const u64 q = u64(1) << e_, n = q - 1;
        u64 gstd = 0;
        for (unsigned i = 0; i <= e_; ++i) gstd |= u64(g[i] & 1) << i;
        auto mul_std = [&](u64 x, u64 y) {  // producto en F_2[t]/(g), bit i = c_i
            u64 r = 0;
            for (unsigned i = 0; i < e_; ++i)
                if ((y >> i) & 1) r ^= x << i;
            for (unsigned k = 2 * e_ - 1; k-- > e_;)
                if ((r >> k) & 1) r ^= gstd << (k - e_);
            return r;
        };
        // Generador del grupo multiplicativo (orden q - 1).
        u64 gen = 0;
        for (u64 c = 2; c < q && !gen; ++c) {
            u64 x = c, ord = 1;
            while (x != 1) { x = mul_std(x, c); ++ord; }
            if (ord == n) gen = c;
        }
        if (q == 2) gen = 1;
        exp_.assign(2 * n + 1, 0);
        log_.assign(q, 0);
        u64 x = 1;
        for (u64 k = 0; k < n; ++k) {
            exp_[k] = exp_[k + n] = rev(x);
            log_[rev(x)] = static_cast<std::uint32_t>(k);
            x = mul_std(x, gen);
        }
        exp_[2 * n] = exp_[0];
        one_ = rev(1);
        for (int c : g) g_.push_back(c);
    }

    // Codificación: bit (e-1-i) = c_i (c_0 en el bit más alto).
    u64 rev(u64 x) const {
        u64 r = 0;
        for (unsigned i = 0; i < e_; ++i) r |= ((x >> i) & 1) << (e_ - 1 - i);
        return r;
    }
    u64 pack(const std::vector<u64>& c) const {
        u64 r = 0;
        for (unsigned i = 0; i < e_; ++i) r |= (c[i] & 1) << (e_ - 1 - i);
        return r;
    }

    // --- interfaz FiniteField ---
    Elem zero() const { return 0; }
    Elem one() const { return one_; }
    Elem from_int(const BigInt& a) const { return mpz_odd_p(a.get_mpz_t()) ? one_ : 0; }
    Elem random() const {  // las mismas e llamadas a rand_below(2) que FqT::random
        u64 r = 0;
        for (unsigned i = 0; i < e_; ++i) r |= base_.random() << (e_ - 1 - i);
        return r;
    }
    Elem add(Elem a, Elem b) const { return a ^ b; }
    Elem sub(Elem a, Elem b) const { return a ^ b; }
    Elem neg(Elem a) const { return a; }
    Elem mul(Elem a, Elem b) const {
        if (a == 0 || b == 0) return 0;
        return exp_[log_[a] + log_[b]];
    }
    Elem inv(Elem a) const {
        if (a == 0) throw division_by_zero("0 is not invertible in F_q");
        return exp_[(exp_.size() - 1) / 2 - log_[a]];
    }
    Elem div(Elem a, Elem b) const { return mul(a, inv(b)); }
    Elem pow(Elem a, const BigInt& k) const {
        const u64 n = (u64(1) << e_) - 1;
        if (a == 0) {
            if (sgn(k) < 0) throw division_by_zero("0 is not invertible in F_q");
            return sgn(k) == 0 ? one_ : 0;
        }
        BigInt r = k % BigInt(static_cast<unsigned long>(n));
        if (sgn(r) < 0) r += static_cast<unsigned long>(n);
        const u64 kk = r.get_ui();
        return exp_[(static_cast<u64>(log_[a]) * kk) % n];
    }
    Elem pth_root(Elem a) const { return pow(a, q_ / 2); }
    bool eq(Elem a, Elem b) const { return a == b; }
    bool is_zero(Elem a) const { return a == 0; }
    bool is_one(Elem a) const { return a == one_; }
    BigInt cardinality() const { return q_; }
    BigInt characteristic() const { return BigInt(2); }
    std::string to_string(Elem a) const { return std::to_string(a); }
    unsigned long fft_capacity() const { return 0; }  // q - 1 impar
    const Elem& fft_base_root() const { return one_; }
    mutable FFTCache<Elem> fft_cache_;
    bool operator==(const FqTab& o) const { return g_ == o.g_; }

    // acumulación: la suma es XOR, no hay nada que diferir
    void acc_clear(Acc& s) const { s = 0; }
    void acc_set(Acc& s, Elem a) const { s = a; }
    void acc_fma(Acc& s, Elem a, Elem b) const { s ^= mul(a, b); }
    Elem acc_reduce(Acc& s) const { return s; }

private:
    Fp64 base_;
    unsigned e_;
    BigInt q_;
    u64 one_ = 1;
    std::vector<u64> exp_;
    std::vector<std::uint32_t> log_;
    std::vector<int> g_;
};
static_assert(FiniteField<FqTab>);

// ---------------------------------------------------------------------------
static bool mismo_factor_tab(const Poly<Fq64>& a, const Poly<FqTab>& b, const FqTab& T) {
    set_seed(4321);
    auto ra = a.factor();
    const BigInt ta = rand_below(BigInt(1) << 64);
    set_seed(4321);
    auto rb = b.factor();
    const BigInt tb = rand_below(BigInt(1) << 64);
    if (ta != tb || T.pack(ra.leading) != rb.leading || ra.factors.size() != rb.factors.size()) return false;
    for (std::size_t j = 0; j < ra.factors.size(); ++j) {
        const auto& ca = ra.factors[j].first.coeffs();
        const auto& cb = rb.factors[j].first.coeffs();
        if (ca.size() != cb.size() || ra.factors[j].second != rb.factors[j].second) return false;
        for (std::size_t i = 0; i < ca.size(); ++i)
            if (T.pack(ca[i]) != cb[i]) return false;
    }
    return true;
}

int main() {
    auto ent = leer_entradas("bench/results/fase6/factor_entradas.txt");
    if (ent.empty()) {
        std::printf("No se encuentra bench/results/fase6/factor_entradas.txt (ejecute flint_ref.py)\n");
        return 1;
    }
    set_seed(6);
    Fq64 k3(Fp64(3), a_big(MOD_3_30)), k256(Fp64(2), a_big(MOD_2_8));
    FqPeq p3(Fp64(3), a_big(MOD_3_30));
    FqSwar s3(Fp64(3), a_big(MOD_3_30));
    FqTab t256(MOD_2_8);

    std::printf("Prototipos fase 6 (2), tema factor (un hilo, mejor de 3 tandas de >= 30 ms, 16 entradas rotando)\n");
    std::printf("\n== Elementos (16 pares rotando)\n");
    elem("F_(3^30) librería", k3);
    elem("F_(3^30) FqPeq", p3);
    elem("F_(3^30) FqSwar", s3);
    {
        bool ok = true;
        for (int i = 0; i < 3000 && ok; ++i) {
            auto a = k3.random(), b = k3.random();
            ok = ok && k3.mul(a, b) == s3.mul(a, b);
            if (!k3.is_zero(b)) ok = ok && k3.inv(b) == s3.inv(b);
            FqSwar::Acc acc;
            Fq64::Acc accl;
            s3.acc_set(acc, a);
            k3.acc_set(accl, a);
            for (int t = 0; t < 40; ++t) {
                auto x = k3.random(), y = k3.random();
                s3.acc_fma(acc, x, y);
                k3.acc_fma(accl, x, y);
            }
            ok = ok && s3.acc_reduce(acc) == k3.acc_reduce(accl);
        }
        std::printf("  FqSwar mul/inv/acc == librería en 3000 casos: %s\n", ok ? "sí" : "NO");
    }
    {  // F_256 por tablas
        std::vector<Fq64::Elem> a, b;
        std::vector<u64> ta, tb;
        for (int i = 0; i < 16; ++i) {
            a.push_back(k256.random());
            auto y = k256.random();
            if (k256.is_zero(y)) y = k256.one();
            b.push_back(y);
            ta.push_back(t256.pack(a.back()));
            tb.push_back(t256.pack(b.back()));
        }
        const double tm = medir_us([&](int i) { sumidero += t256.mul(ta[i], tb[i]); });
        const double ti = medir_us([&](int i) { sumidero += t256.inv(tb[i]); });
        std::printf("  %-24s mul %6.1f ns | inv %6.1f ns\n", "F_256 FqTab", tm * 1e3, ti * 1e3);
        bool ok = true;
        for (int i = 0; i < 256; ++i)
            for (int j = 0; j < 256; ++j) {
                auto x = k256.from_index(BigInt(i)), y = k256.from_index(BigInt(j));
                ok = ok && t256.pack(k256.mul(x, y)) == t256.mul(t256.pack(x), t256.pack(y));
                ok = ok && t256.pack(k256.add(x, y)) == t256.add(t256.pack(x), t256.pack(y));
                if (j) ok = ok && t256.pack(k256.inv(y)) == t256.inv(t256.pack(y));
                // orden: lexicográfico de vector == numérico de la palabra
                ok = ok && ((x < y) == (t256.pack(x) < t256.pack(y)));
            }
        std::printf("  FqTab == librería en las 65 536 parejas (mul, add, inv, orden): %s\n", ok ? "sí" : "NO");
    }

    bench_factor("F_(3^30) grado 5", ent["f3e30_d5"], k3, std::pair{"FqPeq", &p3}, std::pair{"FqSwar", &s3});
    {
        const auto& e30 = ent["f2e8_d30"];
        std::vector<Poly<Fq64>> fl;
        std::vector<Poly<FqTab>> ft;
        bool ok = true;
        for (const auto& c : e30) {
            fl.push_back(poly_de(k256, c));
            std::vector<u64> v;
            for (const auto& x : fl.back().coeffs()) v.push_back(t256.pack(x));
            ft.emplace_back(t256, std::move(v));
            ok = ok && mismo_factor_tab(fl.back(), ft.back(), t256);
        }
        const double tl = medir_us([&](int i) { sumidero += fl[i].factor().factors.size(); });
        const double tt = medir_us([&](int i) { sumidero += ft[i].factor().factors.size(); });
        std::printf("\n== factor F_256 grado 30 (16 entradas)\n  librería              %10.1f us\n"
                    "  FqTab                %10.1f us  (%.2fx)  idéntico a la librería (factores, orden y RNG): %s\n",
                    tl, tt, tl / tt, ok ? "sí" : "NO");
    }
    bench_mul("F_(3^30) grado 200", ent["f3e30_mul_d200"], k3, std::pair{"FqSwar (+ KS rápida)", &s3});
    {
        const auto& m = ent["f2e8_mul_d1000"];
        std::vector<Poly<Fq64>> A, B;
        std::vector<Poly<FqTab>> TA, TB;
        bool ok = true;
        for (std::size_t i = 0; i + 1 < m.size(); i += 2) {
            A.push_back(poly_de(k256, m[i]));
            B.push_back(poly_de(k256, m[i + 1]));
            std::vector<u64> va, vb;
            for (const auto& x : A.back().coeffs()) va.push_back(t256.pack(x));
            for (const auto& x : B.back().coeffs()) vb.push_back(t256.pack(x));
            TA.emplace_back(t256, std::move(va));
            TB.emplace_back(t256, std::move(vb));
        }
        for (std::size_t i = 0; i < A.size(); ++i) {
            auto c = (A[i] * B[i]).coeffs();
            auto d = (TA[i] * TB[i]).coeffs();
            ok = ok && c.size() == d.size();
            for (std::size_t j = 0; ok && j < c.size(); ++j) ok = t256.pack(c[j]) == d[j];
        }
        const double tl = medir_us([&](int i) { sumidero += (A[i] * B[i]).degree(); });
        const double tt = medir_us([&](int i) { sumidero += (TA[i] * TB[i]).degree(); });
        std::printf("\n== producto F_256 grado 1000 (16 pares)\n  librería              %10.1f us\n"
                    "  FqTab (Karatsuba)    %10.1f us  (%.2fx)  idéntico: %s\n",
                    tl, tt, tl / tt, ok ? "sí" : "NO");
    }
    std::printf("\n(sumidero %llu)\n", static_cast<unsigned long long>(sumidero));
    return 0;
}
