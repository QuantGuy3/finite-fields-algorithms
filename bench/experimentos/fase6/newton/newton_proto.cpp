// newton_proto.cpp — Fase 6, tema «newton»: perfil y prototipos de la división de
// Newton, del Reducer y de pow_mod en F_998244353 (Fp64, NTT en palabras).
//
// PROTOTIPO FUERA DE LA LIBRERÍA: no modifica include/.  Usa los internos públicos
// (fft::tables, fft::word::transform, Poly::inverse_series, Poly::Reducer,
// Poly::Composer...) y compara SIEMPRE con las funciones de la librería.
//
// Qué mide (un hilo; poner ALCP_THREADS=1):
//   1) Perfil: coste de una transformada por tamaño; divmod 4000/2000 con el
//      cuerpo «frío» (sin tablas FFT: así lo mide la comparativa desde Python,
//      porque divmod solo usa mul_no_rng y nunca pide la raíz) y «caliente»
//      (tablas ya construidas), con su reparto inversión / cociente / resto;
//      pow_mod(a, q, m) con deg m = 300 y 1000 y su reparto cuadrado / producto /
//      reducción.
//   2) Prototipos (todos con las tablas YA CONSTRUIDAS; en la librería harían
//      falta tablas sin RNG para el caso frío, ver el informe):
//      a) inversión de Newton con producto medio (transformadas de tamaño n en
//         lugar de 2^ceil(log2(n + n/2))) y resto f - g·q calculado módulo
//         x^N - 1 con N >= deg g (la mitad de tamaño);
//      b) Reducer con las transformadas de rev(m)^-1 y de m precalculadas
//         (FFTRep de NTL): 2 transformadas de tamaño 2n + 2 de tamaño n por
//         reducción, en lugar de 6 de tamaño 2n;
//      c) cuadrado con una sola transformada directa (la librería YA lo hace si
//         los dos operandos son el mismo objeto: se comprueba midiendo a*a frente
//         a a*b);
//      d) Composer (Brent-Kung, Frobenius de is_irreducible / DDF de factor) con
//         la transformada del paso gigante y de h precalculadas + Reducer b.
//   Cada variante se comprueba contra la librería (resultado idéntico).
//
// Metodología: 16 entradas distintas rotando, mejor de 3 tandas de >= 30 ms,
// checksum acumulado para que el compilador no elimine el trabajo.
//
// Compilación (portable, como el módulo):
//   g++ -std=c++23 -O3 -Iinclude bench/experimentos/fase6/newton/newton_proto.cpp -lgmpxx -lgmp -static -o build/fase6/newton/newton_proto.exe
// Ejecución (desde la raíz del repo):
//   ALCP_THREADS=1 build/fase6/newton/newton_proto.exe
#include "alcp/alcp.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace alcp;
using P = Poly<Fp64>;
using u64 = std::uint64_t;
using u128 = unsigned __int128;
using vec = std::vector<u64>;
using Tab = FFTTables<u64>;
using reloj = std::chrono::steady_clock;

static u64 g_chk = 0;  // checksum global (evita que se elimine el trabajo)
static void usa(const vec& v) { g_chk = g_chk * 1000003u + v.size() + (v.empty() ? 0 : v[0] + v.back()); }
static void usa(const P& v) { usa(v.coeffs()); }

// Mejor de 3 tandas de >= 30 ms, rotando 16 entradas (fn recibe el índice).
template <class Fn>
double medir_us(Fn&& fn) {
    double mejor = 1e300;
    for (int tanda = 0; tanda < 3; ++tanda) {
        int veces = 0;
        auto t0 = reloj::now();
        double dt = 0;
        do {
            fn(veces % 16);
            ++veces;
            dt = std::chrono::duration<double, std::micro>(reloj::now() - t0).count();
        } while (dt < 30000);
        mejor = std::min(mejor, dt / veces);
    }
    return mejor;
}

static int g_fallos = 0;
static void comprueba(bool ok, const char* que) {
    if (!ok) { ++g_fallos; std::printf("  *** DIFERENTE: %s\n", que); }
}

// ---------------------------------------------------------------------------
// Primitivas en palabras
// ---------------------------------------------------------------------------
struct Ctx {
    const Fp64& K;
    u64 p;
    u64 mu;            // floor(2^64 / p), para Barrett (p < 2^32)
    bool barrett;      // producto punto a punto: Barrett (true) o '%' como la librería
    explicit Ctx(const Fp64& k, bool b) : K(k), p(k.modulus64()), mu(~u64(0) / k.modulus64()), barrett(b) {}
    // Las tablas se suponen YA construidas (el cuerpo ya tiene raíz): no se
    // consume el generador aquí.
    const Tab& tab(unsigned k) const { return fft::tables(K, k); }
    u64 mul(u64 a, u64 b) const {
        if (!barrett) return fft::word::mul_mod(a, b, p);
        const u64 x = a * b;  // < p^2 < 2^60
        const u64 q = static_cast<u64>((static_cast<u128>(x) * mu) >> 64);
        u64 r = x - q * p;     // en [0, 2p)
        r -= p & (0 - static_cast<u64>(r >= p));
        return r;
    }
    u64 sub(u64 a, u64 b) const { u64 r = a - b; return r + (p & (0 - static_cast<u64>(a < b))); }
    u64 add(u64 a, u64 b) const { u64 r = a + b; return r - (p & (0 - static_cast<u64>(r >= p))); }
};

static void fwd(const Tab& t, vec& a) { fft::word::transform(t, a.data(), a.size(), false); }
static void inv(const Tab& t, vec& a) { fft::word::transform(t, a.data(), a.size(), true); }
static void trim(vec& v) { while (!v.empty() && v.back() == 0) v.pop_back(); }
static unsigned lg2(std::size_t n) { return static_cast<unsigned>(fft::ceil_log2(n)); }

// Producto punto a punto por un operando FIJO con cociente de Shoup precalculado.
struct Fijo {
    vec w, ws;  // transformada y sus cocientes de Shoup
    void prepara(const Ctx& C, vec t) {
        w = std::move(t);
        ws.resize(w.size());
        for (std::size_t i = 0; i < w.size(); ++i) ws[i] = fft::word::shoup(w[i], C.p);
    }
    void aplica(const Ctx& C, vec& a) const {
        for (std::size_t i = 0; i < a.size(); ++i) a[i] = fft::word::mul_shoup(w[i], ws[i], a[i], C.p);
    }
};

// Copia la transformada de v (relleno con ceros hasta N, plegado módulo x^N - 1).
static vec fold(const Ctx& C, const vec& v, std::size_t N) {
    vec r(N, 0);
    for (std::size_t i = 0; i < v.size(); ++i) r[i % N] = (i < N) ? v[i] : C.add(r[i % N], v[i]);
    return r;
}

// ---------------------------------------------------------------------------
// a) Inversión de la serie con producto medio
// ---------------------------------------------------------------------------
// x_n = x_k - x_k·((a·x_k) div x^k)·x^k mod x^n, k = ceil(n/2) (misma recursión
// que Poly::inverse_series).  (a mod x^n)·x_k tiene n + k - 1 coeficientes; con
// una convolución CÍCLICA de tamaño N >= n el solapamiento cae en los índices
// < k - 1, y los coeficientes k..n-1 (los únicos que hacen falta) salen limpios:
// es el producto medio.  El segundo producto reutiliza la transformada de x_k.
// 5 transformadas de tamaño 2^ceil(log2 n) por paso, frente a 6 de tamaño
// 2^ceil(log2(n + k - 1)) en la librería.  El inverso módulo x^n es único: el
// resultado es idéntico.
static vec inv_serie_mp(const Ctx& C, const vec& a, std::size_t n) {
    if (n <= 32) {  // caso base: la propia librería (escolar a este tamaño)
        vec h(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(std::min(n, a.size())));
        h.resize(n, 0);
        return P::inverse_series(C.K, h, n);
    }
    const std::size_t k = (n + 1) / 2;
    vec x = inv_serie_mp(C, a, k);
    const unsigned l = lg2(n);
    const std::size_t N = std::size_t(1) << l;
    const Tab& t = C.tab(l);
    vec fa(N, 0), fx(N, 0);
    std::copy(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(std::min(n, a.size())), fa.begin());
    std::copy(x.begin(), x.end(), fx.begin());
    fwd(t, fa);
    fwd(t, fx);
    for (std::size_t i = 0; i < N; ++i) fa[i] = C.mul(fa[i], fx[i]);
    inv(t, fa);
    vec e(N, 0);  // e = coeficientes k..n-1 de a·x_k
    for (std::size_t i = k; i < n; ++i) e[i - k] = fa[i];
    fwd(t, e);
    for (std::size_t i = 0; i < N; ++i) e[i] = C.mul(e[i], fx[i]);
    inv(t, e);  // x_k·e: n - 1 coeficientes <= N, sin solapamiento
    x.resize(n);
    for (std::size_t i = 0; i < n - k; ++i) x[k + i] = C.sub(0, e[i]);
    return x;
}

// divmod f / g (deg f >= deg g >= 1) con a) y el resto módulo x^N2 - 1.
static std::pair<vec, vec> divmod_mp(const Ctx& C, const vec& f, const vec& g) {
    const std::size_t m1 = f.size(), n1 = g.size(), n = n1 - 1, d1 = m1 - n1 + 1;
    vec grev(d1);
    for (std::size_t j = 0; j < d1; ++j) grev[j] = j < n1 ? g[n1 - 1 - j] : 0;
    const vec ginv = inv_serie_mp(C, grev, d1);
    // rev(q) = rev(f)·ginv mod x^d1: producto corto (hace falta tamaño >= 2·d1 - 1).
    const unsigned l1 = lg2(2 * d1 - 1);
    const std::size_t N1 = std::size_t(1) << l1;
    const Tab& t1 = C.tab(l1);
    vec b1(N1, 0), bg(N1, 0);
    for (std::size_t i = 0; i < d1; ++i) b1[i] = f[m1 - 1 - i];
    std::copy(ginv.begin(), ginv.end(), bg.begin());
    fwd(t1, b1);
    fwd(t1, bg);
    for (std::size_t i = 0; i < N1; ++i) b1[i] = C.mul(b1[i], bg[i]);
    inv(t1, b1);
    vec q(d1);
    for (std::size_t i = 0; i < d1; ++i) q[i] = b1[d1 - 1 - i];
    // r = f - g·q tiene grado < n: basta calcularlo módulo x^N2 - 1 con N2 >= n.
    const unsigned l2 = lg2(n);
    const std::size_t N2 = std::size_t(1) << l2;
    const Tab& t2 = C.tab(l2);
    vec x = fold(C, g, N2), y = fold(C, q, N2);
    fwd(t2, x);
    fwd(t2, y);
    for (std::size_t i = 0; i < N2; ++i) x[i] = C.mul(x[i], y[i]);
    inv(t2, x);
    const vec fp = fold(C, f, N2);
    vec r(n);
    for (std::size_t i = 0; i < n; ++i) r[i] = C.sub(fp[i], x[i]);
    trim(q);
    trim(r);
    return {std::move(q), std::move(r)};
}

// ---------------------------------------------------------------------------
// b) Reducer con transformadas precalculadas (FFTRep de NTL)
// ---------------------------------------------------------------------------
// Precalcula ginv = rev(m)^-1 mod x^n (con Poly::inverse_series, una vez), su
// transformada de tamaño N1 = 2^ceil(log2(2n - 1)) y la de m plegado módulo
// x^N2 - 1, N2 = 2^ceil(log2 n).  reduce(a), deg a <= 2n - 1:
//   rev(q) = rev(a)[0..d1)·ginv mod x^d1   (1 directa + 1 inversa de N1; usar
//            ginv completo no cambia los d1 primeros coeficientes y d1 + n - 1 <= N1)
//   r = a - m·q mod (x^N2 - 1)            (1 directa + 1 inversa de N2)
struct ReducerProto {
    const Ctx& C;
    vec m;
    std::size_t n = 0, N1 = 0, N2 = 0;
    const Tab* t1 = nullptr;
    const Tab* t2 = nullptr;
    Fijo G, M;
    ReducerProto(const Ctx& c, const P& mod) : C(c), m(mod.coeffs()) {
        n = m.size() - 1;
        vec rh(n);
        for (std::size_t j = 0; j < n; ++j) rh[j] = m[n - j];
        const vec ginv = P::inverse_series(C.K, rh, n);
        const unsigned l1 = lg2(2 * n - 1), l2 = lg2(n);
        N1 = std::size_t(1) << l1;
        N2 = std::size_t(1) << l2;
        t1 = &C.tab(l1);
        t2 = &C.tab(l2);
        vec g(N1, 0);
        std::copy(ginv.begin(), ginv.end(), g.begin());
        fwd(*t1, g);
        G.prepara(C, std::move(g));
        vec mm = fold(C, m, N2);
        fwd(*t2, mm);
        M.prepara(C, std::move(mm));
    }
    // a con coeficientes recortados, deg a <= 2n - 1.
    vec reduce(const vec& a) const {
        if (a.size() <= n) return a;
        const std::size_t da = a.size() - 1, d1 = da - n + 1;
        vec b(N1, 0);
        for (std::size_t i = 0; i < d1; ++i) b[i] = a[da - i];
        fwd(*t1, b);
        G.aplica(C, b);
        inv(*t1, b);
        vec y(N2, 0);
        for (std::size_t i = 0; i < d1; ++i) y[i] = b[d1 - 1 - i];  // q (d1 <= n <= N2)
        fwd(*t2, y);
        M.aplica(C, y);
        inv(*t2, y);
        vec r = fold(C, a, N2);
        r.resize(n);
        for (std::size_t i = 0; i < n; ++i) r[i] = C.sub(r[i], y[i]);
        trim(r);
        return r;
    }
    // Producto propio (3 transformadas de N1) y reducción.
    vec mulmod(const vec& a, const vec& b) const {
        if (a.empty() || b.empty()) return {};
        vec x(N1, 0), y(N1, 0);
        std::copy(a.begin(), a.end(), x.begin());
        std::copy(b.begin(), b.end(), y.begin());
        fwd(*t1, x);
        fwd(*t1, y);
        for (std::size_t i = 0; i < N1; ++i) x[i] = C.mul(x[i], y[i]);
        inv(*t1, x);
        x.resize(a.size() + b.size() - 1);
        return reduce(x);
    }
    // c) Cuadrado: una sola transformada directa.
    vec sqrmod(const vec& a) const {
        if (a.empty()) return {};
        vec x(N1, 0);
        std::copy(a.begin(), a.end(), x.begin());
        fwd(*t1, x);
        for (std::size_t i = 0; i < N1; ++i) x[i] = C.mul(x[i], x[i]);
        inv(*t1, x);
        x.resize(2 * a.size() - 1);
        return reduce(x);
    }
    // Producto por un operando fijo ya transformado (tamaño N1).
    vec mulmod_fijo(const vec& a, const Fijo& F) const {
        if (a.empty()) return {};
        vec x(N1, 0);
        std::copy(a.begin(), a.end(), x.begin());
        fwd(*t1, x);
        F.aplica(C, x);
        inv(*t1, x);
        trim(x);  // F no es cero, pero el producto puede tener ceros finales tras recortar a 2n - 1
        if (x.size() > 2 * n - 1) x.resize(2 * n - 1);
        trim(x);
        return reduce(x);
    }
    Fijo transformada(const vec& a) const {
        vec x(N1, 0);
        std::copy(a.begin(), a.end(), x.begin());
        fwd(*t1, x);
        Fijo F;
        F.prepara(C, std::move(x));
        return F;
    }
};

// pow_mod con el mismo orden de operaciones que Poly::pow_mod (a ya reducido).
// Variante 1: productos de la librería (Poly, regla Auto) + ReducerProto.
static vec powmod_v1(const ReducerProto& R, const P& a0, BigInt k) {
    const Fp64& K = R.C.K;
    P result = P::one(K), a = a0;
    while (k > 0) {
        if (mpz_odd_p(k.get_mpz_t())) result = P(K, R.reduce((result * a).coeffs()));
        k >>= 1;
        if (k > 0) a = P(K, R.reduce((a * a).coeffs()));
    }
    return result.coeffs();
}
// Variante 2: todo en palabras (productos y cuadrados propios + ReducerProto).
static vec powmod_v2(const ReducerProto& R, vec a, BigInt k) {
    vec result{1};
    while (k > 0) {
        if (mpz_odd_p(k.get_mpz_t())) result = R.mulmod(result, a);
        k >>= 1;
        if (k > 0) a = R.sqrmod(a);
    }
    return result;
}

// ---------------------------------------------------------------------------
// d) Composer (Brent-Kung) con transformadas precalculadas
// ---------------------------------------------------------------------------
struct ComposerProto {
    const ReducerProto& R;
    std::size_t k;
    std::vector<vec> pot;
    Fijo gig;
    ComposerProto(const ReducerProto& r, const vec& h, std::size_t kk) : R(r) {
        k = std::max<std::size_t>(1, std::min(kk, R.n));
        const Fijo H = R.transformada(h);  // h ya reducido (deg < n)
        pot.reserve(k);
        pot.push_back(vec{1});
        for (std::size_t i = 1; i < k; ++i) pot.push_back(R.mulmod_fijo(pot.back(), H));
        gig = R.transformada(R.mulmod_fijo(pot.back(), H));
    }
    vec combina(const vec& g, std::size_t ini, std::size_t cnt) const {
        const u64 p = R.C.p;
        std::size_t n = 0;
        for (std::size_t i = 0; i < cnt; ++i) n = std::max(n, pot[i].size());
        std::vector<u128> s(n, 0);
        for (std::size_t i = 0; i < cnt; ++i) {
            const u64 c = g[ini + i];
            if (c == 0) continue;
            const vec& Pt = pot[i];
            for (std::size_t t = 0; t < Pt.size(); ++t) s[t] += static_cast<u128>(c) * Pt[t];
        }
        vec r(n);
        for (std::size_t t = 0; t < n; ++t) r[t] = static_cast<u64>(s[t] % p);
        trim(r);
        return r;
    }
    vec compose(const vec& g) const {
        if (g.empty()) return {};
        const std::size_t L = g.size(), bloques = (L + k - 1) / k;
        vec acc;
        for (std::size_t j = bloques; j-- > 0;) {
            vec B = combina(g, j * k, std::min(k, L - j * k));
            if (j + 1 == bloques) { acc = std::move(B); continue; }
            acc = R.mulmod_fijo(acc, gig);
            if (acc.size() < B.size()) acc.resize(B.size(), 0);
            for (std::size_t t = 0; t < B.size(); ++t) acc[t] = R.C.add(acc[t], B[t]);
            trim(acc);
        }
        return acc;
    }
};

// ---------------------------------------------------------------------------
// Programa
// ---------------------------------------------------------------------------
static vec aleatorio(const Fp64& K, int d) { return P::random(K, d).coeffs(); }

int main(int argc, char** argv) {
    const bool barrett = argc > 1 && std::string(argv[1]) == "barrett";
    const BigInt p("998244353");
    Fp64 Kw(p);  // cuerpo «caliente»: se le construyen las tablas
    Fp64 Kc(p);  // cuerpo «frío»: nunca pide la raíz (como la comparativa desde Python)
    (void)(P::random(Kw, 300) * P::random(Kw, 300));  // tablas FFT (consume el RNG: da igual aquí)
    for (unsigned k = 1; k <= 14; ++k) (void)fft::tables(Kw, k);
    Ctx C(Kw, barrett);
    std::printf("newton_proto: F_998244353, punto a punto en prototipos: %s\n",
                barrett ? "Barrett" : "'%' (como la librería)");
    std::printf("ALCP_THREADS=%s\n", std::getenv("ALCP_THREADS") ? std::getenv("ALCP_THREADS") : "(sin definir)");

    // ---- 1) coste de una transformada
    std::printf("\n== Transformada en palabras (directa), us\n");
    for (unsigned l = 9; l <= 13; ++l) {
        const std::size_t N = std::size_t(1) << l;
        std::vector<vec> xs;
        for (int i = 0; i < 16; ++i) { vec v = aleatorio(Kw, static_cast<int>(N) - 1); v.resize(N, 0); xs.push_back(v); }
        const Tab& t = C.tab(l);
        double tt = medir_us([&](int i) { fwd(t, xs[i]); g_chk += xs[i][1]; });
        // punto a punto (lib '%') y Barrett
        double tp = medir_us([&](int i) {
            vec& a = xs[i]; const vec& b = xs[(i + 1) % 16];
            for (std::size_t j = 0; j < N; ++j) a[j] = fft::word::mul_mod(a[j], b[j], C.p);
            g_chk += a[3];
        });
        Ctx Cb(Kw, true);
        double tb = medir_us([&](int i) {
            vec& a = xs[i]; const vec& b = xs[(i + 1) % 16];
            for (std::size_t j = 0; j < N; ++j) a[j] = Cb.mul(a[j], b[j]);
            g_chk += a[3];
        });
        std::printf("  N=%5zu  transformada %8.1f   punto a punto %%: %6.1f   Barrett: %6.1f\n", N, tt, tp, tb);
    }

    // ---- 2) divmod 4000 / 2000
    {
        std::printf("\n== divmod grado 4000 / grado 2000 (us)\n");
        std::vector<vec> fv, gv;
        for (int i = 0; i < 16; ++i) { fv.push_back(aleatorio(Kw, 4000)); gv.push_back(aleatorio(Kw, 2000)); }
        std::vector<P> fw, gw, fc, gc;
        for (int i = 0; i < 16; ++i) {
            fw.emplace_back(Kw, fv[i]); gw.emplace_back(Kw, gv[i]);
            fc.emplace_back(Kc, fv[i]); gc.emplace_back(Kc, gv[i]);
        }
        const std::size_t d1 = 2001;
        auto revhead = [&](const vec& g) { vec r(d1); for (std::size_t j = 0; j < d1; ++j) r[j] = j < g.size() ? g[g.size() - 1 - j] : 0; return r; };
        std::vector<vec> ginv;
        for (int i = 0; i < 16; ++i) ginv.push_back(P::inverse_series(Kw, revhead(gv[i]), d1));

        double frio = medir_us([&](int i) { auto r = fc[i].divmod(gc[i]); usa(r.second); });
        comprueba(Kc.fft_cache_.empty(), "el cuerpo frío no debe tener tablas");
        double frio_inv = medir_us([&](int i) { usa(P::inverse_series(Kc, revhead(gv[i]), d1)); });
        double frio_resto = medir_us([&](int i) { usa(fc[i].divmod_with_inverse(gc[i], ginv[i]).second); });
        double cal = medir_us([&](int i) { auto r = fw[i].divmod(gw[i]); usa(r.second); });
        double cal_inv = medir_us([&](int i) { usa(P::inverse_series(Kw, revhead(gv[i]), d1)); });
        double cal_resto = medir_us([&](int i) { usa(fw[i].divmod_with_inverse(gw[i], ginv[i]).second); });
        // productos sueltos (lo que hace divmod_with_inverse por dentro)
        double cal_mulFG = medir_us([&](int i) { usa(fft::convolve(Kw, revhead(fv[i]), ginv[i])); });
        double cal_mulOQ = medir_us([&](int i) { usa(fft::convolve(Kw, gv[i], ginv[(i + 1) % 16])); });
        double frio_mulOQ = medir_us([&](int i) { usa(gc[i].mul_kronecker(P(Kc, ginv[(i + 1) % 16]))); });
        std::printf("  librería, cuerpo FRÍO (Kronecker)  : total %8.1f | inversión %8.1f | cociente+resto %8.1f | 1 producto 2000x2000 Kronecker %8.1f\n",
                    frio, frio_inv, frio_resto, frio_mulOQ);
        std::printf("  librería, cuerpo CALIENTE (NTT)    : total %8.1f | inversión %8.1f | cociente+resto %8.1f | producto 2001x2001 NTT %8.1f, 2001x2001 %8.1f\n",
                    cal, cal_inv, cal_resto, cal_mulFG, cal_mulOQ);
        // prototipo a)
        for (int i = 0; i < 16; ++i) {
            comprueba(inv_serie_mp(C, revhead(gv[i]), d1) == ginv[i], "inversión producto medio");
            auto [q, r] = divmod_mp(C, fv[i], gv[i]);
            auto [ql, rl] = fw[i].divmod(gw[i]);
            comprueba(q == ql.coeffs() && r == rl.coeffs(), "divmod_mp");
        }
        double pr_inv = medir_us([&](int i) { usa(inv_serie_mp(C, revhead(gv[i]), d1)); });
        double pr_tot = medir_us([&](int i) { usa(divmod_mp(C, fv[i], gv[i]).second); });
        std::printf("  prototipo a) (tablas ya hechas)    : total %8.1f | inversión %8.1f | cociente+resto %8.1f\n",
                    pr_tot, pr_inv, pr_tot - pr_inv);
        std::printf("  -> ganancia frente a caliente: total %.2fx, inversión %.2fx;  frente a frío: %.2fx\n",
                    cal / pr_tot, cal_inv / pr_inv, frio / pr_tot);
    }

    // ---- 3) pow_mod(a, q, m), Reducer y Composer
    for (int n : {300, 1000}) {
        std::printf("\n== pow_mod(a, q, m), deg m = %d, q = p (29 cuadrados + 7 productos)\n", n);
        std::vector<P> ms, as, bs;
        std::vector<P::Reducer> reds;
        std::vector<ReducerProto> rps;
        for (int i = 0; i < 16; ++i) {
            ms.push_back(P::random(Kw, n));
            as.push_back(P::random(Kw, n - 1));
            bs.push_back(P::random(Kw, n - 1));
            reds.emplace_back(ms.back());
            (void)reds.back().reduce(as.back() * bs.back());  // inverso perezoso ya calculado
            rps.emplace_back(C, ms.back());
        }
        const BigInt q = Kw.cardinality();
        std::vector<P> prods, sqs;
        for (int i = 0; i < 16; ++i) { prods.push_back(as[i] * bs[i]); sqs.push_back(as[i] * as[i]); }
        // comprobaciones
        for (int i = 0; i < 16; ++i) {
            comprueba(rps[i].reduce(prods[i].coeffs()) == reds[i].reduce(prods[i]).coeffs(), "ReducerProto::reduce");
            comprueba(rps[i].mulmod(as[i].coeffs(), bs[i].coeffs()) == reds[i].reduce(prods[i]).coeffs(), "mulmod");
            comprueba(rps[i].sqrmod(as[i].coeffs()) == reds[i].reduce(sqs[i]).coeffs(), "sqrmod");
        }
        for (int i = 0; i < 4; ++i) {
            const P lib = P::pow_mod(as[i], q, reds[i]);
            comprueba(powmod_v1(rps[i], as[i], q) == lib.coeffs(), "powmod_v1");
            comprueba(powmod_v2(rps[i], as[i].coeffs(), q) == lib.coeffs(), "powmod_v2");
            comprueba(P::pow_mod(as[i], q, ms[i]) == lib, "pow_mod con m");
        }
        double t_pow = medir_us([&](int i) { usa(P::pow_mod(as[i], q, reds[i])); });
        double t_sq = medir_us([&](int i) { usa(as[i] * as[i]); });
        double t_mul = medir_us([&](int i) { usa(as[i] * bs[i]); });
        double t_red = medir_us([&](int i) { usa(reds[i].reduce(prods[i])); });
        // interior del reduce de la librería: los dos productos (rev(a)·ginv y m·q)
        const std::size_t d1 = static_cast<std::size_t>(n - 1);
        std::vector<vec> fr, gi, qs;
        for (int i = 0; i < 16; ++i) {
            const vec& a = prods[i].coeffs();
            vec r(d1); for (std::size_t j = 0; j < d1; ++j) r[j] = a[a.size() - 1 - j];
            fr.push_back(r);
            vec rh(n); for (int j = 0; j < n; ++j) rh[j] = ms[i].coeffs()[n - j];
            vec g = P::inverse_series(Kw, rh, n); g.resize(d1); gi.push_back(g);
            qs.push_back(P::random(Kw, n - 2).coeffs());
        }
        double t_red_fg = medir_us([&](int i) { usa(fft::convolve(Kw, fr[i], gi[i])); });
        double t_red_mq = medir_us([&](int i) { usa(fft::convolve(Kw, ms[i].coeffs(), qs[i])); });
        const double modelo = 29 * (t_sq + t_red) + 7 * (t_mul + t_red);
        std::printf("  librería: pow_mod %9.1f | a*a %7.1f  a*b %7.1f (a*a/a*b = %.2f) | reduce %7.1f (rev(a)*ginv %6.1f + m*q %6.1f + resto %6.1f)\n",
                    t_pow, t_sq, t_mul, t_sq / t_mul, t_red, t_red_fg, t_red_mq, t_red - t_red_fg - t_red_mq);
        std::printf("            modelo 29(sq+red)+7(mul+red) = %.1f  (cuadrados %.0f %%, productos %.0f %%, reducciones %.0f %%)\n",
                    modelo, 100 * 29 * t_sq / modelo, 100 * 7 * t_mul / modelo, 100 * 36 * t_red / modelo);
        double p_red = medir_us([&](int i) { usa(rps[i].reduce(prods[i].coeffs())); });
        double p_sq = medir_us([&](int i) { usa(rps[i].sqrmod(as[i].coeffs())); });
        double p_mul = medir_us([&](int i) { usa(rps[i].mulmod(as[i].coeffs(), bs[i].coeffs())); });
        double p_v1 = medir_us([&](int i) { usa(powmod_v1(rps[i], as[i], q)); });
        double p_v2 = medir_us([&](int i) { usa(powmod_v2(rps[i], as[i].coeffs(), q)); });
        std::printf("  prototipo b) reduce %7.1f (%.2fx) | sqr+reduce %7.1f (lib %7.1f, %.2fx) | mul+reduce %7.1f (lib %7.1f, %.2fx)\n",
                    p_red, t_red / p_red, p_sq, t_sq + t_red, (t_sq + t_red) / p_sq, p_mul, t_mul + t_red, (t_mul + t_red) / p_mul);
        std::printf("  pow_mod: lib %9.1f | v1 (productos lib + Reducer b) %9.1f (%.2fx) | v2 (todo en palabras) %9.1f (%.2fx)\n",
                    t_pow, p_v1, t_pow / p_v1, p_v2, t_pow / p_v2);

        // Composer: k = ceil(sqrt(n)) (DDF de factor) y construcción con k = n (is_irreducible)
        std::size_t k = 1; while (k * k < static_cast<std::size_t>(n)) ++k;
        std::vector<P> hs;
        for (int i = 0; i < 16; ++i) hs.push_back(P::pow_mod(P::variable(Kw), q, reds[i]));  // H_1 = x^q mod m
        std::vector<P::Composer> cl;
        std::vector<ComposerProto> cp;
        for (int i = 0; i < 16; ++i) { cl.emplace_back(hs[i], reds[i], k); cp.emplace_back(rps[i], hs[i].coeffs(), k); }
        for (int i = 0; i < 16; ++i) comprueba(cp[i].compose(as[i].coeffs()) == cl[i].compose(as[i]).coeffs(), "ComposerProto::compose");
        double c_lib = medir_us([&](int i) { usa(cl[i].compose(as[i])); });
        double c_pro = medir_us([&](int i) { usa(cp[i].compose(as[i].coeffs())); });
        double cc_lib = medir_us([&](int i) { P::Composer c(hs[i], reds[i], k); usa(c.compose(as[i])); });
        double cc_pro = medir_us([&](int i) { ComposerProto c(rps[i], hs[i].coeffs(), k); usa(c.compose(as[i].coeffs())); });
        std::printf("  Composer k=%zu: compose lib %8.1f | proto %8.1f (%.2fx) ;  construir+compose lib %8.1f | proto %8.1f (%.2fx)\n",
                    k, c_lib, c_pro, c_lib / c_pro, cc_lib, cc_pro, cc_lib / cc_pro);
    }

    // ---- 4) referencia de extremo a extremo: is_irreducible (grado 400) y factor (grado 300)
    {
        std::printf("\n== Extremo a extremo en la librería (referencia para estimar)\n");
        std::vector<P> irr, rnd;
        while (irr.size() < 4) { P f = P::random(Kw, 400); if (f.is_irreducible()) irr.push_back(f.monic()); }
        for (int i = 0; i < 16; ++i) rnd.push_back(P::random(Kw, 300).monic());
        double ti = medir_us([&](int i) { g_chk += irr[i % 4].is_irreducible(); });
        double tf = medir_us([&](int i) { auto r = rnd[i].factor(); g_chk += r.factors.size(); });
        std::printf("  is_irreducible grado 400 (irreducible): %10.1f us\n  factor grado 300 (aleatorio)          : %10.1f us\n", ti, tf);
        // Coste del producto modular a grado 400 (Composer de is_irreducible, k = n)
        std::vector<P> ms, as, bs; std::vector<P::Reducer> reds; std::vector<ReducerProto> rps;
        for (int i = 0; i < 16; ++i) {
            ms.push_back(irr[i % 4]); as.push_back(P::random(Kw, 399)); bs.push_back(P::random(Kw, 399));
            reds.emplace_back(ms.back()); (void)reds.back().reduce(as.back() * bs.back());
            rps.emplace_back(C, ms.back());
        }
        auto fh = rps[0].transformada(bs[0].coeffs());
        for (int i = 0; i < 16; ++i)
            comprueba(rps[i].mulmod_fijo(as[i].coeffs(), rps[i].transformada(bs[i].coeffs())) == reds[i].reduce(as[i] * bs[i]).coeffs(), "mulmod_fijo");
        double lm = medir_us([&](int i) { usa(reds[i].reduce(as[i].mul(bs[i], MulAlgo::FFT))); });
        double pm = medir_us([&](int i) { usa(rps[i].mulmod(as[i].coeffs(), bs[i].coeffs())); });
        double pf = medir_us([&](int i) { usa(rps[0].mulmod_fijo(as[i].coeffs(), fh)); });
        std::printf("  producto modular grado 400: lib %7.1f | proto %7.1f (%.2fx) | proto con operando fijo transformado %7.1f (%.2fx)\n",
                    lm, pm, lm / pm, pf, lm / pf);
    }

    std::printf("\nfallos de comprobación: %d   checksum %llx\n", g_fallos, static_cast<unsigned long long>(g_chk));
    return g_fallos ? 1 : 0;
}
