// test_fase8_euclides.cpp — Fase 8, propuesta 2: filas del Euclides clásico
// en sitio de Fp64 en AVX2 con p < 2^31 (euclides::avx2::fila2 / axpy).
//
//   1. Diferencial de las piezas: euclides::resto y euclides::menos_producto
//      frente a una copia literal de la versión escalar de la fase 7
//      (namespace ref), con db de 1 a 600 (todas las longitudes hasta 40 y
//      varias no múltiplo de 4 después), cocientes de grado 0, 1 y > 1,
//      coeficientes p - 1, a[da] = 0 y a[db] = 0 tras la primera fila.  Mismo
//      vector resto (entero, incluidos los ceros de arriba), mismo cociente y
//      mismo grado devuelto.
//   2. gcd, gcd_classic, gcd_ext y gcd_ext_classic (g, x, y) frente al
//      Euclides con divmod_school (genérico, no pasa por euclides::), en frío
//      y con la caché FFT llena (half-GCD y su caso base, que usa
//      menos_producto).
//   3. factor: los mismos factores y multiplicidades que en Fp (GMP).
// Primos: 2, 3, 5, 97, 7681, 65537, 998244353, 1073741789, 2147483629 y
// 2147483647 (camino AVX2 si la CPU lo tiene y se compila con optimización) y,
// de control, 2147483659 y 2^61 - 1 (camino escalar).
#include "alcp/alcp.hpp"
#include "check.hpp"

#include <algorithm>
#include <cstdint>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace alcp;
using P = Poly<Fp64>;
using W = std::uint64_t;
using W128 = unsigned __int128;

namespace {

// ------------------------------------------- copia de la versión escalar
// Copia literal de euclides::axpy / resto / menos_producto en 04c1909 (fase 7).
namespace ref {

inline W shoup(W c, W p) { return static_cast<W>((static_cast<W128>(c) << 64) / p); }

inline void axpy(W* r, const W* b, std::size_t len, W c, W cs, W p) {
    for (std::size_t j = 0; j < len; ++j) {
        const W q = static_cast<W>((static_cast<W128>(cs) * b[j]) >> 64);
        W v = c * b[j] - q * p;
        v -= p & (0 - static_cast<W>(v >= p));
        W s = r[j] + v - p;
        s += p & (0 - (s >> 63));
        r[j] = s;
    }
}

inline int grado(const std::vector<W>& a, int hasta) {
    while (hasta >= 0 && a[static_cast<std::size_t>(hasta)] == 0) --hasta;
    return hasta;
}

inline int resto(std::vector<W>& a, int da, const std::vector<W>& b, int db, W inv_lead, W p, std::vector<W>* q) {
    if (q) q->assign(da >= db ? static_cast<std::size_t>(da - db + 1) : 0, 0);
    if (da < db) return da;
    const std::size_t ldb = static_cast<std::size_t>(db);
    if (da == db + 1 && db >= 1 && p < (W(1) << 31)) {
        const auto mulmod = [p](W x, W y) { return x * y % p; };
        const W l1 = a[static_cast<std::size_t>(da)];
        const W f1 = l1 == 0 ? 0 : mulmod(l1, inv_lead);
        const W c1 = f1 == 0 ? 0 : p - f1;
        W adb = a[ldb];
        if (c1 != 0) adb = (c1 * b[ldb - 1] + adb) % p;
        const W f0 = adb == 0 ? 0 : mulmod(adb, inv_lead);
        const W c0 = f0 == 0 ? 0 : p - f0;
        if (q) {
            (*q)[0] = f0;
            (*q)[1] = f1;
        }
        W* r = a.data();
        const W* B = b.data();
        const W mu = ~W(0) / p;
        W prev = 0;
        for (std::size_t j = 0; j < ldb; ++j) {
            const W s = r[j] + c1 * prev + c0 * B[j];
            prev = B[j];
            const W qq = static_cast<W>((static_cast<W128>(s) * mu) >> 64);
            W v = s - qq * p;
            v -= p & (0 - static_cast<W>(v >= p));
            r[j] = v;
        }
        a[static_cast<std::size_t>(da)] = 0;
        a[ldb] = 0;
        return grado(a, db - 1);
    }
    for (int dr = da; dr >= db; --dr) {
        const W lr = a[static_cast<std::size_t>(dr)];
        if (lr == 0) continue;
        const W factor = static_cast<W>(static_cast<W128>(lr) * inv_lead % p);
        if (q) (*q)[static_cast<std::size_t>(dr - db)] = factor;
        const W c = p - factor;
        axpy(a.data() + (dr - db), b.data(), ldb, c, shoup(c, p), p);
        a[static_cast<std::size_t>(dr)] = 0;
    }
    return grado(a, db - 1);
}

inline std::vector<W> menos_producto(const std::vector<W>& x0, const std::vector<W>& q, const std::vector<W>& x1,
                                     W p) {
    std::size_t n = x0.size();
    if (!x1.empty() && !q.empty()) n = std::max(n, q.size() + x1.size() - 1);
    std::vector<W> r(n, 0);
    std::copy(x0.begin(), x0.end(), r.begin());
    if (!x1.empty())
        for (std::size_t i = 0; i < q.size(); ++i)
            if (q[i] != 0) {
                const W c = p - q[i];
                axpy(r.data() + i, x1.data(), x1.size(), c, shoup(c, p), p);
            }
    while (!r.empty() && r.back() == 0) r.pop_back();
    return r;
}

}  // namespace ref

const char* const kPrimos[] = {"2",          "3",          "5",          "97",         "7681",       "65537",
                               "998244353",  "1073741789", "2147483629", "2147483647", "2147483659",
                               "2305843009213693951"};

W powmod(W x, W e, W p) {
    W r = 1 % p;
    for (x %= p; e; e >>= 1, x = static_cast<W>(static_cast<W128>(x) * x % p))
        if (e & 1) r = static_cast<W>(static_cast<W128>(r) * x % p);
    return r;
}
W invmod(W x, W p) { return powmod(x, p - 2, p); }
W mulmod(W x, W y, W p) { return static_cast<W>(static_cast<W128>(x) * y % p); }

// Longitudes de db: todas hasta 40 y varias no múltiplo de 4 (y algunas sí).
std::vector<int> longitudes() {
    std::vector<int> v;
    for (int d = 1; d <= 40; ++d) v.push_back(d);
    for (int d : {63, 64, 65, 66, 127, 128, 129, 255, 257, 301, 511, 513, 598, 599, 600}) v.push_back(d);
    return v;
}

// Un caso de resto: compara la versión nueva con la copia escalar.
bool un_resto(std::vector<W> a, int da, const std::vector<W>& b, int db, W p) {
    const W il = invmod(b[static_cast<std::size_t>(db)], p);
    std::vector<W> a2 = a, q1, q2;
    const int g1 = euclides::resto(a, da, b, db, il, p, &q1);
    const int g2 = ref::resto(a2, da, b, db, il, p, &q2);
    return g1 == g2 && a == a2 && q1 == q2;
}

// ------------------------------------------------------------ 1. piezas
void test_piezas(const char* ps) {
    const std::string tag = std::string(" p = ") + ps;
    const W p = static_cast<W>(std::stoull(ps));
    std::mt19937_64 g(12345 + p);
    auto rnd = [&] { return g() % p; };
    auto nz = [&] { return 1 + g() % (p - 1 == 0 ? 1 : p - 1); };  // en [1, p)
    bool ok_r = true, ok_m = true;
    for (int db : longitudes()) {
        for (int extra : {0, 1, 2, 3, 5}) {  // grado del cociente
            const int da = db + extra;
            std::vector<W> b(static_cast<std::size_t>(db) + 1), a(static_cast<std::size_t>(da) + 1);
            // Aleatorios.
            for (auto& x : b) x = rnd();
            b.back() = nz();
            for (auto& x : a) x = rnd();
            ok_r = ok_r && un_resto(a, da, b, db, p);
            // Sin cociente: también la llamada con q = nullptr.
            {
                std::vector<W> a1 = a, a2 = a;
                const W il = invmod(b.back(), p);
                ok_r = ok_r && euclides::resto(a1, da, b, db, il, p, nullptr) ==
                                   ref::resto(a2, da, b, db, il, p, nullptr) &&
                       a1 == a2;
            }
            // Coeficientes p - 1 (cotas de las sumas sin reducir).
            std::vector<W> bm(b.size(), p - 1), am(a.size(), p - 1);
            ok_r = ok_r && un_resto(am, da, bm, db, p);
            // a[da] = 0 (el factor de la primera fila es 0).
            std::vector<W> a0 = a;
            a0.back() = 0;
            ok_r = ok_r && un_resto(a0, da, b, db, p);
            // a[db] = 0 tras la primera fila (cociente de grado 1: f0 = 0).
            if (extra == 1) {
                std::vector<W> az = a;
                const W f1 = mulmod(az.back(), invmod(b.back(), p), p);
                az[static_cast<std::size_t>(db)] = mulmod(f1, b[static_cast<std::size_t>(db) - 1], p);
                ok_r = ok_r && un_resto(az, da, b, db, p);
                // Y con a[da] = 0 a la vez (las dos filas nulas).
                az.back() = 0;
                az[static_cast<std::size_t>(db)] = 0;
                ok_r = ok_r && un_resto(az, da, b, db, p);
            }
        }
        // menos_producto: x0 - q·x1 con cocientes de 1 a 4 coeficientes.
        for (int lq : {1, 2, 3, 4}) {
            std::vector<W> x0(static_cast<std::size_t>(db) + static_cast<std::size_t>(g() % 3)),
                x1(static_cast<std::size_t>(db) + 1), q(static_cast<std::size_t>(lq));
            for (auto& x : x0) x = rnd();
            for (auto& x : x1) x = rnd();
            for (auto& x : q) x = rnd();
            x1.back() = nz();
            q.back() = nz();
            ok_m = ok_m && euclides::menos_producto(x0, q, x1, p) == ref::menos_producto(x0, q, x1, p);
            std::vector<W> m0(x0.size(), p - 1), m1(x1.size(), p - 1), mq(q.size(), p - 1);
            ok_m = ok_m && euclides::menos_producto(m0, mq, m1, p) == ref::menos_producto(m0, mq, m1, p);
            // x0 = q·x1 exactamente: el resultado se anula y se recorta.
            const std::vector<W> vacio;
            const std::vector<W> qx = ref::menos_producto(vacio, q, x1, p);  // -q·x1
            std::vector<W> mq2(q.size());
            for (std::size_t i = 0; i < q.size(); ++i) mq2[i] = q[i] == 0 ? 0 : p - q[i];
            ok_m = ok_m && euclides::menos_producto(qx, mq2, x1, p).empty();
        }
    }
    CHECK_MSG(ok_r, ("euclides::resto = versión escalar" + tag).c_str());
    CHECK_MSG(ok_m, ("euclides::menos_producto = versión escalar" + tag).c_str());
}

// ------------------------------------------------------------ 2. gcd
// Referencias con divmod_school (genérico, no pasa por euclides::).
P gcd_ref(P a, P b) {
    while (!b.is_zero()) {
        P r = a.divmod_school(b).second;
        a = std::move(b);
        b = std::move(r);
    }
    return a.is_zero() ? a : a.monic();
}

P::ExtGcd gcd_ext_ref(const P& a, const P& b) {
    const Fp64& K = a.field();
    P r0 = a, r1 = b, x0 = P::one(K), x1 = P::zero(K), y0 = P::zero(K), y1 = P::one(K);
    while (!r1.is_zero()) {
        auto [q, r] = r0.divmod_school(r1);
        P x2 = x0 - q.mul_school(x1), y2 = y0 - q.mul_school(y1);
        r0 = std::move(r1);
        r1 = std::move(r);
        x0 = std::move(x1);
        x1 = std::move(x2);
        y0 = std::move(y1);
        y1 = std::move(y2);
    }
    if (r0.is_zero()) return {r0, x0, y0};
    const auto il = K.inv(r0.lead());
    return {r0.scalar_mul(il), x0.scalar_mul(il), y0.scalar_mul(il)};
}

bool igual(const P::ExtGcd& u, const P::ExtGcd& v) { return u.g == v.g && u.x == v.x && u.y == v.y; }

void test_gcd(const char* ps) {
    const std::string tag = std::string(" p = ") + ps;
    set_seed(8);
    const Fp64 K{BigInt(ps)};
    bool ok = true;
    for (int caliente = 0; caliente < 2; ++caliente) {
        if (caliente) (void)(P::random(K, 300) * P::random(K, 300));  // llena la caché FFT si hay NTT
        for (int d : {1, 2, 5, 13, 31, 64, 97, 130, 301, 600}) {
            const P a = P::random(K, d + 1), b = P::random(K, d);
            const P gr = gcd_ref(a, b);
            ok = ok && P::gcd_classic(a, b) == gr && P::gcd(a, b) == gr;
            const auto er = gcd_ext_ref(a, b);
            ok = ok && igual(P::gcd_ext_classic(a, b), er) && igual(P::gcd_ext(a, b), er);
            // Con factor común (mcd no trivial) y cocientes de grado > 1.
            const P g = P::random(K, d / 3 + 1);
            const P x = g.mul_school(P::random(K, d)), y = g.mul_school(P::random(K, std::max(0, d - 3)));
            const P gr2 = gcd_ref(x, y);
            ok = ok && P::gcd_classic(x, y) == gr2 && P::gcd(x, y) == gr2;
            const auto er2 = gcd_ext_ref(x, y);
            ok = ok && igual(P::gcd_ext_classic(x, y), er2) && igual(P::gcd_ext(x, y), er2);
        }
        // Coeficientes p - 1.
        const W pm1 = K.modulus64() - 1;
        const P A(K, std::vector<W>(402, pm1)), B(K, std::vector<W>(401, pm1));
        std::vector<W> vb(399, pm1);
        vb[100] = 0;
        const P C(K, vb);
        ok = ok && P::gcd_classic(A, C) == gcd_ref(A, C) && igual(P::gcd_ext_classic(A, C), gcd_ext_ref(A, C));
        ok = ok && P::gcd_classic(A, B) == gcd_ref(A, B) && igual(P::gcd_ext(A, B), gcd_ext_ref(A, B));
    }
    CHECK_MSG(ok, ("gcd y gcd_ext (g, x, y) = referencia" + tag).c_str());
}

// ------------------------------------------------------------ 3. factor
std::string texto(W x) { return std::to_string(x); }
std::string texto(const BigInt& x) { return x.get_str(); }

std::vector<std::pair<std::vector<std::string>, int>> forma(const auto& F) {
    std::vector<std::pair<std::vector<std::string>, int>> v;
    for (const auto& [f, m] : F.factors) {
        std::vector<std::string> c;
        for (const auto& x : f.coeffs()) c.push_back(texto(x));
        v.emplace_back(std::move(c), m);
    }
    std::sort(v.begin(), v.end());
    return v;
}

void test_factor(const char* ps) {
    const std::string tag = std::string(" p = ") + ps;
    set_seed(21);
    const Fp64 K{BigInt(ps)};
    const Fp G{BigInt(ps)};
    bool ok = true;
    for (int d : {9, 40, 90}) {
        const P u = P::random(K, d), v = P::random(K, d / 3 + 1), w = P::random(K, 3);
        const P f = u.mul_school(v).mul_school(v).mul_school(w);
        std::vector<BigInt> cg;
        for (const auto& x : f.coeffs()) cg.emplace_back(texto(x));
        const Poly<Fp> fg(G, cg);
        ok = ok && forma(f.factor()) == forma(fg.factor());
    }
    CHECK_MSG(ok, ("factor = factor en Fp (GMP)" + tag).c_str());
}

}  // namespace

int main() {
    for (const char* ps : kPrimos) {
        test_piezas(ps);
        test_gcd(ps);
    }
    for (const char* ps : {"2", "3", "97", "65537", "998244353", "2147483647", "2147483659"}) test_factor(ps);
    return check::summary("test_fase8_euclides");
}
