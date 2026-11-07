// perfil.cpp — Fase 6, tema «factor»: dónde se va el tiempo de factor() con q
// pequeño (F_101) y en extensiones (F_(3^30), F_256), y cuánto cuesta la
// aritmética de los elementos de F_q frente a la de polinomios.
//
// Mide (un hilo, entradas variadas rotando 16, mejor de 3 tandas de >= 30 ms):
//   1) primitivas de elementos de F_q (Fq64 = FqT<Fp64>): mul, inv, add y el
//      camino perezoso acc_clear + acc_fma + acc_reduce;
//   2) factor() desglosado: radical, DDF (Frobenius / MCD, o las piezas de
//      Kaltofen-Shoup), EDF (pow_mod / MCD) y multiplicidad;
//   3) recuento de operaciones de elementos de F_q dentro de factor() con un
//      cuerpo «contador» derivado de Fq64 (mismo resultado: se comprueba).
//
// Compilación (desde la raíz del repo, portable como el módulo):
//   g++ -std=c++23 -O3 -Iinclude bench/experimentos/fase6/factor/perfil.cpp -lgmpxx -lgmp -static -o build/fase6/factor/perfil.exe
// Ejecución:  ALCP_THREADS=1 build/fase6/factor/perfil.exe
#include "alcp/alcp.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <vector>

using namespace alcp;
using reloj = std::chrono::steady_clock;

// Módulos fijos de bench/compare_libs.py (orden creciente).
static const std::vector<int> MOD_2_8 = {1, 1, 0, 1, 1, 0, 0, 0, 1};
static const std::vector<int> MOD_3_30 = {2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0,
                                          2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1};
static std::vector<BigInt> a_big(const std::vector<int>& v) {
    std::vector<BigInt> r;
    for (int x : v) r.emplace_back(x);
    return r;
}

static volatile std::uint64_t sumidero = 0;  // evita que el compilador elimine el trabajo

// Mejor de 3 tandas de >= 30 ms; fn recibe el índice de la entrada (rotando nent).
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
// Cuerpo contador: Fq64 que cuenta las operaciones de elementos.  Las
// funciones ocultan (no sobrescriben) las de la base; Poly<Contador> llama a
// las del tipo derivado.  El resultado debe ser idéntico al de Fq64.
struct Cuentas {
    long mul = 0, inv = 0, add = 0, fma = 0, red = 0, neg = 0;
};
static Cuentas cnt;
struct Contador : Fq64 {
    using Fq64::Fq64;
    Elem mul(const Elem& a, const Elem& b) const { ++cnt.mul; return Fq64::mul(a, b); }
    Elem inv(const Elem& a) const { ++cnt.inv; return Fq64::inv(a); }
    Elem div(const Elem& a, const Elem& b) const { return mul(a, inv(b)); }
    Elem add(const Elem& a, const Elem& b) const { ++cnt.add; return Fq64::add(a, b); }
    Elem sub(const Elem& a, const Elem& b) const { ++cnt.add; return Fq64::sub(a, b); }
    Elem neg(const Elem& a) const { ++cnt.neg; return Fq64::neg(a); }
    void acc_fma(Acc& s, const Elem& a, const Elem& b) const { ++cnt.fma; Fq64::acc_fma(s, a, b); }
    Elem acc_reduce(Acc& s) const { ++cnt.red; return Fq64::acc_reduce(s); }
};

// ---------------------------------------------------------------------------
// 1) Primitivas de elementos
template <class K>
void primitivas_elem(const K& F, const char* nombre) {
    using E = typename K::Elem;
    std::vector<E> a, b;
    for (int i = 0; i < 16; ++i) {
        E x = F.random(), y = F.random();
        if (F.is_zero(y)) y = F.one();
        a.push_back(x);
        b.push_back(y);
    }
    const double tm = medir_us([&](int i) { sumidero += F.mul(a[i], b[i])[0]; });
    const double ti = medir_us([&](int i) { sumidero += F.inv(b[i])[0]; });
    const double ta = medir_us([&](int i) { sumidero += F.add(a[i], b[i])[0]; });
    typename K::Acc acc;
    const double tf = medir_us([&](int i) {
        F.acc_clear(acc);
        F.acc_fma(acc, a[i], b[i]);
        sumidero += F.acc_reduce(acc)[0];
    });
    std::printf("== %s: elemento (ns): mul %.0f | inv %.0f | add %.0f | acc_clear+fma+reduce %.0f\n", nombre,
                tm * 1e3, ti * 1e3, ta * 1e3, tf * 1e3);
}

// ---------------------------------------------------------------------------
// 2) factor() desglosado (réplica de factor() con cronómetros).  La DDF se
// replica como en la librería: clásica (Frobenius + MCD por grado) por debajo
// de 40 y la total (distinct_degree) por encima, cuya mitad KS se desglosa en
// ddf_ks_cronometrada.
struct Tiempos {
    double rad = 0, ddf = 0, frob = 0, gcd = 0, edf = 0, edf_pow = 0, edf_gcd = 0, edf_div = 0, mult = 0;
    int veces = 0, n_gcd = 0, n_split = 0;
};

template <class P>
static double us(reloj::time_point a, reloj::time_point b) {
    return std::chrono::duration<double, std::micro>(b - a).count();
}

// Réplica de equal_degree / equal_degree_split (misma secuencia aleatoria)
// separando pow_mod (o traza), MCD y la división g / d.
template <class P>
std::vector<P> edf_cronometrada(const P& h, int r, Tiempos& T) {
    const auto& K = h.field();
    const BigInt q = K.cardinality();
    const bool par = (K.characteristic() == 2);
    const unsigned long s = par ? (bit_length(q) - 1) : 0;
    std::vector<P> fac, pila{h.monic()};
    while (!pila.empty()) {
        P g = std::move(pila.back());
        pila.pop_back();
        if (g.degree() <= 0) continue;
        if (g.degree() == r) { fac.push_back(g.monic()); continue; }
        P d(K);
        for (;;) {
            ++T.n_split;
            auto t0 = reloj::now();
            P a = P::random(K, g.degree() - 1) % g;
            if (a.degree() <= 0) continue;
            P b(K);
            if (!par) {
                BigInt e = (ipow(q, static_cast<unsigned long>(r)) - 1) / 2;
                b = P::pow_mod(a, e, g) - P::one(K);
            } else {
                P tr = P::zero(K), cur = a;
                for (unsigned long i = 0; i < s * static_cast<unsigned long>(r); ++i) {
                    tr = (tr + cur) % g;
                    cur = (cur * cur) % g;
                }
                b = tr;
            }
            auto t1 = reloj::now();
            d = P::gcd(b, g);
            auto t2 = reloj::now();
            T.edf_pow += us<P>(t0, t1);
            T.edf_gcd += us<P>(t1, t2);
            if (d.degree() > 0 && d.degree() < g.degree()) break;
        }
        auto t3 = reloj::now();
        P otro = g / d;
        T.edf_div += us<P>(t3, reloj::now());
        pila.push_back(d);
        pila.push_back(otro);
    }
    return fac;
}

template <class P>
void desglose(const std::vector<P>& f, const char* nombre) {
    Tiempos T;
    for (int rep = 0; rep < 3; ++rep)
        for (std::size_t i = 0; i < f.size(); ++i) {
            const auto& K = f[i].field();
            auto t0 = reloj::now();
            P g = f[i].radical();
            auto t1 = reloj::now();
            std::vector<std::pair<int, P>> dd;
            if (g.degree() >= 40) {
                dd = g.distinct_degree();  // KS: sin desglose aquí
            } else {
                P gg = g.monic(), x = P::variable(K), H = x % gg;
                typename P::Frobenius frob(gg);
                int d = 1;
                while (gg.degree() >= 2 * d) {
                    auto a = reloj::now();
                    H = frob.siguiente(H);
                    auto b = reloj::now();
                    P hd = P::gcd(H - x, gg);
                    auto c = reloj::now();
                    T.frob += us<P>(a, b);
                    T.gcd += us<P>(b, c);
                    ++T.n_gcd;
                    if (hd.degree() > 0) {
                        dd.emplace_back(d, hd.monic());
                        gg = gg / hd;
                        H = H % gg;
                        frob.cambiar_modulo(gg);
                    }
                    ++d;
                }
                if (gg.degree() > 0) dd.emplace_back(gg.degree(), gg.monic());
            }
            auto t2 = reloj::now();
            std::vector<P> irr;
            for (auto& [d, prod] : dd) {
                auto part = edf_cronometrada(prod, d, T);
                irr.insert(irr.end(), part.begin(), part.end());
            }
            auto t3 = reloj::now();
            for (auto& u : irr) sumidero += f[i].multiplicity(u);
            auto t4 = reloj::now();
            T.rad += us<P>(t0, t1);
            T.ddf += us<P>(t1, t2);
            T.edf += us<P>(t2, t3);
            T.mult += us<P>(t3, t4);
            ++T.veces;
        }
    const double tot = T.rad + T.ddf + T.edf + T.mult, v = T.veces * 1e3;
    std::printf("\n== %s: factor desglosado (media de %d llamadas, ms)\n", nombre, T.veces);
    std::printf("  total %.3f | radical %.3f (%.0f%%) | DDF %.3f (%.0f%%; Frobenius %.3f, MCD %.3f, %d MCD/llamada)\n",
                tot / v, T.rad / v, 100 * T.rad / tot, T.ddf / v, 100 * T.ddf / tot, T.frob / v, T.gcd / v,
                T.n_gcd / T.veces);
    std::printf("  EDF %.3f (%.0f%%; pow_mod/traza %.3f, MCD %.3f, división %.3f, %.1f intentos/llamada) | "
                "multiplicidad %.3f (%.0f%%)\n",
                T.edf / v, 100 * T.edf / tot, T.edf_pow / v, T.edf_gcd / v, T.edf_div / v,
                double(T.n_split) / T.veces, T.mult / v, 100 * T.mult / tot);
}

// factor() completo de la librería, mejor de 3 tandas.
template <class P>
double factor_us(const std::vector<P>& f) {
    return medir_us([&](int i) { sumidero += f[i].factor().factors.size(); }, static_cast<int>(f.size()));
}

// 3) Recuento de operaciones de F_q en factor() (cuerpo contador) y
// comprobación de que da lo mismo que Fq64.
template <class PF>
void recuento(const Contador& KC, const std::vector<PF>& f, const char* nombre) {
    using PC = Poly<Contador>;
    cnt = Cuentas{};
    bool ok = true;
    for (const auto& fi : f) {
        PC fc(KC, fi.coeffs());
        set_seed(77);
        auto rc = fc.factor();
        set_seed(77);
        auto rf = fi.factor();
        ok = ok && rc.factors.size() == rf.factors.size();
        for (std::size_t j = 0; ok && j < rc.factors.size(); ++j)
            ok = rc.factors[j].first.coeffs() == rf.factors[j].first.coeffs() &&
                 rc.factors[j].second == rf.factors[j].second;
    }
    const double n = static_cast<double>(f.size());
    std::printf("  recuento por factor() en %s: mul %.0f, inv %.0f, add/sub %.0f, neg %.0f, acc_fma %.0f, "
                "acc_reduce %.0f  [contador == Fq64: %s]\n",
                nombre, cnt.mul / n, cnt.inv / n, cnt.add / n, cnt.neg / n, cnt.fma / n, cnt.red / n,
                ok ? "sí" : "NO");
}

// Producto mónico de factores aleatorios (como _producto_de_factores de compare_libs.py).
template <class K>
Poly<K> producto_de_factores(const K& F, std::initializer_list<std::pair<int, int>> grados) {
    Poly<K> f = Poly<K>::one(F);
    for (auto [d, m] : grados) {
        Poly<K> u = Poly<K>::random(F, d).monic();
        for (int t = 0; t < m; ++t) f = f * u;
    }
    return f;
}

int main() {
    set_seed(6);
    Fp64 f101(101);
    Fq64 k3(Fp64(3), a_big(MOD_3_30));
    Fq64 k256(Fp64(2), a_big(MOD_2_8));
    Contador c3(Fp64(3), a_big(MOD_3_30));
    Contador c256(Fp64(2), a_big(MOD_2_8));

    std::printf("Perfil de factor() — fase 6, tema factor (un hilo, mejor de 3 tandas de >= 30 ms)\n\n");
    primitivas_elem(k3, "F_(3^30) [Fq64]");
    primitivas_elem(k256, "F_256 [Fq64]");

    using P101 = Poly<Fp64>;
    using PQ = Poly<Fq64>;
    std::vector<P101> f40, f100;
    std::vector<PQ> f5, f30;
    for (int i = 0; i < 16; ++i) {
        f40.push_back(producto_de_factores(f101, {{1, 2}, {2, 1}, {3, 2}, {5, 1}, {8, 1}, {15, 1}}));
        f100.push_back(P101::random(f101, 100).monic());
        f5.push_back(PQ::random(k3, 5).monic());
        f30.push_back(PQ::random(k256, 30).monic());
    }
    std::printf("\nfactor() completo (us): F_101 d38 %.1f | F_101 d100 %.1f | F_(3^30) d5 %.1f | F_256 d30 %.1f\n",
                factor_us(f40), factor_us(f100), factor_us(f5), factor_us(f30));

    desglose(f40, "F_101, producto de grado 38 con multiplicidades");
    desglose(f100, "F_101, grado 100 aleatorio");
    desglose(f5, "F_(3^30), grado 5 aleatorio");
    recuento(c3, f5, "F_(3^30) d5");
    desglose(f30, "F_256, grado 30 aleatorio");
    recuento(c256, f30, "F_256 d30");
    std::printf("\n(sumidero %llu)\n", static_cast<unsigned long long>(sumidero));
    return 0;
}
