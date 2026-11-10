// fpn_bench.cpp — Fase 6, tema fpn: aritmética de tamaño fijo para F_p con p
// de 2 palabras frente al Fp actual (mpz).  PROTOTIPO FUERA DE LA LIBRERÍA.
//
// Secciones (argumento opcional: "pool" sustituye el asignador de GMP por
// listas libres por tamaño, para medir cuánto pesan las reservas de mpz):
//   1. Elementos: mul, add, sub, inv y el núcleo «20 acc_fma + acc_reduce»
//      (un coeficiente de mul_school) en Fp (mpz) y FpN2 (fpn2.hpp).
//   2. Polinomios de grado 20 en F_(2^127-1): producto escolar, divmod 38/20 y
//      pow_mod con exponente q = 2^127-1 (paso de Frobenius) en Poly<Fp> y en
//      Poly<FpN2>, comprobando resultado idéntico (convertido a BigInt).
//   3. factor() de 16 polinomios mónicos aleatorios de grado 20 (escenario
//      fp127.factor.d20): tiempo total con Fp y con FpN2 (factores idénticos),
//      desglose por fases (radical, DDF: pow_mod/composición/MCD, EDF:
//      pow_mod/resto, multiplicidad), recuento de operaciones de elemento
//      (cuerpo contador FpCont) y de reservas (GMP y operator new) por fase.
//   4. Variante de la DDF que compone desde el primer paso tras cada
//      cambio de módulo (observación lateral, mismo resultado).
//   5. La entrada exacta de fp127.factor.d20 (build/fase6/fpn/fp127_factor_d20.txt,
//      generada por flint_fp127.py con la semilla de compare_libs.py); con
//      una sola entrada no se puede rotar (se avisa en la salida).
//
// Metodología: entradas variadas (16, rotando), mejor de 3 tandas de >= 30 ms
// (y un múltiplo de 16 llamadas), ALCP_THREADS=1, set_seed(1000 + i) antes de
// cada factor() para que las dos variantes vean el mismo generador.
//
//   g++ -std=c++23 -O3 -Iinclude -Ibench/experimentos/fase6/fpn bench/experimentos/fase6/fpn/fpn_bench.cpp -lgmpxx -lgmp -static -o build/fase6/fpn/fpn_bench.exe
//
// Ejecución (desde la raíz del repositorio):
//   ALCP_THREADS=1 ./build/fase6/fpn/fpn_bench.exe          (asignador normal)
//   ALCP_THREADS=1 ./build/fase6/fpn/fpn_bench.exe pool     (GMP con listas libres)
#include "fpn2.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <string>
#include <vector>

using namespace alcp;
using fpn::E2;
using fpn::FpN2;
using reloj = std::chrono::steady_clock;

// ---------------------------------------------------------------------------
// Recuento de reservas: operator new global y funciones de memoria de GMP.
// ---------------------------------------------------------------------------
struct ContMem {
    long long news = 0, gmp_alloc = 0, gmp_realloc = 0;
};
static ContMem g_mem;

void* operator new(std::size_t n) {
    ++g_mem.news;
    if (void* p = std::malloc(n ? n : 1)) return p;
    throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

// Asignador de GMP: modo normal (malloc) o listas libres por clases de 16 bytes.
static bool g_pool = false;
static void* g_libres[65] = {};  // clases 1..64 (16..1024 bytes)
static inline std::size_t clase(std::size_t n) { return (n + 15) / 16; }
static void* gmp_alloc(std::size_t n) {
    ++g_mem.gmp_alloc;
    if (g_pool) {
        const std::size_t c = clase(n);
        if (c <= 64) {
            if (void* p = g_libres[c]) {
                g_libres[c] = *static_cast<void**>(p);
                return p;
            }
            return std::malloc(c * 16);
        }
    }
    return std::malloc(n);
}
static void gmp_free(void* p, std::size_t n) {
    if (g_pool) {
        const std::size_t c = clase(n);
        if (c <= 64) {
            *static_cast<void**>(p) = g_libres[c];
            g_libres[c] = p;
            return;
        }
    }
    std::free(p);
}
static void* gmp_realloc(void* p, std::size_t viejo, std::size_t nuevo) {
    ++g_mem.gmp_realloc;
    if (g_pool) {
        if (clase(viejo) == clase(nuevo) && clase(nuevo) <= 64) return p;
        void* q = gmp_alloc(nuevo);
        --g_mem.gmp_alloc;
        std::memcpy(q, p, std::min(viejo, nuevo));
        gmp_free(p, viejo);
        return q;
    }
    return std::realloc(p, nuevo);
}

// ---------------------------------------------------------------------------
// Cuerpo contador: Fp que cuenta las operaciones de elemento.
// ---------------------------------------------------------------------------
struct ContOps {
    long long mul = 0, add = 0, sub = 0, neg = 0, inv = 0, fma = 0, red = 0, set = 0;
};
static ContOps g_ops;

struct FpCont : Fp {
    using Fp::Fp;
    Elem add(const Elem& a, const Elem& b) const { ++g_ops.add; return Fp::add(a, b); }
    Elem sub(const Elem& a, const Elem& b) const { ++g_ops.sub; return Fp::sub(a, b); }
    Elem neg(const Elem& a) const { ++g_ops.neg; return Fp::neg(a); }
    Elem mul(const Elem& a, const Elem& b) const { ++g_ops.mul; return Fp::mul(a, b); }
    Elem inv(const Elem& a) const { ++g_ops.inv; return Fp::inv(a); }
    Elem div(const Elem& a, const Elem& b) const { return mul(a, inv(b)); }
    void acc_set(Acc& s, const Elem& a) const { ++g_ops.set; Fp::acc_set(s, a); }
    void acc_fma(Acc& s, const Elem& a, const Elem& b) const { ++g_ops.fma; Fp::acc_fma(s, a, b); }
    Elem acc_reduce(Acc& s) const { ++g_ops.red; return Fp::acc_reduce(s); }
};

// ---------------------------------------------------------------------------
// Utilidades de medida y conversión.
// ---------------------------------------------------------------------------
static double ahora_us() {
    return std::chrono::duration<double, std::micro>(reloj::now().time_since_epoch()).count();
}

// Mejor de 3 tandas de >= 30 ms y un múltiplo de n llamadas; fn(i) con i en [0, n).
template <class Fn>
double medir_us(int n, Fn&& fn) {
    double mejor = 1e300;
    for (int tanda = 0; tanda < 3; ++tanda) {
        long veces = 0;
        const double t0 = ahora_us();
        double dt = 0;
        do {
            for (int i = 0; i < n; ++i) fn(i);
            veces += n;
            dt = ahora_us() - t0;
        } while (dt < 30000);
        mejor = std::min(mejor, dt / static_cast<double>(veces));
    }
    return mejor;
}

static std::uint64_t g_chk = 0;  // suma de control: impide eliminar el trabajo
static void usa(const BigInt& x) { g_chk = g_chk * 31 + (mpz_size(x.get_mpz_t()) ? mpz_getlimbn(x.get_mpz_t(), 0) : 0); }
static void usa(const E2& x) { g_chk = g_chk * 31 + x.lo; }
static void usa(std::size_t x) { g_chk = g_chk * 31 + x; }

static BigInt big(const BigInt& x) { return x; }
static BigInt big(const E2& x) { return fpn::to_big(x); }

template <class F>
std::vector<BigInt> coef_big(const Poly<F>& a) {
    std::vector<BigInt> r;
    for (const auto& c : a.coeffs()) r.push_back(big(c));
    return r;
}

template <class F>
Poly<F> convierte(const F& K, const Poly<Fp>& a) {
    std::vector<typename F::Elem> c;
    for (const auto& x : a.coeffs()) c.push_back(K.from_int(x));
    return Poly<F>(K, std::move(c));
}

static int g_fallos = 0;
static void comprueba(bool ok, const char* que) {
    if (!ok) {
        ++g_fallos;
        std::printf("  *** DISTINTO: %s\n", que);
    }
}

// Factorización como lista de (coeficientes en BigInt, multiplicidad).
template <class F>
std::vector<std::pair<std::vector<BigInt>, int>> fact_big(const typename Poly<F>::Factorization& fz) {
    std::vector<std::pair<std::vector<BigInt>, int>> r;
    for (const auto& [u, e] : fz.factors) r.emplace_back(coef_big(u), e);
    return r;
}

// Costes medidos en la sección 1 para F_(2^127-1) con Fp (ns), para estimar
// cuánto de factor() es aritmética de elementos.
struct Costes {
    double mul = 0, add = 0, inv = 0, inv1 = 0, fma = 0, red = 0, set = 0, malloc = 0;
};
static Costes g_cost;

// ---------------------------------------------------------------------------
// 1. Elementos
// ---------------------------------------------------------------------------
static void seccion_elementos(const Fp& K, const FpN2& KN) {
    constexpr int M = 1024;
    std::vector<BigInt> x(M), y(M);
    std::vector<E2> xn(M), yn(M);
    for (int i = 0; i < M; ++i) {
        do x[i] = K.random(); while (x[i] == 0);
        y[i] = K.random();
        xn[i] = KN.from_int(x[i]);
        yn[i] = KN.from_int(y[i]);
    }
    // Comprobación de resultados idénticos en las 1024 entradas.
    for (int i = 0; i < M; ++i) {
        const int j = (i * 7 + 3) & (M - 1);
        comprueba(big(KN.mul(xn[i], yn[j])) == K.mul(x[i], y[j]), "mul");
        comprueba(big(KN.add(xn[i], yn[j])) == K.add(x[i], y[j]), "add");
        comprueba(big(KN.sub(xn[i], yn[j])) == K.sub(x[i], y[j]), "sub");
        comprueba(big(KN.neg(yn[j])) == K.neg(y[j]), "neg");
        comprueba(big(KN.inv(xn[i])) == K.inv(x[i]), "inv");
        Fp::Acc s;
        FpN2::Acc sn;
        K.acc_clear(s);
        KN.acc_clear(sn);
        for (int t = 0; t < 20; ++t) {
            K.acc_fma(s, x[(i + t) & (M - 1)], y[t]);
            KN.acc_fma(sn, xn[(i + t) & (M - 1)], yn[t]);
        }
        comprueba(big(KN.acc_reduce(sn)) == K.acc_reduce(s), "fma20");
    }
    // Caso extremo de la acumulación: 2^20 productos (p-1)^2.
    {
        const BigInt pm1 = K.characteristic() - 1;
        const E2 pm1n = KN.from_int(pm1);
        Fp::Acc s;
        FpN2::Acc sn;
        K.acc_clear(s);
        KN.acc_clear(sn);
        for (int t = 0; t < (1 << 20); ++t) {
            K.acc_fma(s, pm1, pm1);
            KN.acc_fma(sn, pm1n, pm1n);
        }
        comprueba(big(KN.acc_reduce(sn)) == K.acc_reduce(s), "fma 2^20 extremo");
    }

    std::printf("\n== 1. Elementos de F_p, p = %s (ns por operación, 1024 entradas rotando)\n",
                K.characteristic().get_str().c_str());
    std::printf("%-34s %10s %10s %8s\n", "operación", "Fp (mpz)", "FpN2", "factor");
    auto fila = [](const char* nom, double a, double b) { std::printf("%-34s %10.1f %10.1f %7.1fx\n", nom, a, b, a / b); };
    auto por_op = [&](auto&& fn) { return medir_us(1, [&](int) { for (int i = 0; i < M; ++i) fn(i); }) * 1000.0 / M; };

    g_cost.mul = por_op([&](int i) { usa(K.mul(x[i], y[(i * 7 + 3) & (M - 1)])); });
    fila("mul (devuelve elemento nuevo)", g_cost.mul, por_op([&](int i) { usa(KN.mul(xn[i], yn[(i * 7 + 3) & (M - 1)])); }));
    g_cost.add = por_op([&](int i) { usa(K.add(x[i], y[(i * 7 + 3) & (M - 1)])); });
    fila("add", g_cost.add, por_op([&](int i) { usa(KN.add(xn[i], yn[(i * 7 + 3) & (M - 1)])); }));
    fila("sub", por_op([&](int i) { usa(K.sub(x[i], y[(i * 7 + 3) & (M - 1)])); }),
         por_op([&](int i) { usa(KN.sub(xn[i], yn[(i * 7 + 3) & (M - 1)])); }));
    g_cost.inv = por_op([&](int i) { usa(K.inv(x[i])); });
    fila("inv", g_cost.inv, por_op([&](int i) { usa(KN.inv(xn[i])); }));
    // Cadena dependiente (latencia): r = r·x[i]
    {
        BigInt r = 1;
        E2 rn = KN.one();
        const double a = por_op([&](int i) { r = K.mul(r, x[i]); });
        const double b = por_op([&](int i) { rn = KN.mul(rn, xn[i]); });
        usa(r);
        usa(rn);
        fila("mul en cadena (latencia)", a, b);
    }
    // Núcleo de mul_school: un coeficiente = 20 acc_fma + 1 acc_reduce, con el
    // acumulador nuevo en cada coeficiente (como el vector acc de mul_school).
    const double tf = por_op([&](int i) {
        Fp::Acc s;
        K.acc_clear(s);
        for (int t = 0; t < 20; ++t) K.acc_fma(s, x[(i + t) & (M - 1)], y[t]);
        usa(K.acc_reduce(s));
    });
    const double tfn = por_op([&](int i) {
        FpN2::Acc s;
        KN.acc_clear(s);
        for (int t = 0; t < 20; ++t) KN.acc_fma(s, xn[(i + t) & (M - 1)], yn[t]);
        usa(KN.acc_reduce(s));
    });
    fila("20 acc_fma + acc_reduce", tf, tfn);
    fila("  por acc_fma (aprox. /20)", tf / 20, tfn / 20);
    // Piezas sueltas: acc_fma sobre un acumulador ya reservado, acc_set y
    // acc_reduce (que devuelve un elemento nuevo), inv(1) (el divisor mónico
    // de los restos de pow_mod).
    {
        Fp::Acc s;
        FpN2::Acc sn;
        K.acc_clear(s);
        KN.acc_clear(sn);
        const double a = por_op([&](int i) { K.acc_fma(s, x[i], y[(i * 7 + 3) & (M - 1)]); });
        const double b = por_op([&](int i) { KN.acc_fma(sn, xn[i], yn[(i * 7 + 3) & (M - 1)]); });
        usa(K.acc_reduce(s));
        usa(KN.acc_reduce(sn));
        fila("acc_fma (acumulador ya reservado)", a, b);
        g_cost.fma = a;
    }
    {
        std::vector<Fp::Acc> s(M);
        std::vector<FpN2::Acc> sn(M);
        for (int i = 0; i < M; ++i) {
            K.acc_clear(s[i]);
            KN.acc_clear(sn[i]);
            for (int t = 0; t < 20; ++t) {
                K.acc_fma(s[i], x[(i + t) & (M - 1)], y[t]);
                KN.acc_fma(sn[i], xn[(i + t) & (M - 1)], yn[t]);
            }
        }
        // acc_reduce no consume el acumulador en ninguno de los dos cuerpos.
        const double a = por_op([&](int i) { usa(K.acc_reduce(s[i])); });
        const double b = por_op([&](int i) { usa(KN.acc_reduce(sn[i])); });
        fila("acc_reduce (suma de 20 productos)", a, b);
        g_cost.red = a;
        std::vector<Fp::Acc> s2(M);
        std::vector<FpN2::Acc> sn2(M);
        const double c = por_op([&](int i) { K.acc_set(s2[i], x[i]); });
        const double d = por_op([&](int i) { KN.acc_set(sn2[i], xn[i]); });
        fila("acc_set", c, d);
        g_cost.set = c;
    }
    {
        const BigInt uno = 1;
        const double a = por_op([&](int) { usa(K.inv(uno)); });
        const double b = por_op([&](int) { usa(KN.inv(KN.one())); });
        fila("inv(1)", a, b);
        g_cost.inv1 = a;
    }
    // malloc/free de 16-32 bytes: referencia para el coste de una reserva.
    static void* volatile sumidero;
    const double tm = por_op([&](int i) {
        void* p = std::malloc(16 + (i & 16));
        sumidero = p;
        std::free(sumidero);
    });
    std::printf("%-34s %10.1f\n", "malloc+free 16-32 B (referencia)", tm);
    g_cost.malloc = tm;
}

// ---------------------------------------------------------------------------
// 2. Polinomios de grado 20
// ---------------------------------------------------------------------------
static void seccion_polinomios(const Fp& K, const FpN2& KN) {
    using P = Poly<Fp>;
    using PN = Poly<FpN2>;
    constexpr int N = 16;
    std::vector<P> a, b, c, m;
    for (int i = 0; i < N; ++i) {
        a.push_back(P::random(K, 19));
        b.push_back(P::random(K, 19));
        c.push_back(P::random(K, 38));
        m.push_back(P::random(K, 20).monic());
    }
    std::vector<PN> an, bn, cn, mn;
    for (int i = 0; i < N; ++i) {
        an.push_back(convierte(KN, a[i]));
        bn.push_back(convierte(KN, b[i]));
        cn.push_back(convierte(KN, c[i]));
        mn.push_back(convierte(KN, m[i]));
    }
    const BigInt q = K.cardinality();
    const BigInt e1 = (q - 1) / 2;  // exponente de la EDF con r = 1
    // Comprobación
    for (int i = 0; i < N; ++i) {
        comprueba(coef_big(an[i] * bn[i]) == coef_big(a[i] * b[i]), "mul d19");
        comprueba(coef_big(cn[i] % mn[i]) == coef_big(c[i] % m[i]), "divmod 38/20");
        comprueba(coef_big(PN::pow_mod(an[i], q, mn[i])) == coef_big(P::pow_mod(a[i], q, m[i])), "pow_mod q");
        comprueba(coef_big(PN::gcd(an[i], mn[i])) == coef_big(P::gcd(a[i], m[i])), "gcd 20");
    }
    std::printf("\n== 2. Polinomios de grado 20 en F_(2^127-1) (us, 16 entradas rotando)\n");
    std::printf("%-34s %10s %10s %8s\n", "operación", "Poly<Fp>", "Poly<FpN2>", "factor");
    auto fila = [](const char* nom, double x, double y) { std::printf("%-34s %10.2f %10.2f %7.1fx\n", nom, x, y, x / y); };
    fila("producto 19x19 (escolar)", medir_us(N, [&](int i) { usa((a[i] * b[i]).coeffs()[0]); }),
         medir_us(N, [&](int i) { usa((an[i] * bn[i]).coeffs()[0]); }));
    fila("resto 38 / 20 (escolar)", medir_us(N, [&](int i) { usa((c[i] % m[i]).coeffs()[0]); }),
         medir_us(N, [&](int i) { usa((cn[i] % mn[i]).coeffs()[0]); }));
    fila("gcd grado 19 / 20", medir_us(N, [&](int i) { usa(P::gcd(a[i], m[i]).coeffs()[0]); }),
         medir_us(N, [&](int i) { usa(PN::gcd(an[i], mn[i]).coeffs()[0]); }));
    std::vector<P::Reducer> red;
    std::vector<PN::Reducer> redn;
    for (int i = 0; i < N; ++i) {
        red.emplace_back(m[i]);
        redn.emplace_back(mn[i]);
    }
    fila("pow_mod exponente q (127 bits)", medir_us(N, [&](int i) { usa(P::pow_mod(a[i], q, red[i]).coeffs()[0]); }),
         medir_us(N, [&](int i) { usa(PN::pow_mod(an[i], q, redn[i]).coeffs()[0]); }));
    fila("pow_mod exponente (q-1)/2", medir_us(N, [&](int i) { usa(P::pow_mod(a[i], e1, red[i]).coeffs()[0]); }),
         medir_us(N, [&](int i) { usa(PN::pow_mod(an[i], e1, redn[i]).coeffs()[0]); }));
}

// ---------------------------------------------------------------------------
// 3. Desglose de factor()
// ---------------------------------------------------------------------------
struct Fases {
    double rad = 0, ddf_pow = 0, ddf_comp = 0, ddf_gcd = 0, ddf_otro = 0, edf_pow = 0, edf_otro = 0, mult = 0;
    long n_pow = 0, n_comp = 0, n_gcd = 0, n_split = 0, n_edf_pow = 0;
    double total() const { return rad + ddf_pow + ddf_comp + ddf_gcd + ddf_otro + edf_pow + edf_otro + mult; }
};

// Instantánea de contadores para el recuento por fase.
struct Inst {
    ContMem m;
    ContOps o;
};
static Inst inst() { return {g_mem, g_ops}; }
struct RecFase {
    long long news = 0, gmp = 0, mul = 0, addsub = 0, inv = 0, fma = 0, red = 0, set = 0;
};
static void suma(RecFase& r, const Inst& a, const Inst& b) {
    r.news += b.m.news - a.m.news;
    r.gmp += (b.m.gmp_alloc - a.m.gmp_alloc) + (b.m.gmp_realloc - a.m.gmp_realloc);
    r.mul += b.o.mul - a.o.mul;
    r.addsub += (b.o.add - a.o.add) + (b.o.sub - a.o.sub) + (b.o.neg - a.o.neg);
    r.inv += b.o.inv - a.o.inv;
    r.fma += b.o.fma - a.o.fma;
    r.red += b.o.red - a.o.red;
    r.set += b.o.set - a.o.set;
}
enum { F_RAD, F_DDF_POW, F_DDF_COMP, F_DDF_GCD, F_DDF_OTRO, F_EDF_POW, F_EDF_OTRO, F_MULT, F_N };
static const char* kNomFase[F_N] = {"radical", "DDF pow_mod", "DDF composición", "DDF MCD", "DDF otros",
                                    "EDF pow_mod", "EDF resto", "multiplicidad"};

// Reproduce factor() (camino de grado < 40: DDF clásica con Frobenius, EDF de
// Cantor-Zassenhaus) separando fases.  Devuelve la factorización ordenada.
template <class F>
std::vector<std::pair<std::vector<BigInt>, int>> factor_desglosado(const Poly<F>& f, Fases& T, RecFase* rec) {
    using P = Poly<F>;
    const F& K = f.field();
    double t = ahora_us();
    Inst s = inst();
    auto marca = [&](double& acc, int fase) {
        const double u = ahora_us();
        acc += u - t;
        t = u;
        const Inst s2 = inst();
        if (rec) suma(rec[fase], s, s2);
        s = s2;
    };
    P g = f.radical();
    marca(T.rad, F_RAD);
    // --- DDF (distinct_degree_impl sin Kaltofen-Shoup: grado < 40) ---
    std::vector<std::pair<int, P>> dd;
    {
        P gg = g.monic(), x = P::variable(K), H = x % gg;
        typename P::Frobenius frob(gg);
        bool recien = true;  // el primer paso tras construir o cambiar el módulo es pow_mod
        int d = 1;
        marca(T.ddf_otro, F_DDF_OTRO);
        while (gg.degree() >= 2 * d) {
            H = frob.siguiente(H);
            if (recien) { marca(T.ddf_pow, F_DDF_POW); ++T.n_pow; }
            else { marca(T.ddf_comp, F_DDF_COMP); ++T.n_comp; }
            recien = false;
            P hd = P::gcd(H - x, gg);
            marca(T.ddf_gcd, F_DDF_GCD);
            ++T.n_gcd;
            if (hd.degree() > 0) {
                dd.emplace_back(d, hd.monic());
                gg = gg / hd;
                H = H % gg;
                frob.cambiar_modulo(gg);
                recien = true;
            }
            ++d;
            marca(T.ddf_otro, F_DDF_OTRO);
        }
        if (gg.degree() > 0) dd.emplace_back(gg.degree(), gg.monic());
    }
    marca(T.ddf_otro, F_DDF_OTRO);
    // --- EDF (equal_degree + equal_degree_split, q impar) ---
    std::vector<P> irr;
    const BigInt q = K.cardinality();
    for (auto& [r, prod] : dd) {
        std::vector<P> pila{prod.monic()};
        while (!pila.empty()) {
            P gq = std::move(pila.back());
            pila.pop_back();
            if (gq.degree() <= 0) continue;
            if (gq.degree() == r) { irr.push_back(gq.monic()); continue; }
            P dv(K);
            for (;;) {
                P a = P::random(K, gq.degree() - 1) % gq;
                if (a.degree() <= 0) continue;
                const BigInt e = (ipow(q, static_cast<unsigned long>(r)) - 1) / 2;
                marca(T.edf_otro, F_EDF_OTRO);
                P b = P::pow_mod(a, e, gq) - P::one(K);
                marca(T.edf_pow, F_EDF_POW);
                ++T.n_edf_pow;
                dv = P::gcd(b, gq);
                if (dv.degree() > 0 && dv.degree() < gq.degree()) break;
            }
            ++T.n_split;
            pila.push_back(dv);
            pila.push_back(gq / dv);
        }
    }
    marca(T.edf_otro, F_EDF_OTRO);
    // --- multiplicidades y orden canónico ---
    std::vector<std::pair<P, int>> fs;
    for (auto& u : irr) fs.emplace_back(u, f.multiplicity(u));
    std::sort(fs.begin(), fs.end(), [](const auto& a, const auto& b) { return P::canonical_less(a.first, b.first); });
    marca(T.mult, F_MULT);
    std::vector<std::pair<std::vector<BigInt>, int>> r;
    for (auto& [u, e] : fs) r.emplace_back(coef_big(u), e);
    return r;
}

// Variante (observación lateral): DDF clásica que, tras cada cambio de
// módulo, compone con H_1 mod g en lugar de volver a elevar a q.  Mismos H_d.
template <class F>
std::vector<std::pair<int, std::vector<BigInt>>> ddf_compone(const Poly<F>& g0) {
    using P = Poly<F>;
    const F& K = g0.field();
    P gg = g0.monic(), x = P::variable(K);
    typename P::Reducer red(gg);
    P H1 = P::pow_mod(x, K.cardinality(), red);
    auto kk = [](const P& m) {
        const std::size_t n = static_cast<std::size_t>(std::max(1, m.degree()));
        return std::max<std::size_t>(1, std::min(n, (std::size_t{1} << 18) / n));
    };
    std::optional<typename P::Composer> comp;
    comp.emplace(H1, red, kk(gg));
    P H = H1;
    std::vector<std::pair<int, std::vector<BigInt>>> out;
    int d = 1;
    while (gg.degree() >= 2 * d) {
        if (d > 1) H = comp->compose(H);
        P hd = P::gcd(H - x, gg);
        if (hd.degree() > 0) {
            out.emplace_back(d, coef_big(hd.monic()));
            gg = gg / hd;
            H = H % gg;
            H1 = H1 % gg;
            red = typename P::Reducer(gg);
            comp.emplace(H1, red, kk(gg));
        }
        ++d;
    }
    if (gg.degree() > 0) out.emplace_back(gg.degree(), coef_big(gg.monic()));
    return out;
}

template <class F>
std::vector<std::pair<int, std::vector<BigInt>>> ddf_lib(const Poly<F>& g) {
    std::vector<std::pair<int, std::vector<BigInt>>> r;
    for (auto& [d, u] : g.distinct_degree()) r.emplace_back(d, coef_big(u));
    return r;
}

static void imprime_fases(const char* nom, const Fases& T, int llamadas) {
    const double tot = T.total();
    std::printf("  %-10s total %8.3f ms |", nom, tot / llamadas / 1000);
    const double v[F_N] = {T.rad, T.ddf_pow, T.ddf_comp, T.ddf_gcd, T.ddf_otro, T.edf_pow, T.edf_otro, T.mult};
    for (int k = 0; k < F_N; ++k) std::printf(" %s %.3f (%.0f%%) |", kNomFase[k], v[k] / llamadas / 1000, 100 * v[k] / tot);
    std::printf("\n             por llamada: %.2f pow_mod DDF, %.2f composiciones, %.2f MCD DDF, %.2f pow_mod EDF (%.2f divisiones)\n",
                double(T.n_pow) / llamadas, double(T.n_comp) / llamadas, double(T.n_gcd) / llamadas,
                double(T.n_edf_pow) / llamadas, double(T.n_split) / llamadas);
}

static void seccion_factor(const Fp& K, const FpN2& KN, const FpCont& KC) {
    using P = Poly<Fp>;
    using PN = Poly<FpN2>;
    using PC = Poly<FpCont>;
    constexpr int N = 16;
    set_seed(20261107);
    std::vector<P> f;
    for (int i = 0; i < N; ++i) f.push_back(P::random(K, 20).monic());
    std::vector<PN> fn;
    std::vector<PC> fc;
    for (int i = 0; i < N; ++i) {
        fn.push_back(convierte(KN, f[i]));
        fc.push_back(convierte(KC, f[i]));
    }
    // Comprobación: factor() de la librería con Fp, con FpN2 y el desglose.
    std::printf("\n== 3. factor() de 16 polinomios mónicos aleatorios de grado 20 en F_(2^127-1)\n");
    std::printf("  grados de los factores:");
    for (int i = 0; i < N; ++i) {
        set_seed(1000 + i);
        auto a = fact_big<Fp>(f[i].factor());
        set_seed(1000 + i);
        auto b = fact_big<FpN2>(fn[i].factor());
        comprueba(a == b, "factor Fp frente a FpN2");
        Fases T;
        set_seed(1000 + i);
        auto c = factor_desglosado(f[i], T, nullptr);
        comprueba(a == c, "factor frente al desglose");
        std::printf(" [");
        for (auto& [u, e] : a) std::printf("%zu%s", u.size() - 1, e > 1 ? "^" : "");
        std::printf("]");
        // DDF que compone desde el principio: mismo resultado que distinct_degree.
        comprueba(ddf_lib(f[i].radical()) == ddf_compone(f[i].radical()), "ddf_compone");
    }
    std::printf("\n");

    const double tf = medir_us(N, [&](int i) { set_seed(1000 + i); usa(f[i].factor().factors.size()); });
    const double tfn = medir_us(N, [&](int i) { set_seed(1000 + i); usa(fn[i].factor().factors.size()); });
    std::printf("  factor() total (us/llamada): Poly<Fp> %.1f   Poly<FpN2> %.1f   factor %.2fx\n", tf, tfn, tf / tfn);

    // Desglose: 3 pasadas por las 16 entradas.
    Fases T, TN;
    for (int rep = 0; rep < 3; ++rep)
        for (int i = 0; i < N; ++i) {
            set_seed(1000 + i);
            (void)factor_desglosado(f[i], T, nullptr);
            set_seed(1000 + i);
            (void)factor_desglosado(fn[i], TN, nullptr);
        }
    std::printf("  Desglose (ms por llamada, media de %d):\n", 3 * N);
    imprime_fases("Fp", T, 3 * N);
    imprime_fases("FpN2", TN, 3 * N);

    // Recuento de operaciones de elemento (FpCont) y de reservas, por fase.
    RecFase rec[F_N];
    Fases TC;
    for (int i = 0; i < N; ++i) {
        set_seed(1000 + i);
        (void)factor_desglosado(fc[i], TC, rec);
    }
    RecFase recN[F_N];
    Fases TN2;
    for (int i = 0; i < N; ++i) {
        set_seed(1000 + i);
        (void)factor_desglosado(fn[i], TN2, recN);
    }
    std::printf("  Recuento por llamada (media de %d; Fp contado con FpCont):\n", N);
    std::printf("  %-16s %8s %8s %6s %8s %7s %7s | %9s %8s | %9s | %9s %9s\n", "fase", "mul", "add/sub", "inv",
                "acc_fma", "reduce", "set", "res. GMP", "new Fp", "new FpN2", "medido ms", "arit.est");
    const double v[F_N] = {T.rad, T.ddf_pow, T.ddf_comp, T.ddf_gcd, T.ddf_otro, T.edf_pow, T.edf_otro, T.mult};
    RecFase tot, totN;
    double tot_est = 0;
    for (int k = 0; k < F_N; ++k) {
        const RecFase& r = rec[k];
        // Estimación de la aritmética de elementos con los costes de la sección 1
        // (inv con el coste de un elemento aleatorio: cota superior).
        const double est = (r.mul * g_cost.mul + r.addsub * g_cost.add + r.inv * g_cost.inv + r.fma * g_cost.fma +
                            r.red * g_cost.red + r.set * g_cost.set) / N / 1e6;
        tot_est += est;
        std::printf("  %-16s %8.0f %8.0f %6.0f %8.0f %7.0f %7.0f | %9.0f %8.0f | %9.0f | %9.3f %9.3f\n", kNomFase[k],
                    r.mul / double(N), r.addsub / double(N), r.inv / double(N), r.fma / double(N), r.red / double(N),
                    r.set / double(N), r.gmp / double(N), r.news / double(N), recN[k].news / double(N),
                    v[k] / (3 * N) / 1000, est);
        tot.mul += r.mul; tot.addsub += r.addsub; tot.inv += r.inv; tot.fma += r.fma; tot.red += r.red;
        tot.set += r.set; tot.gmp += r.gmp; tot.news += r.news; totN.news += recN[k].news; totN.gmp += recN[k].gmp;
    }
    std::printf("  %-16s %8.0f %8.0f %6.0f %8.0f %7.0f %7.0f | %9.0f %8.0f | %9.0f | %9.3f %9.3f\n", "TOTAL",
                tot.mul / double(N), tot.addsub / double(N), tot.inv / double(N), tot.fma / double(N),
                tot.red / double(N), tot.set / double(N), tot.gmp / double(N), tot.news / double(N),
                totN.news / double(N), T.total() / (3 * N) / 1000, tot_est);
    std::printf("  Reservas de GMP en Poly<FpN2>: %.0f por llamada.  Estimación de reservas en Poly<Fp>: %.0f x "
                "%.1f ns (malloc+free) = %.3f ms; operator new: %.0f x %.1f ns = %.3f ms\n",
                totN.gmp / double(N), tot.gmp / double(N), g_cost.malloc, tot.gmp / double(N) * g_cost.malloc / 1e6,
                tot.news / double(N), g_cost.malloc, tot.news / double(N) * g_cost.malloc / 1e6);

    // Variante de DDF: tiempo de distinct_degree frente a ddf_compone.
    std::vector<P> rad;
    for (int i = 0; i < N; ++i) rad.push_back(f[i].radical());
    std::vector<PN> radn;
    for (int i = 0; i < N; ++i) radn.push_back(fn[i].radical());
    const double d0 = medir_us(N, [&](int i) { usa(rad[i].distinct_degree().size()); });
    const double d1 = medir_us(N, [&](int i) { usa(ddf_compone(rad[i]).size()); });
    const double d0n = medir_us(N, [&](int i) { usa(radn[i].distinct_degree().size()); });
    const double d1n = medir_us(N, [&](int i) { usa(ddf_compone(radn[i]).size()); });
    std::printf("\n== 4. DDF (us/llamada): distinct_degree de la librería frente a componer desde el primer paso\n");
    std::printf("  Poly<Fp>:   librería %.1f   compone %.1f   (%.2fx)\n", d0, d1, d0 / d1);
    std::printf("  Poly<FpN2>: librería %.1f   compone %.1f   (%.2fx)\n", d0n, d1n, d0n / d1n);
}

// ---------------------------------------------------------------------------
// 5. La entrada exacta del escenario fp127.factor.d20 de compare_libs.py
//    (build/fase6/fpn/fp127_factor_d20.txt, generada con su misma semilla).
// ---------------------------------------------------------------------------
static void seccion_escenario(const Fp& K, const FpN2& KN) {
    using P = Poly<Fp>;
    using PN = Poly<FpN2>;
    std::FILE* fh = std::fopen("build/fase6/fpn/fp127_factor_d20.txt", "r");
    if (!fh) {
        std::printf("\n== 5. (no se encuentra build/fase6/fpn/fp127_factor_d20.txt; se omite)\n");
        return;
    }
    std::vector<BigInt> c;
    char buf[64];
    while (std::fscanf(fh, "%63s", buf) == 1) c.emplace_back(buf);
    std::fclose(fh);
    const P f(K, c);
    const PN fn = convierte(KN, f);
    set_seed(1000);
    const auto a = fact_big<Fp>(f.factor());
    set_seed(1000);
    const auto b = fact_big<FpN2>(fn.factor());
    comprueba(a == b, "escenario: factor Fp frente a FpN2");
    const P g = f.radical();
    const PN gn = fn.radical();
    comprueba(ddf_lib(g) == ddf_compone(g), "escenario: ddf_compone");
    std::printf("\n== 5. Escenario fp127.factor.d20 (la entrada de compare_libs.py; grados:");
    for (auto& [u, e] : a) std::printf(" %zu", u.size() - 1);
    std::printf(")\n");
    // Una sola entrada: aquí no se puede rotar; se repite la misma llamada.
    const double t0 = medir_us(1, [&](int) { set_seed(1000); usa(f.factor().factors.size()); });
    const double t1 = medir_us(1, [&](int) { set_seed(1000); usa(fn.factor().factors.size()); });
    const double d0 = medir_us(1, [&](int) { usa(g.distinct_degree().size()); });
    const double d1 = medir_us(1, [&](int) { usa(ddf_compone(g).size()); });
    const double d0n = medir_us(1, [&](int) { usa(gn.distinct_degree().size()); });
    const double d1n = medir_us(1, [&](int) { usa(ddf_compone(gn).size()); });
    std::printf("  factor() total (ms): Poly<Fp> %.2f   Poly<FpN2> %.2f   (%.2fx)\n", t0 / 1000, t1 / 1000, t0 / t1);
    std::printf("  DDF (ms): Poly<Fp> librería %.2f / compone %.2f;  Poly<FpN2> librería %.2f / compone %.2f\n",
                d0 / 1000, d1 / 1000, d0n / 1000, d1n / 1000);
    std::printf("  Estimación FpN2 + DDF que compone: %.2f ms (factor FpN2 - DDF librería + DDF compone)\n",
                (t1 - d0n + d1n) / 1000);
}

int main(int argc, char** argv) {
    g_pool = (argc > 1 && std::string(argv[1]) == "pool");
    mp_set_memory_functions(gmp_alloc, gmp_realloc, gmp_free);
    std::printf("fpn_bench (fase 6, tema fpn) — asignador de GMP: %s; ALCP_THREADS=%s\n",
                g_pool ? "listas libres (pool)" : "malloc", std::getenv("ALCP_THREADS") ? std::getenv("ALCP_THREADS") : "(sin definir)");
    {
        const BigInt p = (BigInt(1) << 127) - 1;
        Fp K(p);
        FpN2 KN(p);
        FpCont KC(p);
        set_seed(7);
        seccion_elementos(K, KN);
        seccion_polinomios(K, KN);
        seccion_factor(K, KN, KC);
        seccion_escenario(K, KN);
        // Otro primo de 2 palabras que no es de Mersenne (comprobación de la REDC
        // genérica): p = 2^89 - 1 es Mersenne; se usa el siguiente primo a 2^100.
        BigInt p2;
        const BigInt b100 = BigInt(1) << 100;
        mpz_nextprime(p2.get_mpz_t(), b100.get_mpz_t());
        Fp K2(p2);
        FpN2 KN2(p2);
        seccion_elementos(K2, KN2);
        BigInt p3;
        const BigInt b64 = BigInt(1) << 64;
        mpz_nextprime(p3.get_mpz_t(), b64.get_mpz_t());
        Fp K3(p3);
        FpN2 KN3(p3);
        seccion_elementos(K3, KN3);
    }
    std::printf("\nComprobaciones fallidas: %d   (suma de control %llu)\n", g_fallos,
                static_cast<unsigned long long>(g_chk & 0xffff));
    return g_fallos ? 1 : 0;
}
