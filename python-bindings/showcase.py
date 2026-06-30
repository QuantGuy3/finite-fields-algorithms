#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
================================================================================
 showcase.py  --  Recorrido pedagógico y expositivo por la librería alcp_native
================================================================================

Este archivo NO es un test (para eso está test_module.py); es una demostración
COMENTADA, pensada para leerse de arriba a abajo.  Recorre:

  Parte 1 .  Las funciones de la librería, con su significado matemático.
  Parte 2 .  Tiempos de multiplicación: schoolbook vs Karatsuba vs FFT (NTT).
  Parte 3 .  Comparación con otras implementaciones (numpy, scipy, sympy)
             y, comentado, el equivalente en Julia/Nemo.
  Parte 4 .  Un caso donde la implementación Python original fallaba.

Ejecutar:   python python/showcase.py
"""
import os, sys, time, random, shutil, subprocess

sys.path.insert(0, os.path.dirname(__file__))
import alcp_native as al

# Salida UTF-8 en la consola de Windows (para que los acentos salgan bien).
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

random.seed(0)

# ---------------------------------------------------------------- utilidades --
def title(s):
    print("\n" + "=" * 78)
    print("  " + s)
    print("=" * 78)

def sub(s):
    print("\n--- " + s + " " + "-" * max(0, 70 - len(s)))

def bench(fn, repeat=1):
    """Devuelve el mejor tiempo en milisegundos sobre `repeat` ejecuciones."""
    best = float("inf")
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best * 1000.0


# =============================================================================
title("PARTE 1 — Recorrido por la librería")
# =============================================================================

sub("1.1  Enteros: primalidad (Miller–Rabin de GMP) y potencias de primo")
# is_prime usa el test probabilístico de GMP (error < 4^-40): instantáneo
# incluso para enteros gigantes, a diferencia de la criba O(sqrt(p)) original.
print("is_prime(2**127 - 1)          =", al.is_prime(2**127 - 1))   # primo de Mersenne
print("is_prime(2**127 - 3)          =", al.is_prime(2**127 - 3))
# prime_power_factor(q) -> (p, n) con q = p^n  (lanza error si no es potencia)
print("prime_power_factor(3**30)     =", al.prime_power_factor(3**30))
print("prime_power_factor(1024)      =", al.prime_power_factor(1024))

sub("1.2  Cuerpo primo F_p  (los elementos son enteros en [0, p))")
F = al.Fp(7)
print("F = Fp(7);  p =", F.p)
# La aritmética de cuerpo se hace con métodos (los elementos son int normales):
print("3 + 5 (mod 7)  =", F.add(3, 5))
print("3 * 5 (mod 7)  =", F.mul(3, 5))
print("inv(3) en F_7  =", F.inv(3), " ->  3 *", F.inv(3), "=", F.mul(3, F.inv(3)))
print("3^6 (mod 7)    =", F.pow(3, 6), " (pequeño Fermat: a^(p-1)=1)")

sub("1.3  Polinomios F_p[x] con SOBRECARGA DE OPERADORES")
# Los coeficientes se dan en orden ascendente (little-endian): [c0, c1, c2,...]
f = F.poly([1, 0, 1])   # 1 + x^2
g = F.poly([2, 1])      # 2 + x
print("f =", f.to_string(), "   g =", g.to_string())
print("f + g  =>", (f + g).coeffs(), "      (operador __add__)")
print("f * g  =>", (f * g).coeffs(), "  (operador __mul__, auto-despacha el algoritmo)")
print("f**3   =>", (f ** 3).coeffs())
q, r = (f * g).divmod(g)
print("(f*g) // g =>", q.coeffs(), "   (f*g) %% g =>", r.coeffs(), "  -> recupera f, resto 0")

sub("1.4  MCD, identidad de Bézout e inverso modular")
a = F.poly([6, 0, 1])   # x^2 - 1 = (x-1)(x+1)
b = F.poly([6, 1])      # x - 1
print("gcd(x^2-1, x-1) =", a.gcd(b).coeffs(), " (mónico)")
gg, xx, yy = a.gcd_ext(b)               # g = x*a + y*b
print("Bézout: g =", gg.coeffs(), " x =", xx.coeffs(), " y =", yy.coeffs(),
      " ->  x*a+y*b == g:", (xx * a + yy * b) == gg)

sub("1.5  Irreducibilidad (test de Rabin) y factorización (Cantor–Zassenhaus)")
print("x^2+1 irreducible sobre F_7?  ", F.poly([1, 0, 1]).is_irreducible(),
      " (sí: -1 no es residuo cuadrático mod 7)")
print("x^2-1 irreducible sobre F_7?  ", F.poly([6, 0, 1]).is_irreducible(),
      " (no: = (x-1)(x+1))")
comp = F.poly([6, 1]) ** 2 * F.poly([1, 0, 1])     # (x-1)^2 (x^2+1)
factores, lider = comp.factor()
print("factor((x-1)^2 (x^2+1)):  coef.líder =", lider)
for u, e in factores:
    print("     (", u.to_string(), ")^", e, "   irreducible:", u.is_irreducible())
# Comprobación de reconstrucción (lo que el Python original NO comprobaba):
recon = F.poly([lider])
for u, e in factores:
    recon = recon * u ** e
print("   reconstruye el original?", recon == comp)

sub("1.6  Cuerpo de extensión F_q = F_p[x]/(g)")
# Coeficientes de un PolyFq = índices base-p del elemento de F_q en [0, q).
Fq = al.Fq(3, [1, 0, 1])        # F_9 = F_3[x]/(x^2+1);  alternativa: al.Fq(9)
print("F_9: p =", Fq.p, " n =", Fq.n, " q =", Fq.q, " modulo =", Fq.modulus)
h = Fq.poly([1, 2, 5]) * Fq.poly([3, 1])
print("producto sobre F_9[x] -> coef (índices):", h.coeffs())
facs9, lead9 = (Fq.poly([1, 1]) ** 2 * Fq.poly([2, 0, 1])).factor()
print("factorización sobre F_9:", [(u.coeffs(), e) for u, e in facs9])

sub("1.7  Transformada de Fourier (Cooley–Tukey) con raíces precomputadas")
# DFT_{n,g}(a) con n = 2^k y g raíz 2^k-ésima de la unidad en F_p.
F17 = al.Fp(17)                  # 16 = 2^4 divide a 17-1, raíz primitiva g = 3
a = [random.randrange(17) for _ in range(16)]
A = al.dft(F17, 3, 4, a)         # transformada directa
back = al.idft(F17, 3, 4, A)     # inversa: n^{-1} * DFT_{g^{-1}}
print("a            =", a)
print("idft(dft(a)) =", back, " ->  ¿recupera a?", back == a)

sub("1.8  Productos Toeplitz e inversa triangular (vía multiplicación rápida)")
n = 6
col = [random.randrange(1, 101) for _ in range(n)]   # 1ª columna (Toeplitz inf.)
vec = [random.randrange(101) for _ in range(n)]
F101 = al.Fp(101)
y = al.toeplitz_lower_vec(F101, n, col, vec)
inv = al.toeplitz_lower_inverse(F101, n, col)         # 1ª columna de T^{-1}
chk = al.toeplitz_lower_vec(F101, n, col, inv)        # T * T^{-1}e_0 = e_0
print("T * vec      =", y)
print("T * inv(T)   =", chk, " ->  = e_0?", chk == [1] + [0] * (n - 1))


# =============================================================================
title("PARTE 2 — Tiempos:  schoolbook  vs  Karatsuba  vs  FFT (NTT)")
# =============================================================================
# Multiplicamos sobre F_998244353, un primo "NTT-friendly" (2^23 | p-1), de modo
# que la FFT exacta (NTT) es aplicable.  '.mul(other, algo)' fuerza el método;
# 'a * b' elige solo según el tamaño.
P = al.Fp(998244353)
print("\n  grado |  schoolbook |  Karatsuba |     FFT  |  auto (a*b elige)")
print("  ------+-------------+------------+----------+------------------")
for d in [256, 512, 1024, 2048]:
    A = P.poly([random.randrange(998244353) for _ in range(d + 1)])
    B = P.poly([random.randrange(998244353) for _ in range(d + 1)])
    t_s = bench(lambda: A.mul(B, "school")) if d <= 2048 else float("nan")
    t_k = bench(lambda: A.mul(B, "karatsuba"))
    t_f = bench(lambda: A.mul(B, "fft"))
    t_a = bench(lambda: A * B)
    # verificación de que los tres coinciden:
    assert A.mul(B, "school") == A.mul(B, "karatsuba") == A.mul(B, "fft") == (A * B)
    print(f"  {d:5d} | {t_s:9.1f}ms | {t_k:8.1f}ms | {t_f:6.1f}ms | {t_a:6.1f}ms")
print("\n  Lectura: schoolbook es O(n^2); Karatsuba O(n^1.585); FFT O(n log n).")
print("  El cruce Karatsuba->FFT ocurre pronto; por eso 'auto' elige FFT en grados altos.")


# =============================================================================
title("PARTE 3 — Comparación con otras implementaciones")
# =============================================================================

sub("3.1  Multiplicación con primo PEQUEÑO (p = 10007):  alcp vs numpy vs scipy")
# Con un primo pequeño, los coeficientes caben en una palabra de máquina, así que
# la convolución entera de numpy (int64) y la FFT en float de scipy son EXACTAS.
try:
    import numpy as np
    from scipy.signal import fftconvolve

    p = 10007
    Fp = al.Fp(p)
    d = 1500
    aa = [random.randrange(p) for _ in range(d + 1)]
    bb = [random.randrange(p) for _ in range(d + 1)]

    def alcp_mul():
        return (Fp.poly(aa) * Fp.poly(bb)).coeffs()
    def numpy_mul():
        r = np.convolve(np.array(aa, dtype=np.int64), np.array(bb, dtype=np.int64))
        return [int(v) % p for v in r]
    def scipy_mul():
        r = fftconvolve(np.array(aa, float), np.array(bb, float))
        return [int(round(v)) % p for v in r]

    print("¿coinciden alcp == numpy == scipy?  ", alcp_mul() == numpy_mul() == scipy_mul())
    print(f"   alcp  (GMP BigInt por coef.) : {bench(alcp_mul, 3):7.1f} ms")
    print(f"   numpy (int64 nativo, en C)   : {bench(numpy_mul, 3):7.1f} ms")
    print(f"   scipy (FFT float64)          : {bench(scipy_mul, 3):7.1f} ms")
    print("   Lección honesta: para primos pequeños numpy/scipy ganan de calle porque")
    print("   usan enteros/floats NATIVOS; alcp paga el coste de un mpz de GMP por")
    print("   coeficiente.  La ventaja de alcp aparece cuando eso ya no cabe -> 3.1bis.")

    sub("3.1bis  Primo GRANDE (p = 998244353):  cuando numpy/scipy se EQUIVOCAN")
    # Coeficientes ~10^9; sus productos ~10^18 y las sumas desbordan int64 (numpy)
    # y exceden la mantisa de 53 bits (scipy float).  alcp es exacto siempre.
    pbig = 998244353
    Fbig = al.Fp(pbig)
    d = 600
    aa = [random.randrange(pbig) for _ in range(d + 1)]
    bb = [random.randrange(pbig) for _ in range(d + 1)]
    # referencia EXACTA con enteros grandes de Python:
    exact = [0] * (2 * d + 1)
    for i, ai in enumerate(aa):
        for j, bj in enumerate(bb):
            exact[i + j] += ai * bj
    exact = [v % pbig for v in exact]
    r_alcp = (Fbig.poly(aa) * Fbig.poly(bb)).coeffs()
    r_alcp += [0] * (len(exact) - len(r_alcp))
    r_np = [int(v) % pbig for v in np.convolve(np.array(aa, dtype=np.int64), np.array(bb, dtype=np.int64))]
    r_sp = [int(round(v)) % pbig for v in fftconvolve(np.array(aa, float), np.array(bb, float))]
    bad_np = sum(1 for x, y in zip(r_np, exact) if x != y)
    bad_sp = sum(1 for x, y in zip(r_sp, exact) if x != y)
    bad_al = sum(1 for x, y in zip(r_alcp, exact) if x != y)
    print(f"   coeficientes ERRÓNEOS frente a la referencia exacta (de {len(exact)}):")
    print(f"      numpy int64   : {bad_np:5d}   (DESBORDA la palabra de 64 bits)")
    print(f"      scipy float64 : {bad_sp:5d}   (excede la mantisa de 53 bits)")
    print(f"      alcp (GMP)    : {bad_al:5d}   <- exacto: aquí está su razón de ser")
    print("   Moraleja: la NTT exacta de alcp exige un primo con 2^k | p-1; a cambio")
    print("   da el resultado correcto donde las bibliotecas numéricas fallan en silencio.")
except Exception as ex:
    print("   (numpy/scipy no disponibles:", ex, ")")

sub("3.2  Factorización en F_p:  alcp  vs  sympy")
# sympy factoriza sobre GF(p) con su propio algoritmo (Zassenhaus/Shoup).
try:
    from sympy import symbols, Poly
    x = symbols("x")
    p = 101
    Fp = al.Fp(p)

    # construimos un compuesto conocido de grado ~30
    comp = Fp.poly([1])
    for _ in range(4):
        deg = random.randint(2, 4)
        e = random.randint(1, 2)
        base = Fp.poly([random.randrange(p) for _ in range(deg)] + [random.randrange(1, p)]).monic()
        comp = comp * base ** e
    coeffs_le = comp.coeffs()

    def norm_alcp():
        facs, _ = Fp.poly(coeffs_le).factor()
        return sorted((tuple(u.monic().coeffs()), e) for u, e in facs)

    def norm_sympy():
        expr = sum(int(c) * x ** i for i, c in enumerate(coeffs_le))
        _, facs = Poly(expr, x, modulus=p).factor_list()
        out = []
        for poly, mult in facs:
            c = [int(v) % p for v in reversed(poly.all_coeffs())]   # a little-endian
            il = pow(c[-1], p - 2, p)                                # monizar
            out.append((tuple((v * il) % p for v in c), mult))
        return sorted(out)

    fa, fs = norm_alcp(), norm_sympy()
    print("grado del polinomio:", comp.degree())
    print("¿misma factorización (alcp == sympy)?  ", fa == fs)
    print(f"   alcp  : {bench(norm_alcp, 3):7.2f} ms")
    print(f"   sympy : {bench(norm_sympy, 3):7.2f} ms")
except Exception as ex:
    print("   (sympy no disponible:", ex, ")")

sub("3.3  Julia / Nemo (FLINT) — la referencia de oro, comparación REAL")
# Nemo.jl envuelve FLINT (C, sobre GMP): décadas de optimización; es el patrón
# contra el que medirse.  Lanzamos un proceso Julia con LOS MISMOS datos.
def find_julia():
    j = shutil.which("julia")
    if j:
        return j
    cand = os.path.join(os.environ.get("LOCALAPPDATA", ""), "Microsoft", "WindowsApps", "julia.exe")
    return cand if os.path.exists(cand) else None

def parse_factor_set(s):
    out = set()
    for item in filter(None, s.split(";")):
        cs, _, m = item.partition(":")
        out.add((tuple(int(x) for x in cs.split(",")) if cs else (), int(m)))
    return out

try:
    jpath = find_julia()
    if not jpath:
        raise RuntimeError("julia no encontrado (instalar: winget install 9NJNWW8PVKMN)")
    here = os.path.dirname(os.path.abspath(__file__))

    pmul = 998244353
    degm = 2000
    A = [random.randrange(pmul) for _ in range(degm + 1)]
    B = [random.randrange(pmul) for _ in range(degm + 1)]
    Pm = al.Fp(pmul)

    pfac = 101
    Fp101 = al.Fp(pfac)
    comp = Fp101.poly([1])
    for _ in range(5):
        d = random.randint(2, 4); e = random.randint(1, 2)
        base = Fp101.poly([random.randrange(pfac) for _ in range(d)] + [random.randrange(1, pfac)]).monic()
        comp = comp * base ** e
    Fco = comp.coeffs()

    with open(os.path.join(here, "nemo_in.txt"), "w") as fh:
        fh.write(f"{pmul}\n{','.join(map(str,A))}\n{','.join(map(str,B))}\n"
                 f"{pfac}\n{','.join(map(str,Fco))}\n")
    print("   (lanzando Julia/Nemo... el arranque + carga de FLINT tarda unos segundos)")
    subprocess.run([jpath, os.path.join(here, "compare_nemo.jl"),
                    os.path.join(here, "nemo_in.txt"), os.path.join(here, "nemo_out.txt")],
                   check=True, capture_output=True, text=True)
    nemo = {}
    for line in open(os.path.join(here, "nemo_out.txt"), encoding="utf-8"):
        k, _, v = line.partition(" ")
        nemo[k.strip()] = v.strip()

    # multiplicación
    al_prod = (Pm.poly(A) * Pm.poly(B)).coeffs()
    nemo_prod = [int(x) for x in nemo["MULT"].split(",")]
    t_al = bench(lambda: Pm.poly(A) * Pm.poly(B), 3)
    print(f"\n   Multiplicación grado-{degm} sobre F_998244353  (ambos exactos, NTT):")
    print(f"      ¿alcp == Nemo?  {al_prod == nemo_prod}")
    print(f"      alcp (NTT propia)  : {t_al:8.2f} ms")
    print(f"      Nemo/FLINT         : {float(nemo['MULT_MS']):8.2f} ms")

    # factorización
    al_set = set((tuple(u.monic().coeffs()), e) for u, e in Fp101.poly(Fco).factor()[0])
    nemo_set = parse_factor_set(nemo["FACTOR"])
    t_alf = bench(lambda: Fp101.poly(Fco).factor(), 3)
    print(f"\n   Factorización grado-{comp.degree()} sobre F_101:")
    print(f"      ¿mismo conjunto de factores (alcp == Nemo)?  {al_set == nemo_set}")
    print(f"      alcp               : {t_alf:8.2f} ms")
    print(f"      Nemo/FLINT         : {float(nemo['FACTOR_MS']):8.2f} ms")
    print("\n   Interpretación: FLINT es la implementación de referencia (C, muy")
    print("   optimizada); que alcp dé EXACTAMENTE el mismo resultado valida el port,")
    print("   y los tiempos lo sitúan en el orden de magnitud correcto para una")
    print("   implementación didáctica.")
except Exception as ex:
    print("   (Nemo/Julia no disponible:", ex, ")")
    print(r"""   # Equivalente conceptual en Nemo.jl:
   #   using Nemo
   #   R, x = polynomial_ring(Nemo.Native.GF(101), "x")
   #   factor(R([6,1])^2 * R([1,0,1]))""")


# =============================================================================
title("PARTE 4 — Un caso donde la implementación Python ORIGINAL fallaba")
# =============================================================================
# f = x^6 (x+1)^3 sobre F_2.  Los factorizadores Python originales devolvían un
# resultado que NI SIQUIERA reconstruye f; alcp da el correcto.
F2 = al.Fp(2)
f = F2.poly([0, 0, 0, 0, 0, 0, 1, 1, 1, 1])     # x^6 + x^7 + x^8 + x^9
facs, lead = f.factor()
print("f = x^6 (x+1)^3 sobre F_2")
print("alcp:", [(u.coeffs(), e) for u, e in facs], " líder =", lead)
recon = F2.poly([lead])
for u, e in facs:
    recon = recon * u ** e
print("¿alcp reconstruye f?", recon == f, "  (el Python original NO: ver README)")

print("\n" + "=" * 78)
print("  Fin del recorrido.  Núcleo en C++/GMP; aquí sólo se ha llamado a Python.")
print("=" * 78)
