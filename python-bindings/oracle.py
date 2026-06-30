"""
Cross-check oracle.

Runs the ORIGINAL Python implementations (cuerpos_finitos.py, factorizacion.py,
algoritmos_rapidos.py — Prácticas 1/2/3) on seeded pseudo-random inputs and
emits, one record per line, both the inputs and the expected outputs.  The C++
program tests/crosscheck.cpp reads this file, recomputes every result with the
C++ port, and verifies it matches the Python reference exactly (for randomized
factorization it compares the canonical factor multiset).

Record grammar (fields separated by '|'):
    OP | field | in... | out...
where
    field  : "p"            for F_p
             "p;g0,g1,..."  for F_q = F_p[x]/(g), g little-endian monic
    poly   : comma-separated coefficients (little-endian); F_p coeffs are
             integers, F_q coeffs are base-p indices (== conv_a_int)
    factor : "coeffs:mult;coeffs:mult;..." sorted canonically
"""
import sys, os, random

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "ref"))
import cuerpos_finitos as cf
import factorizacion          # patches anillo_*.factorizar (Cantor-Zassenhaus)
import algoritmos_rapidos as ar  # patches mult_fast / divmod_fast, FFT, Toeplitz

random.seed(20260629)
OUT = []


def emit(*tokens):
    OUT.append("|".join(str(t) for t in tokens))


# ---- F_p helpers -----------------------------------------------------------
def poly_fp(fpx, ints):
    return fpx.elem_de_tuple(tuple(ints))

def coeffs_fp(fpx, poly):
    return [int(c) for c in fpx.conv_a_tuple(poly)]

def s_poly(ints):
    return ",".join(str(i) for i in ints) if ints else "0"

def rand_poly_ints(p, deg):
    if deg < 0:
        return [0]
    c = [random.randrange(p) for _ in range(deg)]
    c.append(random.randrange(1, p))   # nonzero leading
    return c

def canon_factors_fp(fpx, faclist):
    items = []
    for (u, e) in faclist:
        items.append((tuple(coeffs_fp(fpx, u)), e))
    items.sort()
    return ";".join(",".join(str(x) for x in cs) + ":" + str(e) for (cs, e) in items)


# ---- F_q helpers -----------------------------------------------------------
def poly_fq(fqx, fq, idxs):
    return fqx.elem_de_tuple(tuple(fq.elem_de_int(i) for i in idxs))

def coeffs_fq(fqx, fq, poly):
    return [fq.conv_a_int(c) for c in fqx.conv_a_tuple(poly)]

def rand_poly_fq_idxs(q, deg):
    if deg < 0:
        return [0]
    c = [random.randrange(q) for _ in range(deg)]
    c.append(random.randrange(1, q))
    return c

def canon_factors_fq(fqx, fq, faclist):
    items = []
    for (u, e) in faclist:
        items.append((tuple(coeffs_fq(fqx, fq, u)), e))
    items.sort()
    return ";".join(",".join(str(x) for x in cs) + ":" + str(e) for (cs, e) in items)


# ===========================================================================
# F_p[x] cases
# ===========================================================================
FP_PRIMES = [2, 3, 5, 7, 13, 101]
for p in FP_PRIMES:
    fp = cf.cuerpo_fp(p)
    fpx = cf.anillo_fp_x(fp)
    for _ in range(12):
        fI = rand_poly_ints(p, random.randint(0, 7))
        gI = rand_poly_ints(p, random.randint(0, 5))
        f, g = poly_fp(fpx, fI), poly_fp(fpx, gI)

        emit("FPX_MUL", p, s_poly(fI), s_poly(gI), s_poly(coeffs_fp(fpx, fpx.mult(f, g))))
        emit("FPX_KARA", p, s_poly(fI), s_poly(gI), s_poly(coeffs_fp(fpx, fpx.mult_fast(f, g))))

        if cf.PolyK(fpx, gI).grado >= 0 and any(gI):
            q, r = fpx.divmod(f, g)
            emit("FPX_DIVMOD", p, s_poly(fI), s_poly(gI),
                 s_poly(coeffs_fp(fpx, q)), s_poly(coeffs_fp(fpx, r)))
            q2, r2 = fpx.divmod_fast(f, g)
            emit("FPX_DIVMOD_FAST", p, s_poly(fI), s_poly(gI),
                 s_poly(coeffs_fp(fpx, q2)), s_poly(coeffs_fp(fpx, r2)))

        emit("FPX_GCD", p, s_poly(fI), s_poly(gI), s_poly(coeffs_fp(fpx, fpx.gcd(f, g))))
        gg, xx, yy = fpx.gcd_ext(f, g)
        emit("FPX_GCDEXT", p, s_poly(fI), s_poly(gI),
             s_poly(coeffs_fp(fpx, gg)), s_poly(coeffs_fp(fpx, xx)), s_poly(coeffs_fp(fpx, yy)))

    # irreducibility + factorization
    for _ in range(10):
        fI = rand_poly_ints(p, random.randint(1, 5))
        f = poly_fp(fpx, fI)
        emit("FPX_IRRED", p, s_poly(fI), 1 if f.es_irreducible else 0)
    for _ in range(8):
        # build a known composite from random factors
        comp = poly_fp(fpx, [1])
        for _ in range(random.randint(1, 3)):
            d = random.randint(1, 3)
            e = random.randint(1, 2)
            base = poly_fp(fpx, rand_poly_ints(p, d))
            comp = fpx.mult(comp, base ** e)
        cI = coeffs_fp(fpx, comp)
        if fpx.grado(comp) >= 1:
            faclist = fpx.factorizar(comp)
            rec = poly_fp(fpx, [1])
            for (u, e) in faclist:
                rec = fpx.mult(rec, u ** e)
            py_ok = 1 if (coeffs_fp(fpx, rec) == cI and
                          all(u.es_irreducible for (u, _) in faclist)) else 0
            emit("FPX_FACTOR", p, s_poly(cI), canon_factors_fp(fpx, faclist), py_ok)

# pot_mod over a fixed prime
fp = cf.cuerpo_fp(7); fpx = cf.anillo_fp_x(fp)
for _ in range(8):
    aI = rand_poly_ints(7, random.randint(1, 4))
    mI = rand_poly_ints(7, random.randint(1, 3))
    k = random.randint(0, 25)
    a, m = poly_fp(fpx, aI), poly_fp(fpx, mI)
    if fpx.grado(m) >= 1:
        emit("FPX_POWMOD", 7, s_poly(aI), s_poly(mI), k, s_poly(coeffs_fp(fpx, fpx.pot_mod(a, k, m))))

# DFT / IDFT over F_17 (16 = 2^4)
fp17 = cf.cuerpo_fp(17); k = 4; n = 1 << k; groot = 3   # 3 has order 16 in F_17*
for _ in range(6):
    aI = [random.randrange(17) for _ in range(n)]
    a_elems = tuple(fp17.elem(v) for v in aI)
    A = ar.fp_fft(fp17, fp17.elem(groot), k, a_elems)
    emit("FP_DFT", 17, groot, k, s_poly(aI), s_poly([int(v) for v in A]))
    B = ar.fp_ifft(fp17, fp17.elem(groot), k, A)
    emit("FP_IDFT", 17, groot, k, s_poly([int(v) for v in A]), s_poly([int(v) for v in B]))

# Toeplitz lower mat-vec over F_101
fp101 = cf.cuerpo_fp(101)
for _ in range(6):
    nn = random.randint(2, 10)
    a = tuple(fp101.elem(random.randrange(101)) for _ in range(nn))
    b = tuple(fp101.elem(random.randrange(101)) for _ in range(nn))
    y = ar.fp_toep_inf_vec(fp101, nn, a, b)
    emit("FP_TOEP_LOWER", 101, nn, s_poly([int(v) for v in a]),
         s_poly([int(v) for v in b]), s_poly([int(v) for v in y]))


# ===========================================================================
# F_q[x] cases (fixed irreducible moduli shared with the C++ side)
# ===========================================================================
FQ_FIELDS = [
    (3, [1, 0, 1]),       # F_9  = F_3[x]/(x^2+1)
    (3, [1, 2, 0, 1]),    # F_27 = F_3[x]/(x^3+2x+1)
    (5, [2, 0, 1]),       # F_25 = F_5[x]/(x^2+2)
]
for (p, gmod) in FQ_FIELDS:
    fp = cf.cuerpo_fp(p)
    fq = cf.cuerpo_fq(fp, tuple(gmod))
    fqx = cf.anillo_fq_x(fq)
    q = fq.fq
    fld = f"{p};{','.join(str(c) for c in gmod)}"
    for _ in range(10):
        fI = rand_poly_fq_idxs(q, random.randint(0, 5))
        gI = rand_poly_fq_idxs(q, random.randint(0, 4))
        f, g = poly_fq(fqx, fq, fI), poly_fq(fqx, fq, gI)
        emit("FQX_MUL", fld, s_poly(fI), s_poly(gI), s_poly(coeffs_fq(fqx, fq, fqx.mult(f, g))))
        if fqx.grado(g) >= 0 and any(gI):
            qq, rr = fqx.divmod(f, g)
            emit("FQX_DIVMOD", fld, s_poly(fI), s_poly(gI),
                 s_poly(coeffs_fq(fqx, fq, qq)), s_poly(coeffs_fq(fqx, fq, rr)))
        emit("FQX_GCD", fld, s_poly(fI), s_poly(gI), s_poly(coeffs_fq(fqx, fq, fqx.gcd(f, g))))
    for _ in range(6):
        fI = rand_poly_fq_idxs(q, random.randint(1, 4))
        f = poly_fq(fqx, fq, fI)
        emit("FQX_IRRED", fld, s_poly(fI), 1 if f.es_irreducible else 0)
    for _ in range(5):
        comp = poly_fq(fqx, fq, [1])
        for _ in range(random.randint(1, 2)):
            d = random.randint(1, 2)
            e = random.randint(1, 2)
            base = poly_fq(fqx, fq, rand_poly_fq_idxs(q, d))
            comp = fqx.mult(comp, base ** e)
        if fqx.grado(comp) >= 1:
            cI = coeffs_fq(fqx, fq, comp)
            faclist = fqx.factorizar(comp)
            rec = poly_fq(fqx, fq, [1])
            for (u, e) in faclist:
                rec = fqx.mult(rec, u ** e)
            py_ok = 1 if (coeffs_fq(fqx, fq, rec) == cI and
                          all(u.es_irreducible for (u, _) in faclist)) else 0
            emit("FQX_FACTOR", fld, s_poly(cI), canon_factors_fq(fqx, fq, faclist), py_ok)


out_path = sys.argv[1] if len(sys.argv) > 1 else "oracle_vectors.txt"
with open(out_path, "w", encoding="utf-8") as fh:
    fh.write("\n".join(OUT) + "\n")
print(f"wrote {len(OUT)} cross-check records to {out_path}")
