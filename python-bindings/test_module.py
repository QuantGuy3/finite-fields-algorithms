"""
End-to-end test of the compiled pybind11 module `alcp_native`.

Exercises F_p and F_q polynomial algebra, factorization round-trips, the
multiplication-algorithm dispatch, FFT/IFFT and Toeplitz routines, and spot-
checks agreement with the ORIGINAL reference Python (where the reference is
self-consistent).  Run:  python python/test_module.py
"""
import os, sys, random

sys.path.insert(0, os.path.dirname(__file__))
import alcp_native as al

fails = 0
def check(cond, msg):
    global fails
    if not cond:
        fails += 1
        print("  FAIL:", msg)

random.seed(7)

# ---- integers --------------------------------------------------------------
check(al.is_prime(998244353), "998244353 prime")
check(not al.is_prime(998244353 * 7), "composite")
check(al.prime_power_factor(3 ** 10) == (3, 10), "prime power 3^10")

# ---- F_p polynomial algebra + factorization round-trip ---------------------
F = al.Fp(101)
for _ in range(200):
    comp = F.poly([1])
    for _ in range(random.randint(1, 3)):
        d = random.randint(1, 3)
        e = random.randint(1, 3)
        base = F.poly([random.randrange(101) for _ in range(d)] + [random.randrange(1, 101)]).monic()
        comp = comp * (base ** e)
    if comp.degree() < 1:
        continue
    facs, lead = comp.factor()
    recon = F.poly([lead])
    for u, e in facs:
        check(u.is_irreducible(), "Fp factor irreducible")
        recon = recon * (u ** e)
    check(recon == comp, "Fp factorization reconstructs")

# multiplication algorithms agree
a = F.poly([random.randrange(101) for _ in range(40)])
b = F.poly([random.randrange(101) for _ in range(35)])
check(a.mul(b, "school") == a.mul(b, "karatsuba") == (a * b), "Fp mul algos agree")

# gcd / Bezout
a = F.poly([random.randrange(101) for _ in range(6)])
b = F.poly([random.randrange(101) for _ in range(4)])
g, x, y = a.gcd_ext(b)
check(g == x * a + y * b, "Fp Bezout identity")

# ---- F_q algebra -----------------------------------------------------------
Fq = al.Fq(3, [1, 0, 1])  # F_9 = F_3[x]/(x^2+1)
check(Fq.q == 9 and Fq.p == 3 and Fq.n == 2, "F_9 params")
for _ in range(100):
    comp = Fq.poly([1])
    for _ in range(random.randint(1, 2)):
        d = random.randint(1, 2)
        e = random.randint(1, 2)
        base = Fq.poly([random.randrange(9) for _ in range(d)] + [random.randrange(1, 9)]).monic()
        comp = comp * (base ** e)
    if comp.degree() < 1:
        continue
    facs, lead = comp.factor()
    recon = Fq.poly([lead])
    for u, e in facs:
        check(u.is_irreducible(), "Fq factor irreducible")
        recon = recon * (u ** e)
    check(recon == comp, "Fq factorization reconstructs")

# auto-built F_q from a prime power
Fq2 = al.Fq(27)
check(Fq2.q == 27, "auto F_27")

# ---- FFT / IFFT over F_17 (16 = 2^4) ---------------------------------------
F17 = al.Fp(17)
k, n, root = 4, 16, 3  # 3 has order 16 in F_17*
a = [random.randrange(17) for _ in range(n)]
A = al.dft(F17, root, k, a)
back = al.idft(F17, root, k, A)
check(back == a, "FFT/IFFT round-trip")

# ---- Toeplitz lower mat-vec vs naive ---------------------------------------
nn = 8
col = [random.randrange(101) for _ in range(nn)]
vec = [random.randrange(101) for _ in range(nn)]
y = al.toeplitz_lower_vec(F, nn, col, vec)
naive = [sum(col[i - j] * vec[j] for j in range(i + 1)) % 101 for i in range(nn)]
check(y == naive, "Toeplitz lower mat-vec")

# inverse: convolving col with its inverse gives e_0
col[0] = col[0] or 1
inv = al.toeplitz_lower_inverse(F, nn, col)
prod = al.toeplitz_lower_vec(F, nn, col, inv)
check(prod == [1] + [0] * (nn - 1), "Toeplitz lower inverse")

# ---- agreement with the ORIGINAL reference Python --------------------------
try:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "ref"))
    import cuerpos_finitos as cf
    fp = cf.cuerpo_fp(101); fpx = cf.anillo_fp_x(fp)
    for _ in range(50):
        fI = [random.randrange(101) for _ in range(5)] + [random.randrange(1, 101)]
        gI = [random.randrange(101) for _ in range(3)] + [random.randrange(1, 101)]
        ref = [int(c) for c in fpx.conv_a_tuple(fpx.mult(fpx.elem_de_tuple(tuple(fI)),
                                                         fpx.elem_de_tuple(tuple(gI))))]
        mine = (F.poly(fI) * F.poly(gI)).coeffs()
        check(mine == ref, "module mul == reference Python mul")
    print("  (reference-Python agreement checks ran)")
except Exception as ex:  # reference modules optional
    print("  (skipped reference comparison:", ex, ")")

if fails == 0:
    print("test_module: ALL PASS")
    sys.exit(0)
else:
    print(f"test_module: {fails} FAILURES")
    sys.exit(1)
