# flint_gcd.py — Referencia: python-flint (FLINT) nmod_poly gcd / xgcd en F_998244353.
# Mismo protocolo que hgcd_proto.cpp: a de grado n, b de grado n-1 aleatorios,
# rotando 8 pares distintos, mejor de 3 tandas de >= 30 ms (us por llamada).
# Uso (desde la raíz del repo):  .venv/Scripts/python.exe bench/experimentos/fase6/hgcd/flint_gcd.py
import random
import time
import flint

P = 998244353
random.seed(1)


def aleatorio(n):
    c = [random.randrange(P) for _ in range(n)] + [random.randrange(1, P)]
    return flint.nmod_poly(c, P)


def medir_us(fn, ne):
    mejor = float("inf")
    for _ in range(3):
        veces = 0
        t0 = time.perf_counter()
        while True:
            fn(veces % ne)
            veces += 1
            dt = time.perf_counter() - t0
            if dt >= 0.03:
                break
        mejor = min(mejor, dt / veces * 1e6)
    return mejor


print("python-flint", flint.__version__, "- F_998244353, us por llamada (incluye la sobrecarga de Python)")
print("%6s %10s %10s %10s" % ("n", "gcd", "xgcd", "mul n*n"))
for n in (256, 384, 512, 768, 1000, 1024, 1536, 2048, 3072, 4096, 8192):
    ne = 8
    A = [aleatorio(n) for _ in range(ne)]
    B = [aleatorio(n - 1) for _ in range(ne)]
    tg = medir_us(lambda i: A[i].gcd(B[i]), ne)
    tx = medir_us(lambda i: A[i].xgcd(B[i]), ne)
    tm = medir_us(lambda i: A[i] * B[i], ne)
    print("%6d %10.1f %10.1f %10.1f" % (n, tg, tx, tm), flush=True)
