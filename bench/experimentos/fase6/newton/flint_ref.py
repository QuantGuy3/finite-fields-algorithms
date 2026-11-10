# flint_ref.py — Fase 6, tema «newton»: referencia python-flint (nmod_poly) en
# F_998244353 con los mismos tamaños que newton_proto.cpp.
#   * divmod grado 4000 / grado 2000
#   * pow_mod(a, q, m) con deg m = 300 y 1000, q = p, a de grado deg m - 1
#   * (si está el módulo de alcp) los mismos casos desde el binding de alcp, para
#     ver el estado «frío» (cuerpo recién creado, sin tablas FFT)
# Metodología: 16 entradas distintas rotando, mejor de 3 tandas de >= 30 ms.
# Ejecución (desde la raíz del repo):
#   ALCP_THREADS=1 .venv/Scripts/python.exe bench/experimentos/fase6/newton/flint_ref.py
import random
import time

import flint

P = 998244353
rng = random.Random(12345)


def aleatorio(d):
    c = [rng.randrange(P) for _ in range(d)] + [rng.randrange(1, P)]
    return c


def medir_us(fn):
    mejor = float("inf")
    for _ in range(3):
        veces = 0
        t0 = time.perf_counter()
        while True:
            fn(veces % 16)
            veces += 1
            dt = time.perf_counter() - t0
            if dt >= 0.03:
                break
        mejor = min(mejor, dt / veces)
    return mejor * 1e6


casos_div = [(aleatorio(4000), aleatorio(2000)) for _ in range(16)]
casos_pow = {n: [(aleatorio(n - 1), aleatorio(n)) for _ in range(16)] for n in (300, 1000)}

print("python-flint", flint.__version__, "— F_998244353, us (mejor de 3 tandas, 16 entradas)")
F = [(flint.nmod_poly(f, P), flint.nmod_poly(g, P)) for f, g in casos_div]
print(f"  flint divmod 4000/2000     : {medir_us(lambda i: divmod(F[i][0], F[i][1])):9.1f}")
for n, casos in casos_pow.items():
    A = [(flint.nmod_poly(a, P), flint.nmod_poly(m, P)) for a, m in casos]
    print(f"  flint pow_mod deg m = {n:4d} : {medir_us(lambda i: A[i][0].pow_mod(P, A[i][1])):9.1f}")
    # Con el inverso de rev(m) precalculado (lo que hace un Reducer)
    try:
        R = [m.reverse().inverse_series_trunc(n) for _, m in A]
        t = medir_us(lambda i: A[i][0].pow_mod(P, A[i][1], R[i]))
        print(f"  flint pow_mod deg m = {n:4d} (con mod_rev_inv precalculado): {t:9.1f}")
    except Exception as e:  # la firma puede variar entre versiones
        print("  (pow_mod con mod_rev_inv no disponible:", e, ")")

try:
    from alcp import alcp_native as nat

    K = nat.Fp(P)
    Fa = [(K.poly(f), K.poly(g)) for f, g in casos_div]
    print(f"  alcp (binding) divmod 4000/2000, cuerpo frío: {medir_us(lambda i: Fa[i][0].divmod(Fa[i][1])):9.1f}")
    (Fa[0][0] * Fa[0][1])  # construye tablas FFT (pide la raíz)
    print(f"  alcp (binding) divmod 4000/2000, cuerpo caliente: {medir_us(lambda i: Fa[i][0].divmod(Fa[i][1])):9.1f}")
    for n, casos in casos_pow.items():
        A = [(K.poly(a), K.poly(m)) for a, m in casos]
        if hasattr(A[0][0], "pow_mod"):
            print(f"  alcp (binding) pow_mod deg m = {n:4d}: {medir_us(lambda i: A[i][0].pow_mod(P, A[i][1])):9.1f}")
except Exception as e:
    print("  (binding de alcp no disponible:", e, ")")
