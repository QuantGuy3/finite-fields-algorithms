"""flint_ref.py — Fase 6, tema «factor»: referencia python-flint con entradas VARIADAS.

Genera 16 entradas deterministas por caso (random.Random con semilla fija), las
guarda en bench/results/fase6/factor_entradas.txt (para que los prototipos C++
midan EXACTAMENTE las mismas entradas) y mide python-flint con la misma
metodología que los prototipos: rotando las 16 entradas, mejor de 3 tandas de
>= 30 ms.  Un hilo (python-flint no usa hilos en estas operaciones).

Formato del fichero: una línea por polinomio, «caso índice c0 c1 ... cn», con
los coeficientes como índices en base p en orden creciente (en F_q, el
elemento sum d_i a^i tiene índice sum d_i p^i).  En los productos, la entrada i
son las líneas 2i (a) y 2i+1 (b).

Uso (desde la raíz del repo):
    .venv/Scripts/python.exe bench/experimentos/fase6/factor/flint_ref.py
"""
from __future__ import annotations

import random
import time
from pathlib import Path

import flint

REPO = Path(__file__).resolve().parents[4]
SALIDA = REPO / "bench" / "results" / "fase6" / "factor_entradas.txt"
MOD_2_8 = [1, 1, 0, 1, 1, 0, 0, 0, 1]
MOD_3_30 = [2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1]
N = 16


def aleatorio(rng, q, grado, monico=False):
    c = [rng.randrange(q) for _ in range(grado)]
    c.append(1 if monico else rng.randrange(1, q))
    return c


def digitos(i, p, n):
    d = []
    for _ in range(n):
        i, r = divmod(i, p)
        d.append(r)
    return d


def medir(fn, n):
    """Mejor de 3 tandas de >= 30 ms rotando n entradas; µs por llamada."""
    mejor = float("inf")
    for _ in range(3):
        veces, t0 = 0, time.perf_counter()
        while True:
            fn(veces % n)
            veces += 1
            dt = time.perf_counter() - t0
            if dt >= 0.03:
                break
        mejor = min(mejor, dt / veces)
    return mejor * 1e6


def main():
    rng = random.Random(20261107)
    lineas = []
    res = []

    # ---- F_101 (nmod_poly) ----
    p = 101
    d40 = []
    for _ in range(N):
        f = flint.nmod_poly([1], p)
        for d, m in ((1, 2), (2, 1), (3, 2), (5, 1), (8, 1), (15, 1)):
            u = flint.nmod_poly(aleatorio(rng, p, d, True), p)
            f = f * u ** m
        d40.append(f)
    d100 = [flint.nmod_poly(aleatorio(rng, p, 100, True), p) for _ in range(N)]
    for nombre, lst in (("fp101_d40", d40), ("fp101_d100", d100)):
        for i, f in enumerate(lst):
            lineas.append(f"{nombre} {i} " + " ".join(str(int(c)) for c in f.coeffs()))
        res.append((nombre + " factor", medir(lambda i, l=lst: l[i].factor(), N)))

    # ---- extensiones (fq_default_poly) ----
    for nombre, pp, n, mod, casos in (
        ("f3e30", 3, 30, MOD_3_30, (("d5", "factor", 5), ("mul_d200", "mul", 200))),
        ("f2e8", 2, 8, MOD_2_8, (("d30", "factor", 30), ("mul_d1000", "mul", 1000))),
    ):
        q = pp ** n
        ctx = flint.fq_default_ctx(pp, n, "a", modulus=flint.fmpz_mod_poly_ctx(pp)(mod))
        R = flint.fq_default_poly_ctx(ctx)
        conv = lambda c: R([ctx(digitos(i, pp, n)) for i in c])
        for sufijo, op, grado in casos:
            caso = f"{nombre}_{sufijo}"
            if op == "factor":
                ent = [aleatorio(rng, q, grado, True) for _ in range(N)]
                for i, c in enumerate(ent):
                    lineas.append(f"{caso} {i} " + " ".join(map(str, c)))
                P = [conv(c) for c in ent]
                res.append((caso + " factor", medir(lambda i, P=P: P[i].factor(), N)))
            else:
                A, B = [], []
                for i in range(N):
                    a, b = aleatorio(rng, q, grado), aleatorio(rng, q, grado)
                    lineas.append(f"{caso} {i} " + " ".join(map(str, a)))
                    lineas.append(f"{caso} {i} " + " ".join(map(str, b)))
                    A.append(conv(a))
                    B.append(conv(b))
                res.append((caso + " mul", medir(lambda i, A=A, B=B: A[i] * B[i], N)))
        print(f"  ({nombre}: tipo FLINT {ctx.type if hasattr(ctx, 'type') else '?'})")

    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    SALIDA.write_text("\n".join(lineas) + "\n")
    print(f"python-flint {flint.__version__}: mejor de 3 tandas de >= 30 ms, rotando {N} entradas")
    for nombre, t in res:
        print(f"  {nombre:24s} {t:10.1f} us")
    print(f"Entradas guardadas en {SALIDA.relative_to(REPO)}")


if __name__ == "__main__":
    main()
