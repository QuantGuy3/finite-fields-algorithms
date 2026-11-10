# Fase 6, tema ntt: referencia de python-flint (nmod_poly.mul) y del módulo
# alcp desde Python para productos de grado 1000 y 10000 en F_998244353.
# Ejecución (desde la raíz del repositorio):
#   ALCP_THREADS=1 .venv/Scripts/python.exe bench/experimentos/fase6/ntt/flint_ref.py
# Metodología: 16 entradas distintas rotando, mejor de 3 tandas de >= 30 ms.
import random
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "python-bindings"))

import flint

P = 998244353
N_ENTRADAS = 16


def medir(f, args):
    """Mejor de 3 tandas de >= 30 ms; f(*args[i]) rotando i; segundos por llamada."""
    mejor = float("inf")
    for _ in range(3):
        reps = 0
        t0 = time.perf_counter()
        while True:
            for a in args:
                f(*a)
            reps += len(args)
            dt = time.perf_counter() - t0
            if dt >= 0.03:
                break
        mejor = min(mejor, dt / reps)
    return mejor


def main():
    rnd = random.Random(20261107)
    print("# referencia python-flint", flint.__version__, "y alcp desde Python, F_998244353")
    try:
        import alcp
        K = alcp.GF(P)
        R = K["x"]
    except Exception as e:  # el módulo puede no estar disponible
        alcp = None
        print("# alcp no disponible:", e)
    for d in (1000, 10000):
        coefs = [([rnd.randrange(P) for _ in range(d + 1)], [rnd.randrange(P) for _ in range(d + 1)])
                 for _ in range(N_ENTRADAS)]
        fa = [(flint.nmod_poly(a, P), flint.nmod_poly(b, P)) for a, b in coefs]
        tf = medir(lambda x, y: x * y, fa)
        linea = f"grado {d:5d}: flint nmod_poly mul {tf * 1e6:9.1f} us"
        if alcp is not None:
            aa = [(R([K(c) for c in a]), R([K(c) for c in b])) for a, b in coefs]
            # Comprobación: mismo producto que FLINT.
            x, y = aa[0]
            fx, fy = fa[0]
            assert [int(c) for c in (fx * fy).coeffs()] == [K.to_index(c) if hasattr(K, "to_index") else int(c)
                                                            for c in (x * y).coeffs], "producto distinto"
            ta = medir(lambda x, y: x * y, aa)
            linea += f"   alcp (Python) {ta * 1e6:9.1f} us   alcp/flint {ta / tf:5.2f}x"
        print(linea)


if __name__ == "__main__":
    main()
