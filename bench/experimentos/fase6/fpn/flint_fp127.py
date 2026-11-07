# flint_fp127.py — Fase 6, tema fpn.  Genera la entrada exacta del escenario
# fp127.factor.d20 de bench/compare_libs.py (misma semilla: sha256 del nombre)
# en build/fase6/fpn/fp127_factor_d20.txt y mide factor() con python-flint
# (referencia; mejor de 3 tandas de >= 30 ms).
#
#   .venv/Scripts/python.exe bench/experimentos/fase6/fpn/flint_fp127.py
import hashlib
import random
import time

import flint

q = 2**127 - 1
semilla = int(hashlib.sha256("fp127.factor.d20".encode()).hexdigest()[:8], 16)
rng = random.Random(semilla)
c = [rng.randrange(q) for _ in range(20)] + [1]  # mónico, como _aleatorio(..., monico=True)
with open("build/fase6/fpn/fp127_factor_d20.txt", "w") as fh:
    fh.write(" ".join(map(str, c)) + "\n")

ctx = flint.fmpz_mod_poly_ctx(q)
f = ctx(c)
grados = sorted(g.degree() for g, _ in f.factor()[1])
mejor = float("inf")
for _ in range(3):
    n, t0 = 0, time.perf_counter()
    while True:
        f.factor()
        n += 1
        dt = time.perf_counter() - t0
        if dt >= 0.03:
            break
    mejor = min(mejor, dt / n)
print(f"python-flint {flint.__version__}: factor() fp127.factor.d20 (grados {grados}): {mejor * 1e3:.2f} ms")
