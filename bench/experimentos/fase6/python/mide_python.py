"""mide_python.py — Fase 6, tema «python»: coste fijo de las llamadas desde Python.

Compara, en F_998244353 y con las MISMAS entradas que bench_cpp.cpp (splitmix64,
16 pares rotando):
  * alcp_native (el módulo real, como lo usa bench/compare_libs.py:
    ``K = nat.Fp(p); f = K.poly(c); f * g``),
  * el prototipo alcp_proto (build/fase6/python/, ver proto_modulo.cpp):
    PolyP (pybind11, objeto opaco sin variant) y RawPoly (tipo CPython a mano),
  * python-flint (nmod_poly / nmod).
Desglosa: producto, suma, construcción desde lista (conversión de entrada),
coeffs() (conversión de salida), degree() (despacho de un método), producto que
no crea objeto (mul_grado) y operaciones de elementos.

Comprueba que el prototipo da EXACTAMENTE los mismos resultados que alcp_native
(productos, sumas, coeffs, elementos con enteros pequeños, grandes, negativos,
bool y float rechazado) y que el hash de los productos coincide con el de
bench_cpp.exe.

Metodología: mejor de 3 tandas de >= 30 ms, 16 entradas rotando, gc desactivado;
se resta el coste del bucle vacío solo en la columna «neto».

Ejecución (desde la raíz del repositorio):
  ALCP_THREADS=1 .venv/Scripts/python.exe bench/experimentos/fase6/python/mide_python.py
"""
from __future__ import annotations

import gc
import sys
import time
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(RAIZ / "python-bindings"))
sys.path.insert(0, str(RAIZ / "build" / "fase6" / "python"))

from alcp import alcp_native as nat  # noqa: E402  (como compare_libs.py)
import alcp_proto as pro  # noqa: E402
import flint  # noqa: E402

P = 998244353
M64 = (1 << 64) - 1


class SplitMix:
    """splitmix64, idéntico al de bench_cpp.cpp."""

    def __init__(self, s: int):
        self.s = s

    def __call__(self) -> int:
        self.s = (self.s + 0x9E3779B97F4A7C15) & M64
        z = self.s
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
        return z ^ (z >> 31)


def coefs(r: SplitMix, d: int) -> list[int]:
    c = [r() % P for _ in range(d)]
    c.append(r() % (P - 1) + 1)
    return c


def entradas(d: int):
    r = SplitMix(d)
    F, G = [], []
    for _ in range(16):
        F.append(coefs(r, d))
        G.append(coefs(r, d))
    return F, G


def medir(fn, args) -> float:
    """ns por llamada de fn(*args[i % 16]); mejor de 3 tandas de >= 30 ms."""
    mejor = float("inf")
    n_args = len(args)
    for _ in range(3):
        n = 0
        t0 = time.perf_counter()
        while True:
            for i in range(256):
                fn(*args[i % n_args])
            n += 256
            dt = time.perf_counter() - t0
            if dt >= 0.030:
                break
        mejor = min(mejor, dt / n)
    return mejor * 1e9


def vacio(*_a):
    return None


def main() -> None:
    gc.disable()
    Kn = nat.Fp(P)
    Kp = pro.Cuerpo(P)
    pro.fijar_p(P)
    lineas: list[str] = []
    out = lambda s="": (print(s), lineas.append(s))  # noqa: E731
    out(f"Python {sys.version.split()[0]}, python-flint {flint.__version__}, "
        f"F_{P}; ns por llamada (mejor de 3 tandas >= 30 ms, 16 entradas rotando)")
    base = medir(vacio, [(1, 2)] * 16)
    out(f"bucle vacío (llamada a una función Python que no hace nada): {base:.0f} ns "
        "(no se resta salvo en la columna «neto»)")
    out()

    # ------------------------------------------------------------------ polinomios
    for d in (10, 100, 1000):
        F, G = entradas(d)
        # Objetos de cada librería.
        nf = [Kn.poly(c) for c in F]; ng = [Kn.poly(c) for c in G]
        pf = [Kp.poly(c) for c in F]; pg = [Kp.poly(c) for c in G]
        rf = [pro.raw_poly(Kp, c) for c in F]; rg = [pro.raw_poly(Kp, c) for c in G]
        hf = [Kp.poly_h(c) for c in F]; hg = [Kp.poly_h(c) for c in G]
        ff = [flint.nmod_poly(c, P) for c in F]; fg = [flint.nmod_poly(c, P) for c in G]

        # ---- comprobación de resultados idénticos (frente a la librería) ----
        h = 0
        for i in range(16):
            ref = (nf[i] * ng[i]).coeffs()
            assert (pf[i] * pg[i]).coeffs() == ref, "PolyP.__mul__"
            assert (rf[i] * rg[i]).coeffs() == ref, "RawPoly.__mul__"
            assert (hf[i] * hg[i]).coeffs() == ref, "PolyH.__mul__"
            assert [int(c) for c in (ff[i] * fg[i]).coeffs()] == ref, "flint"
            assert (pf[i] + pg[i]).coeffs() == (nf[i] + ng[i]).coeffs() == (rf[i] + rg[i]).coeffs()
            assert pf[i].mul_grado(pg[i]) == len(ref) - 1
            assert Kp.poly(F[i]).coeffs() == Kn.poly(F[i]).coeffs() == Kp.poly_stl(F[i]).coeffs()
            h += sum(c * (j + 1) for j, c in enumerate(ref))
        out(f"d={d}: resultados idénticos (16 productos y sumas); hash productos = {h & M64}")

        filas = [
            ("mul (f*g)", [
                ("alcp_native", lambda a, b: a * b, list(zip(nf, ng))),
                ("proto PolyP", lambda a, b: a * b, list(zip(pf, pg))),
                ("proto PolyH (ranura nb_multiply)", lambda a, b: a * b, list(zip(hf, hg))),
                ("proto RawPoly", lambda a, b: a * b, list(zip(rf, rg))),
                ("flint", lambda a, b: a * b, list(zip(ff, fg))),
            ]),
            ("mul sin crear objeto", [
                ("proto PolyP.mul_grado", lambda a, b: a.mul_grado(b), list(zip(pf, pg))),
            ]),
            ("add (f+g)", [
                ("alcp_native", lambda a, b: a + b, list(zip(nf, ng))),
                ("proto PolyP", lambda a, b: a + b, list(zip(pf, pg))),
                ("proto RawPoly", lambda a, b: a + b, list(zip(rf, rg))),
                ("flint", lambda a, b: a + b, list(zip(ff, fg))),
            ]),
            ("lista -> polinomio", [
                ("alcp_native K.poly", lambda c: Kn.poly(c), [(c,) for c in F]),
                ("proto K.poly (directa)", lambda c: Kp.poly(c), [(c,) for c in F]),
                ("proto K.poly_stl (vector<u64>)", lambda c: Kp.poly_stl(c), [(c,) for c in F]),
                ("proto raw_poly", lambda c: pro.raw_poly(Kp, c), [(c,) for c in F]),
                ("flint nmod_poly(c, p)", lambda c: flint.nmod_poly(c, P), [(c,) for c in F]),
            ]),
            ("polinomio -> lista (coeffs)", [
                ("alcp_native", lambda a: a.coeffs(), [(x,) for x in nf]),
                ("proto PolyP", lambda a: a.coeffs(), [(x,) for x in pf]),
                ("proto RawPoly", lambda a: a.coeffs(), [(x,) for x in rf]),
                ("flint (coeffs, nmod)", lambda a: a.coeffs(), [(x,) for x in ff]),
                ("flint [int(c) ...]", lambda a: [int(c) for c in a.coeffs()], [(x,) for x in ff]),
            ]),
            ("degree() (despacho de un método)", [
                ("alcp_native", lambda a: a.degree(), [(x,) for x in nf]),
                ("proto PolyP", lambda a: a.degree(), [(x,) for x in pf]),
                ("proto RawPoly", lambda a: a.degree(), [(x,) for x in rf]),
                ("flint", lambda a: a.degree(), [(x,) for x in ff]),
            ]),
        ]
        for titulo, variantes in filas:
            out(f"  {titulo}")
            for nombre, fn, args in variantes:
                t = medir(fn, args)
                out(f"    {nombre:34s} {t:10.0f} ns   neto {t - base:10.0f} ns")
        out()

    # ------------------------------------------------------------------ elementos
    r = SplitMix(7)
    A = [r() % P for _ in range(16)]
    B = [r() % P for _ in range(16)]
    pares = [(A[i], B[(i + 3) % 16]) for i in range(16)]
    # Comprobación: pequeños, grandes, negativos, bool, frontera de p.
    raros = [0, 1, -1, P - 1, P, P + 1, -P, 2**30 - 1, 2**30, -(2**30), 2**62, 2**63 - 1, 2**63,
             -(2**63), 2**64, 2**100 + 12345, -(2**100) - 7, True, False]
    for a in A + raros:
        for b in B[:4] + raros:
            assert Kp.mul(a, b) == Kn.mul(a, b) == pro.rmul(a, b), (a, b)
            assert Kp.add(a, b) == Kn.add(a, b) == pro.radd(a, b), (a, b)
    for malo in (2.0, "3", None):
        for f in (Kn.mul, Kp.mul, pro.rmul):
            try:
                f(malo, 1)
            except TypeError:
                pass
            else:
                raise AssertionError(f"{f} aceptó {malo!r}")
    out("elementos: resultados idénticos (16 pares + enteros extremos, negativos, bool; float/str rechazados)")
    fa = [flint.nmod(a, P) for a in A]
    fb = [flint.nmod(b, P) for b in B]
    fpares = [(fa[i], fb[(i + 3) % 16]) for i in range(16)]
    elem = [
        ("mul de elementos", [
            ("alcp_native F.mul(a, b)", lambda a, b: Kn.mul(a, b), pares),
            ("proto Cuerpo.mul (directa)", lambda a, b: Kp.mul(a, b), pares),
            ("proto Cuerpo.mul_u64 (caster u64)", lambda a, b: Kp.mul_u64(a, b), pares),
            ("proto rmul (CPython FASTCALL)", lambda a, b: pro.rmul(a, b), pares),
            ("flint nmod * nmod", lambda a, b: a * b, fpares),
            ("Python a * b % p", lambda a, b: a * b % P, pares),
        ]),
        ("add de elementos", [
            ("alcp_native F.add(a, b)", lambda a, b: Kn.add(a, b), pares),
            ("proto Cuerpo.add (directa)", lambda a, b: Kp.add(a, b), pares),
            ("proto radd (CPython FASTCALL)", lambda a, b: pro.radd(a, b), pares),
            ("flint nmod + nmod", lambda a, b: a + b, fpares),
        ]),
        ("llamada vacía", [
            ("proto noop (pybind11)", lambda: pro.noop(), [()]),
            ("proto rnoop (CPython)", lambda: pro.rnoop(), [()]),
        ]),
    ]
    for titulo, variantes in elem:
        out(f"  {titulo}")
        for nombre, fn, args in variantes:
            t = medir(fn, args)
            out(f"    {nombre:34s} {t:10.0f} ns   neto {t - base:10.0f} ns")
    # Todas las variantes van envueltas en una lambda, como el bucle vacío: «neto»
    # es el coste propio de la operación (sin la llamada a la lambda ni el bucle).


if __name__ == "__main__":
    main()
