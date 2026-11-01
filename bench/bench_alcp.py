#!/usr/bin/env python3
"""bench_alcp.py — Sobrecoste de la capa Python (paquete ``alcp``) frente al módulo nativo.

Cada escenario mide el MISMO cálculo dos veces:

* ``alcp.<nombre>``: con la API del paquete (``GF``, elementos, ``K["x"]``...);
* ``nat.<nombre>``: con el módulo nativo ``alcp.alcp_native`` directamente
  (enteros y ``PolyFp``/``PolyFq``), cuando hay un equivalente.

El cociente ``alcp/nat`` cuantifica el sobrecoste de la capa.  Si los dos
cálculos son equivalentes, sus checksums coinciden (se comprueba y se avisa si
no).  El formato del JSON es el de ``bench_python.py``/``run_bench.py``.

Uso::

    python bench/bench_alcp.py --pkg-dir python-bindings --json salida.json
    python bench/bench_alcp.py --pkg-dir DIR --filter "fp\\." --reps 5 --json s.json
    python bench/bench_alcp.py --list
    # A/B intercalado (cada ronda ejecuta A y B en procesos nuevos, alternando):
    python bench/bench_alcp.py --interleave base=DIR_A nuevo=DIR_B --rounds 5 --out-dir DIR

``--pkg-dir`` es el directorio que CONTIENE el paquete ``alcp/`` (con su módulo
nativo compilado dentro).  Para comparar dos versiones de la capa Python con el
mismo binario basta con copiar el paquete (con el ``.pyd``) a otro directorio.
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import pickle
import random
import re
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Callable, Optional

BENCH_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(BENCH_DIR))
# Utilidades compartidas: checksum, semillas, metadatos y salida UTF-8.
from bench_python import _digest, _seed  # noqa: E402
from run_bench import collect_machine_metadata, scrub_paths, utf8_stdout  # noqa: E402

P_NTT = 998244353
P101 = 101


class Escenario:
    """Un cálculo medido con el paquete y (opcionalmente) con el nativo.

    ``setup(alcp, nativo, rng)`` devuelve ``(fn_paquete, fn_nativo | None,
    normalizar_paquete, normalizar_nativo)``; los normalizadores convierten el
    resultado en enteros/listas para el checksum.
    """

    def __init__(self, nombre: str, desc: str, inner: int, setup: Callable):
        self.nombre, self.desc, self.inner, self.setup = nombre, desc, inner, setup


def _rep(fn: Callable[[], Any], n: int) -> Callable[[], Any]:
    if n == 1:
        return fn

    def run():
        r = None
        for _ in range(n):
            r = fn()
        return r
    return run


def _coefs(rng: random.Random, cota: int, grado: int) -> list:
    c = [rng.randrange(cota) for _ in range(grado)]
    c.append(rng.randrange(1, cota))
    return c


def _ind_poli(f) -> list:
    """Índices de los coeficientes de un polinomio del paquete."""
    return [c.index for c in f.coeffs]


def _ind_fac(fac) -> tuple:
    return (sorted((tuple(_ind_poli(u)), e) for u, e in fac), fac.unit.index)


def _ind_fac_nat(fac) -> tuple:
    factores, unidad = fac
    return (sorted((tuple(u.coeffs()), int(e)) for u, e in factores), unidad)


def construir() -> list:
    E: list = []

    def add(nombre: str, desc: str, inner: int = 1):
        def deco(setup):
            E.append(Escenario(nombre, desc, inner, setup))
            return setup
        return deco

    ident = (lambda x: x)

    # ---- aritmética de elementos en bucles Python -------------------------
    def bucle_binario(nombre, desc, n, q_pkg, q_nat, op_pkg, op_nat):
        @add(nombre, desc, n)
        def _(alcp, nat, rng):
            K = q_pkg(alcp)
            Kn = q_nat(nat, K)
            q = K.order
            a0, x0 = rng.randrange(1, q), rng.randrange(1, q)
            a, x = K.from_index(a0), K.from_index(x0)

            def pkg():
                y = x
                for _ in range(n):
                    y = op_pkg(y, a)
                return y

            def nativo():
                y = x0
                for _ in range(n):
                    y = op_nat(Kn, y, a0)
                return y
            return pkg, nativo, (lambda e: e.index), ident

    fp = (lambda alcp: alcp.GF(P_NTT)), (lambda nat, K: nat.Fp(P_NTT))
    f9 = (lambda alcp: alcp.GF(9)), (lambda nat, K: nat.Fq(3, [1, 0, 1]))
    f256 = ((lambda alcp: alcp.GF(256)),
            (lambda nat, K: nat.Fq(2, [c.index for c in K.modulus.coeffs])))

    bucle_binario("fp.add_loop.n20000", "y = y + a, 20000 veces en F_998244353", 20000, *fp,
                  lambda y, a: y + a, lambda K, y, a: K.add(y, a))
    bucle_binario("fp.mul_loop.n20000", "y = y * a, 20000 veces en F_998244353", 20000, *fp,
                  lambda y, a: y * a, lambda K, y, a: K.mul(y, a))
    bucle_binario("fp.sub_loop.n20000", "y = y - a, 20000 veces en F_998244353", 20000, *fp,
                  lambda y, a: y - a, lambda K, y, a: K.sub(y, a))
    bucle_binario("fp.div_loop.n5000", "y = y / a, 5000 veces en F_998244353", 5000, *fp,
                  lambda y, a: y / a, lambda K, y, a: K.div(y, a))
    bucle_binario("fp.mixed_int.n10000", "y = 3*y + 1 (enteros mezclados), 10000 veces en F_998244353",
                  10000, *fp, lambda y, a: 3 * y + 1, lambda K, y, a: K.add(K.mul(3, y), 1))
    bucle_binario("fp.pow_loop.n5000", "y = y**3, 5000 veces en F_998244353", 5000, *fp,
                  lambda y, a: y ** 3, lambda K, y, a: K.pow(y, 3))
    bucle_binario("fp.eq_loop.n20000", "y == a (False), 20000 veces en F_998244353", 20000, *fp,
                  lambda y, a: (y == a) or y, lambda K, y, a: (y == a) or y)
    bucle_binario("f9.mul_add.n10000", "y = y*a + a, 10000 veces en F_9", 10000, *f9,
                  lambda y, a: y * a + a, lambda K, y, a: K.add(K.mul(y, a), a))
    bucle_binario("gf256.mul_loop.n20000", "y = y * a, 20000 veces en GF(2^8)", 20000, *f256,
                  lambda y, a: y * a, lambda K, y, a: K.mul(y, a))
    bucle_binario("gf256.add_int.n10000", "y = y + 1 (entero), 10000 veces en GF(2^8)", 10000, *f256,
                  lambda y, a: y + 1, lambda K, y, a: K.add(y, 1))

    @add("fp.convert_int.n20000", "K(i) para 20000 enteros en F_998244353", 20000)
    def _(alcp, nat, rng):
        K, Kn = alcp.GF(P_NTT), nat.Fp(P_NTT)
        xs = [rng.randrange(-10**12, 10**12) for _ in range(20000)]
        return ((lambda: [K(i) for i in xs]), (lambda: [Kn.from_int(i) for i in xs]),
                lambda r: [e.index for e in r], ident)

    @add("fp.sum_elems.n20000", "sum() de 20000 elementos de F_998244353", 20000)
    def _(alcp, nat, rng):
        K, Kn = alcp.GF(P_NTT), nat.Fp(P_NTT)
        xs = [rng.randrange(P_NTT) for _ in range(20000)]
        es = [K(i) for i in xs]

        def nativo():
            s = 0
            for x in xs:
                s = Kn.add(s, x)
            return s
        return (lambda: sum(es)), nativo, (lambda e: e.index), ident

    # ---- construcción de polinomios ----------------------------------------
    @add("poly.from_list.d200", "R(lista) de grado 200 en F_101 (x20)", 20)
    def _(alcp, nat, rng):
        R, Fn = alcp.GF(P101)["x"], nat.Fp(P101)
        lst = _coefs(rng, P101, 200)
        return _rep(lambda: R(lst), 20), _rep(lambda: Fn.poly(lst), 20), _ind_poli, lambda P: P.coeffs()

    @add("poly.from_elems.d200", "R(lista de elementos) de grado 200 en F_101 (x20)", 20)
    def _(alcp, nat, rng):
        K = alcp.GF(P101)
        R, Fn = K["x"], nat.Fp(P101)
        lst = _coefs(rng, P101, 200)
        es = [K(c) for c in lst]
        return _rep(lambda: R(es), 20), _rep(lambda: Fn.poly(lst), 20), _ind_poli, lambda P: P.coeffs()

    @add("poly.from_dict.d1000", "R({exp: coef}) con 40 términos y grado 1000 en F_101 (x20)", 20)
    def _(alcp, nat, rng):
        R, Fn = alcp.GF(P101)["x"], nat.Fp(P101)
        d = {rng.randrange(1000): rng.randrange(1, P101) for _ in range(40)}
        d[1000] = 1
        dense = [0] * 1001
        for k, v in d.items():
            dense[k] = v
        return _rep(lambda: R(d), 20), _rep(lambda: Fn.poly(dense), 20), _ind_poli, lambda P: P.coeffs()

    @add("poly.from_str.d40", "R(cadena unicode) de grado 40 en F_101 (x20)", 20)
    def _(alcp, nat, rng):
        R, Fn = alcp.GF(P101)["x"], nat.Fp(P101)
        lst = _coefs(rng, P101, 40)
        texto = format(R(lst), "unicode")
        return _rep(lambda: R(texto), 20), _rep(lambda: Fn.poly(lst), 20), _ind_poli, lambda P: P.coeffs()

    @add("poly.from_str_f9.d40", "R(cadena ascii) de grado 40 en F_9 (x10)", 10)
    def _(alcp, nat, rng):
        R, Kn = alcp.GF(9)["x"], nat.Fq(3, [1, 0, 1])
        lst = _coefs(rng, 9, 40)
        texto = format(R.field["x"]([R.field.from_index(i) for i in lst]), "ascii")
        return _rep(lambda: R(texto), 10), _rep(lambda: Kn.poly(lst), 10), _ind_poli, lambda P: P.coeffs()

    @add("poly.parse_long.d500", "R(cadena ascii larga) de grado 500 en F_101", 1)
    def _(alcp, nat, rng):
        R, Fn = alcp.GF(P101)["x"], nat.Fp(P101)
        lst = _coefs(rng, P101, 500)
        texto = format(R(lst), "ascii")
        return (lambda: R(texto)), (lambda: Fn.poly(lst)), _ind_poli, lambda P: P.coeffs()

    @add("poly.from_str_repeat.n200", "R(cadena corta) repetida 200 veces en F_101", 200)
    def _(alcp, nat, rng):
        R, Fn = alcp.GF(P101)["x"], nat.Fp(P101)
        texto = "x^3 + 2*x + 1"
        return ((lambda: [R(texto) for _ in range(200)][-1]),
                (lambda: [Fn.poly([1, 2, 0, 1]) for _ in range(200)][-1]),
                _ind_poli, lambda P: P.coeffs())

    # ---- operaciones entre polinomios --------------------------------------
    def bin_poli(nombre, desc, n, op, grado=100, divisor=50):
        @add(nombre, desc, n)
        def _(alcp, nat, rng):
            R, Fn = alcp.GF(P101)["x"], nat.Fp(P101)
            la, lb = _coefs(rng, P101, grado), _coefs(rng, P101, divisor)
            a, b = R(la), R(lb)
            an, bn = Fn.poly(la), Fn.poly(lb)
            return _rep(lambda: op(a, b), n), _rep(lambda: op(an, bn), n), _ind_poli, lambda P: P.coeffs()

    bin_poli("poly.add.d100", "f + g con grados 100 y 50 en F_101 (x200)", 200, lambda a, b: a + b)
    bin_poli("poly.sub.d100", "f - g con grados 100 y 50 en F_101 (x200)", 200, lambda a, b: a - b)
    bin_poli("poly.mul.d100", "f * g con grados 100 y 50 en F_101 (x100)", 100, lambda a, b: a * b)
    bin_poli("poly.floordiv.d100", "f // g con grados 100 y 50 en F_101 (x100)", 100, lambda a, b: a // b)
    bin_poli("poly.mod.d100", "f % g con grados 100 y 50 en F_101 (x100)", 100, lambda a, b: a % b)
    bin_poli("poly.small_mul.d5", "f * g con grados 5 y 3 en F_101 (x500)", 500, lambda a, b: a * b, 5, 3)

    @add("poly.add_int.d100", "f + 1 con grado 100 en F_101 (x200)", 200)
    def _(alcp, nat, rng):
        R, Fn = alcp.GF(P101)["x"], nat.Fp(P101)
        la = _coefs(rng, P101, 100)
        a, an = R(la), Fn.poly(la)
        return (_rep(lambda: a + 1, 200), _rep(lambda: an + Fn.poly([1]), 200),
                _ind_poli, lambda P: P.coeffs())

    @add("poly.rmul_int.d100", "3 * f con grado 100 en F_101 (x200)", 200)
    def _(alcp, nat, rng):
        R, Fn = alcp.GF(P101)["x"], nat.Fp(P101)
        la = _coefs(rng, P101, 100)
        a, an = R(la), Fn.poly(la)
        # Equivalente nativo de la ruta del paquete: producto por el polinomio constante.
        tres = Fn.poly([3])
        return (_rep(lambda: 3 * a, 200), _rep(lambda: tres * an, 200),
                _ind_poli, lambda P: P.coeffs())

    @add("poly.mul_elem_f9.d100", "f * c (c elemento) con grado 100 en F_9 (x200)", 200)
    def _(alcp, nat, rng):
        K = alcp.GF(9)
        R, Kn = K["x"], nat.Fq(3, [1, 0, 1])
        la = _coefs(rng, 9, 100)
        # Los enteros en una extensión son n·1: se construye desde los índices.
        a = R([K.from_index(i) for i in la])
        an, c, cn = Kn.poly(la), K.from_index(5), Kn.poly([5])
        return (_rep(lambda: a * c, 200), _rep(lambda: an * cn, 200),
                _ind_poli, lambda P: P.coeffs())

    @add("poly.eq.d100", "f == g (iguales) con grado 100 en F_101 (x500)", 500)
    def _(alcp, nat, rng):
        R, Fn = alcp.GF(P101)["x"], nat.Fp(P101)
        la = _coefs(rng, P101, 100)
        a, b, an, bn = R(la), R(list(la)), Fn.poly(la), Fn.poly(list(la))
        return _rep(lambda: a == b, 500), _rep(lambda: an == bn, 500), ident, ident

    # ---- evaluación --------------------------------------------------------
    @add("poly.eval.d100.n500", "f(v) con grado 100 en F_101 para 500 enteros", 500)
    def _(alcp, nat, rng):
        R, Fn = alcp.GF(P101)["x"], nat.Fp(P101)
        la = _coefs(rng, P101, 100)
        a, an = R(la), Fn.poly(la)
        pts = [rng.randrange(P101) for _ in range(500)]
        return ((lambda: [a(v) for v in pts]), (lambda: [an.evaluate(v) for v in pts]),
                lambda r: [e.index for e in r], ident)

    @add("poly.eval_elem.d100.n500", "f(e) con grado 100 en F_101 para 500 elementos", 500)
    def _(alcp, nat, rng):
        K = alcp.GF(P101)
        R, Fn = K["x"], nat.Fp(P101)
        la = _coefs(rng, P101, 100)
        a, an = R(la), Fn.poly(la)
        pts = [rng.randrange(P101) for _ in range(500)]
        es = [K(v) for v in pts]
        return ((lambda: [a(e) for e in es]), (lambda: [an.evaluate(v) for v in pts]),
                lambda r: [e.index for e in r], ident)

    # ---- factorización y raíces --------------------------------------------
    @add("poly.factor.d40", "factor() de un polinomio de grado 40 en F_101", 1)
    def _(alcp, nat, rng):
        R, Fn = alcp.GF(P101)["x"], nat.Fp(P101)
        f = Fn.poly([rng.randrange(1, P101)])
        for d, mult in ((1, 2), (2, 1), (3, 2), (5, 1), (8, 1), (15, 1)):
            u = Fn.poly([rng.randrange(P101) for _ in range(d)] + [1])
            for _ in range(mult):
                f = f * u
        fp = R(f.coeffs())
        return fp.factor, f.factor, _ind_fac, _ind_fac_nat

    @add("poly.roots.d40", "roots() de un polinomio de grado 40 con 20 raíces en F_101", 1)
    def _(alcp, nat, rng):
        R, Fn = alcp.GF(P101)["x"], nat.Fp(P101)
        f = Fn.poly(_coefs(rng, P101, 20))
        for r in rng.sample(range(P101), 20):
            f = f * Fn.poly([(-r) % P101, 1])
        fp = R(f.coeffs())

        def nativo():
            x = Fn.poly([0, 1])
            g = f.gcd(x.pow_mod(P101, f) - x)
            lineales = g.equal_degree(1) if g.degree() > 1 else [g]
            raices = sorted(Fn.neg(u.coeff(0)) for u in lineales)
            return [(r, f.multiplicity(Fn.poly([Fn.neg(r), 1]))) for r in raices]
        return fp.roots, nativo, lambda r: [(e.index, m) for e, m in r], ident

    # ---- representación textual ---------------------------------------------
    def texto(nombre, desc, n, q, grado, espec, nat_q):
        @add(nombre, desc, n)
        def _(alcp, nat, rng):
            K = q(alcp)
            R = K["x"]
            la = _coefs(rng, K.order, grado)
            f = R([K.from_index(i) for i in la])
            fn = nat_q(nat).poly(la)
            return (_rep(lambda: format(f, espec), n), _rep(lambda: fn.to_string("x"), n),
                    ident, None)

    texto("poly.str_unicode.d500", "format(f, 'unicode') con grado 500 en F_101 (x5)", 5,
          lambda alcp: alcp.GF(P101), 500, "unicode", lambda nat: nat.Fp(P101))
    texto("poly.str_ascii.d500", "format(f, 'ascii') con grado 500 en F_101 (x5)", 5,
          lambda alcp: alcp.GF(P101), 500, "ascii", lambda nat: nat.Fp(P101))
    texto("poly.str_latex.d500", "format(f, 'latex') con grado 500 en F_101 (x5)", 5,
          lambda alcp: alcp.GF(P101), 500, "latex", lambda nat: nat.Fp(P101))
    texto("poly.str_f9.d200", "format(f, 'unicode') con grado 200 en F_9 (x5)", 5,
          lambda alcp: alcp.GF(9), 200, "unicode", lambda nat: nat.Fq(3, [1, 0, 1]))

    @add("elem.str_gf256.n2000", "str() de 2000 elementos de GF(2^8)", 2000)
    def _(alcp, nat, rng):
        K = alcp.GF(256)
        es = [K.from_index(rng.randrange(256)) for _ in range(2000)]
        return (lambda: [str(e) for e in es]), None, ident, None

    # ---- pickle -------------------------------------------------------------
    @add("pickle.poly.d500", "pickle.loads(pickle.dumps(f)) con grado 500 en F_101 (x20)", 20)
    def _(alcp, nat, rng):
        R, Fn = alcp.GF(P101)["x"], nat.Fp(P101)
        la = _coefs(rng, P101, 500)
        f, fn = R(la), Fn.poly(la)
        return (_rep(lambda: pickle.loads(pickle.dumps(f)), 20),
                _rep(lambda: Fn.poly(pickle.loads(pickle.dumps(fn.coeffs()))), 20),
                _ind_poli, lambda P: P.coeffs())

    @add("pickle.elems.n1000", "pickle de una lista de 1000 elementos de F_101 (x5)", 5)
    def _(alcp, nat, rng):
        K = alcp.GF(P101)
        xs = [rng.randrange(P101) for _ in range(1000)]
        es = [K(x) for x in xs]
        return (_rep(lambda: pickle.loads(pickle.dumps(es)), 5),
                _rep(lambda: pickle.loads(pickle.dumps(xs)), 5),
                lambda r: [e.index for e in r], ident)

    return E


# ---------------------------------------------------------------------------
# Ejecución
# ---------------------------------------------------------------------------
def medir(fn: Callable[[], Any], reps: int, min_time: float, warmup: int, max_reps: int):
    muestras, calent, out = [], [], None
    gc_was = gc.isenabled()
    try:
        for i in range(warmup + max_reps):
            gc.collect()
            gc.disable()
            t0 = time.perf_counter_ns()
            out = fn()
            dt = (time.perf_counter_ns() - t0) * 1e-9
            if gc_was:
                gc.enable()
            if i < warmup:
                calent.append(dt)
                continue
            muestras.append(dt)
            if len(muestras) >= reps and sum(muestras) >= min_time:
                break
    finally:
        if gc_was:
            gc.enable()
    return muestras, calent, out


def _resultado(nombre: str, desc: str, inner: int) -> dict:
    grupo = nombre.split(".")[0] + "." + nombre.split(".")[1]
    return {"name": nombre, "group": grupo, "op": nombre.split(".")[2], "desc": desc,
            "inner": inner, "samples_s": [], "warmup_s": [], "checksum": "", "error": ""}


def ejecutar(alcp, nat, esc: Escenario, a) -> list:
    rp, rn = (_resultado("alcp." + esc.nombre, esc.desc, esc.inner),
              _resultado("nat." + esc.nombre, esc.desc + " [módulo nativo]", esc.inner))
    try:
        fp, fn, normp, normn = esc.setup(alcp, nat, random.Random(_seed(esc.nombre)))
    except Exception as e:  # se registra y se sigue
        rp["error"] = f"{type(e).__name__}: {e}"
        return [rp]
    res = [rp]
    for r, f, norm in ((rp, fp, normp), (rn, fn, normn)):
        if f is None:
            continue
        if r is rn:
            res.append(rn)
        try:
            r["samples_s"], r["warmup_s"], out = medir(f, a.reps, a.min_time, a.warmup, a.max_reps)
            r["checksum"] = _digest(norm(out) if norm else out)
        except Exception as e:
            r["error"] = f"{type(e).__name__}: {e}"
        if r["samples_s"]:
            r["median_s"] = statistics.median(r["samples_s"])
            r["min_s"] = min(r["samples_s"])
    # Equivalencia: si hay normalizador nativo, los checksums deben coincidir.
    rp["equivalente"] = bool(fn is not None and normn is not None)
    if rp["equivalente"] and not rp["error"] and not rn["error"] and rp["checksum"] != rn["checksum"]:
        rp["error"] = f"checksum distinto del nativo ({rp['checksum']} != {rn['checksum']})"
    return res


def cargar_paquete(pkg_dir: Path):
    sys.path.insert(0, str(pkg_dir))
    import alcp  # noqa: E402
    real = Path(alcp.__file__).resolve().parent.parent
    if real != pkg_dir:
        sys.exit(f"se importó alcp desde {real}, no desde {pkg_dir}")
    return alcp, alcp.alcp_native


def correr(a) -> int:
    pkg_dir = Path(a.pkg_dir).resolve()
    alcp, nat = cargar_paquete(pkg_dir)
    escenarios = filtrar(construir(), a.filter)
    label = a.label or ("alcp-" + pkg_dir.name)
    if not a.quiet:
        print(f"# bench_alcp | {label} | {pkg_dir} | {len(escenarios)} escenarios | "
              f"reps>={a.reps}, min-time {a.min_time:.2f} s", flush=True)
    t0 = time.perf_counter()
    resultados, error = [], False
    for esc in escenarios:
        rs = ejecutar(alcp, nat, esc, a)
        resultados.extend(rs)
        for r in rs:
            error |= bool(r["error"])
            if a.quiet:
                continue
            if r["error"]:
                print(f"{r['name']:<44} ERROR: {r['error']}", flush=True)
            else:
                print(f"{r['name']:<44} {r['median_s'] * 1e3:10.3f} ms  n={len(r['samples_s']):<4} "
                      f"{r['checksum']}", flush=True)
    total = time.perf_counter() - t0
    if not a.quiet:
        imprimir_sobrecoste({label: resultados})
    out = {
        "format": "alcp-bench-result", "format_version": 1, "label": label, "kind": "python",
        "mode": "single", "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "metadata": {**collect_machine_metadata(),
                     "suite": {"suite": "alcp-bench-alcp", "suite_version": 1},
                     "package": {"dir": str(pkg_dir)}},
        "note": a.note or "",
        "runs": [{"label": label, "round": 0, "total_s": total, "reps": a.reps,
                  "min_time_s": a.min_time, "warmup": a.warmup, "filter": a.filter or ""}],
        "scenarios": resultados,
    }
    p = Path(a.json)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(scrub_paths(out), indent=1, ensure_ascii=False), encoding="utf-8")
    if not a.quiet:
        print(f"-> {p}")
    return 1 if error else 0


def filtrar(escenarios: list, filtro: Optional[str]) -> list:
    if not filtro:
        return escenarios
    rx = re.compile(filtro)
    return [e for e in escenarios if rx.search(e.nombre)]


def _medianas(resultados: list, stat: str = "median_s") -> dict:
    return {r["name"]: r.get(stat) for r in resultados if r.get(stat)}


def imprimir_sobrecoste(por_label: dict, stat: str = "median_s") -> None:
    """Tabla del sobrecoste alcp/nat (y, con dos variantes, del cociente B/A)."""
    labels = list(por_label)
    med = {lb: _medianas(rs, stat) for lb, rs in por_label.items()}
    nombres = [r["name"][5:] for r in por_label[labels[0]] if r["name"].startswith("alcp.")]
    cab = "".join(f"{lb[:14]:>16}{'x nat':>8}" for lb in labels)
    extra = f"{'B/A':>8}{'sobrecoste':>12}" if len(labels) == 2 else ""
    print(f"\n{'escenario':<30}{cab}{extra}   (ms; 'x nat' = alcp / nativo; {stat[:-2]})")
    for n in nombres:
        fila, tiempos, nat_t = f"{n:<30}", [], None
        for lb in labels:
            t = med[lb].get("alcp." + n)
            tn = med[lb].get("nat." + n)
            nat_t = tn if nat_t is None else nat_t
            tiempos.append(t)
            ratio = f"{t / tn:7.2f}x" if t and tn else f"{'-':>8}"
            fila += (f"{t * 1e3:16.3f}" if t else f"{'ERR':>16}") + ratio
        if len(labels) == 2 and all(tiempos):
            fila += f"{tiempos[1] / tiempos[0]:8.3f}"
            # Reducción del sobrecoste (t_alcp - t_nat) de A a B, con el nativo medio.
            tn_a, tn_b = med[labels[0]].get("nat." + n), med[labels[1]].get("nat." + n)
            if tn_a and tn_b:
                tn = (tn_a + tn_b) / 2
                sa, sb = tiempos[0] - tn, tiempos[1] - tn
                fila += f"{(1 - sb / sa) * 100:+11.0f}%" if sa > 0 else f"{'-':>12}"
        print(fila)


def intercalar(a) -> int:
    """Ejecuta las variantes alternando en procesos nuevos y agrega las muestras."""
    variantes = []
    for v in a.interleave:
        if "=" not in v:
            sys.exit(f"--interleave espera ETIQUETA=DIRECTORIO, no {v!r}")
        lb, d = v.split("=", 1)
        variantes.append((lb, str(Path(d).resolve())))
    out_dir = Path(a.out_dir or tempfile.mkdtemp(prefix="bench_alcp_"))
    out_dir.mkdir(parents=True, exist_ok=True)
    agregados: dict = {}
    for ronda in range(a.rounds):
        orden = variantes if ronda % 2 == 0 else variantes[::-1]
        for lb, d in orden:
            js = out_dir / f"{lb}-r{ronda}.json"
            cmd = [sys.executable, str(Path(__file__).resolve()), "--pkg-dir", d, "--json", str(js),
                   "--label", lb, "--reps", str(a.reps), "--min-time", str(a.min_time),
                   "--warmup", str(a.warmup), "--quiet"]
            if a.filter:
                cmd += ["--filter", a.filter]
            print(f"# ronda {ronda + 1}/{a.rounds}: {lb}", flush=True)
            subprocess.run(cmd, check=False, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
            datos = json.loads(js.read_text(encoding="utf-8"))
            dest = agregados.setdefault(lb, {})
            for r in datos["scenarios"]:
                acum = dest.setdefault(r["name"], {**r, "samples_s": [], "checksums": set()})
                acum["samples_s"].extend(r["samples_s"])
                acum["checksums"].add(r["checksum"])
                if r["error"]:
                    acum["error"] = r["error"]
    por_label = {}
    error = False
    for lb, escs in agregados.items():
        rs = []
        for r in escs.values():
            if len(r["checksums"]) > 1:
                r["error"] = r["error"] or "checksum distinto entre rondas"
            r["checksum"] = sorted(r.pop("checksums"))[0]
            if r["samples_s"]:
                r["median_s"] = statistics.median(r["samples_s"])
                r["min_s"] = min(r["samples_s"])
            error |= bool(r["error"])
            rs.append(r)
        por_label[lb] = rs
        (out_dir / f"{lb}.json").write_text(json.dumps(scrub_paths({
            "format": "alcp-bench-result", "format_version": 1, "label": lb, "kind": "python",
            "mode": "interleave", "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "metadata": {**collect_machine_metadata(),
                         "suite": {"suite": "alcp-bench-alcp", "suite_version": 1}},
            "note": a.note or "", "runs": [], "scenarios": rs}), indent=1, ensure_ascii=False),
            encoding="utf-8")
    # Checksums de alcp.* deben coincidir entre variantes (mismos resultados).
    if len(por_label) == 2:
        (la, ra), (lb2, rb) = por_label.items()
        cb = {r["name"]: r["checksum"] for r in rb}
        for r in ra:
            if r["name"] in cb and cb[r["name"]] != r["checksum"]:
                print(f"AVISO: {r['name']} da resultados distintos en {la} y {lb2}")
                error = True
    for rs in por_label.values():
        for r in rs:
            if r["error"]:
                print(f"ERROR {r['name']}: {r['error']}")
    imprimir_sobrecoste(por_label, "min_s" if a.stat == "min" else "median_s")
    print(f"-> {out_dir}")
    return 1 if error else 0


def main(argv=None) -> int:
    utf8_stdout()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pkg-dir", default=str(BENCH_DIR.parent / "python-bindings"),
                    help="directorio que contiene el paquete alcp/ (por defecto python-bindings)")
    ap.add_argument("--json", default=None, help="ruta del JSON de salida")
    ap.add_argument("--label", default=None, help="etiqueta de la variante")
    ap.add_argument("--filter", default=None, help="expresión regular sobre el nombre del escenario")
    ap.add_argument("--reps", type=int, default=5, help="repeticiones mínimas (por defecto 5)")
    ap.add_argument("--min-time", type=float, default=0.2, help="tiempo mínimo acumulado en s")
    ap.add_argument("--warmup", type=int, default=1, help="calentamientos (por defecto 1)")
    ap.add_argument("--max-reps", type=int, default=10000, help="tope de repeticiones")
    ap.add_argument("--list", action="store_true", help="lista los escenarios y sale")
    ap.add_argument("--note", default=None, help="nota libre guardada en el JSON")
    ap.add_argument("--quiet", action="store_true", help="sin salida por pantalla")
    ap.add_argument("--interleave", nargs="+", metavar="ETIQUETA=DIR",
                    help="variantes a intercalar (cada una, un directorio con alcp/)")
    ap.add_argument("--rounds", type=int, default=3, help="rondas de --interleave (por defecto 3)")
    ap.add_argument("--out-dir", default=None, help="directorio de los JSON de --interleave")
    ap.add_argument("--stat", choices=("median", "min"), default="median",
                    help="estadístico de la tabla de --interleave (min es más robusto con ruido)")
    a = ap.parse_args(argv)
    if a.list:
        for e in construir():
            print(f"{e.nombre}\t{e.desc}")
        return 0
    if a.interleave:
        return intercalar(a)
    if not a.json:
        sys.exit("falta --json RUTA")
    return correr(a)


if __name__ == "__main__":
    sys.exit(main())
