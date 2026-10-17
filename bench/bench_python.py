#!/usr/bin/env python3
"""bench_python.py — Benchmarks a nivel Python del módulo nativo ``alcp_native``.

Mide el coste visto desde Python (conversión de tipos, sobrecarga por llamada y
operaciones completas) usando SOLO la API del binding ORIGINAL, de modo que el
mismo script sirve para el módulo de la línea base y para el mejorado::

    python bench/bench_python.py --module-dir build/bench/base/module --json bench/results/py-base.json
    python bench/bench_python.py --module-dir DIR --filter "ntt" --reps 3 --json salida.json
    python bench/bench_python.py --module-dir DIR --list

Metodología (igual que la suite C++): entradas deterministas generadas con
``random.Random(semilla fija)`` antes de cronometrar, 1 calentamiento, repeticiones
hasta alcanzar a la vez ``--reps`` y ``--min-time``, recolector de basura
desactivado durante cada medición, todas las muestras guardadas y un checksum
del resultado por escenario.  El JSON tiene el mismo formato que el de
``run_bench.py``, así que se puede comparar con ``run_bench.py compare``
pasando las rutas de los JSON.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import importlib
import json
import random
import re
import statistics
import sys
import time
from pathlib import Path
from typing import Any, Callable

BENCH_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(BENCH_DIR))
# Metadatos de máquina, anonimización de rutas y salida UTF-8 compartidos con run_bench.py.
from run_bench import collect_machine_metadata, scrub_paths, utf8_stdout  # noqa: E402

# Módulo irreducible fijo de F_(3^30) (el mismo que usa bench_suite.cpp).
MOD_3_30 = [2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1]
P_NTT = 998244353
P_M521 = 2**521 - 1


def _digest(obj: Any) -> str:
    """Checksum estable (16 hex) de una estructura de enteros/listas/tuplas."""
    return hashlib.sha256(repr(obj).encode()).hexdigest()[:16]


def _rand_coeffs(rng: random.Random, bound: int, deg: int) -> list[int]:
    c = [rng.randrange(bound) for _ in range(deg)]
    c.append(rng.randrange(1, bound))  # coeficiente principal no nulo
    return c


def _factor_digest(fac: Any) -> str:
    factors, leading = fac
    items = sorted((tuple(u.coeffs()), int(e)) for u, e in factors)
    return _digest((items, leading))


class Scenario:
    def __init__(self, name: str, group: str, op: str, desc: str, inner: int,
                 setup: Callable[[Any, random.Random], tuple[Callable[[], Any], Callable[[Any], str]]]):
        self.name, self.group, self.op, self.desc, self.inner, self.setup = name, group, op, desc, inner, setup


def _repeat(fn: Callable[[], Any], inner: int) -> Callable[[], Any]:
    if inner == 1:
        return fn

    def run() -> Any:
        r = None
        for _ in range(inner):
            r = fn()
        return r

    return run


def build_scenarios() -> list[Scenario]:
    S: list[Scenario] = []

    def add(name: str, group: str, op: str, desc: str, inner: int = 1):
        def deco(setup):
            S.append(Scenario(name, group, op, desc, inner, setup))
            return setup
        return deco

    # ---- construcción y conversión de polinomios en F_998244353 ----
    @add("py.ntt.poly_from_list.d1000", "py.ntt", "poly_from_list",
         "Fp(998244353).poly(lista) con grado 1000 (x20)", 20)
    def _(m, rng):
        F = m.Fp(P_NTT)
        lst = _rand_coeffs(rng, P_NTT, 1000)
        return _repeat(lambda: F.poly(lst), 20), lambda P: _digest(P.coeffs())

    @add("py.ntt.coeffs.d1000", "py.ntt", "coeffs", "PolyFp.coeffs() con grado 1000 (x40)", 40)
    def _(m, rng):
        F = m.Fp(P_NTT)
        P = F.poly(_rand_coeffs(rng, P_NTT, 1000))
        return _repeat(P.coeffs, 40), _digest

    @add("py.ntt.roundtrip.d1000", "py.ntt", "roundtrip",
         "ida y vuelta F.poly(P.coeffs()) con grado 1000 (x20)", 20)
    def _(m, rng):
        F = m.Fp(P_NTT)
        P = F.poly(_rand_coeffs(rng, P_NTT, 1000))
        return _repeat(lambda: F.poly(P.coeffs()), 20), lambda Q: _digest(Q.coeffs())

    # ---- multiplicación de grado 1000 en F_998244353 ----
    for algo, inner in (("auto", 5), ("fft", 5), ("karatsuba", 2), ("school", 1)):
        def mk(algo=algo, inner=inner):
            @add(f"py.ntt.mul_{algo}.d1000", "py.ntt", f"mul_{algo}",
                 "producto de grado 1000 en F_998244353 con " + ("a * b" if algo == "auto" else f".mul(b, '{algo}')")
                 + (f" (x{inner})" if inner > 1 else ""), inner)
            def _(m, rng):
                F = m.Fp(P_NTT)
                a = F.poly(_rand_coeffs(rng, P_NTT, 1000))
                b = F.poly(_rand_coeffs(rng, P_NTT, 1000))
                fn = (lambda: a * b) if algo == "auto" else (lambda: a.mul(b, algo))
                return _repeat(fn, inner), lambda P: _digest(P.coeffs())
        mk()

    # ---- factorización de grado 40 en F_101 ----
    @add("py.p101.factor.d40", "py.p101", "factor", "factor() de un polinomio de grado 40 en F_101")
    def _(m, rng):
        F = m.Fp(101)
        f = F.poly([rng.randrange(1, 101)])
        for d, mult in ((1, 2), (2, 1), (3, 2), (5, 1), (8, 1), (15, 1)):
            u = F.poly([rng.randrange(101) for _ in range(d)] + [1])
            for _ in range(mult):
                f = f * u
        return f.factor, _factor_digest

    # ---- sobrecarga por llamada: operaciones de elemento en bucle Python ----
    @add("py.ntt.fp_mul_loop.n20000", "py.ntt", "fp_mul_loop",
         "20000 llamadas Fp.mul(x, a) en un bucle Python (p = 998244353)", 20000)
    def _(m, rng):
        F = m.Fp(P_NTT)
        a, x0 = rng.randrange(1, P_NTT), rng.randrange(1, P_NTT)

        def run():
            x = x0
            for _ in range(20000):
                x = F.mul(x, a)
            return x
        return run, _digest

    @add("py.m521.fp_mul_loop.n5000", "py.m521", "fp_mul_loop",
         "5000 llamadas Fp.mul(x, a) en un bucle Python (p = 2^521 - 1)", 5000)
    def _(m, rng):
        F = m.Fp(P_M521)
        a, x0 = rng.randrange(1, P_M521), rng.randrange(1, P_M521)

        def run():
            x = x0
            for _ in range(5000):
                x = F.mul(x, a)
            return x
        return run, _digest

    # ---- conversión de enteros grandes: grado 200 con p = 2^521 - 1 ----
    @add("py.m521.poly_from_list.d200", "py.m521", "poly_from_list",
         "Fp(2^521-1).poly(lista) con grado 200 (x20)", 20)
    def _(m, rng):
        F = m.Fp(P_M521)
        lst = _rand_coeffs(rng, P_M521, 200)
        return _repeat(lambda: F.poly(lst), 20), lambda P: _digest(P.coeffs())

    @add("py.m521.coeffs.d200", "py.m521", "coeffs", "PolyFp.coeffs() con grado 200 y p = 2^521-1 (x40)", 40)
    def _(m, rng):
        F = m.Fp(P_M521)
        P = F.poly(_rand_coeffs(rng, P_M521, 200))
        return _repeat(P.coeffs, 40), _digest

    @add("py.m521.mul_auto.d200", "py.m521", "mul_auto", "producto a * b de grado 200 con p = 2^521-1")
    def _(m, rng):
        F = m.Fp(P_M521)
        a = F.poly(_rand_coeffs(rng, P_M521, 200))
        b = F.poly(_rand_coeffs(rng, P_M521, 200))
        return (lambda: a * b), lambda P: _digest(P.coeffs())

    # ---- primalidad ----
    @add("py.is_prime.n2000", "py.misc", "is_prime", "is_prime sobre 2000 impares aleatorios de 256 bits", 2000)
    def _(m, rng):
        xs = [rng.getrandbits(256) | (1 << 255) | 1 for _ in range(2000)]
        return (lambda: [m.is_prime(x) for x in xs]), _digest

    # ---- extensiones F_q ----
    @add("py.fq9.poly_from_list.d500", "py.fq", "poly_from_list",
         "Fq(3, [1,0,1]).poly(índices) con grado 500 (x40)", 40)
    def _(m, rng):
        K = m.Fq(3, [1, 0, 1])
        lst = _rand_coeffs(rng, 9, 500)
        return _repeat(lambda: K.poly(lst), 40), lambda P: _digest(P.coeffs())

    @add("py.fq9.mul_auto.d500", "py.fq", "mul_auto", "producto a * b de grado 500 en F_9 = F_3[x]/(x^2+1)")
    def _(m, rng):
        K = m.Fq(3, [1, 0, 1])
        a = K.poly(_rand_coeffs(rng, 9, 500))
        b = K.poly(_rand_coeffs(rng, 9, 500))
        return (lambda: a * b), lambda P: _digest(P.coeffs())

    @add("py.fq.q3e30.factor.d5", "py.fq", "factor",
         "factor() de un polinomio aleatorio de grado 5 en F_(3^30) (escenario histórico)")
    def _(m, rng):
        K = m.Fq(3, MOD_3_30)
        f = K.poly(_rand_coeffs(rng, 3**30, 5))
        return f.factor, _factor_digest

    # ---- algoritmos rápidos ----
    @add("py.ntt.dft.k10", "py.ntt", "dft", "alcp_native.dft de tamaño 2^10 en F_998244353 (x5)", 5)
    def _(m, rng):
        F = m.Fp(P_NTT)
        root = pow(3, (P_NTT - 1) >> 10, P_NTT)  # raíz FIJA de orden 2^10
        a = [rng.randrange(P_NTT) for _ in range(1 << 10)]
        return _repeat(lambda: m.dft(F, root, 10, a), 5), _digest

    @add("py.ntt.toeplitz_lower_inverse.n1024", "py.ntt", "toeplitz_lower_inverse",
         "alcp_native.toeplitz_lower_inverse con n = 1024 en F_998244353")
    def _(m, rng):
        F = m.Fp(P_NTT)
        a = [rng.randrange(1, P_NTT)] + [rng.randrange(P_NTT) for _ in range(1023)]
        return (lambda: m.toeplitz_lower_inverse(F, 1024, a)), _digest

    return S


def _seed(name: str) -> int:
    return int.from_bytes(hashlib.sha256(name.encode()).digest()[:8], "little")


def run_scenario(m: Any, sc: Scenario, reps: int, min_time: float, warmup: int, max_reps: int) -> dict:
    res: dict[str, Any] = {"name": sc.name, "group": sc.group, "op": sc.op, "desc": sc.desc,
                           "inner": sc.inner, "samples_s": [], "warmup_s": [], "checksum": "", "error": ""}
    try:
        t0 = time.perf_counter()
        fn, digest = sc.setup(m, random.Random(_seed(sc.name)))
        res["setup_s"] = time.perf_counter() - t0
        out = None
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
                    res["warmup_s"].append(dt)
                    continue
                res["samples_s"].append(dt)
                if len(res["samples_s"]) >= reps and sum(res["samples_s"]) >= min_time:
                    break
        finally:
            if gc_was:
                gc.enable()
        res["checksum"] = digest(out)
    except Exception as e:  # se registra y se sigue con el resto
        res["error"] = f"{type(e).__name__}: {e}"
    if res["samples_s"]:
        res["median_s"] = statistics.median(res["samples_s"])
        res["min_s"] = min(res["samples_s"])
    return res


def main(argv: list[str] | None = None) -> int:
    utf8_stdout()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--module-dir", required=True, help="directorio que contiene alcp_native*.pyd|so")
    ap.add_argument("--json", default=None, help="ruta del JSON de salida")
    ap.add_argument("--label", default=None, help="etiqueta de la variante (por defecto py-<directorio>)")
    ap.add_argument("--filter", default=None, help="expresión regular sobre el nombre del escenario")
    ap.add_argument("--reps", type=int, default=7, help="repeticiones mínimas (por defecto 7)")
    ap.add_argument("--min-time", type=float, default=0.3, help="tiempo acumulado mínimo en s (por defecto 0.3)")
    ap.add_argument("--warmup", type=int, default=1, help="calentamientos (por defecto 1)")
    ap.add_argument("--max-reps", type=int, default=10000, help="tope de repeticiones por escenario")
    ap.add_argument("--list", action="store_true", help="lista los escenarios y sale")
    ap.add_argument("--note", default=None, help="nota libre guardada en el JSON")
    a = ap.parse_args(argv)

    scenarios = build_scenarios()
    if a.filter:
        rx = re.compile(a.filter)
        scenarios = [s for s in scenarios if rx.search(s.name)]
    if a.list:
        for s in scenarios:
            print(f"{s.name}\t{s.group}\t{s.op}\t{s.desc}")
        return 0
    if not scenarios:
        sys.exit(f"ningún escenario coincide con el filtro {a.filter!r}")
    if not a.json:
        sys.exit("falta --json RUTA")

    mdir = Path(a.module_dir).resolve()
    if not mdir.is_dir():
        sys.exit(f"no existe el directorio {mdir}")
    sys.path.insert(0, str(mdir))
    m = importlib.import_module("alcp_native")
    mfile = Path(m.__file__).resolve()
    if mfile.parent != mdir:
        sys.exit(f"se importó alcp_native desde {mfile}, no desde {mdir}")

    label = a.label or ("py-" + (mdir.parent.name if mdir.name == "module" else mdir.name))

    build_info = None
    bi = mdir / "build_info.json"
    if bi.is_file():
        build_info = json.loads(bi.read_text(encoding="utf-8"))
    print(f"# bench_python | {label} | {mfile.name} | {len(scenarios)} escenarios | "
          f"reps>={a.reps}, min-time {a.min_time:.2f} s", flush=True)

    t_all = time.perf_counter()
    results = []
    any_error = False
    for sc in scenarios:
        r = run_scenario(m, sc, max(1, a.reps), max(0.0, a.min_time), max(0, a.warmup), max(1, a.max_reps))
        results.append(r)
        if r["error"]:
            any_error = True
            print(f"{sc.name:<42} ERROR: {r['error']}", flush=True)
        else:
            print(f"{sc.name:<42} {r['median_s'] * 1e3:11.3f} ms  min {r['min_s'] * 1e3:11.3f} ms  "
                  f"n={len(r['samples_s']):<5} {r['checksum']}", flush=True)
    total = time.perf_counter() - t_all
    print(f"# total {total:.2f} s")

    out = {
        "format": "alcp-bench-result", "format_version": 1, "label": label, "kind": "python", "mode": "single",
        "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "metadata": {**collect_machine_metadata(), "build": build_info,
                     "module": {"file": str(mfile),
                                "sha256": hashlib.sha256(mfile.read_bytes()).hexdigest()},
                     "suite": {"suite": "alcp-bench-python", "suite_version": 1}},
        "note": a.note or "",
        "runs": [{"label": label, "round": 0, "total_s": total, "reps": a.reps, "min_time_s": a.min_time,
                  "warmup": a.warmup, "filter": a.filter or ""}],
        "scenarios": results,
    }
    p = Path(a.json)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(scrub_paths(out), indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"-> {p}")
    return 1 if any_error else 0


if __name__ == "__main__":
    sys.exit(main())
