#!/usr/bin/env python3
"""compare_libs.py — Comparativa de velocidad de alcp con otras librerías de cuerpos finitos.

Mide las MISMAS operaciones sobre las MISMAS entradas en cada librería y comprueba
que todas devuelven el mismo resultado (forma canónica: coeficientes como índices
en base p, factores mónicos ordenados)::

    python bench/compare_libs.py --json bench/results/comparativa.json
    python bench/compare_libs.py --filter fp101 --libs alcp_native,flint
    python bench/compare_libs.py --list

Librerías (las que no estén instaladas se omiten y se indica):

* ``alcp_original``: el binding del código ORIGINAL (etiqueta ``baseline-original``),
  para ver cuánto ha avanzado la fase 3.
* ``alcp_nivel_a``: el binding de la fase 4, nivel A (rama ``fase4-nivel-a``), con un hilo.
* ``alcp_nivel_b``: el binding de la fase 4, nivel B (commit ``7880a88``), con un hilo.
* ``alcp_fase3``: el binding de la fase 3 (rama ``mejoras-cpp-pybind`` antes de la
  fase 4), con un hilo, para medir lo que aporta la fase 4 frente a ``alcp_native_1h``.
* ``alcp_native``: el binding pybind11 directo del núcleo C++ de este repositorio.
* ``alcp``: el paquete Python de alto nivel (``GF``, ``K["x"]``) sobre el anterior.
* ``flint``: python-flint (FLINT 3, en C).  ``nmod_poly`` cuando p cabe en 64 bits
  (salvo el test de irreducibilidad, que solo existe en ``fmpz_mod_poly``),
  ``fmpz_mod_poly`` para p grande y ``fq_default_poly`` para extensiones.
* ``galois``: galois (NumPy + numba).
* ``sympy``: SymPy con tipos base de Python puro (``SYMPY_GROUND_TYPES=python``).
* ``sympy_flint``: SymPy cuando usa python-flint por debajo (``DUP_Flint``).

Metodología: cada par (librería, escenario) se ejecuta en un PROCESO NUEVO, con
tiempo límite.  La construcción del cuerpo y la conversión de las entradas quedan
FUERA del cronómetro (galois, por ejemplo, compila con numba al crear el cuerpo);
solo se mide la operación.  1 calentamiento y repeticiones hasta reunir a la vez
``--reps`` muestras y ``--min-time`` segundos (máximo ``--max-reps``); se informa
del mínimo y de la mediana.  Recolector de basura desactivado durante la medida.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import random
import statistics
import subprocess
import sys
import time
from pathlib import Path

BENCH_DIR = Path(__file__).resolve().parent
REPO = BENCH_DIR.parent
sys.path.insert(0, str(BENCH_DIR))
sys.path.insert(0, str(REPO / "python-bindings"))

# Módulos de versiones anteriores, compilados con ``run_bench.py module``:
#   final-mod-ref: cabeceras y bindings.cpp de baseline-original (bench/results/fase3/scripts/final.sh)
#   fase3-mod:     cabeceras y bindings.cpp de mejoras-cpp-pybind antes de la fase 4
#   nivelA-mod:    cabeceras y bindings.cpp de la rama fase4-nivel-a
MODULOS_EXTERNOS = {
    "alcp_original": (REPO / "build" / "bench" / "final-mod-ref" / "module").resolve(),
    "alcp_fase3": (REPO / "build" / "bench" / "fase3-mod" / "module").resolve(),
    "alcp_nivel_a": (REPO / "build" / "bench" / "nivelA-mod" / "module").resolve(),
    "alcp_nivel_b": (REPO / "build" / "bench" / "nivelB-mod" / "module").resolve(),
}
P_NTT = 998244353
P_M127 = 2**127 - 1
# Módulo de AES para F_(2^8) y módulo fijo de F_(3^30) (el de bench_suite.cpp).
MOD_2_8 = [1, 1, 0, 1, 1, 0, 0, 0, 1]
MOD_3_30 = [2, 0, 2, 0, 1, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 2, 2, 1, 0, 1, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1]

LIBRERIAS = ["alcp_original", "alcp_fase3", "alcp_nivel_a", "alcp_nivel_b", "alcp_native", "alcp_native_1h", "alcp", "flint", "galois", "sympy", "sympy_flint"]
NOMBRES = {
    "alcp_original": "alcp original (antes de la fase 3)",
    "alcp_fase3": "alcp fase 3, 1 hilo",
    "alcp_nivel_a": "alcp fase 4 nivel A, 1 hilo",
    "alcp_nivel_b": "alcp fase 4 nivel B, 1 hilo",
    "alcp_native": "alcp (binding C++)",
    "alcp_native_1h": "alcp (binding C++, 1 hilo)",
    "alcp": "alcp (paquete Python)",
    "flint": "python-flint",
    "galois": "galois",
    "sympy": "SymPy (Python puro)",
    "sympy_flint": "SymPy (sobre FLINT)",
}


# --------------------------------------------------------------------------- entradas
def _aleatorio(rng: random.Random, q: int, grado: int, monico: bool = False) -> list[int]:
    """Coeficientes (índices en base p, orden creciente) de un polinomio de grado dado."""
    c = [rng.randrange(q) for _ in range(grado)]
    c.append(1 if monico else rng.randrange(1, q))
    return c


def _semilla(nombre: str) -> int:
    return int(hashlib.sha256(nombre.encode()).hexdigest()[:8], 16)


class Escenario:
    """Una operación sobre un cuerpo concreto con entradas deterministas."""

    def __init__(self, nombre, grupo, desc, p, n, modulo, op, entradas):
        self.nombre, self.grupo, self.desc = nombre, grupo, desc
        self.p, self.n, self.modulo, self.op = p, n, modulo, op
        self._entradas = entradas

    @property
    def q(self) -> int:
        return self.p ** self.n

    def entradas(self) -> list[list[int]]:
        return self._entradas(random.Random(_semilla(self.nombre)), self.q)


def _producto_de_factores(nat_poly, rng, q, grados):
    """Producto mónico de factores aleatorios con multiplicidades (entrada de factor())."""
    f = nat_poly([1])
    for d, mult in grados:
        u = nat_poly(_aleatorio(rng, q, d, monico=True))
        for _ in range(mult):
            f = f * u
    return f.coeffs()


def construir_escenarios() -> list[Escenario]:
    E: list[Escenario] = []

    def add(nombre, grupo, desc, p, n, modulo, op):
        def deco(fn):
            E.append(Escenario(nombre, grupo, desc, p, n, modulo, op, fn))
            return fn
        return deco

    # ---- F_998244353 (primo de 30 bits, apto para NTT) ----
    for d in (100, 1000, 10000):
        @add(f"fp30.mul.d{d}", "F_998244353", f"producto de dos polinomios de grado {d}", P_NTT, 1, None, "mul")
        def _(rng, q, d=d):
            return [_aleatorio(rng, q, d), _aleatorio(rng, q, d)]

    @add("fp30.divmod.d4000", "F_998244353", "división con resto de grado 4000 entre grado 2000", P_NTT, 1, None, "divmod")
    def _(rng, q):
        return [_aleatorio(rng, q, 4000), _aleatorio(rng, q, 2000)]

    @add("fp30.gcd.d1000", "F_998244353", "mcd de dos polinomios de grado 1000", P_NTT, 1, None, "gcd")
    def _(rng, q):
        return [_aleatorio(rng, q, 1000), _aleatorio(rng, q, 1000)]

    # Grado alto (fase 4): aquí actúa la división de Newton con el inverso amortizado.
    @add("fp30.irred.d400", "F_998244353", "test de irreducibilidad de un polinomio IRREDUCIBLE de grado 400",
         P_NTT, 1, None, "irred")
    def _(rng, q):
        from alcp import alcp_native as nat
        F = nat.Fp(P_NTT)
        while True:
            c = _aleatorio(rng, q, 400, monico=True)
            if F.poly(c).is_irreducible():
                return [c]

    @add("fp30.factor.d300", "F_998244353", "factor() de un polinomio mónico aleatorio de grado 300",
         P_NTT, 1, None, "factor")
    def _(rng, q):
        return [_aleatorio(rng, q, 300, monico=True)]

    # ---- F_101 (primo pequeño): factorización e irreducibilidad ----
    @add("fp101.factor.d40", "F_101", "factor() de un producto conocido de grado 40 con multiplicidades",
         101, 1, None, "factor")
    def _(rng, q):
        from alcp import alcp_native as nat
        F = nat.Fp(101)
        return [_producto_de_factores(F.poly, rng, q, ((1, 2), (2, 1), (3, 2), (5, 1), (8, 1), (15, 1)))]

    @add("fp101.factor.d100", "F_101", "factor() de un polinomio mónico aleatorio de grado 100", 101, 1, None, "factor")
    def _(rng, q):
        return [_aleatorio(rng, q, 100, monico=True)]

    @add("fp101.irred.d200", "F_101", "test de irreducibilidad de un polinomio IRREDUCIBLE de grado 200",
         101, 1, None, "irred")
    def _(rng, q):
        from alcp import alcp_native as nat
        F = nat.Fp(101)
        while True:          # el primero irreducible de la sucesión determinista
            c = _aleatorio(rng, q, 200, monico=True)
            if F.poly(c).is_irreducible():
                return [c]

    # ---- p = 2^127 - 1 (enteros de varias palabras) ----
    @add("fp127.mul.d500", "F_(2^127-1)", "producto de dos polinomios de grado 500", P_M127, 1, None, "mul")
    def _(rng, q):
        return [_aleatorio(rng, q, 500), _aleatorio(rng, q, 500)]

    @add("fp127.factor.d20", "F_(2^127-1)", "factor() de un polinomio mónico aleatorio de grado 20",
         P_M127, 1, None, "factor")
    def _(rng, q):
        return [_aleatorio(rng, q, 20, monico=True)]

    # ---- extensiones ----
    @add("f2e8.mul.d1000", "F_(2^8)", "producto de grado 1000 en F_256 (módulo de AES)", 2, 8, MOD_2_8, "mul")
    def _(rng, q):
        return [_aleatorio(rng, q, 1000), _aleatorio(rng, q, 1000)]

    @add("f2e8.factor.d30", "F_(2^8)", "factor() de un polinomio mónico aleatorio de grado 30 en F_256",
         2, 8, MOD_2_8, "factor")
    def _(rng, q):
        return [_aleatorio(rng, q, 30, monico=True)]

    @add("f3e30.mul.d200", "F_(3^30)", "producto de grado 200 en F_(3^30)", 3, 30, MOD_3_30, "mul")
    def _(rng, q):
        return [_aleatorio(rng, q, 200), _aleatorio(rng, q, 200)]

    @add("f3e30.factor.d5", "F_(3^30)", "factor() de un polinomio mónico aleatorio de grado 5 en F_(3^30)",
         3, 30, MOD_3_30, "factor")
    def _(rng, q):
        return [_aleatorio(rng, q, 5, monico=True)]

    return E


# --------------------------------------------------------------------------- adaptadores
# Cada adaptador devuelve (función a cronometrar, función que normaliza su resultado).
# Forma canónica: polinomio = lista de índices en orden creciente; divmod = (q, r);
# factor = lista ordenada de (factor mónico, exponente); irred = bool.

def _indice(coefs, p: int) -> int:
    return sum(int(c) * p**i for i, c in enumerate(coefs))


def _digitos(i: int, p: int, n: int) -> list[int]:
    d = []
    for _ in range(n):
        i, r = divmod(i, p)
        d.append(r)
    return d


def _factores(pares):
    return sorted((list(f), int(e)) for f, e in pares)


def adaptador_alcp_native(sc: Escenario, args, nat=None):
    if nat is None:
        from alcp import alcp_native as nat
    K = nat.Fp(sc.p) if sc.n == 1 else nat.Fq(sc.p, sc.modulo)
    a = [K.poly(c) for c in args]
    canon = lambda P: P.coeffs()
    return _operar(sc.op, a, canon,
                   factor=lambda f: f.factor(),
                   norm_factor=lambda r: _factores((u.coeffs(), e) for u, e in r[0]),
                   irred=lambda f: f.is_irreducible(),
                   gcd=lambda f, g: f.gcd(g), divmod_=lambda f, g: f.divmod(g))


def adaptador_externo(lib: str):
    """Adaptador para un módulo de otra versión (MODULOS_EXTERNOS), con la misma API básica."""
    def adaptador(sc: Escenario, args):
        import importlib
        carpeta = MODULOS_EXTERNOS[lib]
        sys.path.insert(0, str(carpeta))
        nat = importlib.import_module("alcp_native")
        if not Path(nat.__file__).resolve().is_relative_to(carpeta):
            raise ImportError("no se encontró el módulo en " + str(carpeta))
        return adaptador_alcp_native(sc, args, nat)
    return adaptador


def adaptador_alcp(sc: Escenario, args):
    import alcp
    K = alcp.GF(sc.p) if sc.n == 1 else alcp.GF(sc.q, modulus=sc.modulo)
    R = K["x"]
    a = [R([K.from_index(i) for i in c]) for c in args]
    canon = lambda P: [K.to_index(c) for c in P.coeffs]   # propiedad en el paquete
    return _operar(sc.op, a, canon,
                   factor=lambda f: f.factor(),
                   norm_factor=lambda r: _factores((canon(u), e) for u, e in r),
                   irred=lambda f: f.is_irreducible(),
                   gcd=lambda f, g: f.gcd(g), divmod_=divmod)


def adaptador_flint(sc: Escenario, args):
    import flint
    p, n = sc.p, sc.n
    if n == 1:
        if p < 2**63 and sc.op != "irred":
            a = [flint.nmod_poly(c, p) for c in args]
        else:
            R = flint.fmpz_mod_poly_ctx(p)
            a = [R(c) for c in args]
        canon = lambda P: [int(c) for c in P.coeffs()]
    else:
        ctx = flint.fq_default_ctx(p, n, "a", modulus=flint.fmpz_mod_poly_ctx(p)(sc.modulo))
        R = flint.fq_default_poly_ctx(ctx)
        a = [R([ctx(_digitos(i, p, n)) for i in c]) for c in args]
        canon = lambda P: [_indice([int(t) for t in c.to_list()], p) for c in P.coeffs()]
    return _operar(sc.op, a, canon,
                   factor=lambda f: f.factor(),
                   norm_factor=lambda r: _factores((canon(u), e) for u, e in r[1]),
                   irred=lambda f: f.is_irreducible(),
                   gcd=lambda f, g: f.gcd(g), divmod_=divmod)


def adaptador_galois(sc: Escenario, args):
    import galois
    if sc.n == 1:
        GF = galois.GF(sc.p)
    else:
        GF = galois.GF(sc.q, irreducible_poly=galois.Poly(sc.modulo[::-1], field=galois.GF(sc.p)))
    a = [galois.Poly(GF(c[::-1])) for c in args]    # galois: orden decreciente
    canon = lambda P: [int(c) for c in P.coeffs[::-1]] if P.degree > 0 or int(P.coeffs[0]) else []
    return _operar(sc.op, a, canon,
                   factor=lambda f: f.factors(),
                   norm_factor=lambda r: _factores((canon(u), e) for u, e in zip(*r)),
                   # galois guarda is_irreducible en una functools.lru_cache: sin __wrapped__
                   # solo se mediría la consulta a la caché tras el calentamiento.
                   irred=lambda f: galois.Poly.is_irreducible.__wrapped__(f),
                   gcd=lambda f, g: galois.gcd(f, g), divmod_=divmod)


def adaptador_sympy(sc: Escenario, args):
    if sc.n != 1:
        raise NoSoportado("SymPy no tiene polinomios sobre extensiones F_(p^n)")
    import sympy
    x = sympy.Symbol("x")
    p = sc.p
    a = [sympy.Poly(c[::-1], x, modulus=p) for c in args]
    canon = lambda P: [int(c) % p for c in P.all_coeffs()[::-1]] if not P.is_zero else []
    return _operar(sc.op, a, canon,
                   factor=lambda f: f.factor_list(),
                   norm_factor=lambda r: _factores((canon(u), e) for u, e in r[1]),
                   irred=lambda f: f.is_irreducible,
                   gcd=lambda f, g: f.gcd(g), divmod_=divmod)


ADAPTADORES = {
    "alcp_original": adaptador_externo("alcp_original"),
    "alcp_fase3": adaptador_externo("alcp_fase3"),
    "alcp_nivel_a": adaptador_externo("alcp_nivel_a"),
    "alcp_nivel_b": adaptador_externo("alcp_nivel_b"),
    "alcp_native": adaptador_alcp_native,
    "alcp_native_1h": adaptador_alcp_native,
    "alcp": adaptador_alcp,
    "flint": adaptador_flint,
    "galois": adaptador_galois,
    "sympy": adaptador_sympy,
    "sympy_flint": adaptador_sympy,
}


class NoSoportado(Exception):
    pass


def _operar(op, a, canon, *, factor, norm_factor, irred, gcd, divmod_):
    """Traduce la operación del escenario a la API de la librería."""
    if op == "mul":
        f, g = a
        return (lambda: f * g), canon
    if op == "divmod":
        f, g = a
        return (lambda: divmod_(f, g)), lambda r: (canon(r[0]), canon(r[1]))
    if op == "gcd":
        f, g = a
        return (lambda: gcd(f, g)), canon
    if op == "factor":
        (f,) = a
        return (lambda: factor(f)), norm_factor
    if op == "irred":
        (f,) = a
        return (lambda: irred(f)), bool
    raise ValueError(op)


# --------------------------------------------------------------------------- medición
def medir(lib: str, sc: Escenario, reps: int, min_time: float, max_reps: int) -> dict:
    args = sc.entradas()
    t0 = time.perf_counter()
    fn, canon = ADAPTADORES[lib](sc, args)
    preparacion = time.perf_counter() - t0
    r = fn()                                   # calentamiento (y JIT, si lo hay)
    resultado = canon(r)
    muestras: list[float] = []
    gc_activo = gc.isenabled()
    gc.disable()
    try:
        total = 0.0
        while len(muestras) < max_reps and (len(muestras) < reps or total < min_time):
            t = time.perf_counter()
            fn()
            dt = time.perf_counter() - t
            muestras.append(dt)
            total += dt
            if dt > 20 and len(muestras) >= 1:    # operación muy lenta: no insistir
                break
    finally:
        if gc_activo:
            gc.enable()
    return {
        "min_s": min(muestras),
        "mediana_s": statistics.median(muestras),
        "reps": len(muestras),
        "preparacion_s": preparacion,
        "checksum": hashlib.sha256(repr(resultado).encode()).hexdigest()[:16],
    }


def version(lib: str) -> str:
    try:
        if lib == "alcp_original":
            return "baseline-original (67331a2)"
        if lib == "alcp_fase3":
            return "mejoras-cpp-pybind antes de la fase 4 (ALCP_THREADS=1)"
        if lib == "alcp_nivel_a":
            return "fase4-nivel-a (ALCP_THREADS=1)"
        if lib == "alcp_nivel_b":
            return "fase 4 nivel B, 7880a88 (ALCP_THREADS=1)"
        if lib.startswith("alcp"):
            return "rama local" + (" (ALCP_THREADS=1)" if lib.endswith("_1h") else "")
        mod = {"flint": "flint", "galois": "galois", "sympy": "sympy", "sympy_flint": "sympy"}[lib]
        m = __import__(mod)
        v = getattr(m, "__version__", "?")
        if lib == "sympy_flint":
            import flint
            v += f" + python-flint {flint.__version__}"
        return v
    except Exception as e:      # librería ausente
        return f"no disponible ({type(e).__name__})"


def entorno(lib: str) -> dict:
    env = dict(os.environ)
    # SymPy elige sus tipos base al importarse: se fija en cada proceso.
    env["SYMPY_GROUND_TYPES"] = "flint" if lib == "sympy_flint" else "python"
    env["PYTHONIOENCODING"] = "utf-8"
    # alcp paraleliza Karatsuba y la NTT; las demás librerías usan un solo hilo.
    if lib in ("alcp_native_1h", "alcp_fase3", "alcp_nivel_a", "alcp_nivel_b"):
        env["ALCP_THREADS"] = "1"
    else:
        env.pop("ALCP_THREADS", None)
    return env


def trabajador(lib: str, nombre: str, a) -> int:
    sc = next(s for s in construir_escenarios() if s.nombre == nombre)
    try:
        if lib == "sympy_flint":
            from sympy.external.gmpy import GROUND_TYPES
            if GROUND_TYPES != "flint":
                raise NoSoportado("SymPy no está usando python-flint")
        res = {"estado": "ok", **medir(lib, sc, a.reps, a.min_time, a.max_reps)}
    except NoSoportado as e:
        res = {"estado": "no_soportado", "motivo": str(e)}
    except ImportError as e:
        res = {"estado": "no_instalada", "motivo": str(e)}
    except Exception as e:
        res = {"estado": "error", "motivo": f"{type(e).__name__}: {e}"}
    print("RESULTADO " + json.dumps(res))
    return 0


def _fmt(s: float | None) -> str:
    if s is None:
        return "—"
    if s < 1e-3:
        return f"{s * 1e6:.1f} µs"
    if s < 1:
        return f"{s * 1e3:.2f} ms"
    return f"{s:.2f} s"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--libs", default=",".join(LIBRERIAS))
    ap.add_argument("--filter", default="", help="subcadena del nombre del escenario")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--min-time", type=float, default=0.5)
    ap.add_argument("--max-reps", type=int, default=500)
    ap.add_argument("--timeout", type=float, default=180, help="segundos por (librería, escenario)")
    ap.add_argument("--json", type=Path)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--worker", nargs=2, metavar=("LIB", "ESCENARIO"), help=argparse.SUPPRESS)
    a = ap.parse_args(argv)
    if a.worker:
        return trabajador(*a.worker, a)

    escenarios = [s for s in construir_escenarios() if a.filter in s.nombre]
    libs = [l for l in a.libs.split(",") if l]
    if a.list:
        for s in escenarios:
            print(f"{s.nombre:22} {s.desc}")
        return 0

    filas = []
    for s in escenarios:
        print(f"\n== {s.nombre}: {s.desc}", flush=True)
        celdas = {}
        for lib in libs:
            cmd = [sys.executable, str(Path(__file__).resolve()), "--worker", lib, s.nombre,
                   "--reps", str(a.reps), "--min-time", str(a.min_time), "--max-reps", str(a.max_reps)]
            for _intento in range(4):
                try:
                    out = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                                         env=entorno(lib), timeout=a.timeout).stdout
                    linea = [l for l in out.splitlines() if l.startswith("RESULTADO ")]
                    res = json.loads(linea[-1][10:]) if linea else {"estado": "error", "motivo": out[-300:]}
                except subprocess.TimeoutExpired:
                    res = {"estado": "tiempo_agotado", "motivo": f"más de {a.timeout:.0f} s"}
                # Smart App Control de Windows bloquea a veces una DLL recién vista: se reintenta.
                if "Application Control" not in res.get("motivo", ""):
                    break
            celdas[lib] = res
            extra = res.get("motivo", "")
            print(f"   {NOMBRES[lib]:27} {_fmt(res.get('min_s')):>11}  {res['estado']:14} "
                  f"{res.get('checksum', '')} {extra[:80]}", flush=True)
        sumas = {r["checksum"] for r in celdas.values() if r["estado"] == "ok"}
        if len(sumas) > 1:
            print(f"   ¡ATENCIÓN! resultados distintos: {sumas}")
        filas.append({"escenario": s.nombre, "grupo": s.grupo, "descripcion": s.desc,
                      "coinciden": len(sumas) <= 1, "resultados": celdas})

    if a.json:
        from run_bench import collect_machine_metadata, scrub_paths
        datos = {"maquina": scrub_paths(collect_machine_metadata()),
                 "versiones": {l: version(l) for l in libs},
                 "parametros": {"reps": a.reps, "min_time": a.min_time, "timeout": a.timeout},
                 "escenarios": filas}
        a.json.parent.mkdir(parents=True, exist_ok=True)
        a.json.write_text(json.dumps(datos, indent=1, ensure_ascii=False), encoding="utf-8")
        print(f"\nJSON escrito en {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
