#!/usr/bin/env python3
"""informe_comparativa.py — Genera docs/optimizacion/comparativa-tabla.md desde el JSON de compare_libs.py.

    python bench/informe_comparativa.py bench/results/comparativa/comparativa.json docs/optimizacion/comparativa-tabla.md
"""
import json
import math
import sys
from pathlib import Path

ORDEN = ["alcp_original", "alcp_fase3", "alcp_native", "alcp_native_1h", "alcp", "flint", "galois", "sympy", "sympy_flint"]
CORTO = {"alcp_original": "alcp original", "alcp_fase3": "alcp fase 3, 1 hilo", "alcp_native": "alcp fase 4", "alcp_native_1h": "alcp fase 4, 1 hilo",
         "alcp": "alcp paquete", "flint": "python-flint", "galois": "galois", "sympy": "SymPy puro",
         "sympy_flint": "SymPy + FLINT"}


def fmt(s):
    if s is None:
        return "—"
    return f"{s * 1e6:.1f} µs" if s < 1e-3 else (f"{s * 1e3:.2f} ms" if s < 1 else f"{s:.2f} s")


def geo(xs):
    xs = [x for x in xs if x]
    return math.exp(sum(map(math.log, xs)) / len(xs)) if xs else None


def main(entrada, salida):
    d = json.loads(Path(entrada).read_text(encoding="utf-8"))
    E = d["escenarios"]
    t = lambda e, l: e["resultados"][l].get("min_s") if e["resultados"][l]["estado"] == "ok" else None
    o = ["# Comparativa con otras librerías: tabla completa", "",
         f"Generado por `bench/informe_comparativa.py` desde `{Path(entrada).as_posix()}`. "
         "Tiempo mínimo por operación (la construcción del cuerpo y de las entradas no se cronometra).", "",
         "| Escenario | " + " | ".join(CORTO[l] for l in ORDEN) + " | Resultados iguales |",
         "|---|" + "---:|" * len(ORDEN) + ":---:|"]
    for e in E:
        celdas = []
        for l in ORDEN:
            r = e["resultados"][l]
            celdas.append(fmt(r.get("min_s")) if r["estado"] == "ok" else
                          {"no_soportado": "n/s", "error": "error", "tiempo_agotado": "> límite"}.get(r["estado"], r["estado"]))
        o.append(f"| `{e['escenario']}` | " + " | ".join(celdas) + f" | {'sí' if e['coinciden'] else '**NO**'} |")
    o += ["", "n/s: la librería no ofrece la operación (SymPy no tiene polinomios sobre F_(p^n)). "
          "error: SymPy 1.14 sobre python-flint 0.9 lanza `TypeError` al factorizar con p de 127 bits.", "",
          "## Cociente frente a python-flint (media geométrica; > 1 = más lento que FLINT)", "",
          "| Librería | Escenarios | Media geométrica | Mejor | Peor |", "|---|---:|---:|---:|---:|"]
    for l in ORDEN:
        if l == "flint":
            continue
        q = [(t(e, l) / t(e, "flint"), e["escenario"]) for e in E if t(e, l) and t(e, "flint")]
        g = geo([x for x, _ in q])
        mejor, peor = min(q), max(q)
        o.append(f"| {CORTO[l]} | {len(q)} | {g:.2f}x | {mejor[0]:.2f}x (`{mejor[1]}`) | {peor[0]:.1f}x (`{peor[1]}`) |")
    o += ["", "## Versiones y máquina", ""]
    for l in ORDEN:
        o.append(f"- {CORTO[l]}: {d['versiones'][l]}")
    m = d["maquina"]
    o.append(f"- Máquina: {m.get('cpu')} ({m.get('physical_cores')} núcleos / {m.get('logical_cpus')} hilos), "
             f"{m.get('os')}, Python {m.get('python')}, "
             + ("enchufado" if (m.get("power") or {}).get("ac_online") else "con batería")
             + f", plan de energía «{(m.get('power_plan') or '').split('(')[-1].rstrip(') ')}»")
    Path(salida).write_text("\n".join(o) + "\n", encoding="utf-8")
    print("\n".join(o))


if __name__ == "__main__":
    main(*sys.argv[1:3])
