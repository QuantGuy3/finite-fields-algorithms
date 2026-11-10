# Análisis de los resultados de la fase 3: genera tablas Markdown (mediana, IQR,
# speedup paso a paso y acumulado, checksums) a partir de los JSON de interleave.
import json, math, re, statistics, sys
from pathlib import Path

R = Path(r"<REPO>/bench/results/fase3")


def stats(xs):
    xs = sorted(xs)
    if not xs:
        return None
    if len(xs) >= 4:
        q = statistics.quantiles(xs, n=4)
        iqr = q[2] - q[0]
    else:
        iqr = xs[-1] - xs[0]
    m = statistics.median(xs)
    return {"med": m, "iqr": iqr, "rel": iqr / m if m else 0.0, "min": xs[0], "n": len(xs)}


def load_dir(d: Path, label: str):
    """Devuelve {escenario: {'st':stats,'ck':checksum,'err':error}}; fusiona rondas sueltas."""
    f = d / f"{label}.json"
    if (d / "limpio" / f"{label}.json").is_file():  # rondas sin suspensión
        f = d / "limpio" / f"{label}.json"
    files = [f] if f.is_file() else sorted(d.glob(f"{label}-r*.json"))
    extra = Path(str(d) + "_extra") / f"{label}.json"  # rondas adicionales de la misma condición
    if extra.is_file():
        files.append(extra)
    out = {}
    for p in files:
        js = json.loads(p.read_text(encoding="utf-8"))
        for s in js["scenarios"]:
            e = out.setdefault(s["name"], {"samples": [], "cks": set(), "err": None})
            e["samples"] += s["samples_s"]
            e["cks"].add(s["checksum"])
            if s.get("error"):
                e["err"] = s["error"]
    for n, e in out.items():
        e["st"] = stats(e["samples"])
        e["ck"] = sorted(e["cks"])[0] if len(e["cks"]) == 1 else "VARIABLE"
        if e["ck"] == "VARIABLE" or str(e["ck"]).startswith(("VARIABLE", "INVALIDO")):
            e["err"] = e["err"] or "checksum variable/invalido"
    return out


def gm(xs):
    xs = [x for x in xs if x and x > 0]
    return math.exp(sum(map(math.log, xs)) / len(xs)) if xs else None


def fx(x):
    return "—" if x is None else (f"{x:.2f}x" if x < 100 else f"{x:.0f}x")


def ft(s):
    if s is None:
        return "—"
    s = s * 1e3
    return f"{s:.3f}" if s < 1 else (f"{s:.2f}" if s < 100 else f"{s:.0f}")


def checks(ref, var, names=None):
    """(n comparados, n distintos, n con error) de var frente a ref."""
    names = names or [n for n in var if n in ref]
    dif = [n for n in names if n in ref and n in var and ref[n]["ck"] != var[n]["ck"]]
    err = [n for n in names if n in var and var[n]["err"]]
    return len(names), dif, err


def grupo(n):
    p = n.split(".")
    return ".".join(p[:2])


def tabla_paso(prev, cur, names, lp, lc, ref_factor=None):
    """Tabla por escenario: mediana e IQR de prev y cur, speedup paso y acumulado."""
    rows = ["| escenario | " + lp + " mediana ms (IQR %) | " + lc + " mediana ms (IQR %) | paso | vs ref |",
            "|---|---:|---:|---:|---:|"]
    sp = []
    for n in names:
        a, b = prev.get(n), cur.get(n)
        if not a or not b:
            continue
        s = a["st"]["med"] / b["st"]["med"]
        sp.append((n, s))
        vr = ref_factor(n, b) if ref_factor else None
        mark = " (!)" if max(a["st"]["rel"], b["st"]["rel"]) > 0.10 else ""
        rows.append(f"| `{n}` | {ft(a['st']['med'])} ({a['st']['rel']*100:.0f}) | {ft(b['st']['med'])} "
                    f"({b['st']['rel']*100:.0f}){mark} | {fx(s)} | {fx(vr)} |")
    return rows, sp


if __name__ == "__main__":
    print("uso: importar desde informe.py")
