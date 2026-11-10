# Genera las tablas Markdown de la fase 3 y un resumen JSON (para el informe final).
import json, re, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from analiza import R, load_dir, gm, fx, ft, checks, grupo, tabla_paso  # noqa: E402

OUT = []
RES = {}


def p(*l):
    OUT.extend(l)


# ---------------------------------------------------------------- extremos
import os
AC = "" if os.environ.get("BAT") else "_ac"
def pick(nombre):
    d = R / (nombre + AC)
    return d if d.is_dir() and any(d.glob("*.json")) else R / nombre
E = pick("extremos")
OUT.append(f"<!-- datos: {E.name} -->")
ext = {lb: load_dir(E, lb) for lb in ("b0", "h0", "a6", "a6gmp", "m9", "m9t1")}
b0, h0x = ext["b0"], ext["h0"]
names_all = list(b0)


def vs_ref_ext(n, h0_chain_med, cur_med):
    """Speedup acumulado frente a la referencia, encadenado: (b0/h0)_extremos * (h0/cur)_cadena."""
    if n not in b0 or n not in h0x:
        return None
    return (b0[n]["st"]["med"] / h0x[n]["st"]["med"]) * (h0_chain_med / cur_med)


p("## 1. Extremos: suite completa (128 escenarios)", "")
p("| variante | escenarios | checksums = ref | errores | media geom. vs ref | escenarios con IQR > 10 % |",
  "|---|---:|---:|---:|---:|---:|")
for lb, d in ext.items():
    n, dif, err = checks(b0, d)
    sp = [b0[k]["st"]["med"] / d[k]["st"]["med"] for k in d if k in b0]
    ruid = sum(1 for k in d if d[k]["st"]["rel"] > 0.10)
    RES.setdefault("ext", {})[lb] = {"gm": gm(sp), "dif": dif, "err": err, "ruido": ruid}
    p(f"| `{lb}` | {n} | {n - len(dif)} | {len(err)} | {fx(gm(sp)) if lb != 'b0' else '1.00x (ref)'} | {ruid} |")
p("")
# Por grupo
grupos = sorted({grupo(n) for n in names_all}, key=names_all.index if False else None)
p("### Media geométrica del speedup frente a la referencia, por grupo", "")
labs = ["h0", "a6", "a6gmp", "m9", "m9t1"]
p("| grupo | n | " + " | ".join(labs) + " |", "|---|---:|" + "---:|" * len(labs))
for g in sorted({grupo(n) for n in names_all}):
    ns = [n for n in names_all if grupo(n) == g]
    cells = []
    for lb in labs:
        cells.append(fx(gm([b0[n]["st"]["med"] / ext[lb][n]["st"]["med"] for n in ns if n in ext[lb]])))
    p(f"| `{g}` | {len(ns)} | " + " | ".join(cells) + " |")
p("")
p("### Por operación", "")
ops = sorted({(n.split(".")[2] if n.count(".") >= 3 else "build") for n in names_all})
p("| operación | n | " + " | ".join(labs) + " |", "|---|---:|" + "---:|" * len(labs))
for o in ops:
    ns = [n for n in names_all if (n.split(".")[2] if n.count(".") >= 3 else "build") == o]
    cells = [fx(gm([b0[n]["st"]["med"] / ext[lb][n]["st"]["med"] for n in ns if n in ext[lb]])) for lb in labs]
    p(f"| `{o}` | {len(ns)} | " + " | ".join(cells) + " |")
p("")
# Coste de robustez: escenarios más afectados
rob = sorted(((b0[n]["st"]["med"] / h0x[n]["st"]["med"], n) for n in names_all if n in h0x))
p("### Coste de robustez (5f49217 frente a la referencia): los 10 escenarios más lentos y los 5 más rápidos", "")
p("| escenario | ref mediana ms (IQR %) | 5f49217 mediana ms (IQR %) | speedup |", "|---|---:|---:|---:|")
for s, n in rob[:10] + rob[-5:]:
    a, b = b0[n]["st"], h0x[n]["st"]
    p(f"| `{n}` | {ft(a['med'])} ({a['rel']*100:.0f}) | {ft(b['med'])} ({b['rel']*100:.0f}) | {fx(s)} |")
p("")
RES["robustez"] = {"gm": gm([s for s, _ in rob]), "min": rob[0], "max": rob[-1],
                   "n_lentos_5": sum(1 for s, _ in rob if s < 0.95), "n_rapidos_5": sum(1 for s, _ in rob if s > 1.05)}

# ---------------------------------------------------------------- cadenas
TEC = {
    "aritmetica": [
        ("h0", "a1", "489728b", "GMP en sitio en Fp y Fq", r"^fq\.|^fp\.(m127|p256)"),
        ("a1", "a2", "0de8cf9", "Reducción perezosa en mul_school/divmod (LazyAccumulation)",
         r"mul_school|divmod|pow_mod|factor|mul_karatsuba|gcd|is_irreducible|toeplitz_divmod"),
        ("a2", "a3", "dcb905c", "Suma/resta de Poly sin relleno ni copias", r"mul_karatsuba|gcd_ext|gcd"),
        ("a3", "a4", "65a5f1c (+5c444c0, 4fb4f1a)", "Backend Fp64/Fq64 (auto)", r"^fp\.(p7|p101|ntt|m61)|^fq\.|^fast\."),
        ("a3", "a4gmp", "65a5f1c, --backend gmp", "FqT<Base> con backend GMP forzado (coste de la plantilla)", r"."),
        ("a4", "a5", "2297387", "Fp64: Barrett sin desplazamientos variables", r"elem_pow|elem_muladd|^fp\.(p7|p101|ntt|m61)"),
        ("a5", "a6", "8f1d896", "Fp64::pow en forma de Montgomery", r"elem_pow"),
    ],
    "multiplicacion": [
        ("h0", "m1", "d4580af", "FFT: tablas directas, bit-reversal precalculado, mariposa sin temporal", r"mul_fft|mul_auto|fast\.ntt"),
        ("m1", "m2", "f68a18b", "Cuadrado con una sola transformada", r"pow_mod|is_irreducible"),
        ("m2", "m3", "9bfc1e7", "Karatsuba sobre vistas", r"mul_karatsuba|mul_auto"),
        ("m3", "m4", "1ad8c02", "NTT en uint64_t + Shoup (p < 2^62)", r"fp\.ntt\.|fast\.ntt\."),
        ("m4", "m5", "70e3453", "NTT multi-limb mpn/Montgomery (p >= 2^62)", r"fp\.p256\."),
        ("m5", "m6", "3660189", "Umbrales de despacho Auto", r"gcd|factor|mul_auto"),
        ("m6", "m7", "06e46b5", "fast::dft/idft en palabras de máquina", r"fast\.ntt\.(dft|idft)"),
        ("m7", "m8", "e8e07e6", "Karatsuba paralelo (std::async)", r"mul_karatsuba|mul_auto|fq\..*(gcd|factor)"),
        ("m8", "m9", "b0a9fe7", "NTT paralela", r"mul_fft|mul_auto|fq\.m61e2"),
        ("m7", "m8t1", "e8e07e6 con ALCP_THREADS=1", "Karatsuba paralelo limitado a 1 hilo (sobrecoste)", r"mul_karatsuba|mul_auto|fq\..*(gcd|factor)"),
        ("m8t1", "m9t1", "b0a9fe7 con ALCP_THREADS=1", "NTT paralela limitada a 1 hilo (sobrecoste)", r"mul_fft|mul_auto|fq\.m61e2"),
    ],
}
for area, lista in TEC.items():
    D = pick(area)
    OUT.append(f"<!-- datos: {D.name} -->")
    if not D.is_dir():
        continue
    labs = sorted({x for t in lista for x in t[:2]})
    data = {lb: load_dir(D, lb) for lb in labs}
    if any(not v for v in data.values()):
        p(f"## {area}: datos incompletos", "")
        continue
    h0c = data["h0"]
    names = list(h0c)
    p(f"## 2.{list(TEC).index(area) + 1}. Cadena {area}", "")
    # Checksums y ruido
    p("| variante | escenarios | checksums = ref (b0) | = 5f49217 | errores | IQR > 10 % |", "|---|---:|---:|---:|---:|---:|")
    for lb in labs:
        n, dif, err = checks(b0, data[lb])
        _, dif2, _ = checks(h0c, data[lb])
        ruid = sum(1 for k in data[lb] if data[lb][k]["st"]["rel"] > 0.10)
        p(f"| `{lb}` | {n} | {n - len(dif)} | {n - len(dif2)} | {len(err)} | {ruid} |")
        RES.setdefault("ck", {})[f"{area}/{lb}"] = (n, dif, dif2, err)
    p("")
    # Resumen por técnica
    p("| commit | técnica | n afect. | media geom. paso (afectados) | rango paso | media geom. paso (control) | media geom. acumulada vs ref (afectados) |",
      "|---|---|---:|---:|---:|---:|---:|")
    detalle = []
    for prev, cur, com, desc, rx in lista:
        A, B = data[prev], data[cur]
        af = [n for n in names if re.search(rx, n) and n in A and n in B]
        ct = [n for n in names if not re.search(rx, n) and n in A and n in B]
        sa = [A[n]["st"]["med"] / B[n]["st"]["med"] for n in af]
        sc = [A[n]["st"]["med"] / B[n]["st"]["med"] for n in ct]
        vr = [vs_ref_ext(n, h0c[n]["st"]["med"], B[n]["st"]["med"]) for n in af]
        vr = [x for x in vr if x]
        rng = f"{min(sa):.2f}-{max(sa):.2f}" if sa else "—"
        RES.setdefault("tec", {})[f"{area}/{cur}"] = {"commit": com, "desc": desc, "gm": gm(sa), "gm_ctl": gm(sc),
                                                     "vsref": gm(vr), "rng": rng, "n": len(af),
                                                     "slow": [(n, s) for n, s in zip(af, sa) if s < 0.95]}
        p(f"| {com} | {desc} | {len(af)} | {fx(gm(sa))} | {rng} | {fx(gm(sc))} | {fx(gm(vr))} |")
        rows, _ = tabla_paso(A, B, af, prev, cur,
                             lambda n, b: vs_ref_ext(n, h0c[n]["st"]["med"], b["st"]["med"]))
        detalle.append((com, desc, rows))
    p("")
    for com, desc, rows in detalle:
        p(f"<details><summary>{com} — {desc}: detalle por escenario</summary>", "", *rows, "", "</details>", "")

# ---------------------------------------------------------------- toolchain
D = R / "toolchain"
if D.is_dir() and (D / "th0-O3.json").is_file():
    labs = ["tb0-O3", "tb0-v3", "tb0-nat", "tb0-lto", "tb0-alloc", "th0-O3", "th0-v3", "th0-nat", "th0-lto", "th0-alloc"]
    data = {lb: load_dir(D, lb) for lb in labs}
    p("## 3. Toolchain (rama wip, perfiles de flags)", "")
    p("| cabeceras | perfil | media geom. vs -O3 | rango | checksums = ref | IQR > 10 % |", "|---|---|---:|---:|---:|---:|")
    for pre, nom in (("tb0", "baseline-original"), ("th0", "5f49217")):
        base = data[f"{pre}-O3"]
        for pf in ("v3", "nat", "lto", "alloc"):
            v = data[f"{pre}-{pf}"]
            s = [base[n]["st"]["med"] / v[n]["st"]["med"] for n in v if n in base]
            n, dif, err = checks(b0, v)
            ruid = sum(1 for k in v if v[k]["st"]["rel"] > 0.10)
            RES.setdefault("tool", {})[f"{pre}-{pf}"] = (gm(s), min(s), max(s), len(dif) + len(err))
            p(f"| {nom} | {pf} | {fx(gm(s))} | {min(s):.2f}-{max(s):.2f} | {n - len(dif)}/{n} | {ruid} |")
    p("")

# ---------------------------------------------------------------- módulo
D = R / "modulo"
if D.is_dir() and list(D.glob("py-modm9-r*.json")):
    labs = ["py-modb0", "py-modh0", "py-moda6", "py-modm9"]
    data = {lb: load_dir(D, lb) for lb in labs}
    base = data["py-modb0"]
    p("## 4. Módulo nativo visto desde Python (bench_python.py)", "")
    p("| escenario | " + " | ".join(f"{lb[3:]} mediana ms (IQR %)" for lb in labs) + " | 5f49217/ref | a6/ref | a6/5f49217 | m9/5f49217 |",
      "|---|" + "---:|" * (len(labs) + 4))
    for n in base:
        cells = []
        for lb in labs:
            e = data[lb].get(n)
            cells.append(f"{ft(e['st']['med'])} ({e['st']['rel']*100:.0f})" if e else "—")
        def sp(x, y):
            return data[x][n]["st"]["med"] / data[y][n]["st"]["med"] if n in data[x] and n in data[y] else None
        p(f"| `{n}` | " + " | ".join(cells) + f" | {fx(sp('py-modb0','py-modh0'))} | {fx(sp('py-modb0','py-moda6'))} | {fx(sp('py-modh0','py-moda6'))} | {fx(sp('py-modh0','py-modm9'))} |")
    for lb in labs[1:]:
        n, dif, err = checks(base, data[lb])
        s = [base[k]["st"]["med"] / data[lb][k]["st"]["med"] for k in data[lb] if k in base]
        RES.setdefault("mod", {})[lb] = (gm(s), len(dif), len(err), n)
    s2 = [data["py-modh0"][k]["st"]["med"] / data["py-moda6"][k]["st"]["med"] for k in base]
    RES["mod_a6_h0"] = gm(s2)
    p("", "| variante | media geom. vs módulo original | checksums distintos | errores |", "|---|---:|---:|---:|")
    for lb, (g, d, e, n) in RES["mod"].items():
        p(f"| `{lb}` | {fx(g)} | {d}/{n} | {e} |")
    p(f"| `py-moda6` frente a `py-modh0` | {fx(RES['mod_a6_h0'])} | | |", "")

Path(sys.argv[1] if len(sys.argv) > 1 else "tablas.md").write_text("\n".join(OUT), encoding="utf-8")
Path(Path(__file__).parent / ("res%s.json" % AC)).write_text(json.dumps(RES, indent=1, default=str), encoding="utf-8")
print("ok", len(OUT))

# ---------------------------------------------------------------- capa Python
D = R / "capa_python"
CAPA = [("q00", "q01", "8fcc151", "métodos nativos ligados en __slots__", r"^alcp\.(fp|f9|gf256)\."),
        ("q01", "q02", "abc2f4d", "camino rápido int en + - *", r"^alcp\.(fp\.mixed_int|gf256\.add_int|fp\.sum_elems)"),
        ("q02", "q03", "27873ea", "analizador: x^k como monomio", r"^alcp\.poly\.(parse_long|from_str)"),
        ("q03", "q04", "29d5628", "caché de cadenas analizadas", r"^alcp\.poly\.from_str"),
        ("q04", "q05", "0ab679a", "_pretty: cada coeficiente formateado una vez", r"^alcp\.poly\.str_"),
        ("q05", "q06", "a3399d9", "caché de textos de elementos", r"^alcp\.(elem\.str_gf256|poly\.str_f9)"),
        ("q06", "q07", "cab2c2a", "pickle de elementos más directo", r"^alcp\.pickle\.elems"),
        ("q07", "q08", "a833bcb", "caminos rápidos int/elemento en conversión y f(v)", r"^alcp\.poly\.(eval|mul_elem|add_int|rmul_int)"),
        ("q08", "q09", "e683a88", "listas de coeficientes con caminos rápidos", r"^alcp\.poly\.from_(list|elems)"),
        ("q09", "q10", "06d922a", "caminos rápidos // % divmod", r"^alcp\.poly\.(floordiv|mod)"),
        ("q10", "q11", "0438b74", "resultado en línea en + - * / (DUDOSA)", r"^alcp\.(fp|f9|gf256)\."),
        ("q11", "q12", "c70bed5", "_pretty: potencias en línea", r"^alcp\.poly\.str_"),
        ("q03", "q13", "ff9ad21 (exp/)", "EXPERIMENTO tokenizador regex", r"^alcp\.poly\.(parse_long|from_str)")]
if D.is_dir() and (D / "q12.json").is_file():
    labs = sorted({x for t in CAPA for x in t[:2]})
    data = {lb: load_dir(D, lb) for lb in labs}
    names = [n for n in data["q00"] if n.startswith("alcp.")]
    p("## 5. Capa Python (bench_alcp.py, mismo .pyd de 5f49217)", "")
    p("| commit | técnica | n afect. | media geom. paso (afectados) | rango | media geom. resto alcp.* |",
      "|---|---|---:|---:|---:|---:|")
    det = []
    for prev, cur, com, desc, rx in CAPA:
        A, B = data[prev], data[cur]
        af = [n for n in names if re.search(rx, n) and n in B]
        ct = [n for n in names if not re.search(rx, n) and n in B]
        sa = [A[n]["st"]["med"] / B[n]["st"]["med"] for n in af]
        sc = [A[n]["st"]["med"] / B[n]["st"]["med"] for n in ct]
        rng = f"{min(sa):.2f}-{max(sa):.2f}" if sa else "—"
        RES.setdefault("capa", {})[cur] = {"commit": com, "gm": gm(sa), "gm_ctl": gm(sc), "rng": rng, "n": len(af)}
        p(f"| {com} | {desc} | {len(af)} | {fx(gm(sa))} | {rng} | {fx(gm(sc))} |")
        rows, _ = tabla_paso(A, B, af, prev, cur)
        det.append((com, desc, rows))
    p("")
    # Checksums: alcp.* debe coincidir entre variantes
    malos = [(lb, n) for lb in labs for n in data[lb] if n in data["q00"] and data[lb][n]["ck"] != data["q00"][n]["ck"]]
    errs = [(lb, n) for lb in labs for n in data[lb] if data[lb][n]["err"]]
    RES["capa_ck"] = (malos, errs)
    p(f"Checksums: {len(malos)} discrepancias frente a 5f49217 y {len(errs)} errores en {len(labs)} variantes.", "")
    p("### Antes y después (5f49217 frente a c70bed5): tiempo y cociente paquete/nativo", "")
    p("| escenario | 5f49217 ms | c70bed5 ms | speedup | cociente alcp/nat antes | después |", "|---|---:|---:|---:|---:|---:|")
    A, B = data["q00"], data["q12"]
    tot = []
    for n in names:
        if n not in B:
            continue
        nn = "nat." + n[5:]
        s = A[n]["st"]["med"] / B[n]["st"]["med"]
        tot.append(s)
        ra = A[n]["st"]["med"] / A[nn]["st"]["med"] if nn in A else None
        rb_ = B[n]["st"]["med"] / B[nn]["st"]["med"] if nn in B else None
        p(f"| `{n}` | {ft(A[n]['st']['med'])} | {ft(B[n]['st']['med'])} | {fx(s)} | {fx(ra)} | {fx(rb_)} |")
    RES["capa_total"] = gm(tot)
    p("", f"Media geométrica del speedup de la capa (41 escenarios alcp.*): {fx(gm(tot))}", "")
    for com, desc, rows in det:
        p(f"<details><summary>{com} — {desc}</summary>", "", *rows, "", "</details>", "")

Path(sys.argv[1] if len(sys.argv) > 1 else "tablas.md").write_text("\n".join(OUT), encoding="utf-8")
Path(Path(__file__).parent / ("res%s.json" % AC)).write_text(json.dumps(RES, indent=1, default=str), encoding="utf-8")
print("ok", len(OUT))

# ---------------------------------------------------------------- focalizadas
FOC = [("barrett", "a4", "a5", "2297387 Barrett sin desplazamientos variables"),
       ("sumaresta", "a2", "a3", "dcb905c suma/resta de Poly"),
       ("karatsuba_vistas", "m2", "m3", "9bfc1e7 Karatsuba sobre vistas"),
       ("ntt_par", "m8", "m9", "b0a9fe7 NTT paralela"),
       ("umbrales_par", "m5", "m6", "3660189 umbrales de despacho Auto"),
       ("umbrales_par", "m7", "m8", "e8e07e6 Karatsuba paralelo"),
       ("umbrales_par", "m7", "m8t1", "e8e07e6 con ALCP_THREADS=1 (coste de la infraestructura paralela)")]
if (R / "focal").is_dir():
    p("## 6. Mediciones focalizadas de las técnicas dudosas (6 rondas x >= 5 repeticiones, con corriente)", "")
    for d, a, b, desc in FOC:
        D = R / "focal" / d
        if not (D / f"{b}.json").is_file():
            continue
        A, B = load_dir(D, a), load_dir(D, b)
        nm = list(A)
        rows, sp = tabla_paso(A, B, nm, a, b)
        s = [x for _, x in sp]
        _, dif, err = checks(b0, B)
        RES.setdefault("focal", {})[f"{d}/{b}"] = {"gm": gm(s), "min": min(s), "max": max(s), "n": len(s),
                                                  "dif": dif, "err": err}
        p(f"### {desc} ({a} -> {b})", "", f"Media geométrica {fx(gm(s))} en {len(s)} escenarios (rango {min(s):.2f}-{max(s):.2f}); "
          f"checksums distintos de la referencia: {len(dif)}; errores: {len(err)}.", "", *rows, "")

Path(sys.argv[1] if len(sys.argv) > 1 else "tablas.md").write_text("\n".join(OUT), encoding="utf-8")
Path(Path(__file__).parent / ("res%s.json" % AC)).write_text(json.dumps(RES, indent=1, default=str), encoding="utf-8")
print("ok", len(OUT))
