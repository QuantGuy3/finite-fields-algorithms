"""une_rondas.py — Une varios JSON de una misma variante (una ronda cada uno) en uno solo.

Uso: python une_rondas.py SALIDA.json ETIQUETA RONDA1.json RONDA2.json ...

Se usa con bench_python.py y bench_alcp.py, que no tienen modo intercalado propio
para dos módulos: las rondas se alternan a mano (A B A B ...) y después se juntan
aquí.  Para cada escenario se concatenan todas las muestras, se recalculan la
mediana y el mínimo, y se guarda la mediana de cada ronda en round_medians_s
(mismo formato que run_bench.py interleave).  Si el checksum cambia entre rondas,
se marca como VARIABLE (igual que run_bench.py), porque indicaría no determinismo.
"""
import json
import statistics
import sys


def main() -> int:
    salida, etiqueta, *entradas = sys.argv[1:]
    datos = [json.load(open(f, encoding="utf-8")) for f in entradas]
    base = datos[0]
    escenarios = []
    for i, esc in enumerate(base["scenarios"]):
        muestras, medianas, checks = [], [], set()
        for d in datos:
            e = d["scenarios"][i]
            assert e["name"] == esc["name"], "las rondas no tienen los mismos escenarios"
            muestras += e["samples_s"]
            medianas.append(e["median_s"])
            checks.add(e["checksum"])
        nuevo = dict(esc)
        nuevo["samples_s"] = muestras
        nuevo["median_s"] = statistics.median(muestras)
        nuevo["min_s"] = min(muestras)
        nuevo["round_medians_s"] = medianas
        if len(checks) > 1:
            nuevo["checksum"] = "VARIABLE:" + ",".join(sorted(str(c) for c in checks))
            nuevo["error"] = "checksum distinto entre rondas"
        escenarios.append(nuevo)
    res = dict(base)
    res["label"] = etiqueta
    res["mode"] = "interleave"
    res["runs"] = [{"round": k, "source": f.replace("\\", "/").split("/")[-1]} for k, f in enumerate(entradas)]
    res["scenarios"] = escenarios
    with open(salida, "w", encoding="utf-8") as fh:
        json.dump(res, fh, ensure_ascii=False, indent=1)
    print(f"{salida}: {len(escenarios)} escenarios, {len(entradas)} rondas")
    return 0


if __name__ == "__main__":
    sys.exit(main())
