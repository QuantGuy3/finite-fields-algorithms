# Elimina de un JSON de interleave las rondas que se solapan con un intervalo de
# Modern Standby (tapa cerrada) y escribe la versión filtrada en <dir>/limpio/.
# Uso: python limpia.py DIR [DIR ...]
import json, subprocess, sys
from datetime import datetime, timedelta
from pathlib import Path


def intervalos():
    ps = ("Get-WinEvent -FilterHashtable @{LogName='System';ProviderName='Microsoft-Windows-Kernel-Power';"
          "Id=506,507;StartTime=(Get-Date).AddHours(-12)} | Sort TimeCreated | "
          "% { '{0:yyyy-MM-ddTHH:mm:ss} {1}' -f $_.TimeCreated,$_.Id }")
    out = subprocess.run(["powershell", "-NoProfile", "-Command", ps], capture_output=True, text=True).stdout
    iv, ini = [], None
    for ln in out.split():
        pass
    toks = out.split()
    for t, i in zip(toks[0::2], toks[1::2]):
        ts = datetime.fromisoformat(t)
        if i == "506":
            ini = ts
        elif i == "507":
            iv.append((ini or ts, ts))
            ini = None
    if ini:
        iv.append((ini, datetime.now()))
    return iv


IV = intervalos()
MARGEN = timedelta(seconds=1)


def contaminada(run):
    fin = datetime.fromisoformat(run["finished"]).replace(tzinfo=None)
    ini = fin - timedelta(seconds=run["total_s"])
    return any(a - MARGEN <= fin and b + MARGEN >= ini for a, b in IV)


for d in map(Path, sys.argv[1:]):
    (d / "limpio").mkdir(exist_ok=True)
    for f in sorted(d.glob("*.json")):
        js = json.loads(f.read_text(encoding="utf-8"))
        if js.get("mode") != "interleave" or not js.get("runs"):
            continue
        malas = {r["round"] for r in js["runs"] if contaminada(r)}
        rondas = sorted({r["round"] for r in js["runs"]})
        buenas = [r for r in rondas if r not in malas]
        for s in js["scenarios"]:
            k = len(s["samples_s"]) // max(1, len(rondas))  # muestras por ronda (reps fijas, min_time 0)
            if k * len(rondas) != len(s["samples_s"]):
                continue  # número variable de muestras: no se puede partir por ronda
            s["samples_s"] = [x for i, x in enumerate(s["samples_s"]) if (i // k) in buenas] or s["samples_s"]
            s["round_medians_s"] = [x for i, x in enumerate(s["round_medians_s"]) if i in buenas] or s["round_medians_s"]
        js["rondas_descartadas_standby"] = sorted(malas)
        (d / "limpio" / f.name).write_text(json.dumps(js, indent=1, ensure_ascii=False), encoding="utf-8")
        print(f"{d.name}/{f.name}: rondas {rondas}, descartadas por suspensión {sorted(malas)}")
