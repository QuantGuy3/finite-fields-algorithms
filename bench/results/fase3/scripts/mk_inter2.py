# Crea inter2.py: inter.py + detección de Modern Standby y repetición de la ejecución afectada.
from pathlib import Path
d = Path(__file__).parent
s = (d / "inter.py").read_text(encoding="utf-8")
old_run = '''    import subprocess
    t0 = time.perf_counter()
    cp = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
'''
new_run = '''    import subprocess
    for intento in range(3):
        ini = time.time()
        t0 = time.perf_counter()
        cp = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
        fin = time.time()
        susp = hubo_suspension(ini, fin)
        with open(S / "med-logs" / "ejecuciones.log", "a", encoding="utf-8") as fl:
            fl.write(f"{label}\\t{time.strftime('%H:%M:%S', time.localtime(ini))}\\t"
                     f"{time.strftime('%H:%M:%S', time.localtime(fin))}\\t{fin - ini:.1f}\\t"
                     f"{'SUSPENSION' if susp else 'ok'}\\t{filt or ''}\\n")
        if not susp:
            break
        print(f"  [{label}] suspension (Modern Standby) durante la ejecucion: se repite", flush=True)
'''
old_def = "def _label_dir(label: str) -> Path:"
new_def = '''def hubo_suspension(ini: float, fin: float) -> bool:
    """True si el registro del sistema muestra Modern Standby (Kernel-Power 506/507) entre ini y fin."""
    import subprocess
    ps = ("$i=[DateTimeOffset]::FromUnixTimeSeconds(%d).LocalDateTime;"
          "$f=[DateTimeOffset]::FromUnixTimeSeconds(%d).LocalDateTime;"
          "@(Get-WinEvent -FilterHashtable @{LogName='System';ProviderName='Microsoft-Windows-Kernel-Power';"
          "Id=506,507;StartTime=$i;EndTime=$f} -ErrorAction SilentlyContinue).Count") % (int(ini) - 2, int(fin) + 2)
    try:
        out = subprocess.run(["powershell", "-NoProfile", "-Command", ps], capture_output=True, text=True, timeout=60)
        return int((out.stdout.strip() or "0").splitlines()[-1]) > 0
    except Exception:
        return False


''' + old_def
assert old_run in s and old_def in s
s = s.replace(old_run, new_run).replace(old_def, new_def, 1)
(d / "inter2.py").write_text(s, encoding="utf-8")
print("ok")
