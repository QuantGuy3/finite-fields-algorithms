#!/usr/bin/env python3
"""run_bench.py — Orquestador del harness de benchmarks de alcp.

Compila ``bench/bench_suite.cpp`` contra un directorio de cabeceras (la línea
base original o la versión mejorada), ejecuta la suite guardando metadatos de
la máquina, intercala variantes para reducir la deriva térmica y compara
resultados (mediana, mínimo, IQR, speedup y verificación de checksums).

Subcomandos::

    build      --include DIR --label NOMBRE [--cxxflags "..."] [--ldflags "..."] [--cxx g++]
    run        --label NOMBRE [--filter RE] [--reps N] [--min-time S] [--warmup N] [--out RUTA]
    interleave --labels A B C [--rounds R] [--filter RE] [--reps N] [--min-time S] [--per-scenario]
    compare    --base A --labels B C [--md salida.md] [--html salida.html]
    module     --source bindings.cpp --include DIR --label NOMBRE [--cxxflags "..."]

Las etiquetas (``--label``) identifican cada variante: el ejecutable va a
``build/bench/<label>/bench_suite.exe`` y los resultados a
``bench/results/<label>.json``.  En ``compare`` se puede pasar, en lugar de una
etiqueta, la ruta a un JSON (también los producidos por ``bench_python.py``).

Ejecutar SIEMPRE con el intérprete del entorno virtual del proyecto.
"""
from __future__ import annotations

import argparse
import ctypes
import datetime as _dt
import hashlib
import html
import json
import math
import os
import platform
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import sysconfig
import time
from pathlib import Path

BENCH_DIR = Path(__file__).resolve().parent
REPO = BENCH_DIR.parent
BUILD_ROOT = REPO / "build" / "bench"
RESULTS_DIR = BENCH_DIR / "results"
SUITE_SRC = BENCH_DIR / "bench_suite.cpp"
IS_WIN = os.name == "nt"
EXE_NAME = "bench_suite.exe" if IS_WIN else "bench_suite"
RESULT_FORMAT = "alcp-bench-result"

# Opciones cuyo valor suele empezar por "-" (p. ej. --cxxflags -O3): se
# reescriben como --opcion=valor antes de argparse para que no se confundan con
# otra opción.
_FLAG_VALUE_OPTS = ("--cxxflags", "--ldflags")


# ---------------------------------------------------------------------------
# Utilidades generales
# ---------------------------------------------------------------------------
def _now_iso() -> str:
    return _dt.datetime.now().astimezone().isoformat(timespec="seconds")


def _run_text(cmd: list[str], timeout: float = 30, encoding: str | None = None) -> str | None:
    """Ejecuta un comando de solo lectura y devuelve su salida, o None si falla."""
    try:
        cp = subprocess.run(cmd, capture_output=True, timeout=timeout, check=False)
    except (OSError, subprocess.SubprocessError):
        return None
    if cp.returncode != 0:
        return None
    raw = cp.stdout
    for enc in ([encoding] if encoding else []) + ["utf-8", "oem" if IS_WIN else "latin-1", "latin-1"]:
        try:
            return raw.decode(enc).strip()
        except (LookupError, UnicodeDecodeError):
            continue
    return raw.decode("utf-8", errors="replace").strip()


def _label_dir(label: str) -> Path:
    if not re.fullmatch(r"[A-Za-z0-9._+-]+", label):
        sys.exit(f"etiqueta inválida '{label}': use solo letras, dígitos y . _ + -")
    return BUILD_ROOT / label


def _sha256_files(paths: list[Path]) -> str:
    h = hashlib.sha256()
    for p in sorted(paths):
        h.update(p.name.encode())
        h.update(p.read_bytes())
    return h.hexdigest()


def _home_variants() -> list[str]:
    home = str(Path.home())
    return sorted({home, home.replace("\\", "/"), home.replace("/", "\\")}, key=len, reverse=True)


def scrub_paths(obj):
    """Sustituye el directorio personal por "~" en todas las cadenas de `obj`.

    Los metadatos contienen rutas absolutas (include, compilador, Python...) que
    incluirían el nombre de usuario; se anonimizan antes de escribir resultados
    que pueden acabar versionados (minimización de datos personales).
    """
    if isinstance(obj, str):
        out = obj
        flags = re.IGNORECASE if IS_WIN else 0
        for h in _home_variants():
            out = re.sub(re.escape(h), "~", out, flags=flags)
        # El nombre de usuario también aparece codificado en algunas rutas
        # temporales (p. ej. "C--Users-<usuario>-Documents"): se sustituye cuando
        # va delimitado por separadores de ruta.
        user = Path.home().name
        if len(user) >= 3:
            out = re.sub(r"(?<=[\\/\-_.])" + re.escape(user) + r"(?=[\\/\-_.]|$)", "<usuario>", out, flags=flags)
        return out
    if isinstance(obj, list):
        return [scrub_paths(x) for x in obj]
    if isinstance(obj, dict):
        return {k: scrub_paths(v) for k, v in obj.items()}
    return obj


def utf8_stdout() -> None:
    """Fuerza UTF-8 en stdout/stderr cuando no son una consola (tuberías, archivos)."""
    for stream in (sys.stdout, sys.stderr):
        try:
            if not stream.isatty():
                stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError, OSError):
            pass


def _split_flags(s: str) -> list[str]:
    # posix=True: comillas al estilo shell; en Windows use / en las rutas.
    return shlex.split(s or "")


# ---------------------------------------------------------------------------
# Metadatos de la máquina y del código
# ---------------------------------------------------------------------------
def _cpu_name() -> str:
    if IS_WIN:
        try:
            import winreg

            with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, r"HARDWARE\DESCRIPTION\System\CentralProcessor\0") as k:
                return str(winreg.QueryValueEx(k, "ProcessorNameString")[0]).strip()
        except OSError:
            pass
    else:
        try:
            for line in Path("/proc/cpuinfo").read_text().splitlines():
                if line.lower().startswith("model name"):
                    return line.split(":", 1)[1].strip()
        except OSError:
            pass
        if sys.platform == "darwin":
            out = _run_text(["sysctl", "-n", "machdep.cpu.brand_string"])
            if out:
                return out
    return platform.processor() or "desconocida"


def _physical_cores() -> int | None:
    try:
        import psutil  # opcional

        return psutil.cpu_count(logical=False)
    except ImportError:
        pass
    if IS_WIN:
        out = _run_text(["powershell", "-NoProfile", "-NonInteractive", "-Command",
                         "(Get-CimInstance Win32_Processor | Measure-Object -Property NumberOfCores -Sum).Sum"],
                        timeout=60)
        if out and out.strip().isdigit():
            return int(out.strip())
        return None
    try:
        cores = set()
        phys = core = None
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("physical id"):
                phys = line.split(":", 1)[1].strip()
            elif line.startswith("core id"):
                core = line.split(":", 1)[1].strip()
            elif not line.strip():
                if core is not None:
                    cores.add((phys, core))
                phys = core = None
        return len(cores) or None
    except OSError:
        return None


def _power_plan() -> str | None:
    if IS_WIN:
        # Solo consulta (no modifica nada).
        return _run_text(["powercfg", "/getactivescheme"], encoding="oem")
    gov = Path("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor")
    try:
        return "governor: " + gov.read_text().strip()
    except OSError:
        return None


def _battery() -> dict | None:
    if IS_WIN:
        class SYSTEM_POWER_STATUS(ctypes.Structure):  # noqa: N801 (nombre de la API Win32)
            _fields_ = [("ACLineStatus", ctypes.c_ubyte), ("BatteryFlag", ctypes.c_ubyte),
                        ("BatteryLifePercent", ctypes.c_ubyte), ("SystemStatusFlag", ctypes.c_ubyte),
                        ("BatteryLifeTime", ctypes.c_ulong), ("BatteryFullLifeTime", ctypes.c_ulong)]

        st = SYSTEM_POWER_STATUS()
        try:
            if not ctypes.windll.kernel32.GetSystemPowerStatus(ctypes.byref(st)):
                return None
        except (AttributeError, OSError):
            return None
        ac = {0: False, 1: True}.get(st.ACLineStatus)
        pct = None if st.BatteryLifePercent == 255 else int(st.BatteryLifePercent)
        return {"ac_online": ac, "battery_percent": pct, "battery_saver": bool(st.SystemStatusFlag & 1),
                "no_battery": bool(st.BatteryFlag & 128)}
    base = Path("/sys/class/power_supply")
    if not base.is_dir():
        return None
    info: dict = {}
    for d in base.iterdir():
        try:
            typ = (d / "type").read_text().strip()
            if typ == "Mains":
                info["ac_online"] = (d / "online").read_text().strip() == "1"
            elif typ == "Battery":
                info["battery_percent"] = int((d / "capacity").read_text().strip())
        except (OSError, ValueError):
            continue
    return info or None


def _compiler_info(cxx: str) -> dict:
    first = _run_text([cxx, "--version"])
    return {"cxx": cxx, "cxx_path": shutil.which(cxx), "version": first.splitlines()[0] if first else None,
            "dumpfullversion": _run_text([cxx, "-dumpfullversion"])}


def _git_info(path: Path) -> dict:
    """Commit y estado del repositorio que contiene `path` (solo lectura)."""
    base = ["git", "--no-optional-locks", "-C", str(path)]
    info = {
        "commit": _run_text(base + ["rev-parse", "HEAD"]),
        "describe": _run_text(base + ["describe", "--tags", "--always"]),
        "branch": _run_text(base + ["rev-parse", "--abbrev-ref", "HEAD"]),
    }
    st = _run_text(base + ["status", "--porcelain", "--", "."])
    info["dirty"] = None if st is None else bool(st.strip())
    return info


_MACHINE_CACHE: dict | None = None


def collect_machine_metadata() -> dict:
    """Metadatos de la máquina (CPU, núcleos, SO, energía).  Se cachean por proceso."""
    global _MACHINE_CACHE
    if _MACHINE_CACHE is None:
        _MACHINE_CACHE = {
            "cpu": _cpu_name(),
            "logical_cpus": os.cpu_count(),
            "physical_cores": _physical_cores(),
            "os": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version.split()[0],
            "python_executable": sys.executable,
        }
    md = dict(_MACHINE_CACHE)
    # Energía y batería pueden cambiar durante la sesión: se leen cada vez.
    md["power_plan"] = _power_plan()
    md["power"] = _battery()
    md["date"] = _now_iso()
    return md


# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------
def default_ldflags() -> str:
    # En Windows se enlaza todo estáticamente (ejecutable autocontenido, sin
    # depender de las DLL de MSYS2); en Linux, GMP dinámica del sistema.
    return "-static -lgmpxx -lgmp" if IS_WIN else "-lgmpxx -lgmp"


def cmd_build(a: argparse.Namespace) -> int:
    include = Path(a.include).resolve()
    if not (include / "alcp" / "alcp.hpp").is_file():
        sys.exit(f"{include} no contiene alcp/alcp.hpp")
    out_dir = _label_dir(a.label)
    out_dir.mkdir(parents=True, exist_ok=True)
    exe = out_dir / EXE_NAME
    ldflags = a.ldflags if a.ldflags is not None else default_ldflags()
    cmd = [a.cxx, "-std=c++23", "-Wall", "-Wextra", *_split_flags(a.cxxflags), f"-I{include}",
           str(SUITE_SRC), *_split_flags(ldflags), "-o", str(exe)]
    print("$ " + " ".join(shlex.quote(c) for c in cmd), flush=True)
    t0 = time.perf_counter()
    cp = subprocess.run(cmd, check=False)
    dt = time.perf_counter() - t0
    if cp.returncode != 0:
        print(f"ERROR: la compilación falló (código {cp.returncode})", file=sys.stderr)
        return cp.returncode
    headers = sorted((include / "alcp").glob("*.hpp"))
    info = {
        "label": a.label,
        "kind": "cpp",
        "include": str(include),
        "include_git": _git_info(include),
        "include_sha256": _sha256_files(headers),
        "suite_source_sha256": _sha256_files([SUITE_SRC]),
        "compiler": _compiler_info(a.cxx),
        "cxxflags": a.cxxflags,
        "ldflags": ldflags,
        "command": cmd,
        "build_seconds": round(dt, 2),
        "built": _now_iso(),
    }
    (out_dir / "build_info.json").write_text(json.dumps(scrub_paths(info), indent=2, ensure_ascii=False),
                                             encoding="utf-8")
    print(f"-> {exe}  ({dt:.1f} s)")
    return 0


# ---------------------------------------------------------------------------
# Ejecución de la suite
# ---------------------------------------------------------------------------
def _load_build_info(label: str) -> dict:
    p = _label_dir(label) / "build_info.json"
    if not p.is_file():
        sys.exit(f"no existe {p}: compile antes con  run_bench.py build --label {label} ...")
    info = json.loads(p.read_text(encoding="utf-8"))
    if info.get("suite_source_sha256") != _sha256_files([SUITE_SRC]):
        print(f"aviso: bench_suite.cpp ha cambiado desde que se compiló '{label}'; "
              f"recompile para medir la versión actual de la suite", file=sys.stderr)
    return info


def _exe_for(label: str) -> Path:
    exe = _label_dir(label) / EXE_NAME
    if not exe.is_file():
        sys.exit(f"no existe {exe}: compile antes con  run_bench.py build --label {label} ...")
    return exe


def _popen_retry(cmd: list[str]) -> subprocess.Popen:
    # En Windows un antivirus puede bloquear durante un instante un .exe recién
    # creado: se reintenta unas pocas veces antes de rendirse.
    last: OSError | None = None
    for _ in range(6):
        try:
            return subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        except OSError as e:  # PermissionError, WinError 32...
            last = e
            time.sleep(0.5)
    raise last  # type: ignore[misc]


def run_suite(label: str, *, filt: str | None, reps: int | None, min_time: float | None,
              warmup: int | None, echo: bool = True) -> dict:
    """Ejecuta bench_suite de `label` y devuelve el JSON que produce."""
    exe = _exe_for(label)
    tmp = _label_dir(label) / "last_run.json"
    if tmp.exists():
        tmp.unlink()
    cmd = [str(exe), "--json", str(tmp)]
    if filt:
        cmd += ["--filter", filt]
    if reps is not None:
        cmd += ["--reps", str(reps)]
    if min_time is not None:
        cmd += ["--min-time", str(min_time)]
    if warmup is not None:
        cmd += ["--warmup", str(warmup)]
    proc = _popen_retry(cmd)
    assert proc.stdout is not None
    for raw in proc.stdout:
        line = raw.decode("utf-8", errors="replace").rstrip("\r\n")
        if echo:
            print(f"[{label}] {line}", flush=True)
    rc = proc.wait()
    if not tmp.is_file():
        sys.exit(f"[{label}] la suite terminó con código {rc} sin escribir resultados")
    data = json.loads(tmp.read_text(encoding="utf-8"))
    data["_returncode"] = rc
    return data


def _suite_meta(d: dict) -> dict:
    keys = ("suite", "suite_version", "compiler", "cplusplus", "gmp_version", "ndebug")
    return {k: d.get(k) for k in keys}


def _run_record(d: dict, label: str, round_no: int) -> dict:
    return {"label": label, "round": round_no, "finished": _now_iso(), "total_s": d.get("total_s"),
            "reps": d.get("reps"), "min_time_s": d.get("min_time_s"), "warmup": d.get("warmup"),
            "filter": d.get("filter"), "returncode": d.get("_returncode")}


def _new_result(label: str, mode: str, build: dict, suite: dict) -> dict:
    return {"format": RESULT_FORMAT, "format_version": 1, "label": label, "kind": "cpp", "mode": mode,
            "created": _now_iso(), "metadata": {**collect_machine_metadata(), "build": build, "suite": suite},
            "runs": [], "scenarios": []}


def _merge_scenarios(res: dict, d: dict) -> None:
    """Añade las muestras de una ejecución `d` al resultado acumulado `res`."""
    index = {s["name"]: s for s in res["scenarios"]}
    for s in d["scenarios"]:
        tgt = index.get(s["name"])
        if tgt is None:
            tgt = {k: s.get(k) for k in ("name", "group", "op", "desc", "inner")}
            tgt.update({"samples_s": [], "round_medians_s": [], "setup_s": [], "checksums": [], "errors": []})
            res["scenarios"].append(tgt)
            index[s["name"]] = tgt
        tgt["samples_s"].extend(s.get("samples_s") or [])
        if s.get("samples_s"):
            tgt["round_medians_s"].append(statistics.median(s["samples_s"]))
        tgt["setup_s"].append(s.get("setup_s"))
        if s.get("checksum"):
            tgt["checksums"].append(s["checksum"])
        if s.get("error"):
            tgt["errors"].append(s["error"])


def _finalize(res: dict) -> dict:
    """Resume checksums y errores de todas las rondas en campos únicos."""
    for s in res["scenarios"]:
        uniq = sorted(set(s["checksums"]))
        if len(uniq) == 1:
            s["checksum"] = uniq[0]
        elif not uniq:
            s["checksum"] = ""
        else:  # el resultado cambió entre rondas: no es determinista
            s["checksum"] = "VARIABLE:" + "|".join(uniq)
        s["error"] = "; ".join(sorted(set(s["errors"])))
        s["median_s"] = statistics.median(s["samples_s"]) if s["samples_s"] else None
        s["min_s"] = min(s["samples_s"]) if s["samples_s"] else None
    return res


def _write_result(res: dict, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(scrub_paths(res), indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"-> {out}")


def cmd_run(a: argparse.Namespace) -> int:
    build = _load_build_info(a.label)
    d = run_suite(a.label, filt=a.filter, reps=a.reps, min_time=a.min_time, warmup=a.warmup)
    res = _new_result(a.label, "single", build, _suite_meta(d))
    if a.note:
        res["note"] = a.note
    res["runs"].append(_run_record(d, a.label, 0))
    _merge_scenarios(res, d)
    _finalize(res)
    out = Path(a.out) if a.out else RESULTS_DIR / f"{a.label}.json"
    _write_result(res, out)
    return 0 if d["_returncode"] == 0 else 1


def _list_scenarios(label: str, filt: str | None) -> list[str]:
    cmd = [str(_exe_for(label)), "--list"] + (["--filter", filt] if filt else [])
    cp = subprocess.run(cmd, capture_output=True, check=False)
    if cp.returncode != 0:
        sys.exit(f"[{label}] --list falló: {cp.stdout.decode('utf-8', 'replace')}")
    return [ln.split("\t", 1)[0] for ln in cp.stdout.decode("utf-8", "replace").splitlines() if ln.strip()]


def cmd_interleave(a: argparse.Namespace) -> int:
    labels: list[str] = a.labels
    if len(set(labels)) != len(labels):
        sys.exit("las etiquetas deben ser distintas")
    builds = {lb: _load_build_info(lb) for lb in labels}
    results: dict[str, dict] = {}
    rc_total = 0

    def one(label: str, filt: str | None, round_no: int, echo: bool) -> None:
        nonlocal rc_total
        d = run_suite(label, filt=filt, reps=a.reps, min_time=a.min_time, warmup=a.warmup, echo=echo)
        if label not in results:
            results[label] = _new_result(label, "interleave", builds[label], _suite_meta(d))
            if a.note:
                results[label]["note"] = a.note
            results[label]["interleave"] = {"labels": labels, "rounds": a.rounds,
                                            "per_scenario": a.per_scenario, "reps_per_round": a.reps,
                                            "min_time_per_round": a.min_time}
        results[label]["runs"].append(_run_record(d, label, round_no))
        _merge_scenarios(results[label], d)
        rc_total |= d["_returncode"]

    t0 = time.perf_counter()
    if a.per_scenario:
        # Intercalado fino: para cada escenario, A B C A B C ... (más procesos,
        # pero la deriva térmica entre variantes de un mismo escenario es mínima).
        names = _list_scenarios(labels[0], a.filter)
        for lb in labels[1:]:
            common = set(_list_scenarios(lb, a.filter))
            names = [n for n in names if n in common]
        for i, name in enumerate(names, 1):
            print(f"== [{i}/{len(names)}] {name}", flush=True)
            for r in range(a.rounds):
                for lb in labels:
                    one(lb, "^" + re.escape(name) + "$", r, echo=False)
            for lb in labels:
                s = next(x for x in results[lb]["scenarios"] if x["name"] == name)
                med = statistics.median(s["samples_s"]) * 1e3 if s["samples_s"] else float("nan")
                print(f"   {lb:<24} mediana {med:11.3f} ms  ({len(s['samples_s'])} muestras)", flush=True)
    else:
        for r in range(a.rounds):
            for lb in labels:
                print(f"== ronda {r + 1}/{a.rounds}: {lb}", flush=True)
                one(lb, a.filter, r, echo=True)
    print(f"intercalado completado en {time.perf_counter() - t0:.1f} s")
    out_dir = Path(a.out_dir) if a.out_dir else RESULTS_DIR
    for lb in labels:
        _write_result(_finalize(results[lb]), out_dir / f"{lb}.json")
    return 0 if rc_total == 0 else 1


# ---------------------------------------------------------------------------
# compare
# ---------------------------------------------------------------------------
def _load_result(ref: str, results_dir: Path) -> tuple[str, dict]:
    p = Path(ref)
    if not (p.suffix.lower() == ".json" and p.is_file()):
        p = results_dir / f"{ref}.json"
    if not p.is_file():
        sys.exit(f"no se encuentra el resultado '{ref}' ({p})")
    d = json.loads(p.read_text(encoding="utf-8"))
    label = d.get("label") or p.stem
    if "metadata" not in d:  # salida directa de bench_suite --json
        d = {"label": label, "metadata": {}, "scenarios": d.get("scenarios", [])}
    return label, d


def _stats(samples: list[float]) -> dict:
    if not samples:
        return {"n": 0, "median": None, "min": None, "iqr": None}
    if len(samples) >= 2:
        q = statistics.quantiles(samples, n=4, method="inclusive")
        iqr = q[2] - q[0]
    else:
        iqr = 0.0
    return {"n": len(samples), "median": statistics.median(samples), "min": min(samples), "iqr": iqr}


def _geomean(xs: list[float]) -> float | None:
    xs = [x for x in xs if x and x > 0 and math.isfinite(x)]
    if not xs:
        return None
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


def _fmt_t(s: float | None) -> str:
    if s is None:
        return "-"
    if s >= 1:
        return f"{s:.3f} s"
    if s >= 1e-3:
        return f"{s * 1e3:.3f} ms"
    return f"{s * 1e6:.1f} us"


def _fmt_x(x: float | None) -> str:
    return "-" if x is None else f"{x:.3f}x"


def build_comparison(base_ref: str, refs: list[str], results_dir: Path) -> dict:
    base_label, base = _load_result(base_ref, results_dir)
    variants = [_load_result(r, results_dir) for r in refs]
    labels = [lb for lb, _ in variants]
    vmaps = [{s["name"]: s for s in d["scenarios"]} for _, d in variants]
    rows = []
    for s in base["scenarios"]:
        bst = _stats(s.get("samples_s") or [])
        bsum = s.get("checksum") or ""
        berr = s.get("error") or ""
        row = {"name": s["name"], "group": s.get("group") or "", "op": s.get("op") or "",
               "desc": s.get("desc") or "", "base": bst, "base_checksum": bsum, "base_error": berr,
               "variants": []}
        for vm in vmaps:
            v = vm.get(s["name"])
            if v is None:
                row["variants"].append(None)
                continue
            vst = _stats(v.get("samples_s") or [])
            vsum = v.get("checksum") or ""
            problems = []
            if berr:
                problems.append(f"error en base: {berr}")
            if v.get("error"):
                problems.append(f"error: {v['error']}")
            for tag, cs in (("base", bsum), ("variante", vsum)):
                if cs.startswith("INVALIDO") or cs.startswith("VARIABLE"):
                    problems.append(f"checksum {tag} {cs}")
            if bsum != vsum:
                problems.append(f"checksum distinto ({bsum} != {vsum})")
            speed = (bst["median"] / vst["median"]) if (bst["median"] and vst["median"]) else None
            row["variants"].append({"stats": vst, "checksum": vsum, "speedup": speed,
                                    "error": "; ".join(problems)})
        rows.append(row)

    # Medias geométricas por grupo, por operación y global (se excluyen los
    # escenarios con error de corrección: su tiempo no es comparable).
    def geo(key: str | None) -> dict:
        out: dict = {}
        for j, lb in enumerate(labels):
            buckets: dict[str, list[float]] = {}
            for r in rows:
                v = r["variants"][j]
                if v is None or v["error"] or v["speedup"] is None:
                    continue
                k = r[key] if key else "global"
                buckets.setdefault(k, []).append(v["speedup"])
            for k, xs in buckets.items():
                out.setdefault(k, {})[lb] = {"geomean": _geomean(xs), "n": len(xs)}
        return out

    n_err = sum(1 for r in rows for v in r["variants"] if v is not None and v["error"])
    # Aviso si las variantes se compilaron con versiones distintas de la suite.
    shas = {lb: ((d.get("metadata") or {}).get("build") or {}).get("suite_source_sha256")
            for lb, d in [(base_label, base)] + variants}
    if len({v for v in shas.values() if v}) > 1:
        print("aviso: las variantes se compilaron con versiones distintas de bench_suite.cpp: "
              + ", ".join(f"{k}={str(v)[:12]}" for k, v in shas.items()), file=sys.stderr)
    missing = [(lb, sum(1 for r in rows if r["variants"][j] is None)) for j, lb in enumerate(labels)]
    return {"base_label": base_label, "base": base, "labels": labels, "variants": [d for _, d in variants],
            "rows": rows, "geo_group": geo("group"), "geo_op": geo("op"), "geo_global": geo(None),
            "n_errors": n_err, "missing": missing}


def _enable_ansi() -> bool:
    if not sys.stdout.isatty():
        return False
    if IS_WIN:
        try:  # activa las secuencias VT en la consola de este proceso
            k32 = ctypes.windll.kernel32
            h = k32.GetStdHandle(-11)
            mode = ctypes.c_ulong()
            if k32.GetConsoleMode(h, ctypes.byref(mode)):
                k32.SetConsoleMode(h, mode.value | 0x0004)
        except (AttributeError, OSError):
            return False
    return True


def _print_comparison(c: dict) -> None:
    ansi = _enable_ansi()
    red = (lambda s: f"\033[1;31m{s}\033[0m") if ansi else (lambda s: s)
    green = (lambda s: f"\033[32m{s}\033[0m") if ansi else (lambda s: s)
    labels = c["labels"]
    w = max([len(r["name"]) for r in c["rows"]] + [10])
    head = f"{'escenario':<{w}}  {'base ' + c['base_label']:>18}"
    for lb in labels:
        head += f"  {lb[:18]:>18}  {'speedup':>8}"
    print(head)
    print("-" * len(head))
    for r in c["rows"]:
        line = f"{r['name']:<{w}}  {_fmt_t(r['base']['median']):>18}"
        errs = []
        for lb, v in zip(labels, r["variants"]):
            if v is None:
                line += f"  {'(falta)':>18}  {'':>8}"
                continue
            sp = _fmt_x(v["speedup"])
            if v["speedup"] and v["speedup"] >= 1.05:
                sp = green(f"{sp:>8}")
            else:
                sp = f"{sp:>8}"
            line += f"  {_fmt_t(v['stats']['median']):>18}  {sp}"
            if v["error"]:
                errs.append(f"{lb}: {v['error']}")
        print(line)
        for e in errs:
            print("    " + red("ERROR de corrección -> " + e))
    print()
    for title, geo in (("por grupo", c["geo_group"]), ("por operación", c["geo_op"]), ("global", c["geo_global"])):
        print(f"Media geométrica del speedup ({title}):")
        for k in sorted(geo):
            parts = [f"{lb}={_fmt_x(geo[k].get(lb, {}).get('geomean'))} (n={geo[k].get(lb, {}).get('n', 0)})"
                     for lb in labels]
            print(f"  {k:<28} " + "  ".join(parts))
    if c["n_errors"]:
        print(red(f"\n{c['n_errors']} ERRORES de corrección (checksums distintos o errores): revíselos."))
    else:
        print(green("\nChecksums: todos coinciden con la base."))
    for lb, n in c["missing"]:
        if n:
            print(f"aviso: a '{lb}' le faltan {n} escenarios de la base")


def _meta_lines(label: str, d: dict) -> list[str]:
    md = d.get("metadata") or {}
    b = md.get("build") or {}
    g = b.get("include_git") or {}
    suite = md.get("suite") or {}
    pw = md.get("power") or {}
    ac = pw.get("ac_online")
    out = [f"**{label}**: modo {d.get('mode', '?')}, {md.get('date', '?')}"
           + (f" (nota: {d['note']})" if d.get("note") else ""),
           f"CPU {md.get('cpu', '?')} ({md.get('physical_cores', '?')} núcleos / {md.get('logical_cpus', '?')} hilos), "
           f"{md.get('os', '?')}",
           f"compilador {suite.get('compiler') or (b.get('compiler') or {}).get('version', '?')}, "
           f"GMP {suite.get('gmp_version', md.get('gmp_version', '?'))}, flags `{b.get('cxxflags', '?')}`",
           f"include `{b.get('include', '?')}` @ {g.get('describe') or g.get('commit') or '?'}"
           f"{' (con cambios sin commit)' if g.get('dirty') else ''}, sha256 {str(b.get('include_sha256', '?'))[:12]}",
           f"energía: {md.get('power_plan') or '?'}; corriente: "
           f"{'sí' if ac else ('no (batería)' if ac is False else '?')}"]
    return out


def _render_md(c: dict) -> str:
    labels = c["labels"]
    L = ["# Comparativa de benchmarks alcp", "",
         f"Base: **{c['base_label']}**. Speedup = mediana(base) / mediana(variante); > 1 significa más rápido.",
         "IQR = rango intercuartílico de las muestras. Las medias geométricas excluyen escenarios con error.", ""]
    for lb, d in [(c["base_label"], c["base"])] + list(zip(labels, c["variants"])):
        L += ["- " + "  \n  ".join(_meta_lines(lb, d))]
    L.append("")
    if c["n_errors"]:
        L += [f'<span style="color:red">**ERROR: {c["n_errors"]} escenarios con checksums distintos o errores.**</span>', ""]
    else:
        L += ["Checksums: todos los escenarios coinciden con la base.", ""]
    L += ["## Media geométrica del speedup", "", "| ámbito | " + " | ".join(labels) + " |",
          "|---|" + "---|" * len(labels)]
    for title, geo in (("global", c["geo_global"]), ("grupo", c["geo_group"]), ("operación", c["geo_op"])):
        for k in sorted(geo):
            name = "**global**" if title == "global" else f"{title} `{k}`"
            L.append(f"| {name} | " + " | ".join(
                f"{_fmt_x(geo[k].get(lb, {}).get('geomean'))} (n={geo[k].get(lb, {}).get('n', 0)})" for lb in labels) + " |")
    L += ["", "## Detalle por escenario", ""]
    hdr = "| escenario | base mediana | base IQR |"
    sep = "|---|---:|---:|"
    for lb in labels:
        hdr += f" {lb} mediana | {lb} mín | {lb} IQR | speedup |"
        sep += "---:|---:|---:|---:|"
    L += [hdr + " checksum |", sep + "---|"]
    for r in c["rows"]:
        line = f"| `{r['name']}` | {_fmt_t(r['base']['median'])} | {_fmt_t(r['base']['iqr'])} |"
        errs = []
        for lb, v in zip(labels, r["variants"]):
            if v is None:
                line += " - | - | - | - |"
                continue
            st = v["stats"]
            line += f" {_fmt_t(st['median'])} | {_fmt_t(st['min'])} | {_fmt_t(st['iqr'])} | {_fmt_x(v['speedup'])} |"
            if v["error"]:
                errs.append(f"{lb}: {v['error']}")
        if errs:
            line += ' <span style="color:red">**ERROR**</span> ' + html.escape("; ".join(errs)) + " |"
        else:
            line += f" ok `{r['base_checksum'][:16]}` |"
        L.append(line)
    return "\n".join(L) + "\n"


def _render_html(c: dict) -> str:
    labels = c["labels"]
    e = html.escape

    def sp_cell(x: float | None) -> str:
        if x is None:
            return "<td>-</td>"
        cls = "up" if x >= 1.05 else ("down" if x <= 0.95 else "")
        return f'<td class="num {cls}">{x:.3f}&times;</td>'

    meta = "".join("<li>" + "<br>".join(e(s).replace("**", "") for s in _meta_lines(lb, d)) + "</li>"
                   for lb, d in [(c["base_label"], c["base"])] + list(zip(labels, c["variants"])))
    geo_rows = []
    for title, geo in (("global", c["geo_global"]), ("grupo", c["geo_group"]), ("operación", c["geo_op"])):
        for k in sorted(geo):
            name = "<b>global</b>" if title == "global" else f"{e(title)} <code>{e(k)}</code>"
            cells = "".join(sp_cell(geo[k].get(lb, {}).get("geomean")) for lb in labels)
            geo_rows.append(f"<tr><td>{name}</td>{cells}</tr>")
    head = "<th>escenario</th><th>base mediana</th><th>base IQR</th>" + "".join(
        f"<th>{e(lb)} mediana</th><th>{e(lb)} mín</th><th>{e(lb)} IQR</th><th>speedup</th>" for lb in labels
    ) + "<th>checksum</th>"
    body = []
    for r in c["rows"]:
        cells = [f'<td><code title="{e(r["desc"])}">{e(r["name"])}</code></td>',
                 f'<td class="num">{_fmt_t(r["base"]["median"])}</td>', f'<td class="num">{_fmt_t(r["base"]["iqr"])}</td>']
        errs = []
        for lb, v in zip(labels, r["variants"]):
            if v is None:
                cells += ["<td>-</td>"] * 4
                continue
            st = v["stats"]
            cells += [f'<td class="num">{_fmt_t(st["median"])}</td>', f'<td class="num">{_fmt_t(st["min"])}</td>',
                      f'<td class="num">{_fmt_t(st["iqr"])}</td>', sp_cell(v["speedup"])]
            if v["error"]:
                errs.append(f"{lb}: {v['error']}")
        if errs:
            cells.append(f'<td class="err">ERROR de corrección: {e("; ".join(errs))}</td>')
            body.append(f'<tr class="errrow">{"".join(cells)}</tr>')
        else:
            cells.append(f'<td><code>{e(r["base_checksum"][:16])}</code></td>')
            body.append(f"<tr>{''.join(cells)}</tr>")
    status = (f'<p class="err big">ERROR: {c["n_errors"]} escenarios con checksums distintos o errores.</p>'
              if c["n_errors"] else '<p class="ok big">Checksums: todos los escenarios coinciden con la base.</p>')
    return f"""<!doctype html>
<html lang="es"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Comparativa alcp</title>
<style>
:root {{ --bg:#ffffff; --fg:#1d232b; --muted:#5b6672; --line:#d9dee4; --head:#f2f4f7;
        --up:#1a7f37; --down:#b42318; --errbg:#fde8e7; }}
@media (prefers-color-scheme: dark) {{
  :root {{ --bg:#14181d; --fg:#e6e9ee; --muted:#9aa5b1; --line:#2c333b; --head:#1c2228;
          --up:#4ac26b; --down:#ff7b72; --errbg:#3a1d1d; }} }}
body {{ background:var(--bg); color:var(--fg); font:14px/1.45 system-ui, sans-serif; margin:0 16px 40px; }}
h1 {{ font-size:22px; margin:20px 0 6px; }} h2 {{ font-size:17px; margin:26px 0 8px; }}
p.note, ul.meta {{ color:var(--muted); }} ul.meta li {{ margin-bottom:8px; }}
.wrap {{ overflow-x:auto; }}
table {{ border-collapse:collapse; font-variant-numeric:tabular-nums; }}
th, td {{ border-bottom:1px solid var(--line); padding:4px 10px; text-align:left; white-space:nowrap; }}
th {{ background:var(--head); position:sticky; top:0; }}
td.num {{ text-align:right; }} .up {{ color:var(--up); font-weight:600; }} .down {{ color:var(--down); font-weight:600; }}
tr.errrow td {{ background:var(--errbg); }} .err {{ color:var(--down); font-weight:700; }}
.ok {{ color:var(--up); font-weight:600; }} .big {{ font-size:15px; }}
</style></head><body>
<h1>Comparativa de benchmarks alcp</h1>
<p class="note">Base: <b>{e(c['base_label'])}</b>. Speedup = mediana(base) / mediana(variante); &gt; 1 es más rápido.
IQR = rango intercuartílico. Las medias geométricas excluyen los escenarios con error.</p>
{status}
<ul class="meta">{meta}</ul>
<h2>Media geométrica del speedup</h2>
<div class="wrap"><table><tr><th>ámbito</th>{''.join(f'<th>{e(lb)}</th>' for lb in labels)}</tr>
{''.join(geo_rows)}</table></div>
<h2>Detalle por escenario</h2>
<div class="wrap"><table><tr>{head}</tr>
{''.join(body)}</table></div>
</body></html>
"""


def cmd_compare(a: argparse.Namespace) -> int:
    rdir = Path(a.results_dir) if a.results_dir else RESULTS_DIR
    c = build_comparison(a.base, a.labels, rdir)
    _print_comparison(c)
    if a.md:
        Path(a.md).parent.mkdir(parents=True, exist_ok=True)
        Path(a.md).write_text(_render_md(c), encoding="utf-8")
        print(f"-> {a.md}")
    if a.html:
        Path(a.html).parent.mkdir(parents=True, exist_ok=True)
        Path(a.html).write_text(_render_html(c), encoding="utf-8")
        print(f"-> {a.html}")
    return 1 if c["n_errors"] else 0


# ---------------------------------------------------------------------------
# module: compila el binding pybind11 alcp_native con unas cabeceras dadas
# ---------------------------------------------------------------------------
def cmd_module(a: argparse.Namespace) -> int:
    try:
        import pybind11
    except ImportError:
        sys.exit("pybind11 no está instalado en este intérprete (use el .venv del proyecto)")
    src = Path(a.source).resolve()
    include = Path(a.include).resolve()
    if not src.is_file():
        sys.exit(f"no existe {src}")
    if not (include / "alcp" / "alcp.hpp").is_file():
        sys.exit(f"{include} no contiene alcp/alcp.hpp")
    out_dir = _label_dir(a.label) / "module"
    out_dir.mkdir(parents=True, exist_ok=True)
    ext = sysconfig.get_config_var("EXT_SUFFIX") or (".pyd" if IS_WIN else ".so")
    out = out_dir / f"alcp_native{ext}"
    py_inc = sysconfig.get_path("include")
    cmd = [a.cxx, "-shared", "-std=c++23", *_split_flags(a.cxxflags), f"-I{include}",
           f"-I{pybind11.get_include()}", f"-I{py_inc}", str(src)]
    if IS_WIN:
        # g++ (MinGW) + CPython de MSVC: se enlaza contra pythonXY de la
        # instalación base y se enlaza estáticamente GMP, libstdc++ y libgcc
        # para que el .pyd solo dependa de DLL del sistema y de pythonXY.dll.
        libs = Path(sys.base_prefix) / "libs"
        ver = f"{sys.version_info.major}{sys.version_info.minor}"
        cmd += [f"-L{libs}", f"-lpython{ver}", "-static", "-lgmpxx", "-lgmp"]
    else:
        cmd = cmd[:1] + ["-fPIC"] + cmd[1:] + ["-lgmpxx", "-lgmp"]
    cmd += ["-o", str(out)]
    print("$ " + " ".join(shlex.quote(c) for c in cmd), flush=True)
    t0 = time.perf_counter()
    cp = subprocess.run(cmd, check=False)
    dt = time.perf_counter() - t0
    if cp.returncode != 0:
        print(f"ERROR: la compilación del módulo falló (código {cp.returncode})", file=sys.stderr)
        return cp.returncode
    info = {"label": a.label, "kind": "module", "source": str(src), "source_sha256": _sha256_files([src]),
            "source_git": _git_info(src.parent), "include": str(include), "include_git": _git_info(include),
            "include_sha256": _sha256_files(sorted((include / "alcp").glob("*.hpp"))),
            "compiler": _compiler_info(a.cxx), "cxxflags": a.cxxflags, "command": cmd,
            "python": sys.version, "ext_suffix": ext, "build_seconds": round(dt, 2), "built": _now_iso()}
    (out_dir / "build_info.json").write_text(json.dumps(scrub_paths(info), indent=2, ensure_ascii=False),
                                             encoding="utf-8")
    # Comprobación rápida: importar el módulo en un proceso aparte.
    probe = ("import sys; sys.path.insert(0, sys.argv[1]); import alcp_native as m; "
             "print('import ok:', m.__file__, '| is_prime(97) =', m.is_prime(97))")
    cp = subprocess.run([sys.executable, "-c", probe, str(out_dir)], check=False)
    print(f"-> {out}  ({dt:.1f} s)")
    return 0 if cp.returncode == 0 else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _normalize_argv(argv: list[str]) -> list[str]:
    out = []
    i = 0
    while i < len(argv):
        if argv[i] in _FLAG_VALUE_OPTS and i + 1 < len(argv):
            out.append(f"{argv[i]}={argv[i + 1]}")
            i += 2
        else:
            out.append(argv[i])
            i += 1
    return out


def main(argv: list[str] | None = None) -> int:
    utf8_stdout()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    default_cxx = os.environ.get("CXX", "g++")

    b = sub.add_parser("build", help="compila bench_suite.cpp contra un directorio de cabeceras")
    b.add_argument("--include", required=True, help="directorio que contiene alcp/alcp.hpp")
    b.add_argument("--label", required=True, help="nombre de la variante (build/bench/<label>/)")
    b.add_argument("--cxxflags", default="-O3", help='flags de compilación (por defecto "-O3")')
    b.add_argument("--ldflags", default=None, help=f'flags de enlace (por defecto "{default_ldflags()}")')
    b.add_argument("--cxx", default=default_cxx, help="compilador (por defecto $CXX o g++)")
    b.set_defaults(func=cmd_build)

    def add_run_opts(p: argparse.ArgumentParser, reps: int | None, min_time: float | None) -> None:
        p.add_argument("--filter", default=None, help="expresión regular sobre el nombre del escenario")
        p.add_argument("--reps", type=int, default=reps, help="repeticiones mínimas por escenario")
        p.add_argument("--min-time", type=float, default=min_time, help="tiempo acumulado mínimo por escenario (s)")
        p.add_argument("--warmup", type=int, default=None, help="calentamientos (por defecto 1)")

    r = sub.add_parser("run", help="ejecuta la suite de una variante y guarda bench/results/<label>.json")
    r.add_argument("--label", required=True)
    add_run_opts(r, None, None)
    r.add_argument("--out", default=None, help="ruta del JSON (por defecto bench/results/<label>.json)")
    r.add_argument("--note", default=None, help='nota libre guardada en el JSON (p. ej. "smoke: no es una medición")')
    r.set_defaults(func=cmd_run)

    i = sub.add_parser("interleave", help="ejecuta variantes alternadas (A B C A B C ...) y agrega")
    i.add_argument("--labels", nargs="+", required=True)
    i.add_argument("--rounds", type=int, default=5)
    add_run_opts(i, 3, 0.0)
    i.add_argument("--per-scenario", action="store_true",
                   help="intercala escenario a escenario (más fino, más lento) en vez de suite completa")
    i.add_argument("--out-dir", default=None, help="directorio de salida (por defecto bench/results)")
    i.add_argument("--note", default=None, help="nota libre guardada en cada JSON")
    i.set_defaults(func=cmd_interleave)

    c = sub.add_parser("compare", help="compara resultados contra una base")
    c.add_argument("--base", required=True, help="etiqueta o ruta .json de la base")
    c.add_argument("--labels", nargs="+", required=True, help="etiquetas o rutas .json de las variantes")
    c.add_argument("--md", default=None, help="escribe la tabla en Markdown")
    c.add_argument("--html", default=None, help="escribe un informe HTML autocontenido")
    c.add_argument("--results-dir", default=None, help="directorio de resultados (por defecto bench/results)")
    c.set_defaults(func=cmd_compare)

    m = sub.add_parser("module", help="compila el módulo pybind11 alcp_native con unas cabeceras dadas")
    m.add_argument("--source", required=True, help="bindings.cpp a compilar")
    m.add_argument("--include", required=True, help="directorio que contiene alcp/alcp.hpp")
    m.add_argument("--label", required=True, help="el módulo va a build/bench/<label>/module/")
    m.add_argument("--cxxflags", default="-O3")
    m.add_argument("--cxx", default=default_cxx)
    m.set_defaults(func=cmd_module)

    args = ap.parse_args(_normalize_argv(sys.argv[1:] if argv is None else argv))
    if getattr(args, "rounds", 1) is not None and getattr(args, "rounds", 1) < 1:
        sys.exit("--rounds debe ser >= 1")
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
