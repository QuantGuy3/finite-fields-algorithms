# Driver de intercalado para la etapa de medición serial (fase 3).
#
# Reutiliza run_bench.py (harness de med-a6, idéntico al de 5f49217 salvo
# bench_suite.cpp) y añade "variantes virtuales":
#   ETIQUETA=BUILD[:backend=gmp|auto][:threads=N]
# p. ej. a6gmp=a6:backend=gmp  o  m9t1=m9:threads=1
# BUILD se busca en med-a6/build/bench y, si no está, en med-t0/build/bench.
#
# Uso: python inter.py --out-dir DIR --rounds N --reps R --min-time T [--filter RE] VAR1 VAR2 ...
import argparse, os, sys, json, time
from pathlib import Path

S = Path(r"<SCRATCH>")
sys.path.insert(0, str(S / "med-a6" / "bench"))
import run_bench as rb  # noqa: E402

ROOTS = [S / "med-a6" / "build" / "bench", S / "med-t0" / "build" / "bench"]
VARS: dict[str, dict] = {}


def parse_var(txt: str) -> str:
    label, _, spec = txt.partition("=")
    parts = (spec or label).split(":")
    build = parts[0]
    opts = dict(p.split("=", 1) for p in parts[1:])
    root = next((r for r in ROOTS if (r / build / rb.EXE_NAME).is_file()), None)
    if root is None:
        sys.exit(f"no existe el build {build}")
    VARS[label] = {"dir": root / build, "backend": opts.get("backend"), "threads": opts.get("threads"),
                   "build": build}
    return label


def _label_dir(label: str) -> Path:
    return VARS[label]["dir"]


def run_suite(label, *, filt, reps, min_time, warmup, echo=True):
    v = VARS[label]
    exe = v["dir"] / rb.EXE_NAME
    tmp = S / "med-logs" / f"last_{label}.json"
    if tmp.exists():
        tmp.unlink()
    cmd = [str(exe), "--json", str(tmp)]
    if v["backend"]:
        cmd += ["--backend", v["backend"]]
    if filt:
        cmd += ["--filter", filt]
    if reps is not None:
        cmd += ["--reps", str(reps)]
    if min_time is not None:
        cmd += ["--min-time", str(min_time)]
    if warmup is not None:
        cmd += ["--warmup", str(warmup)]
    env = dict(os.environ)
    env.pop("ALCP_THREADS", None)
    if v["threads"]:
        env["ALCP_THREADS"] = v["threads"]
    import subprocess
    t0 = time.perf_counter()
    cp = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
    if echo:
        print(f"  [{label}] rc={cp.returncode} {time.perf_counter() - t0:.1f} s", flush=True)
    if not tmp.is_file():
        print(cp.stdout.decode("utf-8", "replace")[-2000:])
        sys.exit(f"[{label}] sin resultados")
    data = json.loads(tmp.read_text(encoding="utf-8"))
    data["_returncode"] = cp.returncode
    return data


def _load_build_info(label):
    p = VARS[label]["dir"] / "build_info.json"
    bi = json.loads(p.read_text(encoding="utf-8")) if p.is_file() else {}
    bi["variante"] = {k: VARS[label][k] for k in ("build", "backend", "threads")}
    return bi


rb._label_dir = _label_dir
rb.run_suite = run_suite
rb._load_build_info = _load_build_info

ap = argparse.ArgumentParser()
ap.add_argument("vars", nargs="+")
ap.add_argument("--rounds", type=int, default=3)
ap.add_argument("--reps", type=int, default=3)
ap.add_argument("--min-time", type=float, default=0.1)
ap.add_argument("--warmup", type=int, default=None)
ap.add_argument("--filter", default=None)
ap.add_argument("--out-dir", required=True)
ap.add_argument("--note", default=None)
a = ap.parse_args()
labels = [parse_var(t) for t in a.vars]
ns = argparse.Namespace(labels=labels, rounds=a.rounds, reps=a.reps, min_time=a.min_time, warmup=a.warmup,
                        filter=a.filter, per_scenario=False, out_dir=a.out_dir, note=a.note)
sys.exit(rb.cmd_interleave(ns))
