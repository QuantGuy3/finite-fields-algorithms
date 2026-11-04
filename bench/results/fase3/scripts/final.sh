#!/usr/bin/env bash
# final.sh — Medición final de la fase 3 (documentación: así se obtuvieron los
# bench/results/fase3/final-*.json).  <REPO> es la raíz del repositorio y <SCRATCH>
# un directorio temporal; hay que sustituirlos antes de ejecutarlo.
# Si Smart App Control bloquea un ejecutable recién compilado, se vuelve a compilar
# (el nuevo enlace cambia su hash).
set -e
cd <REPO>
PY=.venv/Scripts/python.exe
# Cabeceras originales (sin tocar el worktree de la línea base).
git archive -o <SCRATCH>/orig.tar baseline-original include python-bindings/bindings.cpp
tar -xf <SCRATCH>/orig.tar -C <SCRATCH>
# Suite C++: referencia, final portable y final nativa, intercaladas.
$PY bench/run_bench.py build --include <SCRATCH>/include --label final-ref-O3 --cxxflags -O3
$PY bench/run_bench.py build --include include --label final-O3 --cxxflags -O3
$PY bench/run_bench.py build --include include --label final-native --cxxflags "-O3 -march=native"
$PY bench/run_bench.py interleave --labels final-ref-O3 final-O3 final-native --rounds 4 --reps 5 \
    --min-time 0.05 --out-dir bench/results/fase3
# Variantes solo GMP y un hilo (3 rondas alternadas entre sí).
for r in 1 2 3; do
  build/bench/final-O3/bench_suite.exe --reps 5 --min-time 0.05 --backend gmp --json <SCRATCH>/gmp-r$r.json
  ALCP_THREADS=1 build/bench/final-O3/bench_suite.exe --reps 5 --min-time 0.05 --json <SCRATCH>/t1-r$r.json
done
$PY bench/results/fase3/scripts/une_rondas.py bench/results/fase3/final-O3-gmp.json final-O3-gmp <SCRATCH>/gmp-r{1,2,3}.json
$PY bench/results/fase3/scripts/une_rondas.py bench/results/fase3/final-O3-t1.json final-O3-t1 <SCRATCH>/t1-r{1,2,3}.json
$PY bench/run_bench.py compare --base bench/results/fase3/final-ref-O3.json \
    --labels bench/results/fase3/final-O3.json bench/results/fase3/final-native.json \
    --md docs/optimizacion/final.md --html docs/optimizacion/final.html
# Módulo Python: original frente a final, 4 rondas alternadas.
$PY bench/run_bench.py module --source <SCRATCH>/python-bindings/bindings.cpp --include <SCRATCH>/include --label final-mod-ref
$PY bench/run_bench.py module --source python-bindings/bindings.cpp --include include --label final-mod
for r in 1 2 3 4; do
  $PY bench/bench_python.py --module-dir build/bench/final-mod-ref/module --json <SCRATCH>/ref-r$r.json --label py-final-ref
  $PY bench/bench_python.py --module-dir build/bench/final-mod/module --json <SCRATCH>/fin-r$r.json --label py-final
done
$PY bench/results/fase3/scripts/une_rondas.py bench/results/fase3/final-py-ref.json py-final-ref <SCRATCH>/ref-r{1,2,3,4}.json
$PY bench/results/fase3/scripts/une_rondas.py bench/results/fase3/final-py.json py-final <SCRATCH>/fin-r{1,2,3,4}.json
$PY bench/run_bench.py compare --base bench/results/fase3/final-py-ref.json --labels bench/results/fase3/final-py.json \
    --md docs/optimizacion/final-python.md
# Capa Python con el módulo final.
$PY bench/bench_alcp.py --pkg-dir python-bindings --json bench/results/fase3/final-alcp.json --label final-alcp
# Casos históricos (bench_historico.cpp), 2 rondas alternadas de 5 repeticiones.
g++ -std=c++23 -Wall -Wextra -O3 -I<SCRATCH>/include bench/bench_historico.cpp -static -lgmpxx -lgmp -o <SCRATCH>/hist_ref.exe
g++ -std=c++23 -Wall -Wextra -O3 -Iinclude bench/bench_historico.cpp -static -lgmpxx -lgmp -o <SCRATCH>/hist_fin.exe
for r in 1 2; do
  <SCRATCH>/hist_ref.exe 5 > <SCRATCH>/ref-r$r.txt
  <SCRATCH>/hist_fin.exe 5 > <SCRATCH>/fin-r$r.txt
  ALCP_THREADS=1 <SCRATCH>/hist_fin.exe 5 > <SCRATCH>/fin1-r$r.txt
done
