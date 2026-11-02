#!/bin/bash
# Campaña serial de la fase 3 (se lanza tras el intercalado de extremos).
S=<SCRATCH>
R=<REPO>
PY=$R/.venv/Scripts/python.exe
O=$R/bench/results/fase3
L=$S/med-logs
cd $L
bat() { powershell -NoProfile -Command "(Get-CimInstance Win32_Battery).EstimatedChargeRemaining" | tr -d '\r'; }
echo "INICIO campaña $(date +%T) bateria=$(bat)%"

# Filtro aritmética (subconjunto representativo de los escenarios afectados + control)
FA='^(fp\.p101\.(elem_muladd|elem_pow|mul_school\.d512|mul_karatsuba\.d512|divmod|gcd\.|gcd_ext|factor)|fp\.m61\.(elem_muladd|elem_inv|elem_pow|mul_school\.d512|divmod|factor)|fp\.ntt\.(elem_pow|mul_school\.d512|mul_fft\.d4096|factor)|fp\.m127\.(elem_muladd|mul_school\.d512|divmod)|fp\.p256\.(elem_muladd|mul_school\.d512|divmod)|fq\.q3e30\.(elem_muladd|mul_school|divmod)|fq\.q2e8\.(elem_muladd|mul_school|factor)|fq\.q7e4\.(mul_school|gcd)|fq\.m61e2\.(elem_muladd|mul_school|divmod|mul_karatsuba)|fast\.p7\.toeplitz_divmod)'
mkdir -p $O/aritmetica
$PY inter2.py --out-dir $O/aritmetica --rounds 3 --reps 3 --min-time 0 --filter "$FA" --note "fase3 cadena aritmetica" \
  h0 a1 a2 a3 a4 a4gmp=a4:backend=gmp a5 a6 > $L/run-arit.log 2>&1
echo "FIN aritmetica rc=$? $(date +%T) bateria=$(bat)%"

# Filtro multiplicación
FM='^(fp\.ntt\.(mul_karatsuba\.d512|mul_fft|mul_auto|gcd_ext|pow_mod|factor)|fp\.p101\.(mul_karatsuba|mul_auto|gcd_ext|mul_school\.d512)|fp\.p256\.(mul_fft\.d4096|mul_auto|gcd_ext|pow_mod)|fq\.q3e30\.mul_karatsuba|fq\.q2e8\.mul_karatsuba|fq\.m61e2\.(mul_karatsuba|mul_fft|mul_auto|gcd|divmod)|fq\.q7e4\.(gcd|factor)|fast\.ntt\.|fp\.m61\.elem_muladd)'
mkdir -p $O/multiplicacion
$PY inter2.py --out-dir $O/multiplicacion --rounds 3 --reps 3 --min-time 0 --filter "$FM" --note "fase3 cadena multiplicacion" \
  h0 m1 m2 m3 m4 m5 m6 m7 m8 m9 m8t1=m8:threads=1 m9t1=m9:threads=1 > $L/run-mult.log 2>&1
echo "FIN multiplicacion rc=$? $(date +%T) bateria=$(bat)%"

# Módulo Python (bench_python.py, rondas alternadas)
mkdir -p $O/modulo
for r in 1 2 3; do for m in modb0 modh0 moda6 modm9; do
  ti=$(date +%T); $PY $S/med-a6/bench/bench_python.py --module-dir $S/med-a6/build/bench/$m/module --label py-$m --reps 5 --min-time 0.1 \
     --json $O/modulo/py-$m-r$r.json --note "fase3 ronda $r" > $L/run-mod-$m-$r.log 2>&1; echo "  mod $m r$r rc=$? $(date +%T)"; echo "mod $m r$r ini=$ti fin=$(date +%T)" >> $L/ejecuciones.log
done; done
echo "FIN modulo $(date +%T) bateria=$(bat)%"

# Capa Python (bench_alcp.py, mismo .pyd de 5f49217)
mkdir -p $O/capa_python
P=$S/med-pkg2
echo "capa ini=$(date +%T)" >> $L/ejecuciones.log; $PY $S/med-pc/bench/bench_alcp.py --interleave $(for d in $P/q*; do b=$(basename $d); echo "${b%%-*}=$d"; done) \
   --rounds 3 --reps 3 --min-time 0.05 --out-dir $O/capa_python > $L/run-capa.log 2>&1
echo "FIN capa rc=$? $(date +%T) bateria=$(bat)%"

# Toolchain (perfiles de flags de la rama wip, sobre la base original y 5f49217)
FT='^(fp\.p101\.(elem_muladd|mul_school\.d512|mul_karatsuba\.d512|divmod|gcd\.|factor)|fp\.ntt\.(mul_fft\.d4096|pow_mod)|fp\.m61\.(elem_muladd|mul_school\.d512)|fp\.p256\.(elem_muladd|mul_fft\.d4096)|fq\.q2e8\.(elem_muladd|mul_school|factor)|fq\.m61e2\.(mul_school|mul_fft\.d1024)|fast\.ntt\.(dft\.k12|toeplitz_lower_inverse\.n1024))'
mkdir -p $O/toolchain
$PY inter2.py --out-dir $O/toolchain --rounds 3 --reps 3 --min-time 0 --filter "$FT" --note "fase3 toolchain" \
  tb0-O3 tb0-v3 tb0-nat tb0-lto tb0-alloc th0-O3 th0-v3 th0-nat th0-lto th0-alloc > $L/run-tool.log 2>&1
echo "FIN toolchain rc=$? $(date +%T) bateria=$(bat)%"
