#!/bin/bash
# Tercera tanda (corriente): mediciones focalizadas de las técnicas dudosas, más rondas y repeticiones.
S=<SCRATCH>
R=<REPO>
PY=$R/.venv/Scripts/python.exe
O=$R/bench/results/fase3/focal
L=$S/med-logs
cd $L
mkdir -p $O/barrett $O/sumaresta $O/karatsuba_vistas $O/ntt_par
echo "INICIO focal $(date +%T)"
# 2297387 (Barrett sin desplazamientos variables): elementos y operaciones de Fp64
$PY inter2.py --out-dir $O/barrett --rounds 6 --reps 5 --min-time 0.2 --filter '^fp\.(p7|p101|ntt|m61)\.(elem_|mul_school\.d512|divmod|factor)' a4 a5 > $L/run-f-barrett.log 2>&1
echo "FIN barrett rc=$? $(date +%T)"
# dcb905c (suma/resta de Poly): Karatsuba y gcd en ambos backends
$PY inter2.py --out-dir $O/sumaresta --rounds 6 --reps 5 --min-time 0.2 --filter 'mul_karatsuba|gcd' a2 a3 > $L/run-f-sumaresta.log 2>&1
echo "FIN sumaresta rc=$? $(date +%T)"
# 9bfc1e7 (Karatsuba sobre vistas) en F_q
$PY inter2.py --out-dir $O/karatsuba_vistas --rounds 6 --reps 5 --min-time 0.2 --filter '^fq\..*mul_karatsuba|^fp\.(p7|p101|m61|m127|p256)\.mul_karatsuba' m2 m3 > $L/run-f-kvistas.log 2>&1
echo "FIN karatsuba_vistas rc=$? $(date +%T)"
# b0a9fe7 (NTT paralela) en escenarios FFT pequeños y Karatsuba de m61e2
$PY inter2.py --out-dir $O/ntt_par --rounds 6 --reps 5 --min-time 0.2 --filter 'fast\.ntt\.toeplitz|fq\.m61e2\.(mul_karatsuba|gcd|mul_fft\.d1024)|fp\.ntt\.mul_fft\.d4096' m8 m9 m9t1=m9:threads=1 > $L/run-f-nttpar.log 2>&1
echo "FIN ntt_par rc=$? $(date +%T)"
# Rondas adicionales con corriente para los extremos y la cadena de multiplicación
O2=$R/bench/results/fase3
mkdir -p $O2/extremos_ac_extra
$PY inter2.py --out-dir $O2/extremos_ac_extra --rounds 1 --reps 3 --min-time 0 --note "fase3 extremos, ronda adicional (corriente)" b0 h0 a6 a6gmp=a6:backend=gmp m9 m9t1=m9:threads=1 > $L/run-extremos-ac-extra.log 2>&1
echo "FIN extremos_ac_extra rc=$? $(date +%T)"
# 3660189 (umbrales) y e8e07e6 (Karatsuba paralelo): gcd y mul_auto, donde la cadena con corriente dio pérdidas
mkdir -p $O/umbrales_par
$PY inter2.py --out-dir $O/umbrales_par --rounds 6 --reps 5 --min-time 0.2 --filter 'gcd|mul_auto|fq\.q7e4\.factor' m5 m6 m7 m8 m8t1=m8:threads=1 > $L/run-f-umbrales.log 2>&1
echo "FIN umbrales_par rc=$? $(date +%T)"
