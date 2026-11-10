#!/bin/bash
# Segunda tanda, con el portátil ya enchufado (desde 17:06:29).
S=<SCRATCH>
R=<REPO>
PY=$R/.venv/Scripts/python.exe
O=$R/bench/results/fase3
L=$S/med-logs
cd $L
FA='^(fp\.p101\.(elem_muladd|elem_pow|mul_school\.d512|mul_karatsuba\.d512|divmod|gcd\.|gcd_ext|factor)|fp\.m61\.(elem_muladd|elem_inv|elem_pow|mul_school\.d512|divmod|factor)|fp\.ntt\.(elem_pow|mul_school\.d512|mul_fft\.d4096|factor)|fp\.m127\.(elem_muladd|mul_school\.d512|divmod)|fp\.p256\.(elem_muladd|mul_school\.d512|divmod)|fq\.q3e30\.(elem_muladd|mul_school|divmod)|fq\.q2e8\.(elem_muladd|mul_school|factor)|fq\.q7e4\.(mul_school|gcd)|fq\.m61e2\.(elem_muladd|mul_school|divmod|mul_karatsuba)|fast\.p7\.toeplitz_divmod)'
FM='^(fp\.ntt\.(mul_karatsuba\.d512|mul_fft|mul_auto|gcd_ext|pow_mod|factor)|fp\.p101\.(mul_karatsuba|mul_auto|gcd_ext|mul_school\.d512)|fp\.p256\.(mul_fft\.d4096|mul_auto|gcd_ext|pow_mod)|fq\.q3e30\.mul_karatsuba|fq\.q2e8\.mul_karatsuba|fq\.m61e2\.(mul_karatsuba|mul_fft|mul_auto|gcd|divmod)|fq\.q7e4\.(gcd|factor)|fast\.ntt\.|fp\.m61\.elem_muladd)'
echo "INICIO ac $(date +%T)"
mkdir -p $O/aritmetica_ac $O/extremos_ac $O/multiplicacion_ac
$PY inter2.py --out-dir $O/aritmetica_ac --rounds 3 --reps 3 --min-time 0 --filter "$FA" --note "fase3 cadena aritmetica (corriente)" h0 a1 a2 a3 a4 a4gmp=a4:backend=gmp a5 a6 > $L/run-arit-ac.log 2>&1
echo "FIN aritmetica_ac rc=$? $(date +%T)"
for r in 4 5 6; do for m in modb0 modh0 moda6 modm9; do
  ti=$(date +%T); $PY $S/med-a6/bench/bench_python.py --module-dir $S/med-a6/build/bench/$m/module --label py-$m --reps 5 --min-time 0.1 --json $O/modulo/py-$m-r$r.json --note "fase3 ronda $r (corriente)" > $L/run-mod-$m-$r.log 2>&1; echo "mod $m r$r ini=$ti fin=$(date +%T)" >> $L/ejecuciones.log
done; done
echo "FIN modulo r4-6 $(date +%T)"
$PY inter2.py --out-dir $O/extremos_ac --rounds 2 --reps 3 --min-time 0 --note "fase3 extremos suite completa (corriente)" b0 h0 a6 a6gmp=a6:backend=gmp m9 m9t1=m9:threads=1 > $L/run-extremos-ac.log 2>&1
echo "FIN extremos_ac rc=$? $(date +%T)"
$PY inter2.py --out-dir $O/multiplicacion_ac --rounds 2 --reps 3 --min-time 0 --filter "$FM" --note "fase3 cadena multiplicacion (corriente)" h0 m1 m2 m3 m4 m5 m6 m7 m8 m9 m8t1=m8:threads=1 m9t1=m9:threads=1 > $L/run-mult-ac.log 2>&1
echo "FIN multiplicacion_ac rc=$? $(date +%T)"
