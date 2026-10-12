#!/bin/bash
# Strang V — v59: CF-Koeffizienten (6 Linsen) + adjungierte Köpfe als Vektoren (kubisch, Kneser) + Gewichte auf den Linsen
cd ~/v_pruefstand/sg/a/b
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONIOENCODING=utf-8
PY=~/halbexp-venv/bin/python
for S in 0.0126 0.0137 0.0147 0.0166 0.022 0.03; do
  $PY v59c_koeff.py 200,400,600,800,1000,1500,2000 $S > v59c_S$S.log 2>&1 &
done
for M in kubisch kneser; do
  ( $PY v59_elim.py $M 2001 > v59_elim_$M.log 2>&1
    for S in 0.0137 0.0126 0.0147 0.0166 0.022 0.03; do $PY v59_eval.py $M 2001 $S > v59_eval_${M}_S$S.log 2>&1; done ) &
done
wait
echo ENDE > v59_ende.txt
