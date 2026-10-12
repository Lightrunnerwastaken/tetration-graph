#!/bin/bash
# Strang V — (K2) Neulauf mit Ks Korrektur F1 (Mittelpunkt-Scheibe ⊂ Zellscheibe), 8 Arbeiter, Threads begrenzt.
cd ~/v_pruefstand/k2b || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p log
date > log/START
cat k2_jobs_F1.txt | xargs -P 8 -L 1 sh -c 'timeout 7200 ~/halbexp-venv/bin/python -u v06_k2_rand.py "$0" "$1" "$2" "$3" "$4" "$5" "$6" "$7" > log/k2_$0_$3.log 2>&1'
~/halbexp-venv/bin/python -u v07_k2_zusammen.py > log/zusammen.log 2>&1
date > log/FERTIG
