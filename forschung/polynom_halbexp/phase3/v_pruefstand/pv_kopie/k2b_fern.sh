#!/bin/bash
cd ~/v_pruefstand/k2b || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cat k2_jobs_F1_fern.txt | xargs -P 6 -L 1 sh -c 'timeout 7200 ~/halbexp-venv/bin/python -u v06_k2_rand.py "$0" "$1" "$2" "$3" > log/k2_$0_$3.log 2>&1'
date > log/FERN_FERTIG
