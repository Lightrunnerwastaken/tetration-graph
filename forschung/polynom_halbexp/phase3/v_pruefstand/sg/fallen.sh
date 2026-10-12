#!/bin/bash
# Strang V — Fallen-Tests zu (S_γ): Abbildungsgenauigkeit, gestörte Abbildung, anderes Linsenprofil, mehr Präzision. Ein Prozess, nacheinander.
cd ~/v_pruefstand/sg/a/b || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=~/halbexp-venv/bin/python
V12_MFFT=4096 $PY -u v12b.py kubisch 400,800 > f1_mfft4096.log 2>&1
V12_STOER=1e-6 $PY -u v12b.py kubisch 400,800 > f2_stoer1e-6.log 2>&1
V12_STOER=1e-4 $PY -u v12b.py kubisch 400,800 > f3_stoer1e-4.log 2>&1
$PY -u v12b.py kubisch 400,800 0.0166 12 1.3 sin > f4_sinlinse.log 2>&1
$PY -u v12b.py kubisch 800 0.0166 12 1.3 K 250 > f5_bits_plus250.log 2>&1
date > fallen_FERTIG
