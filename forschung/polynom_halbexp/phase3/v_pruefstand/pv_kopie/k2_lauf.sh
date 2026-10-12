#!/bin/bash
# Strang V — (K2): ∫|∂_θ⁶ f_K| auf |x| = ρ in 19 Abschnitten, 8 Arbeiter (Zertifikat 193→400 belegt 12 Kerne, 8 bleiben frei).
cd ~/v_pruefstand/k2 || exit 1
mkdir -p log
date > log/START
cat k2_jobs.txt | xargs -P 8 -L 1 sh -c 'timeout 10800 ~/halbexp-venv/bin/python -u v06_k2_rand.py "$0" "$1" "$2" "$3" > log/k2_$0_$3.log 2>&1'
date > log/FERTIG
