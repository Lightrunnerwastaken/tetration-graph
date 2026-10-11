#!/bin/bash
# Strang X — parallele Verlängerung des Zweig-Zertifikats 193 → 400 (Auftrag Leitung). 12 Arbeiter, längste Übergänge zuerst.
cd ~/tetration-zert || exit 1
mkdir -p log boxen erg
date > log/START
seq 399 -1 193 | xargs -P 12 -I{} sh -c 'timeout 10800 ~/tetration-venv/bin/python zert_job.py {} > log/job_{}.log 2>&1'
date > log/FERTIG
~/tetration-venv/bin/python zert_kette.py 193 400 > log/kette.log 2>&1
date >> log/FERTIG
