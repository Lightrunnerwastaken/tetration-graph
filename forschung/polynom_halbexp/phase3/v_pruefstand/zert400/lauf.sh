#!/bin/bash
# Strang V (Phase 3) — Zweig-Zertifikat 193 -> 400, Kopie von Strang X (x_alternativen/zert/lauf.sh), eigener Ordner und venv halbexp-venv.
# 12 Arbeiter (Server 28 Kerne, >= 16 bleiben frei), laengste Uebergaenge zuerst; danach Kettenpruefung.
cd ~/v_pruefstand/zert400 || exit 1
mkdir -p log boxen erg
date > log/START
seq 399 -1 193 | xargs -P 12 -I{} sh -c 'timeout 10800 ~/halbexp-venv/bin/python zert_job.py {} > log/job_{}.log 2>&1'
date > log/FERTIG
~/halbexp-venv/bin/python zert_kette.py 193 400 > log/kette.log 2>&1
date >> log/FERTIG
