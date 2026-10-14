"""Führender Term n·a_n ≈ −h'(c − ln n) für astronomische n (Kettenregel h'(v) = h'(log v)/v, Hauptzweige).
Zeigt: n·ln n·|a_n| → 0 (wie 1/(ln ln n · ...)), arg a_n hat keinen Grenzwert, sondern wandert mit slog(ln ln n).
Aufruf: python p6_arg.py
"""
import json
from pathlib import Path
import mpmath as mp
from kern import Fix
mp.mp.dps = 50
fx = Fix(mp.e, K=40)
kor = json.loads(Path("daten/korrektur_e.json").read_text())
A = mp.mpf(kor["A"])
c = mp.log(1j * A / (2 * mp.pi))

def hprime(v):
    acc = mp.mpc(1)
    while abs(v) > 30:
        acc /= v
        v = mp.log(v)
    return acc * mp.diff(fx.h, v)

rows = []
for lnln in ["2.12", "3", "5", "10", "15.15", "30", "100", "1000", "3814279.1"]:
    t = mp.exp(mp.mpf(lnln))           # t = ln n
    T = -hprime(c - t)                 # ≈ n·a_n
    rows.append([lnln, mp.nstr(abs(T) * t, 8), mp.nstr(abs(T) * t * mp.log(t), 8), mp.nstr(mp.arg(T), 8)])
    print(f"ln ln n = {lnln:>10}:  n ln n|a_n| = {rows[-1][1]:>12}  n ln n lnln n|a_n| = {rows[-1][2]:>12}  arg = {rows[-1][3]}")
Path("daten/arg_gross_n.json").write_text(json.dumps(rows, indent=1))
