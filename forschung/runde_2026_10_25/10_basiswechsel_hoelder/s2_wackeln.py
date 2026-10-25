"""S2: Wackel-Lemma (Satz A) numerisch: Differenzenquotient von H auf der Skala h = s/U_n'(y) folgt
       D(y, h) = Pi_n(y) * (1 + W(U_n, s)) * (1 + o(1)),   W(u, s) = (P4(u+s) - P4(u)) / s,
mit Pi_n(y) unabhängig von s. Gemessen wird Q(s) = D(y, s/U_n') / (1 + W(U_n, s)); Q muss in s konstant sein.
Vergleich: D(s)/D(s_gross) - 1 gegen W(U_n, s) (ohne Normierung durch die Vorhersage).

Aufruf: python s2_wackeln.py <b> <d> <n> <dps> th1 th2 ...      z. B.  python s2_wackeln.py exp(1) 2 3 80 0.05 0.3
Ergebnis: daten/s2_<b>_<d>_n<n>_dps<dps>.json
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp

import lib

b, d, n, dps = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
ths = sys.argv[5:]
mp.mp.dps = dps
Tb, Td = lib.Tet(b), lib.Tet(d)
SS = ["0.001", "0.01", "0.1", "0.25", "0.5", "0.75", "1", "1.25", "1.5", "2", "3", "5", "10", "100", "1000"]


def W(u, s):
    k = mp.floor(u)
    return (lib.phi4(Tb, Td, u + s - mp.floor(u + s)) - lib.phi4(Tb, Td, u - k)) / s


aus = []
for th in ths:
    t0 = time.time()
    th = mp.mpf(th)
    y = Td.P(th)
    U, dU = y, mp.mpf(1)
    for _ in range(n):
        dU *= Td.dT(U)
        U = Td.T(U)
    H0 = lib.H(Tb, Td, y)
    zeilen = []
    for s in SS:
        s = mp.mpf(s)
        h = s / dU
        D = (lib.H(Tb, Td, y + h) - H0) / h
        w = W(U, s)
        zeilen.append((s, D, w, D / (1 + w)))
    Dg = zeilen[-1][1]
    Qg = zeilen[-1][3]
    print(f"th={mp.nstr(th, 4)}  U_{n}={mp.nstr(U, 8)}  U_{n}'={mp.nstr(dU, 6)}   [{time.time()-t0:.0f}s]")
    print(f"   {'s':>7} {'D/D(1000)-1':>16} {'W(U_n,s)':>16} {'Q(s)/Q(1000)-1':>16}")
    for s, D, w, Q in zeilen:
        print(f"   {mp.nstr(s, 4):>7} {mp.nstr(D / Dg - 1, 8):>16} {mp.nstr(w, 8):>16} {mp.nstr(Q / Qg - 1, 3):>16}")
    rest = max(abs(Q / Qg - 1) for _, _, _, Q in zeilen)
    amp = max(abs(D / Dg - 1) for _, D, _, _ in zeilen)
    print(f"   Wackel-Amplitude max|D/D(1000)-1| = {mp.nstr(amp, 4)},  Rest nach Vorhersage max|Q/Q(1000)-1| = {mp.nstr(rest, 3)}")
    aus.append({"th": mp.nstr(th, 10), "U": mp.nstr(U, 20), "dU": mp.nstr(dU, 20), "amp": mp.nstr(amp, 8), "rest": mp.nstr(rest, 5),
                "zeilen": [[mp.nstr(v, 20) for v in z] for z in zeilen]})
out = Path(__file__).resolve().parent / "daten" / f"s2_{lib._tag(b)}_{lib._tag(d)}_n{n}_dps{dps}.json"
out.write_text(json.dumps(aus, indent=1), encoding="utf-8")
