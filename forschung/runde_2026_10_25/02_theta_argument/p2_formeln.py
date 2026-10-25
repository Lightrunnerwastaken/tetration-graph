"""Gebiet 02, Teil D: explizite Formeln für arg a_n über θ' bzw. h' auf der reellen Achse.

F0 (Etappe 2, führend):  arg a_n = π + arg h'(c − ln n)                       , c = ln(A/2π) + iπ/2   (Fehler O(1/ln² n))
F1 (neu, Achse bei σ_n): arg a_n ≈ π + Φ(σ_n) + π/(2 ln n)                   , Φ(s) = arg(1+θ'(s+i0)) = arg h'(sexp(s−2)+i0),
                         σ_n = 1 + slog(A/(2πn)) ≈ 1/(2πn), also Φ(σ_n) = arg h'(ln(A/2πn) + i0)
F2 (neu, Gerade Im = π): arg a_n ≈ P(ln(ln n − Re c)) + π/(2 ln n)          , P(Y) = arg h'(Y + iπ)
F3 (neu, selbstähnlich): arg a_n ≈ Φ(slog n) − π/ln ln n                      , Φ(slog n) = arg h'(n + i0) = arg h'(ln ln ln n + i0)
Zusätzlich Modul: |a_n| ≈ |1+θ'(σ_n)| / (n · sexp'(σ_n − 2)) und a_n / h'(n+i0) → 1.
Aufruf: python p2_formeln.py [e|2]
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

from kern import Fix

B = sys.argv[1] if len(sys.argv) > 1 else "e"
E2 = Path(r"D:\Users\janis\Documents\Matura\tetration-graph\forschung\identitaeten_2026_10_14\etappe2\daten")
HIER = Path(__file__).resolve().parent
mp.mp.dps = 30
fx = Fix(mp.e if B == "e" else mp.mpf(B), K=40)
lb, lam = fx.lb, fx.lam
kor = json.loads((E2 / f"korrektur_{B}.json").read_text())
A = mp.mpf(kor["A"])
c = mp.log(1j * A / (2 * mp.pi * lb))     # Basis b: ℓ2·ln b = c + log(1−q)
d = json.loads((E2 / f"a_{B}_80_0.003_16384.json").read_text())
a = [mp.mpc(mp.mpf(r), mp.mpf(i)) for r, i in d["a"]]
EPS = mp.mpf("1e-22")


def hp(w):
    """h'(w) per Kettenregel h'(w) = h'(log_b w)/(w ln b), Hauptzweige, bis nahe L; dann Schröder."""
    w = mp.mpc(w)
    acc = mp.mpc(1)
    n = 0
    while abs(w - fx.L) > mp.mpf("1e-3") or n < 5:
        acc /= w * lb
        w = mp.log(w) / lb
        n += 1
    u = w - fx.L
    tol = mp.mpf(10) ** (-mp.mp.dps)
    while abs(u) ** (fx.K + 1) > tol * abs(u):
        acc /= (fx.L + u) * lb
        u = mp.log(fx.L + u) / lb - fx.L
    S = fx.S
    su = mp.polyval(S[::-1], u)
    dsu = mp.polyval([k * S[k] for k in range(len(S) - 1, 0, -1)], u)
    return acc * dsu / su / mp.log(lam)


def arg_mod(x):
    x = float(x)
    return (x + np.pi) % (2 * np.pi) - np.pi


def formeln(t):
    """t = ln n (darf riesig sein). Gibt F0, F1, F3 (Argumente) und die Modulquotienten zurück."""
    v = (c - t) / lb                               # Argument von h im führenden Term (Basis b: h((c−ln n)/ln b))
    F0 = mp.pi + mp.arg(hp(v))
    vr = mp.mpc(mp.re(v), EPS)                      # dieselbe Stelle, aber auf der Achse (von oben)
    F1 = mp.pi + mp.arg(hp(vr)) + mp.atan(mp.im(v) / -mp.re(v))
    # F3: h'(n + i0) über drei Logarithmen: h'(n) = h'(lnlnln n)/(n ln n lnln n), Argument nur vom innersten Faktor
    x3 = mp.log(mp.log(t / lb) / lb) / lb if B == "e" else None
    if B == "e":
        F3 = mp.arg(hp(mp.mpc(x3, EPS))) - mp.pi / mp.log(t)
    else:
        F3 = mp.nan
    return F0, F1, F3


res = {"basis": B, "daten": [], "gross": []}
print("== n ≤ 4096: Daten gegen Formeln (Argumente, mod 2π) ==")
print("    n | arg a_n   | F0        | F1 (σ_n)  | F3 (slog n) | F1−F0     | F3−F0")
for n in [4, 16, 64, 256, 1024, 2048, 4096]:
    F0, F1, F3 = formeln(mp.log(n))
    an = mp.arg(a[n])
    row = [n, float(an), arg_mod(F0), arg_mod(F1), arg_mod(F3)]
    res["daten"].append(row)
    print(f"{n:5d} | {row[1]:+.6f} | {row[2]:+.6f} | {row[3]:+.6f} | {row[4]:+.6f}   | {arg_mod(F1-F0):+.2e} | {arg_mod(F3-F0):+.3f}")

print("\n== grosse n (nur Formeln; F0 ist bis O(1/ln² n) exakt) ==")
print(" lnln n   | F0        | F1−F0      | (F1−F0)·ln n·lnln n | F3−F0     | (F3−F0)·lnln n·lnlnln n")
for L2 in ["2.12", "3", "5", "10", "15.154", "30", "100", "1000", "1e4", "3814279.1", "1e8"]:
    t = mp.exp(mp.mpf(L2))
    F0, F1, F3 = formeln(t)
    e1 = arg_mod(F1 - F0)
    e3 = arg_mod(F3 - F0)
    l2 = mp.mpf(L2)
    row = [L2, arg_mod(F0), e1, float(e1 * t * mp.log(t)), e3, float(e3 * l2 * mp.log(l2))]
    res["gross"].append(row)
    s4 = f"{row[3]:+.4f}" if row[3] is not None else "   –   "
    print(f"{L2:>9} | {row[1]:+.6f} | {e1:+.3e} | {s4:>12} | {e3:+.5f} | {row[5]:+.4f}")
(HIER / "daten").mkdir(exist_ok=True)
(HIER / "daten" / f"p2_{B}.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
