"""Strang P, Prüfung 4 (Kneser-Motor, nur Strang P): unabhängige Nachrechnung von Daten, auf denen Beweise stehen.

M1  Kneser-Koeffizientendatei (4000) gegen den Motor: f_K(x) = sexp(slog(x) + 1/2) an Punkten in D (auch nahe ∂D).
M2  b35b′: f_K(e^z) auf dem linken Bogen von ∂D (Motor) — weg von (−∞, 0]?
M3  Fortsetzung auf f_K(D̄) (für das Kompaktum K' = D̄ ∪ f_K(D̄)): f_K(w) = exp(h(w)) für w = f_K(x), x ∈ Γout.
"""
import sys
import time
from pathlib import Path

import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
EXP = Path(r"D:\Users\janis\Documents\Matura\tetration-graph\experimente")
sys.path.insert(0, str(EXP))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import gemeinsam  # noqa: E402
from reihen import horner, kneser_koeff, fixpunkt_L  # noqa: E402

mp.mp.dps = 45
A = kneser_koeff(4000)
L = fixpunkt_L()
rho = abs(L)


def fK_motor(ws):
    s = gemeinsam.kneser("slog", "exp(1)", ws, dps=40)
    t = gemeinsam.kneser("sexp", "exp(1)", [x + mp.mpf("0.5") for x in s], dps=40)
    return t


t0 = time.time()
print("M1: Reihe (4000 Koeff.) gegen Motor")
pts = [mp.mpc("0.3", "0.5"), mp.mpc("-1.0", "0.6"), mp.mpc("0", "1.2"), mp.mpc("-1.2", "-0.3"),
       0.98 * rho * mp.expj(mp.radians(76)), 0.98 * rho * mp.expj(mp.radians(150)), 0.95 * rho * mp.expj(mp.radians(-20))]
mot = fK_motor(pts)
for p, m in zip(pts, mot):
    s = horner(A, p)
    print(f"  x = {mp.nstr(p, 6):>24}  |Reihe − Motor| = {mp.nstr(abs(s - m), 3)}   f = {mp.nstr(m, 10)}")
print(f"  ({time.time() - t0:.1f} s)")

print("M2: f_K(e^z) auf dem linken Bogen (Motor)")
thL = mp.arg(L)
zs = []
for i in range(13):
    th = thL + mp.mpf("0.0146") + (2 * mp.pi - 2 * thL - 2 * mp.mpf("0.0146")) * i / 12
    zs.append(rho * mp.expj(th))
mot = fK_motor([mp.exp(z) for z in zs])
mina, minarg = 10, 10
for z, m in zip(zs, mot):
    s = horner(A, mp.exp(z))
    mina = min(mina, abs(m))
    minarg = min(minarg, mp.pi - abs(mp.arg(m)))
    print(f"  arg z = {mp.nstr(mp.degrees(mp.arg(z)), 6):>9}°  f_K(e^z) = {mp.nstr(m, 8):>26}  |Reihe − Motor| = {mp.nstr(abs(s - m), 3)}")
print(f"  min |f_K(e^z)| = {mp.nstr(mina, 5)}, min Abstand arg zu ±π = {mp.nstr(minarg, 5)} rad   ({time.time() - t0:.1f} s)")

print("M3: f_K auf f_K(Γout): Motor sexp(slog(w)+1/2) gegen exp(h(w)), w = f_K(x), x = ρ e^{iφ}")
xs = [rho * mp.expj(mp.radians(ph)) for ph in [70, 50, 30, 10, -40]]
ws = [horner(A, x) for x in xs]
mot = fK_motor(ws)
for x, w, m in zip(xs, ws, mot):
    print(f"  φ = {mp.nstr(mp.degrees(mp.arg(x)), 4):>5}°  w = {mp.nstr(w, 8):>26}  |Motor − exp(x)| = {mp.nstr(abs(m - mp.exp(x)), 3)}")
# etwas ausserhalb von f_K(Γout): w = f_K(x·1.03)? nicht im Konvergenzbereich der Reihe; stattdessen w(1+0.05)
ws2 = [w * mp.mpf("1.05") for w in ws[:3]]
mot2 = fK_motor(ws2)
for w, m in zip(ws2, mot2):
    print(f"  w·1.05 = {mp.nstr(w, 8):>26}  f_K = {mp.nstr(m, 10)}  (endlich, Motor-Hauptzweig)")
print(f"Gesamt {time.time() - t0:.1f} s")
