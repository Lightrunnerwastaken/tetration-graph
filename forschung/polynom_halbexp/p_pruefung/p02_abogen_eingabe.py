"""Strang P, Prüfung 2: Eingabe w_eta = [J_inf k_eta]_{<=N} für Pole eta auf dem A-Bogen.

Behauptung ARCHITEKTUR §3/§5.2: (a') ⇒ (E_R) über Satz D8. Dafür muss |w_eta|_{l2_rho} auf ganz dW polynomial sein.
Hier: Norm von [1/(eta − f_K(x))]_{<=N} und von [m(x)/(eta − x)]_{<=N} auf dem A-Bogen und dem D-Bogen.
"""
import json
import sys
from pathlib import Path

import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 60
GEM = Path(__file__).resolve().parent.parent / "gemeinsam"
a = [mp.mpf(s) for s in json.loads((GEM / "kneser_halbexp_koeff_4000.json").read_text())["a"]]
L = mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))
rho, a0 = abs(L), a[0]
R0 = abs(L - a0)


def inv_reihe(eta, N):
    """Koeffizienten von 1/(eta − f(x)) bis x^N."""
    c = eta - a0
    b = [1 / c]
    for k in range(1, N + 1):
        s = mp.mpc(0)
        for j in range(1, k + 1):
            s += a[j] * b[k - j]
        b.append(s / c)
    return b


def l2(b):
    return mp.sqrt(sum(abs(x) ** 2 * rho ** (2 * k) for k, x in enumerate(b)))


faelle = [("A-Scheitel a0+R0", a0 + R0),
          ("A-Bogen psi=60°", a0 + R0 * mp.expj(mp.radians(60))),
          ("A-Bogen psi=90°", a0 + R0 * mp.expj(mp.radians(90))),
          ("D-Bogen 120°", rho * mp.expj(mp.radians(120))),
          ("D-Bogen 180°", -rho)]
Ns = [50, 100, 200, 300]
print("|[1/(eta − f(x))]_{<=N}|_{l2_rho}:")
print("Fall".ljust(20) + "".join(f"N={N:<10d}" for N in Ns) + "Rate/Grad (200→300)")
for name, eta in faelle:
    b = inv_reihe(eta, max(Ns))
    v = [l2(b[:N + 1]) for N in Ns]
    rate = (v[3] / v[2]) ** (mp.mpf(1) / 100)
    print(name.ljust(20) + "".join(f"{mp.nstr(x, 4):<12}" for x in v) + mp.nstr(rate, 6))
