"""E2: (a) groesstes Zentrum x_c mit schlichter Vollscheibe D(x0, x0+2);
(b) rho_pi(c) = Abstand von c zur Menge {|Im sexp| = pi}, Gegenprobe zu Satz W: r_u(c+1) = min(r_u(c), rho_pi(c));
(c) a_pi = sup{a : |Im sexp| < pi auf Re z <= a}: max_y Im sexp(a+iy).
"""
import json
import numpy as np
import lib21 as S
from e1_radius import einfach_A  # noqa  (laeuft e1 nicht erneut: Import nur der Funktion)

erg = {}
# (a)
for delta in (0.004, 0.0015):
    lo, hi = 0.0, 0.5
    for _ in range(18):
        mid = 0.5 * (lo + hi)
        if einfach_A(mid, mid + 2 - delta):
            lo = mid
        else:
            hi = mid
    erg[f"x_c(delta={delta})"] = 0.5 * (lo + hi)
    print("x_c", delta, 0.5 * (lo + hi), flush=True)


# (b) rho_pi(c): kleinster Radius, auf dem max |Im sexp| auf dem Kreis = pi erreicht
def maxim(c, r, n=1500):
    t = np.linspace(0, np.pi, n)
    return np.abs(S.sexp(c + r * np.exp(1j * t)).imag).max()


def rho_pi(c, rmax):
    if maxim(c, rmax) < np.pi:
        return rmax
    lo, hi = 1e-12, rmax
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        if maxim(c, mid) < np.pi:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


for c in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
    rp = rho_pi(c, min(c + 2 - 0.004, 3.0))
    erg[f"rho_pi({c})"] = rp
    print("rho_pi", c, rp, flush=True)

# (c) max_y Im sexp(a+iy), y in (0, 12]
ys = np.concatenate([np.linspace(1e-4, 0.2, 400), np.linspace(0.2, 12, 1200)])
tab = []
for a in np.arange(-3.0, 2.01, 0.1):
    a = round(a, 3)
    v = S.sexp(a + 1j * ys)
    im = v.imag
    k = int(np.argmax(np.abs(im)))
    tab.append((a, float(np.abs(im).max()), float(ys[k])))
    print(f"a={a:5.2f} max|Im sexp| = {np.abs(im).max():.6f} bei y={ys[k]:.4f}", flush=True)
erg["max_im_tabelle"] = tab
# Nullstelle von max|Im| - pi in a
f = lambda a: np.abs(S.sexp(a + 1j * ys).imag).max() - np.pi
lo, hi = None, None
for a, m, _ in tab:
    if m < np.pi:
        lo = a
    elif lo is not None and hi is None:
        hi = a
if lo is not None and hi is not None:
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    erg["a_pi"] = 0.5 * (lo + hi)
    print("a_pi =", 0.5 * (lo + hi))
json.dump(erg, open("daten/e2_pi_band.json", "w", encoding="utf-8"), indent=1)
