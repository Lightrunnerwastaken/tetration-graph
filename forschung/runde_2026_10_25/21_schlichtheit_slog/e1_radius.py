"""E1: Schlichtheitsradius r_u(x0) = sup{ r <= x0+2 : sexp schlicht auf D(x0, r) } fuer reelle Zentren.

Weg A (direkt): Bild des Kreises |z-x0| = r doppelpunktfrei  <=>  schlicht auf der abgeschlossenen Scheibe
         (Darboux; sexp' != 0 nach Etappe 3, Satz A).
Weg B (Reduktion, Satz W): schlicht auf D(x0,r)  <=>  schlicht auf D(x0-1,r) und
         Bildkurve(x0-1,r) trifft Bildkurve(x0-1,r) + 2 pi i k fuer kein k != 0.
"""
import json
import sys
import numpy as np
import lib21 as S

N = int(sys.argv[1]) if len(sys.argv) > 1 else 2500


def kurve(x0, r, n=N):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return S.sexp(x0 + r * np.exp(1j * t))


def einfach_A(x0, r):
    return len(S.doppelpunkte(kurve(x0, r), True, 1)) == 0


def trifft_verschoben(w):
    wz = np.append(w, w[0])
    span = w.imag.max() - w.imag.min()
    kmax = int(span // (2 * np.pi)) + 1
    for k in range(1, kmax + 1):
        if S.schnitt_zweier(wz, wz + 2j * np.pi * k, 1):
            return True
    return False


def einfach_B(x0, r):
    # rekursiv: Scheibe um x0-1 schlicht (direkt pruefen, solange Werte klein) und keine 2 pi i-Paare
    w = kurve(x0 - 1, r)
    if not np.all(np.isfinite(w)):
        return False
    if len(S.doppelpunkte(w, True, 1)) > 0:
        return False
    return not trifft_verschoben(w)


def bisekt(test, x0, lo, hi, it=22):
    if test(x0, hi):
        return hi, True
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if test(x0, mid):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), False


def haupt():
  erg = {}
  for x0 in [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]:
    rmax = x0 + 2 - 0.004
    zeile = {}
    if x0 <= 2.5:
        r, voll = bisekt(einfach_A, x0, 1e-3, rmax)
        zeile["A"] = r; zeile["A_voll"] = voll
    if x0 >= 1.0:
        lo = 1e-12
        r, voll = bisekt(einfach_B, x0, lo, min(rmax, 3.0 if x0 < 2 else 1.5), it=40)
        zeile["B"] = r; zeile["B_voll"] = voll
    d1 = S.dsexp(np.array([x0 - 1.0 + 0j]), h=1e-6)[0].real if x0 - 1 > -2 else float("nan")
    zeile["sexp'(x0-1)"] = d1
    rr = zeile.get("B", zeile.get("A"))
    zeile["r*sexp'(x0-1)/pi"] = rr * d1 / np.pi
    erg[x0] = zeile
    print(x0, zeile, flush=True)

  with open("daten/e1_radius.json", "w", encoding="utf-8") as f:
    json.dump({str(k): v for k, v in erg.items()}, f, indent=1)


if __name__ == "__main__":
  haupt()
