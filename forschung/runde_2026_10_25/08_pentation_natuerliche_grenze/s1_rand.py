"""Schritt 1: Grundpruefungen und oberster Punkt y_max des Randes von Omega = Komponente von {z: H5(z) im Becken von L4}, die i*inf enthaelt.
(a) L4, mu; Bahn von L3 (Fixpunkt von 2^z) unter sexp.
(b) Funktionalgleichung H5(z+1) = H4(H5(z)) auf einem Raster nahe dem Rand (y in [0.02, 0.5]).
(c) y*(x) = hoechstes y mit Klasse != 1 auf senkrechten Linien (grob, dann fein um das Maximum).
Aufruf: python s1_rand.py NX DY"""
import sys, json, time
import numpy as np
from lib import H4, H5, L, L3, mu, klasse, HIER
t0 = time.time()
NX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
DY = float(sys.argv[2]) if len(sys.argv) > 2 else 0.005
print("L4 =", L, " mu =", mu, " |mu| =", abs(mu), " L3 =", L3)
w = np.array([L3]); bahn = []
for k in range(12):
    w = H4(w); bahn.append(w[0])
print("Bahn von L3 unter sexp:", " ".join("%.4f%+.4fi" % (v.real, v.imag) for v in bahn[:6]), "...", bahn[-1])
# (b)
xs = np.linspace(0, 1, 41)[:-1]; ys = np.linspace(0.02, 0.5, 25)
Z = (xs[None, :] + 1j * ys[:, None]).ravel()
fe = np.abs(H5(Z + 1) - H4(H5(Z)))
print("Funktionalgleichung |H5(z+1)-H4(H5 z)|: max %.2e, median %.2e;  min Abstand H5(z) zum Schnitt (-inf,-2]: %.3f" %
      (fe.max(), np.median(fe), np.min(np.where(np.real(H5(Z)) < -2, np.abs(np.imag(H5(Z))), np.abs(H5(Z) + 2)))))
# (c) grob
xg = np.arange(NX) / NX
ygr = np.arange(0.40, 0.15, -DY)
Z = (xg[None, :] + 1j * ygr[:, None]).ravel()
kl, n = klasse(H5(Z), nmax=150)
KL = kl.reshape(len(ygr), NX)
print("Klassen grob:", {k: int((kl == k).sum()) for k in range(4)}, " max Schritte im Becken:", int(n[kl == 1].max()))
ystern = np.array([ygr[np.argmax(KL[:, j] != 1)] if (KL[:, j] != 1).any() else np.nan for j in range(NX)])
jm = np.nanargmax(ystern)
print("grob: y_max ~ %.4f bei x = %.4f" % (ystern[jm], xg[jm]))
top = np.argsort(-np.nan_to_num(ystern, nan=-1))[:8]
print("  hoechste Stellen:", ", ".join("x=%.3f y=%.3f" % (xg[j], ystern[j]) for j in top))
# fein um die drei hoechsten getrennten Stellen
erg = {"grob": {"x": xg.tolist(), "ystern": [None if np.isnan(v) else float(v) for v in ystern]}}
gesehen = []
fein = []
for j in top:
    if any(abs(xg[j] - g) < 0.03 for g in gesehen): continue
    gesehen.append(xg[j])
    x0, y0 = xg[j], ystern[j]
    xf = x0 + np.linspace(-1.5 / NX, 1.5 / NX, 41)
    yf = np.arange(y0 + 1.2 * DY, y0 - 0.2 * DY, -DY / 20)
    Z = (xf[None, :] + 1j * yf[:, None]).ravel()
    kl, n = klasse(H5(Z), nmax=250)
    KF = kl.reshape(len(yf), len(xf))
    ys2 = np.array([yf[np.argmax(KF[:, i] != 1)] if (KF[:, i] != 1).any() else np.nan for i in range(len(xf))])
    i = np.nanargmax(ys2)
    fein.append((float(xf[i]), float(ys2[i]), int((kl == 0).sum())))
    print("fein um x=%.3f: y_max = %.5f bei x = %.5f (Raster %.1e x %.1e), unentschieden: %d" % (x0, ys2[i], xf[i], xf[1] - xf[0], DY / 20, (kl == 0).sum()))
    if len(gesehen) >= 2: break
erg["fein"] = fein
json.dump(erg, open(HIER / "cache" / "s1_rand.json", "w"))
print("Zeit %.1f s" % (time.time() - t0))
