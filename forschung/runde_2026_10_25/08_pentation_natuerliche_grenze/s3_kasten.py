"""Schritt 3: Kastenzaehl-Dimension des Beckenrands (Klasse 1 gegen Rest) im Fenster [X0, X0+W] x [Y0, Y0+H], Raster eps.
Kasten der Kantenlaenge k*eps zaehlt, wenn er Punkte beider Sorten enthaelt; Steigung von log N gegen log(1/(k eps)).
Aufruf: python s3_kasten.py X0 Y0 W H eps"""
import sys, json, time
import numpy as np
from lib import H5, klasse, HIER
t0 = time.time()
X0, Y0, W, Hh, eps = map(float, sys.argv[1:6])
nx, ny = int(round(W / eps)), int(round(Hh / eps))
xs = X0 + eps * np.arange(nx); ys = Y0 + eps * np.arange(ny)
Z = (xs[None, :] + 1j * ys[:, None]).ravel()
kl, n = klasse(H5(Z), nmax=120)
B = (kl.reshape(ny, nx) == 1)
print("Punkte %d, Klassen %s, unentschieden %d" % (nx * ny, np.bincount(kl, minlength=4).tolist(), int((kl == 0).sum())))
erg = []
for k in (1, 2, 4, 8, 16):
    mx, my = nx // (k + 1), ny // (k + 1)
    cnt = 0
    for i in range(my):
        for j in range(mx):
            blk = B[i * k:(i * k) + k + 1, j * k:(j * k) + k + 1]   # Kasten mit Ecken (k+1 Punkte je Kante)
            if blk.any() and not blk.all(): cnt += 1
    # auf gleiche Flaeche normieren: Kastenzahl pro Fenster
    erg.append((k * eps, cnt * (nx * ny) / (mx * my * (k + 1) ** 2) * k * k / 1.0))
    print("Kasten %.4f: N = %d (von %d)" % (k * eps, cnt, mx * my))
# einfache Zaehlung ohne Ueberlappung: Kaesten k x k Punkte, Grenze zwischen Nachbarn mitzaehlen
erg2 = []
for k in (1, 2, 4, 8, 16):
    mx, my = (nx - 1) // k, (ny - 1) // k
    cnt = 0
    for i in range(my):
        for j in range(mx):
            blk = B[i * k:i * k + k + 1, j * k:j * k + k + 1]
            if blk.any() and not blk.all(): cnt += 1
    erg2.append((k * eps, cnt))
le = np.log([1 / e for e, _ in erg2]); lN = np.log([max(c, 1) for _, c in erg2])
print("Kasten (Seite, N):", erg2)
print("lokale Steigungen:", " ".join("%.3f" % v for v in np.diff(lN) / np.diff(le)), "  Gesamtsteigung %.3f" % np.polyfit(le, lN, 1)[0])
json.dump({"fenster": [X0, Y0, W, Hh, eps], "kaesten": erg2}, open(HIER / "cache" / "s3_kasten.json", "w"))
print("Zeit %.1f s" % (time.time() - t0))
