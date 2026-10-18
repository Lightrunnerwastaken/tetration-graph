"""Rang 5: in welchem Becken liegt H5(z)?  Raster z = x + iy, x in [0,1), y in [0.02, 0.5]; Vorwaerts-Iteration von sexp (H4).
Ergebnis: Grenze des Beckens von L4 (dort endet die Koenigs-Darstellung von theta5). Ausgabe daten/becken5.json + Text-Karte.
Aufruf: python s4_becken5.py NX NY"""
import sys, json, warnings, numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER.parents[1] / "eta_hyperoperationen" / "halb_rang_basis2_2026_10_13"))
import e6_lade
from e6_superfunktion import fixpunkt
R = e6_lade.lade_raenge(5); H4, H5 = R[4], R[5]
L = fixpunkt(H4, 1.47 + 1.66j)
NX, NY = int(sys.argv[1]), int(sys.argv[2])
xs = np.arange(NX) / NX; ys = np.linspace(0.5, 0.02, NY)
Z = (xs[None, :] + 1j * ys[:, None]).ravel()
w = H5(Z); klasse = np.zeros(len(w), int)   # 1: L4, 2: conj(L4), 3: entkommt, 0: offen
for k in range(150):
    m = klasse == 0
    if not m.any(): break
    w[m] = H4(w[m])
    klasse[m & (np.abs(w - L) < 1e-6)] = 1
    klasse[m & (np.abs(w - np.conj(L)) < 1e-6)] = 2
    klasse[m & ~(np.abs(w) < 60)] = 3
K = klasse.reshape(NY, NX)
zeichen = {0: "?", 1: ".", 2: "o", 3: "#"}
for i in range(NY):
    print("%5.3f " % ys[i] + "".join(zeichen[v] for v in K[i]))
# hoechster Punkt, der nicht im Becken von L4 liegt
nicht = [(ys[i], xs[j]) for i in range(NY) for j in range(NX) if K[i, j] != 1]
print("hoechster Nicht-L4-Punkt:", max(nicht) if nicht else None)
json.dump({"x": xs.tolist(), "y": ys.tolist(), "klasse": K.tolist()}, open(HIER / "daten" / "becken5.json", "w"))
