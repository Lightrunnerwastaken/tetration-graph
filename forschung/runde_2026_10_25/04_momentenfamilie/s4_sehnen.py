"""Waagrechte Sehnen der Kneser-Kurve Γ_K = θ((0,1)): Anzahl Schnittpunkte je Höhe und Konstante
m* = inf über Sehnen (t<s, gleiche Höhe) von  μ/d,  μ = Re θ(t) − Re θ(s),  d = min(s−t, 1−(s−t)).
Kriterium (BERICHT, Satz J): Lip(φ⁻¹ − id) < m*  ⇒  θ_T doppelpunktfrei."""
import json, sys
from pathlib import Path
H = Path(__file__).resolve().parent
d = json.loads((H / "daten" / sys.argv[1]).read_text())
X = d["x"]; P = [complex(a, b) for a, b in zip(d["re"], d["im"])]
ims = [p.imag for p in P]
lo, hi = min(ims), max(ims)
ny = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
anz = {}; best = (1e9, None)
for k in range(1, ny):
    y = lo + (hi - lo) * k / ny
    cr = []
    for i in range(len(P) - 1):
        a, b = P[i], P[i + 1]
        if (a.imag - y) * (b.imag - y) < 0:
            lam = (y - a.imag) / (b.imag - a.imag)
            cr.append((X[i] + lam * (X[i + 1] - X[i]), a.real + lam * (b.real - a.real)))
    anz[len(cr)] = anz.get(len(cr), 0) + 1
    for i in range(len(cr)):
        for j in range(i + 1, len(cr)):
            (t, rt), (s, rs) = cr[i], cr[j]
            mu = rt - rs; dd = min(s - t, 1 - (s - t))
            r = mu / dd
            if r < best[0]:
                best = (r, (y, t, s, mu, dd))
print(sys.argv[1], "Schnittzahlen je Höhe:", anz)
print("  m* =", f"{best[0]:.4f}", " bei Höhe/t/s/μ/d =", tuple(round(v, 5) for v in best[1]))
