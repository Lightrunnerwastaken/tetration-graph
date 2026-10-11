"""X19b — Gestalt der kohärenten Störung η_p = p_N − T_N f auf |x| = ρ: L¹-Masse relativ zu sup (= R), Breite des Buckels um L,
Mittelwert über Bögen der Länge τ um L (Grad). Aufruf: python x19b_masse.py N1 N2 ..."""
import json, sys
from pathlib import Path
import numpy as np
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
A = json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"]
PD = json.loads((HIER.parent / "gemeinsam" / "pN_N100_bis_1000.json").read_text())["p"]
L = complex(mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337)))
rho = abs(L); aL = np.angle(L)
for N in (int(a) for a in sys.argv[1:]):
    mp.mp.dps = 40
    e = [mp.mpf(a) - mp.mpf(b) for a, b in zip(PD[str(N)], A[: N + 1])]
    phis = aL + np.radians(np.linspace(-180, 180, 7201)[:-1])
    vals = []
    for ph in phis:
        z = mp.mpf(rho) * mp.expj(ph); s = mp.mpc(0)
        for c in reversed(e):
            s = s * z + c
        vals.append(complex(s))
    v = np.array(vals); R = np.abs(v).max()
    d = np.degrees(phis - aL)
    l1 = np.mean(np.abs(v)) / R
    mean0 = np.abs(np.mean(v)) / R
    stats = []
    for tau in (1, 3, 10, 30):
        m = np.abs(d) <= tau
        stats.append(f"τ={tau}°: |Mittel| {np.abs(np.mean(v[m])) / R:.3f}, Mittel|η| {np.mean(np.abs(v[m])) / R:.3f}")
    half = np.abs(d[np.abs(v) >= 0.5 * R])
    print(f"N = {N}: R = {R:.3e}; L¹-Masse/R = {l1:.4f} (Kreismittel |∫η|/R = {mean0:.2e}); |η| ≥ R/2 auf |φ − arg L| ≤ {half.max():.2f}° "
          f"(Anteil {np.mean(np.abs(v) >= 0.5 * R):.4f}); " + "; ".join(stats), flush=True)
