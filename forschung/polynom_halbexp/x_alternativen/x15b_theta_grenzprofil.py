"""X15b — Θ_∞ im K'-Rahmen für ein KONSERVATIVES Grenzprofil des Eckteils:
Γin-Seite: Ws asymptotisches Eckteil-Profil (BRETT W 01:06, W̄-Messung; Ratenzuwachs gegenüber L): τ = 0/2/4/8/30/103.38° → 0/0.0009/0.0024/0.0047/0.0090/0.0115,
dazu Aufschlag +20 % (Fehler); f_K(Γout)-Seite: (a) λ = 0 (konservativ), (b) gefittetes Profil aus x15 (−0.0047 bei 1°, −0.053 bei 10°).
"""
import json, sys
from pathlib import Path
import numpy as np
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER.parent / "p_pruefung"))
import kprime as K
import x12b_theta_eckteil as T
sys.stdout.reconfigure(encoding="utf-8")
pts = json.loads((HIER / "x12_N200_d3.json").read_text())["punkte"][1:]
tw = np.array([0, 2, 4, 8, 30, 103.38]); lw = 1.2 * np.array([0, 0.0009, 0.0024, 0.0047, 0.0090, 0.0115])
to = np.array([0.001, 0.01, 0.1, 1, 10, 76.618]); lo = np.array([-0.000005, -0.00005, -0.00044, -0.00475, -0.0527, -0.25])
for variante in ("aussen 0", "aussen Fit"):
    neg, pos = [], []
    for p in pts:
        eta = complex(*p["eta"])
        if p["rand"] == "in":
            neg.append((K.s_Gin(eta), float(np.interp(p["par"], tw, lw))))
        else:
            lam = 0.0 if variante == "aussen 0" else float(np.interp(np.log(p["par"]), np.log(to), lo))
            pos.append((K.s_Gout_x(K.Mw(eta).real), lam))
    neg = np.array(sorted(neg, key=lambda x: -x[0])); pos = np.array(sorted(pos, key=lambda x: x[0]))
    rs = [T.omega(t, neg, pos) / gv for (u, t, gv) in T.slit]
    out = []
    for dlev in (0.05, 0.1):
        kk = np.exp(dlev); c = (kk**2 * T.T0 - np.conj(T.T0)) / (kk**2 - 1); r = kk * abs(T.T0 - np.conj(T.T0)) / (kk**2 - 1)
        out.append(max(T.omega(c + r*np.exp(1j*ph), neg, pos) / dlev for ph in np.linspace(0, 2*np.pi, 181)[:-1] if (c + r*np.exp(1j*ph)).imag > 1e-9))
    ecke = (np.sum(T.lam_at(T.SP, neg, pos) / T.SP**2 * T.DS) + np.sum(T.lam_at(-T.SP, neg, pos) / T.SP**2 * T.DS)) / K.kappa0
    print(f"{variante}: Schlitz max {max(rs):+.4f} (u→0 {rs[0]:+.4f}); Γ'_0.05 {out[0]:+.4f}; Γ'_0.1 {out[1]:+.4f}; Ecke {ecke:+.3f}"
          f"  ⇒ Θ_∞(0.05) = {max(max(rs), out[0]):+.4f}, Θ_∞(0.1) = {max(max(rs), out[1]):+.4f}")
