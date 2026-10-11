"""X15 — Asymptotisches Θ^c_∞ für den kanonischen Eckteil (Kneser, kleines d): Fit je Randpunkt
  log(Rc_N(η)/Rc_N(L)) = a(η) + b(η) log N + λ(η) N   über N = 200, 300, 400, 600 (x12-Daten),
dann Θ_∞ = sup_{Schlitz ∪ Γ'_δ} ω'(λ)/g' mit dem reinen Ratenprofil λ (polynomiale Faktoren gehen in κ).
Zusätzlich robust: λ_max(η) := max(λ, λ aus dem Fit ohne log-Term) und Fehlerabschätzung über Teilfenster (200–400 gegen 300–600).
Logik und Gitter wie x12b (Ps kprime).
"""
import json
import sys
from pathlib import Path

import numpy as np

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER.parent / "p_pruefung"))
import kprime as K  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
import x12b_theta_eckteil as T  # noqa: E402  (nutzt dessen lam_at/omega/gp; x12b führt beim Import nichts aus, wenn ohne Argumente)

RUNS = [(200, 3), (300, 4), (400, 6), (600, 12)]
daten = {N: json.loads((HIER / f"x12_N{N}_d{d}.json").read_text())["punkte"] for N, d in RUNS}
Ns = np.array([N for N, _ in RUNS], dtype=float)
npkt = len(daten[200])


def fit(idx, Nsel):
    y = np.array([np.log(daten[N][idx]["Rc"] / daten[N][0]["Rc"]) for N in Nsel])
    X = np.vstack([np.ones(len(Nsel)), np.log(Nsel), Nsel]).T
    if len(Nsel) >= 3:
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        return coef[2]
    return (y[-1] - y[0]) / (Nsel[-1] - Nsel[0])


def profil(Nsel):
    pts = []
    for i in range(1, npkt):
        p = daten[200][i]
        pts.append({"rand": p["rand"], "par": p["par"], "eta": p["eta"], "lam": fit(i, Nsel)})
    return pts


def theta(pts):
    neg, pos = [], []
    for p in pts:
        eta = complex(*p["eta"])
        if p["rand"] == "in":
            neg.append((K.s_Gin(eta), p["lam"]))
        else:
            pos.append((K.s_Gout_x(K.Mw(eta).real), p["lam"]))
    neg = np.array(sorted(neg, key=lambda x: -x[0]))
    pos = np.array(sorted(pos, key=lambda x: x[0]))
    rs = [T.omega(t, neg, pos) / gv for (u, t, gv) in T.slit]
    gam = {}
    for dlev in (0.05, 0.1):
        kk = np.exp(dlev)
        c = (kk ** 2 * T.T0 - np.conj(T.T0)) / (kk ** 2 - 1)
        r = kk * abs(T.T0 - np.conj(T.T0)) / (kk ** 2 - 1)
        vals = [T.omega(c + r * np.exp(1j * ph), neg, pos) / dlev for ph in np.linspace(0, 2 * np.pi, 181)[:-1]
                if (c + r * np.exp(1j * ph)).imag > 1e-9]
        gam[dlev] = max(vals)
    ecke = (np.sum(T.lam_at(T.SP, neg, pos) / T.SP ** 2 * T.DS) + np.sum(T.lam_at(-T.SP, neg, pos) / T.SP ** 2 * T.DS)) / K.kappa0
    return max(rs), gam, ecke, neg, pos


for name, Nsel in (("Fit N = 200–600", Ns), ("Fenster 200–400", Ns[:3]), ("Fenster 300–600", Ns[1:])):
    pts = profil(Nsel)
    smax, gam, ecke, neg, pos = theta(pts)
    print(f"{name}: λ Γin bei τ = 0.01/0.1/1/10/103°: "
          + ", ".join(f"{np.interp(np.log(v), np.log([p['par'] for p in pts if p['rand'] == 'in']), [p['lam'] for p in pts if p['rand'] == 'in']):+.5f}"
                      for v in (0.01, 0.1, 1, 10, 103.38))
          + "; f(Γout) bei t = 0.01/0.1/1/10°: "
          + ", ".join(f"{np.interp(np.log(v), np.log([p['par'] for p in pts if p['rand'] == 'out']), [p['lam'] for p in pts if p['rand'] == 'out']):+.5f}"
                      for v in (0.01, 0.1, 1, 10)))
    print(f"   Θ_∞: Schlitz max {smax:+.4f}; Γ'_0.05 {gam[0.05]:+.4f}; Γ'_0.1 {gam[0.1]:+.4f}; Ecken-Grenzwert {ecke:+.3f}"
          f"  ⇒ Θ_∞(δ=0.05) = {max(smax, gam[0.05]):+.4f}, Θ_∞(δ=0.1) = {max(smax, gam[0.1]):+.4f}")
