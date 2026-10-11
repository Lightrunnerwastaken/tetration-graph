"""X12b — Einteilige Poisson–Jensen-Bedingung im K'-Rahmen für den eckfreien Teil R̂ und den kanonischen Eckteil Rc (Kneser, x12-Daten).

λ(η) = (1/N) log(Q(η)/Q(L)) auf ∂K' für Q ∈ {R, R̂, Rc};  Θ = sup_{Schlitz ∪ Γ'_δ} ω'_ζ(λ)/g'(ζ)  (C Satz 13.4 / P-I(ii) mit K').
Halbebenen-Koordinate und Green-Funktion aus Ps Modul p_pruefung/kprime.py (nur gelesen/importiert, Logik wie Ps p15_theta.py).
Aufruf: python x12b_theta_eckteil.py N,d [N,d …]
"""
import json
import sys
from pathlib import Path

import numpy as np

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER.parent / "p_pruefung"))
import kprime as K  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
T0 = K.tau0
T2 = abs(T0) ** 2
U = np.linspace(np.log(1e-12), np.log(1e12), 60001)
SP = np.exp(U)
DS = SP * (U[1] - U[0])


def lam_daten(pts, key, N):
    QL = pts[0][key]
    neg, pos = [], []
    for p in pts[1:]:
        eta = complex(*p["eta"])
        lam = np.log(max(p[key], 1e-300) / QL) / N
        if p["rand"] == "in":
            neg.append((K.s_Gin(eta), lam))
        else:
            pos.append((K.s_Gout_x(K.Mw(eta).real), lam))
    neg = np.array(sorted(neg, key=lambda x: -x[0]))
    pos = np.array(sorted(pos, key=lambda x: x[0]))
    return QL, neg, pos


def lam_at(s, neg, pos):
    s = np.asarray(s, dtype=float)
    a = np.abs(s)
    a = np.where(a > np.sqrt(T2), T2 / a, a)
    out = np.zeros_like(s)
    for side, data in [(s < 0, neg), (s > 0, pos)]:
        ls, lv = np.log(np.abs(data[:, 0])), data[:, 1]
        o = np.argsort(ls)
        ls, lv = ls[o], lv[o]
        la = np.log(a[side])
        v = np.interp(la, ls, lv)
        small = la < ls[0]
        v[small] = lv[0] * np.exp(la[small] - ls[0])
        out[side] = v
    return out


def omega(tz, neg, pos):
    P = lambda s: (tz.imag / np.pi) / ((s - tz.real) ** 2 + tz.imag ** 2)
    return np.sum((lam_at(SP, neg, pos) * P(SP) + lam_at(-SP, neg, pos) * P(-SP)) * DS)


def gp(t):
    return np.log(abs((t - np.conj(T0)) / (t - T0)))


slit = []
for u in np.logspace(-5, np.log10(0.2043), 30):
    t = K.tau(K.L * (1 + u))
    slit.append((u, t, gp(t)))
for arg in sys.argv[1:]:
    pref = "x12"
    if ":" in arg:
        pref, arg = arg.split(":")
    N, d = (int(x) for x in arg.split(","))
    pts = json.loads((HIER / f"{pref}_N{N}_d{d}.json").read_text())["punkte"]
    for key in ("R", "Rhat", "Rc"):
        QL, neg, pos = lam_daten(pts, key, N)
        rs = [omega(t, neg, pos) / gv for (u, t, gv) in slit]
        ecke = (np.sum(lam_at(SP, neg, pos) / SP ** 2 * DS) + np.sum(lam_at(-SP, neg, pos) / SP ** 2 * DS)) / K.kappa0
        gam = []
        for dlev in (0.05, 0.1):
            kk = np.exp(dlev)
            c = (kk ** 2 * T0 - np.conj(T0)) / (kk ** 2 - 1)
            r = kk * abs(T0 - np.conj(T0)) / (kk ** 2 - 1)
            vals = [omega(c + r * np.exp(1j * ph), neg, pos) / dlev for ph in np.linspace(0, 2 * np.pi, 181)[:-1]
                    if (c + r * np.exp(1j * ph)).imag > 1e-9]
            gam.append(max(vals))
        print(f"N={N} d={d} {key:4s}: Q(L) = {QL:.4g}; λ max Γin {neg[:, 1].max():+.5f}, f(Γout) max {pos[:, 1].max():+.5f} min {pos[:, 1].min():+.4f};"
              f" Schlitz ω'/g' max {max(rs):+.4f} (u→0: {rs[0]:+.4f}); Ecke {ecke:+.4f}; Γ'_0.05 {gam[0]:+.4f}, Γ'_0.1 {gam[1]:+.4f}"
              f"  ⇒ Θ = {max(max(rs), gam[0]):+.4f}", flush=True)
