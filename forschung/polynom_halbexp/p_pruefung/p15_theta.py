"""Strang P, Prüfung 15: einteilige Bedingung im K'-Rahmen (C Satz 13.4 mit K' statt W̄):
  λ_N(η) = (1/N) log(R_N(η)/R_N(L)) auf ∂K',   Θ'_N = sup_{ζ ∈ S_L(d) ∪ Γ'_d} ω'_ζ(λ_N)/g'(ζ)  ≤ 1 − ε ?
Randdaten aus p14_RN_N{N}.json, Halbebenen-Koordinate aus kprime (E' → obere Halbebene, Γin ↦ s<0, f_K(Γout) ↦ s>0).
Aufruf: python p15_theta.py 100,200,300
"""
import json
import sys
from pathlib import Path

import numpy as np

import kprime as K

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
T0 = K.tau0
T2 = abs(T0) ** 2


def lam_funktion(N):
    d = json.loads((HIER / f"p14_RN_N{N}.json").read_text())["punkte"]
    RL = d[0]["R"]
    neg, pos = [], []
    for p in d[1:]:
        eta = complex(*p["eta"])
        lam = np.log(p["R"] / RL) / N
        if p["rand"] == "in":
            neg.append((K.s_Gin(eta), lam))
        else:
            pos.append((K.s_Gout_x(K.Mw(eta).real), lam))
    neg = np.array(sorted(neg, key=lambda x: -x[0]))        # s von 0⁻ nach −|τ₀|
    pos = np.array(sorted(pos, key=lambda x: x[0]))
    return RL, neg, pos


def lam_at(s, neg, pos):
    """λ(s) für s ∈ ℝ∖{0}; Spiegelung s ↦ |τ₀|²/s (Konjugation) für |s| > |τ₀|; Interpolation in log|s|."""
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
        v[small] = lv[0] * np.exp(la[small] - ls[0])       # Fortsetzung → 0 linear in |s|
        out[side] = v
    return out


# Integrationsgitter in s (beidseitig, log-verteilt)
U = np.linspace(np.log(1e-12), np.log(1e12), 60001)
SP = np.exp(U)
DS = SP * (U[1] - U[0])


def omega(tz, neg, pos):
    P = lambda s: (tz.imag / np.pi) / ((s - tz.real) ** 2 + tz.imag ** 2)
    lp = lam_at(SP, neg, pos)
    ln = lam_at(-SP, neg, pos)
    return np.sum((lp * P(SP) + ln * P(-SP)) * DS)


def gprime_tau(t):
    return np.log(abs((t - np.conj(T0)) / (t - T0)))


Ns = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [100]
# Schlitzpunkte im τ'-Bild (einmal)
slit = []
for u in np.logspace(-5, np.log10(0.2043), 30):
    t = K.tau(K.L * (1 + u))
    slit.append((u, t, gprime_tau(t)))
for N in Ns:
    RL, neg, pos = lam_funktion(N)
    print(f"N = {N}: R_N(L) = {RL:.3f}; λ_N auf Γin: max {neg[:, 1].max():+.5f} (bei s = {neg[np.argmax(neg[:, 1]), 0]:.3e}),"
          f" auf f(Γout): max {pos[:, 1].max():+.5f}, min {pos[:, 1].min():+.4f}")
    ecke = (np.sum(lam_at(SP, neg, pos) / SP ** 2 * DS) + np.sum(lam_at(-SP, neg, pos) / SP ** 2 * DS)) / K.kappa0
    rs = [omega(t, neg, pos) / gv for (u, t, gv) in slit]
    print(f"   Schlitz: ω'(λ)/g' = " + ", ".join(f"{r:+.4f}" for r in rs[::5]) + f"  (u = 1e-5 … 0.204);  Ecken-Grenzwert (1/κ₀')∫λ/s² = {ecke:+.4f}")
    for dlev in [0.05, 0.1]:
        kk = np.exp(dlev)
        c = (kk ** 2 * T0 - np.conj(T0)) / (kk ** 2 - 1)
        r = kk * abs(T0 - np.conj(T0)) / (kk ** 2 - 1)
        vals = []
        for ph in np.linspace(0, 2 * np.pi, 361)[:-1]:
            t = c + r * np.exp(1j * ph)
            if t.imag <= 1e-9:
                continue
            vals.append(omega(t, neg, pos) / dlev)
        print(f"   Γ'_d (d = {dlev}): max ω'(λ)/d = {max(vals):+.4f}, min = {min(vals):+.4f}")
    print(f"   ⇒ Θ'_N ≈ {max(max(rs), max(vals)):+.4f}")
