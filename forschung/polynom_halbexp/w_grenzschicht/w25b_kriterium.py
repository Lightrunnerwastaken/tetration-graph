"""W25 — Entscheidungsrechnung: vorzeichenbehaftetes Randprofil auf ∂W̄ (D- und A-Bogen) und Poisson–Jensen-Kriterium
Θ = sup_{Schlitz ∪ Γ_d} ω_ζ(λ)/g(ζ) (C Lemma 13.2 / Satz 13.4, P-I(ii)) — aus den genauen Daten w24_N*.json.
Grössen: X ∈ {R_N (ganze kopffreie Antwort, einteilig), K_N (Eckteil, alle r(N) Moden), K_1 (nur u₁)}.
(i) endliches N: λ_N(η) = (1/N) ln(X(η)/X(L));  (ii) asymptotisch: Fit ln X = a + κ ln N + λN je Punkt, λ_rel = λ(η) − λ(L), Fehler aus Teilfenstern;
(iii) „termweise“ Variante: λ := max(λ, 0) auf dem D-Bogen, λ := 0 auf dem A-Bogen ((E_L⁺)-Weg ohne Auslöschung).
Aufruf: python w25_kriterium.py N1,N2,...  [Grösse: R | K | K1]
"""
import json
import sys

import numpy as np

import w24_rand_genau as M

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
GR = sys.argv[2] if len(sys.argv) > 2 else "K"
PRE = sys.argv[3] if len(sys.argv) > 3 else "w24"
if PRE == "w26":
    KO = json.load(open(M.W.GEMEINSAM + "/kneser_halbexp_koeff_4000.json"))["a"]
    M.a0 = float(KO[0]); M.R0 = abs(M.L - M.a0); M.psiL = np.angle(M.L - M.a0)
D = {N: json.load(open(M.W.HIER + f"/{PRE}_N{N}.json")) for N in NS}
names = list(D[NS[0]]["pts"].keys())
bd = [nm for nm in names if nm[0] in "DA"]           # Randpunkte (obere Hälfte)


def lnX(N, nm):
    p = D[N]["pts"][nm]
    if GR == "R":
        return p["lnR"]
    if GR == "K1":
        return p["lnK"][0]
    S = D[N]["S"]
    r = sum(1 for s in S[:len(p["lnK"])] if s > 1.40)
    r = max(r, 1)
    return 0.5 * np.log(sum(np.exp(2 * v) for v in p["lnK"][:r]))


# Randkoordinate s = τ(η) (reell) für Punkt und konjugierten Punkt
def s_of(z):
    t, _, _ = M.tau_map(z)
    return t.real


_, TAU0, _ = M.tau_map(3.0 + 0j)
pts_z = {nm: complex(*D[NS[0]]["pts"][nm]["z"]) for nm in names}
S_pts = []
for nm in bd:
    z = pts_z[nm]
    S_pts.append((nm, s_of(z), +1))
    S_pts.append((nm, s_of(np.conj(z)), -1))


def omega_over_g(zeta, lam_of_s):
    """ω_ζ(λ)/g(ζ) mit λ als Funktion von s (Poisson-Kern der oberen Halbebene), Integration in t = ln|s| beidseitig."""
    tz, t0, _ = M.tau_map(zeta)
    x, y = tz.real, tz.imag
    g = np.log(abs((tz - np.conj(t0)) / (tz - t0)))
    tt = np.linspace(-30, 30, 24001)
    tot = 0.0
    for sg in (+1, -1):
        s = sg * np.exp(tt)
        lam = lam_of_s(s)
        ker = y / ((s - x) ** 2 + y ** 2) / np.pi
        tot += np.trapezoid(lam * ker * np.abs(s), tt)
    return tot / g, g


def baue_profil(werte):
    """werte: dict nm → λ (obere Hälfte). Interpolation in u = sign(s)·ln|s| auf jeder Halbachse, λ → 0 zur Ecke (s → 0 bzw. ±∞)."""
    arr = sorted([(s, werte[nm]) for nm, s, _ in S_pts])
    sp = np.array([a for a, _ in arr])
    lp = np.array([b for _, b in arr])

    def f(s):
        out = np.zeros_like(s)
        for sg in (+1, -1):
            sel = np.sign(sp) == sg
            if not sel.any():
                continue
            xs = np.log(np.abs(sp[sel]))
            ls = lp[sel]
            o = np.argsort(xs)
            xs, ls = xs[o], ls[o]
            m = np.sign(s) == sg
            xv = np.log(np.abs(s[m]))
            v = np.interp(xv, xs, ls)
            # zur Ecke L (s→0, xv < xs[0]) und L̄ (s→±∞, xv > xs[-1]): linear in |s| bzw. 1/|s| auf 0
            lo = xv < xs[0]
            v[lo] = ls[0] * np.exp(xv[lo] - xs[0])
            hi = xv > xs[-1]
            v[hi] = ls[-1] * np.exp(xs[-1] - xv[hi])
            out[m] = v
        return out
    return f


def theta(werte, label):
    f = baue_profil(werte)
    us = [1e-8, 1e-6, 1e-4, 1e-3, 3e-3, 0.01, 0.03, 0.06, 0.1, 0.151]
    sl = [omega_over_g(M.schlitz_punkt(u), f)[0] for u in us]
    gd = {}
    for dd in (0.05, 0.1):
        gd[dd] = max(omega_over_g(z, f)[0] for z in M.gamma_d_punkte(dd, 37))
    print(f"  {label}: Schlitz ω/g bei u = " + " ".join(f"{u:g}:{v:+.3f}" for u, v in zip(us, sl))
          + f" | Γ_0.05: {gd[0.05]:+.3f}, Γ_0.1: {gd[0.1]:+.3f} | Θ = {max(max(sl), gd[0.05]):+.3f} (d=0.05), {max(max(sl), gd[0.1]):+.3f} (d=0.1)")
    return sl, gd


print(f"Grösse {GR}, N = {NS}")
# (i) endliches N
for N in NS:
    l0 = lnX(N, "L")
    w = {nm: (lnX(N, nm) - l0) / N for nm in bd}
    theta(w, f"N={N} endlich")
# (ii) asymptotisch
Ns = np.array(NS, float)
A = np.vstack([np.ones_like(Ns), np.log(Ns), Ns]).T


def fit(nm, sel=None):
    idx = range(len(NS)) if sel is None else sel
    yy = np.array([lnX(NS[i], nm) for i in idx])
    co, *_ = np.linalg.lstsq(A[list(idx)], yy, rcond=None)
    return co


lamL = fit("L")[2]
asym, err = {}, {}
half = len(NS) // 2
for nm in bd + [n for n in names if n[0] == "S"]:
    c = fit(nm)
    c1 = fit(nm, list(range(0, len(NS) - 1)))
    c2 = fit(nm, list(range(1, len(NS))))
    asym[nm] = c[2] - lamL
    err[nm] = max(abs(c1[2] - fit("L", list(range(0, len(NS) - 1)))[2] - asym[nm]),
                  abs(c2[2] - fit("L", list(range(1, len(NS))))[2] - asym[nm]))
print("\nAsymptotisches Profil λ(η) − λ(L) (Fit a + κ ln N + λN; ± = Unterschied der Teilfenster):")
print("  D-Bogen: " + " ".join(f"{nm[1:]}°:{asym[nm]:+.4f}±{err[nm]:.4f}" for nm in bd if nm[0] == "D"))
print("  A-Bogen: " + " ".join(f"{nm[1:]}:{asym[nm]:+.4f}±{err[nm]:.4f}" for nm in bd if nm[0] == "A"))
print("  Schlitz (gemessen, zum Vergleich −g): " + " ".join(f"u={nm[1:]}:{asym[nm]:+.4f}(g={M.green(pts_z[nm]):.4f})" for nm in names if nm[0] == "S"))
theta({nm: asym[nm] for nm in bd}, "asymptotisch")
theta({nm: asym[nm] + err[nm] for nm in bd}, "asympt. + Fehler")
# Exponent auf dem D-Bogen
dpos = [(float(nm[1:]), asym[nm]) for nm in bd if nm[0] == "D" and asym[nm] > 0]
print("  lokale Exponenten D-Bogen: " + " ".join(f"{t1}→{t2}:{np.log(l2/l1)/np.log(t2/t1):.2f}" for (t1, l1), (t2, l2) in zip(dpos[:-1], dpos[1:])))
# (iii) termweise: A-Bogen := 0, D-Bogen := max(λ, 0)
theta({nm: (max(asym[nm] + err[nm], 0) if nm[0] == "D" else 0.0) for nm in bd}, "termweise (λ_A := 0, D + Fehler)")
