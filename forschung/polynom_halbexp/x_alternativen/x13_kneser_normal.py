"""X13 — wie x12 (Kneser oder Fast-Lösung P = p_N auf Ps ∂K'-Gitter), aber Kleinste Quadrate über die NORMALGLEICHUNGEN (für grosse d, z. B. d = 4N),
plus σ_min des Rechteck-Abschnitts (a″)_d.

Operator: J_∞(F)v = v∘F + m_F·v mit F = f_K (m = e^x/f_K' = f_K'∘f_K) oder F = p_N (Tangente an der echten Lösung: m = p_N'∘p_N).
R = T_{N+d}J_∞(F)|_{P_N} (ℓ²_ρ-skaliert), Π̂ = (RᵀR)⁻¹Rᵀ T_{N+d}J_∞(F), Π_N = J_N⁻¹T_N J_∞(F); Rc = |Π_Nk_η − Π̂k_η|₁, R̂ = |Π̂k_η − T_Nk_η|₁.
Ausgabe x13_{F}_N{N}_d{d}.json (Format wie x12, für x12b-artige Θ-Auswertung).
Aufruf: python x13_kneser_normal.py N d [kneser|pN] [prec-Faktor]
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from flint import acb, acb_mat, acb_poly, arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
N, d = int(sys.argv[1]), int(sys.argv[2])
ART = sys.argv[3] if len(sys.argv) > 3 else "kneser"
ctx.prec = int(float(sys.argv[4]) * N) + 400 if len(sys.argv) > 4 else int(1.8 * N) + 400
M = N + d
mp.mp.dps = 60
A_str = json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"]
A_mp = [mp.mpf(s) for s in A_str]
Lm = mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))
rho_m = abs(Lm)
argL = mp.arg(Lm)


def fK_mp(z):
    s = mp.mpc(0)
    for c in reversed(A_mp):
        s = s * z + c
    return s


def to_acb(z):
    return acb(arb(mp.nstr(mp.re(z), 55)), arb(mp.nstr(mp.im(z), 55)))


def randpunkte():
    taus = sorted(set(np.round(np.concatenate([np.logspace(-3, np.log10(103.38), 60), [103.38]]), 6)))
    ts = sorted(set(np.round(np.concatenate([np.logspace(-3, np.log10(76.618), 50), [76.618]]), 6)))
    pts = [("L", 0.0, Lm)]
    for t in taus:
        pts.append(("in", float(t), rho_m * mp.expj(argL + mp.radians(t))))
    for t in ts:
        pts.append(("out", float(t), fK_mp(rho_m * mp.expj(argL - mp.radians(t)))))
    return pts


def mullow(a, b, n):
    return (a * b).truncate(n)


def inv_series(a, n, P):
    g = P([1 / a[0]])
    k = 1
    while k < n:
        k = min(2 * k, n)
        g = mullow(g, P([2]) - mullow(a.truncate(k), g, k), k)
    return g


t0 = time.time()
if ART == "kneser":
    F = arb_poly([arb(s) for s in A_str[: M + 1]])
    Fp = F.derivative()
    ex = arb_poly([1 / arb.fac_ui(k) for k in range(M + 1)])
    m = mullow(ex, inv_series(Fp, M + 1, arb_poly), M + 1)
else:
    pd = json.loads((HIER.parent / "gemeinsam" / "pN_N100_bis_1000.json").read_text())["p"][str(N)]
    F = arb_poly([arb(s) for s in pd])
    Fp = F.derivative()
    # m = p'∘p, abgeschnitten bei M (Horner mit Abschneiden)
    m = arb_poly([0])
    for c in reversed([Fp[k] for k in range(Fp.degree() + 1)]):
        m = mullow(m, F, M + 1) + arb_poly([c])
rho = arb(mp.nstr(rho_m, 55))
rk = [arb(1)]
for k in range(1, M + 1):
    rk.append(rk[-1] * rho)
A = arb_mat(M + 1, N + 1)
pw = arb_poly([1])
for j in range(N + 1):
    for k in range(M + 1):
        A[k, j] = (pw[k] + (m[k - j] if k >= j else 0)) * rk[k] / rk[j]
    pw = mullow(pw, F, M + 1)
JN = arb_mat(N + 1, N + 1)
for i in range(N + 1):
    for j in range(N + 1):
        JN[i, j] = A[i, j]
G = A.transpose() * A                           # RᵀR
Gi = G.inv()
Gf = np.array([[float(Gi[i, j].mid()) for j in range(N + 1)] for i in range(N + 1)])
ew = np.sort(np.linalg.eigvalsh((Gf + Gf.T) / 2))[::-1]
smin = 1 / np.sqrt(ew[:3])
print(f"{ART} N={N} d={d} prec={ctx.prec}: σ_min(R) = {smin[0]:.5f}, {smin[1]:.5f}, {smin[2]:.5f}  ({time.time() - t0:.0f}s)", flush=True)

PTS = randpunkte()
npk = len(PTS)
Fa = acb_poly([acb(F[k]) for k in range(min(F.degree(), M) + 1)])
ma = acb_poly([acb(m[k]) for k in range(m.degree() + 1)])
W = acb_mat(M + 1, npk)
KZ = []
for p, (name, par, eta) in enumerate(PTS):
    e = to_acb(eta)
    inv_f = inv_series(acb_poly([e]) - Fa, M + 1, acb_poly)
    kz = acb_poly([1 / e ** (k + 1) for k in range(M + 1)])
    w = inv_f + mullow(ma, kz, M + 1)
    for k in range(M + 1):
        W[k, p] = w[k] * rk[k]
    KZ.append([kz[k] * rk[k] for k in range(N + 1)])
Wt = acb_mat(N + 1, npk)
for k in range(N + 1):
    for p in range(npk):
        Wt[k, p] = W[k, p]
X = acb_mat(JN).solve(Wt)
Xh = acb_mat(Gi) * (acb_mat(A.transpose()) * W)
print(f"   Lösungen ({time.time() - t0:.0f}s)", flush=True)
out = []
for p, (name, par, eta) in enumerate(PTS):
    R = Rh = Rc = 0.0
    for k in range(N + 1):
        y = X[k, p] - KZ[p][k]
        cc = X[k, p] - Xh[k, p]
        R += abs(complex(y.mid()))
        Rc += abs(complex(cc.mid()))
        Rh += abs(complex((Xh[k, p] - KZ[p][k]).mid()))
    out.append({"rand": name, "par": par, "eta": [float(mp.re(eta)), float(mp.im(eta))], "R": R, "Rhat": Rh, "Rc": Rc})
(HIER / f"x13_{ART}_N{N}_d{d}.json").write_text(json.dumps({"N": N, "d": d, "art": ART, "smin": list(smin), "punkte": out}), encoding="utf-8")
L0 = out[0]
far = max((o for o in out if o['rand'] == 'in'), key=lambda o: o['par'])
print(f"   L: R = {L0['R']:.4f}, R̂ = {L0['Rhat']:.4f}, Rc = {L0['Rc']:.4f};  −ρ: R = {far['R']:.4f}, R̂ = {far['Rhat']:.4f}, Rc = {far['Rc']:.4f}"
      f"  [{time.time() - t0:.0f}s]", flush=True)
