"""X12 — KNESER: Pol-Antwort R, eckfreier Teil R̂ (Kleinste Quadrate, d Zusatzzeilen) und kanonischer Eckteil Rc auf ∂K' (Ps Gitter, p14).

J̃ = ℓ²_ρ-skaliertes J_∞(f_K) (Kneser-Koeffizienten der 4000er-Datei als exakte Daten), m = e^x/f_K' = f_K'∘f_K, Zeilen bis N + d.
Formeln (BERICHT §3.5/3.6): Π̂g = Π_N g − c(g), c(g) = J_N⁻¹Zᵀ(I + ZZᵀ)⁻¹o(g), Z = BJ_N⁻¹, o(g) = BΠ_N g − v(g).
Ausgabe: x12_N{N}_d{d}.json im Format von Ps p14 (Felder rand, par, eta, R, Rhat, Rc) für die Θ-Auswertung (x12b).
Aufruf: python x12_kneser_kprime.py N d
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
M = N + d
ctx.prec = int(float(sys.argv[3]) * N) + 400 if len(sys.argv) > 3 else int(1.8 * N) + 400
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


def randpunkte():   # identisch zu Ps p14
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
f = arb_poly([arb(s) for s in A_str[: M + 1]])
fp = f.derivative()
ex = arb_poly([1 / arb.fac_ui(k) for k in range(M + 1)])
m = mullow(ex, inv_series(fp, M + 1, arb_poly), M + 1)
rho = arb(mp.nstr(rho_m, 55))
rk = [arb(1)]
for k in range(1, M + 1):
    rk.append(rk[-1] * rho)
A = arb_mat(M + 1, N + 1)
pw = arb_poly([1])
for j in range(N + 1):
    for k in range(M + 1):
        A[k, j] = (pw[k] + (m[k - j] if k >= j else 0)) * rk[k] / rk[j]
    pw = mullow(pw, f, M + 1)
JN = arb_mat(N + 1, N + 1)
for i in range(N + 1):
    for j in range(N + 1):
        JN[i, j] = A[i, j]
B = arb_mat(d, N + 1)
for i in range(d):
    for j in range(N + 1):
        B[i, j] = A[N + 1 + i, j]
print(f"Kneser N={N} d={d} prec={ctx.prec}: Matrix ({time.time() - t0:.0f}s)", flush=True)

PTS = randpunkte()
npk = len(PTS)
fa = acb_poly([acb(c) for c in f.coeffs()])
ma = acb_poly([acb(c) for c in m.coeffs()])
W = acb_mat(N + 1, npk)
V = acb_mat(d, npk)
KZ = []
for p, (name, par, eta) in enumerate(PTS):
    e = to_acb(eta)
    inv_f = inv_series(acb_poly([e]) - fa, M + 1, acb_poly)
    kz = acb_poly([1 / e ** (k + 1) for k in range(M + 1)])
    w = inv_f + mullow(ma, kz, M + 1)
    for k in range(M + 1):
        val = w[k] * rk[k]
        if k <= N:
            W[k, p] = val
        else:
            V[k - N - 1, p] = val
    KZ.append([kz[k] * rk[k] for k in range(N + 1)])
JNc = acb_mat(JN)
X = JNc.solve(W)
Zt = acb_mat(JN.transpose()).solve(acb_mat(B.transpose()))
Z = Zt.transpose()
Y = JNc.solve(Zt)
G = Z * Zt
for i in range(d):
    G[i, i] = G[i, i] + 1
O = acb_mat(B) * X - V
Cc = Y * G.solve(O)
print(f"   Lösungen ({time.time() - t0:.0f}s)", flush=True)
out = []
for p, (name, par, eta) in enumerate(PTS):
    R = Rh = Rc = 0.0
    for k in range(N + 1):
        y = X[k, p] - KZ[p][k]
        cc = Cc[k, p]
        R += abs(complex(y.mid()))
        Rc += abs(complex(cc.mid()))
        Rh += abs(complex((y - cc).mid()))
    out.append({"rand": name, "par": par, "eta": [float(mp.re(eta)), float(mp.im(eta))], "R": R, "Rhat": Rh, "Rc": Rc})
(HIER / f"x12_N{N}_d{d}.json").write_text(json.dumps({"N": N, "d": d, "punkte": out}), encoding="utf-8")
L0 = out[0]
sel = lambda r, v: next(o for o in out if o["rand"] == r and abs(o["par"] - v) / v < 0.08)
print(f"   L: R = {L0['R']:.4f}, R̂ = {L0['Rhat']:.4f}, Rc = {L0['Rc']:.4f}", flush=True)
for r, vs in (("in", [0.01, 0.1, 1, 10, 103.38]), ("out", [0.01, 0.1, 1, 10, 76.618])):
    for v in vs:
        try:
            o = sel(r, v)
            print(f"   {r} {o['par']:9.4f}°: R = {o['R']:.4e}  R̂ = {o['Rhat']:.4e}  Rc = {o['Rc']:.4e}", flush=True)
        except StopIteration:
            pass
rad = max(float(X[k, p].rad()) for k in range(0, N + 1, 41) for p in range(0, npk, 7))
print(f"   max Radius (Stichprobe) {rad:.1e}  [{time.time() - t0:.0f}s]", flush=True)
