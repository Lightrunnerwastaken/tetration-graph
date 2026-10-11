"""W28 — Dipol-Modell (w19) mit GEMESSENEN A-Bogen-Daten λ_A (statt 0): Φ_B = λ_A auf dem A-Bogen, 0 auf h(A-Bogen).
Frage: wird c = Φ_B(−1) dann ≈ 0.012 und das Profil nahe L richtig?"""
import sys, json
import numpy as np
import w_kern as W
sys.stdout.reconfigure(encoding="utf-8")
src = open("w19_dipol.py", encoding="utf-8").read()
src = src.split("ds = abs(t_)")[0]
exec(src)
# gemessene λ_A (Eckteil, asymptotisch) als Funktion des Abstands |η − L| (unskaliert)
NS = [500, 600, 700, 800, 900, 1000]
D = {N: json.load(open(f"w24_N{N}.json")) for N in NS}
Ns = np.array(NS, float); Am = np.vstack([np.ones_like(Ns), np.log(Ns), Ns]).T
def lnK(N, nm):
    p = D[N]["pts"][nm]; S = D[N]["S"]; r = max(1, sum(1 for s in S[:len(p["lnK"])] if s > 1.40))
    return 0.5*np.log(sum(np.exp(2*v) for v in p["lnK"][:r]))
lamL = np.linalg.lstsq(Am, np.array([lnK(N,"L") for N in NS]), rcond=None)[0][2]
da, la = [0.0], [0.0]
for nm in D[NS[0]]["pts"]:
    if nm[0] == "A":
        co = np.linalg.lstsq(Am, np.array([lnK(N,nm) for N in NS]), rcond=None)[0]
        da.append(float(nm[1:])); la.append(co[2] - lamL)
o = np.argsort(da); da = np.array(da)[o]; la = np.array(la)[o]
scale = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
scale2 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
ds = abs(t_)
src2 = cur + nr * np.minimum(5 * ds, 0.3 * abs(cur - Lt) * abs(cur - np.conj(Lt)) + 1e-12)
Mat = np.log(np.abs(cur[:, None] - src2[None, :]))
Mat = np.vstack([np.hstack([Mat, np.ones((len(cur), 1))]), np.hstack([np.ones((1, len(src2))), np.zeros((1, 1))])])
# Randdaten: A-Bogen-Teil von cur sind die ersten K+1 Punkte (A von L̃ nach L̃̄)
dat = np.log(np.abs(cur))
nA = len(A)
distA = np.minimum(abs(A*rho - W.L), abs(A*rho - np.conj(W.L)))
dat[:nA] += scale * np.interp(distA, da, la)
nI = len(cur) - nA
hI = hs[::-1][1:-1]
distI = np.minimum(abs(A[::-1][1:-1]*rho - W.L), abs(A[::-1][1:-1]*rho - np.conj(W.L)))
dat[nA:] += scale2 * np.interp(distI, da, la)
rhs = np.concatenate([dat, [0.0]])
sol, *_ = np.linalg.lstsq(Mat, rhs, rcond=None); q, C = sol[:-1], sol[-1]
def PhiB(t): return -np.log(abs(t)) + np.sum(q*np.log(np.abs(t - src2))) + C
argL = np.degrees(np.angle(Lt))
print(f"Skala {scale}: c_pred = Φ_B(−1) = {PhiB(-1+0j):.5f}")
print(" τ: " + " ".join(f"{tau}:{PhiB(np.exp(1j*np.deg2rad(argL+tau))):.4f}" for tau in (1, 2, 4, 8, 15, 30, 60, 103.38)))
