"""Gegenbeispiel „alle“: a, c, d und θ̂_T(−1), θ̂_T(−2), θ̂_T(−3) zugleich exakt null, T ≠ Kneser.

Sehnenverfahren mit der *theoretischen* Ableitung aus S1 (δC = Σ_k F̂_C(k)·c_k, c_k = conj p̂(k)), gebaut aus den
θ-Koeffizienten a_n (etappe2). Das ist robuster als die numerische 12×12-Jacobi-Matrix.
p = t sin 14πx + Σ_{k=1..6} α_k sin 2πkx + β_k (cos 2πkx − 1),  c_k = (β_k + iα_k)/2.
Aufruf: python s2b_alle.py DPS K T
"""
import json
import sys
import time

import mpmath as mp
import numpy as np

import e4_lib as E

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
k = int(sys.argv[2]) if len(sys.argv) > 2 else 6
t = mp.mpf(sys.argv[3]) if len(sys.argv) > 3 else mp.mpf("1e-4")
E.init(dps)
BED = ["a", "c", "d", "neg1", "neg2", "neg3"]
K = 6

d = json.loads((E.HIER.parent / "etappe2" / "daten" / "a_e_80_0.003_16384.json").read_text(encoding="utf-8"))
N = 60
a = np.array([complex(float(x[0]), float(x[1])) for x in d["a"]])[:N]
mul = lambda x, y: np.convolve(x, y)[:N]  # noqa: E731


def inv(x):
    y = np.zeros(N, complex); y[0] = 1 / x[0]
    for j in range(1, N):
        y[j] = -np.dot(x[1:j + 1], y[j - 1::-1]) / x[0]
    return y


g = a.copy(); g[0] = 0
n = np.arange(N)
w = 2j * np.pi * n * g; w[0] += 1
t2 = (2j * np.pi * n) ** 2 * g
F = {"a": 2 * mul(g, w), "c": -2 * mul(t2, inv(mul(w, w))), "d": mul(t2, inv(w))}
for m in (1, 2, 3):
    q = np.zeros(N, complex); q[m] = 1
    F[f"neg{m}"] = mul(q, w)
A = np.array([F[c][1:K + 1] for c in BED])
Ai = np.linalg.inv(A)
print(f"Theorie-Matrix 6×6, Singulärwerte {np.linalg.svd(A, compute_uv=False)}")


def stoer(c):
    par = [(kk + 1, 2 * mp.im(c[kk]), 2 * mp.re(c[kk])) for kk in range(K)] + [(7, t, mp.mpf(0))]
    return ("trig", par)


c = [mp.mpc(0)] * K
t0 = time.time()
for it in range(150):
    r = E.kenngroessen(stoer(c), k, moden=6, momente=4)
    f = [r[x] for x in BED]
    nf = max(abs(v) for v in f)
    print(f"  it {it}: max|C| = {mp.nstr(nf, 3)}  ({time.time()-t0:.0f}s)", flush=True)
    if nf < mp.mpf(10) ** (-dps + 4):
        break
    # Schritt in mp-Genauigkeit mit der double-Inversen (Sehnenverfahren)
    for kk in range(K):
        c[kk] -= mp.fsum(mp.mpc(complex(Ai[kk, j])) * f[j] for j in range(len(BED)))

st = stoer(c)
print("== Lösung:", [(kk, mp.nstr(al, 10), mp.nstr(be, 10)) for kk, al, be in st[1]])
for x in ("a", "b", "c", "d", "neg1", "neg2", "neg3", "neg4", "neg5", "neg6", "V2", "V3", "V4"):
    print(f"   {x:5s} |.| = {mp.nstr(abs(r[x]), 4)}")
mn = min(1 + E.p_wert(st, mp.mpf(j) / 1000)[1] for j in range(1000))
print(f"   min(1+p') = {mp.nstr(mn, 6)}")
out = {"fall": "alle", "dps": dps, "k": k, "bed": BED,
       "p": [[kk, mp.nstr(al, dps), mp.nstr(be, dps)] for kk, al, be in st[1]]}
(E.HIER / "daten" / f"gegenbeispiel_alle_{dps}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
