"""D17 — Starke Konvergenz J_N^{-1}P_N v → G v bis N = 400, obwohl |J_N^{-1}| (Eckspitze) exponentiell wächst.
J bei q⁽²⁾ (zwei Newton-Schritte von Knesers q_N aus; Abstand zu p_N ~ 1e-20). Eingaben v = e^x und
v = (1 − x/L)^{5.447+1.058i} + konj. Fehler gegen N_ref = 400 in ℓ1_{|L|,s}; dazu sup der Lösung auf |x| = |L| nahe arg L.
"""
import json
import sys
import time

import mpmath as mp
import numpy as np
from flint import arb, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
kn = json.load(open(K.GEMEINSAM + "/kneser_halbexp_koeff_1000.json"))["a"]
L = 1.3745570107437075
mp.mp.dps = 60
Lc = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")


def binom_reihe(sig, n):
    c, out = mp.mpc(1), []
    for k in range(n + 1):
        out.append(2 * mp.re(c))
        c = c * (k - sig) / (k + 1) / Lc
    return out


NMAX = max(NS)
eingaben = {"e^x": [mp.mpf(1) / mp.factorial(k) for k in range(NMAX + 1)],
            "Kneser-Typ": binom_reihe(mp.mpc("5.4468", "1.0578"), NMAX)}
loes = {name: {} for name in eingaben}
t0 = time.time()
for N in NS:
    ctx.prec = 500 + 3 * N
    q = [arb(str(x)) for x in kn[:N + 1]]
    for _ in range(2):
        F, J = K.FJ(q, N)
        d = K.loese(J, F)
        q = [arb((q[i] - d[i]).mid()) for i in range(N + 1)]
    _, J = K.FJ(q, N)
    for name, v in eingaben.items():
        x = K.loese(J, [arb(mp.nstr(v[k], 50)) for k in range(N + 1)])
        loes[name][N] = np.array([float(t.mid()) for t in x])
    print(f"N={N} gelöst ({time.time()-t0:.0f}s)", flush=True)
phiL = float(mp.arg(Lc))
for name in eingaben:
    ref = loes[name][NMAX]
    k = np.arange(NMAX + 1)
    print(f"Eingabe {name}:")
    for N in NS[:-1]:
        d = ref.copy()
        d[:N + 1] -= loes[name][N]
        e = [np.abs(d * L ** k * (1.0 + k) ** s).sum() for s in (0, 2)]
        xs = loes[name][N] * L ** np.arange(N + 1)
        M = 8192
        Fv = np.abs(np.fft.ifft(np.concatenate([xs, np.zeros(M - N - 1)])) * M)
        phi = 2 * np.pi * np.arange(M) / M
        nah = np.abs(phi - phiL) < 0.05
        print(f"   N={N}: Fehler ℓ1 s=0 {e[0]:.2e}, s=2 {e[1]:.2e};  sup|Lösung| auf |x|=|L|: gesamt {Fv.max():.4f}, "
              f"nahe arg L (±3°) {Fv[nah].max():.4f}")
