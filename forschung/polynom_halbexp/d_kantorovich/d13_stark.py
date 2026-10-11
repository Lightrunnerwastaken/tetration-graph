"""D13 — Starke Konvergenz J_N(p_N)^{-1} P_N v → G v für feste Eingaben v (trotz Eckschicht):
v1 = e^x, v2 = Re-Teil von (1 − x/L)^{σ} (σ = 5.447 + 1.058i, Knesers Singularitätstyp), v3 = (1 − x/L)^{2.5} + konj.
Fehler |J_N^{-1}P_N v − J_M^{-1}P_M v|_w (auf Grad ≤ N, plus Schwanz der Referenz) mit M = 250, w = |L|^k (1+k)^s.
"""
import sys

import mpmath as mp
import numpy as np
from flint import arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [40, 60, 80, 100, 130, 160, 200, 230]
M = 250
p = K.loesungen(M, bits=1200)
L = 1.3745570107437075
Lc = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
mp.mp.dps = 60


def binom_reihe(sig, n):
    """Koeffizienten von (1 − x/L)^sig + konj, k = 0..n."""
    c, out = mp.mpc(1), []
    for k in range(n + 1):
        out.append(2 * mp.re(c))
        c = c * (k - sig) / (k + 1) / Lc
    return out


eingaben = {
    "e^x": [mp.mpf(1) / mp.factorial(k) for k in range(M + 1)],
    "(1-x/L)^(5.447+1.058i)": binom_reihe(mp.mpc("5.4468", "1.0578"), M),
    "(1-x/L)^2.5": binom_reihe(mp.mpf("2.5"), M),
}


def loese(N, v):
    ctx.prec = 500 + 3 * N
    a = [arb(x.mid()) for x in p[N]]
    _, J = K.FJ(a, N)
    rhs = [arb(mp.nstr(v[k], 50)) for k in range(N + 1)]
    x = K.loese(J, rhs)
    return np.array([float(t.mid()) for t in x])


for name, v in eingaben.items():
    ref = loese(M, v)
    k = np.arange(M + 1)
    zeile = []
    for N in NS:
        x = loese(N, v)
        d = ref.copy()
        d[:N + 1] -= x
        for s in (0, 2, 4):
            w = L ** k * (1.0 + k) ** s
            zeile.append((N, s, float(np.abs(d * w).sum()), float(np.sqrt(((d * w) ** 2).sum()))))
    print(f"Eingabe {name}: |G v|_(ℓ1, s=0) ≈ {np.abs(ref * L ** k).sum():.4g}")
    for s in (0, 2, 4):
        print(f"   s={s}: ℓ1-Fehler " + ", ".join(f"N={N}: {e1:.2e}" for (N, ss, e1, e2) in zeile if ss == s))
        print(f"        ℓ2-Fehler " + ", ".join(f"N={N}: {e2:.2e}" for (N, ss, e1, e2) in zeile if ss == s), flush=True)
