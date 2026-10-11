"""W44 — Rayleigh-Quotient eines Testvektors aus dem Ω_y-Modell (Auftrag Runde 6, Punkt 1 „obere Schranke“).
Testvektor: P(w) = Π_j (w − 1/z_j) mit N Modellpunkten z_j auf ∂Ω_y (Quantile des harmonischen Masses ν: Anteil 0.4728 auf B, Rest auf der
Linsen-Nulllinie), als skalierter Koeffizientenvektor v (v(t) = Σ v_k t^{−k} = P(1/t)), normiert ‖v‖₂ = 1.
Gemessen (arb): R = ‖J̃ᵀv‖/‖v‖ gegen 1/σ₁ (w01); σ₁ ≥ 1/R ist eine STRENGE untere Schranke (bei Einschluss in arb); Überlapp |⟨v, v₁⟩|.
Aufruf: python w44_testvektor.py N1,N2,...
"""
import sys
import time

import numpy as np
from flint import acb, acb_poly, arb, ctx

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
rho = W.RHO
M = np.load(W.HIER + "/w42_modell_kurven.npz")
B, Z, gB, gZ = M["B"], M["Z"], M["gB"], M["gZ"]
B = np.concatenate([np.conj(B[::-1]), B[1:]])
gB = np.concatenate([gB[::-1], gB[1:]])
Z = np.concatenate([Z, np.conj(Z[::-1])[1:]])
gZ = np.concatenate([gZ, gZ[::-1][1:]])


def quantile_punkte(curve, g, n):
    ds = np.abs(np.diff(curve))
    m = np.concatenate([[0], np.cumsum(0.5 * (g[1:] + g[:-1]) * ds)])
    m /= m[-1]
    q = (np.arange(n) + 0.5) / n
    s = np.interp(q, m, np.arange(len(curve)))
    i = np.floor(s).astype(int)
    fr = s - i
    i = np.minimum(i, len(curve) - 2)
    return curve[i] * (1 - fr) + curve[i + 1] * fr


coef, _ = W.kubik_koeff()
for N in [int(v) for v in sys.argv[1].split(",")]:
    t0 = time.time()
    n = N + 1
    ctx.prec = int(1.5 * N) + 400
    ctx.cap = n + 10
    nB = int(round(0.4728 * N))
    zs = np.concatenate([quantile_punkte(B, gB, nB), quantile_punkte(Z, gZ, N - nB)])
    p = acb_poly([1])
    for z in zs:
        w = 1 / complex(z)
        p = p * acb_poly([acb(-w.real, -w.imag), 1])
    c = [p[k] for k in range(n)]                       # v_k = Koeffizient von w^k
    nrm = sum((abs(ck) ** 2 for ck in c), arb(0)).sqrt()
    v = [ck / nrm for ck in c]
    a = W.kubik_a(N, coef)
    P = W.potenzen(a, N)
    J = W.J_mat(a, N, P)
    Dk = [arb(rho) ** k for k in range(n)]
    Dv = [Dk[k] * v[k] for k in range(n)]
    # (J̃ᵀ v)_j = ρ^{−j} Σ_k J[k][j] ρ^k v_k
    res = arb(0)
    for j in range(n):
        s = acb(0)
        for k in range(n):
            s += J[k, j] * Dv[k]
        res += abs(s / Dk[j]) ** 2
    R = res.sqrt()
    d = np.load(W.HIER + f"/w01_N{N}.npz")
    S, Vt = d["S"], d["Vt"]
    vf = np.array([complex(x.real.mid(), x.imag.mid()) for x in v])
    ov = [abs(np.vdot(Vt[i], vf)) for i in range(min(3, Vt.shape[0]))]
    print(f"N={N}: R = {float(R.mid()):.4e} (Einschluss rad {float(R.rad()):.1e}) ⇒ σ₁ ≥ {1/float(R.mid()):.4g}; wahr σ₁ = {S[0]:.4g}, "
          f"1/(Rσ₁) = {1/(float(R.mid())*S[0]):.3g}; Überlapp |⟨v,v_i⟩| i=1,2,3: " + ", ".join(f"{o:.3f}" for o in ov)
          + f"  ({time.time()-t0:.0f}s)", flush=True)
