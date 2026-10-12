"""Strang V — (S_γ) am echten p_N für N = 1500, 2000 (Auftrag Leitung 23:45): ‖(Ĵ^{p_N})⁻¹‖ (Tangente, p_N'∘p_N exakt),
‖(Ĵ_sek(p_N))⁻¹‖ (Sekante, M_P = Σ_{n≤5} g_n η^{n−1}/n!, g_n = f^{(n)}∘f, η = p_N − T_N f_K) und ‖(Ĵ^{f_K})⁻¹‖, Δ_inv = ‖(Ĵ^{p_N})⁻¹ − (Ĵ^{f_K})⁻¹‖.
p_N aus v33_pN_<N>.json (Newton-Sprung, ≥ 2650 Bit). Linse wie v12 (Ks Profilform, Höhe S). Präzision 1.2N + 400 Bit (wie v12/v20).
Aufruf: python v34_norm_pN_gross.py N S [teile: tan,fk,sek]
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from flint import arb, arb_poly, ctx

ARGS = sys.argv[1:]
N = int(ARGS[0])
S = ARGS[1]
TEILE = (ARGS[2] if len(ARGS) > 2 else "tan,fk,sek").split(",")
sys.argv = ["v20", "100", S]
import v20_pN_S as V20  # noqa: E402
V12 = V20.V12

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent


def tr(p, n):
    return p.truncate(n + 1)


def koeff(p, n):
    return [p[i] if i <= p.degree() else arb(0) for i in range(n + 1)]


def reihen_div(z, nn, n):
    zc, nc = koeff(z, n), koeff(nn, n)
    q = [arb(0)] * (n + 1)
    for i in range(n + 1):
        s_ = zc[i]
        for j in range(i):
            s_ -= q[j] * nc[i - j]
        q[i] = s_ / nc[0]
    return arb_poly(q)


def main():
    F, it, d = V12.theodorsen()
    bits = int(1.2 * N) + 400
    ctx.prec = bits
    rho = abs(V12.L_arb())
    a, mf = V20.koeff_fK(N)
    if "tan" in TEILE or "sek" in TEILE:
        p = [arb(arb(s).mid()) for s in json.loads((HIER / f"v33_pN_{N}.json").read_text(encoding="utf-8"))["p"]]
        assert len(p) == N + 1
        P = arb_poly(p)
        eps = sum(abs(float((p[k] - a[k]).mid())) * float(rho.mid()) ** k for k in range(N + 1))
    else:
        eps = float("nan")
    print(f"Kneser N = {N}, Linse S = {S}: ε_N = ‖p_N − T_N f_K‖_ℓ¹ρ = {eps:.3e}  [{bits} Bit]", flush=True)
    inv = {}

    def norm(name, acoef, mcoef):
        t0 = time.time()
        V12.modell_koeff = lambda n_: (acoef, mcoef)
        J = V12.J_hat(N, F, bits)
        ctx.prec = bits
        X = V20.inverse(J)
        inv[name] = X
        v = float(np.linalg.norm(V20.zu_np(X), 2))
        print(f"   {name}: ‖Ĵ⁻¹‖ = {v:.6f}  [{time.time() - t0:.0f} s]", flush=True)
        return v

    if "tan" in TEILE:
        r = arb_poly([0])
        for c in reversed(koeff(P.derivative(), N - 1)):
            r = tr(r * P, N) + c
        norm("Tangente p_N", p, koeff(r, N))
    if "fk" in TEILE:
        norm("f_K", a, mf)
    if "tan" in TEILE and "fk" in TEILE:
        ctx.prec = bits
        dn = float(np.linalg.norm(V20.zu_np(inv["Tangente p_N"] - inv["f_K"]), 2))
        print(f"   Δ_inv = {dn:.3e}  (Δ_inv/ε_N = {dn / eps:.3e})", flush=True)
    if "sek" in TEILE:
        ctx.prec = bits
        aa = [arb(x) for x in V20.AK[:N + 10]]
        ne = N + 7
        f = arb_poly(aa[:ne + 2])
        fp = f.derivative()
        ex = arb_poly([arb(1) / arb.fac_ui(k) for k in range(ne + 1)])
        g = [None, reihen_div(ex, fp, ne)]
        for k in range(1, 6):
            g.append(reihen_div(g[k].derivative(), fp, ne - k))
        eta = tr(P - arb_poly(aa[:N + 1]), N)
        MP = tr(g[1], N)
        ep = arb_poly([1])
        fak = arb(1)
        for k in range(1, 6):
            ep = tr(ep * eta, N)
            fak = fak * k
            MP = MP + tr(tr(g[k + 1], N) * ep, N) * (1 / (fak * (k + 1)))
        norm("Sekante p_N", p, koeff(MP, N))


if __name__ == "__main__":
    main()
