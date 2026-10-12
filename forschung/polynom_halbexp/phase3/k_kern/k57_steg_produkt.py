"""k57 — Zeitfenster „Produkt primal × adjungiert am Steg“ (BERICHT §39). Kubisch.
Köpfe exakt: J̃_n c = e_n, J̃_nᵀ ℓ = e_n (reell). δ_n(z) = Σ c_k (z/ρ)^k, Λ_n(t) = Σ ℓ_k (ρ/t)^k, P_n = δ_n(L)Λ_n(L)/ℓ_nn.
Je Punkt t = ρ(1+s)e^{iθ} am Steg (Γin-Seite):
  Π_n(t) = δ_n(t)Λ_n(t)/(δ_n(L)Λ_n(L))            (Produkt relativ zur Ecke = ΔK_n(t,t)/P_n)
  A_P(t) = m(t)δ_n(t)(ρ/t)^n,  A_A(t) = m(t)Λ_n(t)(t/ρ)^n   (Abweichung von der Quelle; exakt Π_n = A_P·A_A/(m²·ℓ_nn·P_n))
dazu die skalaren Grössen ln|δ_n(L)|, ln|Λ_n(L)|, ln|ℓ_nn|, ln|P_n| und das radiusunabhängige Kreismittel M_n = Σ c_kℓ_k.
Aufruf: python k57_steg_produkt.py 150 300 450 600   (schreibt k57.json, druckt Raten zwischen aufeinanderfolgenden n)
"""
import json
import math
import os
import sys
import time

import mpmath as mp

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, arb_mat, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM, kubik_koeff  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
HIER = os.path.dirname(os.path.abspath(__file__))
THETAS = [0.0, 40.0, 85.0, 100.0, 120.0, 150.0, 180.0]
SS = [-0.03, -0.01, 0.0, 0.01, 0.03, 0.05]


def lauf(n):
    t0 = time.time()
    bits = int(1.3 * n) + 400
    ctx.prec = bits
    mp.mp.dps = int(bits / 3.32) + 10
    nd = mp.mp.dps
    J = J_tilde("kubisch", n, bits=bits); ctx.prec = bits
    e = arb_mat(n + 1, 1); e[n, 0] = arb(1)
    C = J.solve(e, algorithm="approx")
    Lm = J.transpose().solve(e, algorithm="approx")
    c = [mp.mpf(C[k, 0].mid().str(nd, radius=False)) for k in range(n + 1)]
    l = [mp.mpf(Lm[k, 0].mid().str(nd, radius=False)) for k in range(n + 1)]
    L = mp.mpc(L_RE, L_IM); rho = abs(L)
    ck = [mp.mpf(x.mid().str(nd, radius=False)) for x in kubik_koeff()]

    def f(x): return ck[0] + x * (ck[1] + x * (ck[2] + x * ck[3]))

    def fp(x): return ck[1] + x * (2 * ck[2] + x * 3 * ck[3])

    def delta(z):
        w = z / rho; s_ = mp.mpc(0)
        for k in range(n, -1, -1):
            s_ = s_ * w + c[k]
        return s_

    def Lam(t):
        w = rho / t; s_ = mp.mpc(0)
        for k in range(n, -1, -1):
            s_ = s_ * w + l[k]
        return s_

    dL, LL = delta(L), Lam(L)
    lnn = l[n]
    Pn = dL * LL / lnn
    Mn = mp.fsum(c[k] * l[k] for k in range(n + 1))
    res = {"n": n, "bits": bits, "ln_dL": float(mp.log(abs(dL))), "ln_LL": float(mp.log(abs(LL))), "ln_lnn": float(mp.log(abs(lnn))),
           "lnn_c": float(c[n] / lnn), "ln_P": float(mp.log(abs(Pn))), "ln_M_rel": float(mp.log(abs(Mn / (dL * LL)))), "punkte": []}
    for th in THETAS:
        for s in SS:
            t = rho * (1 + s) * mp.expjpi(th / 180)
            m = fp(f(t))
            dt, Lt = delta(t), Lam(t)
            Pi = dt * Lt / (dL * LL)
            AP = m * dt * (rho / t) ** n
            AA = m * Lt * (t / rho) ** n
            chk = AP * AA / (m ** 2 * lnn * Pn) / Pi
            res["punkte"].append({"th": th, "s": s, "ln_Pi": float(mp.log(abs(Pi))), "ln_AP": float(mp.log(abs(AP))), "ln_AA": float(mp.log(abs(AA))),
                                  "ln_m": float(mp.log(abs(m))), "chk": float(abs(chk - 1))})
    res["zeit"] = time.time() - t0
    return res


def main():
    ns = [int(a) for a in sys.argv[1:]]
    pfad = os.path.join(HIER, "k57.json")
    alle = json.load(open(pfad, encoding="utf-8")) if os.path.exists(pfad) else {}
    for n in ns:
        r = lauf(n)
        alle[str(n)] = r
        json.dump(alle, open(pfad, "w", encoding="utf-8"), indent=1)
        print("n = %d (%.0f s): ln|δ(L)| = %.3f, ln|Λ(L)| = %.3f, ln|ℓ_nn| = %.3f (c_nn/ℓ_nn = %.6f), ln|P_n| = %.4f, ln|M_n/(δ(L)Λ(L))| = %.3f, max Prüfung %.1e"
              % (n, r["zeit"], r["ln_dL"], r["ln_LL"], r["ln_lnn"], r["lnn_c"], r["ln_P"], r["ln_M_rel"], max(p["chk"] for p in r["punkte"])))
    ks = sorted(int(k) for k in alle)
    if len(ks) < 2:
        return
    print("\nRaten je Grad zwischen aufeinanderfolgenden n (" + ", ".join("%d→%d" % (a, b) for a, b in zip(ks, ks[1:])) + ")")
    for nm in ("ln_dL", "ln_LL", "ln_lnn", "ln_P", "ln_M_rel"):
        print("  %-9s " % nm + "  ".join("%+.5f" % ((alle[str(b)][nm] - alle[str(a)][nm]) / (b - a)) for a, b in zip(ks, ks[1:])))
    print("  %-14s %s" % ("θ, s", "   h_Π (Fenster)            |  A_P               |  A_A"))
    for i, p in enumerate(alle[str(ks[0])]["punkte"]):
        zeile = []
        for nm in ("ln_Pi", "ln_AP", "ln_AA"):
            zeile.append(" ".join("%+.4f" % ((alle[str(b)]["punkte"][i][nm] - alle[str(a)]["punkte"][i][nm]) / (b - a)) for a, b in zip(ks, ks[1:])))
        print("  %5.0f° %+.2f  %s | %s | %s" % (p["th"], p["s"], zeile[0], zeile[1], zeile[2]))


if __name__ == "__main__":
    main()
