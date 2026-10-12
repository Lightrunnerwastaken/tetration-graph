"""k57b — Folgemessung zu k57 (BERICHT §39). Kubisch, Köpfe exakt wie k57.
(1) Ist das Diagonal-Produkt am Steg explizit? V_n(t) := ΔK_n(t,t)·m(t) = δ_n(t)Λ_n(t)m(t)/ℓ_nn (Vermutung: → 1 auf der Γin-Seite, „lokales Symbol 1/m“).
(2) Wo wächst das Produkt? Π_n(t) = δ_n(t)Λ_n(t)/(δ_n(L)Λ_n(L)) auf den Spiralen Σ₁ (u_A = 80°) und Σ₂ (u_A = −90°), r = 0.01 … 0.2.
(3) Kreismittel: M_n/ℓ_nn = Σ_k c_kℓ_k/ℓ_nn = tr J̃_n⁻¹ − tr J̃_{n−1}⁻¹ (exakt).
Aufruf: python k57b_produkt_symbol.py 150 300 450 600
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
KREIS = [(th, s) for th in (0.0, 40.0, 90.0, 100.0, 120.0, 150.0, 180.0) for s in (-0.01, 0.0, 0.03)]
SPIRALE = [(uA, r) for uA in (80.0, -90.0) for r in (0.01, 0.03, 0.05, 0.1, 0.2)]


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

    lnn = l[n]
    dL, LL = delta(L), Lam(L)
    res = {"n": n, "trink": float(mp.fsum(c[k] * l[k] for k in range(n + 1)) / lnn), "kreis": [], "spirale": []}
    for th, s in KREIS:
        t = rho * (1 + s) * mp.expjpi(th / 180)
        v = delta(t) * Lam(t) * fp(f(t)) / lnn
        res["kreis"].append({"th": th, "s": s, "re": float(mp.re(v)), "im": float(mp.im(v))})
    for uA, r in SPIRALE:
        th = math.radians(uA + 71.46) + math.radians(1.623) * math.log(r / float(rho))
        t = L + r * mp.mpc(math.cos(th), math.sin(th))
        Pi = delta(t) * Lam(t) / (dL * LL)
        res["spirale"].append({"uA": uA, "r": r, "ln_Pi": float(mp.log(abs(Pi))), "aussen": bool(abs(t) > rho)})
    res["zeit"] = time.time() - t0
    return res


def main():
    ns = [int(a) for a in sys.argv[1:]]
    pfad = os.path.join(HIER, "k57b.json")
    alle = json.load(open(pfad, encoding="utf-8")) if os.path.exists(pfad) else {}
    for n in ns:
        alle[str(n)] = lauf(n)
        json.dump(alle, open(pfad, "w", encoding="utf-8"), indent=1)
        print("n = %d fertig (%.0f s), Δtr J̃⁻¹ = M_n/ℓ_nn = %.6f" % (n, alle[str(n)]["zeit"], alle[str(n)]["trink"]))
    ks = sorted(int(k) for k in alle)
    print("\n(1) V_n(t) = ΔK_n(t,t)·m(t) (komplex) für n = " + " / ".join(str(k) for k in ks))
    for i, p in enumerate(alle[str(ks[0])]["kreis"]):
        print("  θ = %5.0f°, s = %+.2f: " % (p["th"], p["s"]) + "  ".join("%+.4f%+.4fi" % (alle[str(k)]["kreis"][i]["re"], alle[str(k)]["kreis"][i]["im"]) for k in ks))
    print("\n(2) Rate h_Π je Grad auf den Spiralen (Fenster " + ", ".join("%d→%d" % (a, b) for a, b in zip(ks, ks[1:])) + ")")
    for i, p in enumerate(alle[str(ks[0])]["spirale"]):
        print("  u_A = %+5.0f°, r = %.2f (%s): " % (p["uA"], p["r"], "aussen" if p["aussen"] else "innen ") +
              "  ".join("%+.4f" % ((alle[str(b)]["spirale"][i]["ln_Pi"] - alle[str(a)]["spirale"][i]["ln_Pi"]) / (b - a)) for a, b in zip(ks, ks[1:])) +
              "   (ln|Π| bei n = %d: %.2f)" % (ks[-1], alle[str(ks[-1])]["spirale"][i]["ln_Pi"]))


if __name__ == "__main__":
    main()
