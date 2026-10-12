"""k57c — Christoffel-Darboux-Probe am Steg (BERICHT §39.4). Kubisch.
Vergleicht K_n(t,t)·m(t) = ε_tᵀJ̃_n⁻¹ẽ_t·m(t) (eine Lösung je Punkt) mit dem logarithmischen Ableitungs-Zähler t·δ_n′(t)/δ_n(t) des Kopfes.
Vermutung: K_n(t,t)·m(t) ≈ t·δ_n′(t)/δ_n(t) + O(1) auf der Γin-Seite ausserhalb D (dann ist ΔK_n(t,t)·m(t) ≈ κ(t) die Ableitung in n).
Aufruf: python k57c_christoffel.py 300
"""
import os
import sys

import mpmath as mp

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_mat, acb_mat, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM, kubik_koeff  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)


def main():
    n = int(sys.argv[1])
    bits = int(1.3 * n) + 400
    ctx.prec = bits
    mp.mp.dps = int(bits / 3.32) + 10
    nd = mp.mp.dps
    Jr = J_tilde("kubisch", n, bits=bits); ctx.prec = bits
    e = arb_mat(n + 1, 1); e[n, 0] = arb(1)
    C = Jr.solve(e, algorithm="approx")
    c = [mp.mpf(C[k, 0].mid().str(nd, radius=False)) for k in range(n + 1)]
    J = acb_mat(Jr)
    L = mp.mpc(L_RE, L_IM); rho = abs(L)
    ck = [mp.mpf(x.mid().str(nd, radius=False)) for x in kubik_koeff()]

    def f(x): return ck[0] + x * (ck[1] + x * (ck[2] + x * ck[3]))

    def fp(x): return ck[1] + x * (2 * ck[2] + x * 3 * ck[3])

    print("n = %d: θ, s, K_n(t,t)·m(t)/n, t·δ_n′/δ_n /n, Differenz (absolut, nicht durch n geteilt)" % n)
    for th in (100.0, 120.0, 150.0, 180.0):
        for s in (0.0, 0.03):
            t = rho * (1 + s) * mp.expjpi(th / 180)
            et = acb_mat(n + 1, 1)
            for k in range(n + 1):
                v = (rho / t) ** k
                et[k, 0] = acb(arb(mp.nstr(mp.re(v), nd)), arb(mp.nstr(mp.im(v), nd)))
            X = J.solve(et, algorithm="approx")
            w = t / rho; K = mp.mpc(0); p = mp.mpc(1)
            for k in range(n + 1):
                K += mp.mpc(mp.mpf(X[k, 0].real.mid().str(nd, radius=False)), mp.mpf(X[k, 0].imag.mid().str(nd, radius=False))) * p
                p *= w
            d = mp.mpc(0); dd = mp.mpc(0)
            for k in range(n, -1, -1):
                dd = dd * w + d
                d = d * w + c[k]
            logab = w * dd / d          # t·δ′/δ  (δ in w = t/ρ, t d/dt = w d/dw)
            Km = K * fp(f(t))
            print("  θ = %5.0f°, s = %.2f:  %s   %s   Δ = %s" % (th, s, mp.nstr(Km / n, 6), mp.nstr(logab / n, 6), mp.nstr(Km - logab, 5)))


if __name__ == "__main__":
    main()
