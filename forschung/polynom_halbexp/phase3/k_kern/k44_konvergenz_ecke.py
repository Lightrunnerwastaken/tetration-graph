"""k44 — Konvergiert die Abschnittslösung u_N = J_N⁻¹T_N W nahe L auf dem Kreis, wenn die Daten in L singulär sind?
Daten W(x) = (1 − x/L)^σ (1 − x/L̄)^σ, σ = 3.5 (Schlitze nach aussen; Koeffizienten ~ k^{−4.5}). Modelle: kubisch mit μ = √L (Drehung 1.45°/e) und |μ| = 2
(Drehung 6.3°/e; dort soll die primale Γout-Hälfte schon bei kleinem n in den Wachstumssektor wandern, k43). Ausgabe |u_N(z)| auf dem Kreis bei
τ = arg z − arg L ∈ {−0.05, −0.02, −0.008, −0.003, +0.003, +0.008, +0.02} (Γout negativ) und Innenpunkt 0.9L; dazu |u_N − u_{N_max}|.
Aufruf: python k44_konvergenz_ecke.py <|μ| oder 'sqrtL'> 100,200,300,400
"""
import os
import sys
import time

import mpmath as mp

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_mat, acb_mat, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def main():
    muarg = sys.argv[1]
    Ns = [int(s) for s in sys.argv[2].split(",")]
    mp.mp.dps = 40
    Lm = mp.mpc(L_RE, L_IM)
    rho = abs(Lm)
    sigma = mp.mpf("3.5")
    if muarg == "sqrtL":
        mu = None
        fac = 1.3
    else:
        mu_m = float(muarg) * mp.expjpi(mp.arg(mp.sqrt(Lm)) / mp.pi)
        mu = (mp.re(mu_m), mp.im(mu_m))
        fac = 3.0
    taus = [-0.05, -0.02, -0.008, -0.003, 0.003, 0.008, 0.02]
    pts = [rho * mp.expjpi((mp.arg(Lm) + t) / mp.pi) for t in taus] + [0.9 * Lm]
    res = {}
    for N in Ns:
        t0 = time.time()
        bits = int(fac * N) + 500
        ctx.prec = bits
        mua = None if mu is None else acb(arb(str(mu[0])), arb(str(mu[1])))
        J = J_tilde("kubisch", N, bits=bits, mu=mua)
        ctx.prec = bits
        La = acb(arb(L_RE), arb(L_IM))
        rhoa = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
        # Binomialreihe (1 − x/L)^σ in arb
        sg = arb("3.5")
        b = [acb(1)]
        for k in range(1, N + 1):
            b.append(b[-1] * (sg - (k - 1)) / k * (-1) / La)
        bc = [x.conjugate() for x in b]
        W = acb_mat(N + 1, 1)
        for k in range(N + 1):
            s = acb(0)
            for j in range(k + 1):
                s += b[j] * bc[k - j]
            W[k, 0] = s * acb(rhoa) ** k
        Jc = acb_mat(J)
        x = Jc.solve(W, algorithm="approx")
        xs = [x[j, 0] for j in range(N + 1)]
        vals = []
        for z in pts:
            q = acb(arb(str(mp.re(z))), arb(str(mp.im(z)))) / acb(rhoa)
            s = acb(0)
            for j in range(N, -1, -1):
                s = s * q + xs[j]
            vals.append(complex(float(s.real.mid()), float(s.imag.mid())))
        res[N] = vals
        print("N = %d (%.0f s): |u_N(z)| = %s" % (N, time.time() - t0, " ".join("%.4g" % abs(v) for v in vals)), flush=True)
    Nm = max(Ns)
    for N in sorted(Ns)[:-1]:
        print("  |u_%d − u_%d| = %s" % (N, Nm, " ".join("%.2e" % abs(a - b) for a, b in zip(res[N], res[Nm]))), flush=True)
    print("Punkte: τ = %s, Innenpunkt 0.9L" % taus)


if __name__ == "__main__":
    main()
