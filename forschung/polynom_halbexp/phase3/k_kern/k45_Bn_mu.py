"""k45 — Test von Glied (c) der Warnung (BERICHT §22.6): kompensiert die globale Normierung die primale Nicht-Lokalität?
B_n := √n·‖c_n‖_{ℓ²_ρ}·|Λ_n(L)| (Gs Eckamplitude; c_n = J̃_n⁻¹e_n, Λ_n mit ℓ_nn = 1, also ⟨ℓ_n, J̃c_n⟩ = 1) und der primale Lokalitätsquotient
Q_n := ‖c_n‖/(|δ_n(L)|·|d_n|^{1/2}) (d_n = ρ n^{−a}) für kubisch μ = √L und |μ| = 2 (arg √L). Nach Lemma K13: sup_{D̄}|δ_n|·|Λ_n(L)| ≤ B_n, also trägt nur B_n
die primale Seite — wächst B_n, wenn Q_n wächst?
Aufruf: python k45_Bn_mu.py <sqrtL|2.0> 100,150,200,...
"""
import os
import sys
import time

import mpmath as mp

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_mat, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def main():
    which = sys.argv[1]
    Ns = [int(s) for s in sys.argv[2].split(",")]
    mp.mp.dps = 30
    Lm = mp.mpc(L_RE, L_IM)
    rho = float(abs(Lm))
    if which == "sqrtL":
        mu, fac, mu_m = None, 1.3, mp.sqrt(Lm)
    else:
        mu_m = float(which) * mp.expjpi(mp.arg(mp.sqrt(Lm)) / mp.pi)
        mu, fac = acb(arb(str(mp.re(mu_m))), arb(str(mp.im(mu_m)))), 3.0
    a = 1 + 1j * mp.log(mu_m) / (2 * mp.pi)
    for n in Ns:
        t0 = time.time()
        bits = int(fac * n) + 500
        J = J_tilde("kubisch", n, bits=bits, mu=mu)
        ctx.prec = bits
        e = arb_mat(n + 1, 1); e[n, 0] = arb(1)
        c = J.solve(e, algorithm="approx")
        l = J.transpose().solve(e, algorithm="approx")
        lnn = l[n, 0]
        nc = arb(0)
        for k in range(n + 1):
            nc += c[k, 0] ** 2
        nc = nc.sqrt()
        rhoa = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
        La = acb(arb(L_RE), arb(L_IM))
        q = acb(rhoa) / La
        s = acb(0)
        for k in range(n, -1, -1):
            s = s * q + acb(l[k, 0] / lnn)
        x = La / acb(rhoa)
        d = acb(0)
        for k in range(n, -1, -1):
            d = d * x + acb(c[k, 0])
        ncf = float(nc.mid()); lam = abs(complex(float(s.real.mid()), float(s.imag.mid())))
        dL = abs(complex(float(d.real.mid()), float(d.imag.mid())))
        dn = rho * float(abs(mp.power(n, -a)))
        ok = nc.rad() < 1e-3 * nc.mid() if hasattr(nc, "rad") else True
        print("n = %d (%d Bit, %.0f s): B_n = %.5f   Q_n = ‖c‖/(|δ(L)|d^{1/2}) = %.4f   (‖c‖ = %.3e, |Λ(L)| = %.3e, Radius ok: %s)"
              % (n, bits, time.time() - t0, (n ** 0.5) * ncf * lam, ncf / (dL * dn ** 0.5), ncf, lam, ok), flush=True)


if __name__ == "__main__":
    main()
