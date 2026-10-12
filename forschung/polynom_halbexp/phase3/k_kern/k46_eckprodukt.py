"""k46 — Eckprodukt P_n := δ_n(L)Λ_n(L) (normiert ⟨ℓ_n, J̃c_n⟩ = 1) und der Diagonalkern K_n(L,L) = ε_Lᵀ J̃_n⁻¹ ẽ_L = L·(J_n⁻¹T_nk_L)(L).
Exakt: K_n − K_{n−1} = P_n (Rang-1-Ränderung). Fragen: (1) Hat P_n eine schnelle Phase oder ist es glatt in n? (2) Skalierung: P_n·n·d_n/(aL) → C₀?
(Hypothese „Euler-Identität“: auf der Kegelskala K_n ≈ (L/d_n)𝒢 ⇒ P_n ≈ (aL/(n d_n))·E𝒢(0,0).) (3) Vergleich mit Os Eckformel K_n(L,L)(μ + 1/μ) ≈ n + 1 + L·T(L).
Aufruf: python k46_eckprodukt.py 200,201,202,203,300,301,302,400,401,402
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
    Ns = [int(s) for s in sys.argv[1].split(",")]
    mp.mp.dps = 30
    Lm = mp.mpc(L_RE, L_IM)
    a = 1 + 1j * Lm / (4 * mp.pi)
    mu = mp.sqrt(Lm)
    prev = None
    for n in Ns:
        t0 = time.time()
        bits = int(1.3 * n) + 400
        J = J_tilde("kubisch", n, bits=bits)
        ctx.prec = bits
        Jc = acb_mat(J)
        rhoa = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
        La = acb(arb(L_RE), arb(L_IM))
        eL = acb_mat(n + 1, 1)   # ẽ_L = (ρ/L)^k
        epsL = acb_mat(n + 1, 1)  # ε_L = (L/ρ)^k
        q = acb(rhoa) / La
        p = acb(1)
        for k in range(n + 1):
            eL[k, 0] = p
            epsL[k, 0] = 1 / p
            p = p * q
        e = acb_mat(n + 1, 1); e[n, 0] = acb(1)
        X = Jc.solve(eL, algorithm="approx")          # J̃⁻¹ẽ_L
        c = Jc.solve(e, algorithm="approx")           # c_n
        l = Jc.transpose().solve(e, algorithm="approx")  # ℓ_n (unnormiert, ℓ_nn = c_nn)
        K = acb(0); dL = acb(0); lam = acb(0)
        for k in range(n + 1):
            K += epsL[k, 0] * X[k, 0]
            dL += epsL[k, 0] * c[k, 0]
            lam += l[k, 0] * eL[k, 0]
        lnn = l[n, 0]
        P = dL * lam / lnn
        Kc = complex(float(K.real.mid()), float(K.imag.mid()))
        Pc = complex(float(P.real.mid()), float(P.imag.mid()))
        dn = abs(complex(float(mp.re(1.3745570107437074 * mp.power(n, -a))), float(mp.im(1.3745570107437074 * mp.power(n, -a)))))
        scal = Pc * n * dn / (complex(a) * complex(Lm))
        o12 = Kc * complex(mu + 1 / mu) - (n + 1)
        print("n = %d (%.0f s): P_n = %.5f%+.5fi |P_n| = %.5f arg = %+.2f° | P_n·n·|d_n|/(aL) = %.5f%+.5fi | K_n(L,L) = %.4f%+.4fi | K(μ+1/μ) − (n+1) = %.4f%+.4fi"
              % (n, time.time() - t0, Pc.real, Pc.imag, abs(Pc), float(mp.degrees(mp.arg(mp.mpc(Pc)))), scal.real, scal.imag, Kc.real, Kc.imag, o12.real, o12.imag), flush=True)
        if prev is not None and prev[0] == n - 1:
            inc = Kc - prev[1]
            print("      Kontrolle K_n − K_{n−1} = %.5f%+.5fi (gegen P_n)" % (inc.real, inc.imag), flush=True)
        prev = (n, Kc)


if __name__ == "__main__":
    main()
