"""k48 — Integraldarstellung der Eckprofile:  F_β(u) = (1/Γ(β)) ∫₀^∞ t^{β−1} e^{−t − u t^{−iα}} dt
(gliedweise: ∫ t^{β−iαk−1}e^{−t}dt = Γ(β − iαk), Majorante e^{|u| t^{Im α}}·e^{−t}).  Prüfung gegen die Reihe (mpmath, hohe Präzision)
und Watson-Leitterm im Potenzsektor:  F_β(u) ≈ γ Γ(γβ) u^{−γβ}/Γ(β) · Σ_j (−1)^j Γ(γ(β+j))/(j! Γ(γβ)) u^{−γj},  γ = i/α.
Aufruf: python k48_gamma_integral.py
"""
import sys

import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")


def main():
    mp.mp.dps = 60
    L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237",
               "1.33723570143068940890116214319371061253950213846051241687")
    al = L / (4 * mp.pi)
    g = 1j / al
    betas = {"A": 0.5 - 2j * al, "P": 0.5 + 1j * al}

    def F_reihe(b, u, K=4000):
        s = mp.mpc(0); term_fac = mp.mpc(1)
        gb = mp.gamma(b)
        for k in range(K):
            c = mp.gamma(b - 1j * al * k) / (mp.factorial(k) * gb) * (-u) ** k
            s += c
            if k > 50 and abs(c) < mp.mpf(10) ** (-45) * (1 + abs(s)):
                break
        return s

    def F_int(b, u):
        f = lambda t: t ** (b - 1) * mp.exp(-t - u * t ** (-1j * al))
        return mp.quad(f, [0, 0.5, 2, 8, 30, 120, mp.inf]) / mp.gamma(b)

    def F_watson(b, u, J=3):
        s = mp.mpc(0)
        for j in range(J):
            s += (-1) ** j * mp.gamma(g * (b + j)) / mp.factorial(j) * u ** (-g * (b + j))
        return g * s / mp.gamma(b)

    print("γ = i/α = %s ; s_0^A = γβ_A = %s ; s_0^P = γβ_P = %s" % (mp.nstr(g, 8), mp.nstr(g * betas["A"], 8), mp.nstr(g * betas["P"], 8)))
    for name, b in betas.items():
        print("Kopf %s (β = %s)" % (name, mp.nstr(b, 6)))
        for u in [mp.mpc(0.7, 0.3), mp.mpc(-1.5, 2), mp.mpc(3, -1), mp.mpc(6, 0), mp.mpc(-4, -3)]:
            fr = F_reihe(b, u); fi = F_int(b, u)
            print("  u = %-14s Reihe %-28s Integral %-28s |Δ|/|F| = %s" % (mp.nstr(u, 3), mp.nstr(fr, 12), mp.nstr(fi, 12), mp.nstr(abs(fr - fi) / abs(fr), 3)))
        for r in [20, 40, 80]:
            for ph in [-30, 0, 30]:
                u = r * mp.expjpi(ph / 180.0)
                fr = F_reihe(b, u, K=20000); fw = F_watson(b, u)
                print("  Potenzsektor |u| = %d, arg = %+d°: Reihe %-26s Watson(3 Glieder) %-26s rel. Abw. %s" % (r, ph, mp.nstr(fr, 8), mp.nstr(fw, 8), mp.nstr(abs(fr - fw) / abs(fr), 3)))


if __name__ == "__main__":
    main()
