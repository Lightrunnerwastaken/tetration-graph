"""k34 — Asymptotik der Eckfunktionen Y (adjungiert) und P (primal): Potenzentwicklung im homogenen Keil aus den Residuen der Mellin–Barnes-Darstellung,
Kleinheit von E_A Y bzw. E_P P, exponentielle Terme (Sattelpunkt) und Lage der Keile (Spiralen).

Y(u) = (1/2πi)∫Γ(s)Γ(β_A + iαs)u^{−s}ds/Γ(β_A),  β_A = ½ − 2iα;  rechte Pole s_m = 2 + (2m+1)·2πi/L, Residuen-Asymptotik
   Y(u) ~ −Σ_m Res_{s_m} = Σ_m (−1)^{m+1}... (hier numerisch: r_m = −Γ(s_m)·(−1)^m/(m!·iα·Γ(β_A))).
P(u): Γ(s)Γ(β_P + iαs)/Γ(β_P), β_P = ½ + iα, Pole t_m = −1 + (2m+1)·2πi/L.
Sattel: Y ~ c·(s*)^{β−½} e^{−a s*}, s* = (u (iα)^{−iα})^{κ₁} (Zweig je Sektor).
Aufruf: python k34_asymptotik.py
"""
import sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 80
L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
ALPHA = L / (4 * mp.pi)
MU = mp.sqrt(L)
A_EXP = 1 + 1j * ALPHA
KAPPA = 1 / A_EXP
K = 1400


def koeff(beta):
    g0 = mp.gamma(beta)
    return [mp.gamma(beta - 1j * ALPHA * k) / (mp.factorial(k) * g0) for k in range(K)]


def ev(c, u):
    s = mp.mpc(0); p = mp.mpc(1)
    for x in c:
        s += x * p; p *= -u
    return s


def residuen(beta, polfam, M):
    """Pole von Γ(β + iαs): β + iαs = −m ⇒ s = (−m − β)/(iα). Residuum von Γ(s)Γ(β+iαs)u^{−s}/Γ(β) dort: Γ(s_m)(−1)^m/(m! iα Γ(β)) u^{−s_m}.
    Rechts schliessen ⇒ Beitrag −Res."""
    out = []
    for m in range(M):
        sm = (-m - beta) / (1j * ALPHA)
        r = -mp.gamma(sm) * (-1) ** m / (mp.factorial(m) * 1j * ALPHA * mp.gamma(beta))
        out.append((sm, r))
    return out


def main():
    for name, beta, op in [("Y (adjungiert)", 0.5 - 2j * ALPHA, "A"), ("P (primal)", 0.5 + 1j * ALPHA, "P")]:
        c = koeff(beta)
        res = residuen(beta, None, 4)
        print("=== %s: β = %s; rechte Pole: %s" % (name, mp.nstr(beta, 6), ", ".join(mp.nstr(s, 6) for s, r in res[:3])))
        print("    Residuen r_m: %s" % ", ".join(mp.nstr(r, 5) for s, r in res[:3]))
        for R in [20, 40, 80]:
            row = []
            for ang in range(-180, 180, 15):
                u = R * mp.expjpi(ang / 180)
                y = ev(c, u)
                asy = sum(r * u ** (-s) for s, r in res)
                if op == "A":
                    h = MU * y + ev(c, u / MU) / MU
                else:
                    h = ev(c, MU * u) + MU * y
                ra = abs(y - asy) / abs(y)
                rh = abs(h) / abs(y)
                row.append("%d:%s/%s" % (ang, ("%.0e" % float(ra)) if ra < 10 else ">", ("%.0e" % float(rh)) if rh < 10 else ">"))
            print("  R = %d: |Y−Asy|/|Y| / |E y|/|y|: %s" % (R, " ".join(row)))
        # Sattelpunkt-Kontrolle im Wachstumssektor
        print("  Sattel: log|y| gegen Re(−a s*) + Re((β−½) log s*), s* = (u (iα)^{−iα})^{κ₁}, Zweige k = −1, 0, 1 (u e^{2πik}):")
        for ang in [150, 180, -150, -120, 120]:
            for R in [40, 80]:
                u = R * mp.expjpi(ang / 180)
                y = ev(c, u)
                vals = []
                for kz in [-1, 0, 1]:
                    lu = mp.log(u) + 2j * mp.pi * kz
                    ss = mp.exp(KAPPA * (lu - 1j * ALPHA * mp.log(1j * ALPHA)))
                    vals.append(float(mp.re(-A_EXP * ss + (beta - 0.5) * mp.log(ss))))
                print("    arg u = %4d, R = %d: log|y| = %.3f; Sattel k=−1/0/1: %s" % (ang, R, float(mp.log(abs(y))), " / ".join("%.3f" % v for v in vals)))


if __name__ == "__main__":
    main()
