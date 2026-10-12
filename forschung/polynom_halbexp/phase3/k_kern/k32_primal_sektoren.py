"""k32 — Sektorstruktur des primalen Kandidaten P(u) = Σ_k Γ(½ + iα − iαk)(−u)^k/(k! Γ(½ + iα)) gegen die gemessenen primalen Kopf-Defekte (k23).
Primaler Defekt |δ∘f + mδ|/|mδ| lokal = |P(μu) + μP(u)|/|μP(u)|, u = λ_P ζ_x, ζ_x = (x − L)/d_x, d_x = −(L²/ρ) n^{−a}, λ_P = e^{iπa} λ_A (k31).
Aufruf: python k32_primal_sektoren.py
"""
import sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 60
L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
ALPHA = L / (4 * mp.pi)
MU = mp.sqrt(L)
A_EXP = 1 + 1j * ALPHA
LAM_A = mp.mpc(-0.203297, -1.406008)
LAM_P = mp.mpc(0.603424, 1.165020)
BP = 0.5 + 1j * ALPHA
G0 = mp.gamma(BP)
K = 700
CP = [mp.gamma(BP - 1j * ALPHA * k) / (mp.factorial(k) * G0) for k in range(K)]


def P(u):
    s = mp.mpc(0); p = mp.mpc(1)
    for c in CP:
        s += c * p; p *= -u
    return s


def main():
    print("λ_P/λ_A = %s, e^{iπa} = %s" % (mp.nstr(LAM_P / LAM_A, 8), mp.nstr(mp.exp(1j * mp.pi * A_EXP), 8)))
    for n in [300, 600]:
        drift = mp.im(A_EXP) * mp.log(n)
        off = -mp.pi - 2 * mp.arg(L) + drift + mp.arg(LAM_P)
        print("n = %d: arg u = θ %+.2f°" % (n, float(mp.degrees(off))))
        for A in [1.0, 4.0, 16.0]:
            row = []
            for th in range(-60, 330, 15):
                u = abs(LAM_P) * A * mp.expjpi((mp.radians(th) + off) / mp.pi)
                pu = P(u)
                d = float(abs(P(MU * u) + MU * pu) / abs(MU * pu))
                row.append("%d:%s" % (th, ("%.0e" % d) if d < 0.1 else ("%.1f" % d if d < 10 else ">10")))
            print("  A = %4.1f Modell primal: %s" % (A, " ".join(row)))
    print("Lokale Exponenten d log P/d log|u| (Vorhersage −(−1 + 2πi/L) = %s):" % mp.nstr(1 - 2j * mp.pi / L, 6))
    for argu in [120, 150, 180, 210, 240]:
        out = []
        for R in [8, 16, 32, 48]:
            u1 = R * mp.expjpi(argu / 180)
            sl = (mp.log(P(1.05 * u1)) - mp.log(P(u1))) / mp.log(1.05)
            out.append("R=%d: %s" % (R, mp.nstr(sl, 4)))
        print("  arg u = %d°: %s" % (argu, " | ".join(out)))
    print("log10|P(u)| auf Kreisen:")
    for R in [2, 8, 23, 40]:
        print("  R = %d: %s" % (R, " ".join("%d:%.1f" % (a, float(mp.log10(abs(P(R * mp.expjpi(a / 180)))))) for a in range(-180, 180, 20))))


if __name__ == "__main__":
    main()
