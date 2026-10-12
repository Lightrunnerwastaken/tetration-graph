"""k33 — Lokale L²-Normen der Eckprofile längs der Tangente des Kreises |x| = ρ in L (auf der Kegelskala = Kreis und Linse).
Primal: ‖c_n‖ ≈ |δ(L)| n^{−Re a/2} · sqrt(I_P(ψ_n)/(π|λ_P|)),  I_P(ψ) = ∫_ℝ |P(t e^{iψ})|² dt,  ψ_n = −90° + arg λ_P − arg L + Im(a) ln n.
Adjungiert (x-Norm ℓ²_ρ, also |w| = 1): analog mit Y, λ_A und ψ_n^A = arg(i w_L λ_A n^{a}) (Tangente von |w| = 1 in w_L).
Gemessen wird I(ψ; R) = ∫_{−R}^{R} — konvergiert es für R → ∞?  Vergleich mit dem gemessenen Quotienten ‖c_n‖/(|δ(L)|n^{−Re a/2}) aus k12/k30.
Aufruf: python k33_linienintegrale.py
"""
import sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 50
L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
ALPHA = L / (4 * mp.pi)
A_EXP = 1 + 1j * ALPHA
LAM_A = mp.mpc(-0.203297, -1.406008)
LAM_P = LAM_A * mp.exp(1j * mp.pi * A_EXP)
K = 600


def koeff(beta):
    g0 = mp.gamma(beta)
    return [mp.gamma(beta - 1j * ALPHA * k) / (mp.factorial(k) * g0) for k in range(K)]


CP = koeff(0.5 + 1j * ALPHA)
CY = koeff(0.5 - 2j * ALPHA)


def ev(c, u):
    s = mp.mpc(0); p = mp.mpc(1)
    for x in c:
        s += x * p; p *= -u
    return s


def linie(c, psi, R, h=0.25):
    e = mp.expjpi(psi / mp.pi)
    tot = mp.mpf(0)
    prof = []
    t = -R
    while t <= R + 1e-12:
        v = abs(ev(c, t * e)) ** 2
        w = h / 2 if abs(abs(t) - R) < 1e-12 else h
        tot += w * v
        if abs(t - round(t)) < 1e-9 and int(round(t)) % 8 == 0:
            prof.append((int(round(t)), float(mp.log10(abs(ev(c, t * e)) + mp.mpf(10) ** -300))))
        t += h
    return tot, prof


def main():
    for n in [400, 800, 1000, 10 ** 4, 10 ** 6]:
        drift = mp.im(A_EXP) * mp.log(n)
        psiP = -mp.pi / 2 + mp.arg(LAM_P) - mp.arg(L) + drift
        # adjungiert: w = ρ/x, Tangente von |w| = 1 in w_L: w − w_L = i w_L t; ζ_G = (w − w_L)/n^{−a}; u = λ_A ζ_G
        wL = abs(L) / L
        psiA = mp.arg(1j * wL * LAM_A) + drift
        print("n = %g: ψ_P = %.2f°, ψ_A = %.2f°" % (n, float(mp.degrees(psiP)), float(mp.degrees(psiA))))
        for name, c, psi, lam in [("primal", CP, psiP, LAM_P), ("adjung.", CY, psiA, LAM_A)]:
            out = []
            for R in [8, 16, 32]:
                I, prof = linie(c, psi, R)
                out.append("R=%d: I=%.4g, Quot=%.4f" % (R, float(I), float(mp.sqrt(I / (mp.pi * abs(lam))))))
            print("   %s: %s | log10|·| längs t: %s" % (name, " ; ".join(out), prof))


if __name__ == "__main__":
    main()
