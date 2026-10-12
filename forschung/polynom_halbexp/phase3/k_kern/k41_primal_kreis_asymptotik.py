"""k41 — Primales Eckprofil längs des KREISES |x| = ρ (die Unbekannten-Norm ℓ²_ρ ist an den Kreis gebunden — keine Linse hilft) für grosse n.
P(e^{iπa}u) mit u = λ_A ζ_G; Kreis bei L in Gs Koordinate exakt (w = ρ/x, |w| = 1). Beide Hälften: Γin (τ > 0) und Γout (τ < 0).
Frage: wann verlässt eine Hälfte den Abklingsektor von F_P (dann ist ‖c_n‖ nicht mehr kegel-lokal)?  Gleichzeitig adjungiert längs Γout (ohne Linse).
Aufruf: python k41_primal_kreis_asymptotik.py
"""
import math
import sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 160
L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
RHO = abs(L)
AL = L / (4 * mp.pi)
A = 1 + 1j * AL
LAMA = mp.mpc(-0.203297, -1.406008)
LAMP = LAMA * mp.exp(1j * mp.pi * A)
K = 1400


def koeff(beta):
    g0 = mp.gamma(beta)
    return [mp.gamma(beta - 1j * AL * k) / (mp.factorial(k) * g0) for k in range(K)]


CP = koeff(0.5 + 1j * AL)
CY = koeff(0.5 - 2j * AL)


def ev(c, u):
    s = mp.mpc(0); p = mp.mpc(1)
    for x in c:
        s += x * p; p *= -u
    return s


def main():
    wL = RHO / L
    for n in [mp.mpf(10) ** e for e in (3, 4, 6, 8, 9, 10, 12, 14)]:
        dn = mp.power(n, -A)
        scale = RHO * abs(dn)
        out = []
        for side, sg in [("Γin", 1), ("Γout", -1)]:
            vals = []
            for Aa in [4, 8, 16, 32, 64, 100]:
                tau = sg * Aa * scale / RHO
                x = RHO * mp.expjpi((mp.arg(L) + tau) / mp.pi)
                zG = (RHO / x - wL) / dn
                # primal in x-Koordinate: ζ_x = (x−L)/d_x mit d_x = −(L²/ρ) n^{−a}; lineare Näherung ζ_x ≈ ζ_G
                zx = (x - L) / (-(L ** 2 / RHO) * dn)
                vp = float(mp.log10(abs(ev(CP, LAMP * zx))))
                va = float(mp.log10(abs(ev(CY, LAMA * zG))))
                vals.append("A%d:%+.1f/%+.1f" % (Aa, vp, va))
            argp = float(mp.degrees(mp.arg(LAMP * (mp.expjpi((mp.arg(L) + sg * 1e-9) / mp.pi) * RHO - L) / (-(L ** 2 / RHO) * dn))))
            out.append("%s (arg u_P = %.1f°): %s" % (side, argp, " ".join(vals)))
        print("n = 1e%d: log10|P| / log10|F_A| längs des Kreises:\n   %s\n   %s" % (int(mp.log10(n)), out[0], out[1]), flush=True)


if __name__ == "__main__":
    main()
