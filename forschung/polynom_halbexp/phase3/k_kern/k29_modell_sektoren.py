"""k29 — Sektorstruktur des Mellin–Barnes-Kandidaten Y(u) = Σ_k Γ(½ − iα(k+2))(−u)^k/(k! Γ(½ − 2iα)) und Vergleich mit den Kopf-Defekten (k23).

H(u) := μY(u) + μ^{−1}Y(u/μ)  (lokale adjungierte Kopfgleichung; Mellin: μ(1+μ^{s−2})Γ(s)Γ(½+iα(s−2)) = 2πμ^{s/2}Γ(s)/Γ(½−iα(s−2)))
       = (2π/Γ(½−2iα)) Σ_k (−u/√μ)^k / (k! Γ(½ + iα(k+2)))   — beide Reihen werden zur Kontrolle ausgewertet.
Defekt wie k23: |H|/|μY|. Abbildung der x-Richtung θ (x − L = r e^{iθ}) auf u (Gs Koordinate w = ρ/x, ζ_G = (w − w_L)/n^{−a}, u = λζ_G):
  arg u = θ + arg(−λ) − 2 arg L + Im(a)·ln n,  |u| = |λ|·A  (A wie in k23: |x − L| = A ρ n^{−Re a}).
Aufruf: python k29_modell_sektoren.py
"""
import sys

import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 60
L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
ALPHA = L / (4 * mp.pi)
MU = mp.sqrt(L)
LAM = mp.mpc(-0.203297, -1.406007)
A_EXP = 1 / (2j * mp.pi / (2j * mp.pi - L / 2))
G0 = mp.gamma(0.5 - 2j * ALPHA)


def koeff_Y(K):
    return [mp.gamma(0.5 - 1j * ALPHA * (k + 2)) / (mp.factorial(k) * G0) for k in range(K)]


def koeff_H(K):
    return [2 * mp.pi / G0 * (-1 / mp.sqrt(MU)) ** k / (mp.factorial(k) * mp.gamma(0.5 + 1j * ALPHA * (k + 2))) for k in range(K)]


CY = None
CH = None


def Y(u):
    s = mp.mpc(0)
    p = mp.mpc(1)
    for c in CY:
        s += c * p
        p *= -u
    return s


def Hreihe(u):
    s = mp.mpc(0)
    p = mp.mpc(1)
    for c in CH:
        s += c * p
        p *= u
    return s


def main():
    global CY, CH
    K = 700
    CY = koeff_Y(K)
    CH = koeff_H(K)
    # Kontrolle der Identität H = μY + Y(·/μ)/μ
    for u in [mp.mpc(1, 0.5), mp.mpc(-3, 4), mp.mpc(10, -7)]:
        h1 = MU * Y(u) + Y(u / MU) / MU
        h2 = Hreihe(u)
        print("Kontrolle u = %s: |H_def − H_Reihe|/|H| = %.1e" % (mp.nstr(u, 4), float(abs(h1 - h2) / abs(h1))))
    print("|λ| = %.4f, arg(−λ) = %.2f°, arg L = %.3f°" % (float(abs(LAM)), float(mp.degrees(mp.arg(-LAM))), float(mp.degrees(mp.arg(L)))))
    # Defekt-Tabelle wie k23 (n = 600, A = 4 und 16)
    for n in [300, 600]:
        drift = float(mp.im(A_EXP)) * mp.log(n)
        print("\nn = %d: arg u = θ %+.2f°" % (n, float(mp.degrees(mp.arg(-LAM)) - 2 * mp.degrees(mp.arg(L)) + mp.degrees(drift))))
        for A in [1.0, 4.0, 16.0]:
            row = []
            for th in range(-60, 330, 15):
                au = mp.radians(th) + mp.arg(-LAM) - 2 * mp.arg(L) + drift
                u = abs(LAM) * A * mp.expjpi(au / mp.pi)
                y = Y(u)
                h = MU * y + Y(u / MU) / MU
                d = float(abs(h) / abs(MU * y))
                row.append("%d:%s" % (th, ("%.0e" % d) if d < 0.1 else ("%.1f" % d if d < 10 else ">10")))
            print("  A = %4.1f Modell-Defekt |H|/|μY|: %s" % (A, " ".join(row)))
    # Potenzverhalten im homogenen Keil: lokale Steigung d log Y / d log u
    print("\nLokale Exponenten d log Y/d log|u| längs Strahlen (Modell), Vorhersage −s₀ = −(2 + 2πi/L) = %s" % mp.nstr(-(2 + 2j * mp.pi / L), 6))
    for argu in [0, 15, 30, 45, 60, -15]:
        out = []
        for R in [8, 16, 32, 48]:
            u1 = R * mp.expjpi(argu / 180)
            u2 = 1.05 * u1
            sl = (mp.log(Y(u2)) - mp.log(Y(u1))) / mp.log(1.05)
            out.append("R=%d: %s" % (R, mp.nstr(sl, 4)))
        print("  arg u = %+d°: %s" % (argu, " | ".join(out)))
    # Betrag längs Kreisen
    print("\nlog10|Y(u)| auf Kreisen |u| = R (Modell):")
    for R in [2, 8, 23, 40]:
        print("  R = %d: %s" % (R, " ".join("%d:%.1f" % (a, float(mp.log10(abs(Y(R * mp.expjpi(a / 180)))))) for a in range(-180, 180, 20))))


if __name__ == "__main__":
    main()
