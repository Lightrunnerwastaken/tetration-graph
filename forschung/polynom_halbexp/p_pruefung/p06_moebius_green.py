"""Strang P, Prüfung 6: E Lemma M2/M3 (Möbius-Drehung) unabhängig nachgerechnet.
Green-Funktion von C \\ (D̄(0,1) ∪ D̄(τ, √(1+τ²))) direkt (Möbius auf Keil + Potenz), Wert am Pol 1/τ gegen
E's Formel e^{−g} = sin(πc/(2ω)) / sin(π(π−2c)/(2ω)), ω = π − c/2; und Schwelle ζ_p ∈ Ḡ' ⇔ c > π/3.
"""
import sys

import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 30


def green_setup(c):
    tau = mp.tan(c / 2)
    R0 = mp.sqrt(1 + tau ** 2)
    P, Pb = mp.mpc(0, 1), mp.mpc(0, -1)

    def M(z):
        return (z - P) / (z - Pb)
    t1 = mp.arg(M(mp.mpf(-1)))           # Bild Einheitskreis (links)
    t2 = mp.arg(M(tau + R0))             # Bild äusserer Kreis
    beta = t2
    om = t1 - t2
    if om < 0:
        om += 2 * mp.pi

    def W(z):
        w = M(z) * mp.expj(-beta)
        ang = mp.arg(w)
        if ang < -mp.mpf(10) ** -20:
            ang += 2 * mp.pi
        return abs(w) ** (mp.pi / om) * mp.expj(ang * mp.pi / om)
    Winf = W(mp.mpf(10) ** 25)

    def g(z):
        w = W(z)
        return mp.log(abs((w - mp.conj(Winf)) / (w - Winf)))
    return tau, R0, om, g


for c in [mp.mpf("0.4"), mp.mpf("0.6"), mp.mpf("0.8"), mp.mpf("1.0"), mp.pi / 3 - mp.mpf("0.01")]:
    tau, R0, om, g = green_setup(c)
    zp = 1 / tau
    gp = g(zp)
    omega = mp.pi - c / 2
    formel = mp.sin(mp.pi * c / (2 * omega)) / mp.sin(mp.pi * (mp.pi - 2 * c) / (2 * omega))
    inG = abs(zp - tau) < R0
    print(f"c = {mp.nstr(c, 6)}: ω(Keil) = {mp.nstr(om, 10)} (π−c/2 = {mp.nstr(omega, 10)}), e^(−g(ζ_p)) = {mp.nstr(mp.exp(-gp), 10)}, "
          f"Formel M3 = {mp.nstr(formel, 10)}, ζ_p ∈ D(a0,R0): {inG}")
print("Schwelle: (1−τ²)² = τ²(1+τ²) bei τ = 1/√3 ⇔ c = π/3:", mp.nstr(2 * mp.atan(1 / mp.sqrt(3)), 12), "=", mp.nstr(mp.pi / 3, 12))
