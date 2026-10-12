"""k28 — Gs normierte Eckprofile S_n (n = 50…1000, je 40 Taylor-Koeffizienten) extrapoliert nach n → ∞ (Ansatz c + A n^{−a} + B n^{−1} + C n^{−2a},
a = 1/κ₁ komplex, wie G) und Vergleich mit dem Mellin–Barnes-Kandidaten c_k = Γ(β₀ − iαk)(−λ)^k/(k! Γ(β₀)), β₀ = ½ − 2iα, α = L/(4π).
λ aus c_1. Ausgabe: relative Abweichung je k, freier β-Fit am Grenzprofil, Vergleich mit Gs eigener Extrapolation (shape_limit_first_5).
Aufruf: python k28_extrapol_vergleich.py
"""
import json
import os
import sys

import mpmath as mp
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 30
HIER = os.path.dirname(os.path.abspath(__file__))
WURZEL = os.path.join(HIER, "..", "..")
L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
ALPHA = L / (4 * mp.pi)
KAPPA = 2j * mp.pi / (2j * mp.pi - L / 2)
A_EXP = 1 / KAPPA


def lade():
    d = json.load(open(os.path.join(WURZEL, "g_gpt_hilfssaetze", "runde5_szego", "spiral_analysis.json"), encoding="utf-8"))
    m = d["models"]["kubisch"]
    prof = {p["n"]: np.array([complex(float(a), float(b)) for a, b in p["natural_shape_coefficients_first_40"]]) for p in m["profiles"]}
    return prof, m


def extrapol(prof, ns, terme):
    """Kleinste Quadrate je Koeffizient: c_k(n) = c_k∞ + Σ A_t n^{−e_t}."""
    a = complex(A_EXP)
    ex = {"a": a, "1": 1.0, "2a": 2 * a, "a+1": a + 1, "2": 2.0}
    M = np.array([[1.0] + [n ** (-ex[t]) for t in terme] for n in ns], dtype=complex)
    Y = np.array([prof[n] for n in ns])
    sol, *_ = np.linalg.lstsq(M, Y, rcond=None)
    res = Y - M @ sol
    return sol[0], np.abs(res).max()


def modell(beta, lam, K):
    lg0 = mp.loggamma(beta)
    return np.array([complex(mp.exp(mp.loggamma(beta - 1j * ALPHA * k) - lg0 - mp.loggamma(k + 1)) * (-lam) ** k) for k in range(K)])


def main():
    prof, m = lade()
    ns_all = sorted(prof)
    print("verfügbare n:", ns_all)
    beta0 = 0.5 - 2j * ALPHA
    for ns, terme in [([200, 300, 400, 500, 600], ["a", "1", "2a"]),
                      ([400, 500, 600, 700, 800, 900, 1000], ["a", "1", "2a"]),
                      ([300, 400, 500, 600, 700, 800, 900, 1000], ["a", "1", "2a", "a+1"]),
                      ([500, 600, 700, 800, 900, 1000], ["a", "1"])]:
        ns = [n for n in ns if n in prof]
        S, rmax = extrapol(prof, ns, terme)
        lam = -S[1] * complex(mp.gamma(beta0) / mp.gamma(beta0 - 1j * ALPHA))
        Mo = modell(beta0, lam, len(S))
        rel = np.abs(S - Mo) / np.maximum(np.abs(S), 1e-300)
        print("\nExtrapolation n ∈ %s, Terme %s (max. Fit-Rest %.1e):" % (ns, terme, rmax))
        print("  λ (aus c_1) = %.6f%+.6fi, |λ| = %.6f, arg λ = %.3f°" % (lam.real, lam.imag, abs(lam), np.degrees(np.angle(lam))))
        print("  S_∞ (k = 0..6): " + " ".join("%.5f%+.5fi" % (x.real, x.imag) for x in S[:7]))
        print("  Modell (k = 0..6): " + " ".join("%.5f%+.5fi" % (x.real, x.imag) for x in Mo[:7]))
        print("  rel. Abweichung k = 2..39: " + " ".join("%.1e" % x for x in rel[2:40:3]))
        print("  max rel. Abw. k ≤ 10: %.2e, k ≤ 20: %.2e, k ≤ 39: %.2e" % (rel[2:11].max(), rel[2:21].max(), rel[2:40].max()))
        # Wiener-Norm-Abstand bei Radius 1.05
        R = 1.05
        w = R ** np.arange(len(S))
        print("  Wiener-Abstand (R = 1.05) Modell − S_∞: %.2e  (Norm S_∞: %.4f)" % (np.sum(np.abs(S - Mo) * w), np.sum(np.abs(S) * w)))
    fe = m["finite_extrapolations"]
    for e in fe:
        print("\nGs eigene Extrapolation %d–%d, erste 5:" % (e["low"], e["high"]), e["shape_limit_first_5"])


if __name__ == "__main__":
    main()
