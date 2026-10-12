"""k31 — Extrapolation der primalen (und adjungierten) Eckprofile aus k30 (n = 200…800, 30 Taylor-Koeffizienten) nach n → ∞ und Vergleich mit
P(u) = Σ_k Γ(½ + iα − iαk)(−u)^k/(k! Γ(½ + iα))  (primal)  bzw.  Y(u) = Σ_k Γ(½ − 2iα − iαk)(−u)^k/(k! Γ(½ − 2iα))  (adjungiert).
Skala λ aus c_1. Zusätzlich freier β-Fit am Grenzprofil (prüft, ob β_P = ½ + iα wirklich der Wert ist).
Aufruf: python k31_primal_extrapol.py
"""
import glob
import os
import re
import sys

import mpmath as mp
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 30
HIER = os.path.dirname(os.path.abspath(__file__))
L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
ALPHA = L / (4 * mp.pi)
A_EXP = complex(1 / (2j * mp.pi / (2j * mp.pi - L / 2)))


def modell(beta, lam, K):
    lg0 = mp.loggamma(beta)
    return np.array([complex(mp.exp(mp.loggamma(beta - 1j * ALPHA * k) - lg0 - mp.loggamma(k + 1)) * (-lam) ** k) for k in range(K)])


def extrapol(prof, ns, terme):
    ex = {"a": A_EXP, "1": 1.0, "2a": 2 * A_EXP, "a+1": A_EXP + 1}
    M = np.array([[1.0] + [n ** (-ex[t]) for t in terme] for n in ns], dtype=complex)
    Y = np.array([prof[n] for n in ns])
    sol, *_ = np.linalg.lstsq(M, Y, rcond=None)
    return sol[0]


def fit_beta(S, kmax):
    """β frei, λ aus c_1: minimiere Σ_{k=2..kmax} |c_k − Modell_k|² R^{2k}."""
    def res(b):
        lam = -S[1] * complex(mp.gamma(b) / mp.gamma(b - 1j * ALPHA))
        Mo = modell(b, lam, kmax + 1)
        return (Mo[2:] - S[2:kmax + 1]) * 1.05 ** np.arange(2, kmax + 1)
    b = complex(0.5)
    for b0 in [0.5 + 1j * complex(ALPHA), 0.5 - 2j * complex(ALPHA), 0.5]:
        bb = complex(b0)
        for _ in range(30):
            r0 = res(bb)
            h = 1e-7
            J = (res(bb + h) - r0) / h
            st = -np.vdot(J, r0) / np.vdot(J, J)
            bb += st
            if abs(st) < 1e-13:
                break
        if np.linalg.norm(res(bb)) < np.linalg.norm(res(b)):
            b = bb
    return b, np.linalg.norm(res(b))


def main():
    prof_p, prof_a = {}, {}
    for f in glob.glob(os.path.join(HIER, "k30_profile_*.npz")):
        n = int(re.findall(r"_(\d+)\.npz", f)[0])
        z = np.load(f)
        prof_p[n] = z["primal"]
        prof_a[n] = z["adj"]
    ns = sorted(prof_p)
    print("n:", ns)
    for name, prof, beta in [("primal", prof_p, 0.5 + 1j * ALPHA), ("adjungiert", prof_a, 0.5 - 2j * ALPHA)]:
        print("\n=== %s, Modell β = %s ===" % (name, mp.nstr(beta, 8)))
        for nsel, terme in [([300, 400, 500, 600, 800], ["a", "1", "2a"]), ([200, 300, 400, 500, 600, 800], ["a", "1", "2a", "a+1"]),
                            ([400, 500, 600, 800], ["a", "1"])]:
            nsel = [n for n in nsel if n in prof]
            S = extrapol(prof, nsel, terme)
            lam = -S[1] * complex(mp.gamma(beta) / mp.gamma(beta - 1j * ALPHA))
            Mo = modell(beta, lam, len(S))
            w = 1.05 ** np.arange(len(S))
            rel = np.abs(S - Mo) / np.abs(S)
            print(" n ∈ %s, %s: λ = %.6f%+.6fi (|λ| = %.5f, arg %.3f°); Wiener-Abstand %.2e (Norm %.4f); max rel k≤10: %.1e" % (
                nsel, terme, lam.real, lam.imag, abs(lam), np.degrees(np.angle(lam)), np.sum(np.abs(S - Mo) * w), np.sum(np.abs(S) * w), rel[2:11].max()))
            print("   S_∞ k=1..5: " + " ".join("%.5f%+.5fi" % (x.real, x.imag) for x in S[1:6]))
            print("   Mod k=1..5: " + " ".join("%.5f%+.5fi" % (x.real, x.imag) for x in Mo[1:6]))
            b, r = fit_beta(S, 12)
            print("   freier β-Fit (k ≤ 12): β = %.6f%+.6fi  (Modell %.6f%+.6fi), Rest %.1e" % (b.real, b.imag, complex(beta).real, complex(beta).imag, r))
        # Konvergenz bei endlichem n
        S = extrapol(prof, [n for n in [300, 400, 500, 600, 800] if n in prof], ["a", "1", "2a"])
        lam = -S[1] * complex(mp.gamma(beta) / mp.gamma(beta - 1j * ALPHA))
        Mo = modell(beta, lam, len(S))
        w = 1.05 ** np.arange(len(S))
        print("  ‖S_n − Modell‖_W·n^{0.894}: " + " ".join("%d:%.4f" % (n, np.sum(np.abs(prof[n] - Mo) * w) * n ** A_EXP.real) for n in ns))


if __name__ == "__main__":
    main()
