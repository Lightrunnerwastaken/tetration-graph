"""k37 — Modellvorhersage (Eckprofil Y) für dünnere Linsen bei erreichbaren n, um den Spiral-Mechanismus am echten Kopf prüfen zu können (k38).
Linse s(τ) = S(1 − exp(−(τ/τ₁)^p)) mit p ∈ {1.3, 1.6, 2.0, 2.5}, S = 0.0166, τ₁ = 12°.  Ausgabe wie k36 (max|Y| längs γ, A ≤ 60).
Aufruf: python k37_duenne_linse_modell.py
"""
import math
import sys
import mpmath as mp
sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import k36_linse_asymptotik as M  # noqa: E402
sys.stdout.reconfigure(encoding="utf-8")


def main():
    for p in [1.3, 1.6, 2.0, 2.5]:
        print("p = %.1f:" % p)
        for n in [300, 600, 1000, 2000, 5000]:
            pr = M.profil(n, lambda t: M.s_tang(t, p=p), 60.0, nA=50)
            m = max(pr, key=lambda t: t[1])
            l2 = math.sqrt(sum(((pr[i][0] - pr[i - 1][0]) * pr[i][1] ** 2) for i in range(1, len(pr))))
            print("   n = %5d: max|Y| = %.3g bei A = %.1f; (∫|Y|²dA)^{1/2} = %.3g" % (n, m[1], m[0], l2), flush=True)


if __name__ == "__main__":
    main()
