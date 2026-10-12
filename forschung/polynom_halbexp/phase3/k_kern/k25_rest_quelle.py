"""k25 — a2-Vortest: Ist der Abschneiderest T der adjungierten Kopfgleichung auf der Kegelskala klein gegen die Quelle 𝒰?

Exakte Identität (W §5, hier für den Kopf, u = e_n): H(ζ) = Σ_j ζ^{−j−1}⟨Λ, J_∞x^j⟩ = 𝒰(ζ) + T(ζ), 𝒰(ζ) = (ρ/ζ)^n/ζ, T(ζ) = Σ_{j>n} γ̃_j(ρ/ζ)^j/ζ,
γ̃_j = Σ_{k≤n} ℓ_k J̃_∞[k][j] (Überlauf des adjungierten Kopfes). Ist T ≪ 𝒰 im äusseren Keil der Kegelskala, so ist (P1) eine Gleichung mit BEKANNTER
Quelle, und die Mellin-Lösung ŷ ∝ 𝒰̂(s)/(1 + μ^{s−2}) ist explizit (Mellin–Barnes, Gammafunktionen).
Gemessen: |T/𝒰| auf Kreisen |ζ − L| = A·ρn^{−Re a}, A = 1, 4, 16, in 24 Richtungen. Aufruf: python k25_rest_quelle.py 300
"""
import cmath
import math
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from flint import arb, acb, arb_mat, arb_poly, ctx  # noqa: E402
from k_kern import J_tilde, kubik_koeff, Lacb, L_RE, L_IM  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
RE_A = 0.8935861640191741


def main():
    n = int(sys.argv[1])
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 250
    bits = int(1.3 * n) + 400
    J = J_tilde("kubisch", n, bits=bits)
    ctx.prec = bits
    e = arb_mat(n + 1, 1)
    e[n, 0] = arb(1)
    l = J.transpose().solve(e, algorithm="approx")
    lf = [l[k, 0] for k in range(n + 1)]
    # Spalten j > n von J̃_∞ in den Zeilen k ≤ n: [x^k] f^j · ρ^{k−j}  (m x^j trägt für j > n nichts zu k ≤ n bei)
    f = arb_poly([arb(x) for x in kubik_koeff()])
    rho = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
    pw = arb_poly([1])
    for _ in range(n):
        pw = (pw * f).truncate(n + 1)
    gam = []
    for j in range(n + 1, n + M + 1):
        pw = (pw * f).truncate(n + 1)
        s = arb(0)
        for k in range(n + 1):
            ck = pw[k] if k <= pw.degree() else arb(0)
            s += lf[k] * ck * rho ** k / rho ** j
        gam.append(s)
    print("n = %d: |γ̃_j| bei j = n+1, n+50, n+%d: %.2e, %.2e, %.2e" % (n, M, abs(float(gam[0].mid())), abs(float(gam[49].mid())), abs(float(gam[-1].mid()))), flush=True)
    L = Lacb()
    rhoc = acb(rho)
    sc = float(rho.mid()) * n ** (-RE_A)
    for A in [1.0, 4.0, 16.0]:
        row = []
        for th in range(-60, 330, 15):
            z = L + acb(sc * A) * acb(cmath.exp(1j * math.radians(th)))
            q = rhoc / z
            U = q ** n / z
            T = acb(0)
            qj = q ** (n + 1)
            for g in gam:
                T += acb(g) * qj
                qj = qj * q
            T = T / z
            Lz = acb(0)
            for k in range(n, -1, -1):
                Lz = Lz * q + acb(lf[k])
            fpz = f.derivative()
            mz = acb(fpz(f(z)) if False else 0)
            row.append((th, abs(complex(((U + T) / U).mid())), abs(complex((Lz / z / U).mid()))))
        print("  A = %4.1f: |H/𝒰| = |1 + T/𝒰|: " % A + " ".join("%d:%.2g" % (th, v) for th, v, w in row), flush=True)
        print("            |Λ/(ζ𝒰)|:         " + " ".join("%d:%.2g" % (th, w) for th, v, w in row), flush=True)


if __name__ == "__main__":
    main()
