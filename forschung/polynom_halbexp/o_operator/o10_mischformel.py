"""O10 — Prüfung der Mischformel (Satz O3): Π_N g = Y_N T⁽⁰⁾_N g + (I − Y_N) T^{(a₀)}_N g  (exakt, jedes f mit J = U S + M).
Kubisches Modell, g = e^x und g = 1/(ζ − x); Π_N g direkt = J_N^{-1} P_N[g∘f + m g] gegen Formel. Rechnung ganz in arb (unskaliert)."""
import sys
from flint import arb, acb, arb_mat, acb_mat, arb_poly, acb_poly, ctx
import o_kern as K
sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 60
ctx.prec = 3 * N + 300
fk, mk, kp = K.modell("kubisch", N)
f = arb_poly(fk); a0 = fk[0]; u = f - arb_poly([a0]); n = N + 1
pw, upw = [arb_poly([1])], [arb_poly([1])]
for _ in range(N):
    pw.append((pw[-1] * f).truncate(n)); upw.append((upw[-1] * u).truncate(n))
from math import comb, factorial
J = arb_mat(n, n); M = arb_mat(n, n)
for j in range(n):
    for k in range(n):
        J[k, j] = (pw[j][k] if k <= pw[j].degree() else 0) + (mk[k - j] if k >= j else 0)
        M[k, j] = mk[k - j] if k >= j else 0
Y = J.solve(M, algorithm="approx")
Sinv = arb_mat(n, n)      # a₀-Koeffizienten -> 0-Koeffizienten: (x−a₀)^j = Σ_k C(j,k)(−a₀)^{j−k} x^k
for j in range(n):
    for k in range(j + 1):
        Sinv[k, j] = arb(comb(j, k)) * (-a0) ** (j - k)
I = arb_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)])
for name in ("exp", "pol"):
    if name == "exp":
        jet0 = [arb(1) / factorial(k) for k in range(n)]
        jeta = [a0.exp() / factorial(k) for k in range(n)]
        # g∘f + m g bis Grad N: e^{f} Reihe
        gf = arb_poly([0]); term = arb_poly([1]); 
        # exp(f) = e^{a0} exp(u): Reihe über u^j/j!
        gf = sum((upw[j] * (a0.exp() / factorial(j)) for j in range(n)), arb_poly([0]))
    else:
        z = arb("2.3")
        jet0 = [1 / z ** (k + 1) for k in range(n)]
        jeta = [1 / (z - a0) ** (k + 1) for k in range(n)]
        gf = sum((upw[j] * jeta[j] for j in range(n)), arb_poly([0]))
    mg = (arb_poly(mk) * arb_poly(jet0)).truncate(n)
    rhs = arb_mat(n, 1, [(gf[k] if k <= gf.degree() else 0) + (mg[k] if k <= mg.degree() else 0) for k in range(n)])
    direkt = J.solve(rhs, algorithm="approx")
    T0 = arb_mat(n, 1, jet0)
    Ta = Sinv * arb_mat(n, 1, jeta)
    formel = Y * T0 + (I - Y) * Ta
    d = max(abs(float((direkt[k, 0] - formel[k, 0]).mid())) * 1.3745570107437075 ** k for k in range(n))
    s = max(abs(float(direkt[k, 0].mid())) * 1.3745570107437075 ** k for k in range(n))
    print(f"N={N}, g={name}: max_k |Π_N g − Formel|_k ρ^k = {d:.2e}  (Größe {s:.2e})")
