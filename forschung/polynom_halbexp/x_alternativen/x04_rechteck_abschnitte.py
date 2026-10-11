"""X04 — Sind die Eckmoden Artefakte des QUADRATISCHEN Abschnitts? Rechteckige Abschnitte T_M J_∞|_{P_N} (M ≥ N), kubisches Modell.

Idee (Literatur: „uneven/rectangular sections“, Hansen JFA 254 (2008), JAMS 24 (2011); Lindner, Infinite Matrices and their Finite Sections
(2006); Heinemeyer–Lindner–Potthast SINUM 46 (2008)): quadratische Abschnitte erzeugen Fast-Kerne an der Abschneidekante, die der
Operator selbst nicht hat; rechteckige Abschnitte (mehr Gleichungen als Unbekannte) messen den Injektionsmodul von J_∞ auf P_N.
Hier: kubisches Modell (Strang D §9.3), f(x) = x + (x² − 2Re(L)x + |L|²)(α + βx), f'(L) = √L — J_∞u = u∘f + (f'∘f)u bildet P_N in
P_{3N} ab, also ist T_{3N}J_∞|_{P_N} EXAKT die volle Abbildung. Norm ℓ²_ρ (Gewichte ρ^k in Bild und Urbild).
Ausgabe: die kleinsten Singulärwerte σ_min von T_M J_∞|_{P_N} für M = N, N+2, N+5, N+10, N+20, N+50, 1.25N, 3N.
Methode: G_M = ÃᵀÃ in arb (hohe Präzision), G_M⁻¹ in arb, dann float64-Eigenwerte (die GROSSEN Eigenwerte von G⁻¹ = 1/σ_min² sind in
double gut konditioniert).  Aufruf: python x04_rechteck_abschnitte.py N
"""
import sys
import time

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
ctx.prec = int(2.6 * N) + 400
t0 = time.time()

Lr = arb("0.31813150520476413531265425158766451720351761387139986692237")
Li = arb("1.33723570143068940890116214319371061253950213846051241687")
R2 = Lr * Lr + Li * Li
rho = R2.sqrt()
# μ = √L (Hauptzweig); (μ−1)/(L−L̄) = (μ−1)/(2i Li)
modL = rho
mu_r = ((modL + Lr) / 2).sqrt()
mu_i = ((modL - Lr) / 2).sqrt()                  # Im L > 0 ⇒ Im √L > 0
# rL = (μ − 1)/(2i·Li) = (μ_i − i(μ_r − 1)) / (2 Li)  ⇒ Re = μ_i/(2Li), Im = −(μ_r − 1)/(2Li)
rL_re = mu_i / (2 * Li)
rL_im = -(mu_r - 1) / (2 * Li)
beta = rL_im / Li
alpha = rL_re - beta * Lr
c = [R2 * alpha, 1 + R2 * beta - 2 * Lr * alpha, alpha - 2 * Lr * beta, beta]
f = arb_poly(c)
fp = arb_poly([c[1], 2 * c[2], 3 * c[3]])
m = fp(f)                                      # f'∘f, Grad 6
print(f"Koeffizienten f: {[float(x.mid()) for x in c]}, f(0) = {float(c[0].mid()):.5f}", flush=True)

MMAX = 3 * N
# Spalten: x^j∘f + m·x^j, skaliert Ã[k][j] = ρ^{k−j}·[x^k](…)
A = arb_mat(MMAX + 1, N + 1)
pot = arb_poly([1])
rk = [arb(1)]
for k in range(1, MMAX + 1):
    rk.append(rk[-1] * rho)
for j in range(N + 1):
    col = pot + m * arb_poly([0] * j + [1])
    inv_rj = 1 / rk[j]
    for k in range(min(col.degree(), MMAX) + 1):
        A[k, j] = col[k] * rk[k] * inv_rj
    pot = pot * f
print(f"Matrix aufgebaut ({time.time() - t0:.0f}s), prec = {ctx.prec}", flush=True)

Ms = sorted(set([N, N + 2, N + 5, N + 10, N + 20, N + 50, int(1.25 * N), MMAX]))
Ms = [M for M in Ms if M <= MMAX]
G = arb_mat(N + 1, N + 1)
start = 0
for M in Ms:
    blk = arb_mat(M + 1 - start, N + 1)
    for i in range(start, M + 1):
        for j in range(N + 1):
            blk[i - start, j] = A[i, j]
    G = G + blk.transpose() * blk
    start = M + 1
    Gi = G.inv()
    Gf = np.array([[float(Gi[i, j].mid()) for j in range(N + 1)] for i in range(N + 1)])
    Gf = (Gf + Gf.T) / 2
    ew = np.sort(np.linalg.eigvalsh(Gf))[::-1]
    sig = 1 / np.sqrt(np.maximum(ew[:6], 1e-300))
    print(f"N={N} M={M} (M−N={M - N}): σ_min… = " + ", ".join(f"{s:.5g}" for s in sig)
          + f"   (1/σ: {1 / sig[0]:.4g}, {1 / sig[1]:.4g}, {1 / sig[2]:.4g})  [{time.time() - t0:.0f}s]", flush=True)
