"""O25 — Eckmodell der Form (FA): λ_L := inf_σ [a(1+Re z) + b(|z|²+Re z)] / (a + b|z|²),
z = μ^{s−1}, s = −½ + iσ, a = e^{−2σθ₀}, b = e^{−2σ(θ₀+π)}  ⇒  Quotient = 1 + Re z·(1+r)/(1+r|z|²), r = e^{−2πσ}.
Zertifikat mit arb: Bisektion in σ mit Intervallauswertung; Schwänze |σ| ≥ S analytisch (siehe Kommentar).
Ausserdem: sup des Quotienten, die grobe Cauchy–Schwarz-Grösse |z|(1+r)/(1+r|z|²) und das Minimum je Zweig.
μ = √L (Kneser, kubisch), optional anderes μ: python o25_ecke_form.py [mu_re mu_im]
"""
import sys

import mpmath as mp
from flint import arb, acb, ctx

import o_kern as K

sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 200
if len(sys.argv) > 2:
    mu = acb(arb(sys.argv[1]), arb(sys.argv[2]))
else:
    mu = K.Lacb().sqrt()
lnmu = mu.log()
pi = arb.pi()


def quot(sig):
    """sig: arb (Ball). Gibt (q, x, g) zurück: q = 1 + x, x = Re z(1+r)/(1+r|z|²), g = |z|(1+r)/(1+r|z|²)."""
    s = acb(arb(-0.5), sig)
    z = ((s - 1) * lnmu).exp()
    r = (-2 * pi * sig).exp()
    az2 = z.real ** 2 + z.imag ** 2
    den = 1 + r * az2
    x = z.real * (1 + r) / den
    g = az2.sqrt() * (1 + r) / den
    return 1 + x, x, g


# 1) Grobrechnung auf Gitter
mp.mp.dps = 30
best = (10, None)
worst_g = (0, None)
qmax = (-10, None)
for i in range(-60000, 60001):
    sg = i / 1000
    q, x, g = quot(arb(sg))
    qm = float(q.mid())
    if qm < best[0]:
        best = (qm, sg)
    if qm > qmax[0]:
        qmax = (qm, sg)
    gm = float(g.mid())
    if gm > worst_g[0]:
        worst_g = (gm, sg)
print(f"μ = {mu.real.str(12)} + {mu.imag.str(12)} i,  ln μ = {lnmu.real.str(10)} + {lnmu.imag.str(10)} i")
print(f"Gitter |σ| ≤ 60: inf Quotient ≈ {best[0]:.8f} bei σ = {best[1]};  sup ≈ {qmax[0]:.6f} bei σ = {qmax[1]};  "
      f"sup |z|(1+r)/(1+r|z|²) ≈ {worst_g[0]:.6f} bei σ = {worst_g[1]}")

# 2) Zertifikat: auf [−S, S] Bisektion mit Intervallauswertung; Ziel: Quotient ≥ c0
c0 = arb(best[0]) - arb("1e-6")
S = 60


def zert(lo, hi, tiefe=0):
    q, _, _ = quot(arb((lo + hi) / 2, (hi - lo) / 2 * (1 + 2 ** -30)))
    if q > c0:
        return 1
    if tiefe > 40:
        raise RuntimeError(f"keine Entscheidung auf [{lo}, {hi}]: q = {q}")
    m = (lo + hi) / 2
    return zert(lo, m, tiefe + 1) + zert(m, hi, tiefe + 1)


n = 0
for k in range(-S, S):
    n += zert(float(k), float(k + 1))
print(f"Zertifikat: Quotient > {c0.mid().str(10)} auf [−{S}, {S}] ({n} Intervalle)")
# 3) Schwänze. σ ≥ S: |z| = |μ|^{−3/2} e^{−σβ} ≤ |μ|^{−3/2}e^{−Sβ} =: ζ_S, r ≤ e^{−2πS}; |x| ≤ |z|(1+r) ≤ 2ζ_S.
#    σ ≤ −S: r|z|² = |μ|^{−3}e^{2|σ|(π+β)} ≥ 1, |x| ≤ |z|(1+r)/(r|z|²) ≤ (1+1/r)/|z| ≤ 2/|z| ≤ 2|μ|^{3/2}e^{−Sβ}.
beta = lnmu.imag
amu = mu.abs_upper()
t1 = 2 * arb(amu) ** arb(-1.5) * (-S * beta).exp()
t2 = 2 * arb(amu) ** arb(1.5) * (-S * beta).exp()
print(f"Schwänze |σ| ≥ {S}: |x| ≤ {max(t1, t2).upper().str(5)}  ⇒  λ_L ≥ {(c0).str(10)} (zertifiziert, falls Schwanzschranke < Abstand)")
