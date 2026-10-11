"""O09 — Satz O2 (Herglotz-Kriterium) für Strang Es Möbius-exp (Drehung um L, L̄ mit g(0)=1, ρ = |L|):
Δ₁(y,z̄) = ρ²(γy+δ) − (αy+β)z̄ ist hermitesch (ρ²γ = −β); D_H = {Δ₁(y,ȳ) > 0} (Kreisscheibe durch L, L̄, enthält 0);
b(y) = ρ²(γy+δ)(γf(y)+δ)²/det.  K akkretiv in H²(D_ρ) ⟺ Re b ≥ 0 auf D_H.  Hier: min Re b auf ∂D_H (mpmath + arb-Einschluss)."""
import sys
import mpmath as mp
from flint import arb, acb, ctx
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 40
L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
A, B = L.real, L.imag
rho2 = A * A + B * B
c = mp.atan((1 - A) / B) - mp.atan(-A / B)
h = mp.tan(c / 2)
al, be, ga, de = B - A * h, h * rho2, -h, B + h * A
det = al * de - be * ga
print("c =", mp.nstr(c, 12), " Hermitesch: ρ²γ + β =", mp.nstr(rho2 * ga + be, 3))
f = lambda y: (al * y + be) / (ga * y + de)
b = lambda y: rho2 * (ga * y + de) * (ga * f(y) + de) ** 2 / det
# D_H: Δ(y,ȳ) = ρ²δ + ρ²γ(y+ȳ) − α|y|² > 0  →  Kreis: |y − y0|² < R²  (α > 0)
y0 = rho2 * ga / al
R2 = rho2 * de / al + y0 ** 2
R = mp.sqrt(R2)
print("α =", mp.nstr(al, 8), " D_H: Mittelpunkt", mp.nstr(y0, 10), " Radius", mp.nstr(R, 10), " |L − y0| =", mp.nstr(abs(L - y0), 10))
yp = -de / ga
print("Pol von f (und b):", mp.nstr(yp, 10), " Abstand zum Mittelpunkt", mp.nstr(abs(yp - y0), 8), "> R ?", abs(yp - y0) > R)
ws = [y0 + R * mp.expj(p) for p in mp.linspace(0, 2 * mp.pi, 20001)]
vals = [mp.re(b(w)) for w in ws]
i = min(range(len(vals)), key=lambda k: vals[k])
print("min Re b auf ∂D_H =", mp.nstr(vals[i], 10), " bei y =", mp.nstr(ws[i], 8), "; Re b(0) =", mp.nstr(mp.re(b(0)), 8))
# rigoroser Einschluss mit arb über 4096 Teilbögen (Ball-Auswertung)
ctx.prec = 200
y0a, Ra = arb(mp.nstr(y0, 40)), arb(mp.nstr(R, 40))
ala, bea, gaa, dea, rho2a = (arb(mp.nstr(t, 40)) for t in (al, be, ga, de, rho2))
deta = ala * dea - bea * gaa
M = 4096
lo = None
pi = arb.pi()
for k in range(M):
    th = (2 * pi * k + pi * arb(0, 1) + pi) / M      # Ball um Bogenmitte, Radius halber Bogen
    th = 2 * pi * (arb(k) + arb("0.5", "0.5")) / M
    w = acb(y0a) + acb(Ra) * acb(th.cos(), th.sin())
    fw = (acb(ala) * w + acb(bea)) / (acb(gaa) * w + acb(dea))
    bw = acb(rho2a) * (acb(gaa) * w + acb(dea)) * (acb(gaa) * fw + acb(dea)) ** 2 / acb(deta)
    lo = bw.real.lower() if lo is None else min(lo, bw.real.lower())
print("arb: untere Schranke von Re b auf ∂D_H (4096 Bögen):", lo)
