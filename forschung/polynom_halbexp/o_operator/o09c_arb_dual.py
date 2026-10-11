"""O09c — rigorose arb-Schranke für Re b auf ∂D_H im Dualbild (Möbius-exp, K⁻¹ = M_{f'}C_h auf H²(D(a₀,R₀))).
Alle Grössen in arb aus den exakten Definitionen (L als 60-stellige Zeichenkette; Rundung im Ball enthalten)."""
import sys
from flint import arb, acb, ctx
sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 200
A = arb("0.31813150520476413531265425158766451720351761387139986692237")
B = arb("1.33723570143068940890116214319371061253950213846051241687")
L = acb(A, B)
rho2 = A * A + B * B
c = ((1 - A) / B).atan() - (-A / B).atan()
h = (c / 2).tan()
al, be, ga, de = B - A * h, h * rho2, -h, B + h * A
det = al * de - be * ga
def f(y): return (acb(al) * y + acb(be)) / (acb(ga) * y + acb(de))
def fp(y): return acb(det) / (acb(ga) * y + acb(de)) ** 2
def kriterium(pa, pb, pc, pd, afun, r2, name):
    # Δ₁ = r²(pc Y + pd) − (pa Y + pb) Z̄ ; hier alle reell ⇒ κ = 0 wenn pd > 0
    assert (r2 * pd) > 0
    # hermitesch: r² pc = −pb ?
    herm = r2 * pc + pb
    A_ = -pa                      # Koeffizient von |Y|²
    Yc = -(r2 * pc) / A_          # c0 + 2 r² pc Re Y + A_|Y|² > 0, A_ < 0
    R2 = Yc * Yc - (r2 * pd) / A_
    R = R2.sqrt()
    M = 4096
    lo = None
    pi = arb.pi()
    for k in range(M):
        th = 2 * pi * (arb(k) + arb("0.5", "0.5")) / M
        w = acb(Yc) + acb(R) * acb(th.cos(), th.sin())
        bw = acb(r2) * (acb(pc) * w + acb(pd)) * afun(w)
        lo = bw.real.lower() if lo is None else min(lo, bw.real.lower())
    print(f"{name}: hermitesch-Defekt {herm}; D_H Mitte {Yc.mid()} Radius {R.mid()}; untere Schranke Re b auf ∂D_H: {float(lo):.6f}")
kriterium(al, be, ga, de, lambda Y: 1 / fp(f(Y)), rho2, "Bild 0")
a0 = be / de
R0sq = (A - a0) ** 2 + B * B
pa = de + a0 * ga; pb = de * a0 - be - a0 * (-ga * a0 + al); pc = -ga; pd = -ga * a0 + al
kriterium(pa, pb, pc, pd, lambda Y: fp(Y + acb(a0)), R0sq, "Dualbild a₀")
