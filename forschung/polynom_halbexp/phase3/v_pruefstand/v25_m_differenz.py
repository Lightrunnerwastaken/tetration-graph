import sys
sys.path.insert(0, r"D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/phase3/v_pruefstand")
N = int(sys.argv[1])
sys.argv = ["v20", str(N), "0.0166"]
import v20_pN_S as V20
from flint import arb, arb_poly, ctx, acb
ctx.prec = 600
rho = abs(acb(-1).lambertw())
p, mp_ = V20.koeff_pN(N)
a, mf = V20.koeff_fK(N)
f = arb_poly(a)
eta = arb_poly([p[j] - a[j] for j in range(N + 1)])
de = eta.derivative()
comp = arb_poly([0])
for c in reversed(de.coeffs()):
    comp = (comp * f).truncate(N + 1) + c
def l1(v): return sum(abs(float(v[j].mid())) * float(rho.mid()) ** j for j in range(N + 1))
dm_exakt = [mp_[j] - mf[j] for j in range(N + 1)]
print(f"N = {N}: ‖η‖ = {l1([p[j]-a[j] for j in range(N+1)]):.3e}; ‖T_N(η'∘f)‖ = {l1([comp[j] if j <= comp.degree() else arb(0) for j in range(N+1)]):.3e}; "
      f"‖m_pN − m_fK‖ (exakt) = {l1(dm_exakt):.3e}")
