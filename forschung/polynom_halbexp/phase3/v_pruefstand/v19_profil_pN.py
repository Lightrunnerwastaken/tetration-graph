import sys, json
sys.path.insert(0, r"D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/phase3/v_pruefstand")
sys.argv = ["v15", "400"]
import v15_A_kopplung as V
from flint import acb, arb, ctx
ctx.prec = 1400
rho = abs(acb(-1).lambertw())
a = V.lade_a()
for N in (318, 400):
    pb = V.BOXEN / f"B_{N}.json"
    if not pb.exists():
        pb = V.BOXEN.parent / "erg" / f"K_{N}.json"
    p = [V.dekodiere(x) for x in json.loads(pb.read_text(encoding="utf-8"))["box"]]
    d = [float((abs(p[k] - a[k]) * rho ** k).mid()) for k in range(N + 1)]
    r = [float((p[k].rad() * rho ** k)) for k in range(N + 1)]
    bins = [(0, 10), (10, 50), (50, int(0.5 * N)), (int(0.5 * N), int(0.9 * N)), (int(0.9 * N), N + 1)]
    print(f"N = {N}: Summe |p_k − a_k|ρ^k je Gradband (Mitten; Boxradius max {max(r):.1e}):")
    for lo, hi in bins:
        print(f"   k ∈ [{lo}, {hi}): {sum(d[lo:hi]):.3e}   max {max(d[lo:hi]):.3e}")
    tail = sum(float((abs(a[k]) * rho ** k).mid()) for k in range(N + 1, len(a)))
    print(f"   Schwanz Σ_(N<k≤640)|a_k|ρ^k = {tail:.3e}")
