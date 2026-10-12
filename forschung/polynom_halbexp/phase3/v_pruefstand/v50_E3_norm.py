"""Strang V — Ks E3: Gs Eckamplitude B_n = √n·‖c_n‖·|Λ_n(L)| bis n = 2000 gegen Q_n·n^{−0.053} (Ks §24.4: B_n = √n·Q_n·|d_n|^{1/2}·|P_n|).

c_n = J̃_n⁻¹e_n (primaler Kopf, Normierung ⟨ℓ_n, J̃c_n⟩ = 1), ‖c_n‖ = ℓ²-Norm = Kreisnorm auf |x| = ρ, δ_n(L) = ε_Lᵀc_n.
Λ_n(L) aus v48 (y_n, gleiche Normierung) — Kontrolle: δ_n(L) hier = w_n aus v48. Q_n := ‖c_n‖/(|δ_n(L)|·|d_n|^{1/2}), d_n = ρn^{−a}.
Aufruf: python v50_E3_norm.py kubisch|kneser n1,n2,… [fac]      (benötigt v48_RL_<MODELL>_N<Nmax>.json mit Nmax ≥ max n)
"""
import json
import math
import sys
import time
from pathlib import Path

HIER = Path(__file__).resolve().parent
MODELL = sys.argv[1]
NS = [int(x) for x in sys.argv[2].split(",")]
FAC = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
sys.argv = ["v48", MODELL, str(max(NS)), "300"]
sys.path.insert(0, str(HIER))
import v48_RL as V48  # noqa: E402
from flint import arb, acb, arb_mat, arb_poly, ctx  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
rl = None
for cand in sorted(HIER.glob(f"v48_RL_{MODELL}_N*.json"), key=lambda p: -int(p.stem.split("_N")[1])):
    d = json.loads(cand.read_text(encoding="utf-8"))
    if d["N"] >= max(NS):
        rl = d["zeilen"]
        break
Lc = complex(float(V48.L_RE), float(V48.L_IM))
rho = abs(Lc)
a_exp = 1 + 1j * Lc / (4 * math.pi)
for n in NS:
    t0 = time.time()
    bits = int(FAC * n) + 400
    ctx.prec = bits
    V48.N = n
    fk, mk = V48.modell(n)
    f = arb_poly(fk)
    rh = (arb(V48.L_RE) ** 2 + arb(V48.L_IM) ** 2).sqrt()
    L = acb(arb(V48.L_RE), arb(V48.L_IM))
    rp = [rh ** k for k in range(n + 1)]
    J = arb_mat(n + 1, n + 1)
    col = arb_poly([1])
    for j in range(n + 1):
        for k in range(n + 1):
            ck = col[k] if k <= col.degree() else arb(0)
            mm = mk[k - j] if k >= j else arb(0)
            J[k, j] = (ck + mm) * rp[k] / rp[j]
        col = (col * f).truncate(n + 1)
    e = arb_mat(n + 1, 1)
    e[n, 0] = arb(1)
    c = J.solve(e, algorithm="approx")
    cf = [c[k, 0] for k in range(n + 1)]
    norm = math.sqrt(sum(float(x.mid()) ** 2 for x in cf))
    q = L / acb(rh)
    s = acb(0)
    for x in reversed(cf):
        s = s * q + x
    dL = complex(float(s.real.mid()), float(s.imag.mid()))
    dn = abs(rho * n ** (-a_exp))
    Q = norm / (abs(dL) * math.sqrt(dn))
    zus = ""
    if rl is not None:
        _, K_, P_, d_, l_ = rl[n]
        lam = complex(*l_)
        B = math.sqrt(n) * norm * abs(lam)
        zus = (f" | δ_n(L) (v48) = {complex(*d_):.6g}; Λ_n(L) = {lam:.6g}; B_n = {B:.5f}; Q_n·n^−0.053·|P_n|-Form: √n Q_n |d_n|^½ |P_n| = "
               f"{math.sqrt(n) * Q * math.sqrt(dn) * abs(complex(*P_)):.5f}")
    print(f"[{MODELL}] n = {n} ({bits} Bit, {time.time() - t0:.0f}s): ‖c_n‖ = {norm:.6g}, δ_n(L) = {dL:.6g}, Q_n = {Q:.5f}" + zus, flush=True)
