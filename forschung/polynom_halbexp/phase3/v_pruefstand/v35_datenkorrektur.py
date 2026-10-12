"""Strang V — Datenkorrektur in (H-ii)/H1 für h = k_η (Ks Anfrage 23:59): Π^P_N k_η − Π_N k_η für das echte p_N, η auf dem Schlitz und auf Γ_d.

Π_N k_η = J_N⁻¹T_N[J_∞k_η],  J_∞k_η = 1/(η − f) + m/(η − x),  J_N = T_N[δ∘f + mδ]|_{P_N}      (f = f_K, m = eˣ/f')
Π^P_N k_η = (T_N𝒥_P)⁻¹T_N[𝒥_Pk_η],  𝒥_Pk_η = 1/(η − P) + M_P/(η − x),  M_P = (f∘P − f∘f)/(P − f) = Σ_{n≤5} g_nη_P^{n−1}/n!   (Sekante, H1)
Alles als Keime modulo x^{N+1} (Ks §10.1: nur Taylor-Koeffizienten bei 0 gehen ein), x-Koordinate, arb, 1.5N + 400 Bit.
Ausgabe je η: ‖Π_Nk_η‖_ℓ¹ρ (≥ sup_D̄), ‖Π^P_Nk_η − Π_Nk_η‖_ℓ¹ρ, Verhältnis zu ε_N = ‖p_N − T_N f_K‖_ℓ¹ρ.
Aufruf: python v35_datenkorrektur.py N [pN-Quelle: v33 | gemeinsam]
"""
import json
import math
import sys
import time
from pathlib import Path

from flint import acb, arb, arb_mat, arb_poly, acb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
N = int(sys.argv[1])
QU = sys.argv[2] if len(sys.argv) > 2 else "v33"
bits = int(1.5 * N) + 400
ctx.prec = bits
L = (-acb(-1).lambertw()).conjugate()
rho = abs(L)
AK = json.loads((HIER.parents[1] / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text(encoding="utf-8"))["a"]
if QU == "v33":
    p = [arb(arb(s).mid()) for s in json.loads((HIER / f"v33_pN_{N}.json").read_text(encoding="utf-8"))["p"]]
else:
    p = [arb(arb(s).mid()) for s in json.loads((HIER.parents[1] / "gemeinsam" / "pN_N100_bis_1000.json").read_text(encoding="utf-8"))["p"][str(N)]]
n = N


def tr(q, m=n):
    return q.truncate(m + 1)


def koeff(q, m=n):
    return [q[i] if i <= q.degree() else 0 * q[0] for i in range(m + 1)]


def reihen_div(z, nn, m):
    zc, nc = koeff(z, m), koeff(nn, m)
    q = [zc[0] * 0] * (m + 1)
    for i in range(m + 1):
        s_ = zc[i]
        for j in range(i):
            s_ -= q[j] * nc[i - j]
        q[i] = s_ / nc[0]
    return q


a = [arb(x) for x in AK[:n + 10]]
ne = n + 7
fe = arb_poly(a[:ne + 2])
fpe = fe.derivative()
ex = arb_poly([arb(1) / arb.fac_ui(k) for k in range(ne + 1)])
g = [None, arb_poly(reihen_div(ex, fpe, ne))]
for k in range(1, 6):
    g.append(arb_poly(reihen_div(g[k].derivative(), fpe, ne - k)))
f = arb_poly(a[:n + 1])
m = tr(g[1])
P = arb_poly(p)
eta = tr(P - f)
MP = tr(g[1])
ep = arb_poly([1])
fak = arb(1)
for k in range(1, 6):
    ep = tr(ep * eta)
    fak = fak * k
    MP = MP + tr(tr(g[k + 1]) * ep) * (1 / (fak * (k + 1)))
eps = sum(abs(float(eta[i].mid())) * float(rho.mid()) ** i for i in range(n + 1))


rp = [rho ** i for i in range(n + 1)]


def matrix(Q, mult):
    J = arb_mat(n + 1, n + 1)
    pw = arb_poly([1])
    mc = koeff(mult)
    for j in range(n + 1):
        if j > 0:
            pw = tr(pw * Q)
        c = koeff(pw)
        for i in range(n + 1):
            J[i, j] = (c[i] + (mc[i - j] if i >= j else 0)) * rp[i] / rp[j]     # ρ-skaliert: Zeilen ·ρ^i, Unbekannte ξ_j = x_jρ^j
    return J


t0 = time.time()
Jf = matrix(f, m)
JP = matrix(P, MP)
print(f"N = {N}: ε_N = {eps:.3e}; Matrizen [{time.time() - t0:.0f} s, {bits} Bit]", flush=True)
etas = []
th = L.arg()
for u in ("0.01", "0.05", "0.19"):
    etas.append((f"Schlitz u = {u}", L * (1 + arb(u))))
GL = json.loads((HIER / "v30_Kgamma_S0.0166.json").read_text(encoding="utf-8"))["niveaulinien"]["0.09"]["punkte"]
for idx in (0, len(GL) // 8, len(GL) // 4, 3 * len(GL) // 8, len(GL) // 2):
    xr, xi = GL[idx]
    etas.append((f"Γ_0.09 bei {math.degrees(math.atan2(xi, xr)):+.0f}°", acb(arb(xr), arb(xi))))
# rechte Seiten (komplex ⇒ Real- und Imaginärteil getrennt; Matrizen reell)
k = len(etas)
Bf = arb_mat(n + 1, 2 * k)
BP = arb_mat(n + 1, 2 * k)
for c_, (name, e) in enumerate(etas):
    for (Q, mult, Bm) in ((f, m, Bf), (P, MP, BP)):
        nen = acb_poly([e - Q[0]] + [-Q[i] for i in range(1, n + 1)])          # η − Q(x)
        inv = reihen_div(acb_poly([acb(1)]), nen, n)                           # 1/(η − Q)
        kx = [1 / e ** (i + 1) for i in range(n + 1)]                          # 1/(η − x) = Σ x^i/η^{i+1}
        mc = koeff(mult)
        mk = [sum((mc[j] * kx[i - j] for j in range(i + 1)), acb(0)) for i in range(n + 1)]
        for i in range(n + 1):
            v = (inv[i] + mk[i]) * rho ** i
            Bm[i, c_] = v.real
            Bm[i, k + c_] = v.imag
t1 = time.time()
Xf = Jf.solve(Bf, algorithm="approx")
XP = JP.solve(BP, algorithm="approx")
print(f"   gelöst [{time.time() - t1:.0f} s]", flush=True)
for c_, (name, e) in enumerate(etas):
    nf = sum(math.hypot(float(Xf[i, c_].mid()), float(Xf[i, k + c_].mid())) for i in range(n + 1))
    nd = sum(math.hypot(float((XP[i, c_] - Xf[i, c_]).mid()), float((XP[i, k + c_] - Xf[i, k + c_]).mid())) for i in range(n + 1))
    print(f"   {name:24s}: ‖Π_Nk_η‖_ℓ¹ρ = {nf:.3e}, ‖Π^P_Nk_η − Π_Nk_η‖_ℓ¹ρ = {nd:.3e}  (/ε_N = {nd / eps:.3e})", flush=True)
