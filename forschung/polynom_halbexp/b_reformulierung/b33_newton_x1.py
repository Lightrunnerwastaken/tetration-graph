"""b33 (Runde 3): Newton zentriert bei x₁ = Π_N f_K statt q_N.
Exakt: F(x₁) = T_N[x₁∘x₁ − eˣ] = T_N Q(x₁ − f_K) (rein quadratisch, da T_N J_∞(x₁ − f_K) = 0).
Gemessen: Grössen sup_{|z|=|L|}, Kohärenz (Überlapp der gewichteten Eingabe mit dem Spitzen-Eingangsvektor v₁ von J(f)⁻¹),
Lage der Newton-Korrekturen im Indexraum (Schwerpunkt von |δ_k||L|^k), Abstand zu p_N (aus gemeinsam/ bzw. selbst iteriert)."""
import json, sys, time
import numpy as np
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b11_projektion_moden.py", encoding="utf-8").read().split("tests = {}")[0])
ctx.prec = 900
aL = abs(L); aLf = float(aL)
G2 = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
P1 = json.load(open(G2 + "polynome_bis_194.json")); P2 = json.load(open(G2 + "polynome_195_bis_232.json"))
kreis = [aL * mp.expjpi(2 * mp.mpf(j) / 360) for j in range(360)] + [L]


def FJ(a, N, mitJ=True):
    f = arb_poly(a)
    P = [arb_poly([1])]
    for _ in range(N):
        P.append((P[-1] * f).truncate(N + 1))
    g_ = arb_poly([0])
    for m in range(N + 1):
        g_ += a[m] * P[m]
    F = [g_[k] - 1 / arb(mp.nstr(mp.factorial(k), 400)) for k in range(N + 1)]
    if not mitJ:
        return F, None
    fpf = arb_poly([0])
    for m in range(1, N + 1):
        fpf += m * a[m] * P[m - 1]
    J = arb_mat(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            J[k, j] = P[j][k] + (fpf[k - j] if j <= k else 0)
    return F, J


def supL(v):
    c = [mp.mpf(x.mid().str(40, radius=False)) if hasattr(x, "mid") else x for x in v]
    return max(abs(horner(c, z)) for z in kreis)


def schwerpunkt(v, N):
    w = np.array([abs(float(x.mid().str(20, radius=False))) * aLf ** k for k, x in enumerate(v)])
    return float(np.dot(np.arange(N + 1), w) / w.sum()) / N


for N in (150, 200, 232, 300):
    t0 = time.time()
    pw, g, Ji = baue(N, J=N)                      # J(f_K)⁻¹
    rhs = []
    ex_ = [arb(1)]
    for k in range(1, N + 1):
        ex_.append(ex_[-1] / k)
    gf = (g * arb_poly([arb(c) for c in K[:N + 1]])).truncate(N + 1)
    x1 = Ji * arb_mat([[ex_[k] + gf[k]] for k in range(N + 1)])
    x1 = [arb(x1[k, 0].mid()) for k in range(N + 1)]
    # Spitzenvektor von J(f)⁻¹ (gewichtet)
    Jf = np.array([[float(Ji[i, j].mid().str(25, radius=False)) for j in range(N + 1)] for i in range(N + 1)])
    w = aLf ** np.arange(N + 1)
    A = (w[:, None] * Jf) / w[None, :]
    U, sv, Vt = np.linalg.svd(A)
    F1, Jx1 = FJ(x1, N)
    r1 = np.array([float(x.mid().str(25, radius=False)) for x in F1]) * w
    ov = abs(np.dot(Vt[0], r1)) / np.linalg.norm(r1)
    d_f = Ji * arb_mat([[x] for x in F1])
    d_x = Jx1.solve(arb_mat([[x] for x in F1]))
    dfv = [d_f[k, 0] for k in range(N + 1)]; dxv = [d_x[k, 0] for k in range(N + 1)]
    x2 = [arb((x1[k] - dxv[k]).mid()) for k in range(N + 1)]
    F2, Jx2 = FJ(x2, N)
    d2 = Jx2.solve(arb_mat([[x] for x in F2])); d2v = [d2[k, 0] for k in range(N + 1)]
    r2 = np.array([float(x.mid().str(25, radius=False)) for x in F2]) * w
    ov2 = abs(np.dot(Vt[0], r2)) / np.linalg.norm(r2) if np.linalg.norm(r2) > 0 else 0
    x3 = [arb((x2[k] - d2v[k]).mid()) for k in range(N + 1)]
    zeile = (f"N={N} ({time.time()-t0:.0f}s) σ₁={sv[0]:.3f}: |F(x₁)|-Überlapp mit v₁ = {ov:.2e}; sup|J(f)⁻¹F(x₁)| = {mp.nstr(supL(dfv),3)},"
             f" sup|J(x₁)⁻¹F(x₁)| = {mp.nstr(supL(dxv),3)} (Schwerpunkt k/N {schwerpunkt(dxv,N):.2f}); zweiter Schritt sup = {mp.nstr(supL(d2v),3)}"
             f" (Schwerpunkt {schwerpunkt(d2v,N):.2f}, Überlapp F(x₂) mit v₁ {ov2:.2e})")
    pN = P1.get(str(N)) or P2.get(str(N))
    if pN is None:
        pN = [mp.mpf(x.mid().str(60, radius=False)) for x in x3]
        quelle = "x₃ (selbst iteriert)"
    else:
        pN = [mp.mpf(s) for s in pN]; quelle = "gemeinsam/"
    x1m = [mp.mpf(x.mid().str(60, radius=False)) for x in x1]
    x2m = [mp.mpf(x.mid().str(60, radius=False)) for x in x2]
    zeile += f"\n      p_N aus {quelle}: sup|x₁ − p_N| = {mp.nstr(max(abs(horner(x1m, z) - horner(pN, z)) for z in kreis), 3)}, sup|x₂ − p_N| = {mp.nstr(max(abs(horner(x2m, z) - horner(pN, z)) for z in kreis), 3)}"
    print(zeile, flush=True)
