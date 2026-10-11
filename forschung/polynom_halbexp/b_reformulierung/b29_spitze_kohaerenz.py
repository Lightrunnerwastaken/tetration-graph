"""b29 (Runde 2): Warum sehen kohärente Eingaben die exponentiell wachsende Eckspitze von J_N⁻¹ nicht?
Gewichtet (ℓ², w_k = |L|^k): Ã = D_w J_N⁻¹ D_w⁻¹ = Σ σ_i u_i v_iᵀ.  Spitze = (σ_1, u_1, v_1) (Schwerpunkt bei k ≈ 0.93N).
Überlapp |⟨v_1, w̃⟩|/|w̃| für Eingaben w̃ = D_w·(T_N J_∞ u):  u = f_K, u = k_ζ (Pol nahe L bzw. fern), u = z^N (inkohärent),
Zufallsvektor. Kohärent ⇔ Überlapp ≪ typischer Wert 1/√N."""
import sys
import numpy as np
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b22_polmoden.py", encoding="utf-8").read().split("mp.mp.dps = 100\nres = {}")[0])
mp.mp.dps = 100
aL = abs(L); a0 = mp.mpf(K[0])
psi = ((mp.arg(L - a0) - mp.pi / 2) + (mp.arg(L) + mp.pi / 2)) / 2
Kc = [mp.mpf(c) for c in K]
for N in (200, 300):
    pw, g, Ji = baue(N, J=N)
    w = np.array([float(aL ** k) for k in range(N + 1)])
    Jf = np.array([[float(Ji[i, j].mid().str(30, radius=False)) for j in range(N + 1)] for i in range(N + 1)])
    A = (w[:, None] * Jf) / w[None, :]
    U, sv, Vt = np.linalg.svd(A)
    v1 = Vt[0]; u1 = U[:, 0]
    print(f"N={N}: σ1 = {sv[0]:.4f}, σ2 = {sv[1]:.4f}; Schwerpunkt v1 bei k = {np.dot(np.arange(N+1), v1**2):.1f}")
    # Phase von v1 längs k (Drehwinkel pro Index ~ arg L?)
    kk = np.arange(N - 40, N + 1)
    ph = np.unwrap(np.angle(v1[kk] + 1e-300j))
    # reelle Matrix ⇒ v1 reell; stattdessen Vorzeichenwechsel-Frequenz
    sw = np.sum(np.diff(np.sign(v1[N - 60:])) != 0) / 60
    print(f"   Vorzeichenwechsel pro Index (oberste 60) = {sw:.3f}  (arg L/π = {float(mp.arg(L)/mp.pi):.3f})")
    eingaben = {}
    # f_K: rechte Seite T_N[eˣ + g f_K]
    gk = [mp.mpf(g[k].mid().str(80, radius=False)) for k in range(N + 1)]
    ex = [1 / mp.factorial(k) for k in range(N + 1)]
    rf = [ex[n] + mp.fsum(gk[i] * Kc[n - i] for i in range(n + 1)) for n in range(N + 1)]
    eingaben["f_K"] = rf
    for name, zeta in (("Pol L+0.02e^{iψ}", L + mp.mpf("0.02") * mp.expj(psi)), ("Pol L+0.2e^{iψ}", L + mp.mpf("0.2") * mp.expj(psi)),
                       ("Pol -1.6", mp.mpf("-1.6")), ("Pol 2.0", mp.mpf("2.0"))):
        eingaben[name] = rhs_pol(N, pw, g, zeta)
    eingaben["x^N (Spitze)"] = [0] * N + [1]
    rng = np.random.default_rng(1)
    eingaben["Zufall"] = list(rng.standard_normal(N + 1) / w)
    for name, r in eingaben.items():
        rt = np.array([complex(x) for x in r]) * w
        ov = abs(np.dot(v1, rt)) / np.linalg.norm(rt)
        resp = np.linalg.norm(A @ rt) / np.linalg.norm(rt)
        print(f"   {name:18s}: |⟨v1,w̃⟩|/|w̃| = {ov:.2e}   |Ã w̃|/|w̃| = {resp:.3e}", flush=True)
