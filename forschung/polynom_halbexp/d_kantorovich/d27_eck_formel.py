"""D27 — Lemma D6 im kubischen Modell: Eckanteil einer Pol-Antwort aus der adjungierten Eckmode.

f(x) = x + (x² − 2Re(L)x + |L|²)(α + βx) (Fixpunkte L, L̄, f'(L) = √L), m = f'∘f (Grad 6), J_N = J_N(f) exakt.
SVD der |L|-umskalierten Inversen: J⁻¹ = Σ σ_i u_i v_iᵀ (reell). V_i(x) := Σ_k v_{i,k} x^{−k} mit v_{i,k} = ṽ_{i,k}|L|^k.
Lemma D6:  ⟨u_i, J⁻¹ P_N y_ζ⟩ = σ_i [ m(ζ)V_i(ζ)/ζ + Σ_j V_i(x_j)/(x_j f'(x_j)) + ℓ_i(P_N Q_ζ) ],
y_ζ = J∞k_ζ = 1/(ζ − f(x)) + m(x)/(ζ − x),  f(x_j) = ζ,  Q_ζ(x) = (m(x) − m(ζ))/(ζ − x).
Teil 1: Prüfung der Formel gegen die direkt gelöste Antwort. Teil 2: Terme der Formel längs Kurven (D-Bogen, Ecke L).
"""
import json
import sys
import time

import numpy as np
from flint import acb, arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
L = complex(0.31813150520476413, 1.3372357014306895)
Lr = abs(L)
mu = np.sqrt(L)
rL = (mu - 1) / (L - L.conjugate())
beta = rL.imag / L.imag
alpha = rL.real - beta * L.real
coef = [Lr ** 2 * alpha, 1 + Lr ** 2 * beta - 2 * L.real * alpha, alpha - 2 * L.real * beta, beta]
f = np.polynomial.Polynomial(coef)
fp = f.deriv()
m = fp(f)                                   # Polynom Grad 6
t0 = time.time()
n = N + 1
ctx.prec = int(0.9 * N) + 300
a = [arb(float(c)) for c in coef] + [arb(0)] * (n - 4)
_, J = K.FJ(a, N)
Ji = J.solve(arb_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)]), algorithm="approx")
k = np.arange(n)
sc = Lr ** k
Jif = np.array([[float(Ji[i, j].mid()) for j in range(n)] for i in range(n)])
G = Jif * sc[:, None] / sc[None, :]
U, S, Vt = np.linalg.svd(G)
print(f"N={N}: σ = {S[0]:.4f}, {S[1]:.4f}, {S[2]:.4f}, {S[3]:.4f}  ({time.time()-t0:.0f}s)", flush=True)


def V(i, x):
    """V_i(x) = Σ ṽ_k (|L|/x)^k (Horner in w = |L|/x)."""
    w = Lr / x
    s = 0j
    for c in Vt[i][::-1]:
        s = s * w + c
    return s


def ell_poly(i, poly):
    """ℓ_i(P_N poly) = Σ_k ṽ_k |L|^k poly_k."""
    c = poly.coef
    return sum(Vt[i][kk] * Lr ** kk * c[kk] for kk in range(min(len(c), n)))


def terme(i, z):
    roots = (f - z).roots()
    T0 = m(z) * V(i, z) / z
    Tj = [V(i, x) / (x * fp(x)) for x in roots]
    # Q_ζ(x) = (m(x) − m(ζ))/(ζ − x) als Polynom in x
    num = m - m(z)
    Q, rem = divmod(num, np.polynomial.Polynomial([z, -1.0]))
    TQ = ell_poly(i, Q)
    return T0, Tj, TQ, roots


def antwort(z):
    """J⁻¹ P_N y_ζ direkt (umskaliert): y_ζ-Koeffizienten aus Reihen."""
    # 1/(ζ − f(x)) als Potenzreihe bis Grad N, m(x)/(ζ − x) ebenso
    g = np.zeros(n, dtype=complex)
    # Reihe von 1/(ζ − f): rekursiv, f Polynom Grad 3
    c0 = z - coef[0]
    fc = [-coef[1], -coef[2], -coef[3]]
    g[0] = 1 / c0
    for kk in range(1, n):
        s = 0j
        for d in range(1, 4):
            if kk - d >= 0:
                s += fc[d - 1] * g[kk - d]
        g[kk] = -s / c0
    geo = z ** (-(k + 1.0))
    mc = m.coef
    h = np.zeros(n, dtype=complex)
    for kk in range(n):
        h[kk] = sum(mc[d] * geo[kk - d] for d in range(len(mc)) if kk - d >= 0)
    y = (g + h) * sc
    return G @ y


if __name__ == "__main__":
    out = {"N": N, "s": [float(x) for x in S[:4]]}
    print("Teil 1 — Prüfung Lemma D6 (Eck-Komponente direkt vs. Formel):")
    for z in (L + 0.02 * np.exp(1j * np.deg2rad(87.15)), -1.45 + 0j, 1.2 * np.exp(2.5j), L + 0.005 * np.exp(1j * np.deg2rad(87.15))):
        r = antwort(z)
        for i in (0, 1):
            direkt = U[:, i] @ r
            T0, Tj, TQ, roots = terme(i, z)
            formel = S[i] * (T0 + sum(Tj) + TQ)
            print(f"  ζ={z:.4f} i={i+1}: direkt {direkt:.6e}  Formel {formel:.6e}  |rel.Diff| {abs(direkt-formel)/abs(direkt):.1e}")
    print("Teil 2 — Terme längs Kurven (i = 1): |σ T0|, |σ T_j| (j = Urbild nahe ζ bzw. L), Summe = Eckanteil")
    zeilen = []
    for name, zs in (("D-Bogen |ζ| = |L|(1+2/N)", [Lr * (1 + 2 / N) * np.exp(1j * np.deg2rad(th)) for th in (80, 90, 120, 150, 180)]),
                     ("Winkelhalb. L + u e^{iψ}", [L + u * np.exp(1j * np.deg2rad(87.15)) for u in (0.3 / N, 1 / N, 3 / N, N ** -0.883, 0.02, 0.05)])):
        print("  " + name)
        for z in zs:
            T0, Tj, TQ, roots = terme(0, z)
            jmin = int(np.argmin(np.abs(roots - L)))
            tot = S[0] * (T0 + sum(Tj) + TQ)
            print(f"    ζ={z:.4f}: |σT0|={abs(S[0]*T0):.2e}  |σT_near|={abs(S[0]*Tj[jmin]):.2e} (Urbild {roots[jmin]:.3f}, |x|/|L|={abs(roots[jmin])/Lr:.4f})"
                  f"  andere {[f'{abs(S[0]*t):.1e}' for jj, t in enumerate(Tj) if jj != jmin]}  |σTQ|={abs(S[0]*TQ):.1e}  → |Eckanteil|={abs(tot):.3e}")
            zeilen.append({"kurve": name, "z": [z.real, z.imag], "T0": abs(S[0] * T0), "Tnear": abs(S[0] * Tj[jmin]), "eck": abs(tot)})
    out["zeilen"] = zeilen
    json.dump(out, open(K.HIER + f"/d27_eck_formel_{N}.json", "w"), indent=1)
    np.savez_compressed(K.HIER + f"/d27_vek_{N}.npz", U=U[:, :4], Vt=Vt[:4], S=S[:6], G=G)
