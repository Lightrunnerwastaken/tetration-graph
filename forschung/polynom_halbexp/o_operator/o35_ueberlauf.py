"""O35 — Eckmoden des quadratischen Abschnitts J_N = T_N J_∞|P_N und ihr Überlauf (Modell exakt bzw. Kneser-Daten).
Rechnung: Y = J_N⁻¹ in arb (gut skaliert, ‖Y‖ ~ σ₁), Gleitkomma-SVD von Y ⇒ Eckmoden v_i = linke Singulärvektoren von Y
(J_N v_i = u_i/s_i). Dann volles Bild J_∞v_i (Grade 0…3N): Massenverteilung oberhalb N.
Rechteck-Abschnitte: σ_min(R_d) = 1/√λ_max((R_dᵀR_d)⁻¹), (R_dᵀR_d)⁻¹ in arb (gut skaliert).
Aufruf: python o35_ueberlauf.py N [D] [modell]
"""
import sys

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

import o_kern as K

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
D = int(sys.argv[2]) if len(sys.argv) > 2 else 8
name = sys.argv[3] if len(sys.argv) > 3 else "kubisch"
ctx.prec = int(3.0 * N) + 400
NR = 3 * N
fk, mk, _ = K.modell(name, NR)
f = arb_poly(fk)
rho = (arb(K.L_RE) ** 2 + arb(K.L_IM) ** 2).sqrt()
pw = [arb_poly([1])]
for _ in range(N):
    pw.append((pw[-1] * f).truncate(NR + 1))
rp = [rho ** k for k in range(NR + 1)]
A = arb_mat(NR + 1, N + 1)
for j in range(N + 1):
    col = pw[j]
    for k in range(NR + 1):
        v = col[k] if k <= col.degree() else arb(0)
        if k >= j and k - j < len(mk):
            v = v + mk[k - j]
        A[k, j] = v * rp[k] / rp[j]


def sub(M, r0, r1):
    S = arb_mat(r1 - r0, M.ncols())
    for i in range(r0, r1):
        for j in range(M.ncols()):
            S[i - r0, j] = M[i, j]
    return S


def tonp(M):
    return np.array([[float(M[i, j].mid()) for j in range(M.ncols())] for i in range(M.nrows())])


# quadratisch
JN = sub(A, 0, N + 1)
I = arb_mat(N + 1, N + 1)
for i in range(N + 1):
    I[i, i] = arb(1)
Y = tonp(JN.solve(I, algorithm="approx"))
U, s, Vt = np.linalg.svd(Y)
print(f"{name} N={N}: σ(J_N⁻¹) grösste: " + ", ".join(f"{x:.4f}" for x in s[:6]) + f"  ⇒ σ_min(J_N) = {1 / s[0]:.5f}", flush=True)
# Rechteck
for d in sorted(set([1, 2, 3, 4, 6, 8, D])):
    if d > D:
        continue
    R = sub(A, 0, N + 1 + d)
    G = R.transpose() * R
    Gi = tonp(G.solve(I, algorithm="approx"))
    lam = np.linalg.eigvalsh((Gi + Gi.T) / 2)
    print(f"d = {d:2d}: σ_min(T_(N+d)J|P_N) = " + ", ".join(f"{1 / np.sqrt(x):.5f}" for x in lam[::-1][:4]), flush=True)
# Überlauf der Eckmoden (Eingänge v_i = U[:, i], denn Y = U S Vᵀ ⇒ J_N U_i = V_i/s_i)
Af = None
for i in range(3):
    # v_i in voller Genauigkeit: J_N v = s_i·Vt[i] (exakt gelöst), dann volles Bild A·v in arb
    yv = arb_mat(N + 1, 1)
    for j in range(N + 1):
        yv[j, 0] = arb(float(s[i] * Vt[i, j]))
    va = JN.solve(yv, algorithm="approx")
    v = np.array([float(va[j, 0].mid()) for j in range(N + 1)])
    outa = A * va
    out = np.array([abs(float(outa[k, 0].mid())) for k in range(NR + 1)])
    tot = np.sum(out ** 2)
    lo = np.sum(out[:N + 1] ** 2)
    prof = [np.sum(out[N + 1:N + 1 + dd] ** 2) / tot for dd in (1, 2, 3, 5, 8, 12, 20, 50)]
    kschw = np.sum(v ** 2 * np.arange(N + 1)) / np.sum(v ** 2)
    hi = np.sum(out[int(1.2 * N):] ** 2) / tot
    print(f"Mode {i} (s = {s[i]:.3f}): Eingang k/N-Schwerpunkt {kschw / N:.3f}; ‖J∞v‖ = {np.sqrt(tot):.4f}, Anteil ≤ N {lo / tot:.2e}; "
          f"Überlauf N+1…N+d (d = 1/2/3/5/8/12/20/50): " + " ".join(f"{p:.3f}" for p in prof) + f"; > 1.2N: {hi:.2e}", flush=True)
    top = np.argsort(out[N + 1:])[::-1][:5] + N + 1
    print(f"      grösste Überlauf-Grade: {list(top)}  Werte {[round(out[t], 4) for t in top]}", flush=True)
    if i == 0:
        # Profil auf dem Kreis: |v|, |T_N J v|, |J v| (letzteres gross)
        ph = np.linspace(-np.pi, np.pi, 721)
        z = np.exp(1j * ph)
        cv = np.array([float(va[j, 0].mid()) for j in range(N + 1)])
        co = np.array([float(outa[k, 0].mid()) for k in range(NR + 1)])
        Vz = np.polynomial.polynomial.polyval(z, cv)
        Tz = np.polynomial.polynomial.polyval(z, co[:N + 1])
        Tz2 = np.polynomial.polynomial.polyval(z, co[:N + 1 + 4])
        thL = np.degrees(np.angle(K.Lc))
        th = np.angle(K.Lc)
        for yv_ in (-200, -100, -50, -20, -10, -5, -2, 0, 2, 5, 10, 20, 40):
            phv = th - yv_ / N
            zz = np.exp(1j * phv)
            print(f"      y = N(θ_L−φ) = {yv_:+5d}: |v| = {abs(np.polynomial.polynomial.polyval(zz, cv)):.3e}, |T_N Jv| = {abs(np.polynomial.polynomial.polyval(zz, co[:N + 1])):.3e}", flush=True)
        prof = np.log10(np.abs(co) + 1e-300)
        print("      log10|[Jv]_k| bei k = " + " ".join(f"{k}:{prof[k]:.1f}" for k in list(range(0, NR + 1, max(1, N // 10)))), flush=True)
        print("      log10|v_k|  bei k/N = " + " ".join(f"{k / N:.2f}:{np.log10(abs(cv[k]) + 1e-300):.1f}" for k in range(0, N + 1, max(1, N // 20))), flush=True)
        print("      log10|[Jv]_k| bei k/N (Kante) = " + " ".join(f"{k / N:.3f}:{prof[k]:.1f}" for k in range(int(0.9 * N), int(1.3 * N), max(1, N // 50))), flush=True)
        print("      fein: k−N = " + " ".join(f"{k - N}:{prof[k]:.2f}" for k in range(N - 4, N + 13)), flush=True)
