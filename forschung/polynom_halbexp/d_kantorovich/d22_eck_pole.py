"""D22 — Eckmode und kohärente Eingaben bis N = 1000 (J_N bei x0 = Π_N f_K, Strang A).
Ein LU-Lauf löst J X = [e_{N−B+1} … e_N | Eingaben]: daraus σ1, σ2, u1, u2 (letzte B Spalten, |L|-umskaliert) und für jede
Eingabe w der Eck-Anteil ⟨u_i, J⁻¹w⟩ = σ_i⟨v_i, w⟩ gegen die Grösse |J⁻¹w| der ganzen Antwort.
Eingaben (Pol-Antworten, Strang B/C): w_ζ = P_N[ J∞^{(x0)} k_ζ ] = P_N[ 1/(ζ − x0(x)) + x0'(x0(x))/(ζ − x) ],
ζ = L + u e^{iψ}, ψ = 87.15° (Winkelhalbierende des Aussenkeils bei L), u = v·N^{−0.883} (v = 1, 4) und u = 0.05;
ζ = −1.6, ζ = 2.0. Exakt (Lemma D5): ⟨u_i, J⁻¹w_ζ⟩ = ⟨u_i, P_N k_ζ⟩ + σ_i⟨v_i, B_N (k_ζ)_hi⟩ — der zweite Term ist der
Eck-Anteil des Konsistenzfehlers.
"""
import json
import sys
import time

import numpy as np
from flint import acb, acb_poly, acb_series, arb, arb_mat, arb_poly, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
B = int(sys.argv[2]) if len(sys.argv) > 2 else 100
PF = float(sys.argv[3]) if len(sys.argv) > 3 else 0.85
PI = json.load(open(K.GEMEINSAM + "/lin_PiN_fK_bis_1000.json"))
Lr = 1.3745570107437075
datei = K.HIER + "/d22_eck_pole.json"
try:
    aus = json.load(open(datei))
except Exception:
    aus = {}
PSI = np.deg2rad(87.15)
t0 = time.time()
for N in NS:
    n = N + 1
    ctx.prec = int(PF * N) + 300
    a = [arb(s) for s in PI[str(N)]["b"][:n]]
    _, J = K.FJ(a, N)
    x0 = acb_poly([acb(c) for c in a])
    # Reihe von x0'(x0(x)) bis Grad N
    P = [acb_poly([1])]
    for _ in range(N):
        P.append((P[-1] * x0).truncate(n))
    mser = acb_poly([0])
    for j in range(1, n):
        mser += j * acb(a[j]) * P[j - 1]
    mser = mser.truncate(n)
    L = acb(arb(K.L_RE), arb(K.L_IM))
    pole = {}
    for v in (1.0, 4.0):
        u = v * N ** (-0.883)
        pole[f"L+{v:g}N^-0.883"] = L + acb(u * np.cos(PSI), u * np.sin(PSI))
    pole["L+0.05"] = L + acb(0.05 * np.cos(PSI), 0.05 * np.sin(PSI))
    pole["-1.6"] = acb(-1.6)
    pole["2.0"] = acb(2.0)
    rhs = []
    for name, z in pole.items():
        # 1/(ζ − x0(x)) als Reihe: Inverse von (ζ − x0) mod x^{n}
        s1 = acb_series((acb_poly([z]) - x0).coeffs(), prec=n).inv().coeffs()
        s1 = s1 + [acb(0)] * (n - len(s1))
        # m(x)/(ζ − x) = m(x) Σ x^k ζ^{−k−1}
        zi = 1 / z
        gl, cur = [], zi
        for _k in range(n):
            gl.append(cur)
            cur = cur * zi
        geo = acb_poly(gl)
        s2 = (mser * geo).truncate(n).coeffs()
        s2 = s2 + [acb(0)] * (n - len(s2))
        w = [s1[k] + s2[k] for k in range(n)]
        rhs.append((name, w))
    ncol = B + 2 * len(rhs)
    E = arb_mat(n, ncol)
    for c in range(B):
        E[N - B + 1 + c, c] = 1
    for i, (name, w) in enumerate(rhs):
        for k in range(n):
            E[k, B + 2 * i] = w[k].real
            E[k, B + 2 * i + 1] = w[k].imag
    X = J.solve(E, algorithm="approx")
    k = np.arange(n)
    sc = Lr ** k
    cols = np.arange(N - B + 1, n)
    G = np.array([[float(X[i, c].mid()) for c in range(B)] for i in range(n)]) * sc[:, None] / (Lr ** cols)[None, :]
    U, S, Vt = np.linalg.svd(G, full_matrices=False)
    z = {"s1": float(S[0]), "s2": float(S[1]), "s3": float(S[2])}
    for i, (name, w) in enumerate(rhs):
        resp = (np.array([float(X[r, B + 2 * i].mid()) for r in range(n)])
                + 1j * np.array([float(X[r, B + 2 * i + 1].mid()) for r in range(n)])) * sc
        tot = float(np.linalg.norm(resp))
        e1c = U[:, 0] @ resp
        e2c = U[:, 1] @ resp
        # natürlicher Anteil ⟨u_i, P_N k_ζ⟩ (umskaliert): k_ζ = Σ ζ^{-k-1} x^k
        zc = complex(float(pole[name].real.mid()), float(pole[name].imag.mid()))
        kz = np.array([zc ** (-(kk + 1)) * Lr ** kk for kk in range(n)])
        n1c = U[:, 0] @ kz
        n2c = U[:, 1] @ kz
        z[name] = {"antwort": tot, "eck1": float(abs(e1c)), "eck2": float(abs(e2c)),
                   "nat1": float(abs(n1c)), "nat2": float(abs(n2c)),
                   "gefahr1": float(abs(e1c - n1c)), "gefahr2": float(abs(e2c - n2c))}
    aus[str(N)] = z
    json.dump(aus, open(datei, "w"), indent=1)
    print(f"N={N} ({time.time()-t0:.0f}s) σ1={S[0]:.3f} σ2={S[1]:.3f} σ3={S[2]:.4f} | " +
          " | ".join(f"{nm}: |J⁻¹w|={z[nm]['antwort']:.2e} Ecke {z[nm]['eck1']:.2e} (nat {z[nm]['nat1']:.1e}, Konsist. {z[nm]['gefahr1']:.1e})" for nm, _ in rhs),
          flush=True)
