"""D28 — Anregung der Eckmode durch Pole AUF dem D-Bogen von ∂W (Strang C, Satz 12.6 braucht: nahe L polynomial,
sonst Rate < 0.068).  Modell: 'kubisch' (f = kubisches Modell, J_N(f) exakt) oder 'exp' (J_N an Π_N f_K, Strang A).
Pole η = |L| e^{iθ}, θ = arg L + τ (τ in Grad), dazu η = L selbst. Ein LU-Lauf: J X = [e_{N−B+1}…e_N | Re/Im P_N y_η],
y_η = J∞k_η = 1/(η − f(x)) + m(x)/(η − x).  Ausgabe: Eckanteil |⟨u1, J⁻¹P_N y_η⟩| und ganze Antwort |J⁻¹P_N y_η − P_N k_η|
(ℓ², |L|-umskaliert), je θ und N.
"""
import json
import sys
import time

import numpy as np
from flint import acb, acb_poly, acb_series, arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
MODELL = sys.argv[1]
NS = [int(x) for x in sys.argv[2].split(",")]
B = 100
L = complex(0.31813150520476413, 1.3372357014306895)
Lr = abs(L)
phiL = np.angle(L)
TAUS = [0.5, 1, 2, 4, 8, 16, 32, 64, 103.38]
datei = K.HIER + f"/d28_bogen_{MODELL}.json"
try:
    aus = json.load(open(datei))
except Exception:
    aus = {}
if MODELL == "kubisch":
    mu = np.sqrt(L)
    rL = (mu - 1) / (L - L.conjugate())
    beta = rL.imag / L.imag
    alpha = rL.real - beta * L.real
    COEF = [Lr ** 2 * alpha, 1 + Lr ** 2 * beta - 2 * L.real * alpha, alpha - 2 * L.real * beta, beta]
else:
    PI = json.load(open(K.GEMEINSAM + "/lin_PiN_fK_bis_1000.json"))
t0 = time.time()
for N in NS:
    n = N + 1
    ctx.prec = int(0.9 * N) + 300
    if MODELL == "kubisch":
        a = [arb(float(c)) for c in COEF] + [arb(0)] * (n - 4)
    else:
        a = [arb(s) for s in PI[str(N)]["b"][:n]]
    P = K.potenzen(a, N)
    _, J = K.FJ(a, N, P)
    x0 = acb_poly([acb(c) for c in a])
    mser = acb_poly([0])
    for j in range(1, n):
        mser += j * acb(a[j]) * acb_poly(P[j - 1].coeffs())
    mser = mser.truncate(n)
    pole = [("L", L)] + [(f"τ={t}", Lr * np.exp(1j * (phiL + np.deg2rad(t)))) for t in TAUS]
    rhs = []
    for name, z in pole:
        zz = acb(z.real, z.imag)
        s1 = acb_series((acb_poly([zz]) - x0).coeffs(), prec=n).inv().coeffs()
        s1 = s1 + [acb(0)] * (n - len(s1))
        zi = 1 / zz
        gl, cur = [], zi
        for _k in range(n):
            gl.append(cur)
            cur = cur * zi
        s2 = (mser * acb_poly(gl)).truncate(n).coeffs()
        s2 = s2 + [acb(0)] * (n - len(s2))
        rhs.append((name, z, [s1[k] + s2[k] for k in range(n)], gl))
    E = arb_mat(n, B + 2 * len(rhs))
    for c in range(B):
        E[N - B + 1 + c, c] = 1
    for i, (name, z, w, gl) in enumerate(rhs):
        for k in range(n):
            E[k, B + 2 * i] = w[k].real
            E[k, B + 2 * i + 1] = w[k].imag
    X = J.solve(E, algorithm="approx")
    k = np.arange(n)
    sc = Lr ** k
    cols = np.arange(N - B + 1, n)
    G = np.array([[float(X[i, c].mid()) for c in range(B)] for i in range(n)]) * sc[:, None] / (Lr ** cols)[None, :]
    U, S, Vt = np.linalg.svd(G, full_matrices=False)
    z = {"s1": float(S[0]), "s2": float(S[1])}
    for i, (name, zc, w, gl) in enumerate(rhs):
        resp = (np.array([float(X[r, B + 2 * i].mid()) for r in range(n)])
                + 1j * np.array([float(X[r, B + 2 * i + 1].mid()) for r in range(n)])) * sc
        kz = np.array([complex(float(gl[r].real.mid()), float(gl[r].imag.mid())) for r in range(n)]) * sc
        z[name] = {"eck1": float(abs(U[:, 0] @ resp)), "eck2": float(abs(U[:, 1] @ resp)),
                   "antwort": float(np.linalg.norm(resp - kz)), "kopf": float(np.linalg.norm(kz))}
    aus[str(N)] = z
    json.dump(aus, open(datei, "w"), indent=1)
    print(f"{MODELL} N={N} ({time.time()-t0:.0f}s) σ1={S[0]:.3f} | " +
          " ".join(f"{nm}: {z[nm]['eck1']:.2e}/{z[nm]['antwort']:.1e}" for nm, *_ in rhs), flush=True)
