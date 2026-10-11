"""O38 — Adjungierte Eckmoden Y_i, T_i (Ws Objekte, W-BERICHT §23.2) und ihre gemeinsame Majorante, ohne Modenwahl.

Bezeichnungen (umskaliert, ℓ²_ρ):  J̃^{-1} = U S Vᵀ,  u_i = U[:, i] (Vorwärts-Eckmode, Einheitsvektor),  y_i = σ_i V[:, i] = J̃^{-T} u_i.
  Y_i(x) = Σ_k y_{i,k} (ρ/x)^k = u_iᵀ J̃^{-1} ẽ_x,                ẽ_x := ((ρ/x)^k)_k  (= abgeschnittener Cauchy-Kern bei x)
  T_i(η) = u_iᵀ (J̃^{-1} c̃_η − δ̃_η),   c = T_N J_∞ (1/(η − ·)),  δ_η = T_N 1/(η − ·)      (W46: Tabs)
Majoranten über ALLE Einheitsvektoren u (also alle Eckmoden zugleich, ohne SVD):
  G(x) := ‖J̃^{-1} ẽ_x‖ ≥ |Y_u(x)|,     G^T(η) := ‖J̃^{-1} c̃_η − δ̃_η‖ ≥ |T_u(η)|.
Punkte: Ecke L; Kreise |t − L| = A·ρ·N^{−1/Re κ₁} (A = 0.25 … 4, wie W46, ψ von der äusseren Normalen); Innenkreise |t − L| = A'·ρ/N.
Aufruf: python o38_adjungiert.py N1,N2,... [modell]   → o38_<modell>.json (wird je N ergänzt)
"""
import json
import os
import sys
import time

import numpy as np
from flint import acb, acb_poly, acb_series, arb, arb_mat, ctx

import o_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(v) for v in sys.argv[1].split(",")]
name = sys.argv[2] if len(sys.argv) > 2 else "kubisch"
HIER = os.path.dirname(os.path.abspath(__file__))
AUS = os.path.join(HIER, f"o38_{name}.json")
L = K.Lc
RHO = K.RHO
argL = np.angle(L)
REK = 1.1182                       # Re κ₁, κ₁ = 2πi/(2πi − log μ)
AS = [0.25, 0.5, 1, 2, 4]
AIN = [1, 4, 16]
PSIS = list(range(-180, 180, 10))
res = json.load(open(AUS, encoding="utf-8")) if os.path.exists(AUS) else {}


def lg(z):
    a = abs(z)
    return float(np.log(a)) if a > 0 else -999.0


for N in NS:
    t0 = time.time()
    n = N + 1
    bits = int(1.3 * N) + 300
    J, M, komplex = K.matrizen(name, N, bits)
    I = arb_mat(n, n)
    for i in range(n):
        I[i, i] = arb(1)
    Ji = J.solve(I, algorithm="approx")
    Jinv = K.to_np(Ji, komplex)
    t1 = time.time()
    U, s, Vt = np.linalg.svd(Jinv)
    r = max(1, int(sum(1 for x in s if x > 1.40)))
    # Modellkoeffizienten (für c_η) in mässiger Genauigkeit
    ctx.prec = 200 + N // 4
    ctx.cap = n + 10                  # Reihenlänge (sonst kappt flint bei 10 Termen!)
    fk, mk, _ = K.modell(name, N)
    fpol = acb_poly([acb(c) for c in fk])
    mpol = acb_poly([acb(c) for c in mk])
    rp = RHO ** np.arange(n)

    def e_vec(x):
        return (RHO / x) ** np.arange(n)

    def c_vec(eta):
        et = acb(eta.real, eta.imag)
        s1 = list(acb_series((acb_poly([et]) - fpol).coeffs(), prec=n).inv().coeffs())
        s1 += [acb(0)] * (n - len(s1))
        kz = acb_poly([et ** (-(j + 1)) for j in range(n)])
        s2 = list((mpol * kz).coeffs())[:n]
        s2 += [acb(0)] * (n - len(s2))
        c = np.array([complex(float(s1[k].real.mid()) + float(s2[k].real.mid()),
                              float(s1[k].imag.mid()) + float(s2[k].imag.mid())) for k in range(n)])
        d = np.array([complex(eta) ** (-(k + 1)) for k in range(n)])
        return c * rp, d * rp

    def auswerten(t):
        ev = e_vec(t)
        g = Jinv @ ev
        cv, dv = c_vec(t)
        gt = Jinv @ cv - dv
        return {"aussen": bool(abs(t) > RHO),
                "G": lg(np.linalg.norm(g)), "GT": lg(np.linalg.norm(gt)),
                "Y": [lg(U[:, i] @ g) for i in range(3)], "T": [lg(U[:, i] @ gt) for i in range(3)],
                "Yr": 0.5 * float(np.log(sum(abs(U[:, i] @ g) ** 2 for i in range(r)))),
                "Tr": 0.5 * float(np.log(sum(abs(U[:, i] @ gt) ** 2 for i in range(r))))}

    out = {"s": [float(x) for x in s[:8]], "r": r, "ecke": auswerten(L), "kreis": {}, "innen": {}}
    for A_ in AS:
        rad = A_ * RHO * N ** (-1 / REK)
        out["kreis"][str(A_)] = [auswerten(L + rad * np.exp(1j * (argL + np.deg2rad(ps)))) for ps in PSIS]
    for A_ in AIN:
        rad = A_ * RHO / N
        out["innen"][str(A_)] = [auswerten(L + rad * np.exp(1j * (argL + np.deg2rad(ps)))) for ps in PSIS]
    out["zeit"] = [t1 - t0, time.time() - t1]
    res[str(N)] = out
    json.dump(res, open(AUS, "w", encoding="utf-8"))
    e = out["ecke"]

    def kb(A_, feld):   # Ws Konvention: innen Y, aussen T
        return max((p[feld[1]] if p["aussen"] else p[feld[0]]) for p in out["kreis"][str(A_)])
    print(f"{name} N={N} ({t1 - t0:.0f}+{time.time() - t1:.0f}s) r={r} lnσ = " + ", ".join(f"{np.log(x):.2f}" for x in s[:4]), flush=True)
    print(f"   Ecke: ln|Y_1| {e['Y'][0]:.3f}  ln|T_1| {e['T'][0]:.3f}  ln G {e['G']:.3f}  ln G^T {e['GT']:.3f}  Yr {e['Yr']:.3f}", flush=True)
    for A_ in AS:
        print(f"   A={A_:<4}: KB(Y_r/T_r) {kb(A_, ('Yr', 'Tr')):.3f}  KB(Y_1/T_1) "
              f"{max((p['T'][0] if p['aussen'] else p['Y'][0]) for p in out['kreis'][str(A_)]):.3f}  "
              f"max ln|Y_1| {max(p['Y'][0] for p in out['kreis'][str(A_)]):.3f}  max ln G {max(p['G'] for p in out['kreis'][str(A_)]):.3f}  "
              f"max ln G^T {max(p['GT'] for p in out['kreis'][str(A_)]):.3f}", flush=True)
    for A_ in AIN:
        pts = out["innen"][str(A_)]
        print(f"   innen A'={A_:<3}: ln|Y_1| {min(p['Y'][0] for p in pts):.3f} … {max(p['Y'][0] for p in pts):.3f}   "
              f"ln G {min(p['G'] for p in pts):.3f} … {max(p['G'] for p in pts):.3f}", flush=True)
