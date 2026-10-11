"""Strang C, Runde 4, c27: Gegenprobe jenseits N ≈ 420 — Pol-Antwort längs ∂W, am Schlitz und Eckprofil (= c22/c23 gross).

Matrix A_f = P_N J_∞(f_K) P_N (Strang As Form: J[k][j] = [f^j]_k + [eˣ/f']_{k−j}; c21: gleiche Antworten wie J_N(q_N)),
Pol-Eingabe b_η = η^{−N−1}[q^{N+1}/(η − q)]_{≤N} (Lemma 9.3), y_η = A_f⁻¹ b_η, R_N(η) = |y_η|_{ℓ¹_ρ}.
Punkte: D-Bogen t ∈ [0,1] (t = 0: L, t = 1: −ρ), A-Bogen, Schlitz L(1 + v/N^{0.883}).
Eckprofil (c22): Normanteile von y_L in Bändern j/N und Mittel von y_j/v_j, v = [(L−x)⁻¹]_{≤N}.
Zwischenspeicher: c27_cache.json (je N ein Eintrag). Aufruf: python -I c27_rand_gross.py 500 [prec]
"""
import json
import os
import sys
import time

import numpy as np
from flint import acb, arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = os.path.dirname(os.path.abspath(__file__))
AK = json.load(open(os.path.join(HIER, "..", "a_moden", "daten", "kneser_koeff_lang.json"), encoding="utf-8"))
CACHE = os.path.join(HIER, "c27_cache.json")

N = int(sys.argv[1])
PREC = int(sys.argv[2]) if len(sys.argv) > 2 else int(1.4 * N + 300)
ctx.prec = PREC
t0 = time.time()

L = acb("0.31813150520476413", "1.3372357014306895")
RHO = abs(L)
a = [arb(s) for s in AK[: N + 2]]
a0 = a[0]
R0 = abs(L - a0)
argL = float(L.arg())
argA = float((L - a0).arg())
rho_f, a0_f, R0_f = float(RHO), float(a0), float(R0)

# Punkte (Name, η als acb)
punkte = []
for t in [0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 1.0]:
    phi = argL + t * (np.pi - argL)
    punkte.append((f"D{t}", L if t == 0 else RHO * acb(0, arb(repr(phi))).exp()))
for t in [0.02, 0.1, 1.0]:
    punkte.append((f"A{t}", a0 + R0 * acb(0, arb(repr(argA * (1 - t)))).exp()))
for v in [0.3, 1.0, 4.0]:
    punkte.append((f"S{v}", L * (1 + arb(v) / arb(N) ** arb("0.883"))))

# Spalten [q^j]_{≤N}, j = 0..N+1
q = arb_poly(a[: N + 1])
spalten = []
P = arb_poly([1])
for j in range(N + 2):
    c = P.coeffs()
    spalten.append(c + [arb(0)] * (N + 1 - len(c)))
    if j <= N:
        P = (P * q).truncate(N + 1)
# eˣ/f' bis x^N
fp = [(k + 1) * a[k + 1] for k in range(N + 1)]
ex = [arb(1)]
for k in range(1, N + 1):
    ex.append(ex[-1] / k)
E = []
for k in range(N + 1):
    s = ex[k]
    for i in range(1, k + 1):
        s -= fp[i] * E[k - i]
    E.append(s / fp[0])
J = arb_mat([[spalten[j][k] + (E[k - j] if j <= k else 0) for j in range(N + 1)] for k in range(N + 1)])
t1 = time.time()

# rechte Seiten
qN1 = arb_poly(spalten[N + 1])
spalten_rhs = []
for _, eta in punkte:
    d0 = eta - a0
    r = [1 / d0]
    for k in range(1, N + 1):
        s = acb(0)
        for j in range(1, k + 1):
            s += a[j] * r[k - j]
        r.append(s / d0)
    vr = (qN1 * arb_poly([x.real for x in r])).truncate(N + 1).coeffs()
    vi = (qN1 * arb_poly([x.imag for x in r])).truncate(N + 1).coeffs()
    vr += [arb(0)] * (N + 1 - len(vr))
    vi += [arb(0)] * (N + 1 - len(vi))
    spalten_rhs += [vr, vi]
RHS = arb_mat([[spalten_rhs[c][k] for c in range(len(spalten_rhs))] for k in range(N + 1)])
t2 = time.time()
X = J.solve(RHS, algorithm="approx")
t3 = time.time()

erg = {"prec": PREC, "R": {}, "zeit": [round(t1 - t0), round(t2 - t1), round(t3 - t2)]}
ys = {}
for c, (name, eta) in enumerate(punkte):
    skal = 1 / eta ** (N + 1)
    y = [(acb(X[k, 2 * c], X[k, 2 * c + 1]) * skal) for k in range(N + 1)]
    yr = np.array([complex(float(z.real.mid()), float(z.imag.mid())) for z in y])
    ys[name] = yr
    R = float(np.sum(np.abs(yr) * rho_f ** np.arange(N + 1)))
    erg["R"][name] = R

# Eckprofil (c22) für y_L und y_{−ρ}
Lc = complex(float(L.real), float(L.imag))
jj = np.arange(N + 1)
baender = [(0, .5), (.5, .7), (.7, .8), (.8, .853), (.853, .9), (.9, .95), (.95, 1.0001)]
prof = {}
for name in ("D0.0", "D1.0"):
    yw = ys[name] * rho_f ** jj
    tot = np.sum(np.abs(yw))
    zeilen = []
    for lo, hi in baender:
        sel = (jj / N >= lo) & (jj / N < hi)
        anteil = float(np.sum(np.abs(yw[sel])) / tot)
        verh = ys[name][sel] * Lc ** (jj[sel] + 1)        # y_j / v_j mit v_j = L^{−j−1}
        zeilen.append([lo, hi, anteil, float(np.mean(np.abs(verh))), float(np.real(np.mean(verh))), float(np.imag(np.mean(verh)))])
    prof[name] = zeilen
erg["profil"] = prof
# Lage der Spitze im Ausgang: |y_{−ρ}(ρe^{iφ})|
M = 1 << 14
w = np.fft.ifft(np.concatenate([ys["D1.0"] * rho_f ** jj, np.zeros(M - N - 1)])) * M
betr = np.abs(w)
mh = int(np.argmax(betr[: M // 2 + 1]))
erg["spitze"] = [360.0 * mh / M, float(betr[mh]), float(betr[M // 2])]
wL = np.fft.ifft(np.concatenate([ys["D0.0"] * rho_f ** jj, np.zeros(M - N - 1)])) * M
bL = np.abs(wL)
mL = int(np.argmax(bL[: M // 2 + 1]))
erg["antwort_L"] = [360.0 * mL / M, float(bL[mL])]

print(f"N={N} prec={PREC} Zeit (Aufbau, RHS, Lösen) = {erg['zeit']} s")
print("  R_N: " + "  ".join(f"{k}:{v:.4g}" for k, v in erg["R"].items()))
for name in ("D0.0", "D1.0"):
    print(f"  Profil y_{name}: " + "  ".join(f"[{lo},{hi}): {an:.3f}/{vm:.3g}/({vr:+.3f}{vi:+.3f}i)" for lo, hi, an, vm, vr, vi in prof[name]))
print(f"  Spitze y_(−ρ) auf |x|=ρ: max {erg['spitze'][1]:.4g} bei φ = {erg['spitze'][0]:.2f}°, Wert bei 180°: {erg['spitze'][2]:.3g};"
      f"  y_L: max {erg['antwort_L'][1]:.4g} bei {erg['antwort_L'][0]:.2f}°")
sys.stdout.flush()
d = json.load(open(CACHE, encoding="utf-8")) if os.path.isfile(CACHE) else {}
d[f"{N}_{PREC}"] = erg
json.dump(d, open(CACHE, "w", encoding="utf-8"))
