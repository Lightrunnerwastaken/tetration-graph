"""O23 — Hochfrequenz-Infimum c_h(n) := inf_{δ ∈ span{x^{n+1},…,x^N}} ‖T_{3N}J_∞δ‖/‖δ‖ (ℓ²_ρ) für kubisch und Kneser (Mittelwert-Rechnung),
dazu die Lage des Minimalvektors (Ort auf dem Kreis, Frequenzschwerpunkt). Frage: liegt c_h(n) ≥ min_{Γin}|m| für alle n (Reserve für grobe Abschätzungen)?"""
import sys
import numpy as np
from flint import arb, arb_mat, arb_poly, ctx
import o_kern as K
sys.stdout.reconfigure(encoding="utf-8")
name = sys.argv[1]; N = int(sys.argv[2])
ctx.prec = int(3.0 * N) + 300
ctx.cap = 4 * N + 20
fk, mk, _ = K.modell(name, 3 * N)
f = arb_poly(fk)
rho = (arb(K.L_RE) ** 2 + arb(K.L_IM) ** 2).sqrt()
pw = [arb_poly([1])]
for _ in range(N):
    pw.append((pw[-1] * f).truncate(3 * N + 1))
R = np.zeros((3 * N + 1, N + 1))
rp = [rho ** k for k in range(3 * N + 1)]
# exakte Spalten in arb, dann QR in mpmath-Präzision? Hier: Gram in arb, Eigenwerte in mpmath
import mpmath as mp
mp.mp.dps = int(1.0 * N) + 60
G = arb_mat(N + 1, N + 1)
cols = []
for j in range(N + 1):
    col = []
    for k in range(3 * N + 1):
        v = pw[j][k] if k <= pw[j].degree() else arb(0)
        if k >= j:
            v = v + mk[k - j]
        col.append(v * rp[k] / rp[j])
    cols.append(col)
Gm = mp.matrix(N + 1, N + 1)
for i in range(N + 1):
    for j in range(i, N + 1):
        s = arb(0)
        for k in range(3 * N + 1):
            s += cols[i][k] * cols[j][k]
        Gm[i, j] = Gm[j, i] = mp.mpf(s.mid().str(mp.mp.dps, radius=False))
for n in [-1, 0, 5, 10, 20, 40, 60, 80, 100, 150]:
    if n + 1 > N - 5:
        break
    sub = Gm[n + 1:, n + 1:]
    ev, Q = mp.eigsy(sub)
    i = min(range(len(ev)), key=lambda t: ev[t])
    v = np.array([float(Q[r, i]) for r in range(len(ev))])
    ks = np.arange(n + 1, N + 1)
    phi = np.linspace(-np.pi, np.pi, 4001)
    vals = np.abs(np.exp(1j * np.outer(phi, ks)) @ v)
    print(f"{name} N={N}: c_h(n={n}) = {mp.nstr(mp.sqrt(ev[i]), 8)}; Minimalvektor: k-Schwerpunkt {(v**2 * ks).sum()/(v**2).sum():.1f}, Max auf Kreis bei {np.degrees(phi[np.argmax(vals)]):.1f}°", flush=True)
