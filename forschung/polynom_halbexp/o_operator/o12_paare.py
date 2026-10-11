"""O12 — Struktur der Eckmoden-Paare (kubisch): Ausgangsvektor u_i (Lage auf dem Kreis, lokaler Exponent radial bei L),
Eingangsvektor v_i (Lage in k: Masse auf k > N − d, Schwerpunkt), aus o07b_kubisch_N.npz."""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
d = np.load(f"o07b_kubisch_{N}.npz")
s, U, V = d["s"], d["U"], d["V"]
Lc = complex(0.31813150520476413, 1.3372357014306895); R = abs(Lc); thL = np.angle(Lc)
k = np.arange(N + 1)
phi = np.linspace(-np.pi, np.pi, 8001)
for i in range(min(6, U.shape[1])):
    u, v = U[:, i], V[:, i]
    # Polynom im umskalierten Argument x/ρ
    vals = np.polynomial.polynomial.polyval(np.exp(1j * phi), u)
    th_max = np.degrees(phi[np.argmax(np.abs(vals))])
    # radial bei L: |u(L(1−r))|, r = 0.003…0.3
    rs = np.array([0.3, 0.1, 0.03, 0.01, 0.003])
    ur = np.abs(np.polynomial.polynomial.polyval((1 - rs) * np.exp(1j * thL), u))
    expo = np.diff(np.log(ur)) / np.diff(np.log(rs))
    w = np.abs(v) ** 2
    vals_v = np.polynomial.polynomial.polyval(np.exp(-1j * phi), v)
    gin = (np.abs(phi) > thL)
    print(f"i={i+1} σ={s[i]:.4g}: u max bei {th_max:.2f}°, radial lok. Exponent (r=0.3→0.003): {np.array2string(expo, precision=2)}; "
          f"v: Schwerpunkt N−k = {N - (w*k).sum():.1f}, Masse k>N−60: {w[k > N-60].sum():.2f}, |V| auf Γin/Γout max: "
          f"{np.abs(vals_v[gin]).max():.2e}/{np.abs(vals_v[~gin]).max():.2e}")
