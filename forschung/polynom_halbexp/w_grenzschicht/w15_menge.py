"""W15 — Welche Menge K macht die (gespiegelte) adjungierte Eckmode Bernstein–Walsh-extremal?
Daten: Φ_y(t) = (1/N)ln|y(t)| (arb, N=600) innen; Vergleich mit g_K(1/t̄) für Zwei-Kreis-Mengen K = D̄ ∪ D̄(c, R) mit Kreis durch L̃, L̃̄
und rechtem Scheitel b (Parameter). Rescaled (|L̃| = 1).
"""
import sys, time
import numpy as np
from flint import acb, arb, arb_mat, ctx
import w_kern as W
sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1]); n = N + 1; ctx.prec = int(1.5 * N) + 400
coef, _ = W.kubik_koeff(); rho = W.RHO
d = np.load(W.HIER + f"/w01_N{N}.npz"); u1 = d["U"][:, 0]
J = W.J_mat(W.kubik_a(N, coef), N); rr = arb(rho); Dk = [rr ** k for k in range(n)]
Z = J.transpose().solve(arb_mat(n, 1, [Dk[k] * arb(float(u1[k])) for k in range(n)]), algorithm="approx")
y = [Z[k, 0] / Dk[k] for k in range(n)]
def pot(t):
    w = 1 / acb(t.real, t.imag); s = acb(0)
    for v in y[::-1]: s = s * w + v
    m = abs(s).mid(); return float(m.log()) / N
Lt = W.L / rho
def green2(w, b):
    """Green-Fkt. (Pol ∞) von D̄ ∪ D̄(c,R), zweiter Kreis durch L̃, L̃̄ und b (reell > 1)."""
    m = (w - Lt) / (w - np.conj(Lt))
    th1 = np.angle((-1 - Lt) / (-1 - np.conj(Lt))); th2 = np.angle((b - Lt) / (b - np.conj(Lt)))
    def inside(lo, hi, t): return ((t - lo) % (2*np.pi)) < ((hi - lo) % (2*np.pi))
    if inside(th2, th1, 0.0): lo, wid = th2, (th1 - th2) % (2*np.pi)
    else: lo, wid = th1, (th2 - th1) % (2*np.pi)
    ang = (np.angle(m) - lo) % (2*np.pi)
    if ang > wid: return 0.0
    tau = abs(m) ** (np.pi/wid) * np.exp(1j*ang*np.pi/wid); tau0 = np.exp(1j*((0.0 - lo) % (2*np.pi))*np.pi/wid)
    return float(np.log(abs((tau - np.conj(tau0)) / (tau - tau0)))), wid/np.pi
pts = [(r, p) for r in (0.3, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95) for p in (0, 30, 60, 76.6, 90, 120, 180)]
data = {(r, p): pot(r*np.exp(1j*np.deg2rad(p))) for r, p in pts}
floor = pot(3.0 + 0j)
print(f"N={N}, Sockel {floor:.4f}")
for b in (1.346, 1.30, 1.25, 1.227, 1.20, 1.15):
    errs = []
    for (r, p), v in data.items():
        g = green2(1/(r*np.exp(-1j*np.deg2rad(p))), b)
        gv = g[0] if isinstance(g, tuple) else g
        if gv > 0: errs.append(v - gv)
    al = green2(3.0+0j, b)[1]
    print(f"b={b:.3f} (α={al:.4f}): Abweichung Φ_y − g_K: mittel {np.mean(errs):+.4f}, min {np.min(errs):+.4f}, max {np.max(errs):+.4f}")
b = 1.227
print("Einzelwerte für b = 1.227: (r, φ): Φ_y / g_K")
for (r, p), v in data.items():
    g = green2(1/(r*np.exp(-1j*np.deg2rad(p))), b); gv = g[0] if isinstance(g, tuple) else g
    print(f"  ({r:.2f},{p:5.1f}): {v:.4f} / {gv:.4f}")
