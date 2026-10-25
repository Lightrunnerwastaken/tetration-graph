"""E2: Log-Singularität des Siegel-Linearisierers auf dem Randkreis.

Erwartung: v(w) = log(w - w*) + O(1) bei w* = phi(-lambda)/lambda (Urbild des asymptotischen Werts), also
n c_n w*^n -> -1.  Liest cache/c_<name>_<N>.npy (Koeffizienten c_n S^n, S=1.7).
"""
import sys
import numpy as np

name = sys.argv[1] if len(sys.argv) > 1 else "gold"
N = int(sys.argv[2]) if len(sys.argv) > 2 else 8000
S = 1.7
al = {"gold": (3 - 5 ** 0.5) / 2, "sqrt2": 2 ** 0.5 - 1}[name]
lam = np.exp(2j * np.pi * al)
c = np.load(f"cache/c_{name}_{N}.npy")  # c_n S^n
n = np.arange(N + 1)
d = n * c  # n c_n S^n  ~  -(S/w*)^n
# w* aus dem Verhältnis über Blöcke: q = d_n/d_{n+1} = w*/S
for lo, hi in [(500, 1000), (1000, 2000), (2000, 4000), (4000, 7999)]:
    q = d[lo:hi] / d[lo + 1:hi + 1]
    qm = np.median(q.real) + 1j * np.median(q.imag)
    # Kleinste Quadrate: log d_n = log A - n log(w*/S)
    y = np.unwrap(np.angle(d[lo:hi]))
    x = n[lo:hi]
    sl_ph = np.polyfit(x, y, 1)[0]
    sl_ab = np.polyfit(x, np.log(np.abs(d[lo:hi])), 1)[0]
    ws = S * np.exp(-(sl_ab + 1j * sl_ph))
    A = np.mean(d[lo:hi] * (ws / S) ** x)
    print(f"n in [{lo},{hi}): w* (Median-Quot.) = {S*qm:.8f}   w* (Fit) = {ws:.10f}  |w*|={abs(ws):.10f}  A=mean(n c_n w*^n) = {A:.5f}")
print("phi(-lam) = lam*w* =", lam * ws, " |.| =", abs(ws))
np.save(f"cache/wstern_{name}.npy", np.array([ws]))
