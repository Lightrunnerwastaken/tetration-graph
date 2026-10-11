import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
L = 0.31813150520476413 + 1.3372357014306895j
fig, ax = plt.subplots(figsize=(7, 7))
for N, c in ((60, "tab:blue"), (120, "tab:orange"), (194, "tab:green")):
    r = np.load(f"b05_nullstellen_{N}.npy")
    ax.plot(r.real, r.imag, ".", color=c, label=f"N={N}", ms=4)
t = np.linspace(0, 2 * np.pi, 400)
ax.plot(abs(L) * np.cos(t), abs(L) * np.sin(t), "k-", lw=0.6)
ax.plot([L.real, L.real], [L.imag, -L.imag], "r*")
ax.set_xlim(-3, 3); ax.set_ylim(-3, 3); ax.set_aspect("equal"); ax.legend()
fig.savefig("b05_nullstellen.png", dpi=100, bbox_inches="tight")
