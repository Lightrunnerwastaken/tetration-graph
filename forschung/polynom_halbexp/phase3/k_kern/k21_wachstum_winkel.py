"""k21 — Wachstumsrate je Grad der adjungierten w-Köpfe über dem Winkel: Λ(θ) = log|G_n(θ)/G_{n−1}(θ)| (gemittelt über n-Fenster).
Liegt das Maximum an der Ecke θ_c, so konzentriert sich der adjungierte Kopf asymptotisch dort (Laplace) — dann ist (H1) asymptotisch eine Eckaussage.
Daten: k13_koepfe_2d_kubisch_300.npz."""
import sys
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
d = np.load(sys.argv[1] if len(sys.argv) > 1 else "k13_koepfe_2d_kubisch_300.npz"); R = d["r"]; C = d["c"]; Nr = R.shape[0] - 1
M = 2048
th = 2 * np.pi * np.arange(M) / M - np.pi
G = R @ np.exp(-1j * np.outer(np.arange(Nr + 1), th))
for lo, hi in ([(200, 231), (270, 301)] if Nr <= 300 else [(300, 331), (450, 481), (570, 601)]):
    Lam = np.mean([np.log(np.abs(G[n]) / np.abs(G[n - 1])) for n in range(lo, hi)], axis=0)
    amp = np.abs(G[hi - 1]) ** 2; amp /= amp.sum(); sel = amp > 1e-6
    im = np.argmax(np.where(sel, Lam, -np.inf))
    ws = np.average(Lam, weights=amp)
    print("n ∈ [%d,%d]: Λ(θ) auf dem Träger: max %.5f bei θ = %.1f°, gewichteter Mittelwert %.5f, min %.5f; Λ(θ_c≈±81°) = %.5f"
          % (lo, hi - 1, Lam[im], np.degrees(th[im]), ws, Lam[sel].min(), Lam[np.argmin(np.abs(th - np.radians(81.1)))]))
    # Profil um die Ecke
    for deg in [60, 70, 76, 79, 81, 83, 86, 90, 100, 120, 150, 180]:
        i = np.argmin(np.abs(th - np.radians(deg)))
        print("     θ = %3d°: Λ = %+.5f, Massenanteil-Dichte %.2e" % (deg, Lam[i], amp[i]))
nr = np.linalg.norm(R, axis=1); nc = np.linalg.norm(C, axis=1)
print("Raten der Normen (n = 200…300): ‖r̂_n‖: %.5f, ‖c_n‖: %.5f, Summe (Kopf w): %.5f"
      % (np.polyfit(range(Nr - 100, Nr + 1), np.log(nr[Nr - 100:Nr + 1]), 1)[0], np.polyfit(range(Nr - 100, Nr + 1), np.log(nc[Nr - 100:Nr + 1]), 1)[0],
         np.polyfit(range(Nr - 100, Nr + 1), np.log(nr[Nr - 100:Nr + 1] * nc[Nr - 100:Nr + 1]), 1)[0]))
