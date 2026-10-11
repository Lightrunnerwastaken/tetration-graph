"""exp.tpid6 — TPID 6 mit Rate: x_n mit ^n x_n = n geht gegen η, und zwar wie η + 2π²η/(e n²).

 TPID 6 (Tetration-Forum): lim x_n = η. Rate (hier): Für x = η + δ braucht der Turm am Engpass bei e etwa
 N(δ) ≈ π √(2η/(e δ)) Schritte (bw.lem.randgesetz, br.satz.sattel_knoten); danach wächst er in wenigen Schritten über n.
 Also n ≈ π √(2η/(e δ_n)) + O(log* n)  ⇒  n² (x_n − η) → 2π²η/e = 10.4906…
 Rechnung in doppelter Genauigkeit: δ_n ≥ 1e−10, Rundung 2e−16 je Schritt ≪ e²δ/η (Engpass-Vorschub).
"""
import math

import mpmath as mp

from gemeinsam import Lauf, plt

lauf = Lauf("exp.tpid6", __file__)
ETA = math.exp(1 / math.e)
KONST = 2 * math.pi ** 2 * ETA / math.e


def erreicht(x, n):
    """True, wenn ^n x ≥ n (Turm von 1 aus, n Stockwerke)."""
    lx, t, ln_n = math.log(x), 1.0, math.log(n)
    for _ in range(n):
        if lx * t >= ln_n:                     # nächstes Stockwerk ist schon ≥ n (vermeidet Überlauf)
            return True
        t = math.exp(lx * t)
    return t >= n


def x_n(n):
    lo, hi = ETA, 2.0
    if n <= 4:
        lo = 1.2
    for _ in range(80):
        mid = (lo + hi) / 2
        if erreicht(mid, n):
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-16 * hi:
            break
    return (lo + hi) / 2


ns = [3, 5, 10, 30, 100, 300, 1000, 3000, 10_000, 30_000, 100_000, 300_000]
werte = {n: x_n(n) for n in ns}
c = {n: n * n * (werte[n] - ETA) for n in ns}
for n in ns:
    print(f"  n = {n:>7}  x_n = {werte[n]:.15f}  n²(x_n − η) = {c[n]:.6f}")

lauf.pruefe("TPID 6: x_n > η für n ≥ 3 und x_n fällt ab n = 3", all(werte[n] > ETA for n in ns) and all(werte[ns[i]] > werte[ns[i + 1]] for i in range(len(ns) - 1)))
lauf.pruefe("TPID 6: x_n → η (Abstand bei n = 3·10⁵ unter 2·10⁻¹⁰)", werte[300_000] - ETA < 2e-10, f"{werte[300_000] - ETA:.3e}")
# Anpassung c_n = K (1 + β/n) über die letzten Punkte (n ≥ 1000)
gross = [n for n in ns if n >= 1000]
A = [[1.0, 1.0 / n] for n in gross]
yv = [c[n] for n in gross]
s11 = sum(a[0] * a[0] for a in A); s12 = sum(a[0] * a[1] for a in A); s22 = sum(a[1] * a[1] for a in A)
t1 = sum(a[0] * y for a, y in zip(A, yv)); t2 = sum(a[1] * y for a, y in zip(A, yv))
det = s11 * s22 - s12 * s12
K_fit = (t1 * s22 - t2 * s12) / det
KB = (s11 * t2 - s12 * t1) / det
lauf.pruefe("Rate: n²(x_n − η) → 2π²η/e = 10.4906 (Anpassung K + K·β/n, n ≥ 1000)", abs(K_fit / KONST - 1) < 2e-3,
            f"K = {K_fit:.5f} gegen {KONST:.5f}, β = {KB / K_fit:.2f}")
lauf.pruefe("ohne Anpassung: n²(x_n − η) bei n = 3·10⁵ auf 0.1 % an 2π²η/e", abs(c[300_000] / KONST - 1) < 1e-3, f"{c[300_000]:.5f}")

lauf.daten = {"x_n": {n: f"{werte[n]:.15f}" for n in ns}, "n2_delta": {n: round(c[n], 6) for n in ns},
              "K_fit": K_fit, "beta": KB / K_fit, "K_theorie": KONST}

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].loglog(ns, [werte[n] - ETA for n in ns], "o", color="#b5523b", label="x_n − η (gerechnet)")
ax[0].loglog(ns, [KONST / n ** 2 for n in ns], "-", color="#3d6fb6", lw=1, label="2π²η/(e n²)")
ax[0].set_xlabel("n"); ax[0].set_title("Lösung von ⁿx = n: Abstand zu η"); ax[0].legend(fontsize=8)
ax[1].semilogx(ns, [c[n] for n in ns], "o-", color="#b5523b", ms=4, label="n²(x_n − η)")
ax[1].axhline(KONST, color="#3d6fb6", ls="--", label="2π²η/e = 10.4906")
ax[1].set_ylim(8, 11); ax[1].set_xlabel("n"); ax[1].set_title("Die Konstante kommt vom Engpass bei e"); ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"x_n − η ≈ 2π²η/(e n²): angepasst K = {K_fit:.5f} gegen {KONST:.5f} (β = {KB / K_fit:.2f}); "
                        f"bei n = 3·10⁵: n²(x_n − η) = {c[300_000]:.5f}.")
lauf.schreibe()
