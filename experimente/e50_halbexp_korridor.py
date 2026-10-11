"""exp.halbexp_korridor — Geometrie von Ω = slog_e(D), D = {|z| < |L|}, und der Korridor V₀ = {z ∈ D : f_K(z) ∈ D}.

 Der Eindeutigkeitssatz kn.satz.halbexp_beschraenkt braucht von Knesers slog = A nur Geometrie:
   (G1) Re A ≤ c_R auf D (Maximum am Rand, bei z = |L|);
   (G2) jede waagrechte Gerade schneidet Ω in einer Strecke (Ω waagrecht konvex), Länge > 1.5 — dann hat B = Ω ∩ (Ω − ½) überall Breite > 1;
   (G3) der linke Rand fällt höchstens linear: x_links(y) ≥ −2 − 0.24|y| (asymptotisch Steigung −Re(1/L)/Im(1/L)·… = −0.2379).
 Oben (y > 0) läuft Ω gegen L, unten (y < 0) spiegelbildlich gegen L̄. Asymptotisch ist A(z) ≈ log(z − L)/L + c: zwei parallele Geraden,
 waagrechter Abstand π/|L|/cos(…) = 2.3495.
 Für die Rückführung kn.satz.polynom_halbexp_reduktion: die Polynome p_N bilden V₀ (auf kompakten Teilen) in D ab.
"""
import mpmath as mp
import numpy as np

from gemeinsam import Lauf, kneser, plt
from polynom_wurzel import horner, loesungen

mp.mp.dps = 40
lauf = Lauf("exp.halbexp_korridor", __file__)
E = "exp(1)"
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
R, phiL = abs(L), mp.arg(L)

# Rand von Ω (obere Hälfte): rechter Bogen φ ∈ [0, φ_L), linker Bogen φ ∈ (φ_L, π]
rechts_phi = [mp.mpf(i) / 60 * (phiL - mp.mpf("0.05")) for i in range(61)] + [phiL - mp.mpf(10) ** (-mp.mpf(k) / 4) for k in range(6, 33)]
links_phi = [phiL + mp.mpf(10) ** (-mp.mpf(k) / 4) for k in range(32, 5, -1)] + \
            [phiL + mp.mpf("0.05") + mp.mpf(i) / 60 * (mp.pi - phiL - mp.mpf("0.05")) for i in range(61)]
wr = kneser("slog", E, [R * mp.expj(p) for p in rechts_phi])
wl = kneser("slog", E, [R * mp.expj(p) for p in links_phi])
xr, yr = np.array([float(mp.re(w)) for w in wr]), np.array([float(mp.im(w)) for w in wr])
xl, yl = np.array([float(mp.re(w)) for w in wl]), np.array([float(mp.im(w)) for w in wl])

c_R = max(xr.max(), xl.max())
lauf.pruefe("(G1) Re slog auf dem Rand von D höchstens c_R = slog(|L|) = 0.31218 (also auf ganz D: Maximumprinzip)",
            abs(c_R - float(mp.re(wr[0]))) < 1e-12 and abs(c_R - 0.312179552787) < 1e-9, round(c_R, 9))
mono = bool(np.all(np.diff(yr) > 0) and np.all(np.diff(yl) < 0))
lauf.pruefe("(G2a) rechter Bogen steigt, linker fällt in Im: jede Höhe trifft jeden Bogen einmal — Ω ist waagrecht konvex", mono)

ys = np.linspace(0, min(yr.max(), yl.max()), 400)
x_rechts = np.interp(ys, yr, xr)
x_links = np.interp(ys, yl[::-1], xl[::-1])
breite = x_rechts - x_links
lauf.pruefe("(G2b) Breite von Ω auf jeder Höhe ≥ 2.0, also Ω ∩ (Ω − ½) breiter als 1 (bis y = 13.6 gerechnet)",
            breite.min() > 2.0, {"min": round(float(breite.min()), 4), "bei y": round(float(ys[breite.argmin()]), 3), "oben": round(float(breite[-1]), 4)})
asym_breite = float(mp.pi / abs(L) / (mp.im(L) / abs(L)))          # Abstand der Geraden π/|L|, geteilt durch die Steigungskomponente
lauf.pruefe("asymptotisch: Breite → π/Im(L) = 2.3495 (zwei parallele Geraden aus A ≈ log(z − L)/L)", abs(breite[-1] - asym_breite) < 5e-3,
            {"π/Im L": round(asym_breite, 5), "gemessen oben": round(float(breite[-1]), 5)})
steig = float(-mp.re(L) / mp.im(L))
lauf.pruefe("(G3) linker Rand x_links(y) ≥ −2 − 0.24·y; Steigung oben → −Re(L)/Im(L) = −0.2379",
            bool(np.all(x_links >= -2 - 0.24 * ys)) and abs((x_links[-1] - x_links[-50]) / (ys[-1] - ys[-50]) - steig) < 5e-3,
            {"Steigung oben": round(float((x_links[-1] - x_links[-50]) / (ys[-1] - ys[-50])), 4), "−Re L/Im L": round(steig, 4)})

# Korridor B = Ω ∩ (Ω − ½): Gitter, T_K(w) ∈ D und T_K(w + ½) ∈ D; Polynome p_N bilden T_K(w) in D ab
hoehen = np.arange(0, 6.01, 0.5)
gitter = []
for y in hoehen:
    a, b = float(np.interp(y, yl[::-1], xl[::-1])), float(np.interp(y, yr, xr)) - 0.5
    for t in np.linspace(0.08, 0.92, 9):
        gitter.append(mp.mpc(a + t * (b - a), y))
z = kneser("sexp", E, gitter)
fk = kneser("sexp", E, [w + mp.mpf("0.5") for w in gitter])
lauf.pruefe("Korridor: für w ∈ Ω ∩ (Ω − ½) liegen z = T_K(w) und f_K(z) = T_K(w + ½) in D (117 Gitterpunkte, 0 ≤ Im w ≤ 6)",
            max(abs(v) for v in z) < R and max(abs(v) for v in fk) < R,
            {"max|z|/|L|": mp.nstr(max(abs(v) for v in z) / R, 6), "max|f_K(z)|/|L|": mp.nstr(max(abs(v) for v in fk) / R, 6)})
loes = loesungen(100)
korr = {N: max(abs(horner(loes[N], v)) for v in z) / R for N in (20, 40, 60, 100)}
abst = {N: max(abs(horner(loes[N], v) - u) for v, u in zip(z, fk)) for N in (20, 40, 60, 100)}
lauf.pruefe("die Polynome p_N bilden dieselben Punkte in D ab: max |p_N(z)|/|L| < 1 (N = 20, 40, 60, 100), Abstand zu f_K fällt",
            all(v < 1 for v in korr.values()) and abst[100] < abst[60] < abst[20],
            {"max|p_N|/|L|": {N: mp.nstr(v, 8) for N, v in korr.items()}, "max|p_N − f_K|": {N: mp.nstr(v, 3) for N, v in abst.items()}})

lauf.daten = {"c_R": c_R, "breite_min": float(breite.min()), "breite_oben": float(breite[-1]), "steigung_links": steig,
              "korridor_max_pN": {N: mp.nstr(v, 10) for N, v in korr.items()}}
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
for s in (1, -1):
    ax[0].plot(xr, s * yr, color="#1f6feb"); ax[0].plot(xl, s * yl, color="#1f6feb")
    ax[0].plot(xr - 0.5, s * yr, color="#d29922", ls="--"); ax[0].plot(xl - 0.5, s * yl, color="#d29922", ls="--")
ax[0].scatter([float(mp.re(w)) for w in gitter], [float(mp.im(w)) for w in gitter], s=4, color="k")
ax[0].set_title("Ω = slog(D) (blau), Ω − ½ (gelb), Korridor-Gitter"); ax[0].set_xlabel("Re w"); ax[0].set_ylabel("Im w")
ax[0].set_xlim(-6, 1.5)
ax[1].plot(ys, breite, label="Breite von Ω"); ax[1].plot(ys, breite - 0.5, label="Breite von Ω ∩ (Ω − ½)")
ax[1].axhline(1, color="k", lw=0.8, ls=":"); ax[1].set_xlabel("Höhe y"); ax[1].legend(fontsize=8); ax[1].set_title("waagrechte Breite")
lauf.bild(fig)
lauf.zusammenfassung = (f"Ω = slog(D) ist waagrecht konvex, Re ≤ {c_R:.5f}, Breite ≥ {breite.min():.3f} (oben → 2.3495), linker Rand mit Steigung −0.238; "
                        "der Korridor Ω ∩ (Ω − ½) ist überall breiter als 1, und die Polynome p_N bilden ihn in D ab.")
lauf.schreibe()
