"""exp.forum_formeln — vier Formeln aus dem Tetration-Forum, nachgerechnet.

 1. Fixpunkte von e^z (L_k = −W_k(−1)): Σ_k L_k^(−m) ist rational (J. Fox 2007, tid 94), Erzeugende: Σ_k L_k^(−m) = −[z^(m−1)] (e^z − 1)/(e^z − z).
 2. Newton-Reihe der ganzzahligen Iterierten (Woon/Helms, tid 12, 719): T(h) = Σ_n C(h,n) Δ^n T(0) ergibt für 1 < b < η die reguläre Tetration.
 3. Reihe des endlichen Turms in u = ln x (Galidakis): x↑↑n stimmt bis u^n mit h(u) = Σ (k+1)^(k−1) u^k / k! (unendlicher Turm) überein.
 4. Der unendliche Turm i^i^i^… = (2i/π)·W(−iπ/2) = 0.4382829… + 0.3605924…i (Euler/Macintyre).
"""
import mpmath as mp

from gemeinsam import Lauf, plt
from regulaer import Regulaer

mp.mp.dps = 60
lauf = Lauf("exp.forum_formeln", __file__)

# 1. Fixpunktsummen
K = 4000
L = [-mp.lambertw(-1, k) for k in range(-K, K)]          # k = −K … K−1, paarweise konjugiert (W_k und W_{−k−1})
erz = mp.taylor(lambda z: (mp.exp(z) - 1) / (mp.exp(z) - z), 0, 12)
fox = [-1, -1, 2, 9, -6, -155, -232, 3969, 20870, -118779]   # Jay Fox: Zähler über (m−1)!, m = 2 … 11
exakt = [-erz[m - 1] for m in range(2, 12)]
lauf.pruefe("Fox' Zähler = −(m−1)!·[z^(m−1)] (e^z−1)/(e^z−z) für m = 2 … 11",
            all(abs(e * mp.factorial(m - 1) - f) < mp.mpf(10) ** -40 for m, e, f in zip(range(2, 12), exakt, fox)),
            [mp.nstr(e, 8) for e in exakt[:4]])
summen, abw = {}, 0
for m in range(2, 7):
    s1 = mp.fsum(l ** (-m) for l in L)
    s2 = mp.fsum(l ** (-m) for l in L[K // 2: K + K // 2])     # halbes K: Schwanzschätzung
    fehler = abs(s1 - exakt[m - 2])
    summen[m] = (mp.nstr(mp.re(s1), 15), mp.nstr(exakt[m - 2], 15), mp.nstr(fehler, 3))
    abw = max(abw, fehler / max(abs(s1 - s2), mp.mpf(10) ** -50))
lauf.pruefe("Σ L_k^(−m) über 8000 Fixpunkte trifft −1, −1/2, 1/3, 3/8, −1/20 bis auf den Schwanz (Fehler < Schwanzgrösse)", abw < 1, summen)

# 2. Newton-Reihe für √2
b = mp.sqrt(2)
R = Regulaer(b, N=120)
T = [mp.mpf(1)]
for _ in range(90):
    T.append(mp.power(b, T[-1]))
h = mp.mpf("0.5")
ziel = R.tet(h)
newton, werte = mp.mpf(0), {}
for n in range(0, 91):
    delta = mp.fsum((-1) ** (n - j) * mp.binomial(n, j) * T[j] for j in range(n + 1))
    newton += mp.binomial(h, n) * delta
    if n in (10, 20, 40, 60, 90):
        werte[n] = mp.nstr(abs(newton - ziel), 3)
lauf.pruefe("Newton-Reihe der ganzen Iterierten bei h = ½ → reguläre Tetration √2 (Abstand fällt mit N)",
            mp.mpf(werte[90]) < mp.mpf(10) ** -6 and mp.mpf(werte[90]) < mp.mpf(werte[20]), {"Ziel": mp.nstr(ziel, 20), "Abstand": werte})

# 3. Turmreihe um x = 1
def turm_u(u, n):
    t = mp.exp(u)
    for _ in range(n - 1):
        t = mp.exp(u * t)
    return t


hk = [mp.mpf(k + 1) ** (k - 1) / mp.factorial(k) for k in range(10)]
gleich = {}
for n in range(1, 8):
    a = mp.taylor(lambda u: turm_u(u, n), 0, 9)
    erste_abw = next((k for k in range(10) if abs(a[k] - hk[k]) > mp.mpf(10) ** -40), None)
    gleich[n] = erste_abw
lauf.pruefe("x↑↑n und der unendliche Turm h(ln x) stimmen genau bis (ln x)^n überein, ab (ln x)^(n+1) nicht mehr (n = 1 … 7)",
            all(gleich[n] == n + 1 for n in range(1, 8)), gleich)

# 4. i^i^i^…
lnI = mp.mpc(0, mp.pi / 2)
hi = -mp.lambertw(-lnI) / lnI
z = mp.mpc(1)
bahn = [z]
for _ in range(3000):
    z = mp.power(mp.j, z)
    bahn.append(z)
lam = hi * lnI
lauf.pruefe("i^i^i^… = (2i/π)·W(−iπ/2) = 0.4382829… + 0.3605924…i; Multiplikator |λ| < 1, Iteration konvergiert",
            abs(bahn[-1] - hi) < mp.mpf(10) ** -30 and abs(lam) < 1, {"Wert": mp.nstr(hi, 16), "|λ|": mp.nstr(abs(lam), 8), "arg λ/2π": mp.nstr(mp.arg(lam) / (2 * mp.pi), 8)})

lauf.daten = {"fixpunktsummen": {m: v[1] for m, v in summen.items()}, "newton_abstand": werte, "i_turm": mp.nstr(hi, 25)}
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
ax[0].plot([float(mp.re(w)) for w in bahn[:60]], [float(mp.im(w)) for w in bahn[:60]], "o-", ms=3, lw=0.7)
ax[0].plot([float(mp.re(hi))], [float(mp.im(hi))], "r*", ms=10)
ax[0].set_title("Bahn 1, i, i^i, … → i^i^i^…"); ax[0].set_aspect("equal")
ns = sorted(werte)
ax[1].semilogy(ns, [float(mp.mpf(werte[n])) for n in ns], "o-")
ax[1].set_xlabel("Glieder N"); ax[1].set_title("Newton-Reihe √2 bei h = ½: Abstand zur regulären Tetration")
lauf.bild(fig)
lauf.zusammenfassung = ("Fox' rationale Fixpunktsummen folgen aus (e^z−1)/(e^z−z); die Newton-Reihe der ganzen Iterierten trifft die reguläre Tetration; "
                        "x↑↑n stimmt bis (ln x)^n mit dem unendlichen Turm überein; i^i^i^… = 0.43828 + 0.36059i.")
lauf.schreibe()
