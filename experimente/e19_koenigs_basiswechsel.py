"""exp.koenigs_basiswechsel — linearer Basiswechsel über Koenigs-Koordinaten (Rekonstruktion).

Basen mit anziehendem Fixpunkt, Multiplikator λ: Log b = λ e^{−λ}, Fixpunkt e^λ.
Ψ = σ_b^{-1} ∘ (κ σ_d), κ = σ_b(1)/σ_d(1). Dann
  H_b(Ψ(T_d(x))) = (ln λ_d / ln λ_b) · x       (Koenigs-Höhe H_b = log_{λ_b}(σ_b/σ_b(1)))
und für λ_d = λ_b^k:  f_b^{∘k} ∘ Ψ = Ψ ∘ f_d  („gleiche Spirale, jeder k-te Punkt“).
"""
import mpmath as mp

from gemeinsam import Lauf, plt
from koenigs_komplex import KoenigsKomplex

mp.mp.dps = 30
lauf = Lauf("exp.koenigs_basiswechsel", __file__)


def psi(B, D, z):
    kappa = B.sigma(mp.mpc(1)) / D.sigma(mp.mpc(1))
    w = kappa * D.sigma(z)
    return B.sigma_inv(w, B.L + w)


def hoehe(B, z):
    return mp.log(B.sigma(z) / B.sigma(mp.mpc(1))) / B.log_lam


# ------------------------------------------------------------------ reell: √2 und λ_d = (ln 2)^2
lb = mp.log(2)
B = KoenigsKomplex(lb, N=40)
D = KoenigsKomplex(lb ** 2, N=40)
lauf.pruefe("Rezeptur: λ = ln 2 liefert die Basis √2 mit Fixpunkt 2", abs(mp.exp(B.u) - mp.sqrt(2)) < 1e-25 and abs(B.L - 2) < 1e-25,
            mp.nstr(mp.re(mp.exp(B.u)), 12))
d_basis = mp.re(mp.exp(D.u))
res = 0
for n in (6, 9, 12):
    y = mp.mpc(1)
    for _ in range(n):
        y = D.f(y)
    links = B.f(B.f(psi(B, D, y)))
    rechts = psi(B, D, D.f(y))
    res = max(res, abs(links - rechts))
lauf.pruefe(f"Konjugation f_√2 ∘ f_√2 ∘ Ψ = Ψ ∘ f_d (d = {mp.nstr(d_basis, 8)}) auf dem Orbit-Schwanz", res < mp.mpf(10) ** -18, mp.nstr(res, 3))
res = 0
for x in (mp.mpf("0.5"), mp.mpf("1.25"), mp.mpf("3.5")):
    res = max(res, abs(hoehe(B, psi(B, D, D.tet(x))) - 2 * x))
lauf.pruefe("Koenigs-Höhe: H_√2(Ψ(T_d(x))) = 2x exakt (x = 0.5, 1.25, 3.5)", res < mp.mpf(10) ** -15, mp.nstr(res, 3))
D2 = KoenigsKomplex(lb ** mp.mpf("1.7"), N=40)
res = max(abs(hoehe(B, psi(B, D2, D2.tet(x))) - mp.mpf("1.7") * x) for x in (mp.mpf("0.5"), mp.mpf("2.2")))
lauf.pruefe("beliebiges Verhältnis: λ_d = λ_b^{1.7} ⇒ Höhe linear mit Steigung 1.7", res < mp.mpf(10) ** -15, mp.nstr(res, 3))
naiv = [mp.re(hoehe(B, D.tet(mp.mpf(x)))) for x in (5, 20, 40)]
lauf.pruefe("naiver Vergleich der Werte (Höhe_√2 von T_d(x)) läuft gegen eine Konstante — kein Basiswechsel im Phasensinn", abs(naiv[2] - naiv[1]) < 1e-5,
            [mp.nstr(v, 8) for v in naiv])

# ------------------------------------------------------------------ komplex: Spiralen
lb_c = mp.mpc(0.85 * mp.cos(1.1), 0.85 * mp.sin(1.1))
Bc = KoenigsKomplex(lb_c, N=40)
Dc = KoenigsKomplex(lb_c ** 2, N=40, log_lam=2 * mp.log(lb_c))
orbit_b, orbit_d = [mp.mpc(1)], [mp.mpc(1)]
for _ in range(60):
    orbit_b.append(Bc.f(orbit_b[-1]))
for _ in range(30):
    orbit_d.append(Dc.f(orbit_d[-1]))
abw = max(abs(psi(Bc, Dc, orbit_d[n]) - orbit_b[2 * n]) for n in range(10, 26))
lauf.pruefe("komplexe Spiralen: Ψ(z_{d,n}) = z_{b,2n} für n = 10…25", abw < mp.mpf(10) ** -15, mp.nstr(abw, 3))
lauf.daten = {"d_basis_reell": mp.nstr(d_basis, 15), "d_fixpunkt": mp.nstr(mp.re(D.L), 15), "naiv": [mp.nstr(v, 10) for v in naiv],
              "b_komplex": mp.nstr(mp.exp(Bc.u), 10), "d_komplex": mp.nstr(mp.exp(Dc.u), 10)}

fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.6))
ax[0].plot([float(z.real) for z in orbit_b], [float(z.imag) for z in orbit_b], ".-", lw=0.6, ms=4, label=f"Orbit Basis b (λ = 0.85·e^(1.1i))")
ax[0].plot([float(z.real) for z in orbit_d], [float(z.imag) for z in orbit_d], ".-", lw=0.6, ms=4, label="Orbit Basis d (λ_d = λ²)")
mapped = [psi(Bc, Dc, orbit_d[n]) for n in range(6, 30)]
ax[0].plot([float(z.real) for z in mapped], [float(z.imag) for z in mapped], "o", mfc="none", ms=8, color="#b5523b", label="Ψ(Orbit d) = jeder 2. Punkt von b")
ax[0].set_aspect("equal"); ax[0].legend(fontsize=7); ax[0].set_title("Zwei Basen, eine Spirale")
xs = [mp.mpf(i) / 10 for i in range(0, 41)]
ax[1].plot([float(x) for x in xs], [float(mp.re(hoehe(B, psi(B, D, D.tet(x))))) for x in xs], label="Koenigs-Höhe H_√2(Ψ(T_d(x)))")
ax[1].plot([float(x) for x in xs], [float(2 * x) for x in xs], "--", color="#888", label="2x")
ax[1].set_xlabel("Höhe x zur Basis d"); ax[1].set_title("Basiswechsel in Koenigs-Koordinaten: exakt linear")
ax[1].legend(fontsize=8)
lauf.bild(fig)
lauf.zusammenfassung = (f"√2 und d = {mp.nstr(d_basis, 10)} (λ_d = (ln 2)²): Ψ konjugiert f_d zu f_√2∘f_√2; Koenigs-Höhe exakt 2x; "
                        f"beliebiges Verhältnis 1.7 ebenfalls linear. Komplexe Spiralen: Ψ bildet den d-Orbit auf jeden zweiten b-Punkt ab.")
lauf.schreibe()
