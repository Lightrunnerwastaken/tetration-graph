"""Strang V — Reserve-Fenster für die Linse (Auftrag Leitung 21:50), Teil (b): geometrische Obergrenze U(τ).

Für Punkte x = ρ·e^{s}·e^{i(θ_L + τ)} über Γin (τ = Winkel ab arg L, 0 < τ < τ_end = 2π − 2θ_L) wird geprüft, bis zu welcher Höhe s
die von Satz K1 benutzte Geometrie hält:
  (G1) f(ℬ) ⊂ D, d. h. |f(x)| < ρ (dann trifft f(ℬ) die Schlitze S_L, S_L̄ nicht, und K_γ = K' ∪ ℬ̄ behält die K'-Ecke);
  (G2) f'(x) ≠ 0 (sonst hat m = f'∘f … bzw. für K1 m = eˣ/f' Pole) — für Kneser über f_K'(x) = f_K'(eˣ)eˣ/f_K(eˣ) automatisch, solange eˣ ∈ D;
  (G3) Lage zum Schlitz: die Linse liegt über Γin, der Schlitz auf dem Strahl arg L — trifft sich nur in L.
Auswertung NUR ausserhalb von D bzw. auf ∂D (kein float im Innern bei grossen Polynomen):
  kubisch: f exakt (Polynom 3. Grades); Kneser: f_K(x) = Log f_K(eˣ) (f∘exp = exp∘f), f_K(eˣ) mit |eˣ| < ρ aus den zertifizierten
  Koeffizienten (641 Terme, |b_k| ≤ 461 k^{−6} ⇒ Abschneidefehler ≤ 1e-12), Zweig von Log durch Stetigkeit ab s = 0 (Kontrolle gegen die Reihe auf ∂D).
Ausgabe: Tabelle U(τ) für beide Modelle und lokale Näherung U ≈ τ·tan(arg μ) bei L.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
GEM = HIER.parents[1] / "gemeinsam"
L = complex(0.31813150520476413, 1.3372357014306895)
rho = abs(L)
thL = math.atan2(L.imag, L.real)
tau_end = 2 * math.pi - 2 * thL
mu = np.sqrt(L)

# kubisch: reell, f(L) = L, f'(L) = μ
A = np.array([[1, L.real, (L * L).real, (L ** 3).real], [0, L.imag, (L * L).imag, (L ** 3).imag],
              [0, 1, (2 * L).real, (3 * L * L).real], [0, 0, (2 * L).imag, (3 * L * L).imag]])
c = np.linalg.solve(A, np.array([L.real, L.imag, mu.real, mu.imag]))
f_kub = np.polynomial.Polynomial(c)
df_kub = f_kub.deriv()
nullst = df_kub.roots()

d = json.loads((GEM / "kneser_halbexp_koeff_zert_640.json").read_text(encoding="utf-8"))
a = np.array([float(x[0]) for x in d["a"]])


def fK_innen(y):
    """Reihe für |y| ≤ ρ (Horner, komplex)."""
    s = np.zeros_like(y, dtype=complex)
    for k in range(len(a) - 1, -1, -1):
        s = s * y + a[k]
    return s


def fK_aussen(x, ref):
    """f_K(x) = Log f_K(eˣ) + 2πik, k so, dass Stetigkeit zum Referenzwert ref gilt."""
    v = np.log(fK_innen(np.exp(x)))
    k = np.round((ref - v).imag / (2 * np.pi))
    return v + 2j * np.pi * k


def U_profil(modell, taus, smax=0.5, ds=1e-4):
    out = []
    for tau in taus:
        ph = thL + tau
        s = 0.0
        x0 = rho * np.exp(1j * ph)
        if modell == "kneser":
            ref = fK_innen(np.array([x0 * (1 - 1e-12)]))[0]
            if abs(np.exp(x0)) >= rho:
                out.append((tau, float("nan"), "eˣ ∉ D"))
                continue
        grenze = None
        while s < smax:
            x = rho * np.exp(s) * np.exp(1j * ph)
            if modell == "kubisch":
                fx = f_kub(x)
            else:
                if abs(np.exp(x)) >= rho:
                    grenze = (s, "eˣ verlässt D")
                    break
                fx = fK_aussen(np.array([x]), ref)[0]
                ref = fx
            if abs(fx) >= rho:
                grenze = (s, "|f| = ρ")
                break
            s += ds
        if grenze is None:
            grenze = (smax, "smax")
        out.append((tau, grenze[0], grenze[1]))
    return out


def main():
    taus = [math.radians(t) for t in (0.5, 1, 2, 4, 8, 15, 30, 45, 60, 80, 103.38)]
    print(f"kubisch: f'-Nullstellen {np.round(nullst, 4)}, |·| = {np.round(np.abs(nullst), 4)}")
    erg = {}
    for modell in ("kubisch", "kneser"):
        res = U_profil(modell, taus)
        erg[modell] = [(math.degrees(t), u, grund) for t, u, grund in res]
        print(f"\n{modell}: U(τ) = sup s mit |f(ρe^{{s}}e^{{i(θ_L+τ)}})| < ρ")
        for t, u, grund in res:
            print(f"  τ = {math.degrees(t):7.2f}°   U = {u:.4f}   ({grund})   τ·tan(arg μ) = {t * math.tan(np.angle(mu)):.4f}")
    (HIER / "v13_U_profil.json").write_text(json.dumps(erg, indent=1, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
