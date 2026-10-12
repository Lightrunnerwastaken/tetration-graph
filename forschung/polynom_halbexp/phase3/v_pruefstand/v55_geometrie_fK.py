"""Strang V — v55: Ks k54 (Globale Geometrie von E, k_kern §36.2) mit f_K statt der Kubik (Auftrag Leitung 11:40). Kopie von k_kern/k54_geometrie_E.py,
nur f, f′ ersetzt: f_K aus Taylor (A, 4000 Koeff.) in D̄, Ps oberer Struktur (Zweig k₀ = 2, v44) für |x − L| < 0.35 bzw. |x − L̄| < 0.35 (konjugiert),
auf dem Linsenstreifen ausserhalb D fern von L über f(t) = Log f_K(eᵗ) (Zweig nahe f_K(ρe^{iφ})). f′ entsprechend (Taylor-Ableitung, f′(t) = f′(eᵗ)eᵗ/f(eᵗ),
nahe L numerisch, Schrittweite 1e-7).
Ursprünglicher Kopf (K):
𝒞: Linse |t| = ρ(1+s) über Γin (ausserhalb U_{r1}(L, L̄)), Spiralen Σ₁ (u_A = u1, aussen), Σ₂ (u_A = u2, innen) bis r1, Verbindungsbögen, Γout-Bogen |t| = ρ.
E := Abschluss von (int 𝒞) ∪ f(int 𝒞) ∪ (D ∖ U_{r0}(L, L̄)).
Prüft: (1) Γin/Γout-Aufteilung; (2) Univalenz von f auf int 𝒞 (Windungszahl von f(𝒞) um f(z) = 1); (3) E zusammenhängend, ohne Löcher (Raster);
(4) der Spiralschnitt in der Lücke meidet E; (5) Dreipunkt-Bedingung (Ahlfors) auf dem Rand von E; (6) Aussenöffnung bei L numerisch.
Aufruf: python k54_geometrie_E.py [u1 u2 s r1 r0]
"""
import math
import sys

import numpy as np
from matplotlib.path import Path
from scipy import ndimage

sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

L = 0.31813150520476413 + 1.3372357014306895j
RHO = abs(L)
MU = np.sqrt(L)
T = math.radians(1.623)          # Spiralsteigung je e-Faktor
ARG_LAM0 = -71.46


import cmath  # noqa: E402
import json  # noqa: E402
from pathlib import Path as _P  # noqa: E402

HIER = _P(__file__).resolve().parent
sys.path.insert(0, str(HIER / "pv_kopie"))
A4000 = np.array([complex(float(x[0] if isinstance(x, list) else x)) for x in
                  json.loads((HIER.parent.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text(encoding="utf-8"))["a"]])
DA = A4000[1:] * np.arange(1, len(A4000))
from flint import acb, arb, ctx  # noqa: E402
from kc_kneser import Kneser  # noqa: E402

KN = Kneser(160, laut=False)
_R = KN.R
ctx.prec = 160
_Lb = _R.L
_c = KN.c
_TPI = 2 * arb.pi() * acb(0, 1)
_sqL = (_Lb / 2).exp()
_codd = [_c[n] if n % 2 == 1 else acb(0) for n in range(len(_c))]
_tha = _Lb.imag
_pi = arb.pi()


def _mid(z):
    return acb(z.real.mid(), z.imag.mid())


def f_ober(x):
    """f_K nahe L (Schlitzebene, Schnitt S_L), Ps obere Struktur, Zweig k₀ = 2."""
    x = acb(x.real, x.imag)
    w = _mid(x - _Lb)
    for _ in range(60):
        fv, df = _R.P_dP_ball(w)
        dw = _mid((fv - x) / df)
        w = _mid(w - dw)
        if dw.abs_upper() < arb(2) ** (-140):
            break
    A = w.arg()
    while bool(A <= _tha + 4 * _pi):
        A += 2 * _pi
    while bool(A > _tha + 6 * _pi):
        A -= 2 * _pi
    zeta = acb(abs(w).log(), A) / _Lb
    z = zeta - _c[0]
    for _ in range(100):
        zn = _mid(zeta - _R._horner(_c, (_TPI * z).exp()))
        if (zn - z).abs_upper() < arb(2) ** (-140):
            z = zn
            break
        z = zn
    d = -2 * _R._horner(_codd, (_TPI * z).exp())
    v = _R.P_ball(_sqL * w * (_Lb * d).exp())
    return complex(float(v.real.mid()), float(v.imag.mid()))


def _f1(x):
    """f_K an einem Punkt (obere Halbebene oder reell)."""
    if abs(x - L) < 0.35:
        return f_ober(x)
    if abs(x) <= RHO * (1 + 1e-12):
        return complex(np.polyval(A4000[::-1], x))
    # Linsenstreifen ausserhalb D: f(t) = Log f(eᵗ) + 2πik, k so, dass nahe f_K(ρ e^{i arg t})
    e = cmath.exp(x)
    g = cmath.log(complex(np.polyval(A4000[::-1], e)))
    ref = complex(np.polyval(A4000[::-1], RHO * x / abs(x)))
    k = round(((ref - g) / (2j * math.pi)).real)
    return g + 2j * math.pi * k


def f(x):
    x = np.asarray(x, dtype=complex)
    out = np.empty(x.shape, dtype=complex)
    flat, of = x.ravel(), out.ravel()
    inD = (np.abs(flat) <= RHO) & (np.abs(flat - L) >= 0.35) & (np.abs(flat - L.conjugate()) >= 0.35)
    of[inD] = np.polyval(A4000[::-1], flat[inD])
    for i in np.nonzero(~inD)[0]:
        z = flat[i]
        of[i] = _f1(z) if z.imag >= 0 else _f1(z.conjugate()).conjugate()
    return out


def fstrich(x):
    x = np.asarray(x, dtype=complex)
    out = np.empty(x.shape, dtype=complex)
    flat, of = x.ravel(), out.ravel()
    inD = (np.abs(flat) <= RHO) & (np.abs(flat - L) >= 0.35) & (np.abs(flat - L.conjugate()) >= 0.35)
    of[inD] = np.polyval(DA[::-1], flat[inD])
    h = 1e-6
    for i in np.nonzero(~inD)[0]:
        z = flat[i]
        of[i] = (f(np.array([z + h]))[0] - f(np.array([z - h]))[0]) / (2 * h)
    return out


def spiral(uA, rs, Lp=L, conj=False):
    th = np.radians(uA + 71.46) + T * np.log(rs / RHO)
    if conj:
        return (Lp + rs * np.exp(1j * th)).conjugate()
    return Lp + rs * np.exp(1j * th)


def kontur(u1, u2, s, r1, m=500):
    """Geschlossene Kontur, gegen den Uhrzeigersinn um 0: von L über Σ₁ → Linse über Γin → (konj. Σ₁ bei L̄) → L̄ → (konj. Σ₂) → Γout-Bogen → Σ₂ → L."""
    argL = np.angle(L)
    rs = np.geomspace(1e-9, r1, m)
    S1 = spiral(u1, rs)                          # von L nach aussen (r wachsend)
    S2 = spiral(u2, rs)
    # Linse: Kreis |t| = ρ(1+s), Winkel zwischen den Spiral-Endpunkten
    a1 = np.angle(S1[-1]); a1 = a1 if a1 > 0 else a1 + 2 * math.pi
    lens = RHO * (1 + s) * np.exp(1j * np.linspace(a1 + 0.02, 2 * math.pi - a1 - 0.02, 3 * m))
    # Verbindung S1-Ende → Linsenanfang (gerade Strecke), ebenso konjugiert
    con1 = np.linspace(S1[-1], lens[0], 200)
    con1b = np.conj(con1[::-1])
    S1b = np.conj(S1)[::-1]                       # bei L̄: von aussen nach L̄
    # Γout: von L̄-Seite: Σ₂ konjugiert nach aussen, Verbindung zum Kreis, Kreisbogen durch 0°, Verbindung, Σ₂ nach L
    a2 = np.angle(S2[-1])                         # Endwinkel (von 0 aus gesehen) des Σ₂-Endes
    # Punkt auf |t| = ρ nahe dem Σ₂-Ende: radial nach aussen projizieren
    q2 = RHO * np.exp(1j * a2)
    con2 = np.linspace(S2[-1], q2, 200)           # Σ₂-Ende → Kreis
    arc = RHO * np.exp(1j * np.linspace(-a2, a2, 3 * m))   # von L̄-Seite (Winkel −a2) gegen Uhrzeiger … Γout-Bogen durch 0
    S2b_out = np.conj(S2)                         # bei L̄: von L̄ nach aussen
    con2b = np.conj(con2)                         # Σ₂b-Ende → Kreis (konj.)
    parts = [S1, con1, lens, con1b, S1b, S2b_out, con2b, arc, con2[::-1], S2[::-1]]
    names = ["Σ1", "V1", "Linse", "V1b", "Σ1b", "Σ2b", "V2b", "Γout", "V2", "Σ2"]
    return parts, names


def windung(poly, pts):
    path = Path(np.column_stack([poly.real, poly.imag]))
    return path.contains_points(np.column_stack([pts.real, pts.imag]))


def main():
    a = sys.argv[1:]
    u1 = float(a[0]) if len(a) > 0 else 80.0
    u2 = float(a[1]) if len(a) > 1 else -90.0
    s = float(a[2]) if len(a) > 2 else 0.03
    r1 = float(a[3]) if len(a) > 3 else 0.05
    r0 = float(a[4]) if len(a) > 4 else 0.02
    print("Parameter: u1 = %.0f°, u2 = %.0f°, Linse s = %.3f, r1 = %.3f, r0 = %.3f" % (u1, u2, s, r1, r0))
    # (1) Γin/Γout
    ph = np.linspace(0, 2 * math.pi, 3600, endpoint=False)
    fz = np.abs(f(RHO * np.exp(1j * ph)))
    gin = fz < RHO
    wechsel = np.degrees(ph[np.nonzero(np.diff(gin.astype(int)))[0]])
    print("(1) Γin = {|f| < ρ} auf |x| = ρ: Wechsel bei %s Grad; arg L = %.2f°, max|f|/ρ auf Γout = %.3f" % (np.round(wechsel, 2), math.degrees(np.angle(L)), fz[~gin].max() / RHO))
    parts, names = kontur(u1, u2, s, r1)
    Cc = np.concatenate(parts)
    # (2) Univalenz
    xs = np.linspace(-2.2, 2.2, 441)
    X, Y = np.meshgrid(xs, xs)
    Z = X + 1j * Y
    inC = windung(Cc, Z.ravel()).reshape(Z.shape)
    fC = f(Cc)
    zin = Z[inC]
    sel = zin[np.random.default_rng(1).choice(len(zin), 3000, replace=False)]
    # Windungszahl von f(𝒞) um f(z): Summe der Winkeländerungen
    def wind(curve, w):
        d = curve - w
        ang = np.unwrap(np.angle(d))
        return (ang[-1] - ang[0] + np.angle(d[0] / d[-1])) / (2 * math.pi)
    wz = np.array([wind(np.append(fC, fC[0]), f(z)) for z in sel[:600]])
    print("(2) Univalenz: Windungszahl von f(𝒞) um f(z) für 600 innere z: min %.3f, max %.3f; min|f'| auf int 𝒞 (Raster) = %.4f"
          % (wz.min(), wz.max(), np.abs(fstrich(zin)).min()))
    # (3) E auf dem Raster
    inFC = windung(fC, Z.ravel()).reshape(Z.shape)
    inD = (np.abs(Z) <= RHO) & (np.abs(Z - L) >= r0) & (np.abs(Z - L.conjugate()) >= r0)
    E = inC | inFC | inD
    lab, nE = ndimage.label(E)
    lab2, nC = ndimage.label(~E)
    print("(3) Raster 441²: E hat %d Zusammenhangskomponente(n); Komplement hat %d (1 = keine Löcher); Anteil int𝒞 / f(int𝒞) / D-Teil: %d / %d / %d Punkte"
          % (nE, nC, inC.sum(), inFC.sum(), inD.sum()))
    # (4) Schnitt in der Lücke: Spirale in der Mitte zwischen f(Σ₂) und Σ₁
    gap_mid = 0.5 * ((u1) + (u2 + 38.3 + 360))     # in u_A-Grad (Σ₁ bei u1, f(Σ₂) bei u2+38.3)
    gap_mid = gap_mid - 360 if gap_mid > 180 else gap_mid
    rs = np.geomspace(1e-6, 1.0, 4000)
    cut = spiral(gap_mid - 0.0, rs) if False else L + rs * np.exp(1j * (np.radians(0.5 * (u1 + 71.46 + u2 + 71.46 + 38.3)) + T * np.log(rs / RHO)))
    inE_cut = windung(Cc, cut) | windung(fC, cut) | ((np.abs(cut) <= RHO) & (np.abs(cut - L) >= r0))
    hit = rs[inE_cut]
    print("(4) Schnittspirale (Richtung bei r = ρ: %.1f°) trifft E: %s" % (0.5 * (u1 + 71.46 + u2 + 71.46 + 38.3),
          "nein (r ∈ [1e-6, 1])" if len(hit) == 0 else "JA bei r ∈ [%.3g, %.3g] (%d Punkte)" % (hit.min(), hit.max(), len(hit))))
    # (6) Aussenöffnung nahe L: Winkelbereich, der bei Radius r nicht in E liegt
    for r in (0.03, 0.01, 0.001, 1e-4, 1e-6):
        th = np.linspace(-math.pi, math.pi, 7200, endpoint=False)
        pts = L + r * np.exp(1j * th)
        inE = windung(Cc, pts) | windung(fC, pts) | ((np.abs(pts) <= RHO) & (np.abs(pts - L) >= r0))
        frei = np.degrees(th[~inE])
        if len(frei):
            # zusammenhängende Lücke (zyklisch)
            sorted_f = np.sort(np.mod(frei, 360))
            gaps = np.diff(np.concatenate([sorted_f, [sorted_f[0] + 360]]))
            k = np.argmax(gaps)
            start = sorted_f[(k + 1) % len(sorted_f)]; end = sorted_f[k]
            width = (end - start) % 360 + 0.05
            print("(6) r = %.0e: Lücke (nicht in E) von %.1f° bis %.1f°, Breite ≈ %.1f° (Vorhersage Θ = %.1f°)" % (r, start, end % 360, width, u1 - u2 - 38.3))
    # (5) Dreipunkt-Bedingung auf dem Rand von E: Rand = Linse/Σ₁-Teil von 𝒞 und f(Γout-Teil von 𝒞) (+ konj.) — Kontrolle gegen das Raster
    rand = np.concatenate([parts[0][::-1][::5], f(parts[9])[::5], f(parts[8][::-1]), f(parts[7][::-1][::3])])
    print("(5) Dreipunkt-Bedingung (Teilstück nahe L: Σ₁ und f(Σ₂ ∪ V2 ∪ Γout)): siehe Bericht; Stichprobe:", end=" ")
    pts = np.concatenate([parts[0][::40], f(parts[9][::40])])
    worst = 0.0
    seq = np.concatenate([parts[0][::-1], f(parts[9][::-1])[1:], f(parts[8][::-1])[1:], f(parts[7][::-1])[1:]])   # Σ₁ (aussen → L), f(Σ₂) (L → aussen), f(V2), f(Γout)
    keep = np.concatenate([[True], np.abs(np.diff(seq)) > 1e-13]); seq = seq[keep]
    n = len(seq)
    for i in range(0, n, 3):
        d = np.abs(seq[i:] - seq[i])
        cm = np.maximum.accumulate(d)
        ok = d > 1e-14
        if ok.sum() > 1:
            worst = max(worst, np.max(cm[ok] / d[ok]))
    print("max diam(Teilbogen)/Sehne = %.2f" % worst)


if __name__ == "__main__":
    main()
