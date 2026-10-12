"""M23 — wie M18, aber für KNESER (Kantengrössen aus mk_kneser_arb/mk_ecke_kneser; Quelle Vorzertifikat oder Datei von P).
|t| ≤ 13: F, F' auf Kreisbögen (Taylor + Rest), h per parametrischem Krawczyk; |t| > 13: Eckschranken (S2), rechts per Symmetrie.
---- Ursprüngliche Beschreibung (M18):
M18 — Kommutatorterm der IMS-Zerlegung streng (arb): ‖R‖ ≤ Schur-Schranke von |R|, R = K̃ − Σ_j χ_j K̃ χ_j,
Kern R_{ee'}(t,t') = S(t,t')·K̃_{ee'}(t,t'), S = ½Σ_j(χ_j(t) − χ_j(t'))² ≤ min(1, (π²/8W²)(t−t')²), S = 0 im selben flachen Bereich.
K̃_{ee'}(t,t') = √(w_e(t)/w_e'(t')) W_e(t) P_e'(Φ_e(t), t')  (Streifen, unitäres Bild; Φ_0 = G, Φ_1 = H).
Abschneider: χ_L = cos(πs/2), χ_m = sin(πs/2) auf [−T2, −T1] (s = (t+T2)/W), gespiegelt rechts, χ_m = 1 auf [−T1, T1].
Verfahren: Intervallraster der Weite d auf [−Tb, Tb]² (Kugeln für t, t'), Obersummen; Schwänze |t| > Tb oder |t'| > Tb analytisch grob.
Aufruf: python m23_kneser_kommutator.py T1 T2 d Tb [quelle]
"""
import sys
import time

from flint import acb, arb, ctx

sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 90
T1 = arb(sys.argv[1]) if len(sys.argv) > 1 else arb(5)
T2 = arb(sys.argv[2]) if len(sys.argv) > 2 else arb(13)
d = arb(sys.argv[3]) if len(sys.argv) > 3 else arb("0.1")
Tb = arb(sys.argv[4]) if len(sys.argv) > 4 else arb(28)
W = T2 - T1

import numpy as np

from mk_ecke_kneser import eckschranken
from mk_kneser import h_newton
from mk_kneser_arb import KneserArb

quelle = sys.argv[5] if len(sys.argv) > 5 else "vorzertifikat"
KN = KneserArb(quelle)
pi = arb.pi()
th = KN.th
L = KN.eith
Lb = L.conjugate()
eith = L
sinth = th.sin()
ell = KN.ell
mu = KN.mu
Ta = arb(13)
emith = acb(th.cos(), -th.sin())


def xi(M):
    return Lb + (Lb - L) / (M - 1)


def xi_inv(x):
    z = ((x - L) / (x - Lb) * emith).log()
    assert bool(z.imag > 0) and bool(z.imag < pi), "Streifenzweig verletzt"
    return z


def winkel(x):
    """Winkelkugel einer Kugel x nahe dem Einheitskreis, ohne Sprung am Verzweigungsschnitt von arg."""
    x0 = acb(x.real.mid(), x.imag.mid())
    r = x0 / abs(x0)
    a0 = r.arg()
    return (x * r.conjugate()).arg() + a0


def kante_mitte(t, e):
    q = t.exp()
    M = eith * q if e == 0 else -eith * q
    x = xi(M)
    phi = winkel(x)
    w = (sinth / pi) * q / abs(1 - M) ** 2
    if e == 0:
        Phi = xi_inv(KN.F_kreis(phi))
        Wt = KN.Fp_kreis(phi) * (-KN.rho * x).exp()
    else:
        z0 = complex(float(x.real.mid()), float(x.imag.mid()))
        v0 = complex(h_newton(np.array([z0]))[0])
        v = KN.h(acb(phi.cos(), phi.sin()), v0)
        assert v is not None, "Krawczyk gescheitert"
        Phi = xi_inv(v)
        Wt = KN.Fp_kreis(phi)
    return Phi, Wt, w


def kante_links(t, e):
    """t ≤ −Ta: Eckschranken."""
    q = arb(t.upper()).exp()
    E = eckschranken(KN, q)
    Q = t.exp()
    M = eith * Q if e == 0 else -eith * Q
    w = (sinth / pi) * Q / abs(1 - M) ** 2
    if e == 0:
        g = E["gG"]
        Phi = acb(t) + ell + acb(arb(0, g.upper()), arb(0, g.upper()))
        g = E["gWin"]
        Wt = 1 / mu + acb(arb(0, g.upper()), arb(0, g.upper()))
    else:
        g = E["gH"]
        Phi = acb(t, pi) - ell + acb(arb(0, g.upper()), arb(0, g.upper()))
        g = E["gWout"]
        Wt = mu + acb(arb(0, g.upper()), arb(0, g.upper()))
    return Phi, Wt, w


def kante(t, e):
    if bool(t <= -Ta):
        return kante_links(t, e)
    if bool(t >= Ta):
        Ph, Wt, w = kante_links(-t, e)
        return acb(-Ph.real, Ph.imag), Wt.conjugate(), w
    return kante_mitte(t, e)


def poisson(e, p, tp):
    u, y = p.real, p.imag
    den = (tp - u).cosh() - y.cos() if e == 0 else (tp - u).cosh() + y.cos()
    return y.sin() / den / (2 * pi)


def profil(t):
    """Kugel-Einschlüsse (χ_L, χ_m, χ_R) für Kugel t; exakt in den flachen Bereichen."""
    lo, hi = arb(t.lower()), arb(t.upper())
    if bool(hi <= -T2):
        return arb(1), arb(0), arb(0)
    if bool(lo >= -T1) and bool(hi <= T1):
        return arb(0), arb(1), arb(0)
    if bool(lo >= T2):
        return arb(0), arb(0), arb(1)

    def stueck(a, b):   # s ∈ [a, b] ⊂ [0, 1] → Einschlüsse von cos(πs/2), sin(πs/2)
        ca, cb = (pi * a / 2).cos(), (pi * b / 2).cos()
        sa_, sb_ = (pi * a / 2).sin(), (pi * b / 2).sin()
        c = arb((ca.upper() + cb.lower()) / 2, ((ca.upper() - cb.lower()) / 2).upper())
        s_ = arb((sb_.upper() + sa_.lower()) / 2, ((sb_.upper() - sa_.lower()) / 2).upper())
        return c, s_

    def kl(v):
        return arb(min(max(v, arb(0)), arb(1)))

    cL, mL = stueck(kl(((lo + T2) / W).lower()), kl(((hi + T2) / W).upper()))
    cR, mR = stueck(kl(((T2 - hi) / W).lower()), kl(((T2 - lo) / W).upper()))
    if bool(hi < 0):
        return cL, mL, arb(0)
    if bool(lo > 0):
        return arb(0), mR, cR
    return arb(0), arb(1), arb(0)


def S_bound(t, tp):
    a = profil(t)
    b = profil(tp)
    s = sum((x - y) ** 2 for x, y in zip(a, b)) / 2
    lip = (pi ** 2 / (8 * W ** 2)) * (t - tp) ** 2
    up = min(s.upper(), lip.upper(), arb(1).upper())
    return arb(up)


t0 = time.time()
N = int(float((2 * Tb / d).upper()) + 0.5)
mids = [-Tb + d * (arb(i) + arb("0.5")) for i in range(N)]
balls = [arb(m.mid(), (d / 2).upper()) for m in mids]
NSUB = 24


def kante_fein(tb, e):
    """Einschluss von (Φ_e, W_e, w_e) über die Kugel tb durch Unterteilung in NSUB Teilkugeln und Vereinigung (Ecken: direkt)."""
    if bool(arb(tb.upper()) <= -Ta) or bool(arb(tb.lower()) >= Ta):
        return kante(tb, e)
    lo = arb(tb.lower())
    wd = (arb(tb.upper()) - lo) / NSUB
    res = None
    for k in range(NSUB):
        sub = arb((lo + wd * (arb(k) + arb("0.5"))).mid(), (wd / 2).upper() * (1 + 1e-15))
        r = kante(sub, e)
        res = r if res is None else tuple(x.union(y) for x, y in zip(res, r))
    return res


kanten = {e: [kante_fein(tb, e) for tb in balls] for e in (0, 1)}
sqw = {e: [k[2].sqrt() for k in kanten[e]] for e in (0, 1)}
profs = [profil(tb) for tb in balls]


def flach(tb):
    """Index des flachen Bereichs (0 links, 1 Mitte, 2 rechts) oder None."""
    lo, hi = arb(tb.lower()), arb(tb.upper())
    if bool(hi <= -T2):
        return 0
    if bool(lo >= -T1) and bool(hi <= T1):
        return 1
    if bool(lo >= T2):
        return 2
    return None


fl = [flach(tb) for tb in balls]
print(f"Raster: {N} Intervalle je Achse, Weite {d}, [−{Tb}, {Tb}]; flache Kästchen: {sum(x is not None for x in fl)} ({time.time()-t0:.0f}s)",
      flush=True)
zeile = {e: [arb(0)] * N for e in (0, 1)}
spalte = {e: [arb(0)] * N for e in (0, 1)}
lipc = pi ** 2 / (8 * W ** 2)
for i in range(N):
    for j in range(N):
        if fl[i] is not None and fl[i] == fl[j]:
            continue
        # S-Schranke
        s = sum((x - y) ** 2 for x, y in zip(profs[i], profs[j])) / 2
        dt = abs(mids[i] - mids[j]) + d
        Sij = min(s.upper(), (lipc * dt ** 2).upper(), arb(1).upper())
        if Sij == 0:
            continue
        Sij = arb(Sij)
        for e in (0, 1):
            Phi, Wt, w = kanten[e][i]
            for ep in (0, 1):
                k = (sqw[e][i] / sqw[ep][j]) * abs(Wt) * poisson(ep, Phi, balls[j])
                v = arb((Sij * abs(k)).upper()) * d
                zeile[e][i] += v
                spalte[ep][j] += v
    if i % 50 == 0:
        print(f"  Zeile {i}/{N} ({time.time()-t0:.0f}s)", flush=True)
Rz = max(max(z.upper() for z in zeile[e]) for e in (0, 1))
Cs = max(max(z.upper() for z in spalte[e]) for e in (0, 1))
print(f"Zeilensumme ≤ {float(Rz):.5f}, Spaltensumme ≤ {float(Cs):.5f} (ohne Schwänze)")
# Schwänze: Paare mit |t| > Tb oder |t'| > Tb tragen nur bei, wenn der andere Punkt in [−T2, T2] liegt (sonst S = 0),
# also |Δ| ≥ Tb − T2. Grobe Majorante: |K̃| ≤ e^{|Δ|/2}·C_w·Wmax·(1/2π)/(cosh(|Δ| − 1) − 1) für |Δ| ≥ 2, mit
# C_w = sup √(w(t)/w(t'))e^{−|Δ|/2} ≤ sup_t (ρ-Faktoren) ≤ 1/(1−cos θ) grob; Integral ∫_{D}^{∞} e^{Δ/2}·4e^{−(Δ−1)}/(2π) dΔ ·2.
Dm = Tb - T2
Wmax = arb(3)
Cw = (2 / (1 - th.cos())).sqrt()
# Faktoren: 2 Kanten e' × 2 Seiten (t' < −Tb und t' > Tb) × Integral ∫_Dm^∞ e^{−x/2}dx = 2e^{−Dm/2}
schwanz = 2 * 2 * Cw * Wmax * 4 * arb(1).exp() / (2 * pi) * 2 * (-Dm / 2).exp()
# Voraussetzungen der Schwanzmajorante P ≤ (1/2π)·4e·e^{−|Δ|} (|Δ| ≥ 15): |Re Φ_e(t) − t| ≤ 1 und |W_e| ≤ 3 für alle t;
# auf dem Raster hier geprüft, ausserhalb (|t| > Tb) aus den Schwarz-Schranken von m17 (|ReΦ − t| ≤ Re ℓ + γ e^{−Tb}).
vor_ok = True
for e in (0, 1):
    for tb, (Ph, Wt, w) in zip(balls, kanten[e]):
        if not (bool(abs(Ph.real - tb) < 1) and bool(abs(Wt) < 3)):
            vor_ok = False
print(f"Voraussetzungen der Schwanzmajorante auf dem Raster: {'erfüllt' if vor_ok else 'VERLETZT'}")
print(f"Schwanzbeitrag ≤ {float(schwanz.upper()):.2e}")
norm = ((arb(Rz) + schwanz) * (arb(Cs) + schwanz)).sqrt()
print(f"‖R‖ ≤ {float(norm.upper()):.5f}  (T1 = {T1}, T2 = {T2}, W = {W}) ({time.time()-t0:.0f}s)")
