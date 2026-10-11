"""M18 — Kommutatorterm der IMS-Zerlegung streng (arb): ‖R‖ ≤ Schur-Schranke von |R|, R = K̃ − Σ_j χ_j K̃ χ_j,
Kern R_{ee'}(t,t') = S(t,t')·K̃_{ee'}(t,t'), S = ½Σ_j(χ_j(t) − χ_j(t'))² ≤ min(1, (π²/8W²)(t−t')²), S = 0 im selben flachen Bereich.
K̃_{ee'}(t,t') = √(w_e(t)/w_e'(t')) W_e(t) P_e'(Φ_e(t), t')  (Streifen, unitäres Bild; Φ_0 = G, Φ_1 = H).
Abschneider: χ_L = cos(πs/2), χ_m = sin(πs/2) auf [−T2, −T1] (s = (t+T2)/W), gespiegelt rechts, χ_m = 1 auf [−T1, T1].
Verfahren: Intervallraster der Weite d auf [−Tb, Tb]² (Kugeln für t, t'), Obersummen; Schwänze |t| > Tb oder |t'| > Tb analytisch grob.
Aufruf: python m18_kommutator.py T1 T2 d Tb
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

pi = arb.pi()
eps = arb("0.31382")
th = arb("76.618") * pi / 180
c0 = 1 - 2 * eps * th.cos()
L = acb(th.cos(), th.sin())
Lb = L.conjugate()
eith = L
sinth = th.sin()


def f(z):
    return eps + c0 * z + eps * z * z


def fp(z):
    return c0 + 2 * eps * z


def xi(M):
    return (Lb * M - L) / (M - 1)


def h(x):
    return (-c0 + (c0 * c0 - 4 * eps * (eps - x)).sqrt()) / (2 * eps)


def Lam(x):
    return ((c0 + eps * (x + L)) / (c0 + eps * (x + Lb))).log()


def kante(t, e):
    """Für Kugel t: (Φ_e(t) als acb, W_e(t) acb, w_e(t) arb)."""
    q = t.exp()
    M = eith * q if e == 0 else -eith * q
    x = xi(M)
    w = (sinth / pi) * q / abs(1 - M) ** 2
    if e == 0:
        Phi = acb(t) + Lam(x)
        Wt = 1 / fp(f(x))
    else:
        Phi = acb(t, pi) - Lam(h(x))
        Wt = fp(x)
    return Phi, Wt, w


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
    """Einschluss von (Φ_e, W_e, w_e) über die Kugel tb durch Unterteilung in NSUB Teilkugeln und Vereinigung."""
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
