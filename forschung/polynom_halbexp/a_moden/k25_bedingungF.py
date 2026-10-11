"""Strang A, Anfrage Strang C (BRETT 20:56), Bedingung (F) für Satz 9.5:
(1) Wo liegt der Schnitt des Motors (slog) bei L? (slog auf |z − L| = 0.1, Sprünge suchen)
(2) max|f_K| auf der Kontur Γ_d = Rand von W + d (W = D(0,|L|) ∪ D(a₀, R₀)), ohne |ζ − L|, |ζ − L̄| < 0.25, d = 0.05, 0.1, 0.2;
    f_K = sexp(slog(ζ) + ½) vom Motor, Stetigkeit längs der Kontur geprüft (Zweigsprünge melden).
(3) Sprung über den radialen Schlitz ζ = L(1 + u), 0 < u ≤ 0.2: Δf aus der exakten lokalen Formel (Ufer: Log(1 − ζ/L) = ln u ± iπ),
    Δf/u^{5.447} und Vergleich mit dem Leitterm d₁ t^{1+γ}(e^{±2πi(1+γ)} − 1).
(4) Holomorphie im Ring zwischen ∂D und Γ_d: Cauchy-Integrale ∮ f dζ über kleine Kreise (r = 0.1) mit Mittelpunkten im Ring = 0?
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from fatou_backend import FatouGP  # noqa: E402

HIER = Path(__file__).resolve().parent
DATEN = HIER / "daten"
mp.mp.dps = 40
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
G = 2j * mp.pi / L
SL = mp.exp(L / 2)
LML = mp.log(-L)
RL = abs(L)
cx = lambda p: mp.mpc(mp.mpf(p[0]), mp.mpf(p[1]))  # noqa: E731
sth = json.loads((DATEN / "singularitaet_theta.json").read_text(encoding="utf-8"))
phi = [cx(p) for p in sth["phi"]]
psi = [cx(p) for p in sth["psi"]]
eps = [cx(p) for p in sth["eps_q"]]
sing = json.loads((HIER.parent / "gemeinsam" / "kneser_singularitaet_L.json").read_text(encoding="utf-8"))
d1 = cx(sing["d_n"]["1"])


def pv(c, x):
    s = mp.mpc(0)
    for a in reversed(c):
        s = s * x + a
    return s


gp = FatouGP(dps=40, fatou_gp="fork")
S = gp.session("exp(1)")

# (1) Schnitt des Motors um L
n = 360
kreis = [L + mp.mpf("0.1") * mp.expjpi(2 * mp.mpf(j) / n) for j in range(n)]
A = S.slog_batch(kreis)
spr = [(j, abs(A[(j + 1) % n] - A[j])) for j in range(n)]
gross = [(j, x) for j, x in spr if x > 0.3]
print("(1) slog auf |z − L| = 0.1: Sprünge > 0.3 bei Winkel (Grad, gemessen von L aus): "
      + ", ".join(f"{360 * j / n + 0.5:.1f}° (|ΔA| = {mp.nstr(x, 3)})" for j, x in gross))
print(f"    Richtung radial nach aussen = arg L = {float(mp.degrees(mp.arg(L))):.1f}°; radial nach innen = {float(mp.degrees(mp.arg(-L))) % 360:.1f}°")

# (2) Kontur Γ_d
a0 = mp.mpf("0.49856328794111443468")
R0 = abs(L - a0)
for d in (mp.mpf("0.05"), mp.mpf("0.1"), mp.mpf("0.2")):
    pts = []
    m = 2000
    for j in range(m):
        th = 2 * mp.pi * j / m
        # Rand von W + d: in Richtung th der äusserste der beiden Kreise (Strahl von einem inneren Punkt ist ungeeignet);
        # hier: Punkte auf beiden Kreisen, die ausserhalb des jeweils anderen (+d) liegen
        z1 = (RL + d) * mp.expj(th)
        z2 = a0 + (R0 + d) * mp.expj(th)
        for z in (z1, z2):
            if abs(z) >= RL + d - mp.mpf(10) ** -20 and abs(z - a0) >= R0 + d - mp.mpf(10) ** -20 \
                    and abs(z - L) >= mp.mpf("0.25") and abs(z - mp.conj(L)) >= mp.mpf("0.25"):
                pts.append(z)
    pts.sort(key=lambda z: float(mp.arg(z)))
    Ap = S.slog_batch(pts)
    fp = S.sexp_batch([a + mp.mpf("0.5") for a in Ap])
    sprung = max(abs(Ap[i + 1] - Ap[i]) for i in range(len(pts) - 1))
    imax = max(range(len(pts)), key=lambda i: abs(fp[i]))
    print(f"(2) d = {d}: {len(pts)} Punkte, max|f_K| = {mp.nstr(abs(fp[imax]), 6)} bei ζ = {mp.nstr(pts[imax], 5)}; "
          f"grösster Schritt in slog längs der Kontur {mp.nstr(sprung, 3)} (Zweigsprung falls ≫ Schrittweite)")


# (3) Sprung über den radialen Schlitz (exakte lokale Formel)
def f_ufer(u, seite):
    z = L * (1 + u)
    t = z - L
    s = pv(psi, t)
    logs = LML + (mp.log(u) + seite * 1j * mp.pi) + mp.log(s / t)
    qq = mp.exp(G * logs)
    e = mp.fsum(eps[k] * qq ** k for k in range(1, len(eps)))
    return pv(phi, SL * s * mp.exp(L * e)), qq


print("(3) Sprung über ζ = L(1+u):  u, |q| (beide Ufer), |Δf|, |Δf|/u^5.447, Leitterm-Vorhersage |Δf|/u^5.447")
for u in ("0.005", "0.01", "0.02", "0.05", "0.1", "0.15", "0.2"):
    u = mp.mpf(u)
    fo, qo = f_ufer(u, 1)
    fu, qu = f_ufer(u, -1)
    df = fo - fu
    t = L * u
    lead = d1 * (mp.exp((1 + G) * (LML + mp.log(u) + 1j * mp.pi)) - mp.exp((1 + G) * (LML + mp.log(u) - 1j * mp.pi)))
    print(f"    u = {u}: |q| = {mp.nstr(abs(qo), 3)} / {mp.nstr(abs(qu), 3)}, |Δf| = {mp.nstr(abs(df), 4)}, "
          f"|Δf|/u^5.447 = {mp.nstr(abs(df) / u ** mp.mpf('5.447'), 5)}, Leitterm {mp.nstr(abs(lead) / u ** mp.mpf('5.447'), 5)}")
# Kontrolle der Ufer gegen den Motor (Motor-Schnitt liegt anderswo ⇒ beide Ufer sollten Motorwerte nahe dem Schlitz treffen)
for u in ("0.05", "0.2"):
    u = mp.mpf(u)
    nrm = 1j * L / RL
    for seite, sgn in ((1, 1), (-1, -1)):
        z = L * (1 + u) + sgn * mp.mpf(10) ** -12 * nrm
        fm = S.sexp_batch([S.slog_batch([z])[0] + mp.mpf("0.5")])[0]
        fl, _ = f_ufer(u, seite)
        print(f"    Ufer {seite:+d}, u = {u}: |Motor − lokale Formel| = {mp.nstr(abs(fm - fl), 3)}")

# (4) Cauchy-Integrale über kleine Kreise im Ring
print("(4) |∮ f dζ| / (2π r max|f|) über Kreise r = 0.08 um Punkte im Ring |L| < |c| < |L| + 0.15 (ohne Schlitznähe):")
for th in range(0, 360, 30):
    c = (RL + mp.mpf("0.07")) * mp.expj(mp.pi * th / 180)
    if min(abs(c - L), abs(c - mp.conj(L))) < mp.mpf("0.3"):
        continue
    m = 64
    zz = [c + mp.mpf("0.08") * mp.expjpi(2 * mp.mpf(j) / m) for j in range(m)]
    ff = S.sexp_batch([a + mp.mpf("0.5") for a in S.slog_batch(zz)])
    integ = mp.fsum(f * (z - c) for f, z in zip(ff, zz)) * 2j * mp.pi / m
    print(f"    Mittelpunkt arg {th}°: {mp.nstr(abs(integ) / (2 * mp.pi * mp.mpf('0.08') * max(abs(f) for f in ff)), 3)}")
gp.close()
