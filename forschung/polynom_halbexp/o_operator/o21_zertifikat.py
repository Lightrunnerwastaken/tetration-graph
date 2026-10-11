"""O21 — Strenge (arb) Zertifizierung der Pfadbedingungen für die Index-Homotopie (Satz O7b):
f_t = (1−t) f_M + t f_Z, t ∈ [0,1], f_M = Möbius-exp (Strang E), f_Z = kubisch | radial. Alle f_t sind reell und fixieren L, L̄ exakt.
Wegen Symmetrie (f_t reell) genügt φ ∈ [0, π]. Geprüft für ALLE (t, φ) per adaptiver Bisektion in Ball-Arithmetik:
 (P1) Re f_t'(ρe^{iφ}) > 0          ⇒ (Minimumprinzip, Noshiro–Warschawski) f_t univalent auf einer Umgebung von D̄
 (P2) Re m_t(ρe^{iφ}) > 0, m_t = f_t'∘f_t   ⇒ m_t nullstellenfrei auf D̄
 (P3) Re f_t(ρe^{iφ}) < p (Pol von f_M)     ⇒ f_t(D̄) ∌ p, also f_t, m_t holomorph um D̄
 (P4) g_t(φ) = |f_t(ρe^{iφ})|² − ρ²: > 0 auf [0, θ_L − η], < 0 auf [θ_L + η, π], g_t' < 0 auf [θ_L − η, θ_L + η]  (⇒ Γin/Γout-Struktur)
 (P5) X_t := π·Im(Log μ_t)/|Log μ_t|² > 3/2 und 0 < arg μ_t < π, μ_t = f_t'(L)  ⇒ 1 + μ_t^{s−1} ≠ 0 auf Re s = −½ (keine Nullstelle kreuzt)
Aufruf: python o21_zertifikat.py kubisch|radial
"""
import sys
import time

from flint import acb, arb, ctx

import o_kern as K

sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 120
ziel = sys.argv[1] if len(sys.argv) > 1 else "kubisch"
A, B = arb(K.L_RE), arb(K.L_IM)
L = acb(A, B)
rho = (A * A + B * B).sqrt()
thL = arb.atan2(B, A)
# Möbius-exp
al, be, ga, de, cM = K.moeb_ell_par()
pol = -de / ga
zc = K.kubik_koeff() if ziel == "kubisch" else K.radial_koeff()


def fM(x):
    return (acb(al) * x + acb(be)) / (acb(ga) * x + acb(de))


def fMp(x):
    return acb(al * de - be * ga) / (acb(ga) * x + acb(de)) ** 2


def fZ(x):
    s = acb(0)
    for c in reversed(zc):
        s = s * x + acb(c)
    return s


def fZp(x):
    s = acb(0)
    for k in range(len(zc) - 1, 0, -1):
        s = s * x + acb(zc[k] * k)
    return s


def ft(t, x):
    return acb(1 - t) * fM(x) + acb(t) * fZ(x)


def ftp(t, x):
    return acb(1 - t) * fMp(x) + acb(t) * fZp(x)


def ball(lo, hi):
    # Ball, der [lo, hi] sicher enthält (Rundung von Mittelpunkt/Radius durch Zuschlag abgefangen)
    return arb((lo + hi) / 2, (hi - lo) / 2 * (1 + 2.0 ** -40) + 2.0 ** -60)


stat = {"evals": 0}


def pruefe(bed, t0, t1, p0, p1, tiefe=0):
    """Adaptive Bisektion: True, wenn bed(t-Ball, φ-Ball) auf dem ganzen Rechteck zertifiziert."""
    stat["evals"] += 1
    t, ph = ball(t0, t1), ball(p0, p1)
    r = bed(t, ph)
    if r is True:
        return True
    if tiefe > 40:
        raise RuntimeError(f"nicht zertifizierbar bei t∈[{t0},{t1}], φ∈[{p0},{p1}]")
    # Bisektion in der breiteren (relativen) Richtung
    if (t1 - t0) >= (p1 - p0) / 3:
        tm = (t0 + t1) / 2
        return pruefe(bed, t0, tm, p0, p1, tiefe + 1) and pruefe(bed, tm, t1, p0, p1, tiefe + 1)
    pm = (p0 + p1) / 2
    return pruefe(bed, t0, t1, p0, pm, tiefe + 1) and pruefe(bed, t0, t1, pm, p1, tiefe + 1)


def x_of(ph):
    return acb(rho) * acb(ph.cos(), ph.sin())


def P1(t, ph):
    return ftp(t, x_of(ph)).real > 0


def P2(t, ph):
    y = ft(t, x_of(ph))
    return ftp(t, y).real > 0


def P3(t, ph):
    return ft(t, x_of(ph)).real < pol


def g(t, ph):
    y = ft(t, x_of(ph))
    return y.real ** 2 + y.imag ** 2 - rho ** 2


def gprime(t, ph):
    x = x_of(ph)
    y = ft(t, x)
    yp = ftp(t, x) * acb(0, 1) * x
    return 2 * (y.real * yp.real + y.imag * yp.imag)


eta = arb("0.02")
th = thL
t0 = time.time()
ergebnisse = []
for name, bed, lo, hi in (
        ("P1 Re f_t' > 0", P1, arb(0), arb.pi()),
        ("P2 Re m_t > 0", P2, arb(0), arb.pi()),
        ("P3 Re f_t < p", P3, arb(0), arb.pi()),
        ("P4a g_t > 0 auf [0, θ_L−η]", lambda t, ph: g(t, ph) > 0, arb(0), th - eta),
        ("P4b g_t < 0 auf [θ_L+η, π]", lambda t, ph: g(t, ph) < 0, th + eta, arb.pi()),
        ("P4c g_t' < 0 auf [θ_L−η, θ_L+η]", lambda t, ph: gprime(t, ph) < 0, th - eta, th + eta)):
    stat["evals"] = 0
    # Intervallgrenzen als Mittelpunkte (Bälle mit kleinem Radius): wir nehmen die oberen/unteren Schranken konservativ
    lo_f = float(lo.lower()) if lo.rad() > 0 else float(lo)
    hi_f = float(hi.upper()) if hi.rad() > 0 else float(hi)
    if name.startswith("P4b"):
        lo_f = float((th + eta).lower())
    if name.startswith("P4a"):
        hi_f = float((th - eta).upper())
    if name.startswith("P4c"):
        lo_f, hi_f = float((th - eta).lower()), float((th + eta).upper())
    # Startgitter 8 × 8
    ok = True
    for i in range(8):
        for j in range(8):
            ok = ok and pruefe(bed, i / 8, (i + 1) / 8, lo_f + (hi_f - lo_f) * j / 8, lo_f + (hi_f - lo_f) * (j + 1) / 8)
    ergebnisse.append((name, ok, stat["evals"]))
    print(f"{ziel}: {name}: {'ZERTIFIZIERT' if ok else 'FEHLER'} ({stat['evals']} Rechtecke, {time.time()-t0:.0f}s)", flush=True)

# P5: μ_t = f_t'(L) (exakt) für t-Bälle
okP5 = True
minX = None
for i in range(256):
    t = ball(i / 256, (i + 1) / 256)
    mu = ftp(t, L)
    lg = mu.log()
    a_, b_ = lg.real, lg.imag
    X = arb.pi() * b_ / (a_ * a_ + b_ * b_)
    okP5 = okP5 and bool(X > arb("1.5")) and bool(b_ > 0) and bool(b_ < arb.pi())
    minX = X.lower() if minX is None else min(minX, X.lower())
print(f"{ziel}: P5 X_t > 3/2 für alle t: {'ZERTIFIZIERT' if okP5 else 'FEHLER'}; min X_t ≥ {float(minX):.4f}")
print(f"θ_L = {th}, η = {eta}, Pol p = {pol}")
