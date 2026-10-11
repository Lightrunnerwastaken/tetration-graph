"""Strang C, Runde 4, c29: Poisson–Jensen mit gemessenen Randdaten — (i) Gesetz R = Green-Funktion (Gittertest wie c17),
(ii) Eck-Bedingung (E) in Randform für die GESAMTE Pol-Antwort: Θ_N := sup_{Schlitz ∪ {g=d}} ω_ζ(λ_N)/g(ζ),
λ_N(η) := (1/N)·log(R_N(η)/R_N(L)) (vorzeichenbehaftet, = 0 in der Ecke).

Satz 12.6: log R_N(ζ) ≤ ω_ζ(log R_N|∂W) − (N+2) g(ζ) für ζ ∉ W̄. Randdaten auf ∂W (obere Hälften, Spiegelung R_N(η̄) = R_N(η)),
Randparameter s = τ(η) ∈ ℝ (L ↦ 0, L̄ ↦ ∞, −ρ ↦ −1, a0+R0 ↦ +1, Spiegelung s ↦ |τ0|²/s), Poisson-Integral über ℝ mit stückweise
linearen Randdaten exakt integriert. Matrix A_f wie c27.
Aufruf: python -I c29_poisson_rand.py gitter        (N = 80, 160; Gitter von c17)
        python -I c29_poisson_rand.py theta 200,400  (Θ_N)
"""
import json
import os
import sys
import time

import numpy as np
from flint import acb, arb, arb_mat, arb_poly, ctx
from mpmath import mp, mpf, mpc, pi as mpi, exp as mexp, log as mlog, arg as marg, conj as mconj

sys.stdout.reconfigure(encoding="utf-8")
HIER = os.path.dirname(os.path.abspath(__file__))
AK = json.load(open(os.path.join(HIER, "..", "a_moden", "daten", "kneser_koeff_lang.json"), encoding="utf-8"))
CACHE = os.path.join(HIER, "c29_cache.json")

# ---------- Geometrie (wie c24/c28) ----------
mp.dps = 30
Lm = mpc("0.31813150520476413", "1.3372357014306895")
RHOm = abs(Lm)
a0m = mpf(AK[0])
R0m = abs(Lm - a0m)


def Mob(z):
    return (z - Lm) / (z - mconj(Lm))


_ta, _tb = marg(Mob(-RHOm)), marg(Mob(a0m + R0m))
LO, HI = sorted([_ta, _tb])
if not (LO < 0 < HI):
    LO, HI = HI, LO + 2 * mpi
ALPHA = (HI - LO) / mpi


def tau(z):
    w = Mob(z) * mexp(-1j * LO)
    a = marg(w)
    if a < 0:
        a += 2 * mpi
    a = min(max(a, mpf(0)), ALPHA * mpi)
    return abs(w) ** (1 / ALPHA) * mexp(1j * a / ALPHA)


TAU0 = tau(mpf(10) ** 30)
T0B2 = float(abs(TAU0) ** 2)


def g_tau(t):
    return float(mlog(abs((t - mconj(TAU0)) / (t - TAU0))))


argL = float(marg(Lm))
bL = float(marg(Lm - a0m))


def randpunkte():
    """obere Hälften: D-Bogen φ = argL + (π−argL)w, A-Bogen ψ = bL(1−w); w geometrisch an der Ecke verdichtet."""
    ws = sorted(set([0.0] + list(np.geomspace(1e-9, 0.1, 40)) + list(np.linspace(0.1, 1.0, 28)[1:])))
    pts = []
    for w in ws:
        if w == 0.0:
            pts.append(("L", w, Lm))
            continue
        pts.append(("D", w, RHOm * mexp(1j * (argL + (float(mpi) - argL) * w)) * (1 + mpf(10) ** -25)))
        pts.append(("A", w, a0m + R0m * mexp(1j * bL * (1 - w)) * 1 + (mexp(1j * bL * (1 - w)) * mpf(10) ** -25)))
    return pts


def s_von(z):
    t = tau(z)
    return float(t.real)


# ---------- Pol-Antworten ----------
def antworten(N, punkte_z, prec=None):
    prec = prec or int(1.4 * N + 300)
    ctx.prec = prec
    a = [arb(s) for s in AK[: N + 2]]
    q = arb_poly(a[: N + 1])
    spalten = []
    P = arb_poly([1])
    for j in range(N + 2):
        c = P.coeffs()
        spalten.append(c + [arb(0)] * (N + 1 - len(c)))
        if j <= N:
            P = (P * q).truncate(N + 1)
    fp = [(k + 1) * a[k + 1] for k in range(N + 1)]
    ex = [arb(1)]
    for k in range(1, N + 1):
        ex.append(ex[-1] / k)
    E = []
    for k in range(N + 1):
        s = ex[k]
        for i in range(1, k + 1):
            s -= fp[i] * E[k - i]
        E.append(s / fp[0])
    J = arb_mat([[spalten[j][k] + (E[k - j] if j <= k else 0) for j in range(N + 1)] for k in range(N + 1)])
    qN1 = arb_poly(spalten[N + 1])
    cols = []
    etas = []
    for z in punkte_z:
        eta = acb(mp.nstr(z.real, 28), mp.nstr(z.imag, 28))
        etas.append(eta)
        d0 = eta - a[0]
        r = [1 / d0]
        for k in range(1, N + 1):
            s = acb(0)
            for j in range(1, k + 1):
                s += a[j] * r[k - j]
            r.append(s / d0)
        vr = (qN1 * arb_poly([x.real for x in r])).truncate(N + 1).coeffs()
        vi = (qN1 * arb_poly([x.imag for x in r])).truncate(N + 1).coeffs()
        cols += [vr + [arb(0)] * (N + 1 - len(vr)), vi + [arb(0)] * (N + 1 - len(vi))]
    RHS = arb_mat([[cols[c][k] for c in range(len(cols))] for k in range(N + 1)])
    X = J.solve(RHS, algorithm="approx")
    rho = float(RHOm)
    out = []
    for c, eta in enumerate(etas):
        skal = 1 / eta ** (N + 1)
        tot = 0.0
        for k in range(N + 1):
            z = acb(X[k, 2 * c], X[k, 2 * c + 1]) * skal
            tot += abs(complex(float(z.real.mid()), float(z.imag.mid()))) * rho ** k
        out.append(tot)
    return out


# ---------- Poisson-Integral über ℝ, stückweise linear ----------
def poisson(s_arr, F_arr, t):
    """(1/π)∫ F(s) Im t/|s−t|² ds, F stückweise linear zwischen den (sortierten) s, ausserhalb konstant (Randwerte)."""
    x, y = float(t.real), float(t.imag)
    th = np.arctan((s_arr - x) / y)
    ln = np.log((s_arr - x) ** 2 + y ** 2)
    dth = np.diff(th)
    dln = np.diff(ln)
    slope = np.diff(F_arr) / np.diff(s_arr)
    tot = np.sum(F_arr[:-1] * dth + slope * ((x - s_arr[:-1]) * dth + 0.5 * y * dln))
    tot += F_arr[0] * (th[0] + np.pi / 2) + F_arr[-1] * (np.pi / 2 - th[-1])
    return float(tot / np.pi)


def randdaten(N):
    key = str(N)
    d = json.load(open(CACHE, encoding="utf-8")) if os.path.isfile(CACHE) else {}
    pts = randpunkte()
    if key in d and len(d[key]) == len(pts):
        R = d[key]
    else:
        t0 = time.time()
        R = antworten(N, [p[2] for p in pts])
        d = json.load(open(CACHE, encoding="utf-8")) if os.path.isfile(CACHE) else {}
        d[key] = R
        json.dump(d, open(CACHE, "w", encoding="utf-8"))
        print(f"   [N={N}: {len(pts)} Randpunkte in {time.time() - t0:.0f} s]")
    s = np.array([0.0 if p[0] == "L" else s_von(p[2]) for p in pts])
    logR = np.log(np.array(R))
    # Spiegelung: s' = |τ0|²/s (L̄-Seite), gleiche Werte
    s_all = np.concatenate([s, T0B2 / s[s != 0]])
    F_all = np.concatenate([logR, logR[s != 0]])
    o = np.argsort(s_all)
    s_all, F_all = s_all[o], F_all[o]
    keep = np.concatenate([[True], np.diff(s_all) > 1e-12 * (1 + np.abs(s_all[1:]))])
    return s_all[keep], F_all[keep], R[0]


modus = sys.argv[1]
if modus == "gitter":
    gitter = [mpc(x, y) for x in (-2.2, -1.6, -1.0, -0.4, 0.2, 0.8, 1.4, 1.9, 2.4) for y in (0.0, 0.6, 1.2, 1.8, 2.4)
              if not (x == 0.2 and y == 0.0)]

    def in_W(z):
        return abs(z) < RHOm or abs(z - a0m) < R0m

    aussen = [z for z in gitter if not in_W(z)]
    daten = {}
    for N in (80, 160):
        sN, FN, RL = randdaten(N)
        RZ = antworten(N, aussen)
        daten[N] = (sN, FN, RZ)
    print("ζ | e^{−g} | Max-Formel | beob. (R160/R80)^{1/80} | Poisson-Jensen-Vorhersage | Jensen-Lücke N=80 / 160")
    abw_g, abw_max, abw_pj = [], [], []
    for i, z in enumerate(aussen):
        t = tau(z)
        g = g_tau(t)
        mx = float(max(RHOm / abs(z), R0m / abs(z - a0m)))
        beob = (daten[160][2][i] / daten[80][2][i]) ** (1 / 80)
        om = {N: poisson(daten[N][0], daten[N][1], t) for N in (80, 160)}
        pj = np.exp((om[160] - om[80]) / 80 - g)
        luecke = {N: om[N] - (N + 2) * g - np.log(daten[N][2][i]) for N in (80, 160)}
        abw_g.append(np.log(beob) + g)
        abw_max.append(np.log(beob / mx))
        abw_pj.append(np.log(beob / pj))
        print(f"{complex(z):+.1f} | {np.exp(-g):.4f} | {mx:.4f} | {beob:.4f} | {pj:.4f} | {luecke[80]:.2f} / {luecke[160]:.2f}")
    for name, a in (("log(beob) + g", abw_g), ("log(beob/Max)", abw_max), ("log(beob/PJ)", abw_pj)):
        a = np.array(a)
        print(f"{name}: Mittel {a.mean():+.4f}, Streuung {a.std():.4f}, min {a.min():+.4f}, max {a.max():+.4f}, n = {len(a)}")
elif modus == "theta":
    for N in [int(x) for x in sys.argv[2].split(",")]:
        sN, FN, RL = randdaten(N)
        lam = (FN - np.log(RL)) / N
        zeile = []
        # Schlitz
        best, arg_best = -1e9, None
        for u in np.geomspace(1e-7, 0.15, 60):
            t = tau(Lm * (1 + mpf(u)))
            r = poisson(sN, lam, t) / g_tau(t)
            if r > best:
                best, arg_best = r, u
        lim = []
        for u in (1e-6, 1e-3, 1e-2, 0.05):
            t = tau(Lm * (1 + mpf(u)))
            lim.append(f"u={u}: {poisson(sN, lam, t) / g_tau(t):+.3f}")
        # Niveaulinien
        nv = []
        for d in (0.05, 0.1):
            k = np.exp(d)
            c = (k ** 2 * complex(TAU0) - np.conj(complex(TAU0))) / (k ** 2 - 1)
            r = k * abs(complex(TAU0) - np.conj(complex(TAU0))) / (k ** 2 - 1)
            m = max(poisson(sN, lam, c + r * np.exp(2j * np.pi * j / 360)) for j in range(360)) / d
            nv.append(f"d={d}: {m:+.3f}")
        # Randprofil (zur Info): λ_N an ausgewählten Stellen
        print(f"N={N}: R_N(L) = {RL:.4g}; Θ_Schlitz = {best:+.3f} (bei u = {arg_best:.2g}); " + ", ".join(lim)
              + "; sup_{g=d} ω(λ)/d: " + ", ".join(nv))
        sys.stdout.flush()
