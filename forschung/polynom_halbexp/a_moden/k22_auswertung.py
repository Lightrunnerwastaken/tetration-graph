"""Strang A, Runde 2: Auswertung von P_N = Π_N f_K (aus daten/lin_PN.json).

(1) S_a(N) = max_{z∈D} |P_N(z)|·(|z−L||z−L̄|)^a, a = 0, 1, 2, 3, 3.4 (Maximum auf |z| = |L|, da |P_N g^a| subharmonisch;
    grobes FFT-Gitter + feines Gitter um arg L), mit Ort t* = θ* − arg L (in Einheiten 1/N).
(2) P_N(L) − L.  (3) Modenamplitude c_N (Koeffizient von e^{−2πiw}, Linie w ∈ [−1.6, −0.6)), K_N = c_N(δ_{−1})_N/(N u_N).
(4) Eckschicht: E_N(z) = (P_N − f_K)(z) bei z = L(1 − v/N) für v aus einem Gitter (f_K nahe L aus der exakten lokalen Formel).
Ausgabe: daten/lin_auswertung.json, gemeinsam/eckschicht_profil.json.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from flint import acb, acb_poly, arb, ctx  # noqa: E402

from lin_kern import DATEN, GEM, lade  # noqa: E402

mp.mp.dps = 50
ctx.prec = 200
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
sing = json.loads((GEM / "kneser_singularitaet_L.json").read_text(encoding="utf-8"))
d1, k0 = cx(sing["d_n"]["1"]), cx(sing["kappa"]["0"])
D1 = d1 * mp.exp((1 + G) * LML)
Dd = L * SL * mp.exp(-2j * mp.pi * k0) * mp.exp((1 - G) * LML)
q = [mp.mpf(x) for x in json.loads((GEM / "kneser_halbexp_koeff_1000.json").read_text(encoding="utf-8"))["a"]]
qpoly = acb_poly([acb(arb(mp.nstr(x, 60))) for x in q])


def pv(c, x):
    s = mp.mpc(0)
    for a in reversed(c):
        s = s * x + a
    return s


def fK(z):
    """f_K: nahe L exakte lokale Formel (Zweig log(z − L) = Log(−L) + Log(1 − z/L)), sonst Taylor-Reihe um 0."""
    if abs(z - L) < mp.mpf("0.12"):
        t = z - L
        s = pv(psi, t)
        logs = LML + mp.log(1 - z / L) + mp.log(s / t)
        qq = mp.exp(G * logs)
        e = mp.fsum(eps[k] * qq ** k for k in range(1, len(eps)))
        return pv(phi, SL * s * mp.exp(L * e))
    w = qpoly(acb(arb(mp.nstr(mp.re(z), 55)), arb(mp.nstr(mp.im(z), 55))))
    return mp.mpc(mp.mpf(w.real.mid().str(50, radius=False)), mp.mpf(w.imag.mid().str(50, radius=False)))


# Linie für die Modenamplitude (Motor einmal, dann Zwischenspeicher)
lin_pfad = DATEN / "linie48.json"
if lin_pfad.is_file():
    ld = json.loads(lin_pfad.read_text(encoding="utf-8"))
    ws = [mp.mpf(x) for x in ld["w"]]
    zs = [mp.mpf(x) for x in ld["z"]]
    Tp = [mp.mpf(x) for x in ld["Tp"]]
else:
    from fatou_backend import FatouGP
    from flint import acb as _acb
    kd = json.loads((DATEN / "kreis_r1.34_M4096.json").read_text(encoding="utf-8"))
    M = kd["M"]
    Ah = [mp.mpc(*x) for x in kd["A"]]
    Av = Ah + [mp.conj(x) for x in Ah[1:-1][::-1]]
    ctx.prec = 260
    dd = _acb.dft([_acb(arb(mp.nstr(mp.re(x), 58)), arb(mp.nstr(mp.im(x), 58))) for x in Av])
    Ak = [mp.mpf(dd[k].real.mid().str(58, radius=False)) / M / mp.mpf("1.34") ** k for k in range(1200)]
    ctx.prec = 200
    m = 48
    ws = [mp.mpf("-1.6") + mp.mpf(j) / m for j in range(m)]
    gp = FatouGP(dps=57, fatou_gp="fork")
    zs = [mp.re(z) for z in gp.session("exp(1)").sexp_batch(ws)]
    gp.close()
    Tp = [mp.re(pv([k * q[k] for k in range(1, len(q))], z) / pv([k * Ak[k] for k in range(1, len(Ak))], z)) for z in zs]
    lin_pfad.write_text(json.dumps({"w": [mp.nstr(x, 55) for x in ws], "z": [mp.nstr(x, 55) for x in zs],
                                    "Tp": [mp.nstr(x, 55) for x in Tp]}), encoding="utf-8")
fz = [pv(q, z) for z in zs]


def koef(al, N):
    return mp.exp(mp.loggamma(N - al) - mp.loggamma(-al) - mp.loggamma(N + 1)) * L ** (-N)


AS = [0, 1, 2, 3, 3.4]
VS_R = [0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64]
PHIS = [-80, -60, -40, -20, 0, 20, 40, 60, 80]
VS_C = [0.5, 1, 2, 4, 8, 16]
Ldat = lade()
aus_pfad = DATEN / "lin_auswertung.json"
aus = json.loads(aus_pfad.read_text(encoding="utf-8")) if aus_pfad.is_file() else {}
t0 = time.time()
argL = float(mp.arg(L))
Lc = complex(L)
for N in sorted(Ldat):
    if str(N) in aus:
        continue
    b = [mp.mpf(x) for x in Ldat[N]["b"]]
    bp = acb_poly([acb(arb(x)) for x in Ldat[N]["b"]])

    def P(z):
        w = bp(acb(arb(mp.nstr(mp.re(z), 45)), arb(mp.nstr(mp.im(z), 45))))
        return mp.mpc(mp.mpf(w.real.mid().str(45, radius=False)), mp.mpf(w.imag.mid().str(45, radius=False)))

    # (1) Kreis
    c = np.array([float(x * RL ** k) for k, x in enumerate(b)])
    Mg = 1 << int(np.ceil(np.log2(16 * (N + 1))))
    th1 = 2 * np.pi * np.arange(Mg) / Mg
    v1 = np.fft.ifft(c, Mg) * Mg
    th2 = argL + np.linspace(-60.0 / N, 60.0 / N, 6001)
    zz = float(RL) * np.exp(1j * th2)
    v2 = np.polyval(c[::-1], np.exp(1j * th2))
    th = np.concatenate([th1, th2])
    vv = np.concatenate([v1, v2])
    zall = float(RL) * np.exp(1j * th)
    wgt = np.abs(zall - Lc) * np.abs(zall - np.conj(Lc))
    S = {}
    for a in AS:
        y = np.abs(vv) * wgt ** a
        i = int(np.argmax(y))
        tt = (th[i] + np.pi) % (2 * np.pi) - np.pi
        S[str(a)] = [float(y[i]), float((abs(tt) - argL) * N) if a > 0 else float(tt)]
    # (2) P(L) − L
    PL = P(L) - L
    # (3) Amplitude
    th_ = [(P(z) - f) / t for z, f, t in zip(zs, fz, Tp)]
    cN = mp.fsum(t * mp.expjpi(2 * w) for t, w in zip(th_, ws)) / len(ws)
    c0 = mp.fsum(th_) / len(ws)
    KN = cN * Dd * koef(1 - G, N) / (N * D1 * koef(1 + G, N))
    # (4) Eckschicht
    radial = {}
    for v in VS_R:
        z = L * (1 - mp.mpf(v) / N)
        E = P(z) - (L if v == 0 else fK(z))
        radial[str(v)] = [mp.nstr(mp.re(E), 15), mp.nstr(mp.im(E), 15)]
    komplex = {}
    for v in VS_C:
        for ph in PHIS:
            vv_ = mp.mpf(v) * mp.expj(mp.pi * ph / 180)
            z = L * (1 - vv_ / N)
            E = P(z) - fK(z)
            komplex[f"{v}|{ph}"] = [mp.nstr(mp.re(E), 15), mp.nstr(mp.im(E), 15)]
    aus[str(N)] = {"S": S, "PL": [mp.nstr(mp.re(PL), 15), mp.nstr(mp.im(PL), 15)],
                   "c": [mp.nstr(mp.re(cN), 15), mp.nstr(mp.im(cN), 15)], "c0": mp.nstr(abs(c0), 3),
                   "K": [mp.nstr(mp.re(KN), 12), mp.nstr(mp.im(KN), 12)], "radial": radial, "komplex": komplex}
    aus_pfad.write_text(json.dumps(aus, indent=0), encoding="utf-8")
    print(f"N = {N}: S_a = " + ", ".join(f"{a}: {S[str(a)][0]:.4g} (t*N={S[str(a)][1]:.3g})" for a in AS)
          + f";  |P(L)−L| = {mp.nstr(abs(PL), 4)};  |c| = {mp.nstr(abs(cN), 4)}, K = {mp.nstr(KN, 5)}, |c0|/|c| = "
          f"{mp.nstr(abs(c0) / abs(cN), 2)}  ({time.time() - t0:.0f} s)", flush=True)
