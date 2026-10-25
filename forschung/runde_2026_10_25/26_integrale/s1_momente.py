"""Gebiet 26 (Integrale): senkrechte Momente, Raabe-Binet-Form für jede Gerade, Fourier-Laplace-Formel.

Bausteine nur lesend aus forschung/identitaeten_2026_10_14/etappe3 (lib3: Schröder/Abel an L, θ-Reihe aus etappe2;
Taylor-Daten nahe der Achse). Kein Motor, kein Schreiben ausserhalb dieses Ordners.

Geprüft wird (Basis e, s = sexp(x0+it), D(t) = e^s − s, K(x) = ∫_0^∞ (sexp(x+it) − L) dt, Λ(x) = ∫_x^{x+1} ln sexp'):
  RB   i K(x) + Λ(x) − L (x + 3/2) = ln κ   für jedes x   (V3 auf jeder Geraden)
  M1   ∫_ℝ t D(t) dt = 2i (arg κ + (x0 + 2) Im L)
  M1r  ∫_{x0}^{x0+1} x T(x) dx = Re W1 + ReL (x0+1/2) + x0 (R(x0) − ReL) + Im K(x0+1),  W1 = ∫_0^∞ t D dt
  V    ∫_ℝ (sexp(x0+it) − L_{sgn t}) dt = 2 arg κ + (2 x0 + 3) Im L
  FL   ∫_{x0}^{x0+1} e^{2πimx} T(x) dx = e^{2πimx0} (δ_{m0} L − i ∫_0^∞ e^{−2πmt} D(t) dt),  m ≥ 0
  M2   zweites Moment: braucht die neue Konstante w0 = Im J(0) = ∫_0^∞ t Im(sexp(it) − L) dt
Aufruf: python s1_momente.py [DPS]
"""
import json
import sys
from functools import lru_cache
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
E3 = HIER.parents[1] / "identitaeten_2026_10_14" / "etappe3"
sys.path.insert(0, str(E3))
from lib3 import Basis, ThetaReihe, lade_theta  # noqa: E402

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
mp.mp.dps = dps
TAY = json.loads((E3 / "daten" / "taylor_130_256.json").read_text(encoding="utf-8"))
ZEN = {mp.mpf(k): [mp.mpf(x) for x in v] for k, v in TAY["zentren"].items()}
DZEN = {c: [k * a for k, a in enumerate(v)][1:] for c, v in ZEN.items()}
B = Basis(70)
L = B.L
TR = ThetaReihe(lade_theta(), B)
# κ aus etappe3 (BERICHT_ETAPPE3.md, 40 Stellen) — wird hier unabhängig aus RB zurückgewonnen
KAPPA_E3 = mp.mpc("1.717640195943158037022581479756274116225", "-2.384443011730485878838953547660210549918")


def T_nah(z):
    c = min(ZEN, key=lambda c: abs(z - c))
    if abs(z - c) <= mp.mpf("0.1"):
        return mp.polyval(ZEN[c][::-1], z - c)
    return None


@lru_cache(maxsize=None)
def T_c(x0s, ts):
    """sexp(x0 + i t), t >= 0."""
    z = mp.mpc(mp.mpf(x0s), mp.mpf(ts))
    if mp.mpf(ts) <= mp.mpf("0.09"):
        w = T_nah(z)
        if w is not None:
            return w
    return TR.sexp(z)


def T(x0, t):
    return T_c(mp.nstr(x0, dps + 5), mp.nstr(t, dps + 5))


def T_reell(x):
    """sexp und ln sexp' auf [−0.5, 1.5] (rechts über sexp(x) = e^{sexp(x−1)})."""
    x = mp.mpf(x)
    if x <= mp.mpf("0.5"):
        c = min(ZEN, key=lambda c: abs(x - c))
        v = mp.polyval(ZEN[c][::-1], x - c)
        d = mp.polyval(DZEN[c][::-1], x - c)
        return v, mp.log(d)
    v1, l1 = T_reell(x - 1)
    return mp.exp(v1), v1 + l1


KN = [0, mp.mpf("0.09"), mp.mpf("0.5"), 2, 6, 15, 35, 70, 140]


def q0inf(f):
    return mp.quad(f, KN)


def qper(f, x0):
    pts = [x0] + [p for p in (mp.mpf("-0.3"), mp.mpf("-0.1"), mp.mpf("0.1"), mp.mpf("0.3"), mp.mpf("0.5"), mp.mpf("0.7"), mp.mpf("0.9"), mp.mpf("1.1")) if x0 < p < x0 + 1] + [x0 + 1]
    return mp.quad(f, pts)


def zeile(name, ist, soll):
    print(f"{name:34s} ist = {mp.nstr(ist, 30):>70s}\n{'':34s} soll= {mp.nstr(soll, 30):>70s}   |Abw| = {mp.nstr(abs(ist - soll), 3)}", flush=True)


if __name__ == "__main__":
    ImL, ReL = mp.im(L), mp.re(L)
    print(f"dps={dps}  L = {mp.nstr(L, 30)}")
    erg = {}
    lnk = []
    for x0 in (mp.mpf(0), mp.mpf("-0.4"), mp.mpf("0.2")):
        D = lambda t: mp.exp(T(x0, t)) - T(x0, t)
        K = q0inf(lambda t: T(x0, t) - L)
        Kr = q0inf(lambda t: mp.exp(T(x0, t)) - L)            # K(x0+1)
        I0 = q0inf(D)
        W1 = q0inf(lambda t: t * D(t))
        Lam = qper(lambda x: T_reell(x)[1], x0)
        R = qper(lambda x: T_reell(x)[0], x0)
        Rx = qper(lambda x: x * T_reell(x)[0], x0)
        lk = 1j * K + Lam - L * (x0 + mp.mpf(3) / 2)
        lnk.append(lk)
        print(f"\n--- Gerade x0 = {mp.nstr(x0, 3)}")
        zeile("RB: iK+Λ−L(x+3/2)  vs ln κ (E3)", lk, mp.log(KAPPA_E3))
        argk = mp.arg(KAPPA_E3)
        zeile("V1-Hälfte: Re ∫_0^∞ D", mp.re(I0), ImL)
        zeile("V2: Im ∫_0^∞ D vs R − ReL", mp.im(I0), R - ReL)
        zeile("M1: ∫_ℝ t D dt", 2j * mp.im(W1), 2j * (argk + (x0 + 2) * ImL))
        zeile("M1r: ∫ x T dx", Rx, mp.re(W1) + ReL * (x0 + mp.mpf(1) / 2) + x0 * (R - ReL) + mp.im(Kr))
        zeile("V: ∫_ℝ (s − L_sgn) dt", 2 * mp.re(K), 2 * argk + (2 * x0 + 3) * ImL)
        for m in (1, 2):
            cm = qper(lambda x: mp.expjpi(2 * m * x) * T_reell(x)[0], x0)
            lap = q0inf(lambda t: mp.exp(-2 * mp.pi * m * t) * D(t))
            zeile(f"FL m={m}", cm, mp.expjpi(2 * m * x0) * (-1j) * lap)
        erg[mp.nstr(x0, 3)] = dict(K=K, W1=W1, Lam=Lam, R=R, Rx=Rx)
    # M2: zweites Moment auf x0 = 0, neue Konstante w0
    x0 = mp.mpf(0)
    D = lambda t: mp.exp(T(x0, t)) - T(x0, t)
    W2 = q0inf(lambda t: t ** 2 * D(t))
    J0 = q0inf(lambda t: t * (T(x0, t) - L))
    J1 = q0inf(lambda t: t * (mp.exp(T(x0, t)) - L))
    print("\n--- zweites Moment (x0 = 0)")
    K0 = erg["0.0"]["K"]
    zeile("Im J(1) − Im J(0) vs ∫_0^1 Re K", mp.im(J1) - mp.im(J0), mp.re(K0) + ImL / 2)
    Rx2 = qper(lambda x: x ** 2 * T_reell(x)[0], x0)
    W1 = erg["0.0"]["W1"]
    I0 = q0inf(D)
    rhs = Rx2 - L / 3 + 1j * 0 * I0 - 0 * W1 + 1j * erg["0.0"]["K"] * 0  # Hilfsterme unten
    # i W2 = ∫x²(T−L) + i x0² I0 − 2 x0 W1 + i(2x0+1) K(x0+1) − 2 J(x0+1), hier x0 = 0
    Kr = q0inf(lambda t: mp.exp(T(x0, t)) - L)
    zeile("M2: i W2 (Cauchy mit z²)", 1j * W2, Rx2 - L / 3 + 1j * Kr - 2 * J1)
    print(f"w0 = Im J(0) = {mp.nstr(mp.im(J0), 35)}")
    print(f"Re J(0)      = {mp.nstr(mp.re(J0), 35)}")
    print(f"∫_ℝ t² D dt  = {mp.nstr(2 * mp.re(W2), 35)}")
    print(f"ln κ (RB, Mittel der drei Geraden) = {mp.nstr(sum(lnk) / 3, 35)}")
    print(f"Streuung RB = {mp.nstr(max(abs(a - b) for a in lnk for b in lnk), 3)}")
    out = {k: {kk: [mp.nstr(mp.re(v), dps), mp.nstr(mp.im(v), dps)] for kk, v in d.items()} for k, d in erg.items()}
    out["J0"] = [mp.nstr(mp.re(J0), dps), mp.nstr(mp.im(J0), dps)]
    out["W2_0"] = [mp.nstr(mp.re(W2), dps), mp.nstr(mp.im(W2), dps)]
    (HIER / "cache").mkdir(exist_ok=True)
    (HIER / "cache" / f"momente_{dps}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
