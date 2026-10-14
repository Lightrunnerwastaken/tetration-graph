"""Senkrechte Integrale (Cauchy auf dem Halbstreifen) für Kneser-sexp, Basis e.

V1  ∫_{-∞}^{∞} [Φ(sexp(x0+1+it)) − Φ(sexp(x0+it))] dt = (Φ(L) − Φ(L̄))/i     (Φ(w) = w, w², 1/w)
V2  ∫_{x0}^{x0+1} sexp = Re L + ∫_0^∞ Im[sexp(x0+1+it) − sexp(x0+it)] dt         (Raabe-Form)
V3  L·(a_0 + 1/2) ≡ ∫_0^1 ln sexp'(x) dx − L + i ∫_0^∞ (sexp(it) − L) dt   (mod 2πi)
V4  Re ∫_0^∞ (sexp(x0+it) − L) dt − Im L · x0 hängt nicht von x0 ab.
Werte: |t| ≤ 0.1 aus den Taylor-Reihen (daten/taylor_130_256.json), sonst sexp = F(z + θ(z)) mit θ-Reihe (etappe2).
Aufruf: python s3_senkrecht.py DPS
"""
import json
import sys

import mpmath as mp

from lib3 import HIER, Basis, ThetaReihe, lade_theta

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 50
mp.mp.dps = dps
T = json.loads((HIER / "daten" / "taylor_130_256.json").read_text(encoding="utf-8"))
ZEN = {mp.mpf(k): [mp.mpf(x) for x in v] for k, v in T["zentren"].items()}
B = Basis(70)
L = B.L
TR = ThetaReihe(lade_theta(), B)


def sexp(z):
    """sexp für Im z >= 0, Re z in [−0.5, 0.5] (Taylor nahe der Achse) bzw. beliebig (θ-Reihe)."""
    if mp.im(z) <= mp.mpf("0.09"):
        c = min(ZEN, key=lambda c: abs(z - c))
        if abs(z - c) <= mp.mpf("0.1"):
            return mp.polyval(ZEN[c][::-1], z - c)
    return TR.sexp(z)


def quad0inf(f):
    return mp.quad(f, [0, mp.mpf("0.09"), 0.5, 2, 6, 15, 35, 70, 140])


def V1(x0, Phi, PhiL):
    g = lambda t: Phi(mp.exp(sexp(mp.mpc(x0, t)))) - Phi(sexp(mp.mpc(x0, t)))
    oben = quad0inf(g)
    # unten: Konjugation, sexp(x0 − it) = conj sexp(x0 + it), Φ reell auf ℝ
    gesamt = oben + mp.conj(oben)
    soll = (PhiL - mp.conj(PhiL)) / 1j
    return gesamt, soll


if __name__ == "__main__":
    print("L =", mp.nstr(L, 30), " 2 Im L =", mp.nstr(2 * mp.im(L), 30))
    for x0 in (mp.mpf(0), mp.mpf("-0.4"), mp.mpf("0.2")):
        for name, Phi in (("w", lambda w: w), ("w^2", lambda w: w ** 2), ("1/w", lambda w: 1 / w)):
            g, s = V1(x0, Phi, Phi(L))
            print(f"V1 x0={mp.nstr(x0,2):5s} Φ={name:4s}: Integral = {mp.nstr(g, 35)}   Soll = {mp.nstr(s, 20)}   |Abw| = {mp.nstr(abs(g - s), 3)}", flush=True)
    # V2 mit x0 = 0 und x0 = -0.4
    for x0 in (mp.mpf(0), mp.mpf("-0.4")):
        I_v = quad0inf(lambda t: mp.im(mp.exp(sexp(mp.mpc(x0, t))) - sexp(mp.mpc(x0, t))))
        lhs = mp.quad(lambda x: sexp(mp.mpf(x)) if abs(x) <= 0.5 else mp.exp(sexp(mp.mpf(x) - 1)), [x0, x0 + mp.mpf("0.5"), 0.5, x0 + 1] if x0 < 0 else [0, 0.5, 1])
        print(f"V2 x0={mp.nstr(x0,2)}: ∫ sexp = {mp.nstr(lhs, 35)},  Re L + senkrecht = {mp.nstr(mp.re(L) + I_v, 35)},  |Abw| = {mp.nstr(abs(lhs - mp.re(L) - I_v), 3)}", flush=True)
    # V3
    K0 = quad0inf(lambda t: sexp(mp.mpc(0, t)) - L)

    def lnd(x):
        c = min(ZEN, key=lambda c: abs(x - c))
        co = ZEN[c]
        return mp.log(mp.polyval([k * co[k] for k in range(len(co) - 1, 0, -1)], x - c))
    # ∫_0^1 ln sexp' : auf [1/2,1] ln sexp'(x) = sexp(x−1) + ln sexp'(x−1)
    I_ld = mp.quad(lnd, [0, 0.25, 0.5]) + mp.quad(lambda x: sexp(x - 1) + lnd(x - 1), [0.5, 0.75, 1])
    nm = json.loads((HIER / "daten" / "nullmoden_40.json").read_text(encoding="utf-8"))
    a0 = lade_theta()[0]
    lhs = L * (a0 + mp.mpf(1) / 2)
    rhs = I_ld - L + 1j * K0
    k = (lhs - rhs) / (2j * mp.pi)
    print(f"V3: ∫ln sexp' = {mp.nstr(I_ld, 35)},  K0 = ∫_0^∞(sexp(it)−L)dt = {mp.nstr(K0, 35)}")
    print(f"V3: (L(a0+1/2) − [∫ln sexp' − L + iK0])/(2πi) = {mp.nstr(k, 30)}")
    for x0 in (mp.mpf("-0.4"), mp.mpf("0.2")):
        Kx = quad0inf(lambda t: sexp(mp.mpc(x0, t)) - L)
        print(f"V4 x0={mp.nstr(x0,2)}: Re K(x0) − Im L·x0 − Re K(0) = {mp.nstr(mp.re(Kx) - mp.im(L) * x0 - mp.re(K0), 3)}", flush=True)
