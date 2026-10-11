"""Strang P, Satz P-V Schritt 1: Startdaten aus dem Kneser-Motor (gp-tetration-Fork, dps 50; nur lesen/benutzen).

(a) T̃₀: sexp auf |z| = 1 (K = 512 Punkte) ⇒ Taylor-Koeffizienten um 0 per DFT.
(b) θ̃₀: θ(z) = χ(sexp(z)) − z auf der Geraden Im z = h₀ (256 Punkte), χ = log_L Ψ = m + log(w)/L (Ψ(y) = L^m w, kc_basis.Psi),
    stetig gemacht modulo γ = 2πi/L; Fourier-Koeffizienten c_k, k ≥ 0 (negative Moden als Kontrolle).
Ausgabe p26_start.json. Motor: höchstens dieser eine Prozess.
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp
from flint import acb, arb, ctx

import kc_basis as K

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
DPS = 50
KT = 512
KTH = 256
H0 = mp.mpf("0.28")


def main():
    from fatou_backend import FatouGP
    mp.mp.dps = DPS + 10
    t0 = time.time()
    motor = FatouGP(dps=DPS, fatou_gp="fork")
    s = motor.session("exp(1)")
    # (a)
    zs = [mp.expjpi(2 * mp.mpf(j) / KT) for j in range(KT)]
    werte = s.sexp_batch(zs)
    print("sexp auf |z|=1:", round(time.time() - t0, 1), "s")
    sys.stdout.flush()
    t = []
    for n in range(KT):
        acc = mp.fsum(werte[j] * mp.expjpi(-2 * mp.mpf(n * j) / KT) for j in range(KT)) / KT
        t.append(acc)
    print("t0..t3:", [mp.nstr(x, 20) for x in t[:4]], " max|Im t_n|:", mp.nstr(max(abs(mp.im(x)) for x in t), 5))
    # (b)
    xs = [mp.mpf(j) / KTH - mp.mpf("0.5") for j in range(KTH)]
    zl = [mp.mpc(x, H0) for x in xs]
    wl = s.sexp_batch(zl)
    print("sexp auf Im z = 0.28:", round(time.time() - t0, 1), "s")
    sys.stdout.flush()
    K.setup(400)
    Lc = K.L()
    gam = 2 * acb.pi() * acb(0, 1) / Lc
    th = []
    for x, y in zip(zl, wl):
        ya = acb(mp.nstr(mp.re(y), DPS + 5), mp.nstr(mp.im(y), DPS + 5))
        w, m = K.Psi(ya)
        chi = m + w.log() / Lc
        th.append(chi - acb(mp.nstr(mp.re(x), 30), mp.nstr(mp.im(x), 30)))
    # stetig machen modulo γ
    for j in range(1, len(th)):
        best = min(range(-3, 4), key=lambda q: float(abs(th[j] + q * gam - th[j - 1]).mid()))
        th[j] = th[j] + best * gam
    sprung = float(abs(th[-1] - th[0]).mid())
    print("θ: max Sprung benachbart", max(float(abs(th[j] - th[j - 1]).mid()) for j in range(1, len(th))), " θ(Ende) − θ(Anfang):", sprung)
    c = []
    twopi = 2 * acb.pi()
    for k in range(-8, KTH // 2):
        acc = acb(0)
        for j, x in enumerate(xs):
            acc += th[j] * (-(twopi * acb(0, 1)) * k * acb(mp.nstr(x, 40))).exp()
        acc = acc / KTH * (twopi * k * arb(mp.nstr(H0, 30))).exp()
        c.append((k, acc))
    neg = max(float(abs(v).mid()) for k, v in c if k < 0)
    print("negative Moden (sollen ~0):", neg)
    for k, v in c[8:14]:
        print("  c_%d = %s" % (k, v.str(25)))
    ck = [v for k, v in c if k >= 0]
    raus = {"dps": DPS, "h0": "0.28",
            "t": [[mp.nstr(mp.re(x), DPS + 3), mp.nstr(mp.im(x), 8)] for x in t],
            "c": [[v.real.mid().str(DPS + 3, radius=False), v.imag.mid().str(DPS + 3, radius=False)] for v in ck],
            "neg_moden": neg}
    (HIER / "p26_start.json").write_text(json.dumps(raus, indent=0), encoding="utf-8")
    print("fertig", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
