"""Strang A: Prüfung der lokalen Form δ_{−1}(z) ≈ L^{3/2} e^{−2πiκ₀} t^{1−γ} an L (t = z − L, t^α = (−L)^α(1 − z/L)^α).

δ_{−1}(z) = T'(A(z) + ½)·e^{−2πiA(z)}, T' per zentraler Differenz am Motor (dps 57, h = 1e-18), z = φ(e^{Lu}) nahe L.
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
mp.mp.dps = 70
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
G = 2j * mp.pi / L
SL = mp.exp(L / 2)
d = json.loads((DATEN / "singularitaet_theta.json").read_text(encoding="utf-8"))
cx = lambda p: mp.mpc(mp.mpf(p[0]), mp.mpf(p[1]))  # noqa: E731
phi = [cx(p) for p in d["phi"]]
k0 = cx(d["kappa"]["0"])


def pv(c, x):
    s = mp.mpc(0)
    for a in reversed(c):
        s = s * x + a
    return s


gp = FatouGP(dps=57, fatou_gp="fork")
s = gp.session("exp(1)")
h = mp.mpf(10) ** -18
for v in ("1.5", "2.5", "3.5", "4.5"):
    u = mp.mpc((mp.mpf("-1.8") - mp.re(L) * mp.mpf(v)) / mp.im(L), mp.mpf(v))
    z = pv(phi, mp.exp(L * u))
    A = s.slog_batch([z])[0]
    tp, tm = s.sexp_batch([A + mp.mpf("0.5") + h, A + mp.mpf("0.5") - h])
    Tp = (tp - tm) / (2 * h)
    delta = Tp * mp.exp(-2j * mp.pi * A)
    t = z - L
    lokal = L * SL * mp.exp(-2j * mp.pi * k0) * mp.exp((1 - G) * (mp.log(-L) + mp.log(1 - z / L)))
    print(f"|t| = {mp.nstr(abs(t), 3)}: δ_(−1)(z) / (L^(3/2) e^(−2πiκ₀) t^(1−γ)) = {mp.nstr(delta / lokal, 10)}")
gp.close()
