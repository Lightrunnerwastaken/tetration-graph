"""Strang V — Eingangswerte von Xs Lemma X8′ (A_θ) streng (arb, Vs (K2)-Datei): min_{∂D}|f'|, max_{∂D}|m/f'|, M2 = max_{∂D}|f''∘f|.
Auf ∂D (x = ρz, |z| = 1): f'(x) = F'(z), f''(x) = F''(z)/ρ, m = eˣ/f', m/f' = eˣ/f'², f''∘f = m'/f' = eˣ(f' − f'')/f'³ (aus m·f' = eˣ).
Intervalle über φ ∈ [0, π] (Symmetrie), Bisektion bis zur geforderten Schärfe. A_θ = (5/ρ)/min|f'| + max|m/f'|/(ρc₀) + M2/(c₀N), c₀ = 0.641.
"""
import sys
from pathlib import Path

from flint import acb, arb, ctx

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER / "m_kopie"))
from mk_kneser_arb import KneserArb  # noqa: E402
sys.argv = ["v40", "0.0166"]
import v40_lipschitz_zert as V40  # noqa: E402  (F2_kreis)

sys.stdout.reconfigure(encoding="utf-8")
KN = V40.KN
ctx.prec = 120
pi = arb.pi()
rho = KN.rho


def werte(phi):
    z = acb(phi.cos(), phi.sin())
    f1 = KN.Fp_kreis(phi)
    f2 = V40.F2_kreis(phi) / rho
    ex = (rho * z).exp()
    return abs(f1), abs(ex / f1 ** 2), abs(ex * (f1 - f2) / f1 ** 3)


def zerlege(lo, hi, tiefe=0):
    phi = arb(((lo + hi) / 2).mid(), ((hi - lo) / 2).upper() * (1 + 1e-9))
    a, b, c = werte(phi)
    if tiefe >= 16 or (a.rad() < 1e-4 and b.rad() < 1e-4 and c.rad() < 1e-4):
        return a.lower(), b.upper(), c.upper()
    m = (lo + hi) / 2
    r1, r2 = zerlege(lo, m, tiefe + 1), zerlege(m, hi, tiefe + 1)
    return min(r1[0], r2[0]), max(r1[1], r2[1]), max(r1[2], r2[2])


mn, mmf, M2 = zerlege(arb(0), pi)
c0 = arb("0.641")
for N in (318, 400):
    At = (5 / rho) / arb(mn) + arb(mmf) / (rho * c0) + arb(M2) / (c0 * N)
    print(f"N = {N}: A_θ ≤ {float(At.upper()):.4f}", flush=True)
print(f"min|f'| ≥ {float(mn):.5f} (x21: 0.3456), max|m/f'| ≤ {float(mmf):.5f} (2.1176), M2 = max|f''∘f| ≤ {float(M2):.5f} (0.8183)", flush=True)
