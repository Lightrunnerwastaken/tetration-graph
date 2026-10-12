"""Strang V — T1, Teil 2: primaler Kopf δ_n = J_n⁻¹e_n längs ∂D nahe L in Modellen mit schnellerer Spirale (Ks Glied (a)/(c) der Warnung §22.6).

Kubisches Modell wie v45 (f(L) = L, f'(L) = μ frei). Kopf c_n = J̃_n⁻¹e_n, δ_n(x) = Σ c_k (x/ρ)^k. Ausgabe je n:
  log10 max|δ_n(x)/δ_n(L)| auf ∂D in Bändern A = |x − L|/|d_n| (d_n = ρ n^{−a}, a = 1 + i log μ/2π) für Γout (τ < 0) und Γin (τ > 0),
  sowie log10(‖c_n‖_{ℓ²}/|δ_n(L)|) (Kreisnorm relativ zum Eckwert; wächst sie super-polynomial, ist die primale Norm nicht-lokal).
Aufruf: python v47_t1_kopf.py |μ| argμ_grad n1,n2,... [fac]
"""
import cmath
import math
import sys
import time

import numpy as np
from flint import arb, acb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
L_RE = "0.31813150520476413531265425158766451720351761387139986692237"
L_IM = "1.33723570143068940890116214319371061253950213846051241687"
MU_ABS, MU_ARG = sys.argv[1], sys.argv[2]
NS = [int(x) for x in sys.argv[3].split(",")]
FAC = float(sys.argv[4]) if len(sys.argv) > 4 else 1.2


def koeff(prec):
    ctx.prec = prec
    L = acb(arb(L_RE), arb(L_IM))
    phi = arb(MU_ARG) * arb.pi() / 180
    mu = arb(MU_ABS) * acb(phi.cos(), phi.sin())
    rL = (mu - 1) / (L - L.conjugate())
    beta = rL.imag / L.imag
    alpha = rL.real - beta * L.real
    R2 = L.real ** 2 + L.imag ** 2
    import os
    if os.environ.get("V_QUINT"):
        # Quintik f = x + q(x)(α + βx + γx² + δx³), q = x² − 2Re L·x + ρ², f'(L) = μ: α + βL = rL − γL² − δL³ (γ, δ aus V_QUINT)
        g_, d_ = [arb(v) for v in os.environ["V_QUINT"].split(",")]
        Rr = rL - g_ * L ** 2 - d_ * L ** 3
        beta = Rr.imag / L.imag
        alpha = Rr.real - beta * L.real
        f = arb_poly([0, 1]) + arb_poly([R2, -2 * L.real, 1]) * arb_poly([alpha, beta, g_, d_])
        return [f[k] for k in range(f.degree() + 1)], L
    return [R2 * alpha, 1 + R2 * beta - 2 * L.real * alpha, alpha - 2 * L.real * beta, beta], L


def main():
    Lc = complex(float(L_RE), float(L_IM))
    rho, th = abs(Lc), cmath.phase(Lc)
    mu = float(MU_ABS) * cmath.exp(1j * math.radians(float(MU_ARG)))
    a = 1 + 1j * cmath.log(mu) / (2 * math.pi)
    print(f"[kubisch, |μ| = {MU_ABS}, arg μ = {MU_ARG}°] a = {a.real:.4f}{a.imag:+.4f}i (Drehung {math.degrees(a.imag):.2f}°/e-Faktor)", flush=True)
    for n in NS:
        t0 = time.time()
        bits = int(FAC * n) + 300
        c, L = koeff(bits)
        ctx.prec = bits
        rho_a = (L.real ** 2 + L.imag ** 2).sqrt()
        f = arb_poly(c)
        m = f.derivative()(f)
        pw = [arb_poly([1])]
        for _ in range(n):
            pw.append((pw[-1] * f).truncate(n + 1))
        mk = [m[k] if k <= m.degree() else arb(0) for k in range(n + 1)]
        rp = [rho_a ** k for k in range(n + 1)]
        J = arb_mat(n + 1, n + 1)
        for j in range(n + 1):
            col = pw[j]
            for k in range(n + 1):
                ck = col[k] if k <= col.degree() else arb(0)
                mm = mk[k - j] if k >= j else arb(0)
                J[k, j] = (ck + mm) * rp[k] / rp[j]
        e = arb_mat(n + 1, 1)
        e[n, 0] = arb(1)
        x = J.solve(e, algorithm="approx")
        cf = np.array([float(x[k, 0].mid()) for k in range(n + 1)])
        dL = np.polyval(cf[::-1], Lc / rho)
        dn = abs(rho * n ** (-a))
        tg = np.concatenate([-np.geomspace(1.5, 1e-5, 3000), np.geomspace(1e-5, 1.5, 3000)])
        zg = np.exp(1j * (th + tg))
        vals = np.abs(np.polyval(cf[::-1], zg)) / abs(dL)
        A = np.abs(rho * zg - Lc) / dn
        norm = np.linalg.norm(cf) / abs(dL)
        zeile = []
        for seite, msk in (("Γout", tg < 0), ("Γin", tg > 0)):
            for lo, hi in ((0, 4), (4, 16), (16, 64), (64, 256), (256, 1e9)):
                sel = msk & (A >= lo) & (A < hi)
                zeile.append(f"{seite} A∈[{lo},{hi if hi < 1e8 else '∞'}): {math.log10(vals[sel].max()) if sel.any() else float('nan'):+.2f}")
        print(f"n = {n} ({bits} Bit, {time.time() - t0:.0f}s): |d_n| = {dn:.2e}, log10‖c_n‖/|δ_n(L)| = {math.log10(norm):+.3f}; log10 max|δ/δ(L)|: " +
              ", ".join(zeile), flush=True)


if __name__ == "__main__":
    main()
