"""Strang V — Test T1 (Auftrag Leitung 08:58): konvergiert Π_N U → U auf ∂D nahe L, wenn U in L singulär ist und die Spirale schnell ist?

Modell: kubisch f(x) = x + (x² − 2Re L·x + ρ²)(α + βx), reelle Koeffizienten, f(L) = L, f'(L) = μ (frei, |μ| und arg μ), L = Knesers Fixpunkt.
Exakte Lösung U = g^s, g(x) = (1 − x/L)(1 − x/L̄) (Schlitze in L, L̄), Daten W = J_∞U = U∘f + (f'∘f)·U als Taylor-Reihe (G = g∘f Polynom, G^s per
Miller-Rekursion). Π_N U = J_N⁻¹T_N W (dieselbe Abschnittsaufgabe wie für f_K). Ausgabe |Π_N U − U| an Innenpunkten und auf ∂D bei τ = arg z − arg L
(τ < 0: Γout-Seite), dazu die Lage der übrigen Urbilder von L unter f (müssen ausserhalb D̄ liegen) und f(0).
Aufruf: python v45_t1_modell.py |μ| argμ_grad s N1,N2,... [fac]
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
MU_ABS = sys.argv[1]
MU_ARG = sys.argv[2]
S = sys.argv[3]
NS = [int(x) for x in sys.argv[4].split(",")]
FAC = float(sys.argv[5]) if len(sys.argv) > 5 else 1.5
TAUS = [-0.30, -0.20, -0.12, -0.08, -0.05, -0.03, -0.015, -0.007, 0.007, 0.015, 0.03, 0.05, 0.08, 0.12, 0.20]


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


def miller(G, s, N):
    """Koeffizienten von (G/G_0)^s bis Grad N (G Polynom, Liste)."""
    g = [x / G[0] for x in G]
    d = len(g) - 1
    P = [arb(1)]
    for n in range(1, N + 1):
        acc = arb(0)
        for j in range(1, min(n, d) + 1):
            acc += ((s + 1) * j - n) * g[j] * P[n - j]
        P.append(acc / n)
    return P


def main():
    t0 = time.time()
    ctx.prec = 200
    c, L = koeff(200)
    cf = [complex(float(x.mid())) for x in c]
    Lc = complex(float(L.real.mid()), float(L.imag.mid()))
    rho = abs(Lc)
    # Urbilder von L unter f: f(x) − L = 0 (kubisch)
    pc = list(cf)
    pc[0] = pc[0] - Lc
    wurz = np.roots(pc[::-1])
    fdach = lambda x: np.polyval(np.array(cf[::-1]), x)  # noqa: E731
    # max|f| auf dem Kreis, Bild der Kreishälften
    th = cmath.phase(Lc)
    mx = max(abs(fdach(rho * cmath.exp(1j * t))) for t in np.linspace(0, 2 * math.pi, 2000))
    gin = [abs(fdach(rho * cmath.exp(1j * t))) / rho for t in np.linspace(th + 0.01, math.pi, 400)]
    print(f"[{'quintisch ' + __import__('os').environ.get('V_QUINT', '') if __import__('os').environ.get('V_QUINT') else 'kubisch'}, |μ| = {MU_ABS}, arg μ = {MU_ARG}°, s = {S}] f = {', '.join(f'{x.real:.5f}' for x in cf)}; Urbilder von L: "
          f"{', '.join(f'{w.real:+.3f}{w.imag:+.3f}i (|·|/ρ = {abs(w) / rho:.3f})' for w in wurz)}; max_{{|x|=ρ}}|f|/ρ = {mx / rho:.3f}; "
          f"max |f(Γin)|/ρ = {max(gin):.3f}", flush=True)
    pts = [rho * cmath.exp(1j * (th + t)) for t in TAUS] + [0.9 * Lc, 0.5 * Lc, -0.9 * rho, 0.0]
    namen = [f"τ={t:+.3f}" for t in TAUS] + ["0.9L", "0.5L", "−0.9ρ", "0"]
    for N in NS:
        t1 = time.time()
        bits = int(FAC * N) + 300
        c, L = koeff(bits)
        ctx.prec = bits
        s = arb(S)
        rho_a = (L.real ** 2 + L.imag ** 2).sqrt()
        f = arb_poly(c)
        g = arb_poly([arb(1), -2 * L.real / rho_a ** 2, 1 / rho_a ** 2])
        G = g(f)
        Gl = [G[k] for k in range(G.degree() + 1)]
        U = miller([g[0], g[1], g[2]], s, N)
        UF0 = Gl[0] ** s                      # G(0) = g(f(0)) > 0 (f(0) reell)
        UF = [UF0 * x for x in miller(Gl, s, N)]
        fp = f.derivative()
        m = fp(f)
        Up = arb_poly(U)
        mU = (m * Up)
        W = [UF[k] + (mU[k] if k <= mU.degree() else arb(0)) for k in range(N + 1)]
        # J̃ = D J D⁻¹, D = diag ρ^k
        pw = [arb_poly([1])]
        for _ in range(N):
            pw.append((pw[-1] * f).truncate(N + 1))
        mk = [m[k] if k <= m.degree() else arb(0) for k in range(N + 1)]
        rp = [rho_a ** k for k in range(N + 1)]
        J = arb_mat(N + 1, N + 1)
        for j in range(N + 1):
            col = pw[j]
            for k in range(N + 1):
                ck = col[k] if k <= col.degree() else arb(0)
                mm = mk[k - j] if k >= j else arb(0)
                J[k, j] = (ck + mm) * rp[k] / rp[j]
        Wt = arb_mat(N + 1, 1)
        for k in range(N + 1):
            Wt[k, 0] = W[k] * rp[k]
        x = J.solve(Wt, algorithm="approx")
        ut = [x[k, 0] for k in range(N + 1)]
        werte = []
        for z in pts:
            q = z / rho
            sv = 0j
            for k in range(N, -1, -1):
                sv = sv * q + complex(float(ut[k].mid()))
            # exakte U(z) = g(z)^s mit Hauptzweig (z in D̄)
            gz = (1 - z / Lc) * (1 - z / Lc.conjugate())
            Uz = cmath.exp(float(S) * cmath.log(gz))
            werte.append(abs(sv - Uz))
        dl2 = math.sqrt(sum((float(ut[k].mid()) - float((U[k] * rp[k]).mid())) ** 2 for k in range(N + 1)))
        # feines τ-Gitter: Bandmaxima auf ∂D (Γout τ < 0, Γin τ > 0) und Wert in L
        uf = np.array([float(ut[k].mid()) for k in range(N + 1)])
        np.save(f"v45k_{MU_ABS}_{MU_ARG}_s{S}_N{N}.npy", uf)
        tg = np.concatenate([-np.geomspace(0.6, 1e-4, 1500), np.geomspace(1e-4, 0.6, 1500)])
        zg = np.exp(1j * (th + tg))
        pv = np.polyval(uf[::-1], zg)
        zz = rho * zg
        gz = (1 - zz / Lc) * (1 - zz / Lc.conjugate())
        Ug = np.exp(float(S) * np.log(gz))
        eg = np.abs(pv - Ug)
        baender = [(-0.6, -0.1), (-0.1, -0.02), (-0.02, -0.004), (-0.004, 0), (0, 0.004), (0.004, 0.02), (0.02, 0.1), (0.1, 0.6)]
        bm = [float(eg[(tg >= a) & (tg <= b)].max()) for a, b in baender]
        eL = abs(np.polyval(uf[::-1], Lc / rho))
        print(f"   Bänder max|Π_N U − U| auf ∂D: " + " ".join(f"[{a:+.3f},{b:+.3f}] {v:.2e}" for (a, b), v in zip(baender, bm)) +
              f" | in L: {eL:.2e} | sup_∂D {eg.max():.2e} bei τ = {tg[eg.argmax()]:+.4f}", flush=True)
        print(f"N = {N} ({bits} Bit, {time.time() - t1:.0f}s): ‖ũ − T_N Ũ‖_ℓ² = {dl2:.3e}; |Π_N U − U|: " +
              " ".join(f"{n} {v:.2e}" for n, v in zip(namen, werte)), flush=True)


if __name__ == "__main__":
    main()
