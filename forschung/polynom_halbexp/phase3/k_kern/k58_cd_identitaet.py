"""k58 — Christoffel-Darboux-Identität (CD1) für J̃_n prüfen und am Steg zerlegen (BERICHT §40). Kubisch.
(CD1)  (z − t)·K_n(z,t) = ρ·q̃(z)·Λ_n(t) − t·ν(z) + Ξ(z,t)
  K_n(z,t) = ε_zᵀJ̃_n⁻¹ẽ_t;  Λ_n(t) = e_nᵀJ̃_n⁻¹ẽ_t (adjungierter Kopf, roh);  ν = J_n⁻¹·1;
  q̃(z) = (z/ρ)^{n+1} − σ̃(z),  J̃σ̃ = T_n f^{n+1} (skaliert, /ρ^{n+1})  — monischer Kopf der Ordnung n+1 (T_n J q = 0);
  Ξ(z,t) = J_n⁻¹T_n[(f − x)·r_t∘f](z),  r_t = K_n(·,t);  r_t∘f bis Grad n = ẽ_t − M̃r̃ (aus J̃r̃ = ẽ_t).
Diagonale (Ableitung von (CD1) nach z in z = t):  K(t,t) = ρq̃′(t)Λ(t) − tν′(t) + ∂_zΞ(t,t);  Nebenbedingung 0 = ρq̃(t)Λ(t) − tν(t) + Ξ(t,t).
Aufruf: python k58_cd_identitaet.py 150 300
"""
import math
import os
import sys
import time

import mpmath as mp

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_mat, acb_mat, arb_poly, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM, kubik_koeff, modell  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
PUNKTE = [(th, s) for th in (100.0, 120.0, 150.0, 180.0) for s in (0.0, 0.03)]


def mpc_of(x, nd):
    return mp.mpc(mp.mpf(x.real.mid().str(nd, radius=False)), mp.mpf(x.imag.mid().str(nd, radius=False)))


def horner(co, w):
    s_ = mp.mpc(0)
    for c in reversed(co):
        s_ = s_ * w + c
    return s_


def horner_d(co, w):
    d = mp.mpc(0); dd = mp.mpc(0)
    for c in reversed(co):
        dd = dd * w + d
        d = d * w + c
    return d, dd


def lauf(n):
    t0 = time.time()
    bits = int(1.3 * n) + 400
    ctx.prec = bits
    mp.mp.dps = int(bits / 3.32) + 10
    nd = mp.mp.dps
    L = mp.mpc(L_RE, L_IM); rho = abs(L)
    Jr = J_tilde("kubisch", n, bits=bits); ctx.prec = bits
    fk, mk = modell("kubisch", n)
    rho_a = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
    # rechte Seiten: e_0 (ν), b_σ (σ̃)
    f = arb_poly(fk[:4])
    pw = arb_poly([1])
    for _ in range(n + 1):
        pw = arb_poly([pw[k] if k <= pw.degree() else arb(0) for k in range(n + 1)]) * f
    B = arb_mat(n + 1, 2)
    B[0, 0] = arb(1)
    for k in range(n + 1):
        ckk = pw[k] if k <= pw.degree() else arb(0)
        B[k, 1] = ckk * rho_a ** k / rho_a ** (n + 1)
    X = Jr.solve(B, algorithm="approx")
    nu = [mp.mpf(X[k, 0].mid().str(nd, radius=False)) for k in range(n + 1)]
    sig = [mp.mpf(X[k, 1].mid().str(nd, radius=False)) for k in range(n + 1)]
    qt = [-x for x in sig] + [mp.mpf(1)]                       # q̃ in w = z/ρ, Grad n+1
    en = arb_mat(n + 1, 1); en[n, 0] = arb(1)
    Y = Jr.transpose().solve(en, algorithm="approx")
    lam = [mp.mpf(Y[k, 0].mid().str(nd, radius=False)) for k in range(n + 1)]
    Yc = Jr.solve(en, algorithm="approx")
    cn = [mp.mpf(Yc[k, 0].mid().str(nd, radius=False)) for k in range(n + 1)]
    mt = [mp.mpf(mk[k].mid().str(nd, radius=False)) * rho ** k for k in range(n + 1)]   # m skaliert
    ff = [mp.mpf(fk[k].mid().str(nd, radius=False)) for k in range(4)]
    phit = [(ff[i] - (1 if i == 1 else 0)) * rho ** i for i in range(4)]                # (f − x) skaliert
    ck = [mp.mpf(x.mid().str(nd, radius=False)) for x in kubik_koeff()]

    def fz(x): return ck[0] + x * (ck[1] + x * (ck[2] + x * ck[3]))

    def fpz(x): return ck[1] + x * (2 * ck[2] + x * 3 * ck[3])

    # K(·,t) für alle t (eine komplexe Lösung mit mehreren Spalten)
    Jc = acb_mat(Jr)
    E = acb_mat(n + 1, len(PUNKTE))
    ts = []
    for j, (th, s) in enumerate(PUNKTE):
        t = rho * (1 + s) * mp.expjpi(th / 180)
        ts.append(t)
        for k in range(n + 1):
            v = (rho / t) ** k
            E[k, j] = acb(arb(mp.nstr(mp.re(v), nd)), arb(mp.nstr(mp.im(v), nd)))
    R = Jc.solve(E, algorithm="approx")
    rs = [[mpc_of(R[k, j], nd) for k in range(n + 1)] for j in range(len(PUNKTE))]
    # Ξ-Daten: T_n[(f − x)·(ẽ_t − M̃ r̃)]
    G = acb_mat(n + 1, len(PUNKTE))
    for j, t in enumerate(ts):
        et = [(rho / t) ** k for k in range(n + 1)]
        comp = [et[k] - mp.fsum(mt[k - i] * rs[j][i] for i in range(k + 1)) for k in range(n + 1)]
        for k in range(n + 1):
            g = mp.fsum(phit[i] * comp[k - i] for i in range(min(3, k) + 1))
            G[k, j] = acb(arb(mp.nstr(mp.re(g), nd)), arb(mp.nstr(mp.im(g), nd)))
    XI = Jc.solve(G, algorithm="approx")
    xis = [[mpc_of(XI[k, j], nd) for k in range(n + 1)] for j in range(len(PUNKTE))]
    print("n = %d (%.0f s)" % (n, time.time() - t0))
    print("  Prüfung (CD1) abseits der Diagonale (z = t·e^{i·5°}) und Nebenbedingung auf der Diagonale (relativ):")
    zeilen = []
    for j, (th, s) in enumerate(PUNKTE):
        t = ts[j]
        Lam = horner(lam, rho / t)
        z = t * mp.expjpi(5.0 / 180)
        lhs = (z - t) * horner(rs[j], z / rho)
        rhs = rho * horner(qt, z / rho) * Lam - t * horner(nu, z / rho) + horner(xis[j], z / rho)
        con = rho * horner(qt, t / rho) * Lam - t * horner(nu, t / rho) + horner(xis[j], t / rho)
        # Diagonale
        Ktt = horner(rs[j], t / rho)
        q0, q1 = horner_d(qt, t / rho)
        n0, n1 = horner_d(nu, t / rho)
        x0, x1 = horner_d(xis[j], t / rho)
        d0, d1 = horner_d(cn, t / rho)
        w = t / rho
        A1 = rho * q1 / rho * Lam            # ρ·q̃′(t) mit d/dz = (1/ρ) d/dw
        A2 = -t * n1 / rho
        A3 = x1 / rho
        m = fpz(fz(t))
        zeilen.append((th, s, m, Ktt, A1, A2, A3, rho * q0 * Lam, t * n0, x0, w * q1 / q0, w * d1 / d0))
        print("    θ = %5.0f°, s = %.2f: |LHS − RHS|/|LHS| = %s, |Nebenbed.|/|tν| = %s" % (th, s, mp.nstr(abs(lhs - rhs) / abs(lhs), 3), mp.nstr(abs(con) / abs(t * n0), 3)))
        Bt = fpz(fz(t)) * rho * horner(qt, t / rho) * Lam / t
        for dphi in (5.0, 20.0):
            z2 = t * mp.expjpi(dphi / 180)
            reg = horner(xis[j], z2 / rho) + (1 - 1 / Bt) * rho * horner(qt, z2 / rho) * Lam
            print("        z = t·e^{i%.0f°}: Ξ = %s,  ρq̃(z)Λ(t) = %s,  Ξ_reg = Ξ + (1 − 1/B)ρq̃Λ = %s" % (dphi, mp.nstr(horner(xis[j], z2 / rho), 6), mp.nstr(rho * horner(qt, z2 / rho) * Lam, 6), mp.nstr(reg, 6)))
    print("  Diagonale, alles mal m(t):  m·K(t,t) | m·ρq̃′Λ | −m·tν′ | m·∂Ξ | Summe−mK | m·ρq̃Λ/t | m·ν | m·Ξ/t | t(log q̃)′ | t(log δ_n)′")
    for (th, s, m, Ktt, A1, A2, A3, B1, B2, B3, lq, ld) in zeilen:
        t = rho * (1 + s) * mp.expjpi(th / 180)
        f6 = lambda v: mp.nstr(v, 6)
        print("    θ = %5.0f°, s = %.2f: %s | %s | %s | %s | %s | %s | %s | %s | %s | %s" % (
            th, s, f6(m * Ktt), f6(m * A1), f6(m * A2), f6(m * A3), mp.nstr(abs(m * (A1 + A2 + A3 - Ktt)), 2),
            f6(m * B1 / t), f6(m * B2 / t), f6(m * B3 / t), f6(lq), f6(ld)))


def main():
    for a in sys.argv[1:]:
        lauf(int(a))


if __name__ == "__main__":
    main()
