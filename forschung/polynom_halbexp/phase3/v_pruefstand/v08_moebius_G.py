"""Strang V — Plan A2: loxodromisches Möbius-Modell f = M⁻¹(μM), M(x) = (x − L)/(x − L̄), μ = √L (Ds d24, Xs Satz 1.5).
Misst wie Ks k01 (kubisch/Kneser): G_N(t) = ‖J̃_N⁻¹ẽ_t‖ auf Γout, in L, auf einer Linse über Γin und (Kontrolle) auf dem Γin-Kreis,
dazu σ₁(J̃_N⁻¹) (Eckrate c). Frage: trägt das Modell die universelle Ecke (G_N(L) ~ √N, Γout konstant) — dann ist es ein Kandidat für (H*_G).
J_N[k][j] = [x^k]f^j + [x^{k−j}]m, m = f'∘f, komplexe Koeffizienten (acb). Präzision ≈ N + 250 Bit wie k01.
Aufruf: python v08_moebius_G.py 100,200,300 [mu_faktor]
"""
import math
import sys
import time

from flint import acb, acb_mat, acb_poly, arb, ctx

sys.stdout.reconfigure(encoding="utf-8")
Ns = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else "100,200").split(",")]
MUF = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
Lc = complex(0.31813150520476413, 1.3372357014306895)
RHO = abs(Lc)
argL = math.atan2(Lc.imag, Lc.real)


def modell(N):
    L = (-acb(-1).lambertw()).conjugate()
    Lb = L.conjugate()
    mu = L.sqrt() * MUF
    a, b, c, d = L - Lb * mu, -(L * Lb) * (1 - mu), 1 - mu, mu * L - Lb
    # f = (a x + b)/(c x + d) als Potenzreihe
    den = acb_poly([d, c])
    inv = [acb(1) / d]
    for k in range(1, N + 2):
        inv.append(-inv[-1] * c / d)
    invp = acb_poly(inv)
    f = (acb_poly([b, a]) * invp)
    f = acb_poly([f[k] for k in range(N + 2)])
    fp = (acb_poly([a * d - b * c]) * invp * invp)
    fp = acb_poly([fp[k] for k in range(N + 2)])
    # m = f'∘f
    m = acb_poly([0])
    pw = acb_poly([1])
    for k in range(N + 1):
        m += fp[k] * pw
        pw = acb_poly([(pw * f)[i] for i in range(N + 1)])
    return L, mu, f, m


def J_tilde(N):
    ctx.prec = N + 250
    L, mu, f, m = modell(N)
    rho = abs(L)
    pw = [acb_poly([1])]
    for _ in range(N):
        q = pw[-1] * f
        pw.append(acb_poly([q[i] for i in range(N + 1)]))
    J = acb_mat(N + 1, N + 1)
    rp = [rho ** k for k in range(N + 1)]
    for j in range(N + 1):
        for k in range(N + 1):
            ck = pw[j][k]
            mm = m[k - j] if k >= j else acb(0)
            J[k, j] = (ck + mm) * rp[k] / rp[j]
    return J, L, mu, f


def f_wert(f, x, N):
    s = acb(0)
    for k in range(N + 1, -1, -1):
        s = s * x + f[k]
    return s


def punkte(f, N):
    # Γout/Γin über |f(x)| ≶ ρ auf dem Kreis (Möbius: f(Kreis) ist ein Kreis durch L, L̄)
    P = []
    for deg in range(-180, 181, 15):
        x = RHO * complex(math.cos(math.radians(deg)), math.sin(math.radians(deg)))
        fx = complex(f_wert(f, acb(x.real, x.imag), N).mid()) if False else None
        P.append(("Kreis %d°" % deg, x))
    P.append(("L", Lc))
    P.append(("L̄", Lc.conjugate()))
    return P


def main():
    for N in Ns:
        t0 = time.time()
        J, L, mu, f = J_tilde(N)
        # Klassifikation der Kreispunkte
        P = []
        for deg in list(range(-170, 181, 10)):
            ph = math.radians(deg)
            x = RHO * complex(math.cos(ph), math.sin(ph))
            # f(x) exakt über Möbius
            Lb = Lc.conjugate()
            M = (x - Lc) / (x - Lb)
            mm = complex(mu.real.mid(), mu.imag.mid()) * M
            fx = (Lc - Lb * mm) / (1 - mm)
            seite = "in" if abs(fx) < RHO else "out"
            P.append(("Kreis %+4d° (%s)" % (deg, seite), x, seite, deg))
        P.append(("L", Lc, "L", None))
        # Linse über Γin: Radius ρ·(1 + s), s = 0.01·(Abstand zu L in rad)^1.5 begrenzt — grob
        lins = []
        for name, x, seite, deg in P:
            if seite == "in":
                dl = abs(math.radians(deg) - argL)
                dl = min(dl, 2 * math.pi - dl)
                s = 0.012 * min(dl, 1.0) ** 1.5
                lins.append(("Linse %+4d°" % deg, x * math.exp(s), "linse", deg))
        P = P + lins
        B = acb_mat(N + 1, len(P))
        rho = abs(L)
        for cidx, (_, t, _, _) in enumerate(P):
            q = acb(rho) / acb(t.real, t.imag)
            p = acb(1)
            for k in range(N + 1):
                B[k, cidx] = p
                p = p * q
        X = J.solve(B, algorithm="approx")
        G = []
        for cidx in range(X.ncols()):
            s = arb(0)
            for k in range(X.nrows()):
                s += abs(X[k, cidx]) ** 2
            G.append(float(s.sqrt().mid()))
        print(f"N = {N} ({time.time() - t0:.0f}s, prec {ctx.prec}), μ = {mu.str(6)}", flush=True)
        for (name, t, seite, deg), g in zip(P, G):
            print(f"  {name:22s} ln G = {math.log(g):9.4f}   ln G − ½ln N = {math.log(g) - 0.5 * math.log(N):8.4f}", flush=True)


if __name__ == "__main__":
    main()
