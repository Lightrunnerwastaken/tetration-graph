"""b34 (Runde 3): bilineare Pol-Antwort Ξ_N(ζ,ζ') := J_N⁻¹ T_N[ ½ D²Φ(f)[E_ζ, E_ζ'] ],  E_ζ := Π_N k_ζ − k_ζ,
D²Φ(f)[a,b] = (a'∘f) b + (b'∘f) a + (f''∘f) a b   (Φ(p) = p∘p).  Test von (H*₂): dist(ζ)dist(ζ')·max_{D̄}|Ξ_N| ≤ C N^{κ₂}?
Reihen: f'∘f = g = eˣ/f', f''∘f = g'/f'.  (k_ζ)'∘f = 1/(ζ − f)².  Polynom-Teile über Potenzen von f mod x^{N+1}."""
import sys, time
import numpy as np
import mpmath as mp
from flint import arb, acb, arb_poly, acb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b22_polmoden.py", encoding="utf-8").read().split("mp.mp.dps = 100\nres = {}")[0])
mp.mp.dps = 100
aL = abs(L); a0 = mp.mpf(K[0]); R0 = abs(L - a0)
psi = ((mp.arg(L - a0) - mp.pi / 2) + (mp.arg(L) + mp.pi / 2)) / 2
kreis = [aL * mp.expjpi(2 * mp.mpf(j) / 240) for j in range(240)] + [L, mp.conj(L)]


def dist_G(z):
    return min(abs(z) - aL, abs(z - a0) - R0)


def mul(a, b, n):
    return [mp.fsum(a[i] * b[k - i] for i in range(k + 1)) for k in range(n + 1)]


def ser_inv(a, n):
    r = [1 / a[0]]
    for k in range(1, n + 1):
        r.append(-mp.fsum(a[i] * r[k - i] for i in range(1, k + 1)) / a[0])
    return r


for N in (100, 160, 220):
    t0 = time.time()
    pw, g, Ji = baue(N, J=N)
    Kc = [mp.mpf(c) for c in K[:N + 3]]
    gk = [mp.mpf(g[k].mid().str(100, radius=False)) for k in range(N + 1)]
    fp = [(k + 1) * Kc[k + 1] for k in range(N + 1)]
    gp = [(k + 1) * gk[k + 1] for k in range(N)] + [mp.mpf(0)]
    fpp_f = mul(gp, ser_inv(fp, N), N)                       # f''∘f = g'/f'
    PW = [[mp.mpf(pw[j][k].mid().str(100, radius=False)) for k in range(N + 1)] for j in range(N + 1)]   # f^j mod x^{N+1}
    E = {}
    def bereite(zeta):
        P = loese(Ji, rhs_pol(N, pw, g, zeta))             # Π_N k_ζ (Koeffizienten)
        Ek = [P[k] - zeta ** (-(k + 1)) for k in range(N + 1)]            # Taylor von E_ζ bei 0 (bis N)
        dP = [(k + 1) * P[k + 1] for k in range(N)] + [mp.mpf(0)]
        dP_f = [mp.fsum(dP[j] * PW[j][k] for j in range(N + 1)) for k in range(N + 1)]   # P'∘f
        inv = ser_inv([zeta - Kc[0]] + [-c for c in Kc[1:N + 1]], N)       # 1/(ζ − f)
        dk_f = mul(inv, inv, N)                                             # k_ζ'∘f = 1/(ζ − f)²
        dE_f = [dP_f[k] - dk_f[k] for k in range(N + 1)]
        return P, Ek, dE_f
    orte = {"L,t=0.5": L + mp.mpf(0.5) / N * mp.expj(psi), "L,t=2": L + mp.mpf(2) / N * mp.expj(psi),
            "L,t=8": L + mp.mpf(8) / N * mp.expj(psi), "L̄,t=2": mp.conj(L + mp.mpf(2) / N * mp.expj(psi)),
            "links,t=2": -(aL + mp.mpf(2) / N)}
    for k_, z_ in orte.items():
        E[k_] = bereite(z_)
    zeile = []
    for (n1, n2) in (("L,t=0.5", "L,t=0.5"), ("L,t=2", "L,t=2"), ("L,t=8", "L,t=8"), ("L,t=2", "L̄,t=2"), ("L,t=2", "links,t=2"), ("links,t=2", "links,t=2")):
        _, Ea, dEa = E[n1]; _, Eb, dEb = E[n2]
        rhs = [(mul(dEa, Eb, N)[k] + mul(dEb, Ea, N)[k] + mul(fpp_f, mul(Ea, Eb, N), N)[k]) / 2 for k in range(N + 1)]
        Xi = loese(Ji, rhs)
        m = max(abs(horner(Xi, z)) for z in kreis)
        d = dist_G(orte[n1]) * dist_G(orte[n2])
        zeile.append(f"{n1}×{n2}: {mp.nstr(d * m, 3)}")
    print(f"N={N} ({time.time()-t0:.0f}s): dist·dist'·max|Ξ| — " + " | ".join(zeile), flush=True)
