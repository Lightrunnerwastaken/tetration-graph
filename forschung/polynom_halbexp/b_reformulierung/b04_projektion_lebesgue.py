"""b04: Die J-Projektion Π_N u := J_N^{-1} T_N[u∘p + p'∘p·u] (J_N = Linearisierung an p = p_N).
Π_N ist eine Projektion auf P_N (Π_N P = P) und tötet ker J_∞. Gemessen wird ihre 'Lebesgue-Konstante'
   Λ_N(r, R) = max_{|z|=r} (R/2π) ∫ |Π_N[k_t](z)| dφ,  k_t(z) = 1/(t − z), t = R e^{iφ},
also eine Schranke für ||Π_N u||_{D_r} ≤ Λ_N(r,R) ||u||_{D_R}. Frage: wächst Λ_N polynomial oder exponentiell in N,
und hängt das davon ab, ob R < |L| oder R > |L| ist?"""
import json, sys, time
import mpmath as mp
from flint import arb, acb, arb_mat, acb_mat, arb_poly, acb_poly, ctx
sys.stdout.reconfigure(encoding="utf-8")
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
RAW = json.load(open(G + "polynome_bis_194.json"))
ctx.prec = 400


def setup(N):
    a = [arb(s) for s in RAW[str(N)]]
    f = arb_poly(a)
    P = [arb_poly([1])]
    for _ in range(N):
        P.append((P[-1] * f).truncate(N + 1))
    fpf = arb_poly([0])
    for m in range(1, N + 1):
        fpf += m * a[m] * P[m - 1]
    fpf = fpf.truncate(N + 1)
    J = arb_mat(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            J[k, j] = P[j][k] + (fpf[k - j] if j <= k else 0)
    return a, fpf, J.inv()


def Pi_kt(a, fpf, Ji, N, t):
    """Koeffizienten von Π_N[k_t] (komplex)."""
    t = acb(t)
    c0 = 1 / (t - a[0])                      # 1/(t − a0 − h) = Σ b_n x^n, b_n = c0 Σ_{k≥1} a_k b_{n−k}
    s1 = [c0]
    for n in range(1, N + 1):
        s = acb(0)
        for k in range(1, n + 1):
            s += a[k] * s1[n - k]
        s1.append(c0 * s)
    s2 = acb_poly([t ** (-(j + 1)) for j in range(N + 1)])
    s2 = (acb_poly(fpf) * s2)
    b = [s1[k] + s2[k] for k in range(N + 1)]
    br = arb_mat([[x.real] for x in b]); bi = arb_mat([[x.imag] for x in b])
    xr = Ji * br; xi = Ji * bi
    return [acb(xr[k, 0], xi[k, 0]) for k in range(N + 1)]


def horner(c, z):
    s = acb(0)
    for x in reversed(c):
        s = s * z + x
    return s


M = 96
for N in (20, 40, 60, 80):
    t0 = time.time()
    a, fpf, Ji = setup(N)
    out = []
    for R in (1.1, 1.25, 1.35, 1.45, 1.6, 1.8):
        coefs = [Pi_kt(a, fpf, Ji, N, R * complex(mp.expjpi(2 * mp.mpf(j) / M))) for j in range(M)]
        for r in (0.5, 1.0):
            if r >= R:
                continue
            lam = 0.0
            for jz in range(0, M, 4):
                z = acb(r * complex(mp.expjpi(2 * mp.mpf(jz) / M + mp.mpf(1) / M)))
                s = sum(abs(complex(horner(c, z).mid())) for c in coefs) * R / M
                lam = max(lam, s)
            out.append(f"R={R} r={r}: {lam:.3e}")
    print(f"N={N} ({time.time()-t0:.0f}s):", " | ".join(out), flush=True)
