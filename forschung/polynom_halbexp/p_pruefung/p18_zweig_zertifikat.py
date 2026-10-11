"""Strang P, Prüfung 18: rigoroses Zweig-Zertifikat (Auftrag Leitung, Punkt 3).

Für N = 1 … N_max: Box X_N ∋ p_N (der σ-Homotopie-Zweig von p_1 = x + 1/2). Übergang N → M = N+1 (Krawczyk mit Parameter):
  G_{s,c}(a) := T_M[p_a∘p_a − eˣ] − (1 − s)·c·x^M,   s ∈ [0,1],  c ∈ C := [x^M](p∘p) − 1/M! über p ∈ X_N.
Ist K(Z) := z̃ − Y·G_{S,C}(z̃) + (I − Y·J(Z))(Z − z̃) ⊂ int Z (komponentenweise) und pad(X_N) ⊂ Z, so hat G_{s,c} für jedes (s,c) genau eine
Nullstelle a(s,c) ∈ Z; bei s = 0 ist das pad(p_N) (G_{0,c}(pad p_N) = 0 für das wahre c), a(·,c) ist stetig (Eindeutigkeit + Kompaktheit),
also verläuft der σ-Homotopie-Pfad ganz in Z und p_M = a(1) ∈ Z (J(a) regulär auf Z ⇒ der Pfad ist DIE Fortsetzung). Danach wird die Box für
s = 1 per Krawczyk/Intervall-Newton verkleinert. Arb-Intervallrechnung (python-flint), Präzision 4N + 400 Bit.
Aufruf: python p18_zweig_zertifikat.py NMAX   → p18_zertifikat.json (Breiten, Krawczyk-Reserven, Zeiten)
"""
import json
import sys
import time
from pathlib import Path

from flint import arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60


def F_J(a, N, mit_J=True):
    """F_k = [x^k](p∘p) − 1/k!, J[k][j] = [x^k]p^j + [x^{k−j}](p'∘p); a = Liste arb (auch Intervalle)."""
    f = arb_poly(a)
    P = [arb_poly([1])]
    for _ in range(N):
        P.append((P[-1] * f).truncate(N + 1))
    g, fpf = arb_poly([0]), arb_poly([0])
    for m in range(N + 1):
        g += a[m] * P[m]
        if m and mit_J:
            fpf += m * a[m] * P[m - 1]
    F = arb_mat(N + 1, 1)
    fak = arb(1)
    for k in range(N + 1):
        fak = fak * k if k else fak
        F[k, 0] = g[k] - 1 / fak
    if not mit_J:
        return F, None
    J = arb_mat(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            J[k, j] = P[j][k] + (fpf[k - j] if j <= k else 0)
    return F, J


def mitte_mat(M):
    return arb_mat([[arb(M[i, j].mid()) for j in range(M.ncols())] for i in range(M.nrows())])


def newton_mid(a, N, schritte=6):
    """Mittelpunkts-Newton (ohne Radien) für das System N."""
    a = [arb(x.mid()) for x in a]
    for _ in range(schritte):
        F, J = F_J(a, N)
        d = mitte_mat(J).solve(mitte_mat(F))
        a = [arb((a[i] - d[i, 0]).mid()) for i in range(N + 1)]
    return a


def krawczyk(Z, zc, G_zc, Y, N):
    _, JZ = F_J(Z, N)
    eins = arb_mat([[1 if i == j else 0 for j in range(N + 1)] for i in range(N + 1)])
    d = arb_mat([[Z[k] - zc[k]] for k in range(N + 1)])
    K = arb_mat([[zc[k]] for k in range(N + 1)]) - Y * G_zc + (eins - Y * JZ) * d
    return [K[k, 0] for k in range(N + 1)]


def enthalten(K, Z):
    """K ⊂ int Z komponentenweise; Rückgabe (ok, minimale relative Reserve)."""
    ok, res = True, 1e300
    for kk, zz in zip(K, Z):
        lo_ok = (kk.mid() - kk.rad()) > (zz.mid() - zz.rad())
        hi_ok = (kk.mid() + kk.rad()) < (zz.mid() + zz.rad())
        if not (lo_ok and hi_ok):
            ok = False
        r = float(zz.rad())
        if r > 0:
            res = min(res, (r - float(abs(kk.mid() - zz.mid()) + kk.rad())) / r)
    return ok, res


def newton_s(a, N, s, C_mid, schritte=8):
    """Mittelpunkts-Newton für G_s(a) = F(a) − (1−s)c x^N."""
    a = [arb(x.mid()) for x in a]
    for _ in range(schritte):
        F, J = F_J(a, N)
        F[N, 0] = F[N, 0] - (1 - s) * C_mid
        d = mitte_mat(J).solve(mitte_mat(F))
        a = [arb((a[i] - d[i, 0]).mid()) for i in range(N + 1)]
    return a


def segment(Xa, sa, sb, C, M):
    """Krawczyk für s ∈ [sa, sb], Startbox Xa (enthält a(sa)). Rückgabe (ok, Box bei sb, Reserve)."""
    Cm = arb(C.mid())
    am = [arb(x.mid()) for x in Xa]
    ab = newton_s(am, M, sb, Cm)
    zc = [arb(((am[k] + ab[k]) / 2).mid()) for k in range(M + 1)]
    Z = []
    RHO = arb("1.3745570107437075")
    hmax = max(abs(ab[k] - am[k]) * RHO ** k for k in range(M + 1)) / 2
    xmax = max(Xa[k].rad() * RHO ** k for k in range(M + 1))
    for k in range(M + 1):
        half = abs(ab[k] - am[k]) / 2
        r = (half * arb("1.6") + Xa[k].rad() * 2 + (hmax * arb("0.25") + xmax * 2) * RHO ** (-k)
             + arb(2) ** (-(ctx.prec // 2)) + abs(ab[k]) * arb(2) ** (-(ctx.prec // 3)))
        Z.append(arb(zc[k].mid(), r.upper()))
    S = arb((sa + sb) / 2, (sb - sa) / 2)
    Fz, Jz = F_J(zc, M)
    Y = mitte_mat(Jz).inv()
    Gz = arb_mat([[Fz[k, 0] - ((1 - S) * C if k == M else 0)] for k in range(M + 1)])
    K = krawczyk(Z, zc, Gz, Y, M)
    ok, res = enthalten(K, Z)
    a_in = all((p.mid() - p.rad()) > (z.mid() - z.rad()) and (p.mid() + p.rad()) < (z.mid() + z.rad()) for p, z in zip(Xa, Z))
    if not (ok and a_in):
        return False, None, res
    # verkleinern bei s = sb
    Fa, Ja = F_J(ab, M)
    Fa[M, 0] = Fa[M, 0] - (1 - arb(sb)) * C
    Y1 = mitte_mat(Ja).inv()
    Xn = Z
    for _ in range(4):
        K1 = krawczyk(Xn, ab, Fa, Y1, M)
        Xn = [k1 if float(k1.rad()) < float(z.rad()) else z for k1, z in zip(K1, Xn)]
    return True, Xn, res


protokoll = []
ctx.prec = 600
X = [arb("0.5"), arb(1)]                              # p_1 exakt
t0 = time.time()
for N in range(1, NMAX):
    M = N + 1
    ctx.prec = 4 * M + 400
    padX = X + [arb(0)]
    Fp, _ = F_J(padX, M, mit_J=False)
    C = Fp[M, 0]                                       # c-Intervall über X_N
    stapel = [(0.0, 1.0)]
    Xa, sa, nseg, resmin, ok_all = padX, 0.0, 0, 1e9, True
    while sa < 1.0:
        sb = 1.0
        for _ in range(30):
            ok, Xb, res = segment(Xa, sa, sb, C, M)
            if ok:
                break
            sb = sa + (sb - sa) / 2
        if not ok:
            ok_all = False
            break
        Xa, sa = Xb, sb
        nseg += 1
        resmin = min(resmin, res)
    if not ok_all:
        print(f"N = {N} → {M}: Segment ab s = {sa:.3g} nicht zertifizierbar — Abbruch", flush=True)
        protokoll.append({"N": N, "ok": False})
        break
    X = Xa
    breite = max(float(x.rad()) * 1.3745570107 ** k for k, x in enumerate(X))
    protokoll.append({"N": N, "ok": True, "segmente": nseg, "reserve": resmin, "breite_rho": breite, "t": round(time.time() - t0, 1)})
    (HIER / "p18_zertifikat.json").write_text(json.dumps({"protokoll": protokoll, "letztes_N": M,
        "box": [[str(x.mid()), str(x.rad())] for x in X]}), encoding="utf-8")
    if M % 5 == 0 or M < 6 or M >= 38:
        print(f"N = {N} → {M}: Pfad zertifiziert ({nseg} Segmente, min. Reserve {resmin:.3f}), Box-Radius·ρ^k max {breite:.2e}, "
              f"{time.time() - t0:.0f} s", flush=True)
(HIER / "p18_zertifikat.json").write_text(json.dumps({"protokoll": protokoll,
                                                       "letztes_N": protokoll[-1]["N"] + (1 if protokoll[-1]["ok"] else 0),
                                                       "box": [[str(x.mid()), str(x.rad())] for x in X]}), encoding="utf-8")
print("zertifiziert bis N =", protokoll[-1]["N"] + (1 if protokoll[-1]["ok"] else 0))
