"""zert_kern — Krawczyk-Übergang N → N+1 längs der σ-Homotopie, KOPIE von Ps p24_zweig_kohaerent.py (p_pruefung/, nur gelesen),
angepasst von Strang X für die parallele Verlängerung des Zweig-Zertifikats (Auftrag Leitung, Meilenstein 7):
  * NC (Zahl der Tschebyschow-Knoten) als Parameter (P: fest 7) — ab M ≈ 210 nötig, weil der Wrapping-Term ∝ 1.88^M·r wächst;
  * newton_pN: hochgenaue (nicht strenge) Lösung des N-Systems als Mitte der Startbox B_N;
  * sonst wörtlich Ps Funktionen (F_J, komp, majorante_sup, uebergang, kodiere/dekodiere).
Die Strenge liegt wie bei P allein in uebergang(): für JEDEN Startpunkt in der Box B_N (als Intervall) verläuft der Homotopiepfad in der
Röhre und endet in der Ausgabebox K_{N+1}. Die Kette wird danach mit zert_kette.py geprüft (K_{N+1} ⊂ B_{N+1}, P-Box von p_193 ⊂ B_193).
"""
import time

from flint import arb, arb_mat, arb_poly, ctx

RHO = 1.37455701


def ob(x):
    return x.abs_upper()


def potenzen(f, M):
    P = [arb_poly([1])]
    for _ in range(M):
        P.append((P[-1] * f).truncate(M + 1))
    return P


def F_J(a, M, mit_J=True):
    """F_k = [x^k](p∘p) − 1/k!, J[k][j] = [x^k]p^j + [x^{k−j}](p'∘p)."""
    f = arb_poly(a)
    P = potenzen(f, M)
    g, fpf = arb_poly([0]), arb_poly([0])
    for m in range(M + 1):
        g += a[m] * P[m]
        if m and mit_J:
            fpf += m * a[m] * P[m - 1]
    F, fak = [], arb(1)
    for k in range(M + 1):
        if k:
            fak *= k
        F.append(g[k] - 1 / fak)
    if not mit_J:
        return F, None
    J = arb_mat(M + 1, M + 1)
    for k in range(M + 1):
        for j in range(M + 1):
            J[k, j] = P[j][k] + (fpf[k - j] if j <= k else 0)
    return F, J


def komp(h, P0, M):
    c = h.coeffs()
    r = arb_poly([c[-1]]) if c else arb_poly([0])
    for k in range(len(c) - 2, -1, -1):
        r = (r * P0).truncate(M + 1) + c[k]
    return r


def matvec(A, v):
    n = A.nrows()
    x = A * arb_mat([[t] for t in v])
    return [x[i, 0] for i in range(n)]


def punkt(v):
    return [arb(t.mid()) for t in v]


def norm_rho(v, rho):
    return max(float(ob(t)) * rho ** k for k, t in enumerate(v)) if v else 0.0


def majorante_sup(Ahat, M):
    f = arb_poly(Ahat)
    P = potenzen(f, M)
    comp, fpf = arb_poly([0]), arb_poly([0])
    for m in range(M + 1):
        comp += Ahat[m] * P[m]
        if m:
            fpf += m * Ahat[m] * P[m - 1]
    return P, comp, fpf


def amax(a, b):
    return (a + b + abs(a - b)) / 2


def kodiere(t):
    m, r, e = t.mid_rad_10exp(ctx.dps + 20)
    return [str(m), str(r), int(e)]


def dekodiere(x):
    m, r, e = x
    z = arb(10) ** int(e)
    return arb(int(m)) * z + arb(0, 1) * (arb(int(r)) * z)


def newton_pN(start, N, max_it=12):
    """Nicht strenge Newton-Lösung von T_N[a∘a − eˣ] = 0 (Punktarithmetik, Präzision ctx.prec) ab start (Liste von arb/Strings).
    Gibt (Mitte, Protokoll) zurück; Protokoll: ρ-gewichtete Residuen je Schritt."""
    z = [arb(s) if not isinstance(s, arb) else arb(s.mid()) for s in start]
    z = (z + [arb(0)] * (N + 1))[: N + 1]
    prot = []
    for _ in range(max_it):
        F, J = F_J(z, N)
        res = norm_rho(punkt(F), RHO)
        prot.append(res)
        if res < max(2.0 ** (-(ctx.prec - 60)), 1e-290):
            break
        d = J.mid().solve(arb_mat([[f.mid()] for f in F])).entries()
        z = punkt([z[k] - d[k] for k in range(N + 1)])
        if len(prot) >= 2 and prot[-1] < 1e-100 and prot[-1] > 0.01 * prot[-2]:
            break                        # Rundungsboden erreicht
    F, _ = F_J(z, N, mit_J=False)
    prot.append(norm_rho(punkt(F), RHO))
    return z, prot


def uebergang(box, N, NC):
    """Ps uebergang() mit Knotenzahl NC als Parameter (sonst unverändert)."""
    M = N + 1
    t0 = time.time()
    ctx.prec = 2 * M + 320
    rho = RHO
    w = [arb(rho) ** (-k) for k in range(M + 1)]
    eM = [arb(1) if k == M else arb(0) for k in range(M + 1)]
    P0 = box + [arb(0)]
    F0, _ = F_J(P0, M, mit_J=False)
    C = F0[M]
    cm = arb(C.mid())
    P0m = punkt(P0)
    Fm, J0 = F_J(P0m, M)
    Y0 = J0.mid().inv().mid()
    a1 = [-cm * t for t in matvec(Y0, eM)]
    a1p = punkt(a1)
    d1 = arb_poly(a1p)
    P0p = arb_poly(P0m)
    q = (komp(d1.derivative(), P0p, M) * d1).truncate(M + 1) \
        + (komp(P0p.derivative().derivative(), P0p, M) * d1 * d1).truncate(M + 1) * arb(0.5)
    qv = [q[k] for k in range(M + 1)]
    a2 = punkt([-t for t in matvec(Y0, qv)])
    pi = arb.pi()
    knoten = [(1 - ((2 * j + 1) * pi / (2 * NC)).cos()) / 2 for j in range(NC)]
    knoten = [arb(s.mid()) for s in knoten]
    V, G, Jk = [], [], []
    for s in knoten:
        z = [P0m[k] + s * a1p[k] + s * s * a2[k] for k in range(M + 1)]
        z = punkt(z)
        alt = None
        for it in range(8):
            Fz, _ = F_J(z, M, mit_J=False)
            Gz = [Fz[k] - (1 - s) * cm * eM[k] for k in range(M + 1)]
            n = norm_rho(Gz, rho)
            if alt is not None and n > alt / 50:
                break
            alt = n
            z = punkt([z[k] - t for k, t in enumerate(matvec(Y0, Gz))])
        # X-Änderung: Ps Sehnen-Newton (mit Y0) bleibt bei ~1e-51 stehen; für M ≳ 220 muss der Röhrenradius r ∝ 1.88^{−M} sinken,
        # also volle Newton-Schritte an den Knoten bis an die Rechengenauigkeit (Prädiktor, nicht Teil der Strenge).
        for it2 in range(5):
            Fz, Jn = F_J(z, M)
            Gz = [Fz[k] - (1 - s) * cm * eM[k] for k in range(M + 1)]
            if norm_rho(Gz, rho) < max(2.0 ** (-(ctx.prec - 120)), 1e-290):
                break
            dz = Jn.mid().solve(arb_mat([[g.mid()] for g in Gz])).entries()
            z = punkt([z[k] - dz[k] for k in range(M + 1)])
        Fz, Jz = F_J(z, M)
        G.append([Fz[k] - (1 - s) * cm * eM[k] for k in range(M + 1)])
        V.append(z)
        Jk.append(Jz)
    Vm = arb_mat([[knoten[i] ** n for n in range(NC)] for i in range(NC)])
    Vinv = Vm.inv()
    ahat = [[sum((Vinv[n, i] * V[i][k] for i in range(NC)), arb(0)) for k in range(M + 1)] for n in range(NC)]
    z0 = ahat[0]
    z1 = [sum((ahat[n][k] for n in range(NC)), arb(0)) for k in range(M + 1)]
    Y = Jk[NC // 2].mid().inv().mid()
    absY = arb_mat([[ob(Y[i, j]) for j in range(M + 1)] for i in range(M + 1)])
    lam = 2 / pi * arb(NC).log() + 1
    d1n = norm_rho(a1p, rho)
    fak = [arb(1)]
    for k in range(1, M + 1):
        fak.append(fak[-1] * k)
    bester = None
    for f in [0.01, 0.1, 1.0, 10.0]:
        Rk = arb(max(4.0, min(1e30, f / max(d1n * M, 1e-300))))
        A_R = [sum((ob(ahat[n][k]) * Rk ** n for n in range(NC)), arb(0)) for k in range(M + 1)]
        PRk, compRk, fpfRk = majorante_sup(A_R, M)
        faktk = 2 * arb(2) ** (-2 * NC) / (Rk - 1) ** NC
        remGk = [faktk * (compRk[k] + 1 / fak[k] + ((1 + Rk) * ob(cm) if k == M else 0)) for k in range(M + 1)]
        guete = max(float(remGk[k] / w[k]) for k in range(M + 1))
        if bester is None or guete < bester[0]:
            bester = (guete, Rk, PRk, fpfRk, faktk, remGk)
    _, R, PR, fpfR, fakt, remG = bester
    remJ = arb_mat(M + 1, M + 1)
    for k in range(M + 1):
        for j in range(M + 1):
            remJ[k, j] = fakt * (PR[j][k] + (fpfR[k - j] if j <= k else 0))
    Gb = []
    for k in range(M + 1):
        mx = ob(G[0][k])
        for i in range(1, NC):
            mx = amax(mx, ob(G[i][k]))
        Gb.append(lam * mx + remG[k] + (C.rad() if k == M else 0))
    T1 = matvec(absY, Gb)
    eins = arb_mat([[1 if i == j else 0 for j in range(M + 1)] for i in range(M + 1)])
    Eb = arb_mat(M + 1, M + 1)
    for i in range(NC):
        E = eins - Y * Jk[i]
        for k in range(M + 1):
            for j in range(M + 1):
                e = ob(E[k, j])
                Eb[k, j] = e if i == 0 else amax(Eb[k, j], e)
    YremJ = absY * remJ
    for k in range(M + 1):
        for j in range(M + 1):
            Eb[k, j] = lam * Eb[k, j] + YremJ[k, j]
    Ebw = matvec(Eb, w)
    t1 = max(float(T1[k] / w[k]) for k in range(M + 1))
    theta = max(float(Ebw[k] / w[k]) for k in range(M + 1))
    s0 = max(float((ob(P0m[k] - z0[k]) + P0[k].rad()) / w[k]) for k in range(M + 1))
    A1 = [sum((ob(ahat[n][k]) for n in range(NC)), arb(0)) for k in range(M + 1)]
    P1, _, fpf1 = majorante_sup(A1, M)
    ok, r = False, max(2 * t1 / max(1 - theta, 1e-9), 2 * s0, 1e-300)
    info = {}
    for versuch in range(12):
        Bh = [A1[k] + arb(r) * w[k] for k in range(M + 1)]
        PB, _, fpfB = majorante_sup(Bh, M)
        Jb = arb_mat(M + 1, M + 1)
        for k in range(M + 1):
            for j in range(M + 1):
                Jb[k, j] = PB[j][k] - P1[j][k] + ((fpfB[k - j] - fpf1[k - j]) if j <= k else 0)
        T3 = matvec(absY, matvec(Jb, [arb(r) * t for t in w]))
        diff = [arb(r) * w[k] - T1[k] - arb(r) * Ebw[k] - T3[k] for k in range(M + 1)]
        reserve = min(float(diff[k] / (arb(r) * w[k])) for k in range(M + 1))
        streng = all(diff[k] > 0 for k in range(M + 1))
        start_ok = all((ob(P0m[k] - z0[k]) + P0[k].rad()) < arb(r) * w[k] for k in range(M + 1))
        info = dict(r=r, reserve=reserve, start_ok=start_ok,
                    t3=max(float(T3[k] / (arb(r) * w[k])) for k in range(M + 1)))
        info["streng"] = streng
        if streng and start_ok:
            ok = True
            break
        r *= 4
    if not ok:
        return False, None, dict(M=M, NC=NC, c=float(ob(cm)), d1=d1n, t1=t1, theta=theta, **info, t=time.time() - t0)
    z = punkt(z1)
    for _ in range(3):
        Fz, Jz = F_J(z, M)
        z = punkt([z[k] - t for k, t in enumerate(Jz.mid().solve(arb_mat([[f.mid()] for f in Fz])).entries())])
    Fz, Jz = F_J(z, M)
    Y1 = Jz.mid().inv().mid()
    korr = matvec(absY, [ob(f) for f in Fz])
    r1 = [100 * ob(korr[k]) + arb(2) ** (-ctx.prec // 2) * w[k] for k in range(M + 1)]
    Zb = [arb(z[k].mid(), r1[k]) for k in range(M + 1)]
    FZ, JZ = F_J(Zb, M)
    dvec = arb_mat([[Zb[k] - z[k]] for k in range(M + 1)])
    Kb = arb_mat([[z[k]] for k in range(M + 1)]) - Y1 * arb_mat([[f] for f in Fz]) + (eins - Y1 * JZ) * dvec
    Kb = [Kb[k, 0] for k in range(M + 1)]
    innen = all((Kb[k].mid() - Kb[k].rad() > Zb[k].mid() - Zb[k].rad()) and
                (Kb[k].mid() + Kb[k].rad() < Zb[k].mid() + Zb[k].rad()) for k in range(M + 1))
    in_roehre = all(ob(Zb[k] - z1[k]) + Zb[k].rad() < arb(r) * w[k] for k in range(M + 1))
    breite = max(float(Kb[k].rad()) * rho ** k for k in range(M + 1))
    proto = dict(M=M, NC=NC, ok=bool(innen and in_roehre), c=float(ob(cm)), d1=d1n, R=float(R), t1=t1, theta=theta, **info,
                 innen=innen, in_roehre=in_roehre, breite_rho=breite, prec=ctx.prec, t=round(time.time() - t0, 1))
    return proto["ok"], Kb, proto
