"""b41 (Runde 4): Zerlegung des bilinearen Pol-Kerns über die Ableitungsidentität (b40):
  Ξ(ζ,ζ') := −½ J⁻¹T_N D²Φ(f)[E_ζ,E_ζ']  =  ¼ Π_N[(E_ζE_ζ')'/f']  −  (N+1)(E_ζ'(0)ρ_ζ + E_ζ(0)ρ_ζ')/(4f'(0)) · t_N,
  E_ζ := Π_N k_ζ − k_ζ,  ρ_ζ := [x^{N+1}] J_∞E_ζ,  t_N := J_N⁻¹ x^N  (Spitzen-/Kopfantwort).
Gemessen (alles mit python-flint, arb/acb): Übereinstimmung direkt ↔ zerlegt; Grössen beider Teile;
T_N := sup_D̄|t_N|/|L|^N; (H*) auf D̄ und auf dem äusseren Bogen von Ḡ' ((H*_G)); |E_ζ(0)|, |ρ_ζ||L|^{N+1}.
Aufruf: python b41_zerlegung.py N1 N2 ...   (Ausgabe zeilenweise, Zwischenstände in b41_<N>.json)"""
import json, sys, time
import mpmath as mp
from flint import arb, acb, arb_poly, acb_poly, arb_mat, acb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
import os
ctx.prec = int(os.environ.get("PREC", "800"))
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
K = json.load(open(G + "kneser_halbexp_koeff_1000.json"))["a"]
mp.mp.dps = 50
Lm = -mp.lambertw(-1); Lm = mp.conj(Lm) if mp.im(Lm) < 0 else Lm
L = acb(arb(mp.nstr(mp.re(Lm), 45)), arb(mp.nstr(mp.im(Lm), 45)))
aL = abs(L)
a0 = arb(K[0]); R0 = abs(L - a0)
psi = ((mp.arg(Lm - mp.mpf(K[0])) - mp.pi / 2) + (mp.arg(Lm) + mp.pi / 2)) / 2   # Winkelhalbierende Aussenkeil
phiL = float(mp.arg(Lm - mp.mpf(K[0])))
kreis = [aL * acb(arb(mp.nstr(mp.cos(2 * mp.pi * j / 240), 40)), arb(mp.nstr(mp.sin(2 * mp.pi * j / 240), 40))) for j in range(240)] + [L, L.conjugate()]
bogen = []
for j in range(1, 160):
    ph = -phiL + 2 * phiL * j / 160
    z = a0 + R0 * acb(arb(mp.nstr(mp.cos(ph), 40)), arb(mp.nstr(mp.sin(ph), 40)))
    if float(abs(z).mid()) > float(aL.mid()):
        bogen.append(z)


def A(p):
    return acb_poly([acb(c) for c in p.coeffs()])


def tr(p, n):           # mod x^n
    c = p.coeffs()[:n]
    return acb_poly(c) if c else acb_poly([0])


def coef(p, k):
    c = p.coeffs()
    return c[k] if k < len(c) else acb(0)


def inv_ser(p, n):      # 1/p mod x^n (Newton)
    y = acb_poly([1 / coef(p, 0)])
    m = 1
    while m < n:
        m = min(2 * m, n)
        y = tr(y * (acb_poly([2]) - tr(p, m) * y), m)
    return y


def supv(p, pts):
    return max(float(abs(p(z)).mid()) for z in pts)


def dist_G(z):
    return min(float((abs(z) - aL).mid()), float((abs(z - a0) - R0).mid()))


def mv(Ji, vec):        # acb_mat * Spaltenvektor (Liste acb) → Liste
    r = Ji * acb_mat([[c] for c in vec])
    return [r[k, 0] for k in range(r.nrows())]


def lauf(N):
    t0 = time.time()
    M = N + 2                                          # Reihen mod x^{N+2}
    f = arb_poly([arb(c) for c in K[:M + 1]])
    pw = [arb_poly([1])]
    for _ in range(M):
        pw.append((pw[-1] * f).truncate(M))
    ex = [arb(1)]
    for k in range(1, M + 1):
        ex.append(ex[-1] / k)
    fp = A(f.derivative()); fpA = tr(fp, M)
    g = tr(acb_poly([acb(c) for c in ex[:M]]) * inv_ser(fpA, M), M)
    gp = g.derivative()
    ifp = inv_ser(fpA, M)
    fpp_f = tr(gp * ifp, M)
    Jm = arb_mat(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            cpw = pw[j].coeffs()
            Jm[k, j] = (cpw[k] if k < len(cpw) else arb(0)) + (coef(g, k - j).real if j <= k else arb(0))
    Ji = acb_mat(Jm.inv())
    PW = acb_mat(M, M)                                 # PW[k, j] = [x^k] f^j, k,j < M
    for k in range(M):
        for j in range(M):
            cpw = pw[j].coeffs()
            PW[k, j] = acb(cpw[k]) if k < len(cpw) else acb(0)
    fA = A(f)
    tN = mv(Ji, [acb(0)] * N + [acb(1)])
    tNp = acb_poly(tN)
    TN = supv(tNp, kreis) / float((aL ** N).mid())
    w = [float(abs(c).mid()) * float((aL ** k).mid()) for k, c in enumerate(tN)]
    schw = sum(k * x for k, x in enumerate(w)) / sum(w) / N
    print(f"N={N}: Aufbau {time.time()-t0:.0f}s; Kopfantwort T_N = sup|J⁻¹x^N|/|L|^N = {TN:.4g} (Schwerpunkt k/N {schw:.2f})", flush=True)

    def bereite(zeta):
        izf = inv_ser(acb_poly([zeta]) - tr(fA, M), M)                # 1/(ζ − f)
        kz = acb_poly([zeta ** (-(m + 1)) for m in range(M)])          # k_ζ (Jets bei 0)
        rhs = tr(izf + tr(g * kz, M), N + 1)
        P = mv(Ji, [coef(rhs, k) for k in range(N + 1)])
        Pp = acb_poly(P)
        E = tr(Pp - kz, M)                                             # E_ζ mod x^{N+2}
        Pf = PW * acb_mat([[c] for c in P] + [[acb(0)]])               # P∘f mod x^{N+2}
        Ef = tr(acb_poly([Pf[k, 0] for k in range(M)]) - izf, M)        # E_ζ∘f
        r = tr(Ef + g * E, M)                                          # J_∞E_ζ
        rest = max(float(abs(coef(r, k)).mid()) * float((aL ** k).mid()) for k in range(N + 1))
        rho = coef(r, N + 1)
        return dict(P=Pp, E=E, Ef=Ef, rho=rho, rest=rest)

    orte = {"L,t=0.5": 0.5, "L,t=2": 2.0, "L,t=8": 8.0}
    Z = {}
    for nm, t in orte.items():
        Z[nm] = L + acb(arb(mp.nstr(t / N * mp.cos(psi), 40)), arb(mp.nstr(t / N * mp.sin(psi), 40)))
    Z["L̄,t=2"] = Z["L,t=2"].conjugate()
    Z["links,t=2"] = acb(-(aL + arb(2) / N))
    Z["rechts,t=2"] = acb(a0 + R0 + arb(2) / N)
    D = {}
    for nm, z in Z.items():
        D[nm] = d = bereite(z)
        dz = dist_G(z)
        hD = dz * supv(d["P"], kreis); hG = dz * supv(d["P"], bogen) if bogen else float("nan")
        print(f"   ζ={nm}: dist·sup_D̄|Π_N k_ζ| = {hD:.4f}, dist·sup_(äusserer Bogen)|Π_N k_ζ| = {hG:.4f};"
              f" |E_ζ(0)| = {float(abs(coef(d['E'],0)).mid()):.3e}, |ρ_ζ||L|^(N+1) = {float((abs(d['rho'])*aL**(N+1)).mid()):.3e},"
              f" Kontrolle max_k≤N |[J_∞E]_k||L|^k = {d['rest']:.1e}", flush=True)
    f0 = coef(fpA, 0)
    erg = []
    for (n1, n2) in (("L,t=0.5", "L,t=0.5"), ("L,t=2", "L,t=2"), ("L,t=8", "L,t=8"), ("L,t=2", "L̄,t=2"),
                     ("L,t=2", "links,t=2"), ("links,t=2", "links,t=2"), ("rechts,t=2", "L,t=2")):
        a, b = D[n1], D[n2]
        # direkt: D²Φ[a,b] = (a'∘f) b + (b'∘f) a + (f''∘f) a b,   a'∘f = (a∘f)'/f'
        daf = tr(tr(a["Ef"].derivative(), N + 1) * ifp, N + 1)
        dbf = tr(tr(b["Ef"].derivative(), N + 1) * ifp, N + 1)
        Ea, Eb = tr(a["E"], N + 1), tr(b["E"], N + 1)
        D2 = tr(daf * Eb + dbf * Ea + tr(fpp_f * tr(Ea * Eb, N + 1), N + 1), N + 1)
        Xd = acb_poly(mv(Ji, [-coef(D2, k) / 2 for k in range(N + 1)]))
        # zerlegt: v = (E_a E_b)'/f',  v∘f = ((E_a∘f)(E_b∘f))'/(f'·g)
        prod = tr(a["E"] * b["E"], M); prodf = tr(a["Ef"] * b["Ef"], M)
        v = tr(tr(prod.derivative(), N + 1) * ifp, N + 1)
        vf = tr(tr(prodf.derivative(), N + 1) * inv_ser(tr(fpA * g, N + 1), N + 1), N + 1)
        rhs = tr(vf + g * v, N + 1)
        PiV = acb_poly(mv(Ji, [coef(rhs, k) for k in range(N + 1)]))
        teil1 = acb_poly([c / 4 for c in PiV.coeffs()])
        ctop = (N + 1) * (coef(b["E"], 0) * a["rho"] + coef(a["E"], 0) * b["rho"]) / (4 * f0)
        teil2 = acb_poly([-ctop * c for c in tN])
        dd = dist_G(Z[n1]) * dist_G(Z[n2])
        sX = supv(Xd, kreis); s1 = supv(teil1, kreis); s2 = supv(teil2, kreis)
        dif = supv(Xd - teil1 - teil2, kreis)
        erg.append(dict(paar=f"{n1}×{n2}", direkt=dd * sX, teil_PiN=dd * s1, teil_kopf=dd * s2, abweichung=dd * dif))
        print(f"   {n1}×{n2}: dist·dist'·sup|Ξ| = {dd*sX:.3e} | Π_N-Teil {dd*s1:.3e} | Kopfteil {dd*s2:.3e} | Abweichung {dd*dif:.1e}", flush=True)
    json.dump(dict(N=N, T_N=TN, schwerpunkt=schw, paare=erg), open(f"b41_{N}.json", "w"), indent=1)
    print(f"N={N} fertig nach {time.time()-t0:.0f}s", flush=True)


for N in [int(s) for s in sys.argv[1:]]:
    lauf(N)
