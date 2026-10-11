"""b45 (Runde 5, Abgleich von Unstimmigkeiten):
(1) Bs (H*)-Messung (sup_D̄|Π_N k_ζ| mit Kopf 1/dist) gegen Cs R_N (Pol-Antwort OHNE Kopf, y_ζ = Π_N k_ζ − T_N k_ζ):
    für den linken Pol ζ = −(|L| + 2/N) und den Eckpol L + (2/N)e^{iψ}: sup_D̄|y_ζ|, Ort des Maximums, sup_D̄|Π_N k_ζ|.
(2) Cs Befund „p_N konvergiert bei x = 1.8“ (N ≤ 180) gegen b42 (Wachstum am Halbmondrand ab N ≈ 300):
    |Π_N f_K(x) − f_K(x)| bei x = 1.8 und x = 1.6 für N = 100…500 (f_K(x) = exp(f_K(log x)), Knesers Reihe).
Aufruf: python b45_abgleich.py"""
import json, os, sys, time
import mpmath as mp
from flint import arb, acb, arb_poly, acb_poly, arb_mat, acb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
K = json.load(open(G + "kneser_halbexp_koeff_1000.json"))["a"]
d4 = json.load(open(G + "kneser_halbexp_koeff_4000.json"))
lin = json.load(open(G + "lin_PiN_fK_bis_1000.json"))
mp.mp.dps = 200
K4 = [mp.mpf(c) for c in (d4["a"] if isinstance(d4, dict) else d4)]   # erst NACH dps setzen (Fehler im ersten Lauf: 15 Stellen)
Lm = -mp.lambertw(-1); Lm = mp.conj(Lm) if mp.im(Lm) < 0 else Lm
aL = abs(Lm)


def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s


# (2) zuerst (billig)
for x in (mp.mpf("1.6"), mp.mpf("1.8")):
    fx = mp.exp(horner(K4, mp.log(x)))
    zeile = []
    for N in (100, 200, 300, 400, 500):
        b = [mp.mpf(s) for s in lin[str(N)]["b"]]
        zeile.append(f"N={N}: {mp.nstr(abs(horner(b, x) - fx), 3)}")
    print(f"(2) |Π_N f_K({mp.nstr(x,2)}) − f_K| — " + ", ".join(zeile), flush=True)

# (1)
psi = ((mp.arg(Lm - mp.mpf(K[0])) - mp.pi / 2) + (mp.arg(Lm) + mp.pi / 2)) / 2
kreis = [aL * mp.expjpi(2 * mp.mpf(j) / 720) for j in range(720)]


def coef(p, k):
    c = p.coeffs()
    return c[k] if k < len(c) else (acb(0) if isinstance(p, acb_poly) else arb(0))


def tr(p, n):
    c = p.coeffs()[:n]
    return type(p)(c) if c else type(p)([0])


def inv_ser(p, n):
    y = type(p)([1 / coef(p, 0)])
    m = 1
    while m < n:
        m = min(2 * m, n)
        y = tr(y * (type(p)([2]) - tr(p, m) * y), m)
    return y


for N in (200, 300, 400):
    t0 = time.time()
    ctx.prec = 600 + 2 * N
    M = N + 1
    f = arb_poly([arb(c) for c in K[:M + 1]])
    pw = [arb_poly([1])]
    for _ in range(N):
        pw.append((pw[-1] * f).truncate(M))
    ex = [arb(1)]
    for k in range(1, M + 1):
        ex.append(ex[-1] / k)
    fp = tr(f.derivative(), M)
    g = tr(arb_poly(ex[:M]) * inv_ser(fp, M), M)
    Jm = arb_mat(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            Jm[k, j] = coef(pw[j], k) + (coef(g, k - j) if j <= k else arb(0))
    Ji = acb_mat(Jm.inv())
    fA = acb_poly([acb(c) for c in f.coeffs()])
    gA = acb_poly([acb(c) for c in g.coeffs()])
    for nm, zeta in (("links", -(aL + mp.mpf(2) / N)), ("Ecke", Lm + mp.mpf(2) / N * mp.expj(psi))):
        z = acb(arb(mp.nstr(mp.re(zeta), 60)), arb(mp.nstr(mp.im(zeta), 60)))
        izf = inv_ser(acb_poly([z]) - tr(fA, M), M)
        kz = acb_poly([z ** (-(m + 1)) for m in range(M)])
        rhs = tr(izf + tr(gA * kz, M), N + 1)
        r = Ji * acb_mat([[coef(rhs, k)] for k in range(N + 1)])
        P = [mp.mpc(mp.mpf(r[k, 0].real.mid().str(60, radius=False)), mp.mpf(r[k, 0].imag.mid().str(60, radius=False))) for k in range(N + 1)]
        y = [P[k] - zeta ** (-(k + 1)) for k in range(N + 1)]
        werte = [(abs(horner(y, w)), float(mp.arg(w) * 180 / mp.pi)) for w in kreis]
        my, ort = max(werte)
        mP = max(abs(horner(P, w)) for w in kreis)
        l1 = sum(abs(c) * aL ** k for k, c in enumerate(y))
        print(f"(1) N={N} ({time.time()-t0:.0f}s) ζ={nm}: sup_D̄|y_ζ| = {mp.nstr(my,4)} bei arg z = {ort:.1f}°, ℓ¹_ρ|y_ζ| = {mp.nstr(l1,4)},"
              f" sup_D̄|Π_N k_ζ| = {mp.nstr(mP,4)} (·dist = {mp.nstr(mP*2/N,4)})", flush=True)
