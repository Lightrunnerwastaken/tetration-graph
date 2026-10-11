"""b46 (Runde 5, Abgleich 6.5): Ist Strang Es Sekanten-Abbildung T(P) := Π^{(P)}_N f_K auch in Richtung der Eckmode
Lipschitz-klein — oder verstärkt die Eckmode (σ₁ ~ e^{0.012N}) Störungen der festen Kugel?
T(P) = (T_N 𝒥_P|_{P_N})⁻¹ T_N[eˣ + M_P·P],  𝒥_P v = v∘P + M_P v,
M_P = Σ_{j≥0} (f^{(j+1)}∘f)(P − f)^j/(j+1)!  mit  f'∘f = g = eˣ/f',  f^{(j+1)}∘f = (f^{(j)}∘f)'/f'  (keine Division durch P − f).
P = p_N + η·q, q ∈ {u₁ (Ausgangsvektor der Eckmode von J_N(f)⁻¹ in ℓ²_{|L|}), zufällige Richtung}, normiert Σ|q_k|²|L|^{2k} = 1.
Gemessen: sup_{|z|=|L|}|T(P) − T(p_N)| / η  (Lipschitz-Quotient) und sup|T(P) − f_K|.
Aufruf: python b46_sekante_eckrichtung.py N1 N2 …"""
import json, os, sys, time
import numpy as np
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
K = json.load(open(G + "kneser_halbexp_koeff_1000.json"))["a"]
PN = json.load(open(G + "pN_N100_bis_1000.json"))["p"]
mp.mp.dps = 60
d4 = json.load(open(G + "kneser_halbexp_koeff_4000.json"))
K4 = [mp.mpf(c) for c in (d4["a"] if isinstance(d4, dict) else d4)]
Lm = -mp.lambertw(-1); Lm = mp.conj(Lm) if mp.im(Lm) < 0 else Lm
aL = abs(Lm); aLf = float(aL)
kreis = [aL * mp.expjpi(2 * mp.mpf(j) / 240) for j in range(240)] + [Lm, mp.conj(Lm)]


def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s


fKw = [horner(K4, z) for z in kreis]


def coef(p, k):
    c = p.coeffs()
    return c[k] if k < len(c) else arb(0)


def tr(p, n):
    c = p.coeffs()[:n]
    return arb_poly(c) if c else arb_poly([0])


def inv_ser(p, n):
    y = arb_poly([1 / coef(p, 0)])
    m = 1
    while m < n:
        m = min(2 * m, n)
        y = tr(y * (arb_poly([2]) - tr(p, m) * y), m)
    return y


def mpl(x):
    return mp.mpf(x.mid().str(50, radius=False))


def supL(c, ref=None):
    cm = [mpl(x) if hasattr(x, "mid") else mp.mpf(x) for x in c]
    if ref is None:
        return max(abs(horner(cm, z)) for z in kreis)
    return max(abs(horner(cm, z) - r) for z, r in zip(kreis, ref))


for N in [int(s) for s in sys.argv[1:]] or [200, 300]:
    t0 = time.time()
    ctx.prec = 700 + 2 * N
    M = N + 1
    f = arb_poly([arb(c) for c in K[:M + 1]])
    ex = [arb(1)]
    for k in range(1, M + 1):
        ex.append(ex[-1] / k)
    fp = tr(f.derivative(), M); ifp = inv_ser(fp, M)
    g = tr(arb_poly(ex[:M]) * ifp, M)
    # Ableitungen f^{(j+1)}∘f, j = 0..J
    J = 40
    D = [g]
    for j in range(J):
        D.append(tr(tr(D[-1].derivative(), M) * ifp, M))
    fakt = [arb(1)]
    for j in range(1, J + 2):
        fakt.append(fakt[-1] * j)
    # J_N(f) und Eckmode
    pwf = [arb_poly([1])]
    for _ in range(N):
        pwf.append((pwf[-1] * f).truncate(M))
    Jf = arb_mat(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            Jf[k, j] = coef(pwf[j], k) + (coef(g, k - j) if j <= k else arb(0))
    Jfi = Jf.inv()
    A = np.array([[float(Jfi[i, j].mid().str(20, radius=False)) for j in range(N + 1)] for i in range(N + 1)])
    w = aLf ** np.arange(N + 1)
    U, sv, Vt = np.linalg.svd((w[:, None] * A) / w[None, :])
    u1 = U[:, 0] / w                                     # Koeffizienten, Σ|u_k|²|L|^{2k} = 1
    rng = np.random.default_rng(1)
    qr = rng.standard_normal(N + 1); qr = qr / np.linalg.norm(qr) / w
    pN = [arb(s) for s in PN[str(N)]]
    print(f"N={N}: Aufbau {time.time()-t0:.0f}s, σ₁ = {sv[0]:.3f}, σ₂ = {sv[1]:.3f}", flush=True)

    def T_von(Pc):
        P = arb_poly(Pc)
        d = tr(P - f, M)
        Mp = arb_poly([0]); dj = arb_poly([1])
        for j in range(J + 1):
            Mp = tr(Mp + tr(D[j] * dj, M) * (1 / fakt[j + 1]), M)
            dj = tr(dj * d, M)
        pw = [arb_poly([1])]
        for _ in range(N):
            pw.append((pw[-1] * P).truncate(M))
        Jm = arb_mat(N + 1, N + 1)
        for k in range(N + 1):
            for j in range(N + 1):
                Jm[k, j] = coef(pw[j], k) + (coef(Mp, k - j) if j <= k else arb(0))
        rhs = tr(arb_poly(ex[:M]) + tr(Mp * P, M), N + 1)
        sol = Jm.solve(arb_mat([[coef(rhs, k)] for k in range(N + 1)]))
        return [sol[k, 0] for k in range(N + 1)]

    T0 = T_von(pN)
    print(f"   T(p_N): sup|T − f_K| = {mp.nstr(supL(T0, fKw), 4)}, sup|T − p_N| = {mp.nstr(supL([T0[k] - pN[k] for k in range(N + 1)]), 3)} ({time.time()-t0:.0f}s)", flush=True)
    for name, q in (("Eckmode u₁", u1), ("zufällig", qr)):
        qn = [arb(repr(float(x))) for x in q]
        sq = supL(qn)
        for eta in (1e-4, 1e-3, 1e-2):
            Pc = [pN[k] + arb(eta) * qn[k] for k in range(N + 1)]
            Tq = T_von(Pc)
            dT = supL([Tq[k] - T0[k] for k in range(N + 1)])
            print(f"   q = {name} (sup_D̄|q| = {mp.nstr(sq,3)}), η = {eta:g}: sup|T(P) − T(p_N)|/η = {mp.nstr(dT/eta, 3)},"
                  f" sup|T(P) − f_K| = {mp.nstr(supL(Tq, fKw), 3)} ({time.time()-t0:.0f}s)", flush=True)
