"""Runde 2, Schritt 6: Sekanten-Operatoren 𝒥_p u = u∘p + M_{p,f}·u, M_{p,f} = (f∘p − f∘f)/(p − f) = (1+τ²)/((1−τp)(1−τf)),
Π^{(p)} = (T_N𝒥_p|_{P_N})⁻¹T_N𝒥_p.  (a) Identität p_N − f = r·E^{(p_N)}(ζ_p,·), r = (1+τ²)/τ² (Lemma M4) prüfen.
(b) (H*) für 𝒥_p gleichmässig über eine Kugel um f: p = p_N + η·q (q zufällig, ‖q‖_{D̄} = 1), Grösse dist·max|Π^{(p)}k_ζ|."""
import sys, time, random
import mpmath as mp
from flint import arb, acb, arb_mat, arb_poly, ctx
from kern import zweig, g_moebius_ell
from moebius_lin import auswerten

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 40
c = mp.mpf(sys.argv[1])
Ns = [int(x) for x in sys.argv[2].split(",")]
etas = [mp.mpf(x) for x in sys.argv[3].split(",")]
BITS = 900
ctx.prec = BITS
tau = mp.tan(c / 2)
ta = arb(mp.nstr(tau, 120))
loes = zweig(g_moebius_ell(c), max(Ns), bits=BITS)
K = 192
zs = [acb(mp.nstr(mp.cos(2 * mp.pi * j / K), 30), mp.nstr(mp.sin(2 * mp.pi * j / K), 30)) for j in range(K)]


def sek_matrix(pk, N):
    """T_N 𝒥_p auf P_N und die Reihen von p, M (arb)."""
    p = arb_poly(pk)
    fk = [ta] + [(1 + ta * ta) * ta ** (k - 1) for k in range(1, N + 1)]
    f = arb_poly(fk)
    P = [arb_poly([1])]
    for _ in range(N):
        P.append((P[-1] * p).truncate(N + 1))
    # 1/(1 − τp): p = p0 + h, 1/(1 − τp0 − τh) = Σ (τh)^n/(1 − τp0)^{n+1}
    def inv1m(q):
        q0 = q[0]
        h = q - arb_poly([q0])
        s, hp = arb_poly([0]), arb_poly([1])
        for n in range(N + 1):
            s += hp * (ta ** n / (1 - ta * q0) ** (n + 1))
            hp = (hp * h).truncate(N + 1)
        return s.truncate(N + 1)
    M = ((1 + ta * ta) * inv1m(p) * inv1m(f)).truncate(N + 1)
    A = arb_mat(N + 1, N + 1)
    for j in range(N + 1):
        col = P[j] + (M * arb_poly([0] * j + [1])).truncate(N + 1)
        for k in range(N + 1):
            A[k, j] = col[k]
    return A, p, M


def pol_antw(A_inv, p, M, N, zeta):
    zeta = acb(zeta)
    p0 = p[0]
    h = p - arb_poly([p0])
    hp, v = arb_poly([1]), [acb(0)] * (N + 1)
    for n in range(N + 1):
        fac = 1 / (zeta - p0) ** (n + 1)
        for k in range(N + 1):
            v[k] += hp[k] * fac
        hp = (hp * h).truncate(N + 1)
    for k in range(N + 1):
        s = acb(0)
        for i in range(k + 1):
            s += M[i] / zeta ** (k - i + 1)
        v[k] += s
    return [sum((A_inv[k, j] * v[j] for j in range(N + 1)), acb(0)) for k in range(N + 1)]


def maxD(coef, zeta=None):
    best = mp.mpf(0)
    for z in zs:
        val = auswerten(coef, z) - (1 / (acb(zeta) - z) if zeta is not None else 0)
        best = max(best, abs(mp.mpc(val.real.mid().str(20, radius=False), val.imag.mid().str(20, radius=False))))
    return best


R0 = mp.sqrt(1 + tau ** 2)
om = mp.pi - c / 2
alpha = om / mp.pi
random.seed(1)
for N in Ns:
    t0 = time.time()
    pk = [arb(mp.nstr(x, 200)) for x in loes[N]]
    # (a) Identität an p = p_N
    A, p, M = sek_matrix(pk, N)
    Ai = A.inv()
    zp = acb(mp.nstr(1 / tau, 60))
    E = pol_antw(Ai, p, M, N, zp)
    r = (1 + ta * ta) / (ta * ta)
    dmax = mp.mpf(0)
    for z in zs:
        lhs = auswerten([acb(x) for x in pk], z) - (acb(ta) + z) / (1 - acb(ta) * z)
        rhs = r * (auswerten(E, z) - 1 / (zp - z))
        d = lhs - rhs
        dmax = max(dmax, abs(mp.mpc(d.real.mid().str(20, radius=False), d.imag.mid().str(20, radius=False))))
    print(f"N={N}: Lemma M4 geprüft: max|p_N − f − r·E^(p_N)(ζ_p,·)| = {mp.nstr(dmax, 3)}", flush=True)
    # (b) Gleichmässigkeit von (H*) über eine Kugel
    for eta in etas:
        q = [random.gauss(0, 1) / (k + 1) ** 2 for k in range(N + 1)]
        qn = sum(abs(x) for x in q)
        pk2 = [pk[k] + arb(mp.nstr(eta * q[k] / qn, 30)) for k in range(N + 1)]
        A, p, M = sek_matrix(pk2, N)
        Ai = A.inv()
        werte = []
        for (b, nrm, ecke) in ((mp.mpc(0, 1), mp.expj(mp.pi / 2 + c / 4), True), (mp.mpc(-1, 0), mp.mpc(-1, 0), False),
                               (tau + R0, mp.mpc(1, 0), False)):
            for t in (mp.mpf("0.5"), mp.mpf(8)):
                d = t / (N ** alpha if ecke else N)
                zeta = b + d * nrm
                coef = pol_antw(Ai, p, M, N, acb(mp.nstr(mp.re(zeta), 30), mp.nstr(mp.im(zeta), 30)))
                werte.append(mp.nstr(d * maxD(coef), 3))
        Ep = pol_antw(Ai, p, M, N, zp)
        print(f"   η = {mp.nstr(eta, 2)}: dist·max|Π^(p)k_ζ| (Ecke t=.5/8, links, rechts) = {werte};  max|E^(p)(ζ_p,·)| = {mp.nstr(maxD(Ep, zp), 3)}  ({time.time() - t0:.0f}s)", flush=True)
