r"""Strang C, Runde 3, c19: (H_R) über Rationalität in ζ + Bernstein–Walsh.

Lemma 12.1 (exakt): y_ζ := A^{-1}B_N[(ζ−x)^{-1}]_{>N} = ζ^{-N-1} Σ_{m=0}^{N} A^{-1}[q^{N+1} g^m]_{≤N} (ζ−a0)^{-m-1},  g = q − a0.
  ⇒ y_ζ rational in ζ, Pole nur bei 0 und a0, y_ζ = O(ζ^{-N-2}) bei ∞.
Satz 12.2 (Bernstein–Walsh): R_N(ζ) ≤ S_N · exp(−(N+2) g_E(ζ)),  S_N := max_{∂W} R_N,  E = Ĉ \ W̄,  g_E Greensche Funktion
  von E mit Pol ∞. g_E explizit: Möbius M(ζ) = (ζ−L)/(ζ−L̄) bildet E auf einen Sektor ab, Potenz bildet ihn auf eine Halbebene ab.

Dieses Skript: (1) g_E, Öffnungswinkel, Asymptotik bei L; (2) S_N auf ∂W (fein, verdichtet bei L) für mehrere N;
(3) Probe der Ungleichung an Schlitz-, Kontur- und Gitterpunkten; (4) Schranke für α_N aus Satz 9.5 + 12.2 (C_f = 0.95, max|f| = 3.6).
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, mpc, pi, exp, log, arg, sqrt, conj, quad, gamma
from flint import arb_mat, arb, ctx

mp.dps = 130
ctx.prec = 950
HIER = os.path.dirname(os.path.abspath(__file__))
G = os.path.join(HIER, "..", "gemeinsam")
KN = [mpf(x) for x in json.load(open(os.path.join(G, "kneser_halbexp_koeff_1000.json"), encoding="utf-8"))["a"]]
L = mpc("0.31813150520476413", "1.3372357014306895")
Lq = conj(L)
RHO = abs(L)
a0 = KN[0]
R0 = abs(L - a0)
SIGMA = (1 + 2j * pi / L).real
CF, MF = mpf("0.95"), mpf("3.6")

# ---------- (1) Greensche Funktion von E = Ĉ \ (D(0,ρ) ∪ D(a0,R0))
def M(z):
    return (z - L) / (z - Lq)

th_a = arg(M(-RHO))            # Bild des Randbogens von D(0,ρ) (aussen von D(a0,R0))
th_b = arg(M(a0 + R0))         # Bild des Randbogens von D(a0,R0) (aussen von D(0,ρ))
th_inf = arg(mpc(1))           # M(∞) = 1, arg 0
# Sektor, der arg 0 enthält: zwischen lo und hi
lo, hi = sorted([th_a, th_b])
if not (lo < th_inf < hi):
    lo, hi = hi, lo + 2 * pi
OEFF = hi - lo
ALPHA = OEFF / pi


def tau(z):
    w = M(z) * exp(-1j * lo)
    a = arg(w)
    if a < 0:
        a += 2 * pi
    return abs(w) ** (1 / ALPHA) * exp(1j * a / ALPHA)


TAU0 = tau(mpf(10) ** 30)


def gE(z):
    t = tau(z)
    return log(abs((t - conj(TAU0)) / (t - TAU0)))


print(f"Öffnungswinkel von E bei L: {mp.nstr(OEFF * 180 / pi, 8)}°, α = {mp.nstr(ALPHA, 8)}")
print("Probe g_E: auf ∂W (soll 0):", mp.nstr(gE(-RHO), 3), mp.nstr(gE(a0 + R0), 3),
      " | g_E(10^4) − log 10^4 =", mp.nstr(gE(mpf(10) ** 4) - log(mpf(10) ** 4), 8), "(Robin-Konstante)")
for u in (1e-3, 1e-2, 0.05, 0.1, 0.2):
    z = L * (1 + u)
    print(f"  g_E(L(1+{u})) = {mp.nstr(gE(z), 6)},  g_E/u^(1/α) = {mp.nstr(gE(z) / mpf(u) ** (1 / ALPHA), 6)};"
          f"  φ = min(log|ζ|/ρ, log|ζ−a0|/R0) = {mp.nstr(min(log(abs(z) / RHO), log(abs(z - a0) / R0)), 6)}")
sys.stdout.flush()


# ---------- Pol-Antwort
def mul(a, b, n):
    c = [mpf(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[: n + 1 - i]):
            c[i + j] += ai * bj
    return c


def baue(N):
    q = KN[: N + 1]
    pot = [[mpf(1)] + [mpf(0)] * N]
    for _ in range(N + 1):
        pot.append(mul(pot[-1], q, N))
    dq = [k * q[k] for k in range(1, N + 1)]
    dpp = [sum(dq[j] * pot[j][k] for j in range(len(dq))) for k in range(N + 1)]
    A = arb_mat([[arb(mp.nstr(pot[j][k] + (dpp[k - j] if j <= k else 0), 125)) for j in range(N + 1)] for k in range(N + 1)])
    Ai = A.inv()
    # z_m = A^{-1}[q^{N+1} g^m]_{≤N}
    g = [mpf(0)] + q[1:]
    w = pot[N + 1][:]
    Z = []
    Wm = arb_mat([[arb(0)] * (N + 1) for _ in range(N + 1)])
    spalten = []
    for m in range(N + 1):
        spalten.append(w[:])
        w = mul(w, g, N)
    Wm = arb_mat([[arb(mp.nstr(spalten[m][j], 125)) for m in range(N + 1)] for j in range(N + 1)])
    Zm = Ai * Wm
    Z = [[mpf(Zm[j, m].mid().str(125, radius=False)) for j in range(N + 1)] for m in range(N + 1)]
    qq = [sum(q[j] * pot[j][k] for j in range(N + 1)) for k in range(N + 1)]
    F = [qq[k] - mpf(1) / mp.factorial(k) for k in range(N + 1)]
    d1 = Ai * arb_mat([[arb(mp.nstr(-x, 125))] for x in F])
    alpha = sum(abs(mpf(d1[j, 0].mid().str(120, radius=False))) * RHO ** j for j in range(N + 1))
    return Z, alpha


def RN(Z, N, z):
    s = 1 / (z - a0)
    # y = z^{-N-1} Σ_m Z[m] s^{m+1}
    acc = [mpc(0)] * (N + 1)
    sp = s
    for m in range(N + 1):
        zm = Z[m]
        for j in range(N + 1):
            acc[j] += zm[j] * sp
        sp *= s
    return sum(abs(acc[j]) * RHO ** j for j in range(N + 1)) / abs(z) ** (N + 1)


def rand_W(n1=120, n2=120):
    pts = []
    # Bogen von D(0,ρ) ausserhalb D(a0,R0): Winkel von arg L bis 2π − arg L (über π)
    aL = arg(L)
    for k in range(n1 + 1):
        t = aL + (2 * pi - 2 * aL) * k / n1
        pts.append(RHO * exp(1j * t))
    # Bogen von D(a0,R0) ausserhalb D(0,ρ): Winkel um a0 von arg(L−a0) abwärts über 0 bis −arg(L−a0)
    b = arg(L - a0)
    for k in range(n2 + 1):
        t = b - 2 * b * k / n2
        pts.append(a0 + R0 * exp(1j * t))
    # Verdichtung bei L auf beiden Bögen
    for d in (1e-4, 1e-3, 3e-3, 1e-2, 3e-2):
        pts.append(RHO * exp(1j * (aL + d)))
        pts.append(a0 + R0 * exp(1j * (b - d)))
    return pts


NS = [int(x) for x in (sys.argv[1].split(",") if len(sys.argv) > 1 else "40,80,120,160")]
RAND = rand_W()
for N in NS:
    Z, alphaN = baue(N)
    werte = [(RN(Z, N, z), z) for z in RAND]
    SN, zmax = max(werte, key=lambda t: t[0])
    beiL = RN(Z, N, L)
    print(f"\nN={N}: S_N = max_∂W R_N = {mp.nstr(SN, 5)} bei ζ = {mp.nstr(zmax, 5)} (|ζ−L| = {mp.nstr(abs(zmax - L), 3)}); "
          f"R_N(L) = {mp.nstr(beiL, 5)}; α_N = {mp.nstr(alphaN, 5)}")
    # (3) Probe der Ungleichung
    proben = [L * (1 + u) for u in (1e-3, 1e-2, 0.03, 0.1, 0.3)] + [mpc(1.9, 0), mpc(-1.6, 0.6), mpc(0.2, 1.8), mpc(2.4, 2.4)]
    zeile = []
    for z in proben:
        r = RN(Z, N, z)
        schr = SN * exp(-(N + 2) * gE(z))
        zeile.append(f"{mp.nstr(z, 3)}: R/Schranke={mp.nstr(r / schr, 3)}")
    print("   Probe R_N ≤ S_N e^{−(N+2)g_E}: " + "; ".join(zeile))
    # (4) Schranke für α_N: Schlitz (u bis u0) + Kontur ∂(W+0.1) (Länge ~11, min g_E dort)
    u0 = mpf("0.2")
    schlitz = quad(lambda u: CF * u ** SIGMA * exp(-(N + 2) * gE(L * (1 + u))) * RHO, [0, mpf(1) / N, mpf(10) / N, u0])
    schlitz = 2 * SN * schlitz / (2 * pi)
    print(f"   Schranke Schlitzteil 2·(1/2π)∫ C_f u^σ S_N e^(−(N+2)g_E) |L| du = {mp.nstr(schlitz, 4)};  α_N = {mp.nstr(alphaN, 4)};"
          f"  Verhältnis {mp.nstr(schlitz / alphaN, 3)}")
    sys.stdout.flush()
