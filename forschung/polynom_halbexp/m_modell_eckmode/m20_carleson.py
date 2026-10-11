"""M20 — Carleson-Konstante C_h: ∫_Γout |g(h(e^{iφ}))|² dφ/2π ≤ C_h² ‖g‖²_{H²(D)} für g ∈ H²(D), streng (arb).
Streifen: Z := g∘ξ·(ξ'/2π)^{1/2} ∈ H²(S), ‖Z‖ = ‖g‖. Dann ∫_Γout|g∘h|² = ∫|Z(H(t))|² r(t) dt, r(t) = |ξ'(t+iπ)|/|ξ'(H(t))|.
Mit s = Re H(t) (Re H' ≥ d0 > 0) und Höhe y(s) := Im H ∈ [y1, y2]:
  |Z(s+iy)|² ≤ (1+ε)|Z(s+iy1)|² + (1+1/ε)(y2−y1)∫_{y1}^{y2}|Z'(s+iη)|²dη, Fourier: ‖Z(·+iy)‖ ≤ ‖Z‖, ‖Z'(·+iη)‖² ≤ κ‖Z‖²,
  κ = sup_{η∈[y1,y2], σ} σ² e^{−2ση}/(1+e^{−2πσ}).  ⇒  C_h² ≤ (r_max/d0)·[(1+ε) + (1+1/ε)(y2−y1)²κ].
Raster auf [−Tb, 0] (Symmetrie t ↦ −t), für t < −Tb Schwarz-/Cauchy-Schranken aus g3 (m17: max_{|M|=R}|g3| = γ3·R).
Aufruf: python m20_carleson.py [Tb] [Zellweite]
"""
import sys

from flint import acb, arb, ctx

sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 100
Tb = arb(sys.argv[1]) if len(sys.argv) > 1 else arb(6)
dz = arb(sys.argv[2]) if len(sys.argv) > 2 else arb("0.02")
pi = arb.pi()
eps = arb("0.31382")
th = arb("76.618") * pi / 180
c0 = 1 - 2 * eps * th.cos()
L = acb(th.cos(), th.sin())
Lb = L.conjugate()
mu = acb(1, 2 * eps * th.sin())
ell = mu.log()
b = mu.arg()
cc = eps * (L - Lb)
dd = c0 + eps * Lb


def f(z):
    return eps * (z + c0 / (2 * eps)) ** 2 + eps - c0 * c0 / (4 * eps)


def xi(M):
    return Lb + (Lb - L) / (M - 1)


def hinv(x):
    return (-c0 + (c0 * c0 - 4 * eps * eps + 4 * eps * x).sqrt()) / (2 * eps)


def Lam(x):
    return (1 + cc / (dd + eps * x)).log()


def Lam_abl(x):
    w = dd + eps * x
    return -cc * eps / (w * (w + cc))


def absxi_abl(zeta):
    M = L * zeta.exp()
    return 2 * th.sin() * zeta.real.exp() / abs(M - 1) ** 2


def daten(t):
    """Für reelle Kugel t: (H, Re H', r) auf der Oberkante."""
    M = -L * t.exp()
    x = xi(M)
    v = hinv(x)
    H = acb(t, pi) - Lam(v)
    Hp = 1 - Lam_abl(v) / (c0 + 2 * eps * v) * (L - Lb) * M / (M - 1) ** 2
    r = absxi_abl(acb(t, pi)) / absxi_abl(H)
    return H, Hp.real, r


N = int(float((Tb / dz).upper()) + 0.5)
y1, y2 = None, None
d0, rmax = None, arb(0)
for i in range(N):
    m = -Tb + dz * (arb(i) + arb("0.5"))
    tb = arb(m.mid(), (dz / 2).upper() * (1 + 1e-12))
    H, hp, r = daten(tb)
    lo, hi = H.imag.lower(), H.imag.upper()
    y1 = lo if y1 is None else arb(min(y1, lo))
    y2 = hi if y2 is None else arb(max(y2, hi))
    d0 = hp.lower() if d0 is None else arb(min(d0, hp.lower()))
    rmax = arb(max(rmax, r.upper()))
print(f"Raster [−{Tb}, 0]: Im H ∈ [{float(y1):.5f}, {float(y2):.5f}], Re H' ≥ {float(d0):.5f}, r ≤ {float(rmax):.5f}")
# Linkes Ende t < −Tb: |H − (t + iπ − ℓ)| ≤ γ3 e^t, |Re H' − 1| ≤ e^t·A/(R − e^t), A = max_{|M|=R}|g3|
R = arb("0.3")
gam3 = arb("0.34535")      # obere Schranke aus m17 (γ3 ≤ 0.3453442…)
A = gam3 * R
q = (-Tb).exp()
y1e, y2e = pi - b - gam3 * q, pi - b + gam3 * q
d0e = 1 - q * A / (R - q)
# r am Ende: r = e^{t − Re H}·|M_H − 1|²/|M_{t+iπ} − 1|²; e^{t − Re H} ≤ e^{Re ℓ + γ3 q}; |M_H| = e^{Re H} ≤ e^{−Tb}·e^{|Re ℓ|+γ3 q}
qH = q * (ell.real + gam3 * q).exp()
re_ = (ell.real + gam3 * q).exp() * (1 + qH) ** 2 / (1 - q) ** 2
print(f"Ende t < −{Tb}: Im H ∈ [{float(y1e.lower()):.5f}, {float(y2e.upper()):.5f}], Re H' ≥ {float(d0e.lower()):.5f}, r ≤ {float(re_.upper()):.5f}")
Y1 = arb(min(y1, y1e.lower()))
Y2 = arb(max(y2, y2e.upper()))
D0 = arb(min(d0, d0e.lower()))
RM = arb(max(rmax, re_.upper()))
# κ = sup σ² e^{−2ση}/(1+e^{−2πσ}) über η ∈ [Y1, Y2]: σ < 0: ≤ σ² e^{−2|σ|(π−Y2)}·... ≤ (1/(e(π−Y2)))²; σ > 0: ≤ (1/(eY1))²
kap = arb(max((1 / (arb(1).exp() * (pi - Y2))) ** 2, (1 / (arb(1).exp() * Y1)) ** 2).upper())
best = None
for k in range(1, 200):
    e_ = arb(k) / 1000
    val = (RM / D0) * ((1 + e_) + (1 + 1 / e_) * (Y2 - Y1) ** 2 * kap)
    if best is None or val.upper() < best[0].upper():
        best = (val, e_)
Ch = best[0].sqrt()
# m_in := min_{Γin}|m|, m = f'∘f; Γin = {θ ≤ |φ| ≤ π}, Symmetrie φ ↦ −φ; Intervallunterteilung auf [θ, π]
nphi = 20000
m_in = None
for k in range(nphi):
    a_ = th + (pi - th) * arb(k) / nphi
    b_ = th + (pi - th) * arb(k + 1) / nphi
    ph = arb(((a_ + b_) / 2).mid(), ((b_ - a_) / 2).upper() * (1 + 1e-12))
    z = acb(ph.cos(), ph.sin())
    v = abs(c0 + 2 * eps * f(z)).lower()
    m_in = v if m_in is None else arb(min(m_in, v))
print(f"min_Γin |m| ≥ {float(m_in):.6f} (zertifiziert)")
CE = (1 / m_in ** 2 + Ch ** 2).sqrt()
print(f"κ ≤ {float(kap):.4f}, ε = {float(best[1]):.3f}: C_h ≤ {float(Ch.upper()):.5f}; C_E = (1/m_in² + C_h²)^½ ≤ {float(CE.upper()):.5f}")
c = arb("0.8211")
print(f"(a″)_∞: σ_min(J_∞|P_N) ≥ c/C_E ≥ {float((c / CE).lower()):.4f} für alle N (c = 0.8211 aus Satz M6)")

# --- Verfeinerung (Kanten getrennt): ‖𝒯x‖² ≤ ‖Z0‖²/m_in² + (r/d0)∫|Z(s+iy(s))|²ds und für y1 = min Höhe, θy = y1/π:
# ‖Z(·+iy1)‖² = ∫|Ẑ|²e^{−2σy1} ≤ α‖Z0‖² + β‖Z1‖² mit α = (1−θy)s^θy, β = θy s^{θy−1} (Young, x = e^{−2σ}, beliebiges s > 0);
# Ableitungsterm ≤ κ'‖Z‖², κ' = (1+1/ε)(y2−y1)²κ. ⇒ C_E'² = max(1/m_in² + R[(1+ε)α + κ'], R[(1+ε)β + κ']), R = r_max/d0.
def ce_fein(m_in, RM, D0, Y1, Y2, kap):
    R = RM / D0
    th_y = Y1 / pi
    beste = None
    for ke in range(1, 120):
        e_ = arb(ke) / 1000
        kp = (1 + 1 / e_) * (Y2 - Y1) ** 2 * kap
        for js in range(1, 200):
            s_ = arb(js) / 200
            al = (1 - th_y) * s_ ** th_y
            be = th_y * s_ ** (th_y - 1)
            A = 1 / m_in ** 2 + R * ((1 + e_) * al + kp)
            Bv = R * ((1 + e_) * be + kp)
            v = arb(max(A.upper(), Bv.upper()))
            if beste is None or v < beste[0]:
                beste = (v, e_, s_)
    return beste


bf = ce_fein(m_in, RM, D0, Y1, Y2, kap)
CEf = bf[0].sqrt()
print(f"Verfeinert (Kanten getrennt, ε = {float(bf[1]):.3f}, s = {float(bf[2]):.3f}): C_E' ≤ {float(CEf.upper()):.5f} ⇒ σ_min(J_∞|P_N) ≥ "
      f"{float((c / CEf).lower()):.4f}")
