"""O09b — Satz O2 allgemein (M_a C_φ, φ Möbius, H²(D(0,r))) und Anwendung auf das Dualbild bei a₀ für Möbius-exp:
K' = K^{-1} = M_{f'} C_h auf H²(D(a₀,R₀)); verschoben Y = y − a₀: a(Y) = f'(Y+a₀), φ(Y) = h(Y+a₀) − a₀.
Kriterium: Δ₁(Y,Z̄) = r²(γY+δ) − (αY+β)Z̄ = e^{iκ}·hermitesch;  b = r² e^{−iκ}(γY+δ) a(Y);  Re b ≥ 0 auf D_H ⟺ akkretiv."""
import sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 40

def kriterium(al, be, ga, de, afun, r, name):
    r2 = r * r
    # Δ₁ = r²δ + r²γ Y − β Z̄ − α Y Z̄ ; hermitesch bis auf Phase: wähle κ = arg(r²δ)
    kap = mp.arg(r2 * de)
    e = mp.expj(-kap)
    c0, cy, cz, cyz = e * r2 * de, e * r2 * ga, -e * be, -e * al
    herm = abs(mp.im(c0)) + abs(mp.im(cyz)) + abs(cy - mp.conj(cz))
    # D_H = {c0 + cy Y + conj(cy) Ȳ + cyz |Y|² > 0}, cyz reell
    A_ = mp.re(cyz); C0 = mp.re(c0)
    # c0 + 2Re(cy Y) + A_|Y|² > 0
    if A_ < 0:   # Kreisscheibe |Y − Yc|² < R²
        Yc = -mp.conj(cy) / A_
        R2 = abs(Yc) ** 2 - C0 / A_
        kind = "Scheibe"
    else:
        Yc = -mp.conj(cy) / A_
        R2 = abs(Yc) ** 2 - C0 / A_
        kind = "Äusseres"
    R = mp.sqrt(R2)
    b = lambda Y: r2 * e * (ga * Y + de) * afun(Y)
    ws = [Yc + R * mp.expj(p) for p in mp.linspace(0, 2 * mp.pi, 40001)]
    vals = [mp.re(b(w)) for w in ws]
    i = min(range(len(vals)), key=lambda k: vals[k])
    Yp = -de / ga
    print(f"{name}: hermitesch-Defekt {mp.nstr(herm,3)}; D_H = {kind} Mitte {mp.nstr(Yc,8)} Radius {mp.nstr(R,8)}; Pol Möbius {mp.nstr(Yp,8)} "
          f"(Abstand {mp.nstr(abs(Yp-Yc),6)}); min Re b auf ∂D_H = {mp.nstr(vals[i],8)} bei {mp.nstr(ws[i],6)}; Re b(0) = {mp.nstr(mp.re(b(0)),6)}")

L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
A, B = L.real, L.imag
rho2 = A * A + B * B
c = mp.atan((1 - A) / B) - mp.atan(-A / B)
h = mp.tan(c / 2)
al, be, ga, de = B - A * h, h * rho2, -h, B + h * A
det = al * de - be * ga
f = lambda y: (al * y + be) / (ga * y + de)
fp = lambda y: det / (ga * y + de) ** 2
m = lambda y: fp(f(y))
# Bild bei 0: K = M_{1/m} C_f auf H²(D_ρ)
kriterium(al, be, ga, de, lambda Y: 1 / m(Y), mp.sqrt(rho2), "Bild 0 (K)")
# Dualbild: a0 = f(0), R0 = |L − a0|; h = f^{-1} = (δ y − β)/(−γ y + α); φ(Y) = h(Y + a0) − a0
a0 = f(0); R0 = abs(L - a0)
# h(Y+a0) − a0 = (δ(Y+a0) − β − a0(−γ(Y+a0) + α)) / (−γ(Y+a0) + α)
pa = de + a0 * ga; pb = de * a0 - be - a0 * (-ga * a0 + al); pc = -ga; pd = -ga * a0 + al
kriterium(pa, pb, pc, pd, lambda Y: fp(Y + a0), R0, "Dualbild a₀ (K⁻¹)")
