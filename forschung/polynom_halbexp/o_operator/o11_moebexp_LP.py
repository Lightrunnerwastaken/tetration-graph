"""O11 — Satz O5 (unbedingte lineare Konvergenz im Möbius-exp-Modell): Geometrie K' = D̄ ∪ f(D̄), Bedingung ζ_p ∉ K' (⇔ |h(ζ_p)| > ρ),
Green-Funktion von ℂ∖K' (Sektor in M = (x−L)/(x−L̄)), Rate e^{−g_{K'}(ζ_p)} gegen die wahre Rate e^{−g_W(ζ_p)} = 0.4092 (Strang E),
sowie max_{∂D}|1/m| (= Schranke für ‖J_N^{-1}‖ nach Satz O2 + D8)."""
import sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 30
L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
Lb = mp.conj(L); rho = abs(L)
A, B = L.real, L.imag
c = mp.atan((1 - A) / B) - mp.atan(-A / B)
h_ = mp.tan(c / 2)
al, be, ga, de = B - A * h_, h_ * rho ** 2, -h_, B + h_ * A
det = al * de - be * ga
f = lambda y: (al * y + be) / (ga * y + de)
hinv = lambda y: (de * y - be) / (-ga * y + al)
fp = lambda y: det / (ga * y + de) ** 2
m = lambda y: fp(f(y))
zp = -de / ga
print("ζ_p =", mp.nstr(zp, 10), " h(ζ_p) =", mp.nstr(hinv(zp), 10), " |h(ζ_p)|/ρ =", mp.nstr(abs(hinv(zp)) / rho, 8))
M = lambda x: (x - L) / (x - Lb)
lam = M(f(mp.mpf(0.3))) / M(mp.mpf(0.3))
print("λ = M∘f/M =", mp.nstr(lam, 10), " |λ| =", mp.nstr(abs(lam), 10), " arg λ =", mp.nstr(mp.arg(lam), 10), " (c =", mp.nstr(c, 10), ")")
# Bild von ∂D: Gerade durch 0; Richtung aus M(ρ e^{iφ})
d1 = mp.arg(M(rho * mp.expj(0.3)))
d2 = mp.arg(M(rho * mp.expj(3.0)))
w0 = M(0)
print("Randrichtungen von M(∂D):", mp.nstr(d1, 8), mp.nstr(d2, 8), "; M(0) =", mp.nstr(w0, 8), " arg", mp.nstr(mp.arg(w0), 8), "; M(∞) = 1")
# H = Halbebene mit w0; Grenzwinkel θa < arg < θa+π
th = [d1, d2]
# bestimme θa so, dass arg w0 in (θa, θa+π)
for t in th + [t + mp.pi for t in th] + [t - mp.pi for t in th]:
    if 0 < mp.arg(w0 * mp.expj(-t)) < mp.pi:
        tha = t; break
print("H = {", mp.nstr(tha, 8), "< arg w <", mp.nstr(tha + mp.pi, 8), "}")
phi_l = mp.arg(lam)
# K' = H ∪ λH: Sektor θa ≤ arg ≤ θa + π + arg λ (falls arg λ > 0) ; Komplement: θa+π+argλ < arg < θa + 2π
lo, hi = tha + mp.pi + phi_l, tha + 2 * mp.pi
om = hi - lo
def g_sektor(w, lo, hi):
    om = hi - lo
    a = mp.arg(w * mp.expj(-lo)) % (2 * mp.pi)
    v = abs(w) ** (mp.pi / om) * mp.expj(a * mp.pi / om)
    a1 = mp.arg(mp.mpc(1) * mp.expj(-lo)) % (2 * mp.pi)
    v1 = mp.expj(a1 * mp.pi / om)
    return mp.log(abs((v - mp.conj(v1)) / (v - v1)))
wp = M(zp)
gK = g_sektor(wp, lo, hi)
# W-Geometrie: D(a0,R0) ↔ λ^{1/2}H
loW = tha + mp.pi + phi_l / 2
gW = g_sektor(wp, loW, hi)
print("Komplement von K' (Öffnung", mp.nstr(om, 8), "): e^{−g_K'(ζ_p)} =", mp.nstr(mp.exp(-gK), 8),
      "; W-Geometrie (Öffnung", mp.nstr(hi - loW, 8), "): e^{−g_W(ζ_p)} =", mp.nstr(mp.exp(-gW), 8))
print("max_{∂D}|1/m| =", mp.nstr(max(abs(1 / m(rho * mp.expj(t))) for t in mp.linspace(0, 2 * mp.pi, 4001)), 8),
      "; max_{|x|≤ρ} Pole von m: |h(ζ_p)| =", mp.nstr(abs(hinv(zp)), 8))
