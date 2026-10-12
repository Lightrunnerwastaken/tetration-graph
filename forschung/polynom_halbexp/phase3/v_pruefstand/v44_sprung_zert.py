"""Strang V — Sprungschranke |Δf_K(L(1+u))| ≤ C_J·u^σ STRENG aus Ps oberer Struktur (Satz P-V, p32a, kc_rig), u ∈ (0, U_MAX].

Grundlage (P, bewiesen): sexp_K(z) = P(e^{L(z + θ_K(z))}) für Im z ≥ y₁ = 0.36 (alle Re z), θ_K 1-periodisch, |θ_K − θ̃| ≤ δθ = 4.39e-167
(p32a), θ̃ = Σ c̃_n qⁿ (241 Moden, p26_iter6), P Poincaré-Funktion an L (kc_rig: Koeffizienten + Restkonstante, Bälle).
Herleitung: x nahe L, W := P⁻¹(x) (der kleine Zweig, Krawczyk). Jedes z mit z + θ_K(z) = (ln|W| + iA)/L, A ≡ arg W (mod 2π), ist ein Urbild
von x unter sexp_K. Für x ∈ D nahe L ist slog_K(x) eines davon: A ∈ (θ + π/2, θ + 3π/2) + 2πk₀ — k₀ wird unten an x₁ = 1.2e^{iθ} über Ps
slog-Krawczyk (dieselbe Verzweigung, mit der P-V die f_K-Werte auf |x| = 1.2 erzeugt) bestimmt. Dann f_K(x) = sexp_K(z + ½) =
P(e^{L(z + ½ + θ_K(z+½))}) = P(√L·W·e^{Lδ(z)}), δ(z) := θ_K(z+½) − θ_K(z) = −2Σ_{n ungerade} κ_n qⁿ, √L := e^{L/2}.
Die beiden Ufer des Schlitzes S_L = {L(1+u)} (Fortsetzung von D aus über die Γin- bzw. die Γout-Seite) haben A_in = θ + 2πk₀ bzw.
A_out = θ + 2π(k₀+1), also Δf = P(√L W e^{Lδ_in}) − P(√L W e^{Lδ_out}) und (Mittelwertsatz, gemeinsames W)
    |Δf| ≤ sup_H |P'| · |√L W| · |e^{Lδ_in} − e^{Lδ_out}|,   H = √L·W·(1 + D(0, max|e^{Lδ} − 1|)).
Teil 1 (u ≤ u₁, analytisch): |W| ≤ ρu/(1−κ₂), |arg W − θ| ≤ asin κ₂ (Rouché auf D(0, r_s)), |q| ≤ |e^{γ log W}|·e^{2π(Im c̃₀ + δθ + S₁q_max)},
    γ = 2πi/L, |δ| ≤ 2S_ung|q| ⇒ |Δf| ≤ K₀·u^{1+Re γ}, also |Δf| ≤ K₀ u^{SIG} für SIG ≤ 1 + Re γ = 5.44694… (u ≤ 1).
Teil 2 (u₁ ≤ u ≤ U_MAX): Intervalle [u_a, u_b], alles in Bällen (W per Krawczyk, z per Brouwer-Box mit Kontraktion, Im z ≥ y₁ geprüft).
Aufruf: python v44_sprung_zert.py [U_MAX] [u1]
Ausgabe: v44_sprung.json (Tabelle u_a, u_b, Schranke |Δf|, Quotient |Δf|/u_a^SIG; C_J(u_max) kumulativ).
"""
import json
import math
import sys
import time
from pathlib import Path

from flint import acb, arb, ctx

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER / "pv_kopie"))
import kc_rig  # noqa: E402
from kc_kneser import Kneser  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
U_MAX = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
U1 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
PREC = 300
t0 = time.time()
KN = Kneser(PREC)
R = KN.R
ctx.prec = PREC
L = R.L
c = KN.c
pi = arb.pi()
I = acb(0, 1)
TPI = 2 * pi * I
dth = arb(json.loads((HIER / "pv_kopie" / "p32a_theta_K.json").read_text(encoding="utf-8"))["delta_theta"])
dth = arb(dth.upper())
Y1 = arb("0.36")
rho = abs(L)
th = L.imag                    # arg L = Im L, da Log L = L
gam = TPI / L
SIG_WAHR = 1 + gam.real
SIG = arb("5.4469")
assert SIG < SIG_WAHR
sqL = (L / 2).exp()
cob = [x.abs_upper() for x in c]
print(f"L = {L.str(12)}, ρ = {rho.str(10)}, γ = 2πi/L = {gam.str(10)}, σ = 1 + Re γ = {SIG_WAHR.str(10)} (verwendet SIG = 5.4469), δθ = {dth.str(3)}",
      flush=True)


def ob(x):
    return x.abs_upper()


def ball(z, r):
    return kc_rig.disc(z, r)


def theta_box(z):
    q = (TPI * z).exp()
    return R._horner(c, q) + acb(arb(0, dth), arb(0, dth))


def dtheta_box(z):
    """θ̃'(z) auf der Box (für die Kontraktion)."""
    q = (TPI * z).exp()
    dc = [n * cn for n, cn in enumerate(c)]
    return R._horner(dc, q) * TPI


codd = [c[n] if n % 2 == 1 else acb(0) for n in range(len(c))]


def delta_box(z):
    """δ(z) = θ_K(z+½) − θ_K(z) = −2Σ_{ungerade} κ_n qⁿ; Fehler θ_K − θ̃ an beiden Stellen ≤ δθ."""
    q = (TPI * z).exp()
    return -2 * R._horner(codd, q) + acb(arb(0, 2 * dth), arb(0, 2 * dth))


def P_inv(X, w0):
    """Box ⊇ {W : P(W) = x, x ∈ X} mit Eindeutigkeit (Krawczyk). w0: Startwert (acb, Mitte)."""
    w = w0
    Xm = acb(X.real.mid(), X.imag.mid())
    for _ in range(80):
        f, df = R.P_dP_ball(w)
        dw = (f - Xm) / df
        dw = acb(dw.real.mid(), dw.imag.mid())
        w = acb((w - dw).real.mid(), (w - dw).imag.mid())
        if ob(dw) < arb(2) ** (-PREC + 30):
            break
    f0, d0 = R.P_dP_ball(w)
    Y = 1 / acb(d0.real.mid(), d0.imag.mid())
    r = 2 * ob(f0 - X) / abs(d0).lower() + arb(10) ** -80
    for _ in range(12):
        B = ball(w, r)
        _, dB = R.P_dP_ball(B)
        Kc = w - Y * (f0 - X)
        Kr = ob(1 - Y * dB) * r
        if ob(Kc - w) + Kr < r:
            return w, r
        r = r * 2
    raise RuntimeError("Krawczyk für P⁻¹ versagt")


def z_fix(zeta):
    """Box ⊇ Fixpunkt von z ↦ ζ − θ_K(z) (für alle ζ im Ball); Brouwer + Kontraktion; prüft Im z ≥ y₁."""
    zm = acb(zeta.real.mid(), zeta.imag.mid()) - c[0]
    for _ in range(200):
        zn = acb(zeta.real.mid(), zeta.imag.mid()) - R._horner(c, (TPI * zm).exp())
        zn = acb(zn.real.mid(), zn.imag.mid())
        if ob(zn - zm) < arb(2) ** (-PREC + 30):
            zm = zn
            break
        zm = zn
    r = 2 * (ob(zeta - acb(zeta.real.mid(), zeta.imag.mid())) + ob(zeta - theta_box(zm) - zm)) + arb(10) ** -80
    for _ in range(12):
        B = ball(zm, r)
        assert B.imag.lower() > Y1 + arb("0.005"), f"Im z = {B.imag.str(5)} < y₁"
        T = zeta - theta_box(B)
        kontr = ob(dtheta_box(B)) + dth / arb("0.005")       # |θ_K'| ≤ |θ̃'| + δθ/Abstand (Cauchy)
        if ob(T - zm) < r and kontr < 1:
            return B, kontr
        r = r * 2
    raise RuntimeError("Fixpunkt-Box versagt")


def A_bereich(Wc, A_ziel):
    """arg W (Hauptwert) + 2πm so, dass nahe A_ziel."""
    a = Wc.arg()
    m = round(float(((A_ziel - a) / (2 * pi)).mid()))
    return m


# ---------------------------------------------------------------- (0) Verzweigung k₀ an x₁ = 1.2 e^{iθ} (Ps slog-Krawczyk, Fortsetzung längs |x| = 1.2)
x1 = arb("1.2") * acb(th.cos(), th.sin())
w = KN.slog_num(acb(arb("1.2")), acb(0))
nst = 400
for j in range(1, nst + 1):
    phi = th * j / nst
    w = KN.slog_num(arb("1.2") * acb(phi.cos(), phi.sin()), w)
ws, rs = KN.slog_krawczyk(x1, w, r0=arb("1e-30"))
assert ws is not None, "slog-Krawczyk an x₁ versagt"
W1, rW1 = P_inv(acb(x1), x1 - L)
treffer = []
for k0 in range(-1, 5):
    m = A_bereich(W1, th + pi + 2 * pi * k0)
    lw = ball(W1, rW1).log() + TPI * m
    try:
        Bz, _ = z_fix(lw / L)
    except (AssertionError, RuntimeError) as e:
        print(f"  k₀ = {k0}: {e}")
        continue
    d = ob(Bz - ws)
    enth = bool(d + 0 < rs)
    treffer.append((k0, enth))
    print(f"  k₀ = {k0}: z = {Bz.str(8)}, |z − slog(x₁)| ≤ {d.str(3)}  ⇒ {'IN der Krawczyk-Scheibe (Radius ' + rs.str(2) + ')' if enth else 'nicht'}", flush=True)
K0 = [k for k, e in treffer if e]
assert len(K0) == 1, f"Verzweigung nicht eindeutig: {treffer}"
k0 = K0[0]
print(f"(0) slog_K(x₁) = {ws.str(10)} (Krawczyk-Radius {rs.str(2)}) = z_{{k₀}}(x₁) genau für k₀ = {k0}; Ufer A_in = θ + {2 * k0}π, A_out = θ + {2 * k0 + 2}π "
      f"({time.time() - t0:.0f}s)", flush=True)
A_in = th + 2 * pi * k0
A_out = th + 2 * pi * (k0 + 1)

# ---------------------------------------------------------------- (1) analytisch für u ≤ u₁
rs_ = arb(__import__("os").environ.get("V44_RS", "0.1"))
kap2 = sum((ob(R.p[n]) * rs_ ** (n - 1) for n in range(2, R.n0 + 1)), arb(0)) + R.C * rs_ ** R.n0
rp = arb(__import__("os").environ.get("V44_RP", "0.125"))
kap1 = sum((n * ob(R.p[n]) * rp ** (n - 1) for n in range(2, R.n0 + 1)), arb(0)) + R.C * (2 * rp) ** (R.n0 + 1) / rp
assert kap2 < 1 and kap1 < 1
u1 = arb(U1)
assert rho * u1 < rs_ * (1 - kap2), "u₁ zu gross für Rouché auf D(0, r_s)"
Wmax = rho * u1 / (1 - kap2)
a_arg = (kap2).asin()
# Im ζ ≥ ((A − a)·Re L − ln|W|·Im L)/ρ² am ungünstigsten Ufer (A_in) und |W| = Wmax
imzeta_min = ((A_in - a_arg) * L.real - Wmax.log() * L.imag) / rho ** 2
Y = arb("0.45")
qmax = (-2 * pi * Y).exp()
S1 = sum((ob(c[n]) * qmax ** (n - 1) for n in range(1, len(c))), arb(0))
S1 += sum((dth * (2 * pi * n * Y1).exp() * qmax ** (n - 1) for n in range(1, 2000)), arb(0))
S1 += dth * (2 * pi * 2000 * Y1).exp() * qmax ** 1999 / (1 - (2 * pi * Y1).exp() * qmax)
Sodd = sum((ob(c[n]) * qmax ** (n - 1) for n in range(1, len(c), 2)), arb(0)) + (S1 - sum((ob(c[n]) * qmax ** (n - 1) for n in range(1, len(c))), arb(0)))
assert imzeta_min - (c[0].imag + dth) - S1 * qmax > Y, "Halbebene Im z ≥ Y nicht invariant"
kontr = 2 * pi * sum((n * (ob(c[n]) + dth * (2 * pi * n * Y1).exp()) * qmax ** n for n in range(1, len(c))), arb(0))
assert kontr < 1
fq = (2 * pi * (c[0].imag + dth + S1 * qmax)).exp()


def ezeta(A, Wm):
    return (gam.real * Wm.log() - gam.imag * A).exp()


qin = ezeta(A_in - a_arg, Wmax) * fq
qout = ezeta(A_out - a_arg, Wmax) * fq
assert qin < qmax and qout < qmax
dmax = 2 * Sodd * qin
rH = sqL.abs_upper() * Wmax * (rho * dmax).exp()
assert rH < rp
Mp = 1 + kap1
# K₀: |Δf| ≤ Mp·√ρ·|W|·ρ(|δ_in| + |δ_out|)e^{ρ dmax}, |W| ≤ ρu/(1−κ₂), |δ_b| ≤ 2S_ung·fq·(ρu/(1−κ₂))^{Re γ}·e^{−Im γ(A_b − a)}
K0c = Mp * rho.sqrt() * (rho / (1 - kap2)) * rho * 2 * Sodd * fq * (rho / (1 - kap2)) ** gam.real \
    * ((-gam.imag * (A_in - a_arg)).exp() + (-gam.imag * (A_out - a_arg)).exp()) * (rho * dmax).exp()
K0c = arb(K0c.upper())
print(f"(1) u ≤ u₁ = {U1}: κ₂ ≤ {kap2.str(3)} (D(0,{rs_.str(2)})), κ₁ ≤ {kap1.str(3)} (D(0,{rp.str(3)})), Im ζ ≥ {imzeta_min.str(4)}, Halbebene Im z ≥ {Y.str(3)} "
      f"invariant, Kontraktion {kontr.str(3)}, |q_in| ≤ {qin.str(3)}, |q_out| ≤ {qout.str(3)}  ⇒  |Δf| ≤ K₀·u^{{1+Re γ}}, K₀ = {K0c.str(5)}", flush=True)

# ---------------------------------------------------------------- (2) Intervalle u₁ ≤ u ≤ U_MAX
tab = []
ua = U1
wprev = None
fak = float(__import__("os").environ.get("V44_FAK", "1.01"))
CJ_kum = arb(K0c)
Cmax = {}
while ua < U_MAX - 1e-12:
    ub = min(ua * fak, U_MAX)
    U = arb((ua + ub) / 2, (ub - ua) / 2 * 1.0000001)
    X = L * (1 + U)
    w0 = wprev if wprev is not None else L * arb((ua + ub) / 2)
    try:
        Wc, rW = P_inv(X, w0)
        wprev = Wc
        Wb = ball(Wc, rW)
        erg = {}
        for name, Az in (("in", A_in), ("out", A_out)):
            m = A_bereich(Wc, Az)
            lw = Wb.log() + TPI * m
            assert abs(float((lw.imag - Az).mid())) < 1.0
            Bz, kz = z_fix(lw / L)
            erg[name] = (Bz, delta_box(Bz), kz)
    except (AssertionError, RuntimeError) as e:
        print(f"  Abbruch bei u = {ua:.4f}: {e}", flush=True)
        break
    din, dout = erg["in"][1], erg["out"][1]
    ein, eout = (L * din).exp(), (L * dout).exp()
    mm = max(ob(ein - 1), ob(eout - 1))
    H = sqL * Wb * (1 + acb(arb(0, mm), arb(0, mm)))
    _, dP = R.P_dP_ball(H)
    sprung = ob(dP) * ob(sqL * Wb) * ob(ein - eout)
    quot = sprung / arb(ua) ** SIG
    CJ_kum = arb(max(CJ_kum.upper(), quot.upper()))
    tab.append(dict(u_a=ua, u_b=ub, sprung=float(sprung.upper()), quot=float(quot.upper()), CJ_bis_ub=float(CJ_kum.upper()),
                    imz_in=float(erg["in"][0].imag.lower()), imz_out=float(erg["out"][0].imag.lower()), absW=float(ob(Wb))))
    ua = ub
for umax in (0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.451, 0.46, 0.48, 0.5):
    if umax <= U_MAX + 1e-12:
        werte = [K0c.upper()] + [t["quot"] for t in tab if t["u_a"] < umax - 1e-12]
        Cmax[str(umax)] = float(max(werte))
print("(2) Intervalle (Faktor 1.01): " + ", ".join(f"u = {t['u_a']:.3f}: |Δf| ≤ {t['sprung']:.3e}, /u^σ ≤ {t['quot']:.4f}, Im z_in ≥ {t['imz_in']:.3f}"
                                                     for t in tab[::max(1, len(tab) // 8)]), flush=True)
print("⇒ C_J(u ≤ u_max) = " + ", ".join(f"{k}: {v:.4f}" for k, v in Cmax.items()) + f"   (SIG = 5.4469, streng)  ({time.time() - t0:.0f}s)", flush=True)
json.dump(dict(beschreibung="Sprung |Δf_K(L(1+u))| ≤ C_J u^SIG, SIG = 5.4469 (≤ 1 + Re 2πi/L), streng aus P-V/p32a/kc_rig; Verzweigung k₀ über Ps slog-Krawczyk an 1.2e^{iθ}",
               k0=k0, SIG=5.4469, sigma_wahr=float(SIG_WAHR.mid()), K0_u_klein=float(K0c.upper()), u1=U1, C_J=Cmax, tabelle=tab),
          open(HIER / "v44_sprung.json", "w", encoding="utf-8"), ensure_ascii=False, indent=0)
