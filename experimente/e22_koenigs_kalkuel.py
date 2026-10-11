"""exp.koenigs_kalkuel — Koenigs-Spektralreihe, wertseitiger und quasikonformer Basiswechsel (Juli-Synthese, §§7, 15, 17, 19).

 (1) T(x) = a + Σ c_n λ^{nx}, c_n = p_n σ(1)^n; Ableiten ↔ n log λ, Integrieren ↔ 1/(n log λ), Shift-Summe ↔ 1/(1 − λ^n);
     Bruchsumme in Lambert-Form Σ(x;h) = a h + Σ c_n λ^{nx} (1 − λ^{nh})/(1 − λ^n).
 (2) K_{b,d} = T_b ∘ A_d: χ_b ∘ K = χ_d^ρ, ρ = log λ_b / log λ_d; Verzweigung K(y) − a_b ~ C (a_d − y)^ρ;
     Jacobi-Formel K'(T_d(x+h))/K'(T_d(x)) = (ln b/ln d)^h b^{Σ_b(x;h)} / d^{Σ_d(x;h)}.
 (3) Quasikonformer Basiswechsel: K = Ψ_b ∘ exp ∘ M ∘ Log ∘ σ_d, M reell-linear, M(Ω_d) = Ω_b, M(2πi) = 2πi.
 (4) Normalform: H(z) = c z + κ, c = e^κ, konjugiert e^{cz} zu e^w + κ; H(1) = E_κ(κ).
"""
import mpmath as mp

from fixpunkt import AmFixpunkt
from gemeinsam import Lauf, plt
from regulaer import Regulaer

mp.mp.dps = 40
lauf = Lauf("exp.koenigs_kalkuel", __file__)

# ------------------------------------------------------------------ (1) Spektralreihe zur Basis √2
b = mp.sqrt(2)
R = Regulaer(b, N=80, tol=mp.mpf("1e-3"))
a, lam, lb = R.a, R.lam_a, R.lb
F = AmFixpunkt(lb, a, N=80)
w1 = F.sigma(mp.mpf(1))
N = len(F.p) - 1
c = [mp.mpf(0)] + [mp.re(F.p[n] * w1 ** n) for n in range(1, N + 1)]
ll = mp.log(lam)


def spek(x, r=0):
    return (a if r == 0 else 0) + mp.fsum(c[n] * (n * ll) ** r * lam ** (n * x) for n in range(1, N + 1))


def Tstrich(x):
    return 1 / mp.diff(R.slog, R.tet(x))


x0 = mp.mpf(6)
lauf.pruefe("Spektralreihe: a + Σ c_n λ^{nx} = reguläre Tetration (x = 6, 9)",
            max(abs(spek(x) - R.tet(x)) for x in (x0, mp.mpf(9))) < mp.mpf(10) ** -25, mp.nstr(abs(spek(x0) - R.tet(x0)), 3))
lauf.pruefe("Ableiten ↔ Faktor n·log λ: Spektral-T′(6) = 1/slog′(T(6))", abs(spek(x0, 1) - Tstrich(x0)) < mp.mpf(10) ** -15,
            mp.nstr(abs(spek(x0, 1) - Tstrich(x0)), 3))
lauf.pruefe("c₂ = ln b·c₁²/(2(λ − 1)) (Formel der Synthese)", abs(c[2] - lb * c[1] ** 2 / (2 * (lam - 1))) < mp.mpf(10) ** -30)

# Shift-Summe der Ableitung, unabhängig über die Transferformel T'(y+1) = ln b · T(y+1) · T'(y)
y, d, direkt = R.tet(x0), Tstrich(x0), mp.mpf(0)
for _ in range(400):
    direkt += d
    y_neu = mp.power(b, y)
    d, y = lb * y_neu * d, y_neu
lambert_d = mp.fsum(c[n] * n * ll * lam ** (n * x0) / (1 - lam ** n) for n in range(1, N + 1))
lauf.pruefe("Shift-Summe ↔ 1/(1 − λ^n): Σ_k T′(6+k) über Transferformel = Lambert-Reihe", abs(direkt - lambert_d) < mp.mpf(10) ** -15,
            mp.nstr(abs(direkt - lambert_d), 3))

# Bruchsumme: direkt über Iteration gegen Lambert-Form
def sigma_direkt(x, h, M=400):
    s, y1, y2 = a * h, R.tet(x), R.tet(x + h)
    for _ in range(M):
        s += y1 - y2
        y1, y2 = mp.power(b, y1), mp.power(b, y2)
    return s


h = mp.mpf("0.5")
lam_form = a * h + mp.fsum(c[n] * lam ** (n * x0) * (1 - lam ** (n * h)) / (1 - lam ** n) for n in range(1, N + 1))
lauf.pruefe("Bruchsumme in Lambert-Form = direkte Bruchsumme (x = 6, h = 1/2)", abs(lam_form - sigma_direkt(x0, h)) < mp.mpf(10) ** -25,
            mp.nstr(abs(lam_form - sigma_direkt(x0, h)), 3))

# Integral-Schwanz ↔ 1/(n log λ): ∫_6^∞ (T − a) = ∫_6^7 Σ_m (f^m(T(t)) − a) dt
def schwanz_innen(t):
    y, s = R.tet(t), mp.mpf(0)
    for _ in range(400):
        s += y - a
        y = mp.power(b, y)
    return s


mp.mp.dps = 25
integral = mp.quad(schwanz_innen, [x0, x0 + 1])
mp.mp.dps = 40
spek_int = -mp.fsum(c[n] * lam ** (n * x0) / (n * ll) for n in range(1, N + 1))
lauf.pruefe("Integrieren ↔ 1/(n log λ): ∫_6^∞ (T − a) dt = −Σ c_n λ^{6n}/(n log λ)", abs(integral - spek_int) < mp.mpf(10) ** -18,
            mp.nstr(abs(integral - spek_int), 3))

# ------------------------------------------------------------------ (2) wertseitiger Basiswechsel √2 ← 1.3
Rd = Regulaer(mp.mpf("1.3"), N=80, tol=mp.mpf("1e-3"))
Fd = AmFixpunkt(Rd.lb, Rd.a, N=80)
rho = mp.log(lam) / mp.log(Rd.lam_a)


def K(y):
    return R.tet(Rd.slog(y))


def chi(Fx, z):
    return Fx.sigma(z) / Fx.sigma(mp.mpf(1))


res = max(abs(chi(F, K(yy)) - chi(Fd, yy) ** rho) for yy in (mp.mpf("1.05"), mp.mpf("1.2"), mp.mpf("1.35")))
lauf.pruefe(f"wertseitig: χ_√2(K(y)) = χ_1.3(y)^ρ mit ρ = log λ_b/log λ_d = {mp.nstr(rho, 8)}", res < mp.mpf(10) ** -20, mp.nstr(res, 3))
pkt = []
for k in (20, 30):                         # nah am Fixpunkt, damit die Korrektur O(K − a_b) der Zielseite klein ist
    yy = Rd.a - mp.mpf(10) ** -k
    pkt.append((mp.log(mp.mpf(10) ** -k), mp.log(abs(K(yy) - a))))
exponent = (pkt[1][1] - pkt[0][1]) / (pkt[1][0] - pkt[0][0])
lauf.pruefe("Verzweigung am Quellfixpunkt: K(y) − a_b ≈ C·(a_d − y)^ρ", abs(exponent - rho) < mp.mpf(10) ** -6, f"Exponent {mp.nstr(exponent, 10)}")

xj, hj = mp.mpf("0.3"), mp.mpf("0.5")
Td_strich = lambda x: 1 / mp.diff(Rd.slog, Rd.tet(x))
links = mp.log((Tstrich(xj + hj) / Td_strich(xj + hj)) / (Tstrich(xj) / Td_strich(xj)))


def sigma_d(x, h, M=400):
    s, y1, y2 = Rd.a * h, Rd.tet(x), Rd.tet(x + h)
    for _ in range(M):
        s += y1 - y2
        y1, y2 = mp.power(Rd.b, y1), mp.power(Rd.b, y2)
    return s


rechts = hj * mp.log(lb / Rd.lb) + lb * sigma_direkt(xj, hj) - Rd.lb * sigma_d(xj, hj)
lauf.pruefe("Jacobi-Formel: log K′ ändert sich um h·log(ln b/ln d) + ln b·Σ_b − ln d·Σ_d", abs(links - rechts) < mp.mpf(10) ** -15,
            mp.nstr(abs(links - rechts), 3))

# ------------------------------------------------------------------ (3) quasikonformer Basiswechsel (komplexe Basen)
def basis_aus_lam(lmb):
    lmb = mp.mpc(lmb)
    return AmFixpunkt(lmb * mp.exp(-lmb), mp.exp(lmb), N=60)


Bq = basis_aus_lam(0.6 * mp.expj(0.9))
Dq = basis_aus_lam(0.5 * mp.expj(2.0))


def M_abb(v, Om_d, Om_b):
    # v = s·Ω_d + t·2πi (s, t reell) ↦ s·Ω_b + t·2πi
    A = mp.matrix([[mp.re(Om_d), 0], [mp.im(Om_d), 2 * mp.pi]])
    s, t = mp.lu_solve(A, mp.matrix([mp.re(v), mp.im(v)]))
    return s * Om_b + t * 2j * mp.pi


def Kq(z, Bx=Bq, Dx=Dq):
    return Bx.Psi(mp.exp(M_abb(mp.log(Dx.sigma(z)), Dx.log_lam, Bx.log_lam)))


zs = [Dq.L + mp.mpc("0.05", "0.03"), Dq.L + mp.mpc("-0.04", "0.06"), Dq.L + mp.mpc("0.02", "-0.07")]
res = max(abs(Kq(mp.exp(Dq.u * z)) - mp.exp(Bq.u * Kq(z))) for z in zs)
lauf.pruefe("quasikonform: K ∘ f_d = f_b ∘ K (Rückkehrabbildungen konjugiert, Zeit erhalten)", res < mp.mpf(10) ** -25, mp.nstr(res, 3))
eps = mp.mpf(10) ** -12


def dbar(fun, z):
    return (fun(z + eps) - fun(z - eps) + 1j * (fun(z + 1j * eps) - fun(z - 1j * eps))) / (4 * eps)


nicht_holo = abs(dbar(Kq, zs[0]))
Bgleich = basis_aus_lam(0.5 * mp.expj(2.0) * mp.mpf("1.0"))      # gleicher Multiplikator: Ω_b = Ω_d
holo = abs(dbar(lambda z: Kq(z, Bgleich, Dq), zs[0]))
lauf.pruefe("Hindernis: K ist nicht holomorph (∂K/∂z̄ ≠ 0), bei gleichem Multiplikator dagegen holomorph",
            nicht_holo > 1e-3 and holo < mp.mpf(10) ** -10, f"|∂̄K| = {mp.nstr(nicht_holo, 3)} bzw. {mp.nstr(holo, 3)}")

# ------------------------------------------------------------------ (4) Normalform e^w + κ
cc = mp.mpc(0, mp.pi / 2)                                   # Basis i
kappa = mp.log(cc)
H = lambda z: cc * z + kappa
E = lambda w: mp.exp(w) + kappa
res = max(abs(H(mp.exp(cc * z)) - E(H(z))) for z in (mp.mpc(0.3, 0.2), mp.mpc(-1, 2), mp.mpc(2, -0.5)))
lauf.pruefe("Normalform: H(z) = cz + κ konjugiert e^{cz} zu e^w + κ, und H(1) = E_κ(κ) (Turm = singuläre Bahn)",
            res < mp.mpf(10) ** -35 and abs(H(1) - E(kappa)) < mp.mpf(10) ** -35, mp.nstr(res, 3))

lauf.daten = {"sqrt2": {"a": mp.nstr(a, 20), "lambda": mp.nstr(lam, 20), "c": [mp.nstr(c[n], 15) for n in range(1, 7)]},
              "rho_sqrt2_1.3": mp.nstr(rho, 15), "verzweigungsexponent": mp.nstr(exponent, 12),
              "integral_schwanz_6": mp.nstr(integral, 20), "qk_dbar": mp.nstr(nicht_holo, 6)}

fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
xs = [mp.mpf(i) / 10 for i in range(30, 101)]
ax[0].semilogy([float(x) for x in xs], [float(abs(spek(x) - a)) for x in xs], label="|T(x) − a| (Spektralreihe)")
for n in (1, 2, 3):
    ax[0].semilogy([float(x) for x in xs], [float(abs(c[n] * lam ** (n * x))) for x in xs], "--", lw=0.8, label=f"|c_{n} λ^{{{n}x}}|")
ax[0].set_xlabel("x"); ax[0].set_title("Koenigs-Spektrum der √2-Tetration"); ax[0].legend(fontsize=7)
ys = [Rd.a - mp.mpf(10) ** (-mp.mpf(k) / 4) for k in range(4, 60)]
ax[1].loglog([float(Rd.a - yy) for yy in ys], [float(abs(K(yy) - a)) for yy in ys], label="|K(y) − a_√2|")
ax[1].loglog([float(Rd.a - yy) for yy in ys], [float((Rd.a - yy) ** rho) for yy in ys], ":", color="#888", label=f"(a_d − y)^ρ, ρ = {mp.nstr(rho, 4)}")
ax[1].set_xlabel("a_1.3 − y"); ax[1].set_title("Wertseitiger Basiswechsel: Potenz-Verzweigung"); ax[1].legend(fontsize=7)
lauf.bild(fig)
lauf.zusammenfassung = (f"Spektralreihe, Ableitungs-, Integral- und Shift-Kalkül und Lambert-Bruchsumme für √2 bestätigt. "
                        f"K = T_√2∘A_1.3 ist in Koenigs-Koordinaten die Potenz ρ = {mp.nstr(rho, 8)} (Verzweigungsexponent gemessen {mp.nstr(exponent, 8)}); "
                        "Jacobi-Formel über Bruchsummen bestätigt. Quasikonformer Basiswechsel konjugiert, ist aber nur bei gleichem Multiplikator holomorph. "
                        "Normalform e^w + κ bestätigt.")
lauf.schreibe()
