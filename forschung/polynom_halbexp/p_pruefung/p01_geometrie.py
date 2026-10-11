"""Strang P, Prüfung 1: Geometrie von W = D ∪ D(a0,R0), Lage von f_K(D) relativ zum A-Bogen,
explizite Green-Funktion (B Lemma 7.3 / C Satz 12.2), Konstanten alpha, cap, g am Schlitz.

Nur Daten aus gemeinsam/ (Kneser-Koeffizienten 4000), kein Motor.
"""
import json
import sys
from pathlib import Path

import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 40
HIER = Path(__file__).resolve().parent
GEM = HIER.parent / "gemeinsam"

a = [mp.mpf(s) for s in json.loads((GEM / "kneser_halbexp_koeff_4000.json").read_text())["a"]]
da = [k * a[k] for k in range(1, len(a))]


def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s


def f(z):
    return horner(a, z)


def fp(z):
    return horner(da, z)


# Fixpunkt L
L = mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))
rho = abs(L)
a0 = a[0]
R0 = abs(L - a0)
print("L =", mp.nstr(L, 20), " rho =", mp.nstr(rho, 17), " a0 =", mp.nstr(a0, 17), " R0 =", mp.nstr(R0, 17))
print("Kontrolle f(L) - L (Reihe 4000 Terme, auf dem Rand):", mp.nstr(abs(f(L) - L), 5))
mu = fp(L)
print("f'(L) =", mp.nstr(mu, 12), " sqrt(L) = e^{L/2} =", mp.nstr(mp.exp(L / 2), 12), " |mu| =", mp.nstr(abs(mu), 8),
      " arg mu (Grad) =", mp.nstr(mp.degrees(mp.arg(mu)), 8))

# Winkel an der Ecke
argL = mp.degrees(mp.arg(L))
argLa = mp.degrees(mp.arg(L - a0))
aussen = 180 - (argLa - argL)
print(f"arg L = {mp.nstr(argL, 10)}°, arg(L-a0) = {mp.nstr(argLa, 10)}°, Aussenwinkel = {mp.nstr(aussen, 10)}°, "
      f"alpha = {mp.nstr(aussen / 180, 10)}")
# Aussenkeil: Richtungen ausserhalb D: (argL-90, argL+90); ausserhalb D(a0): (argLa-90, argLa+90)
print(f"Aussenkeil von W bei L: ({mp.nstr(argLa - 90, 6)}°, {mp.nstr(argL + 90, 6)}°)")
argmu = mp.degrees(mp.arg(mu))
print(f"f(D) nahe L: Richtungen ({mp.nstr(argL + 90 + argmu, 6)}°, {mp.nstr(argL + 270 + argmu, 6)}°) "
      f"[= D-Halbscheibe gedreht um arg mu]")
print(f"Aussenkeil von W ∪ f(D) bei L: ({mp.nstr(argL - 90 + argmu, 6)}°, {mp.nstr(argL + 90, 6)}°), "
      f"Oeffnung {mp.nstr(180 - argmu, 6)}°, alpha' = {mp.nstr((180 - argmu) / 180, 6)}")

# ---- Green-Funktion von C \ W-bar (Moebius + Potenz)
Lb = mp.conj(L)


def M(z):
    return (z - L) / (z - Lb)


# Strahlen: Bild von dD und dD(a0)
th1 = mp.arg(M(-rho))           # Bild des linken D-Bogens
th2 = mp.arg(M(a0 + R0))        # Bild des aeusseren A-Bogens
print("arg M(-rho) =", mp.nstr(mp.degrees(th1), 10), "  arg M(a0+R0) =", mp.nstr(mp.degrees(th2), 10))
# Aussengebiet: Keil zwischen th2 und th1 (Richtung pruefen mit Punkt 3.0)
t_inf = mp.arg(M(mp.mpf(10) ** 6))
print("arg M(gross) =", mp.nstr(mp.degrees(t_inf), 10))
beta = th2
omega = th1 - th2
if omega < 0:
    omega += 2 * mp.pi
print("Keil: beta =", mp.nstr(mp.degrees(beta), 10), " omega =", mp.nstr(mp.degrees(omega), 12))


def Wmap(z):
    w = M(z) * mp.exp(-1j * beta)
    ang = mp.arg(w)
    if ang < 0:
        ang += 2 * mp.pi
    return mp.power(abs(w), mp.pi / omega) * mp.exp(1j * ang * mp.pi / omega)


W_inf = Wmap(mp.mpf(10) ** 30)


def g(z):
    w = Wmap(z)
    return mp.log(abs((w - mp.conj(W_inf)) / (w - W_inf)))


print("W_inf =", mp.nstr(W_inf, 12))
# Probe: g = 0 auf beiden Boegen
r1 = max(abs(g(rho * mp.expjpi(t))) for t in [0.6, 0.8, 1.0, 1.2, 1.4])
r2 = max(abs(g(a0 + R0 * mp.expj(t))) for t in [-1.0, -0.5, 0, 0.5, 1.0])
print("g auf D-Bogen (max):", mp.nstr(r1, 3), " auf A-Bogen:", mp.nstr(r2, 3))
Z = mp.mpf(10) ** 12
robin = g(Z) - mp.log(Z)
print("g - log|z| bei grossem z:", mp.nstr(robin, 10), " cap =", mp.nstr(mp.exp(-robin), 10))
# Harmonizitaet (diskreter Laplace) an zwei Punkten
for z0 in [mp.mpc(2.2, 0), mp.mpc(-1, 1.5)]:
    hh = mp.mpf("1e-4")
    lap = (g(z0 + hh) + g(z0 - hh) + g(z0 + 1j * hh) + g(z0 - 1j * hh) - 4 * g(z0)) / hh ** 2
    print("Laplace g bei", mp.nstr(z0, 4), ":", mp.nstr(lap, 3))
# g am Schlitz
bis = mp.radians((mp.degrees(beta) + 0))  # unused
psi_b = (argLa - 90 + argL + 90) / 2
print("Winkelhalbierende psi =", mp.nstr(psi_b, 8), "°")
for u in ["1e-4", "1e-3", "1e-2", "0.05", "0.1", "0.2"]:
    u = mp.mpf(u)
    gb = g(L + u * mp.expj(mp.radians(psi_b)))
    gr = g(L * (1 + u))
    print(f"u = {mp.nstr(u, 3):6s}  g(L+u e^(i psi))/u^(1/alpha) = {mp.nstr(gb / u ** (180 / aussen), 6)}   "
          f"g(L(1+u))/u^(1/alpha) = {mp.nstr(gr / u ** (180 / aussen), 6)}")
print("g(2.2) =", mp.nstr(g(mp.mpf(2.2)), 8), " e^-g =", mp.nstr(mp.exp(-g(mp.mpf(2.2))), 6),
      "  Max-Formel:", mp.nstr(max(rho / 2.2, R0 / (2.2 - a0)), 6))
# u0(d): g(L(1+u0)) = d
for d in [0.05, 0.1]:
    u0 = mp.findroot(lambda u: g(L * (1 + u)) - d, 0.1)
    print(f"u0({d}) = {mp.nstr(u0, 6)}  (|zeta-L| = {mp.nstr(u0 * rho, 5)})")

# ---- Lage des A-Bogens relativ zu f(D): Urbild h(eta) per Newton, Fortsetzung von L aus
print("\nA-Bogen eta = a0 + R0 e^{i psi}: Urbild h(eta) unter f (Fortsetzung von L aus)")
psiL = mp.arg(L - a0)
x = L
out = []
for j in range(1, 41):
    psi = psiL - j * psiL / 40
    eta = a0 + R0 * mp.expj(psi)
    # Newton
    for _ in range(60):
        dx = (f(x) - eta) / fp(x)
        x -= dx
        if abs(dx) < mp.mpf(10) ** -30:
            break
    out.append((mp.degrees(psi), abs(x) / rho, mp.degrees(mp.arg(x - L)), abs(eta - L)))
for (p, r, ang, dist) in out[::4] + [out[-1]]:
    print(f"  psi = {mp.nstr(p, 6):>8}°  |eta-L| = {mp.nstr(dist, 4):>8}  |h(eta)|/rho = {mp.nstr(r, 10)}")
print("  => A-Bogen liegt in f(D) genau dort, wo |h|/rho < 1.")

# ---- f(Gamma_out) gegen A-Bogen: radialer Abstand
print("\nf(Gamma_out) = Bild des rechten Kreisbogens: Lage relativ zu D(a0,R0)")
for phi_deg in [76, 70, 60, 45, 30, 15, 0]:
    z = rho * mp.expj(mp.radians(phi_deg))
    w = f(z)
    print(f"  arg x = {phi_deg:3d}°: f(x) = {mp.nstr(w, 8):>28}  |f(x)| / rho = {mp.nstr(abs(w) / rho, 6)}  "
          f"|f(x)-a0|/R0 = {mp.nstr(abs(w - a0) / R0, 8)}")

# ---- Ist der D-Bogen ausserhalb f(D-bar)? Windungszahl von f(dD) um eta
print("\nWindungszahl von f(dD) um Punkte des D-Bogens und des A-Bogens:")
K = 4000
pts = [rho * mp.expj(2 * mp.pi * k / K) for k in range(K)]
img = [f(z) for z in pts]


def wind(eta):
    s = mp.mpf(0)
    for k in range(K):
        w1 = img[k] - eta
        w2 = img[(k + 1) % K] - eta
        s += mp.arg(w2 / w1)
    return s / (2 * mp.pi)


for name, eta in [("D-Bogen 90°", rho * mp.expj(mp.radians(90))), ("D-Bogen 120°", rho * mp.expj(mp.radians(120))),
                  ("D-Bogen 180°", -rho), ("A-Bogen psi=60°", a0 + R0 * mp.expj(mp.radians(60))),
                  ("A-Bogen psi=0°", a0 + R0), ("A-Bogen psi=90°", a0 + R0 * mp.expj(mp.radians(90)))]:
    print(f"  {name:18s} Windungszahl = {mp.nstr(wind(eta), 6)}")
