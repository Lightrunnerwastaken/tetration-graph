"""s5: exakte Monodromie bei L:  slog = u + theta(u), u = regulaerer slog (Schroeder) bei L, theta 1-periodisch
(Fourier-Koeffizienten kappa_m nur lesend aus polynom_halbexp/gemeinsam/kneser_singularitaet_L.json).
Umlauf um L:  u -> u + gamma, also  M_L(s) = (u+gamma) + theta(u+gamma),  u = (id+theta)^(-1)(s).
Test: linker Wert (P-Blatt) an Re w = Re L - 0.001 + iy  ->  vorhergesagter rechter Wert (K-Blatt) gegen Motor (s4)."""
import json
import mpmath as mp

mp.mp.dps = 40
D = json.load(open(r"D:\Users\janis\Documents\Matura\tetration-graph\forschung\polynom_halbexp\gemeinsam\kneser_singularitaet_L.json"))
kap = {int(m): mp.mpc(mp.mpf(v[0]), mp.mpf(v[1])) for m, v in D["kappa"].items()}
L = mp.conj(-mp.lambertw(-1))
g = 2j * mp.pi / L


def th(u):
    return sum(c * mp.exp(2j * mp.pi * m * u) for m, c in kap.items())


def dth(u):
    return sum(c * 2j * mp.pi * m * mp.exp(2j * mp.pi * m * u) for m, c in kap.items())


def u_von_s(s):
    u = s - kap[0]
    for _ in range(60):
        u = u - (u + th(u) - s) / (1 + dth(u))
    return u


# Motorwerte aus s4 (dps 30 und 40 gleich auf allen gedruckten Stellen), Re w = Re L -+ 0.001
daten = {
    2: (mp.mpc("-2.88508868646", "0.560125095136"), mp.mpc("1.56198972697", "1.61553302785")),
    3: (mp.mpc("-2.75580372771", "0.0502174297513"), mp.mpc("1.74990531889", "1.12531692605")),
    5: (mp.mpc("-1.13105531951", "-2.70688257058"), mp.mpc("1.91953863695", "0.797061053461")),
}
for y, (sl, sr) in daten.items():
    u = u_von_s(sl)
    pred = (u + g) + th(u + g)
    groesster = max(abs(kap[m] * mp.exp(2j * mp.pi * m * u)) for m in kap if m > 0)
    print(f"y={y}: links {mp.nstr(sl,12)}  u={mp.nstr(u,10)}  Vorhersage rechts {mp.nstr(pred,12)}  Motor rechts {mp.nstr(sr,12)}  Abw {mp.nstr(abs(pred-sr),3)}  (groesstes Fourier-Glied {mp.nstr(groesster,3)})")
