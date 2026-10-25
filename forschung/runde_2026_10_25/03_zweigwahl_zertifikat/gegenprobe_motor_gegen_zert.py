"""Gegenprobe beim Eintragen (25.10.2026) zu Annahme (E): Motor-θ (Etappe 2, a_0 … a_4096) gegen das zertifizierte θ̃ aus
kn.satz.kneser_e_zertifiziert (forschung/polynom_halbexp/p_pruefung/p26_iter6.json; dort |θ_K − θ̃| ≤ 4.4e-167 für Im z ≥ 0.36, p32a).
Nur Lesen, Gleitkomma (mpmath, dps 80) — ein Beleg für (E), kein Beweis (bei z₀ = 0.05i braucht man θ_K nahe der Achse).

Aufruf (aus diesem Ordner):
  C:/Users/janis/miniconda3/python.exe -I -X utf8 gegenprobe_motor_gegen_zert.py \
      ../../identitaeten_2026_10_14/etappe2/daten/a_e_80_0.003_16384.json ../../polynom_halbexp/p_pruefung/p26_iter6.json
Ergebnis 25.10.2026: θ_motor − θ̃ = 3·(2πi/L) bis auf ≤ 3.2e-79 an 0.36i, 0.4i, 0.3+0.5i, i; max_{1≤n≤60} |a_n − c̃_n| = 5.3e-82.
"""
import json
import sys

import mpmath as mp

mp.mp.dps = 80
motor = json.load(open(sys.argv[1]))
zert = json.load(open(sys.argv[2]))
a = [mp.mpc(mp.mpf(r), mp.mpf(i)) for r, i in motor['a']]


def als_c(x):
    if isinstance(x, (list, tuple)):
        return mp.mpc(mp.mpf(str(x[0])), mp.mpf(str(x[1])))
    return mp.mpc(mp.mpf(str(x)))


c = [als_c(x) for x in zert['c']]
print('Moden: Motor', len(a), ' zertifiziert', len(c))
L = -mp.lambertw(-1)
L = L if mp.im(L) > 0 else mp.conj(L)        # Fixpunkt e^L = L mit Im L > 0, auf 80 Stellen


def th(coef, z):
    q = mp.exp(2j * mp.pi * z)
    return sum(cn * q ** n for n, cn in enumerate(coef))


for z in [mp.mpc(0, '0.36'), mp.mpc(0, '0.4'), mp.mpc('0.3', '0.5'), mp.mpc(0, '1')]:
    d = th(a, z) - th(c, z)
    k = d * L / (2j * mp.pi)    # Verschiebung um 2πi k/L ist erlaubt (mod 2πi/ℓ)
    kr = mp.nint(mp.re(k))
    print('z =', mp.nstr(z, 4), ' θ_motor − θ̃ =', mp.nstr(d, 5), '  in Einheiten 2πi/L:', mp.nstr(k, 5),
          f'  |θ_motor − θ̃ − {int(kr)}·2πi/L| =', mp.nstr(abs(d - kr * 2j * mp.pi / L), 3))
m = max(abs(a[n] - c[n]) for n in range(1, 61))
print('max_{1≤n≤60} |a_n − c̃_n| =', mp.nstr(m, 4), ' (zertifiziert: |κ_n − c̃_n| ≤ 4.4e-167·e^{2πn·0.36}, bei n = 60:',
      mp.nstr(mp.mpf('4.4e-167') * mp.exp(2 * mp.pi * 60 * mp.mpf('0.36')), 3), ')')
print('a_0 Motor =', mp.nstr(a[0], 15), '  c̃_0 =', mp.nstr(c[0], 15))
