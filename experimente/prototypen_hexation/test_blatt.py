"""Tests der Blattverfolgung: (1) reeller Weg = Hauptzweig, (2) geschlossener Weg kehrt zurück, (3) Weg, der den Schnitt der Kette kreuzt,
(4) Weg mit grosser Spannweite (riesige Zwischenwerte), (5) T_pfad mit Riesenwerten gegen direkte Rechnung."""
import sys, time
HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente"
sys.path.insert(0, HIER + "/prototypen_hexation")
sys.path.insert(0, HIER)
import mpmath as mp
from blatt import BlattPentation, verfolge
mp.mp.dps = 40
t0 = time.time()
BP = BlattPentation("2")
# (1) reeller Weg
ws = [mp.mpf("-1.5") + mp.mpf(k) / 10 for k in range(0, 36)]
a, s = BP.P_pfad(ws)
b = BP.Pe.P(ws)
print("(1) reell: max |P_pfad − P_haupt| =", mp.nstr(max(abs(x - y) for x, y in zip(a, b)), 3), " Sprung", mp.nstr(max(s), 3),
      f"[{time.time() - t0:.0f}s]", flush=True)
# (2) geschlossener Weg um 0.5 mit Radius 1.2 (Start und Ende bei 0.5 + 1.2 = 1.7, reell)
kreis = [mp.mpf("0.5") + mp.mpf("1.2") * mp.expjpi(2 * mp.mpf(j) / 400) for j in range(0, 401)]
v, s = BP.P_pfad(kreis)
print("(2) Kreis: Rückkehr P(Ende) − P(Start) =", mp.nstr(v[-1] - v[0], 3), " grösster Sprung:", mp.nstr(max(s), 3),
      " Hauptzweig unterwegs abweichend an", sum(1 for x, y in zip(v, BP.Pe.P(kreis)) if abs(x - y) > 1e-20), "Punkten", f"[{time.time() - t0:.0f}s]", flush=True)
# (3) senkrechter Weg 1 → 1 + 3i (Höhen) — vergleicht verfolgt gegen Hauptzweig (früher: oben 2.0745+0.99i verfolgt, 1.9727−1.3552i Haupt)
senk = [mp.mpf(1) + 1j * mp.mpf(k) / 100 for k in range(0, 301)]
senk, v3, s3 = verfolge(senk, BP.P_pfad)
s3 = [s3]
h3 = BP.Pe.P(senk)
abw = [i for i in range(len(senk)) if abs(v3[i] - h3[i]) > 1e-15]
print("(3) senkrecht 1 → 1+3i: verfolgt ≠ Hauptzweig ab Index", abw[0] if abw else None, "von", len(senk), " grösster Sprung", mp.nstr(max(s3), 3),
      " oben:", mp.nstr(v3[-1], 10), "(Haupt", mp.nstr(h3[-1], 10), ")", f"[{time.time() - t0:.0f}s]", flush=True)
# (4) Weg 1 → 2.3 reell, dann senkrecht hoch (P-Werte bis in die Tausende, Zwischenstufen riesig)
weg = [mp.mpf(1) + mp.mpf(k) / 100 for k in range(0, 131)] + [mp.mpf("2.3") + 1j * mp.mpf(k) / 100 for k in range(1, 61)]
p4, w4, s4 = verfolge(weg, BP.P_pfad)
h4 = BP.Pe.P([w for w in p4 if mp.im(w) == 0])
print("(4) 1 → 2.3 → 2.3+0.6i:", len(p4), "Punkte, grösster Sprung", mp.nstr(s4, 3), " reeller Teil gegen Hauptzweig",
      mp.nstr(max(abs(x - y) / abs(y) for x, y in zip(w4, h4)), 3), " Ende", mp.nstr(w4[-1], 10), f"[{time.time() - t0:.0f}s]", flush=True)
# (5) T_pfad direkt: Weg von −2.5+0.3i (braucht k = 1) nach 3.9+0.1i; Endwert = T(3.9+0.1i) des Hauptzweigs?
us = [mp.mpc("-2.5", "0.3") + (mp.mpc("3.9", "0.1") - mp.mpc("-2.5", "0.3")) * mp.mpf(j) / 640 for j in range(641)]
p5, w5, s5 = verfolge(us, BP.T_pfad)
print("(5) T entlang −2.5+0.3i → 3.9+0.1i:", len(p5), "Punkte, Sprung", mp.nstr(s5, 3), " Ende", mp.nstr(w5[-1], 12),
      " Hauptzweig", mp.nstr(BP.T_haupt([p5[-1]])[0], 12), f"[{time.time() - t0:.0f}s]", flush=True)
