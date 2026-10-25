"""s4: (a) Schnitt des Motor-Blatts auf der Senkrechten Re w = Re L: Werte links/rechts.
(b) s* = slog_K(L + 2 pi i k) auf dem Kneser-Blatt (rechts), Pruefung sexp(s*) = L + 2 pi i k und sexp(s*+1) = L.
Aufruf: python s4_schnitt_und_sstern.py <dps>"""
import sys
import mpmath as mp
import motor

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
mp.mp.dps = dps + 10
L = mp.conj(-mp.lambertw(-1))  # oberer Fixpunkt von exp, L = e^L, Im L > 0 (-W_0(-1) ist der untere!)
g = 2j * mp.pi / L
d = mp.mpf("0.001")
ys = [mp.mpf(2), mp.mpf(3), mp.mpf(5), mp.mpf(7), mp.mpf(9), mp.mpf(12)]
aus = []
for y in ys:
    aus += [f"slog({motor.gpzahl(mp.mpc(L.real - d, y), dps)})", f"slog({motor.gpzahl(mp.mpc(L.real + d, y), dps)})"]
for k in [1, 2, 3, -1]:
    aus.append(f"slog({motor.gpzahl(L + 2j * mp.pi * k + mp.mpf('1e-30'), dps)})")
res = motor.werte("exp(1)", aus, dps)
print("(a) Senkrechte Re w = Re L +- 0.001")
for i, y in enumerate(ys):
    a, b = res[2 * i], res[2 * i + 1]
    print(f"  y={mp.nstr(y,3):>4}: links {mp.nstr(a,12):>34}  rechts {mp.nstr(b,12):>34}  rechts-links {mp.nstr(b-a,10)}  (gamma={mp.nstr(g,8)})")
ss = res[2 * len(ys):]
aus2 = []
for s in ss:
    aus2 += [f"sexp({motor.gpzahl(s, dps)})", f"sexp({motor.gpzahl(s + 1, dps)})"]
r2 = motor.werte("exp(1)", aus2, dps)
print("(b) s*_k = slog(L + 2 pi i k) rechts der Senkrechten")
for j, k in enumerate([1, 2, 3, -1]):
    s = ss[j]
    print(f"  k={k:>2}: s* = {mp.nstr(s, dps-2)}")
    print(f"        |sexp(s*) - (L+2pi i k)| = {mp.nstr(abs(r2[2*j] - L - 2j*mp.pi*k), 3)}   |sexp(s*+1) - L| = {mp.nstr(abs(r2[2*j+1] - L), 3)}")
motor.schliessen()
