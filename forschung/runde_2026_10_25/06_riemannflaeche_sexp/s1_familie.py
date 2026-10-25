"""S1: Gitterfamilien der reinen Log-Verzweigungspunkte.
Loesungen von sexp(w) = c (c = 2 pi i k) laengs w0 + j*gamma, gamma = 2 pi i / L; Newton mit Motor (Hauptblatt).
Aufruf: python s1_familie.py DPS     (zwei Einstellungen: 40 und 30)
"""
import sys, json
import mpmath as mp
from motor import werte, gpzahl, schliessen

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 40
mp.mp.dps = DPS + 10

def sexp(zs):
    return werte("exp(1)", [f"sexp({gpzahl(z, DPS + 10)})" for z in zs], DPS)

L = mp.mpc("0.318131505204764135312654251587", "1.33723570143068940890116214319")
gam = 2j * mp.pi / L

def newton(w, c, it=25):
    h = mp.mpf(10) ** (-(DPS // 3))
    for _ in range(it):
        a, b, f0 = sexp([w + h, w - h, w])
        d = (a - b) / (2 * h)
        st = (f0 - c) / d
        if abs(st) > 0.5:
            st = st / abs(st) * 0.5
        w = w - st
        if abs(st) < mp.mpf(10) ** (-DPS + 6):
            return w, d
    return w, d

aus = {}
for k, w0 in [(1, mp.mpc("2.0043124825946880039", "0.70941354086940362668")),
              (2, mp.mpc("2.15587108746", "0.49064172727")),
              (-1, mp.mpc("2.0043124825946880039", "-0.70941354086940362668"))]:
    c = 2j * mp.pi * k
    ws = []
    for j in range(0, 3):
        g = w0 + j * (gam if k > 0 else mp.conj(gam)) if not ws else ws[-1][0] + (gam if k > 0 else mp.conj(gam))
        w, d = newton(g, c)
        ws.append((w, d))
        print(f"k={k} j={j}: w={mp.nstr(w, 22)}  sexp'(w)={mp.nstr(d, 8)}", flush=True)
    for j in range(len(ws) - 1):
        dd = ws[j + 1][0] - ws[j][0] - (gam if k > 0 else mp.conj(gam))
        print(f"   w_{j+1}-w_{j}-gamma = {mp.nstr(dd, 6)}  |.|={mp.nstr(abs(dd), 4)}", flush=True)
    aus[str(k)] = [[mp.nstr(mp.re(w), DPS - 5), mp.nstr(mp.im(w), DPS - 5)] for w, d in ws]
print("gamma =", mp.nstr(gam, 20), " exp(-2 pi Im gamma) =", mp.nstr(mp.exp(-2 * mp.pi * mp.im(gam)), 6))
json.dump(aus, open(f"daten/familie_{DPS}.json", "w"), indent=1)
schliessen()
