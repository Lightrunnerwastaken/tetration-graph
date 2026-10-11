"""Modell 8: lange Zweig-Läufe zur Schärfe der Schwelle und zur Zweigwahl.
 (a) elliptische Drehung c = 1.05 (knapp ÜBER π/3 = 1.0472, vorhergesagte Rate ρ = sin(c/2)/cos c = 1.0073): Fehler fällt erst,
     muss dann wieder wachsen (Umkehr bei grossem N).
 (b) hyperbolisch c = 2: bleibt der Zweig auf der falschen Familie, oder findet er f_M doch noch?"""
import sys, time
import mpmath as mp
from kern import zweig, g_moebius_ell, g_moebius_hyp

mp.mp.dps = 60
fall = sys.argv[1]
NMAX = int(sys.argv[2])
t0 = time.time()
if fall == "ell":
    c = mp.mpf(sys.argv[3]) if len(sys.argv) > 3 else mp.mpf("1.05")
    ref = mp.tan(c / 2)
    gk = g_moebius_ell(c)
    print(f"elliptisch c = {c}: ρ = {mp.nstr(mp.sin(c / 2) / mp.cos(c), 6)}", flush=True)
else:
    c = mp.mpf(sys.argv[3]) if len(sys.argv) > 3 else mp.mpf(2)
    ref = mp.tanh(c / 2)
    gk = g_moebius_hyp(c)
    print(f"hyperbolisch c = {c}: T = {mp.nstr(ref, 6)}", flush=True)


def aus(N, a):
    if N % 5 == 0:
        print(f"N={N:4d} a0−f(0) = {mp.nstr(a[0] - ref, 5):>12}  ({time.time() - t0:.0f}s)", flush=True)


zweig(gk, NMAX, bits=int(sys.argv[4]) if len(sys.argv) > 4 else 1200, ausgabe=aus, still=False)
