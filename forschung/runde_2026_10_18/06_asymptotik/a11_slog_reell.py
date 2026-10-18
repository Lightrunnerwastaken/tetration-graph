"""A11: slog auf der reellen Achse gegen den linear interpolierten log*.
slog_lin(x) = n - 1 + log^[n] x  mit log^[n] x in (0,1].  Exakt: slog(x) - slog_lin(x) = P(slog x - n + 1),  P(u) = u + 1 - sexp(u), u in [-1,0]
(1-periodische Fortsetzung). Kein Grenzwert; Mittel = 1/2 - int_{-1}^0 sexp; Schwankungsbreite = max P - min P.
"""
import json, mpmath as mp
mp.mp.dps = 40
for f in ["daten/taylor_0_1.5_400_50.json", "daten/taylor_0_1_512_60.json"]:
    C = [mp.mpf(s) for s in json.load(open(f))][:180]
    S = lambda u: mp.polyval(C[::-1], u)
    dS = [k * C[k] for k in range(1, len(C))]
    S1 = lambda u: mp.polyval(dS[::-1], u)
    P = lambda u: u + 1 - S(u)
    # Extrema: sexp'(u) = 1
    us = [mp.findroot(lambda u: S1(u) - 1, u0) for u0 in (-0.8, -0.2)]
    I = mp.quad(S, [-1, 0])
    print(f, "\n  Extremstellen u:", [mp.nstr(u, 20) for u in us], "\n  P dort:", [mp.nstr(P(u), 20) for u in us])
    print("  int_{-1}^0 sexp =", mp.nstr(I, 25), "  Mittel von P = 1/2 - I =", mp.nstr(mp.mpf(1)/2 - I, 20))
    print("  Nullstelle von P in (-1,0):", mp.nstr(mp.findroot(P, -0.5), 20))
    print("  int P^2 =", mp.nstr(mp.quad(lambda u: P(u)**2, [-1, 0]), 15))
