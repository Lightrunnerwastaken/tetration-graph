"""F4: Konkavitaet von slog und die Ableitungs-Welle Q(u) = 1/sexp'(u) = slog'(x)/slog_lin'(x), Basis e, Kneser.
(a) Schranken auf [-1,0], die den Ein-Wendepunkt-Satz fuer sexp_e auf (-2,oo) auf eine Periode zurueckfuehren:
    sexp'' hat genau eine Nullstelle in [-1,0]; sexp''' > 0 auf [-1,0]; min sexp' = m; max|sexp''| = |sexp''(-1)|.
(b) Q: Mittel q0 = int_0^1 slog'(x)^2 dx >= 1; Sprung von Q' an den ganzen Zahlen (= 1 exakt); Bernoulli-Naeherung q0 - B_2(u+1)/2 + ...
Zwei Taylor-Saetze als Gegenprobe (Argument: Dateiname).
"""
import json, sys
from pathlib import Path
import mpmath as mp

HIER = Path(__file__).resolve().parent
mp.mp.dps = 32
satz = sys.argv[1] if len(sys.argv) > 1 else "taylor_0_1_512_60.json"
C0 = [mp.mpf(s) for s in json.load(open(HIER / "daten" / satz))][:180]
def abl(j):
    co = [mp.ff(k, j) * C0[k] for k in range(j, len(C0))]
    return lambda u: mp.polyval(co[::-1], u)
S, S1, S2, S3, S4 = (abl(j) for j in range(5))
G = mp.linspace(-1, 0, 401)
print("Satz:", satz)
print("min sexp' auf [-1,0] (Gitter 401):", mp.nstr(min(S1(u) for u in G), 20))
print("min sexp''' auf [-1,0]:", mp.nstr(min(S3(u) for u in G), 12), "  bei u =", mp.nstr(min(G, key=S3), 6))
print("sexp''(-1) =", mp.nstr(S2(-1), 20), " sexp''(0) =", mp.nstr(S2(0), 20), " max|sexp''| =", mp.nstr(max(abs(S2(u)) for u in G), 20))
print("Vorzeichenwechsel sexp'' auf Gitter:", sum(1 for a, b in zip(G[:-1], G[1:]) if S2(a) * S2(b) < 0))
# Konvexitaet auf [0,1]: sexp''(u+1) = sexp(u+1) (sexp'(u)^2 + sexp''(u)) > 0 ?
print("min (sexp'^2 + sexp'') auf [-1,0]:", mp.nstr(min(S1(u) ** 2 + S2(u) for u in G), 12))
# Konkavitaet auf (-2,-1]: Zaehler sexp''(v) sexp(v) - sexp'(v)^2 < 0 fuer v in (-1,0]
print("max (sexp'' sexp - sexp'^2) auf [-1,0]:", mp.nstr(max(S2(v) * S(v) - S1(v) ** 2 for v in G), 12))
xw = mp.findroot(S2, -0.52)
m = S1(xw)
print("x_w =", mp.nstr(xw, 28), " X_w = sexp(x_w) =", mp.nstr(S(xw), 28))
print("max slog' auf R = 1/m =", mp.nstr(1 / m, 28), "  slog''(X_w)=0;  23/22 =", mp.nstr(mp.mpf(23) / 22, 12))
print("slog'(0) = 1/c =", mp.nstr(1 / C0[1], 25))
# (b) Q
Q = lambda u: 1 / S1(u)
q0 = mp.quad(Q, [-1, xw, 0])
print("q0 = int_{-1}^0 du/sexp'(u) = int_0^1 slog'^2 dx =", mp.nstr(q0, 25))
Q1 = lambda u: -S2(u) / S1(u) ** 2
Q2 = lambda u: (2 * S2(u) ** 2 - S3(u) * S1(u)) / S1(u) ** 3
print("Sprung Q'(0+)-Q'(0-) = Q'(-1)-Q'(0) =", mp.nstr(Q1(-1) - Q1(0), 25))
JQ2 = Q2(-1) - Q2(0)
print("Sprung Q''(0+)-Q''(0-) =", mp.nstr(JQ2, 20))
for name, F in [("q0 - B2/2", lambda u: q0 - mp.bernpoly(2, u + 1) / 2),
                ("q0 - B2/2 - JQ2*B3/6", lambda u: q0 - mp.bernpoly(2, u + 1) / 2 - JQ2 * mp.bernpoly(3, u + 1) / 6)]:
    print(f"  max|Q - ({name})| =", mp.nstr(max(abs(Q(u) - F(u)) for u in G), 6))
print("Spannweite Q:", mp.nstr(1 / m - 1 / C0[1], 12), "  Bernoulli-Spannweite 1/8 =", 0.125)
