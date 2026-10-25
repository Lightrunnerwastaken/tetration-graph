"""Kern-Lemma des Exponenten: H_{s}^{-1}(V) -> z_s = x_s* - 1 fuer Re V -> -oo, mit Annaeherung c * lambda_{s-2}^V
(s = 4: slog, z_s = -2, lambda_2 = b;  s = 6: ihex, lambda_4;  s = 8: ihyper(8), lambda_6).
Daraus alpha_{s+2} = ln lambda_{s-2} / ln lambda_s.  Danach (optional) Rang 9: p_8, lambda_8 -> alpha_10.
Aufruf: python p2_lemma.py DPS [rang9]"""
import sys, json, time
from pathlib import Path
import mpmath as mp
import motor
DPS = int(sys.argv[1]); R9 = len(sys.argv) > 2
B = "2"; mp.mp.dps = DPS + 10
t0 = time.time()
out = {}
def gz(v): return motor.gpzahl(v, DPS + 5)
def zeigen(name, fn, zs, lam, t, ks, cvor):
    V = [mp.mpc(-k, t) for k in ks]
    w = motor.roh(B, [f"iferr({fn}({gz(v)}),E,0)" for v in V], DPS)
    zeilen = []
    for i, k in enumerate(ks):
        if w[i] == 0: print(f"  {name} k={k}: Fehler"); continue
        d = mp.mpc(w[i]) - zs
        c = d / mp.power(lam, V[i])
        q = (mp.mpc(w[i - 1]) - zs) / d if i > 0 and w[i - 1] != 0 else None
        print(f"  {name} V=-{k}+{mp.nstr(t,4)}i: G^-1(V)-z_s = {mp.nstr(d,10)}  /lam^V = {mp.nstr(c,12)}" + (f"  Quotient = {mp.nstr(q,12)}" if q is not None else ""), flush=True)
        zeilen.append([k, str(d), str(c), str(q)])
    print(f"  Vorhersage lambda = {mp.nstr(lam,12)}, Konstante = {mp.nstr(cvor,12)}", flush=True)
    return zeilen
# s = 4: slog
r = motor.roh(B, ["pentx(0.5)", "hexx(0.5)", "lambdat", "pcent", "derivnum(t=-1,sexp(t))"], DPS)[2:]
lam4, pcent, sd = [mp.mpc(v) for v in r]
lam4 = mp.re(lam4); pcent = mp.re(pcent)
out["s4"] = zeigen("slog", "slog", mp.mpf(-2), mp.mpf(2), mp.pi / (2 * mp.log(lam4)), range(1, 13), 1 / sd)
print(f"[{time.time()-t0:.0f}s]", flush=True)
# s = 6: ihexx, z_s = x6* - 1
x6 = mp.mpc(motor.roh(B, ["ihexx(-1.7439091761327370266)"], DPS)[0])   # p_4 aus dem Bericht r18/02; unten mit Motor-Fixpunkt nachgeprueft
p4 = mp.mpf("-1.7")
for _ in range(45): p4 = mp.re(mp.mpc(motor.roh(B, [f"slog({gz(p4)})"], DPS)[0]))
x6 = mp.mpc(motor.roh(B, [f"ihexx({gz(p4)})"], DPS)[0])
h = mp.mpf(10) ** (-6)
A6 = (mp.mpc(motor.roh(B, [f"hexx({gz(x6 + h)})"], DPS)[0]) - mp.mpc(motor.roh(B, [f"hexx({gz(x6 - h)})"], DPS)[0])) / (2 * h)
print(f"p4 = {mp.nstr(p4,18)}  x6* = {mp.nstr(x6,18)}  A6 = {mp.nstr(A6,12)}", flush=True)
# hex(z) = ipent(hex(z+1)) = log(sigma(hex(z+1)-p4)) / ln lam4 - pcent  => lam4^V = sigma(..)*lam4^pcent, sigma'(0)=1
out["s6"] = zeigen("ihex", "ihexx", x6 - 1, lam4, mp.mpf("0.698"), range(1, 13), mp.power(lam4, pcent) / A6)
out["p4"], out["x6"], out["A6"] = str(p4), str(x6), str(A6)
print(f"[{time.time()-t0:.0f}s]", flush=True)
# s = 8: ihyper(8), z_s = x8* - 1
r = motor.roh(B, ["precision(hyperdigits(8),10)", "hyp[7]", "hylam[7]", "hyc[7]"], DPS)
p6, lam6, c7 = mp.re(mp.mpc(r[1])), mp.re(mp.mpc(r[2])), mp.mpc(r[3])
x8 = mp.mpc(motor.roh(B, [f"ihyper(8,{gz(p6)})"], DPS)[0])
h = mp.mpf(10) ** (-5)
A8 = (mp.mpc(motor.roh(B, [f"hyper(8,{gz(x8 + h)})"], DPS)[0]) - mp.mpc(motor.roh(B, [f"hyper(8,{gz(x8 - h)})"], DPS)[0])) / (2 * h)
print(f"Rang 8: {r[0]} Stellen  p6 = {mp.nstr(p6,18)}  lam6 = {mp.nstr(lam6,15)}  hyc7 = {mp.nstr(c7,12)}  x8* = {mp.nstr(x8,18)}  A8 = {mp.nstr(A8,12)}  [{time.time()-t0:.0f}s]", flush=True)
# H8(z) = log(sigma6(H8(z+1)-p6)/hyc7)/ln lam6  =>  z - z_s ~ hyc7 * lam6^V / A8
out["s8"] = zeigen("ihyper8", "ihyper(8,", x8 - 1, lam6, mp.mpf("0.5"), range(1, 11), c7 / A8) if False else None
V = [mp.mpc(-k, mp.mpf("0.5")) for k in range(1, 11)]
w = motor.roh(B, [f"iferr(ihyper(8,{gz(v)}),E,0)" for v in V], DPS)
z8 = []
for i, v in enumerate(V):
    if w[i] == 0: print("  ihyper8 Fehler bei", v); continue
    d = mp.mpc(w[i]) - (x8 - 1); c = d / mp.power(lam6, v)
    q = (mp.mpc(w[i - 1]) - (x8 - 1)) / d if i > 0 and w[i - 1] != 0 else None
    print(f"  ihyper8 V={mp.nstr(v,4)}: G^-1(V)-z_s = {mp.nstr(d,10)}  /lam6^V = {mp.nstr(c,12)}" + (f"  Quotient = {mp.nstr(q,12)}" if q is not None else ""), flush=True)
    z8.append([str(v), str(d), str(c), str(q)])
print(f"  Vorhersage lambda6 = {mp.nstr(lam6,12)}, Konstante hyc7/A8 = {mp.nstr(c7/A8,12)}", flush=True)
out["s8"] = z8
out.update({"p6": str(p6), "lam6": str(lam6), "hyc7": str(c7), "x8": str(x8), "A8": str(A8), "lam4": str(lam4), "pcent": str(pcent)})
P = Path(__file__).resolve().parent / "daten" / f"lemma_{DPS}.json"
json.dump(out, open(P, "w"), indent=1)
print(f"[{time.time()-t0:.0f}s] Lemma fertig", flush=True)
if R9:
    r = motor.roh(B, ["precision(hyperdigits(9),10)", "hyp[9]", "hylam[9]", "hyc[9]"], DPS)
    out.update({"r9_digits": str(r[0]), "p8": str(r[1]), "lam8": str(r[2]), "hyc9": str(r[3])})
    lam8 = mp.re(mp.mpc(r[2]))
    print(f"Rang 9: {r[0]} Stellen  p8 = {r[1]}  lam8 = {r[2]}  alpha10 = ln lam6/ln lam8 = {mp.nstr(mp.log(lam6)/mp.log(lam8),15)}  [{time.time()-t0:.0f}s]", flush=True)
    json.dump(out, open(P, "w"), indent=1)
motor.schliessen()
