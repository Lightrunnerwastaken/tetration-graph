"""Asymptotik-Analyse der Folgen eta_r, e_r, beta_r, iota_r (r=3..8).
Liest eta_values.csv (Kopie aus forschung/eta_hyperoperationen) und cache/jet_*.json.
Aufruf: python -B analyse.py   (Laufzeit < 2 s)
"""
import csv, json, glob
from pathlib import Path
import mpmath as mp
mp.mp.dps = 40
HERE = Path(__file__).absolute().parent
rows = list(csv.DictReader((HERE / "eta_values.csv").open(encoding="utf-8-sig")))
R = [int(z["map_rank"]) for z in rows]
eta = {int(z["map_rank"]): mp.mpf(z["eta"]) for z in rows}
E = {int(z["map_rank"]): mp.mpf(z["e"]) for z in rows}
d = {r: 2 - eta[r] for r in R}; eps = {r: 4 - E[r] for r in R}

def p_loc(s, r, shift=0):
    return mp.log(s[r-1] / s[r]) / mp.log((r + shift) / (r - 1 + shift))

def fit3(s, r1, r2, r3):
    """A*(r+c)^(-p) exakt durch drei Punkte; gibt (A,c,p)."""
    def g(c):
        return p_loc2(s, r1, r2, c) - p_loc2(s, r2, r3, c)
    def p_loc2(s, ra, rb, c):
        return mp.log(s[ra] / s[rb]) / mp.log((rb + c) / (ra + c))
    c = mp.findroot(g, 0.0)
    p = p_loc2(s, r1, r2, c)
    A = s[r3] * (r3 + c) ** p
    return A, c, p

def fit_lim(s, r1, r2, r3):
    """s_r = L - A r^(-q) exakt durch drei Punkte (Aitken-artig ueber q)."""
    def g(q):
        return (s[r2] - s[r1]) / (r1**-q - r2**-q) - (s[r3] - s[r2]) / (r2**-q - r3**-q)
    q = mp.findroot(g, 1.0)
    A = (s[r2] - s[r1]) / (r1**-q - r2**-q)
    L = s[r1] + A * r1**-q
    return L, A, q

print("r  delta  eps  eps/delta  p_loc(delta,r) p_loc(delta,m=r-2) p_loc(eps,r)")
for r in R:
    line = f"{r} {mp.nstr(d[r],12)} {mp.nstr(eps[r],12)} {mp.nstr(eps[r]/d[r],8)}"
    if r > 3:
        line += f" {mp.nstr(p_loc(d,r),8)} {mp.nstr(p_loc(d,r,-2),8) if r>3 else ''} {mp.nstr(p_loc(eps,r),8)}"
    print(line)

print("\n3-Parameter-Modell delta = A (r+c)^-p")
for tri in [(4,5,6),(5,6,7),(6,7,8)]:
    A,c,p = fit3(d,*tri)
    pred = {k: 2 - A*(k+c)**-p for k in (8,9,10)}
    print(tri, "A",mp.nstr(A,8),"c",mp.nstr(c,8),"p",mp.nstr(p,8),
          "eta8_pred",mp.nstr(pred[8],10),"err",mp.nstr(pred[8]-eta[8],3),"eta9_pred",mp.nstr(pred[9],10),"eta10",mp.nstr(pred[10],10))
print("\n3-Parameter-Modell eps = A (r+c)^-p")
for tri in [(4,5,6),(5,6,7),(6,7,8)]:
    A,c,p = fit3(eps,*tri)
    pred = {k: 4 - A*(k+c)**-p for k in (8,9,10)}
    print(tri, "A",mp.nstr(A,8),"c",mp.nstr(c,8),"p",mp.nstr(p,8),
          "e8_pred",mp.nstr(pred[8],10),"err",mp.nstr(pred[8]-E[8],3),"e9_pred",mp.nstr(pred[9],10))

# Grenzwert-Tests ohne Annahme L=2 bzw. 4
print("\nFreier Grenzwert L - A r^-q durch drei Punkte")
for name, s in [("eta", eta), ("e", E)]:
    for tri in [(4,5,6),(5,6,7),(6,7,8)]:
        try:
            L,A,q = fit_lim(s,*tri); print(name, tri, "L",mp.nstr(L,8),"q",mp.nstr(q,6))
        except Exception as ex:
            print(name, tri, "kein Fit", ex)

# Jets
jets = {}
for f in sorted(glob.glob(str(HERE/"cache"/"jet_r*.json"))):
    j = json.load(open(f))
    if j.get("status") == "completed":
        jets.setdefault(j["rank"], []).append(j)
beta = {3: 1/(2*mp.e)}; iota = {3: mp.mpf(1)/3}
print("\nJets (beta = F''/2, iota = 1 - a3/a2^2):")
for r in sorted(jets):
    for j in jets[r]:
        print(r, "M",j["M"],"rho",j["rho"],"dps",j["dps"],"beta",j["beta"][:18],"iota",j["iota"][:18])
    beta[r] = mp.mpf(jets[r][0]["beta"]); iota[r] = mp.mpf(jets[r][0]["iota"])
if 8 not in beta:
    beta[8] = mp.mpf("0.5368242374")  # aus BERICHT der Runde 12.10. (Abschnitt 8, Punkt 4)
print("\nbeta_r:", {r: mp.nstr(beta[r],10) for r in sorted(beta)})
print("iota_r:", {r: mp.nstr(iota[r],10) for r in sorted(iota)})
print("Inkremente beta:", [mp.nstr(beta[r]-beta[r-1],6) for r in sorted(beta) if r-1 in beta])
print("Inkremente iota:", [mp.nstr(iota[r]-iota[r-1],6) for r in sorted(iota) if r-1 in iota])
for tri in [(5,6,7),(6,7,8)]:
    if all(t in beta for t in tri):
        try:
            L,A,q = fit_lim(beta,*tri); print("beta Grenzfit", tri, "L",mp.nstr(L,8),"q",mp.nstr(q,6))
        except Exception as ex: print("beta fit", ex)
# Skaleninvariante Kombinationen
print("\nbeta*eps, beta*eps^2, beta*(e-2), beta*delta^(1/2), beta*(e-2)^2:")
for r in sorted(beta):
    print(r, mp.nstr(beta[r]*eps[r],8), mp.nstr(beta[r]*eps[r]**2,8), mp.nstr(beta[r]*(E[r]-2),8),
          mp.nstr(beta[r]*mp.sqrt(d[r]),8), mp.nstr(beta[r]*(E[r]-2)**2,8))

# Normierungen: v_r = Abstand * N(r); Aitken-Extrapolation der Folge v_r
print("\nNormierungen (Aitken aus je drei aufeinanderfolgenden v_r):")
for name, s in [("delta", d), ("eps", eps)]:
    for lab, f in [("r^1.5", lambda r: r**1.5), ("r ln r", lambda r: r*mp.log(r)),
                   ("(r-1) ln(r-1)", lambda r: (r-1)*mp.log(r-1)), ("r ln(r)^2", lambda r: r*mp.log(r)**2)]:
        v = [s[r]*f(r) for r in range(3, 9)]
        ai = [(v[i+2]*v[i]-v[i+1]**2)/(v[i+2]-2*v[i+1]+v[i]) for i in range(len(v)-2)]
        print(name, lab, [mp.nstr(x, 6) for x in v], "Aitken", [mp.nstr(x, 6) for x in ai])
v = [d[r]*r*mp.log(r) for r in range(3, 9)]
q = (v[5]-v[4])/(v[4]-v[3])
v9 = v[5] + (v[5]-v[4])*q
print("delta*r*ln r: Quotient der letzten Zuwaechse", mp.nstr(q, 5), "-> v9", mp.nstr(v9, 6),
      "eta9(r ln r)", mp.nstr(2 - v9/(9*mp.log(9)), 8),
      "eta9(Potenzfit r=3..6)", mp.nstr(2 - mp.mpf("2.7849974")*mp.mpf(9)**mp.mpf("-1.4671731"), 8))
