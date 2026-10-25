"""Zertifikat (Ballarithmetik, python-flint arb/acb): genau ein Wendepunkt von Knesers sexp_b auf (-2, ∞).

Annahme H (die einzige ungeprüfte Eingabe; Prüfungen dazu in w_pruef.py), für jedes Zentrum c der Liste:
  H1  Die gespeicherten Motorwerte s_j an den M Knoten z~_j des Kreises |z - c| = r erfüllen |s_j - sexp(z~_j)| <= EPS.
  H2  |sexp(z)| <= C für |z - c| <= R  (grobe Schranke; R > r).
Alles Weitere ist streng:
  * a_k (Taylor-Koeffizienten von sexp in c): DFT der Bälle s_j in acb; Fehler EPS r^-k, Knotenverschiebung
    |z~_j - z_j| * C/(R-r-...), Aliasing sum_{m>=1} |a_{k+mM}| r^{mM} <= C R^-k q/(1-q), q = (r/R)^M.
    a_k ist reell (sexp reell auf R), also a_k im Realteil-Ball.
  * Rest der Reihe für die j-te Ableitung (j <= 3) auf |z - c| <= RHO: C * sum_{k>=K} k^3 RHO^(k-3) R^-k.
  * Ausserhalb [-1/2, 1/2] über die Turmgleichung T(x+1) = b^T(x) (Kettenregel in Bällen).
Geprüft wird auf [A, B] (Basis e: [-1.6, 0.6], Basis 2: [-1.1, 1.0]):
  sexp' > 0;  sexp'' < 0 links von J;  sexp''' > 0 auf J und Vorzeichenwechsel von sexp'' in J;  sexp'' > 0 rechts von J.
Daraus mit zwei elementaren Lemmata (siehe Bericht) genau ein Wendepunkt auf (-2, ∞); Intervall-Newton liefert ihn eng.
Aufruf: python w_reell.py e   bzw.   python w_reell.py 2      (je < 1 min); optional EPS und C: python w_reell.py e 1e-8 1000
"""
import json
import sys
import time
from fractions import Fraction

from flint import acb, arb, arb_poly, ctx

from w_daten import HIER, kreis

ctx.prec = 640

KONF = {
    "e": dict(tag="e", dps=130, M=256, zentren=["-0.4", "-0.2", "0", "0.2", "0.4"], r="0.2", R="0.3",
              C=10, EPS="1e-30", A=Fraction(-16, 10), B=Fraction(6, 10), lnb=lambda: arb(1)),
    "2": dict(tag="2", dps=80, M=256, zentren=["-0.4", "-0.2", "0", "0.2", "0.4"], r="0.2", R="0.3",
              C=10, EPS="1e-30", A=Fraction(-11, 10), B=Fraction(1, 1), lnb=lambda: arb(2).log()),
}
RHO = arb("0.11")      # Auswertung nur für |x - c| <= 0.1 (+ Ballbreite)
K = 128                # benutzte Taylor-Glieder (= M/2)


def ball(mid, rad):
    return arb(mid) + arb(0, rad)


def cball(z, rad):
    return z + acb(arb(0, rad), arb(0, rad))


def taylor(cfg, c):
    M = cfg["M"]; r = arb(cfg["r"]); R = arb(cfg["R"]); C = arb(cfg["C"]); eps = arb(cfg["EPS"])
    daten = kreis(cfg["tag"], cfg["dps"], c, cfg["r"], M)
    cc = arb(c)
    # Schranke |sexp'| <= C/(R - r - 1e-3) auf einer Umgebung der Knoten (Cauchy)
    dmax = C / (R - r - arb("1e-3"))
    s = []
    for j, (zr, zi, sr, si) in enumerate(daten):
        zt = acb(arb(zr), arb(zi))
        zj = cc + r * acb(arb(2 * j) / M).exp_pi_i()
        shift = abs(zt - zj).upper()
        if shift > arb("1e-3"):
            raise RuntimeError("Knoten zu weit weg")
        val = acb(arb(sr), arb(si))
        s.append(cball(val, (eps + shift * dmax).upper()))
    q = (r / R) ** M
    alias = C * q / (1 - q)
    a = []
    imax = arb(0)
    for k in range(K):
        w = [acb(arb(-2 * j * k) / M).exp_pi_i() for j in range(M)]
        d = sum((s[j] * w[j] for j in range(M)), acb(0)) / M / r ** k
        ak = d.real + arb(0, (alias / R ** k).upper())
        if not d.imag.contains(0):
            raise RuntimeError(f"c={c}: Im a_{k} enthält 0 nicht ({d.imag})")
        a.append(ak)
    # Rest k >= K für Ableitungen j <= 3 auf |z - c| <= RHO
    t = RHO / R
    rest = C * arb(K) ** 3 * t ** K / (1 - t * (arb(K + 1) / K) ** 3) / RHO ** 3
    if not (t * (arb(K + 1) / K) ** 3 < 1):
        raise RuntimeError("Rest divergiert")
    p0 = arb_poly(a)
    polys = [p0, p0.derivative(), p0.derivative().derivative(), p0.derivative().derivative().derivative()]
    return dict(c=cc, polys=polys, rest=rest.upper(), a=a)


class Sexp:
    def __init__(self, cfg):
        self.cfg = cfg
        self.lnb = cfg["lnb"]()
        self.T = {}
        for c in cfg["zentren"]:
            t0 = time.time()
            self.T[Fraction(c)] = taylor(cfg, c)
            tt = self.T[Fraction(c)]
            print(f"  Zentrum {c}: Taylor-Bälle in {time.time()-t0:.1f}s, Rest {tt['rest'].mid().str(3, radius=False)}, "
                  f"a_2 = {tt['a'][2].str(25)}", flush=True)

    def werte(self, lo, hi):
        """Bälle für (T, T', T'', T''') auf [lo, hi] (Fractions)."""
        if hi <= Fraction(1, 2) and lo >= Fraction(-1, 2):
            mitte = (lo + hi) / 2
            c = min(self.T, key=lambda c: abs(c - mitte))
            if max(abs(lo - c), abs(hi - c)) > Fraction(11, 100):
                raise RuntimeError("zu weit vom Zentrum")
            t = self.T[c]
            x = intervall(lo, hi)
            u = x - t["c"]
            return [p(u) + arb(0, t["rest"]) for p in t["polys"]]
        lam = self.lnb
        if hi <= Fraction(-1, 2):          # T(x) = log_b T(x+1)
            g0, g1, g2, g3 = self.werte(lo + 1, hi + 1)
            if not g0 > 0:
                raise RuntimeError("log eines nicht positiven Balls")
            v1 = g1 / g0; v2 = g2 / g0; v3 = g3 / g0
            return [g0.log() / lam, v1 / lam, (v2 - v1 ** 2) / lam, (v3 - 3 * v1 * v2 + 2 * v1 ** 3) / lam]
        if lo >= Fraction(1, 2):           # T(x) = b^T(x-1)
            g0, g1, g2, g3 = self.werte(lo - 1, hi - 1)
            E = (lam * g0).exp()
            return [E, E * lam * g1, E * (lam * g2 + lam ** 2 * g1 ** 2),
                    E * (lam * g3 + 3 * lam ** 2 * g1 * g2 + lam ** 3 * g1 ** 3)]
        raise RuntimeError(f"Intervall überquert ±1/2: {lo}, {hi}")


def arb_bruch(f):
    return arb(f.numerator) / f.denominator


def intervall(lo, hi):
    """Ball, der [lo, hi] sicher enthält."""
    rad = float((hi - lo) / 2) * (1 + 1e-9) + 1e-300
    return arb_bruch((lo + hi) / 2) + arb(0, rad)


def zertifiziere(basis):
    cfg = KONF[basis]
    print(f"Basis {basis}: EPS = {cfg['EPS']}, C = {cfg['C']}, r = {cfg['r']}, R = {cfg['R']}, Motor dps {cfg['dps']}, M = {cfg['M']}")
    S = Sexp(cfg)
    A, B = cfg["A"], cfg["B"]
    schnitte = sorted({A, B, Fraction(-3, 2), Fraction(-1, 2), Fraction(1, 2)} )
    schnitte = [s for s in schnitte if A <= s <= B]
    # Grundgitter Breite 1/500, dann adaptive Halbierung
    stuecke = []
    for lo, hi in zip(schnitte[:-1], schnitte[1:]):
        n = int((hi - lo) * 500)
        stuecke += [(lo + (hi - lo) * i / n, lo + (hi - lo) * (i + 1) / n) for i in range(n)]
    erg = []   # (lo, hi, vorzeichen T'', vorzeichen T''')
    zaehler = dict(auswertungen=0)

    def vz(b):
        return 1 if b > 0 else (-1 if b < 0 else 0)

    def pruefe(lo, hi, tiefe=0):
        zaehler["auswertungen"] += 1
        w = S.werte(lo, hi)
        if not w[1] > 0:
            if tiefe < 12:
                m = (lo + hi) / 2
                return pruefe(lo, m, tiefe + 1) + pruefe(m, hi, tiefe + 1)
            raise RuntimeError(f"T' > 0 nicht gezeigt auf [{float(lo)}, {float(hi)}]: {w[1]}")
        s2, s3 = vz(w[2]), vz(w[3])
        if s2 == 0 and s3 == 0 and tiefe < 12:
            m = (lo + hi) / 2
            return pruefe(lo, m, tiefe + 1) + pruefe(m, hi, tiefe + 1)
        return [(lo, hi, s2, s3)]

    t0 = time.time()
    for lo, hi in stuecke:
        erg += pruefe(lo, hi)
    print(f"  {len(erg)} Teilintervalle, {zaehler['auswertungen']} Auswertungen, {time.time()-t0:.1f}s")
    # Auswertung der Vorzeichenfolge
    unbest = [i for i, e in enumerate(erg) if e[2] == 0]
    if not unbest:
        raise RuntimeError("kein Kandidat für den Wendepunkt")
    i0, i1 = unbest[0], unbest[-1]
    # J = zusammenhängender Block um die unbestimmten Stücke, auf dem T''' > 0
    if not all(erg[i][3] == 1 for i in range(i0, i1 + 1)):
        raise RuntimeError("T''' > 0 auf J nicht gezeigt")
    links = erg[:i0]; rechts = erg[i1 + 1:]
    assert all(e[2] == -1 for e in links), "T'' < 0 links nicht gezeigt"
    assert all(e[2] == 1 for e in rechts), "T'' > 0 rechts nicht gezeigt"
    Jlo, Jhi = erg[i0][0], erg[i1][1]
    # Vorzeichen an den Rändern von J folgt aus den Nachbarstücken (abgeschlossene Intervalle, gemeinsame Endpunkte)
    print(f"  T' > 0 auf [{float(A)}, {float(B)}]")
    print(f"  T'' < 0 auf [{float(A)}, {float(Jlo)}],  T''' > 0 auf J = [{float(Jlo)}, {float(Jhi)}],  T'' > 0 auf [{float(Jhi)}, {float(B)}]")
    # Intervall-Newton für die Nullstelle von T'' in J
    lo, hi = Jlo, Jhi
    for it in range(60):
        m = (lo + hi) / 2
        fm = S.werte(m, m)[2]
        d = S.werte(lo, hi)[3]
        N = arb_bruch(m) - fm / d
        # N als Bruch-Schranken (nach aussen gerundet)
        lo_n = _frac_unten(N); hi_n = _frac_oben(N)
        neu_lo, neu_hi = max(lo, lo_n), min(hi, hi_n)
        if neu_lo > neu_hi:
            raise RuntimeError("Newton: leerer Schnitt (keine Nullstelle)")
        if (neu_hi - neu_lo) >= (hi - lo) * Fraction(99, 100):
            lo, hi = neu_lo, neu_hi
            break
        lo, hi = neu_lo, neu_hi
    breite = float(hi - lo)
    mitte = intervall(lo, hi)
    print(f"  Wendepunkt x_w in {mitte.str(40)}   (Breite {breite:.2e})")
    w = S.werte(lo, hi)
    print(f"  dort sexp = {w[0].str(30)}, sexp' = {w[1].str(30)}, sexp''' = {w[3].str(20)}")
    out = dict(basis=basis, EPS=cfg["EPS"], C=cfg["C"], r=cfg["r"], R=cfg["R"], A=str(cfg["A"]), B=str(cfg["B"]),
               J=[str(Jlo), str(Jhi)], J_float=[float(Jlo), float(Jhi)], x_w=mitte.str(40), breite=breite,
               teilintervalle=len(erg), sexp_w=w[0].str(30), sexp1_w=w[1].str(30), sexp3_w=w[3].str(20))
    (HIER / f"zertifikat_{basis}_eps{cfg['EPS']}_C{cfg['C']}.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
    return out


def _frac_unten(x):
    s = x.lower().str(80, radius=False)
    f = Fraction(s)
    return f - Fraction(1, 10 ** 70)


def _frac_oben(x):
    s = x.upper().str(80, radius=False)
    f = Fraction(s)
    return f + Fraction(1, 10 ** 70)


if __name__ == "__main__":
    b = sys.argv[1] if len(sys.argv) > 1 else "e"
    if len(sys.argv) > 2:                 # Empfindlichkeit: python w_reell.py e 1e-10 1000
        KONF[b]["EPS"] = sys.argv[2]
        KONF[b]["C"] = int(sys.argv[3]) if len(sys.argv) > 3 else KONF[b]["C"]
    zertifiziere(b)
