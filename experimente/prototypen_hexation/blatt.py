"""Blattverfolgung: analytische Fortsetzung der Kneser-Tetration T_b und der Pentation P_b entlang von Wegen.

 T ist im Hauptzweig auf ℂ ∖ (−∞, −2] analytisch. Fortsetzung über den Schnitt: T(u) = log_b T(u+1) (mit Zweigwahl), also
   T(u) = L_{n_1}( L_{n_2}( … L_{n_k}( T(u + k) ) … ) ),   L_n(v) = (Log v + 2πi n) / ln b,
 mit T(u + k) im Hauptzweig, solange Re(u + k) > −2 entlang des ganzen Weges. Die n_ℓ werden am Wegstart aus dem Hauptzweig bestimmt
 (dort korrekt) und danach nach Stetigkeit gewählt. Jeder Schritt meldet die Winkeländerung |Δ Log v| zwischen Nachbarpunkten;
 ist sie zu gross (> ~1), muss der Weg verfeinert werden (verfolge()).
"""
import sys
HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente"
if HIER not in sys.path:
    sys.path.insert(0, HIER)
import mpmath as mp

from gemeinsam import kneser
import cache_drossel  # noqa: F401  (Zwischenspeicher seltener schreiben)
from stufen import basis_text, psi_wert
from hexation import Pentation

ZWEI_PI_I = 2j * mp.pi


class BlattPentation:
    def __init__(self, b, M=8):
        self.b = b
        self.Pe = Pentation(b, M=M)
        self.M = M
        self.lnb = mp.mpf(1) if b == "exp(1)" else mp.log(mp.mpf(b))

    def T_haupt(self, us):
        return kneser("sexp", basis_text(self.b), us)

    GRENZE = mp.mpf("4.5")       # bis hier rechnet der Motor; darüber exp in mpmath
    XA = mp.mpf("1.5")           # ab hier analytischer Logarithmus log T(x) = ln b · T(x − 1)
    k_extra = 0                  # zusätzliche Verschiebung (nur zur Kontrolle: Ergebnis darf sich nicht ändern)

    def referenzen(self, xs):
        """Hauptzweig T(x) für viele x; None, wenn astronomisch gross oder das Argument nicht mehr genau ist."""
        out = [None] * len(xs)
        klein = [i for i, x in enumerate(xs) if mp.re(x) <= self.GRENZE]
        if klein:
            for i, v in zip(klein, self.T_haupt([xs[i] for i in klein])):
                out[i] = v
        gross = [i for i, x in enumerate(xs) if mp.re(x) > self.GRENZE]
        if gross:
            js = {i: int(mp.ceil(mp.re(xs[i]) - self.GRENZE)) for i in gross}
            schranke = mp.mpf(10) ** (mp.mp.dps - 10)
            for i, v in zip(gross, self.T_haupt([xs[i] - js[i] for i in gross])):
                for _ in range(js[i]):
                    a = self.lnb * v
                    if abs(mp.re(a)) > schranke or abs(mp.im(a)) > schranke:
                        v = None
                        break
                    v = mp.exp(a)
                out[i] = v
        return out

    def T_pfad(self, us):
        """T entlang des Weges us (geordnet); Start us[0] im Hauptzweig korrekt. Rückgabe: Werte, Sprung je Punkt.

        Stufe l trägt den Wert V an u + l, entweder als Zahl ('z') oder als Abstand D zum Hauptzweig, V = T(u + l) + D ('d').
        Für 'd' ist log V = ln b · T(u + l − 1) + log1p(D / T(u + l)) + 2πi m exakt; m ändert sich zwischen zwei 'd'-Nachbarn nicht.
        So werden riesige Werte nie in Winkeln verfolgt, nur die mässigen ('z') nach Stetigkeit."""
        N, lnb = len(us), self.lnb
        sprung = [mp.mpf(0)] * N
        k = max(0, max(int(mp.ceil(mp.mpf("-1.7") - mp.re(u))) for u in us)) + self.k_extra
        if k == 0:
            v = self.referenzen(us)
            if any(w is None for w in v):
                raise ValueError("Wert zu gross")
            return v, sprung
        bedarf = sorted({(i, l) for i in range(N) for l in range(k + 1)
                         if l == k or i == 0 or mp.re(us[i] + l) >= self.XA - 1})
        R = dict(zip(bedarf, self.referenzen([us[i] + l for i, l in bedarf])))
        zust = [("d", mp.mpc(0)) if mp.re(us[i] + k) >= self.XA else ("z", R[(i, k)]) for i in range(N)]
        if any(a == "z" and w is None for a, w in zust):
            raise ValueError("Wert zu gross")
        for l in range(k, 0, -1):
            neu, vorher = [], None                      # vorher = (Art, stetiger Log, m, Korrektur, u)
            for i in range(N):
                art, w = zust[i]
                if art == "d":
                    Rl, Rl1 = R[(i, l)], R[(i, l - 1)]
                    korr = mp.log1p(w / Rl) if Rl is not None else mp.mpc(0)
                    basis = lnb * Rl1 if Rl1 is not None else None
                    if i == 0:
                        m = 0
                    elif vorher[0] == "d":
                        m = vorher[2]
                        sprung[i] = max(sprung[i], abs(korr - vorher[3]) + abs(us[i] - vorher[4]) / 2)
                    elif basis is None or vorher[1] is None:
                        m, sprung[i] = 0, mp.inf
                    else:
                        m = int(mp.nint(mp.im(vorher[1] - basis - korr) / (2 * mp.pi)))
                        sprung[i] = max(sprung[i], abs(basis + korr + ZWEI_PI_I * m - vorher[1]))
                    lw = basis + korr + ZWEI_PI_I * m if basis is not None else None
                    D = (korr + ZWEI_PI_I * m) / lnb
                    vorher = ("d", lw, m, korr, us[i])
                    V = None
                else:
                    lw0 = mp.log(w)
                    if i == 0:
                        n = mp.nint(mp.im(lnb * R[(0, l - 1)] - lw0) / (2 * mp.pi))
                    elif vorher[1] is None:
                        n, sprung[i] = 0, mp.inf
                    else:
                        n = mp.nint(mp.im(vorher[1] - lw0) / (2 * mp.pi))
                    lw = lw0 + ZWEI_PI_I * n
                    if i > 0 and vorher[1] is not None:
                        sprung[i] = max(sprung[i], abs(lw - vorher[1]))
                    vorher = ("z", lw, None, None, us[i])
                    V = lw / lnb
                R1 = R.get((i, l - 1))
                gross = mp.re(us[i] + l - 1) >= self.XA
                if art == "d":
                    if R1 is None or (gross and abs(D) < abs(R1) / 2):
                        neu.append(("d", D))
                    else:
                        neu.append(("z", R1 + D))
                else:
                    if gross and R1 is not None and abs(V - R1) < abs(R1) / 2:
                        neu.append(("d", V - R1))
                    else:
                        neu.append(("z", V))
            zust = neu
        werte = []
        for i, (art, w) in enumerate(zust):
            if art == "z":
                werte.append(w)
            elif R[(i, 0)] is None:
                raise ValueError("Wert zu gross")
            else:
                werte.append(R[(i, 0)] + w)
        return werte, sprung

    def P_pfad(self, ws):
        """Pentation entlang des Weges ws (Höhen, geordnet); Start ws[0] mit korrektem Hauptzweig. Rückgabe: Werte, Sprung je Punkt."""
        Pe = self.Pe
        ys = [psi_wert(Pe.L, Pe.psiT, Pe.c * Pe.lam ** (w - self.M)) for w in ws]
        sprung = [mp.mpf(0)] * len(ws)
        for _ in range(self.M):
            ys, s = self.T_pfad(ys)
            sprung = [max(a, b) for a, b in zip(sprung, s)]
        return ys, sprung


def verfolge(param, auswerten, tol=mp.mpf("0.8"), runden=12, log=None):
    """Adaptive Wegverfeinerung: auswerten(param_liste) → (werte, sprung); fügt Mittelpunkte ein, wo sprung > tol."""
    param = list(param)
    for r in range(runden):
        werte, sprung = auswerten(param)
        schlecht = [i for i in range(1, len(param)) if sprung[i] > tol]
        if log:
            log(f"   Verfeinerung {r}: {len(param)} Punkte, max Sprung {mp.nstr(max(sprung), 3)}, zu verfeinern: {len(schlecht)}"
                + (f" (erste Stelle {mp.nstr(param[schlecht[0]], 6)})" if schlecht else ""))
        if not schlecht or r == runden - 1:
            return param, werte, max(sprung)
        neu = list(param)
        for i in sorted(schlecht, reverse=True):
            neu.insert(i, (param[i - 1] + param[i]) / 2)
        param = neu
