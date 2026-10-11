"""mk_kneser_arb — Kugel-Einschlüsse für die Kneser-Halbexponentialfunktion in z = x/ρ, F(z) = f_K(ρz)/ρ = (1/ρ)Σ b_k z^k.

SCHNITTSTELLE ZU STRANG P (alles Weitere wird hier streng abgeleitet):
  (K1) b_k = a_k ρ^k für k ≤ K_daten als Kugeln (Mitte, Radius),
  (K2) Schwanzschranke |b_k| ≤ C·k^{−p} für alle k > K_daten (p > 3).
Quellen:
  quelle="vorzertifikat": Mitten aus gemeinsam/kneser_halbexp_koeff_4000.json (Strang A, NICHT zertifiziert), Radien konservativ
      1e-15·|b_k| + 1e-40, Schwanz C = 100, p = 6.4 für k > 4000 (beobachtet max_{k≥500} |b_k|k^{6.45} = 74.2) — ANNAHME, nicht bewiesen.
  quelle=<Pfad>: JSON von Strang P mit Feldern {"K": int, "b_mitte": [str], "b_radius": [str], "C": str, "p": str} — dann streng.
Abgeleitete Grössen (streng relativ zu (K1)+(K2)):
  F, F' an Punkten |z0| ≤ 1 (Horner mit K_eval Gliedern + Restsumme), an Kugeln in D̄ (Mittelwertsatz mit S1 = sup|F'|, S2 = sup|F''|),
  auf Kreisbögen (Winkelkugel), h = F⁻¹ (parametrischer Krawczyk in D̄), Fortsetzung nahe Γin über F(z) = log(ρF(w))/ρ, w = e^{ρz}/ρ.
"""
import json
import os

from flint import acb, arb, ctx

HIER = os.path.dirname(os.path.abspath(__file__))
GEM = os.path.join(HIER, "..", "gemeinsam")


def _summe_potenz(C, p, K, j):
    """Obere Schranke für Σ_{k>K} k^j·C·k^{−p} ≤ C·∫_K^∞ x^{j−p}dx = C K^{j−p+1}/(p−j−1) (p − j > 1)."""
    return C * arb(K) ** (j - p + 1) / (p - j - 1)


class KneserArb:
    def __init__(self, quelle="vorzertifikat", K_eval=400, prec=120):
        ctx.prec = prec
        self.prec = prec
        L = (-acb(-1).lambertw()).conjugate()      # obere Fixstelle von exp, streng (arb): L = conj(−W₀(−1))
        self.L = L
        self.rho = abs(L)
        self.th = L.arg()
        self.mu = L.sqrt()
        self.ell = self.mu.log()
        self.eith = acb(self.th.cos(), self.th.sin())
        if quelle == "vorzertifikat":
            d = json.load(open(os.path.join(GEM, "kneser_halbexp_koeff_4000.json"), encoding="utf-8"))
            b = []
            for k, x in enumerate(d["a"]):
                m = arb(x) * self.rho ** k
                b.append(arb(m.mid(), (abs(m).upper() * arb("1e-15") + arb("1e-40")).upper()))
            self.C, self.p = arb(100), arb("6.4")
            self.zertifiziert = False
        elif quelle.startswith("pzert"):
            # Strang P, Satz P-V: gemeinsam/kneser_halbexp_koeff_zert_640.json, Eintrag n = [Mitte, |a_n − Mitte| ≤ e_n, e_n ρ^n].
            # (K1) streng. (K2) noch NICHT geliefert: Schwanz |b_k| ≤ C k^{−p} für k > 640 als Parameter: "pzert" (C = 100, p = 6.4)
            # oder "pzert:C:p". zertifiziert = "K1" (Ergebnis streng bis auf (K2)).
            teile = quelle.split(":")
            self.C = arb(teile[1]) if len(teile) > 1 else arb(100)
            self.p = arb(teile[2]) if len(teile) > 2 else arb("6.4")
            d = json.load(open(os.path.join(GEM, "kneser_halbexp_koeff_zert_640.json"), encoding="utf-8"))
            b = []
            for n, (mi, er, _) in enumerate(d["a"]):
                rn = self.rho ** n
                m = arb(mi) * rn
                # Ps Radien e_n enthalten die Rundung der 140-stelligen Dezimalmitten nicht (bis 5e-141 bei a_0; m27): addiert.
                ziff = len(mi.replace("-", "").replace(".", "").lstrip("0").split("e")[0].split("E")[0])
                dez = abs(m) * arb(10) ** (-(ziff - 2))
                b.append(arb(m.mid(), (m.rad() + arb(er) * rn + dez).upper() * (1 + arb("1e-10")).upper()))
            self.zertifiziert = "K1 streng, K2 angenommen (C = %s, p = %s)" % (self.C, self.p)
        else:
            d = json.load(open(quelle, encoding="utf-8"))
            b = []
            for mi, ra in zip(d["b_mitte"], d["b_radius"]):
                m = arb(mi)
                # Dezimalrundung der Mitte (falls nicht im Radius enthalten) und Aufrundung des Radius-Strings zur Sicherheit
                ziff = len(mi.replace("-", "").replace(".", "").lstrip("0").split("e")[0].split("E")[0])
                dez = abs(m) * arb(10) ** (-(ziff - 2))
                b.append(arb(m.mid(), (m.rad() + arb(ra) * (1 + arb("1e-3")) + dez).upper()))
            self.C, self.p = arb(d["C"]), arb(d["p"])
            self.zertifiziert = True
        self.Kd = len(b) - 1
        self.K_eval = min(K_eval, self.Kd)
        self.b = b
        Ke = self.K_eval
        # Restsummen über K_eval < k ≤ Kd (aus den Daten) + Schwanz k > Kd (aus (K2))
        def rest(j):
            s_ = arb(0)
            for k in range(Ke + 1, self.Kd + 1):
                s_ += abs(b[k]) * arb(k) ** j if j > 0 else abs(b[k])
            return s_ + _summe_potenz(self.C, self.p, self.Kd, j)
        self.R0 = arb(rest(0).upper())
        self.R1 = arb((rest(1)).upper())
        # Supremumsschranken auf D̄
        def summe(j):
            s_ = arb(0)
            for k in range(self.Kd + 1):
                f = arb(1)
                for i in range(j):
                    f *= (k - i)
                s_ += abs(b[k]) * f
            # Schwanz Σ k(k−1)…(k−j+1)|b_k| ≤ Σ k^j |b_k|
            return s_ + _summe_potenz(self.C, self.p, self.Kd, j)
        self.S1 = arb((summe(1) / self.rho).upper())       # sup|F'| auf D̄
        self.S2 = arb((summe(2) / self.rho).upper())       # sup|F''|
        self.S3 = arb((summe(3) / self.rho).upper())       # sup|F'''|
        # exakte Mitten (Radius 0) für die Horner-Auswertung; Koeffizientenradien getrennt als reelle Summen
        self.bk = [acb(arb(x.mid())) for x in b[:Ke + 1]]
        self.bk1 = [acb(arb((k * b[k]).mid())) for k in range(1, Ke + 1)]
        self.Rc0 = arb(sum((x.rad() for x in b[:Ke + 1]), arb(0)).upper())
        self.Rc1 = arb(sum((arb(k) * b[k].rad() for k in range(1, Ke + 1)), arb(0)).upper())
        self.prec_horner = 400        # Rundungsradius 2^-400 × Umwickelfaktor ≤ √2^K_eval bleibt vernachlässigbar

    # --- Punkte und Kugeln in D̄
    @staticmethod
    def _exakt(z0):
        return acb(arb(z0.real.mid()), arb(z0.imag.mid()))

    def _horner(self, koeff, z0):
        """Horner mit exaktem Punkt (Radius 0) und exakten Koeffizienten bei hoher Präzision: die Rechteck-Umwicklung
        (Faktor ≤ √2 je Schritt) wirkt nur auf Rundungsradien ~2^-400."""
        alt = ctx.prec
        ctx.prec = self.prec_horner
        z = self._exakt(z0)
        assert bool(abs(z) <= 1), "Auswertepunkt ausserhalb D̄"
        s = acb(0)
        for c in reversed(koeff):
            s = s * z + c
        ctx.prec = alt
        return s

    def F_punkt(self, z0):
        e = self.Rc0 + self.R0
        s = self._horner(self.bk, z0)
        return (s + acb(arb(0, e.upper()), arb(0, e.upper()))) / self.rho

    def Fp_punkt(self, z0):
        e = self.Rc1 + self.R1
        s = self._horner(self.bk1, z0)
        return (s + acb(arb(0, e.upper()), arb(0, e.upper()))) / self.rho

    @staticmethod
    def _mitte_radius(z):
        z0 = acb(z.real.mid(), z.imag.mid())
        r = (z.real.rad() ** 2 + z.imag.rad() ** 2).sqrt()
        return z0, arb(r.upper())

    def in_scheibe(self, z):
        return bool(abs(z) <= 1)

    def F_kugel(self, z):
        """F über eine Kugel z ⊂ D̄ (geprüft)."""
        assert self.in_scheibe(z), "Kugel nicht in D̄"
        z0, r = self._mitte_radius(z)
        v = self.F_punkt(z0)
        e = r * self.S1
        return v + acb(arb(0, e.upper()), arb(0, e.upper()))

    def Fp_kugel(self, z):
        assert self.in_scheibe(z), "Kugel nicht in D̄"
        z0, r = self._mitte_radius(z)
        v = self.Fp_punkt(z0)
        e = r * self.S2
        return v + acb(arb(0, e.upper()), arb(0, e.upper()))

    def F_kreis(self, phi):
        """F(e^{iφ}) für eine Winkelkugel φ (Bogen ⊂ ∂D): Mittelwertsatz längs des Bogens (Länge ≤ 2 rad φ)."""
        z0, r = self._kreispunkt(phi)
        e = r * self.S1
        return self.F_punkt(z0) + acb(arb(0, e.upper()), arb(0, e.upper()))

    def _kreispunkt(self, phi):
        """Exakter Punkt z0 knapp innerhalb des Kreises nahe e^{i mid φ} und Radius r ≥ |z − z0| für alle z = e^{iφ}, φ ∈ Kugel."""
        p0 = arb(phi.mid())
        z = acb(p0.cos(), p0.sin()) * (1 - arb(2) ** -100)
        z0 = self._exakt(z)
        r = arb(phi.rad()) + arb(2) ** -98
        return z0, arb(r.upper())

    def Fp_kreis(self, phi):
        z0, r = self._kreispunkt(phi)
        e = r * self.S2
        return self.Fp_punkt(z0) + acb(arb(0, e.upper()), arb(0, e.upper()))

    # --- Umkehrzweig h (parametrischer Krawczyk)
    def h(self, w, v0, r0=None):
        """Einschluss von h(w) für eine Kugel w (alle Punkte), Startpunkt v0 (komplex, float-genau).
        Krawczyk: X = v0 + Kugel(r0), K(X) = v0 − Y(F(v0) − w) + (1 − Y F'(X))(X − v0) ⊂ int X ⇒ für jedes w in der Kugel genau eine
        Lösung in X. Verlangt X ⊂ D̄. Rückgabe: K(X) (enthält die Lösungen) oder None."""
        v0 = acb(v0)
        w0, rw = self._mitte_radius(w)
        Fp0 = self.Fp_punkt(v0)
        Y = 1 / acb(Fp0.real.mid(), Fp0.imag.mid())
        F0 = self.F_punkt(v0)
        if r0 is None:
            r0 = arb(4) * (abs(Y) * (abs(F0 - w0) + rw)).upper() + arb("1e-25")
        for _ in range(6):
            X = v0 + acb(arb(0, r0.upper()), arb(0, r0.upper()))
            if not self.in_scheibe(X):
                return None
            KX = v0 - Y * (F0 - w) + (1 - Y * self.Fp_kugel(X)) * (X - v0)
            innen = bool(abs(KX.real - v0.real).upper() < r0) and bool(abs(KX.imag - v0.imag).upper() < r0)
            if innen:
                return KX
            r0 = r0 * 2
        return None

    # --- Fortsetzung nahe Γin (auch ausserhalb D̄): F(z) = log(ρF(w))/ρ, F'(z) = F'(w)·w/F(w), w = e^{ρz}/ρ
    def F_gin(self, z):
        w = (self.rho * z).exp() / self.rho
        Fw = self.F_kugel(w)
        return (self.rho * Fw).log() / self.rho, w, Fw

    def Fp_gin(self, z):
        w = (self.rho * z).exp() / self.rho
        return self.Fp_kugel(w) * w / self.F_kugel(w)


if __name__ == "__main__":
    import sys
    import time
    sys.stdout.reconfigure(encoding="utf-8")
    t0 = time.time()
    K = KneserArb()
    print(f"Quelle: Vorzertifikat (zertifiziert = {K.zertifiziert}); K_daten = {K.Kd}, K_eval = {K.K_eval}; "
          f"R0 ≤ {float(K.R0):.2e}, R1 ≤ {float(K.R1):.2e}; S1 ≤ {float(K.S1):.5f}, S2 ≤ {float(K.S2):.5f}, S3 ≤ {float(K.S3):.5f} "
          f"({time.time()-t0:.1f}s)")
    print("F(e^{iθ}) − e^{iθ} ∈", K.F_punkt(K.eith) - K.eith)
    print("F'(e^{iθ}) − √L ∈", K.Fp_punkt(K.eith) - K.mu)
    w = acb(arb("0.6"), arb("0.5"))
    import numpy as np
    from mk_kneser import h_newton
    v0 = complex(h_newton(np.array([complex(0.6, 0.5)]))[0])
    print("h(0.6+0.5i) ∈", K.h(w, v0), " (float", v0, ")")
    z = acb(arb("-1.05"), arb("0.1"))
    print("F_gin(−1.05+0.1i) ∈", K.F_gin(z)[0])
    from mk_kneser import F as Ff
    print("  float:", Ff(np.array([complex(-1.05, 0.1)]))[0])
