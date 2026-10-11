"""Strang P — Satz P-V als wiederverwendbares Modul (rigoros, arb). Logik 1:1 aus p29_fK_koeffizienten.py übernommen
(dort als Closures in main()); p29 bleibt unverändert.

Kneser(prec): lädt p26_iter6.json (t̃, c̃), p27_zert.json, p28_zert.json, berechnet ε_T, ε_U wie p29 und stellt bereit:
  rep_ball, sexp_ball((zc, r)) → (Wert, Ableitung) als Bälle für sexp_K auf der Scheibe,
  sexp_num, slog_num (Numerik), slog_krawczyk(x, w0, xrad) → (Mitte, Radius) einer Scheibe, die slog_K(X) enthält
  (X = Scheibe um x mit Radius xrad; Eindeutigkeit von sexp_K(z) = x in der Scheibe für jedes x ∈ X).
"""
import json
import time
from pathlib import Path

from flint import acb, arb, ctx

import kc_basis as KB
import kc_rig

HIER = Path(__file__).resolve().parent


def amax(a, b):
    return (a + b + abs(a - b)) / 2


class Kneser:
    def __init__(self, prec=700, n0=300, m_extra=250, laut=True):
        t0 = time.time()
        ctx.prec = prec
        self.prec = prec
        Z7 = json.loads((HIER / "p27_zert.json").read_text(encoding="utf-8"))
        Z8 = json.loads((HIER / "p28_zert.json").read_text(encoding="utf-8"))
        self.Z7, self.Z8 = Z7, Z8
        self.H, self.RHO = arb(Z7["H"]), arb(Z7["rho"])
        d_max = arb(Z7["d_max"]).abs_upper()
        e_max = arb(Z7["e_max"]).abs_upper()
        phi_f = arb(Z8["phi_fehler"]).abs_upper()
        assert Z8["pruefung_A"] and Z8["pruefung_B"]
        R = kc_rig.Rig(prec, n0=n0, m_extra=m_extra)
        self.R = R
        d = json.loads((HIER / "p26_iter6.json").read_text(encoding="utf-8"))
        t = [acb(arb(arb(a).mid())) for a, b in d["t"]]
        c = [acb(arb(arb(a).mid()), arb(arb(b).mid())) for a, b in d["c"]]
        self.t, self.c = t, c
        self.dt = [n * t[n] for n in range(1, len(t))]
        dt = self.dt
        supT1 = arb(0)
        for i in range(12):
            for j in range(10):
                cz = acb(arb(-0.55 + 0.1 * i), arb(-0.45 + 0.1 * j))
                supT1 = amax(supT1, R.T_ball(dt, kc_rig.disc(cz, arb("0.075"))).abs_upper())
        tpi_ = 2 * arb.pi()
        supU1 = arb(0)
        for i in range(12):
            y = 0.35
            while y <= 3.05:
                zc = acb(arb(-0.55 + 0.1 * i), arb(y))
                a, tau = R.U_taylor(c, zc, 10, arb("0.075"), arb("0.2"))
                supU1 = amax(supU1, sum((j * abs(a[j]) * arb("0.075") ** (j - 1) for j in range(1, 11)), arb(0))
                             + tau / arb("0.025"))
                y += 0.1
        w_oben = arb(Z8["w_oben"]).abs_upper()
        assert w_oben < arb("0.1")
        P1 = 1 + sum((n * abs(R.p[n]) * arb("0.1") ** (n - 1) for n in range(2, R.n0 + 1)), arb(0)) \
            + R.C * arb("0.2") ** (R.n0 + 1) / arb("0.1")
        th1 = sum((tpi_ * kk * abs(ck) * (-tpi_ * 3 * kk).exp() for kk, ck in enumerate(c)), arb(0))
        supU1 = amax(supU1, P1 * w_oben * R.absL * (1 + th1))
        lip = amax(supT1, supU1) * arb("1.001")
        tn = abs(t[0] - 1)
        tau_b = 2 * (tn + lip * phi_f)
        self.epsT = lip * (phi_f + tau_b) + supT1 * d_max
        self.epsU = lip * (phi_f + tau_b) + supU1 * (d_max + e_max)
        etaB = arb(Z7["eta_B"]).abs_upper()
        self.eps = amax(self.epsT, self.epsU) + etaB
        self.hf = float(self.H.mid())
        self.w_oben = w_oben
        KB.setup(prec)
        self.tn_num = [KB.mid(x) for x in t]
        self.cn_num = [KB.mid(x) for x in c]
        if laut:
            print(f"Kneser: ε_T ≤ {self.epsT.str(3)}, ε_U ≤ {self.epsU.str(3)}, ε ≤ {self.eps.str(3)} ({time.time() - t0:.0f}s)",
                  flush=True)

    # ---------------------------------------------------------------- wie p29
    def rep_ball(self, S, m_tay=8):
        R, t, dt, c, eps = self.R, self.t, self.dt, self.c, self.eps
        zc, r = S
        im = float(zc.imag.mid())
        if abs(im) <= self.hf:
            Dk = kc_rig.disc(zc, r) if r > 0 else zc
            v = R.T_ball(t, Dk)
            dv = R.T_ball(dt, Dk)
        else:
            spiegel = im < 0
            zz = zc.conjugate() if spiegel else zc
            if r > arb("1e-100"):
                a, tau = R.U_taylor(c, zz, m_tay, r, arb("0.1"))
                v = a[0] + acb(arb(0, 1), arb(0, 1)) * (sum((abs(a[j]) * r ** j for j in range(1, m_tay + 1)), arb(0)) + tau)
                dv = a[1] + acb(arb(0, 1), arb(0, 1)) * (sum((j * abs(a[j]) * r ** (j - 1) for j in range(2, m_tay + 1)), arb(0))
                                                         + 2 * tau / r)
            else:
                v, dv = R.U_dU_ball(c, zz if r == 0 else kc_rig.disc(zz, r))
            if spiegel:
                v, dv = v.conjugate(), dv.conjugate()
        eb = acb(arb(0, eps), arb(0, eps))
        return v + eb, dv + eb * 100

    def sexp_ball(self, S):
        zc, r = S
        k = 0
        while float(zc.real.mid()) >= 0.5:
            zc -= 1
            k += 1
        while float(zc.real.mid()) < -0.5:
            zc += 1
            k -= 1
        v, dv = self.rep_ball((zc, r))
        for _ in range(k):
            v = v.exp()
            dv = dv * v
        for _ in range(-k):
            dv = dv / v
            v = v.log()
        return v, dv

    def sexp_num(self, w):
        k = 0
        while float(w.real.mid()) >= 0.5:
            w -= 1
            k += 1
        while float(w.real.mid()) < -0.5:
            w += 1
            k -= 1
        im = float(w.imag.mid())
        if abs(im) <= self.hf:
            v, dv = KB.T_dT(self.tn_num, w)
        elif im > 0:
            v, dv = KB.U_dU(self.cn_num, w)
        else:
            v, dv = KB.U_dU(self.cn_num, w.conjugate())
            v, dv = v.conjugate(), dv.conjugate()
        for _ in range(k):
            v = v.exp()
            dv = dv * v
        for _ in range(-k):
            dv = dv / v
            v = v.log()
        return KB.mid(v), KB.mid(dv)

    def slog_num(self, x, w):
        for _ in range(100):
            v, dv = self.sexp_num(w)
            dw = KB.mid((v - x) / dv)
            w = w - dw
            if abs(dw).mid() < arb(2) ** (-self.prec + 40):
                break
        return w

    def slog_krawczyk(self, x, w0, xrad=arb(0), r0=None):
        """Scheibe D(w̃, r0) ⊇ slog_K(D(x, xrad)) per 1D-Krawczyk (Eindeutigkeit in der Scheibe). Rückgabe (w̃, r0) oder (None, w̃)."""
        w = self.slog_num(x, w0)
        g0, dg = self.sexp_ball((w, arb(0)))
        g0 = g0 - x + acb(arb(0, xrad), arb(0, xrad))
        if r0 is None:
            r0 = 2 * g0.abs_upper() / dg.abs_lower() + arb(10) ** -300
        for _ in range(6):
            _, Dp = self.sexp_ball((w, r0))
            d0 = KB.mid(Dp)
            Kc = -g0 / d0
            Kr = abs(1 - Dp / d0) * r0
            if abs(Kc) + Kr < r0:
                return w, r0
            r0 = r0 * 2
        return None, w
