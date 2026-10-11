"""Strang P — Satz P-V, RIGOROSE Bausteine (Ballarithmetik, python-flint).

Rig(prec, n0, r0):
  * L als Ball (Lambert W), p_n (n ≤ n0) als Bälle per Rekursion, Restkonstante C mit |P − T_{n0}| ≤ C|w|^{n0+1} auf |w| ≤ r0
    (Fixpunktargument wie p25; Eindeutigkeit der formalen Lösung, L^n ≠ L).
  * P_ball(w): Einschluss von P(w) für einen acb-Ball w (Reduktion w/L^k bis |w| ≤ r_ausw, Horner + Restball, k-mal exp).
  * U_taylor(c, zc, m, r, R): Taylor-Modell von Ũ(zc + h) = P(exp(L(zc + h + θ̃(zc + h)))) auf |h| ≤ r:
        Koeffizienten a_0..a_m (Bälle) und Restschranke τ mit |Ũ(zc+h) − Σ a_j h^j| ≤ τ.
    Rest = Cauchy-Schwanz (grobe Schranke von F₁ auf |h| = R, F₁ = exp^k∘T_{n0}∘(w/L^k)) + Wert-Fehler aus P − T_{n0}.
  * T_taylor(t, zc, r): exakte verschobene Koeffizienten des Polynoms T̃ (endlich).
  * expT_taylor(t, zc, m, r, R): Taylor-Modell von exp(T̃(zc + h)).
Alle Schranken sind obere Schranken als arb (Vergleiche immer mit arb-Semantik).
"""
from flint import acb, acb_series, arb, ctx


def ob(x):
    return x.abs_upper()


def disc(center, rad):
    """acb-Ball, der die Kreisscheibe |z − center| ≤ rad enthält."""
    return center + acb(arb(0, rad), arb(0, rad))


class Rig:
    def __init__(self, prec, n0=400, r0="0.5", r_ausw="0.12", m_extra=300):
        ctx.prec = prec
        self.prec = prec
        L = -acb(-1).lambertw(-1)
        assert abs(L - L.exp()) < arb(10) ** (-(prec // 4))
        self.L = L
        self.absL = abs(L)
        n0 = int(n0)
        # Koeffizienten (Bälle, keine Mittelpunkt-Abschneidung)
        p = [L, acb(1)] + [acb(0)] * (n0 - 1)
        a = [acb(0), acb(1)] + [acb(0)] * (n0 - 1)
        e = [acb(1), acb(1)] + [acb(0)] * (n0 - 1)
        Lp = L
        for n in range(2, n0 + 1):
            Lp = Lp * L
            s = acb(0)
            for k in range(1, n):
                s += k * a[k] * e[n - k]
            q = s / n
            pn = L * q / (Lp - L)
            p[n] = pn
            a[n] = pn
            e[n] = q + pn
        self.p = p
        self.n0 = n0
        self.r0 = arb(r0)
        self.r_ausw = arb(r_ausw)
        self.C = self._restkonstante(m_extra)

    def _exp_reihe(self, a, n):
        e = [acb(1)] + [acb(0)] * n
        for mm in range(1, n + 1):
            s = acb(0)
            for k in range(1, min(mm, len(a) - 1) + 1):
                s += k * a[k] * e[mm - k]
            e[mm] = s / mm
        return e

    def _restkonstante(self, m_extra):
        """|P − T_{n0}| ≤ C|w|^{n0+1} auf |w| ≤ r0 (Kontraktionsargument, siehe p25)."""
        L, p, n0, r0 = self.L, self.p, self.n0, self.r0
        absL = self.absL
        nmax = n0 + m_extra
        b = [acb(0)] + [p[n] / L ** n for n in range(1, n0 + 1)] + [acb(0)] * m_extra
        Eb = self._exp_reihe(b, nmax)
        delta = [L * Eb[k] - (p[k] if k <= n0 else 0) for k in range(nmax + 1)]
        R1 = arb(2)
        S1 = sum((ob(p[n]) * (R1 / absL) ** n for n in range(1, n0 + 1)), arb(0))
        S2 = absL + sum((ob(p[n]) * R1 ** n for n in range(1, n0 + 1)), arb(0))
        B1 = absL * S1.exp() + S2
        D = sum((ob(delta[k]) * r0 ** (k - n0 - 1) for k in range(n0 + 1, nmax + 1)), arb(0))
        D += B1 * R1 ** (-(n0 + 1)) * (r0 / R1) ** m_extra / (1 - r0 / R1)
        S1r = sum((ob(p[n]) * (r0 / absL) ** n for n in range(1, n0 + 1)), arb(0))
        M1 = absL * S1r.exp()
        fakt = M1 * absL ** (-(n0 + 1))
        C = D / (1 - fakt) * 2
        for _ in range(40):
            if D + fakt * C * (C * (r0 / absL) ** (n0 + 1)).exp() < C:
                break
            C = C * 2
        lip = fakt * (C * (r0 / absL) ** (n0 + 1)).exp()
        assert D + fakt * C * (C * (r0 / absL) ** (n0 + 1)).exp() < C and lip < 1
        self.lip = lip
        return C

    # ---------------------------------------------------------------- Punkt-/Ball-Auswertung
    def _horner(self, coeffs, w):
        s = acb(0)
        for c in reversed(coeffs):
            s = s * w + c
        return s

    def P_ball(self, w):
        k = 0
        while not (ob(w) <= self.r_ausw):
            w = w / self.L
            k += 1
            if k > 200:
                raise RuntimeError("P_ball: Reduktion versagt")
        rest = self.C * ob(w) ** (self.n0 + 1)
        s = self._horner(self.p, w) + acb(arb(0, rest), arb(0, rest))
        for _ in range(k):
            s = s.exp()
        return s

    def P_dP_ball(self, w):
        """Einschluss von (P(w), P'(w)); Rest-Ableitung per Cauchy auf dem Kreis vom Radius |w| (2|w| ≤ r0)."""
        k = 0
        while not (ob(w) <= self.r_ausw):
            w = w / self.L
            k += 1
        s, ds = acb(0), acb(0)
        for cf in reversed(self.p):
            ds = ds * w + s
            s = s * w + cf
        aw = ob(w)
        rest = self.C * aw ** (self.n0 + 1)
        drest = self.C * (2 * aw) ** (self.n0 + 1) / aw if aw > 0 else arb(0)
        s = s + acb(arb(0, rest), arb(0, rest))
        ds = ds + acb(arb(0, drest), arb(0, drest))
        fak = acb(1)
        for _ in range(k):
            s = s.exp()
            fak = fak * s
        return s, ds * fak / self.L ** k

    def U_dU_ball(self, c, z):
        I = acb(0, 1)
        tpi = 2 * acb.pi() * I
        q = (tpi * z).exp()
        th, dth = acb(0), acb(0)
        for k in range(len(c) - 1, -1, -1):
            th = th * q + c[k]
            dth = dth * q + k * c[k]
        dth = dth * tpi
        w = (self.L * (z + th)).exp()
        f, df = self.P_dP_ball(w)
        return f, df * w * self.L * (1 + dth)

    def theta_ball(self, c, z):
        q = (2 * acb.pi() * acb(0, 1) * z).exp()
        return self._horner(c, q)

    def U_ball(self, c, z):
        return self.P_ball((self.L * (z + self.theta_ball(c, z))).exp())

    def T_ball(self, t, z):
        return self._horner(t, z)

    # ---------------------------------------------------------------- Taylor-Modelle
    def _grob_sup_kreis(self, f, zc, R, nbox=48):
        """obere Schranke von |f(zc + h)| auf |h| = R (Überdeckung durch nbox Bälle)."""
        pi = arb.pi()
        rb = R * pi / nbox * arb("1.05")
        best = arb(0)
        for j in range(nbox):
            h = R * (2 * pi * acb(0, 1) * j / nbox).exp()
            v = ob(f(disc(zc + h, rb)))
            best = v if j == 0 else (best + v + abs(best - v)) / 2
        return best

    def U_taylor(self, c, zc, m, r, R):
        ctx.cap = m + 1
        L = self.L
        I = acb(0, 1)
        tpi = 2 * acb.pi() * I
        qc = (tpi * zc).exp()
        # θ̃(zc + h) = Σ_k c_k qc^k e^{2πikh}
        th = [acb(0)] * (m + 1)
        qk = acb(1)
        for k, ck in enumerate(c):
            g = ck * qk
            fac = acb(1)
            tk = tpi * k
            for j in range(m + 1):
                th[j] += g * fac
                fac = fac * tk / (j + 1)
            qk = qk * qc
        u = [L * (zc + th[0]), L * (1 + th[1])] + [L * th[j] for j in range(2, m + 1)]
        w = acb_series(u).exp()
        # k so wählen, dass w/L^k auf |h| ≤ R in |·| ≤ r0/2 liegt (grobe Schranke)
        Wsup = self._grob_sup_kreis(lambda z: (L * (z + self.theta_ball(c, z))).exp(), zc, R)
        Wsup0 = ob((L * (zc + self.theta_ball(c, zc))).exp())
        Wmax = (Wsup + Wsup0 + abs(Wsup - Wsup0)) / 2       # max-Prinzip: sup auf der Scheibe ≤ sup auf dem Rand
        k = 0
        while not (Wmax / self.absL ** k <= self.r0 / 2):
            k += 1
        Lk = L ** k
        om = [x / Lk for x in w.coeffs()] + [acb(0)] * (m + 1 - len(w.coeffs()))
        om0 = om[0]
        # verschobenes Polynom T_{n0}(om0 + y) = Σ_j b_j y^j, j ≤ m
        b = [acb(0)] * (m + 1)
        # Horner-Verschiebung: b = Koeffizienten von Σ p_n (om0 + y)^n
        shift = [acb(0)] * (m + 1)
        for n in range(self.n0, -1, -1):
            # shift := shift * (om0 + y) + p_n
            neu = [acb(0)] * (m + 1)
            for j in range(m + 1):
                neu[j] = shift[j] * om0 + (shift[j - 1] if j >= 1 else 0)
            neu[0] += self.p[n]
            shift = neu
        b = shift
        y = acb_series([acb(0)] + om[1:])
        F = acb_series([b[m]])
        for j in range(m - 1, -1, -1):
            F = F * y + b[j]
        for _ in range(k):
            F = F.exp()
        a = F.coeffs() + [acb(0)] * (m + 1 - len(F.coeffs()))
        # Cauchy-Schwanz von F₁ (T_{n0}-Komposition) auf |h| = R
        def F1(z):
            ww = (L * (z + self.theta_ball(c, z))).exp() / Lk
            s = self._horner(self.p, ww)
            for _ in range(k):
                s = s.exp()
            return s
        MR = self._grob_sup_kreis(F1, zc, R)
        q = r / R
        schwanz = MR * q ** (m + 1) / (1 - q)
        # Wert-Fehler P − T_{n0}: |Δ| ≤ C|ω|^{n0+1}, Lipschitz von exp^k über grobe Bälle
        Wr = ob((L * (disc(zc, r) + self.theta_ball(c, disc(zc, r)))).exp() / Lk)
        eps = self.C * Wr ** (self.n0 + 1)
        A = self._horner(self.p, (L * (disc(zc, r) + self.theta_ball(c, disc(zc, r)))).exp() / Lk)
        A = A + acb(arb(0, eps), arb(0, eps))
        lipk = arb(1)
        for _ in range(k):
            Ae = A.exp()
            lipk = lipk * ob(Ae)
            A = Ae
        return a, schwanz + eps * lipk

    def T_taylor(self, t, zc):
        """alle Koeffizienten von T̃(zc + h) (exakt, endlich)."""
        n = len(t)
        shift = [acb(0)] * n
        for i in range(n - 1, -1, -1):
            neu = [acb(0)] * n
            for j in range(n):
                neu[j] = shift[j] * zc + (shift[j - 1] if j >= 1 else 0)
            neu[0] += t[i]
            shift = neu
        return shift

    def expT_taylor(self, t, zc, m, r, R):
        ctx.cap = m + 1
        tc = self.T_taylor(t, zc)
        E = acb_series(tc[:m + 1]).exp()
        a = E.coeffs() + [acb(0)] * (m + 1 - len(E.coeffs()))
        MR = self._grob_sup_kreis(lambda z: self.T_ball(t, z).exp(), zc, R)
        q = r / R
        return a, MR * q ** (m + 1) / (1 - q)
