"""Strang P, Prüfung 32 (Meilenstein 6, Punkt 2): f_K auf ∂D rigoros — θ-Ableitungen bis 6, (K2), (K3).

Ziel: ‖∂_θ⁶ f_K(ρe^{iθ})‖_{L¹(0,2π)} ≤ 2π·C₆  ⇒  |b_k| ≤ C₆·k^{−6} für alle k ≥ 1 (sechsmal partiell integrieren; f_K ∈ C⁵(∂D),
∂_θ⁵f absolut stetig, ∂_θ⁶f ~ |θ − θ_L|^{−0.553} integrierbar). Symmetrie f(x̄) = conj f(x) ⇒ ∫_0^{2π} = 2∫_0^π.

Zellen [θ_a, θ_b] ⊂ [0, π]; Bogenkugel X ⊇ Bogen; Z = Krawczyk-Scheibe ⊇ slog_K(X) (kc_kneser.slog_krawczyk); dann
  T-Gebiet (|Im Z| ≤ 0.42): z-Parameter, x = S(z), f = S(z+½), ∂_θ = (iS/S')∂_z; Jets von S aus T̃ (± ε/δʲ) und exp/log-Verschiebung.
  Strukturgebiet (Im Z ≥ 0.39 > y₁ = 0.36): EXAKTE Form (p32a): sexp_K(z) = P(e^{L(z+θ_K(z))}) ⇒ mit w = e^{L(z+θ_K)}, q = e^{2πiz}:
      x = P(w),  f = P(H(w)),  H(w) = √L·w·R(q),  R(q) = exp(L(θ_K(−q) − θ_K(q))),
      H^{(k)}(w) = √L w^{1−k} G_k(q),  G_0 = R,  G_k = (2 − k)G_{k−1} + 𝒟G_{k−1},  𝒟 = γσ(q)d/dq,  σ = q/(1 + 2πi q θ_K'(q)),
      ∂_θ = (iP(w)/P'(w)) d/dw.  Für k ≥ 2 ist G_k = O(q) (Faktor σ) ⇒ keine Auslöschung, Kugelauswertung trägt bis w → 0.
      θ_K = θ̃ + Δ, |[ηʲ]Δ(q₀+η)| ≤ δθ e^{2πj y₁}/(1−t)^{j+1}, t = e^{2πy₁}|q₀| (p32a).
  Rest |θ − θ_L| < δ_min: |∂_θ⁶f| ≤ B₀ + c₆|H₆|, |H₆| ≤ √L(M₆/q_t)C_q|w|^{Reγ−5}/6! (Schwarz-Lemma für G₆), integriert exakt.
Ausgabe: p32_rand.json (C6, Teilintegrale, Zellen, f/f' an L ± δ, sup|∂⁵f| usw.), p32_zellen.json (Zellenliste).
Aufruf: python p32_rand_ableitungen.py SEGMENT [delta_min_exp] [prec]   SEGMENT ∈ {rechts, links, test}
"""
import json
import math
import sys
import time
from pathlib import Path

from flint import acb, acb_series, arb, ctx

import kc_rig
from kc_kneser import Kneser, amax

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
SEG = sys.argv[1] if len(sys.argv) > 1 else "test"
DMIN_EXP = int(sys.argv[2]) if len(sys.argv) > 2 else 40
PREC = int(sys.argv[3]) if len(sys.argv) > 3 else 700
M = 6                      # höchste θ-Ableitung


def ob(x):
    return x.abs_upper()


def disc(c, r):
    return kc_rig.disc(c, r)


class Rand:
    def __init__(self):
        self.K = Kneser(PREC)
        K = self.K
        R = K.R
        self.R = R
        self.L = R.L
        self.sqrtL = (R.L / 2).exp()                      # L^{1/2} = e^{L/2} (log L = L)
        self.gam = 2 * acb.pi() * acb(0, 1) / R.L
        Z = json.loads((HIER / "p32a_theta_K.json").read_text(encoding="utf-8"))
        self.y1 = arb(Z["y1"])
        self.dth = arb(Z["delta_theta"])
        self.e2y1 = (2 * arb.pi() * self.y1).exp()
        self.c = K.c
        self.t = K.t
        self.eps = K.eps

    # ------------------------------------------------------------ Hilfen
    @staticmethod
    def shift(coeffs, z, m):
        """[h^j] von Σ coeffs[n] (z+h)^n, j ≤ m (Horner-Verschiebung, abgeschnitten)."""
        s = [acb(0)] * (m + 1)
        for n in range(len(coeffs) - 1, -1, -1):
            neu = [s[j] * z + (s[j - 1] if j >= 1 else 0) for j in range(m + 1)]
            neu[0] += coeffs[n]
            s = neu
        return s

    def _P_jet_direkt(self, y, m):
        """[ε^j] P(y + ε), j ≤ m, für |y| ≤ 0.12 (Ball): Taylor T_{n0} verschoben + Rest per Cauchy (Radius 0.12)."""
        R = self.R
        b = self.shift(R.p, y, m)
        rc = arb("0.12")
        rest = R.C * (ob(y) + rc) ** (R.n0 + 1)
        assert ob(y) + rc <= R.r0
        return [b[j] + disc(acb(0), rest / rc ** j) for j in range(m + 1)]

    def _P_jet_punkt(self, y, m):
        """[ε^j] P(y + ε) an einem PUNKT y (Reduktion y/L^k, dann k-mal exp als Reihe; Punktarithmetik ⇒ scharf)."""
        R = self.R
        k = 0
        while not (ob(y) <= R.r_ausw):
            y = y / R.L
            k += 1
            if k > 300:
                raise RuntimeError("P_jet: Reduktion")
        b = self._P_jet_direkt(y, m)
        if k == 0:
            return b
        Lk = R.L ** k
        ctx.cap = m + 1
        s = acb_series([b[j] / Lk ** j for j in range(m + 1)])
        for _ in range(k):
            s = s.exp()
        cf = s.coeffs()
        return cf + [acb(0)] * (m + 1 - len(cf))

    def P_jet(self, y, m):
        """[ε^j] P(y' + ε), j ≤ m, gültig für ALLE y' im Ball y.
        |y| ≤ 0.12: direkt. Sonst Taylor-Modell um den Mittelpunkt y_c: Punkt-Jets bis M2 = m + 48, Rest τ ≤ M_R·Σ_{i>M2}(1/4)^i
        (M_R = grobes sup|P| auf |ε| = R = 16 r_y, Cauchy auf R_P = R/4), Jets über den Ball per Verschiebung + τ/(R_P − r_y)^j."""
        if ob(y) <= arb("0.12"):
            return self._P_jet_direkt(y, m)
        R = self.R
        yc = acb(y.real.mid(), y.imag.mid())
        ry = (y.real.rad() ** 2 + y.imag.rad() ** 2).sqrt() * arb("1.0001") + arb(10) ** (-(PREC // 3))
        M2 = m + 48
        ph = self._P_jet_punkt(yc, M2)
        Rk = amax(16 * ry, abs(yc).lower() / 20)              # nicht winzig (sonst τ/R_P^j riesig)
        MR = R._grob_sup_kreis(R.P_ball, yc, Rk)
        RP = Rk / 4
        tau = MR * (arb(1) / 4) ** (M2 + 1) * arb(4) / 3
        dlt = disc(acb(0), ry)
        aus = []
        for j in range(m + 1):
            s = acb(0)
            for i in range(M2, j - 1, -1):
                s = s * dlt + ph[i] * math.comb(i, j)
            aus.append(s + disc(acb(0), tau / (RP - ry) ** j))
        return aus

    def theta_jets(self, q0, m):
        """[η^j] θ_K(q0 + η), j ≤ m: θ̃ exakt + Δ-Schranke."""
        tj = self.shift(self.c, q0, m)
        t = self.e2y1 * ob(q0)
        assert t < 1
        return [tj[j] + disc(acb(0), self.dth * self.e2y1 ** j / (1 - t) ** (j + 1)) for j in range(m + 1)]

    def G_werte(self, q0, m=M):
        """G_k(q0), k = 0..m (Kugeln)."""
        ctx.cap = m + 2
        tp = self.theta_jets(q0, m + 1)
        tm = self.theta_jets(-q0, m + 1)
        thp = acb_series(tp)
        thm = acb_series([tm[j] * (-1) ** j for j in range(m + 2)])     # θ_K(−q0 − η)
        Rs = (self.L * (thm - thp)).exp()
        dth = thp.derivative()
        qs = acb_series([q0, acb(1)])
        sig = qs / (1 + 2 * acb.pi() * acb(0, 1) * qs * dth)
        G = Rs
        werte = [G.coeffs()[0]]
        for k in range(1, m + 1):
            G = (2 - k) * G + self.gam * sig * G.derivative()
            werte.append(G.coeffs()[0])
        return werte, Rs.coeffs()[0]

    def theta_d(self, X, F, m=M):
        """∂_θ^j f, j = 0..m, aus Reihen X(ε) = x, F(ε) = f (gleiche Variable): ∂_θ = (iX/X') d/dε."""
        I = acb(0, 1)
        A = I * X / X.derivative()
        aus = [F.coeffs()[0]]
        for _ in range(m):
            F = A * F.derivative()
            aus.append(F.coeffs()[0])
        return aus

    # ------------------------------------------------------------ Strukturgebiet
    def struktur(self, Z, m=M, H6_null=False):
        """Z: acb-Ball (z-Werte). Rückgabe [∂_θ^j f] j ≤ m als Kugeln."""
        tpiI = 2 * acb.pi() * acb(0, 1)
        q0 = (tpiI * Z).exp()
        th0 = self.theta_jets(q0, 0)[0]
        w0 = (self.L * (Z + th0)).exp()
        Gw, R0 = self.G_werte(q0, m)
        fak = [1]
        for k in range(1, m + 1):
            fak.append(fak[-1] * k)
        Hk = [self.sqrtL * w0 * R0] + [self.sqrtL * w0 ** (1 - k) * Gw[k] / fak[k] for k in range(1, m + 1)]
        return self._aus_H(w0, Hk, m)

    def _aus_H(self, w0, Hk, m):
        ctx.cap = m + 1
        Pw = self.P_jet(w0, m)
        PH = self.P_jet(Hk[0], m)
        X = acb_series(Pw)
        Hs = acb_series([acb(0)] + Hk[1:])
        F = acb_series(PH)(Hs)
        return self.theta_d(X, F, m), Pw, PH

    # ------------------------------------------------------------ T-Gebiet
    def S_jets(self, Zc, Zr, m):
        """[h^j] sexp_K(z + h) für z ∈ D(Zc, Zr), j ≤ m (Kugeln). Zc Punkt-acb, Zr arb."""
        k = 0
        while float(Zc.real.mid()) >= 0.5:
            Zc -= 1
            k += 1
        while float(Zc.real.mid()) < -0.5:
            Zc += 1
            k -= 1
        im = abs(Zc.imag) + Zr
        re = abs(Zc.real) + Zr
        delta = amax(arb(0), arb("0.6") - re)
        d2 = arb("0.45") - im
        delta = (delta + d2 - abs(delta - d2)) / 2
        assert delta > arb("0.005"), "T-Jets: zu nah am Rand"
        Zb = disc(Zc, Zr) if Zr > 0 else Zc
        b = self.shift(self.t, Zb, m)
        b = [b[j] + disc(acb(0), self.eps / delta ** j) for j in range(m + 1)]
        ctx.cap = m + 1
        s = acb_series(b)
        for _ in range(k):
            s = s.exp()
        for _ in range(-k):
            s = s.log()
        cf = s.coeffs()
        return cf + [acb(0)] * (m + 1 - len(cf))

    def t_gebiet(self, Zc, Zr, m=M):
        ctx.cap = m + 1
        X = acb_series(self.S_jets(Zc, Zr, m))
        F = acb_series(self.S_jets(Zc + arb("0.5"), Zr, m))
        return self.theta_d(X, F, m)


# ---------------------------------------------------------------- Zellen-Schleife
def segment(RD, seite, dmin):
    K = RD.K
    L = RD.L
    rho = RD.R.absL
    thL = float(L.imag.mid())                     # arg L = Im L (L = e^L)
    I = acb(0, 1)
    if seite == "rechts":
        th, richtung, w = 0.0, 1, acb("0.31")
    else:
        th, richtung, w = math.pi, -1, acb("-1.75")
    integral = arb(0)
    sup5 = arb(0)
    zellen = []
    t0 = time.time()
    while True:
        d = abs(thL - th)
        if d <= dmin * 1.0000001:
            break
        h = min(0.004, d / 16)                        # Zellbreite
        while True:
            th_b = th + richtung * h
            if abs(thL - th_b) < dmin:
                th_b = thL - richtung * dmin
            hb = abs(th_b - th) / 2
            thm = (th + th_b) / 2
            xm = rho * (I * arb(thm)).exp()
            xrad = rho * arb(hb) * arb("1.0001") + arb(10) ** (-(PREC // 4))
            zc, zr = K.slog_krawczyk(xm, w, xrad)
            if zc is None:
                h /= 2
                continue
            imlo = float((zc.imag - zr).lower())
            imhi = float((zc.imag + zr).upper())
            # Mittelpunkt: eigene (winzige) Krawczyk-Scheibe für slog_K(x_m)
            zm, zmr = K.slog_krawczyk(xm, zc, arb(0))
            if zm is None:
                raise RuntimeError(f"Punkt-Krawczyk versagt bei θ = {thm}")
            if abs(imhi) <= 0.42 and abs(imlo) <= 0.42:
                art = "T"
                ball = RD.t_gebiet(zc, zr, 7)
                mitte = RD.t_gebiet(zm, zmr, 7)
            elif imlo >= 0.39 and float((zm.imag - zmr).lower()) >= 0.39:
                art = "S"
                ball, _, _ = RD.struktur(disc(zc, zr), 7)
                mitte, _, _ = RD.struktur(disc(zm, zmr), 7)
            else:
                h /= 2
                continue
            c6 = ob(mitte[6])
            s7 = ob(ball[7])
            # Mittelwertform: |∂⁶f(θ)| ≤ |∂⁶f(θ_m)| + |θ − θ_m|·sup|∂⁷f|  ⇒  ∫_Zelle ≤ 2h_b|∂⁶f(θ_m)| + h_b²·sup|∂⁷f|
            if arb(hb) * s7 > c6 * arb("0.2") + 1 and h > 1e-9 * max(d, 1e-300):
                h /= 2
                continue
            sup6 = c6 + arb(hb) * s7
            ref = c6
            break
        integral += c6 * arb(2 * hb) + arb(hb) ** 2 * s7
        sup5 = amax(sup5, ob(mitte[5]) + arb(hb) * ob(ball[6]))
        zellen.append(dict(a=th, b=th_b, art=art, sup6=float(sup6), mitte6=float(ref), sup7=float(s7), z=zm.str(8),
                           f=mitte[0].str(12), d1=mitte[1].str(12)))
        w = zm
        th = th_b
        if len(zellen) % 100 == 0:
            print(f"  {seite}: {len(zellen)} Zellen, θ = {th:.12g}, d = {abs(thL - th):.3e}, Art {art}, "
                  f"|∂⁶f(θ_m)| = {float(c6):.4g}, sup|∂⁷f| = {float(s7):.4g}, ∫ bisher = {integral.str(6)}  ({time.time() - t0:.0f}s)", flush=True)
    return integral, sup5, zellen, w


def rest_schranke(RD, w_letzt, dmin):
    """Rest |θ − θ_L| ≤ δ_min einer Seite (siehe Kopf)."""
    R = RD.R
    L = RD.L
    rho = R.absL
    gam = RD.gam
    # |w| ≤ w_t auf dem Rest: |x − L| ≤ ρ δ_min, |w| ≤ |x−L|/(1 − κ₂(w_t)), κ₂ ≤ 2|p₂|w_t (grob, w_t winzig)
    w_t = rho * arb(dmin) * arb("1.01")
    kap2 = sum((ob(R.p[n]) * w_t ** (n - 1) for n in range(2, 40)), arb(0))
    assert kap2 < arb("0.001")
    w_t = w_t / (1 - kap2)
    # Im u im Rest: u = L(z + θ_K(z)); aus der letzten Zelle (z-Kugel) + Spielraum 0.01 (Variation ≤ δ_min + 2κ₂)
    zc = w_letzt
    q0 = (2 * acb.pi() * acb(0, 1) * zc).exp()
    th0 = RD.theta_jets(q0, 0)[0]
    u = L * (zc + th0)
    imu_lo = u.imag.lower() - arb("0.01")
    imu_hi = u.imag.upper() + arb("0.01")
    # sup Im θ_K auf kleinem |q|: θ_K(q) ∈ θ̃(disc(0, q_klein)) ± δθ
    q_klein = arb("1e-30")
    thq = RD.theta_jets(disc(acb(0), q_klein), 0)[0]
    Cq = ob((-gam.imag * (imu_lo if gam.imag > 0 else imu_hi) + 2 * arb.pi() * thq.imag.upper()).exp())
    q_t = Cq * w_t ** gam.real
    assert q_t < q_klein
    # G_k auf disc(0, q_t): M_k; Schwarz: |G_k(q)| ≤ M_k |q|/q_t (k ≥ 2, G_k(0) = 0)
    Gw, R0 = RD.G_werte(disc(acb(0), q_t))
    Mk = [ob(g) for g in Gw]
    fak = [1]
    for k in range(1, 7):
        fak.append(fak[-1] * k)
    sL = ob(RD.sqrtL)
    Hk = [disc(acb(0), sL * w_t * ob(R0)), RD.sqrtL * Gw[1]]
    for k in range(2, 6):
        Hk.append(disc(acb(0), sL * (Mk[k] / q_t) * Cq * w_t ** (gam.real + 1 - k) / fak[k]))
    Hk.append(acb(0))                                         # H_6 separat
    w0 = disc(acb(0), w_t)
    aus, Pw, PH = RD._aus_H(w0, Hk, M)
    B0 = ob(aus[6])
    A0 = acb(0, 1) * Pw[0] / Pw[1]
    c6 = ob(A0) ** 6 * 720 * ob(PH[1])
    K6 = c6 * sL * (Mk[6] / q_t) * Cq / 720                  # |c6 H_6| ≤ K6 |w|^{Reγ − 5}
    # |w| ≥ c_w |τ|: |x − L| = 2ρ|sin(τ/2)| ≥ 0.99ρ|τ|, |w| ≥ |x−L|/(1+κ₂)
    c_w = arb("0.99") * rho / (1 + kap2)
    ex = gam.real - 5                                         # −0.553
    dm = arb(dmin)
    integral = B0 * dm + K6 * c_w ** ex * dm ** (ex + 1) / (ex + 1)
    return integral, dict(w_t=w_t.str(3), q_t=q_t.str(3), Cq=Cq.str(3), B0=B0.str(3), K6=K6.str(3),
                          M6=Mk[6].str(3), sup5_rest=ob(aus[5]).str(3), f_rest=aus[0].str(6), d1_rest=aus[1].str(6))


def test(RD):
    """Vergleich T-Gebiet ↔ Strukturform im Überlapp (Punktwerte) und Plausibilität."""
    K = RD.K
    rho = RD.R.absL
    w = acb("0.31")
    for th in [0.05, 0.2, 0.36, 0.39, 0.40]:
        x = rho * (acb(0, 1) * arb(th)).exp()
        w = K.slog_num(x, w)
        print(f"θ = {th}: z = {w.str(8)}")
        if abs(float(w.imag.mid())) <= 0.42:
            a = RD.t_gebiet(w, arb(0))
            print("   T:", [x.str(10) for x in a])
        if float(w.imag.mid()) >= 0.39:
            b, _, _ = RD.struktur(w)
            print("   S:", [x.str(10) for x in b])
    for th in [1.0, 1.3, 1.337, 1.3372]:
        x = rho * (acb(0, 1) * arb(th)).exp()
        w = K.slog_num(x, w)
        b, _, _ = RD.struktur(w)
        print(f"θ = {th}: z = {w.str(6)}  S: ∂⁶f = {b[6].str(8)}, f = {b[0].str(10)}, ∂f = {b[1].str(10)}")


def main():
    t0 = time.time()
    RD = Rand()
    if SEG == "test":
        test(RD)
        return
    dmin = 10.0 ** (-DMIN_EXP)
    integral, sup5, zellen, w = segment(RD, SEG, dmin)
    rest, info = rest_schranke(RD, w, dmin)
    print(f"{SEG}: {len(zellen)} Zellen, ∫ = {integral.str(8)}, Rest ≤ {rest.str(3)}, sup|∂⁵f| ≤ {sup5.str(5)}  ({time.time() - t0:.0f}s)")
    print("Rest:", json.dumps(info, ensure_ascii=False))
    (HIER / f"p32_rand_{SEG}.json").write_text(json.dumps({
        "segment": SEG, "delta_min": dmin, "prec": PREC, "integral_zellen": integral.str(10, radius=False),
        "integral_obere": float(integral.upper()), "rest": float(rest.upper()), "rest_info": info,
        "sup5": float(sup5), "zellen": zellen}, indent=0, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
