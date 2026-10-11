"""Strang P, Satz P-V Schritt 5: ZERTIFIZIERTE Taylor-Koeffizienten von f_K = sexp(slog + ½) um 0 (rigoros, arb).

Grundlage: p27/p28 — sexp_K = T̃ ± ε_T auf R_L = [−½,½]×[−H,H], sexp_K = Ũ ± ε_U auf Re ∈ [−½,½], Im ≥ H − ρ/2 (gespiegelt),
Funktionalgleichung (exp/Log, Prüfungen A/B) außerhalb des Streifens.
(1) ε_T, ε_U aus |φ − id|, |τ|, d_max, Y_max und groben Ableitungsschranken.
(2) sexp_ball(S): Ball, der sexp_K auf der (kleinen) Kreisscheibe S enthält (Reduktion in den Streifen, Darstellung, exp/Log).
(3) f_K(x) für x auf |x| = r_f: numerisches w̃ = slog(x), Krawczyk 1D für sexp_K(s) = x auf D(w̃, r₀) ⇒ f_K(x) ∈ sexp_ball(D + ½).
(4) DFT auf |x| = r_f (K Punkte, Konjugationssymmetrie) ⇒ a_n mit Ballradius; Aliasfehler ≤ M r'^{−n}(r_f/r')^K/(1 − (r_f/r')^K),
    M ≥ sup_{|x|=r'}|f_K| über Krawczyk für Bogenfamilien (x ∈ Bogenscheibe) auf |x| = r'.
Aufruf: python p29_fK_koeffizienten.py DATEN.json [prec] [K] [nmax]
"""
import json
import math
import sys
import time
from pathlib import Path

from flint import acb, arb, ctx

import kc_basis as KB
import kc_rig

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
DATEN = HIER / sys.argv[1]
PREC = int(sys.argv[2]) if len(sys.argv) > 2 else 700
KDFT = int(sys.argv[3]) if len(sys.argv) > 3 else 4096
NMAX = int(sys.argv[4]) if len(sys.argv) > 4 else 640
RF = "1.2"
RA = "1.30"


def amax(a, b):
    return (a + b + abs(a - b)) / 2


def main():
    t0 = time.time()
    ctx.prec = PREC
    Z7 = json.loads((HIER / "p27_zert.json").read_text(encoding="utf-8"))
    Z8 = json.loads((HIER / "p28_zert.json").read_text(encoding="utf-8"))
    H, RHO = arb(Z7["H"]), arb(Z7["rho"])
    d_max = arb(Z7["d_max"]).abs_upper()
    e_max = arb(Z7["e_max"]).abs_upper()
    phi_f = arb(Z8["phi_fehler"]).abs_upper()
    assert Z8["pruefung_A"] and Z8["pruefung_B"]
    R = kc_rig.Rig(PREC, n0=300, m_extra=250)
    d = json.loads(DATEN.read_text(encoding="utf-8"))
    t = [acb(arb(arb(a).mid())) for a, b in d["t"]]
    c = [acb(arb(arb(a).mid()), arb(arb(b).mid())) for a, b in d["c"]]
    dt = [n * t[n] for n in range(1, len(t))]
    # ---------------- (1) grobe Ableitungsschranken und ε_T, ε_U
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
            supU1 = amax(supU1, sum((j * abs(a[j]) * arb("0.075") ** (j - 1) for j in range(1, 11)), arb(0)) + tau / arb("0.025"))
            y += 0.1
    # Im ≥ 3: |Ũ'| ≤ |P'(w)|·|w|·|L|·(1+|θ̃'|) mit |w| ≤ w_oben (p28), |P'| ≤ 1 + Σ n|p_n|0.1^{n−1} + Rest auf |w| ≤ 0.1
    w_oben = arb(Z8["w_oben"]).abs_upper()
    assert w_oben < arb("0.1")
    P1 = 1 + sum((n * abs(R.p[n]) * arb("0.1") ** (n - 1) for n in range(2, R.n0 + 1)), arb(0)) + R.C * arb("0.2") ** (R.n0 + 1) / arb("0.1")
    th1 = sum((tpi_ * kk * abs(ck) * (-tpi_ * 3 * kk).exp() for kk, ck in enumerate(c)), arb(0))
    supU1 = amax(supU1, P1 * w_oben * R.absL * (1 + th1))
    lip = amax(supT1, supU1) * arb("1.001")
    tn = abs(t[0] - 1)
    tau_b = 2 * (tn + lip * phi_f)
    epsT = lip * (phi_f + tau_b) + supT1 * d_max
    epsU = lip * (phi_f + tau_b) + supU1 * (d_max + e_max)
    etaB = arb(Z7["eta_B"]).abs_upper()
    eps = amax(epsT, epsU) + etaB          # gilt für beide Darstellungen im ganzen Band |Im| ≤ H + ρ/2 und auf |Re| ≤ 0.6
    print(f"sup|T̃'| ≤ {supT1.str(4)}, sup|Ũ'| ≤ {supU1.str(4)};  |τ| ≤ {tau_b.str(3)};  ε_T ≤ {epsT.str(3)}, ε_U ≤ {epsU.str(3)}  ({time.time() - t0:.0f}s)")
    sys.stdout.flush()
    hf = float(H.mid())

    # ---------------- (2) sexp-Ball auf kleinen Scheiben (und Ableitung)
    def rep_ball(S, m_tay=8):
        """(Wert, Ableitung) als Bälle für die Darstellung im Streifen, S = Scheibe (Mitte, Radius) mit Re ∈ [−½, ½)."""
        zc, r = S
        im = float(zc.imag.mid())
        if abs(im) <= hf:
            Dk = kc_rig.disc(zc, r) if r > 0 else zc
            v = R.T_ball(t, Dk)
            dv = R.T_ball(dt, Dk)
            e = eps
        else:
            spiegel = im < 0
            zz = zc.conjugate() if spiegel else zc
            if r > arb("1e-100"):
                a, tau = R.U_taylor(c, zz, m_tay, r, arb("0.1"))
                v = a[0] + acb(arb(0, 1), arb(0, 1)) * (sum((abs(a[j]) * r ** j for j in range(1, m_tay + 1)), arb(0)) + tau)
                dv = a[1] + acb(arb(0, 1), arb(0, 1)) * (sum((j * abs(a[j]) * r ** (j - 1) for j in range(2, m_tay + 1)), arb(0)) + 2 * tau / r)
            else:
                v, dv = R.U_dU_ball(c, zz if r == 0 else kc_rig.disc(zz, r))
            if spiegel:
                v, dv = v.conjugate(), dv.conjugate()
            e = eps
        # Fehlerbälle: Wert ±ε; Ableitung ±ε/0.01 (Cauchy, ε gilt auf einer 0.01-Umgebung)
        eb = acb(arb(0, e), arb(0, e))
        return v + eb, dv + eb * 100

    def sexp_ball(S):
        zc, r = S
        k = 0
        while float(zc.real.mid()) >= 0.5:
            zc -= 1
            k += 1
        while float(zc.real.mid()) < -0.5:
            zc += 1
            k -= 1
        v, dv = rep_ball((zc, r))
        for _ in range(k):
            v = v.exp()
            dv = dv * v
        for _ in range(-k):
            dv = dv / v
            v = v.log()
        return v, dv

    # numerisches slog (schnell, nur Startwert)
    KB.setup(PREC)
    tn_num = [KB.mid(x) for x in t]
    cn_num = [KB.mid(x) for x in c]

    def sexp_num(w):
        k = 0
        while float(w.real.mid()) >= 0.5:
            w -= 1
            k += 1
        while float(w.real.mid()) < -0.5:
            w += 1
            k -= 1
        im = float(w.imag.mid())
        if abs(im) <= hf:
            v, dv = KB.T_dT(tn_num, w)
        elif im > 0:
            v, dv = KB.U_dU(cn_num, w)
        else:
            v, dv = KB.U_dU(cn_num, w.conjugate())
            v, dv = v.conjugate(), dv.conjugate()
        for _ in range(k):
            v = v.exp()
            dv = dv * v
        for _ in range(-k):
            dv = dv / v
            v = v.log()
        return KB.mid(v), KB.mid(dv)

    def slog_num(x, w):
        for _ in range(80):
            v, dv = sexp_num(w)
            dw = KB.mid((v - x) / dv)
            w = w - dw
            if abs(dw).mid() < arb(2) ** (-PREC + 40):
                break
        return w

    def fK_zert(x, w0, r0=None, xrad=arb(0)):
        """f_K auf der Scheibe D(x, xrad) per Krawczyk auf D(w̃, r0); Rückgabe (Ball von f_K, w̃) oder (None, w̃)."""
        w = slog_num(x, w0)
        g0, dg = sexp_ball((w, arb(0)))
        g0 = g0 - x + acb(arb(0, xrad), arb(0, xrad))
        if r0 is None:
            r0 = 4 * g0.abs_upper() / dg.abs_lower() + arb(10) ** -300
        _, Dp = sexp_ball((w, r0))
        d0 = KB.mid(Dp)
        Kc = -g0 / d0
        Kr = abs(1 - Dp / d0) * r0
        if not (abs(Kc) + Kr < r0):
            return None, w
        v, _ = sexp_ball((w + arb("0.5"), r0))
        return v, w

    # ---------------- (3)+(4) DFT auf |x| = r_f
    rf = arb(RF)
    I = acb(0, 1)
    tpi = 2 * arb.pi()
    werte = [None] * KDFT
    w = acb(0)
    for j in range(KDFT // 2 + 1):
        x = rf * (tpi * I * j / KDFT).exp()
        v, w = fK_zert(x, w)
        if v is None:
            raise RuntimeError(f"Krawczyk versagt bei j = {j}")
        werte[j] = v
        if j % 256 == 0:
            print(f"  f_K auf |x| = {RF}: j = {j}, f_K ≈ {v.str(8)}, Radius {float(v.real.rad() + v.imag.rad()):.2e}  ({time.time() - t0:.0f}s)")
            sys.stdout.flush()
    for j in range(KDFT // 2 + 1, KDFT):
        werte[j] = werte[KDFT - j].conjugate()
    wurzeln = [(-tpi * I * j / KDFT).exp() for j in range(KDFT)]
    koeff = []
    for n in range(NMAX + 1):
        acc = acb(0)
        for j in range(KDFT):
            acc += werte[j] * wurzeln[(n * j) % KDFT]
        koeff.append(acc / KDFT / rf ** n)
    print(f"DFT fertig ({time.time() - t0:.0f}s)")
    # ---------------- M ≥ sup_{|x|=r'}|f_K| über Bogenscheiben
    ra = arb(RA)
    rad0 = arb("0.008")
    hmax = 0.002

    def bogen_h(phi, w):
        x = ra * (I * arb(phi)).exp()
        ww = slog_num(x, w)
        _, dv = sexp_num(ww)
        return min(hmax, 0.25 * float(rad0.mid()) * float(abs(dv).mid())), ww

    M = arb(0)
    w = acb(0)
    phi = 0.0
    h, w = bogen_h(phi, w)
    anzahl = 0
    while True:
        x = ra * (I * arb(phi)).exp()
        v, w = fK_zert(x, w, r0=rad0, xrad=arb(h) * arb("1.0001"))
        if v is None:
            raise RuntimeError(f"Bogen-Krawczyk versagt bei φ = {phi}")
        M = amax(M, v.abs_upper())
        anzahl += 1
        if phi >= math.pi:
            break
        schritt = 1.8 * h / float(ra.mid())
        while True:
            phi_n = min(math.pi, phi + schritt)
            h_n, w_n = bogen_h(phi_n, w)
            # Überdeckung: Bogenabstand ≤ 0.99·(h + h_n)  (Sehne ≤ Bogen; Scheibenradius h deckt den Bogen bis Winkelabstand h/r')
            if (phi_n - phi) * float(ra.mid()) <= 0.99 * (h + h_n):
                break
            schritt *= 0.7
        phi, h = phi_n, h_n
    print(f"Bogenüberdeckung |x| = {RA}: {anzahl} Scheiben (adaptiv)")
    q = (rf / ra) ** KDFT
    print(f"M = sup_|x|={RA} |f_K| ≤ {M.str(4)};  (r_f/r')^K = {q.str(3)}  ({time.time() - t0:.0f}s)")
    rho = R.absL
    aus = []
    for n, a in enumerate(koeff):
        alias = M * ra ** (-n) * q / (1 - q)
        rad = a.real.rad() + alias
        aus.append([a.real.mid().str(140, radius=False), rad.str(5, radius=False),
                    float((rad * rho ** n).mid()), a.imag.str(5)])
    gew = max(x[2] for x in aus)
    print(f"Koeffizienten a_0..a_{NMAX}: max_n Fehler·ρ^n = {gew:.3e}")
    for n in (0, 1, 2, 4, 100, 300, 600):
        if n <= NMAX:
            print(f"  a_{n} = {koeff[n].real.str(30)}  (Fehler {aus[n][1]})")
    (HIER / "p29_fK_koeff.json").write_text(json.dumps({
        "beschreibung": "Zertifizierte Taylor-Koeffizienten von Knesers f_K um 0 (Strang P, Satz P-V). Eintrag: [Mittelpunkt, Fehlerschranke, Fehler·ρ^n, Imaginärteil-Ball]",
        "eps_T": epsT.str(5), "eps_U": epsU.str(5), "M": M.str(5), "r_f": RF, "r_alias": RA, "K": KDFT,
        "a": aus}, indent=0), encoding="utf-8")
    print(f"fertig ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
