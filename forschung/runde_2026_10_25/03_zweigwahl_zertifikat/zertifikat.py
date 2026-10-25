"""Gebiet 03: Intervall-Zertifikat der Zweigwahl (S4) für Basis e und 2, mit python-flint (arb/acb).

Behauptung S4 (Satz Z, luecken/BERICHT_LUECKEN_ANALYTISCH.md):  sigma(sexp z) = exp(l*(z+theta(z))) auf D,
l = Log(lambda), lambda = L ln b. Nach Lemma H genügt EIN Punkt z0. Hier wird gezeigt:

    P := exp(l*(z0 + theta_K(z0)))  liegt in sigma(H)                                  (*)

Daraus folgt S4 an z0 OHNE Psi auszuwerten: P = sigma(w') mit w' in H, also Psi(P) = w' (Lemma U(c)),
und Psi(P) = F(z0+theta_K(z0)) = sexp(z0) (Kneser-Darstellung V1), also sigma(sexp z0) = P.

Werkzeuge (alle rigoros, Bälle von arb):
  L        Intervall-Newton für b^w = w (Einschluss + Eindeutigkeit im Ball).
  sigma    sigma(w) = lambda^n sigma(v_n), v_{k+1} = Log(v_k)/ln b (Hauptzweig; bildet H in H ab),
           und nahe L de Branges: f(t) = sigma(L+R t)/R ist schlicht auf der Einheitsscheibe (B(L,R) in H, R <= Im L),
           f(0) = 0, f'(0) = 1, |a_k| <= k  =>  |sigma(L+u) - u| <= R*(t/(1-t)^2 - t),
                                               |sigma'(L+u) - 1| <= (1+t)/(1-t)^3 - 1,   t = |u|/R.
  Koebe-1/4  sigma schlicht auf B(w1, Im w1)  =>  sigma(B(w1, Im w1)) enthält B(sigma(w1), |sigma'(w1)| Im w1 / 4).
  Testpunkt w1 ist beliebig (keine Annahme): hier per Gleitkomma-Newton so gewählt, dass sigma(w1) ~ P.
Einzige Motor-Annahme (E): |theta_K(z0) - theta_num(z0)| <= tau  (mod 2 pi i / l), theta_num = Summe der Datei-a_n.
Ausgegeben wird tau_max: das grösste tau, für das (*) noch folgt.

Aufruf: C:/Users/janis/miniconda3/python.exe zertifikat.py      (Laufzeit < 1 min)
"""
import json
import sys
from pathlib import Path

import mpmath as mp
from flint import acb, arb, ctx

HIER = Path(__file__).resolve().parent
ID = HIER.parent.parent / "identitaeten_2026_10_14"
E2 = ID / "etappe2" / "daten"
E3 = ID / "etappe3" / "daten"

ctx.prec = 800
mp.mp.dps = 200
N_IT = {"e": 450, "2": 700}   # Log-Schritte: de-Branges-Rest ~ |lambda|^-n ~ 1e-61

DATEN = {
    "e": {"haupt": E2 / "a_e_80_0.003_16384.json", "zweit": E2 / "a_e_60_0.002_16384.json",
          "taylor": E3 / "taylor_130_256.json"},
    "2": {"haupt": E2 / "a_2_80_0.003_16384.json", "zweit": E2 / "a_2_60_0.002_16384.json",
          "taylor": E3 / "taylor_b2_80_256.json"},
}
PUNKTE = [("0.05i", mp.mpc(0, "0.05")), ("0.07e^{i pi/3}", mp.mpf("0.07") * mp.expjpi(mp.mpf(1) / 3)),
          ("0.08e^{3i pi/4}", mp.mpf("0.08") * mp.expjpi(mp.mpf(3) / 4)), ("0.5+0.05i", mp.mpc("0.5", "0.05"))]


def A(x):
    """mpmath -> exakter acb-Punkt (Mittelpunkt, Radius 0)."""
    x = mp.mpc(x)
    return acb(arb(mp.nstr(x.real, 220)), arb(mp.nstr(x.imag, 220)))


def M(x):
    """acb -> mpmath (Mittelpunkt)."""
    m = x.mid()
    return mp.mpc(mp.mpf(m.real.str(230, radius=False)), mp.mpf(m.imag.str(230, radius=False)))


def up(x):
    """obere Schranke einer arb-Zahl als mpf."""
    return mp.mpf((x.mid() + x.rad()).str(30, radius=False)) * (1 + mp.mpf(10) ** -25)


def lo(x):
    return mp.mpf((x.mid() - x.rad()).str(30, radius=False)) * (1 - mp.mpf(10) ** -25)


class Basis:
    def __init__(self, b):
        self.name = b
        self.lb = arb(1) if b == "e" else arb(2).log()
        lbm = mp.mpf(1) if b == "e" else mp.log(2)
        Lc = -mp.lambertw(-lbm) / lbm
        Lm = mp.conj(Lc) if mp.im(Lc) < 0 else Lc
        # Intervall-Newton für f(w) = exp(lb w) - w auf X = B(Lm, 1e-40)
        r = arb("1e-40")
        X = acb(arb(mp.nstr(Lm.real, 220), r), arb(mp.nstr(Lm.imag, 220), r))
        m = A(Lm)
        f = lambda w: (self.lb * w).exp() - w
        fp = lambda w: self.lb * (self.lb * w).exp() - 1
        N = m - f(m) / fp(X)
        assert X.contains(N), "Intervall-Newton für L fehlgeschlagen"
        self.L = N                                   # enthält den eindeutigen Fixpunkt in X
        self.lam = self.L * self.lb
        self.l = self.lam.log()                      # Log lambda (Hauptzweig), F(s) = Psi(e^{l s})
        assert lo(self.L.imag) > 0 and up(self.L.imag * self.lb) < 3.14159
        self.R = arb(1)                              # B(L, 1) liegt in H (Im L > 1.33 für e, > 1.56 für 2)
        assert lo(self.L.imag) > 1.01
        self.Lm, self.lbm = Lm, lbm
        self.lamm = Lm * lbm

    # ---- Gleitkomma (nur für den Testpunkt w1)
    def sigma_f(self, w, n=None):
        n = n or N_IT[self.name]
        d = mp.mpc(1)
        for _ in range(n):
            d /= (w * self.lbm)
            w = mp.log(w) / self.lbm
        return self.lamm ** n * (w - self.Lm), self.lamm ** n * d

    def w1_aus(self, P, start):
        w = mp.mpc(start)
        for _ in range(60):
            s, ds = self.sigma_f(w)
            dw = (s - P) / ds
            w -= dw
            if abs(dw) < mp.mpf(10) ** -150:
                break
        return w

    # ---- rigoros
    def sigma(self, w, n):
        """Einschluss von sigma(w) und sigma'(w) (w: acb in H)."""
        prod = acb(1)
        v = w
        for _ in range(n):
            assert lo(v.imag) > 0, "Iterierte verlässt H"
            prod = prod / (v * self.lb)
            v = v.log() / self.lb
        u = v - self.L
        t = u.abs_upper() / self.R
        assert up(t) < 0.5
        E = self.R * (t / (1 - t) ** 2 - t)
        Ed = (1 + t) / (1 - t) ** 3 - 1
        su = u + kugel(E)
        dsu = 1 + kugel(Ed)
        lamn = self.lam ** n
        return lamn * su, lamn * dsu * prod, t


def kugel(r):
    """komplexer Ball um 0, der die Kreisscheibe |x| <= r enthält (r: arb)."""
    rs = mp.nstr(up(r), 25)
    return acb(arb(0, rs), arb(0, rs))


def lade_a(datei):
    d = json.loads(Path(datei).read_text(encoding="utf-8"))
    return [acb(arb(x[0]), arb(x[1])) for x in d["a"]]


def theta(a, z):
    q = (2 * acb.pi() * acb(0, 1) * z).exp()
    s = acb(0)
    qn = acb(1)
    for ak in a:
        s += ak * qn
        qn *= q
    return s


def lade_taylor(datei):
    d = json.loads(Path(datei).read_text())
    z = d["zentren"] if "zentren" in d else d
    return {float(k): [acb(arb(x)) for x in c] for k, c in z.items()}


def poly(c, z):
    s = acb(0)
    for ck in reversed(c):
        s = s * z + ck
    return s


def zertifikat(B, a, a2, tay, name, z0m):
    z0 = A(z0m)
    th = theta(a, z0)
    th2 = theta(a2, z0)
    X = B.l * (z0 + th)
    P = X.exp()
    Pm = M(P)
    # Testpunkt: Start aus Taylor-Daten (nur als Startwert, keine Annahme)
    zc = 0.0 if abs(z0m) < 0.1 else 0.4
    tz = poly(tay[zc], z0 - A(zc))
    start = M(tz)
    w1m = B.w1_aus(Pm, start)
    w1 = A(w1m)
    n = N_IT[B.name]
    sw, dsw, t = B.sigma(w1, n)
    r1 = w1.imag
    rho = dsw.abs_lower() * r1 / 4
    dist = (P - sw).abs_upper()
    # tau_max: |P - sigma(w1)| + |P| (e^{|l| tau} - 1) < rho
    marge = rho - dist
    assert lo(marge) > 0
    tau = ((marge / P.abs_upper()).log1p() / B.l.abs_upper())
    # Folgerung: sexp(z0) in B(w1, r1 t*), t*/(1+t*)^2 = eps, eps = dist/(r1 |sigma'|)
    eps = dist / (r1 * dsw.abs_lower())
    tstar = ((1 - 2 * eps) - (1 - 4 * eps).sqrt()) / (2 * eps)
    rs = r1 * tstar
    st = (w1 - tz).abs_upper()
    # Zweigzahl k (Konvention lib3/kern.h): theta = 3 - z + h(u3) + 2 pi i k / l
    u = w1m
    for _ in range(3):
        u = mp.log(u) / B.lbm
    nh = 0
    v = u
    while abs(v - B.Lm) > mp.mpf("1e-3") or (B.name == "2" and nh < 5):
        v = mp.log(v) / B.lbm
        nh += 1
    Y = X - (nh + 3) * B.l          # P / lambda^{nh+3} = e^Y, h(u3) = nh + Log(e^Y)/l
    k_lo = (Y.imag / (2 * arb.pi()))
    return {
        "basis": B.name, "z0": name,
        "theta_num(z0)": th.str(20, radius=False), "rad theta_num": float(up(th.rad())),
        "|theta80-theta60|": float(up((th - th2).abs_upper())),
        "|P|": float(up(P.abs_upper())), "rad P": float(up(P.rad())),
        "Im w1": float(lo(r1)), "t_deBranges": float(up(t)),
        "dist_up": float(up(dist)), "rho_lo": float(lo(rho)), "tau_max_lo": float(lo(tau)),
        "sexp(z0) Mitte": mp.nstr(w1m, 25), "sexp(z0) Radius": float(up(rs)),
        "|w1 - sexp_Taylor|": (float(up(st)) if st is not None else None),
        "n_h(lib3)": nh, "Im(Y)/2pi (k=0 iff in (-1/2,1/2))": [float(lo(k_lo)), float(up(k_lo))],
    }


def main():
    alle = []
    for b in ["e", "2"]:
        B = Basis(b)
        print(f"Basis {b}: L = {B.L.str(25)}   l = Log(lambda) = {B.l.str(20)}")
        a = lade_a(DATEN[b]["haupt"])
        a2 = lade_a(DATEN[b]["zweit"])[: len(a)]
        tay = lade_taylor(DATEN[b]["taylor"])
        for name, z0m in PUNKTE:
            r = zertifikat(B, a, a2, tay, name, z0m)
            alle.append(r)
            ok = r["dist_up"] < r["rho_lo"]
            print(f"  z0 = {name:16s} |P - sigma(w1)| <= {r['dist_up']:.2e}  <  rho >= {r['rho_lo']:.4f}  {'OK' if ok else 'FEHLT'}"
                  f"   tau_max >= {r['tau_max_lo']:.4e}   |th80-th60| = {r['|theta80-theta60|']:.1e}"
                  f"   sexp(z0) Radius {r['sexp(z0) Radius']:.1e}  |w1-Taylor| {r['|w1 - sexp_Taylor|']}"
                  f"   n_h={r['n_h(lib3)']} Im Y/2pi in {r['Im(Y)/2pi (k=0 iff in (-1/2,1/2))']}")
    (HIER / "daten").mkdir(exist_ok=True)
    (HIER / "daten" / "zertifikat.json").write_text(json.dumps(alle, indent=1, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
