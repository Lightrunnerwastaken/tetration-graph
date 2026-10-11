"""Strang P, Satz P-V Schritt 2: hochgenaue Kneser-Daten (T̃, θ̃) per Cousin-Newton auf dem Zylinder (NUMERIK, nicht rigoros).

Ein Schritt (um die aktuellen Daten T̃, θ̃; Ũ = P(L^{z+θ̃}) ist exakt):
  e(t): Ũ(t + e) = T̃(t) auf der Schleife Γ₊ = [−½, ½] + iH;  d(t): T̃(t + 1 + d) = e^{T̃(t)} auf der Naht Σ = −½ + i[−H, H];
  Knoten-Konstante c = e(n₊) − e(n₊+1) − d(n₊) (Rampe c·(t − n₊) auf Γ₊ ⇒ Knotenbedingung exakt; am Fixpunkt c = 0);
  Ω(z) = C + Σ_γ (1/2i)∫_γ J(t) cot π(t − z) dt  (Γ₊ links→rechts, Γ₋ = Spiegel rechts→links, Σ aufwärts), W(0) = 0;
  V = Ω oben (Fourier v_k = Σ_γ ∫ J e^{−2πikt} dt), W = Ω im geschnittenen Rechteck;
  θ̃_neu = θ̃∘(id + V) + V (Fourier-Faltung bis Ordnung 4), T̃_neu = DFT von F_neu auf |z| = 0.95
  (F_neu = Ũ(z₀ + V(z₀)) für |Im z₀| ≥ 0.35, sonst T̃(z₀ + W(z₀)); z₀ = z − round(Re z), dann exp bzw. log).
Aufruf: python p26b_newton.py EINGABE.json AUSGABE.json [PREC]
"""
import json
import math
import sys
import time
from pathlib import Path

from flint import acb, arb, ctx

import kc_basis as K

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
EIN = HIER / sys.argv[1]
AUS = HIER / sys.argv[2]
PREC = int(sys.argv[3]) if len(sys.argv) > 3 else 870
NL, NS = 2000, 600
KTH = 240
NT = 520
KC = 1024
RC = "0.95"
H_NUM = "0.25"


def m(x):
    return acb(x.real.mid(), x.imag.mid())


def lade():
    d = json.loads(EIN.read_text(encoding="utf-8"))
    t = [acb(arb(a), arb(b)) for a, b in d["t"]]
    c = [acb(arb(a), arb(b)) for a, b in d["c"]]
    return t, c


def gl(n, a, b):
    xs, ws = [], []
    for k in range(n):
        x, w = arb.legendre_p_root(n, k, weight=True)
        xs.append((a + b) / 2 + (b - a) / 2 * x)
        ws.append((b - a) / 2 * w)
    return xs, ws


def solve_e(t, c, tz):
    T, _ = K.T_dT(t, tz)
    e = acb(0)
    for _ in range(8):
        U, dU = K.U_dU(c, tz + e)
        de = m((U - T) / dU)
        e = e - de
        if abs(de).mid() < arb(2) ** (-PREC + 30):
            break
    return m(e)


def solve_d(t, s):
    E = K.T_dT(t, s)[0].exp()
    d = acb(0)
    for _ in range(8):
        T1, dT1 = K.T_dT(t, s + 1 + d)
        dd = m((T1 - E) / dT1)
        d = d - dd
        if abs(dd).mid() < arb(2) ** (-PREC + 30):
            break
    return m(d)


def faltung(a, b, n):
    out = [acb(0)] * n
    for i, ai in enumerate(a[:n]):
        if ai == 0:
            continue
        for j in range(0, n - i):
            if j < len(b):
                out[i + j] += ai * b[j]
    return out


def main():
    ctx.prec = PREC
    K.setup(PREC)
    t, c = lade()
    t = [m(x) for x in t] + [acb(0)] * max(0, NT - len(t))
    t = t[:NT]
    c = [m(x) for x in c]
    t0 = time.time()
    H = arb(H_NUM)
    I = acb(0, 1)
    pi = arb.pi()
    npl = acb(-0.5, H)
    # Schleife
    xs, ws = gl(NL, arb(-0.5), arb(0.5))
    ts = [acb(x, H) for x in xs]
    es = [solve_e(t, c, tz) for tz in ts]
    e_n0 = solve_e(t, c, npl)
    e_n1 = solve_e(t, c, npl + 1)
    print(f"e auf Γ₊: max|e| = {max(float(abs(x).mid()) for x in es):.3e}  ({time.time() - t0:.0f}s)")
    sys.stdout.flush()
    # Naht
    ys, wsy = gl(NS, -H, H)
    ss = [acb(-0.5, y) for y in ys]
    ds = [solve_d(t, s) for s in ss]
    d_n = solve_d(t, npl)
    cc = m(e_n0 - e_n1 - d_n)
    print(f"d auf Σ: max|d| = {max(float(abs(x).mid()) for x in ds):.3e};  Knotenkonstante c = {float(abs(cc).mid()):.3e}  ({time.time() - t0:.0f}s)")
    sys.stdout.flush()
    # Knoten (Stützstelle, Gewicht·Dichte) der drei Bögen
    knoten = []
    for tz, w, e in zip(ts, ws, es):
        J = e + cc * (tz - npl)
        knoten.append((tz, m(w * J)))
        knoten.append((tz.conjugate(), m(-w * J.conjugate())))
    for s, w, d in zip(ss, wsy, ds):
        knoten.append((s, m(I * w * d)))
    halbi = 1 / (2 * I)

    def Omega0(z):
        acc = acb(0)
        for nd, q in knoten:
            acc += q * (pi * (nd - z)).cot()
        return acc * halbi

    w0 = acb(0)                                   # Normierung: T̃(W(0)) = 1, also F(0) = 1
    for _ in range(6):
        T0, dT0 = K.T_dT(t, w0)
        w0 = m(w0 - (T0 - 1) / dT0)
    C = m(w0 - Omega0(acb(0)))
    # Fourier-Koeffizienten von V
    v = [acb(0)] * (KTH + 1)
    v[0] = C + sum((q for _, q in knoten), acb(0)) / 2
    Es = [m((-2 * pi * I * nd).exp()) for nd, _ in knoten]
    pot = [acb(1)] * len(knoten)
    for k in range(1, KTH + 1):
        acc = acb(0)
        for j, (nd, q) in enumerate(knoten):
            pot[j] = m(pot[j] * Es[j])
            acc += q * pot[j]
        v[k] = m(acc)
    print(f"V: |v_0| = {float(abs(v[0]).mid()):.3e}, |v_1| = {float(abs(v[1]).mid()):.3e}, |v_KTH|·e^(-2πK·H) = "
          f"{float((abs(v[KTH]) * (-2 * pi * KTH * H).exp()).mid()):.3e}  ({time.time() - t0:.0f}s)")
    sys.stdout.flush()

    # Schleife mit Ausbuchtung nach oben (für Auswertepunkte knapp unter Γ₊, |Re| < 0.1)
    bump = []
    for (a, b, im_a, im_b, n) in [(-0.5, -0.3, 0, 0, 400), (-0.3, -0.3, 0, 0.15, 300), (-0.3, 0.3, 0.15, 0.15, 1200),
                                  (0.3, 0.3, 0.15, 0, 300), (0.3, 0.5, 0, 0, 400)]:
        za = acb(arb(a), H + arb(im_a))
        zb = acb(arb(b), H + arb(im_b))
        xs_, ws_ = gl(n, arb(-1), arb(1))
        for x_, w_ in zip(xs_, ws_):
            tz = (za + zb) / 2 + (zb - za) / 2 * x_
            dt = (zb - za) / 2 * w_
            e = solve_e(t, c, tz)
            bump.append((tz, m(dt * (e + cc * (tz - npl)))))
    knoten_b = bump + [kq for kq in knoten if float(kq[0].imag.mid()) < float(H.mid()) - 1e-30]
    # (obere Schleife ersetzt; untere Schleife und Naht bleiben)

    def Omega_b(z):
        acc = acb(0)
        for nd, q in knoten_b:
            acc += q * (pi * (nd - z)).cot()
        return acc * halbi

    print(f"Ausbuchtung: {len(bump)} Knoten ({time.time() - t0:.0f}s)")
    sys.stdout.flush()

    def V_eval(z):
        q = (2 * pi * I * z).exp()
        s = acb(0)
        for k in range(KTH, -1, -1):
            s = s * q + v[k]
        return s

    # θ̃_neu = θ̃∘(id+V) + V
    cn = c + [acb(0)] * (KTH + 1 - len(c))
    cn = cn[:KTH + 1]
    neu = [cn[k] + v[k] for k in range(KTH + 1)]
    Vp = list(v)
    fak = 1
    for ordnung in range(1, 5):
        fak *= ordnung
        abl = [cn[k] * (2 * pi * I * k) ** ordnung for k in range(KTH + 1)]
        term = faltung(abl, Vp, KTH + 1)
        neu = [neu[k] + term[k] / fak for k in range(KTH + 1)]
        Vp = faltung(Vp, v, KTH + 1)
    c_neu = [m(x) for x in neu]
    # T̃_neu per DFT auf |z| = RC
    rc = arb(RC)
    werte = [None] * KC
    hn = float(H.mid())
    for j in range(KC // 2 + 1):
        z = rc * (2 * pi * I * j / KC).exp()
        sh = int(round(float(z.real.mid())))
        z0 = z - sh
        im0 = float(z0.imag.mid())
        if im0 >= hn:
            F = K.U_dU(c, z0 + V_eval(z0))[0]
        elif im0 >= hn - 0.1:
            W = Omega_b(z0) + C
            F = K.T_dT(t, z0 + W)[0]
        else:
            W = Omega0(z0) + C
            F = K.T_dT(t, z0 + W)[0]
        if sh == 1:
            F = F.exp()
        elif sh == -1:
            F = F.log()
        elif sh != 0:
            raise RuntimeError(sh)
        if float(F.real.rad()) > 1e-200 or float(F.imag.rad()) > 1e-200:
            print(f"  WARNUNG F[{j}] z0={complex(z0.real.mid(), z0.imag.mid())} rad {float(F.real.rad()):.2e}")
        werte[j] = m(F)
    for j in range(KC // 2 + 1, KC):
        werte[j] = werte[KC - j].conjugate()
    print(f"F_neu auf |z| = {RC} berechnet ({time.time() - t0:.0f}s)")
    sys.stdout.flush()
    t_neu = []
    wurzeln = [m((-2 * pi * I * j / KC).exp()) for j in range(KC)]   # direkt, ohne Radius-Akkumulation
    for n in range(NT):
        acc = acb(0)
        for j in range(KC):
            acc += werte[j] * wurzeln[(n * j) % KC]
        tn = acc / KC / rc ** n
        if not (tn.real.rad() < abs(tn.real.mid()) * arb(10) ** -20 or tn.real.rad() < arb(10) ** -200):
            print(f"  WARNUNG t_{n}: mid {float(tn.real.mid()):.3e} rad {float(tn.real.rad()):.3e}")
        t_neu.append(m(tn))
    im_max = max(float(abs(x.imag).mid()) for x in t_neu)
    aend_t = max(float((abs(t_neu[n] - t[n]) * rc ** n).mid()) for n in range(NT))
    aend_c = max(float((abs(c_neu[k] - (cn[k])) * (-2 * pi * k * H).exp()).mid()) for k in range(KTH + 1))
    print(f"T̃_neu: t₀ − 1 = {float(abs(t_neu[0] - 1).mid()):.3e}, max|Im t_n| = {im_max:.3e}, Änderung (·0.95^n) {aend_t:.3e}; "
          f"θ̃-Änderung (·e^(-2πk·H)) {aend_c:.3e}  ({time.time() - t0:.0f}s)")
    t_neu = [acb(x.real.mid()) for x in t_neu]
    dig = int(PREC * 0.30103) + 5
    raus = {"prec": PREC, "H": H_NUM, "diagnose": {"max_e": max(float(abs(x).mid()) for x in es),
                                                     "max_d": max(float(abs(x).mid()) for x in ds),
                                                     "c": float(abs(cc).mid()), "aend_t": aend_t, "aend_c": aend_c},
            "t": [[x.real.mid().str(dig, radius=False), "0"] for x in t_neu],
            "c": [[x.real.mid().str(dig, radius=False), x.imag.mid().str(dig, radius=False)] for x in c_neu]}
    AUS.write_text(json.dumps(raus, indent=0), encoding="utf-8")
    print("gespeichert", AUS.name)


if __name__ == "__main__":
    main()
