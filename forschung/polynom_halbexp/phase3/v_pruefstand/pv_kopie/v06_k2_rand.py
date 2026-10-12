"""Strang V — (K2) über ∫|∂_θ⁶ f_K(ρe^{iθ})|dθ, auf Grundlage von Ps p32 (Kopie p32_rand_ableitungen.py, unverändert importiert).

Änderungen gegenüber p32 (Diagnose v05, BERICHT §2.2):
  (1) nan/inf in |∂⁶f(θ_m)| oder sup_Zelle|∂⁷f| gilt als Fehlschlag ⇒ Zelle halbieren (p32 nahm nan-Zellen an, weil `nan > x` falsch ist).
  (2) Annahme einer Zelle, sobald hb·sup|∂⁷f| ≤ max(ACC_ABS, ACC_REL·|∂⁶f(θ_m)|). Die Schranke 2hb|∂⁶f(θ_m)| + hb²·sup|∂⁷f| ist für JEDE
      endliche Zelle streng (Mittelwertsatz); die Annahmeregel steuert nur den Aufwand. Ziel C₆ ≤ 5e5 (v04) erlaubt grobe Zellen.
  (3) Zellbreite adaptiv: Start mit min(1.5·h_vorher, H_MAX, d/16).
  (4) Zweigtreue: am gemeinsamen Endpunkt x* zweier Nachbarzellen Punkt-Krawczyk ⇒ winzige Scheibe D* mit genau einer Lösung von
      sexp_K(s) = x*; liegt D* in BEIDEN Zellscheiben, benutzen beide Zellen denselben Zweig von slog (Lücke L3 aus BERICHT §1.1).
  (5) Variable ist der Abstand d = |θ − θ_L| (Python-float, beliebig klein); θ = θ_L ∓ d wird in arb gebildet (p32 rechnete θ als float —
      unbrauchbar für d < 1e-15). Abschnitte [D_A, D_B] einzeln rechenbar; Startwert numerisch fortgesetzt; Abschnittsgrenzen prüft v07.
Ausgabe: k2_<seite>_<i>.json (Integral als obere Schranke, Rand-z-Scheiben, ggf. Rest |θ − θ_L| < δ_min nach p32.rest_schranke).
Aufruf: python v06_k2_rand.py SEITE D_A D_B INDEX [ACC_ABS] [ACC_REL] [H_MAX]
  SEITE rechts (θ < θ_L) / links (θ > θ_L); D_A > D_B Abstände zu θ_L; D_A = "Rand" heisst θ = 0 bzw. π; D_B = "L" heisst bis δ_min (+ Rest).
"""
import json
import math
import sys
import time
from pathlib import Path

from flint import acb, arb

ARGS = [x for x in sys.argv[1:] if x != ""]          # xargs übergibt fehlende Felder als ""
sys.argv = [sys.argv[0], "test", "40"]          # p32 liest argv[1:3] (Segment, δ_min-Exponent); Präzision 700 wie bei P
import p32_rand_ableitungen as P32

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
SEITE = ARGS[0]
D_A_ARG, D_B_ARG = ARGS[1], ARGS[2]
INDEX = int(ARGS[3])
ACC_ABS = float(ARGS[4]) if len(ARGS) > 4 else 1000.0
ACC_REL = float(ARGS[5]) if len(ARGS) > 5 else 10.0
H_MAX = float(ARGS[6]) if len(ARGS) > 6 else 0.004
DMIN = 10.0 ** (-int(ARGS[7])) if len(ARGS) > 7 else 1e-40     # δ_min (p32.rest_schranke verlangt q_t < 1e-30; 1e-12 genügt)
ob = P32.ob


def endlich(x):
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def main():
    t0 = time.time()
    RD = P32.Rand()
    K = RD.K
    L = RD.L
    rho = RD.R.absL
    thL = L.imag                                   # arg L = Im L exakt (L = e^L), als arb
    thLf = float(thL.mid())
    I = acb(0, 1)
    s = -1 if SEITE == "rechts" else 1             # θ = θ_L + s·d
    d_rand = thLf if SEITE == "rechts" else math.pi - thLf
    D_A = d_rand if D_A_ARG == "Rand" else float(D_A_ARG)
    bis_L = D_B_ARG == "L"
    D_B = DMIN if bis_L else float(D_B_ARG)

    def x_von(d):
        return rho * (I * (thL + s * arb(d))).exp()

    # numerischer Startwert: Fortsetzung von slog ab θ = 0 bzw. π, gleichmässig bis d = max(D_A, 0.05), dann geometrisch
    w = acb("0.31") if SEITE == "rechts" else acb("-1.75")
    d_gl = max(D_A, 0.05)
    n = max(1, int((d_rand - d_gl) / 0.002))
    for i in range(1, n + 1):
        w = K.slog_num(x_von(d_rand - (d_rand - d_gl) * i / n), w)
    dd = d_gl
    while dd > D_A * 1.0000001:
        dd = max(dd * 0.8, D_A)
        w = K.slog_num(x_von(dd), w)
    print(f"[{SEITE} {INDEX}] Start d = {D_A:.6g}, z ≈ {w.str(8)}  ({time.time() - t0:.0f}s)", flush=True)
    integral = arb(0)
    sup5 = arb(0)
    zellen = []
    d = D_A
    h_vor = H_MAX
    z_vorher = None
    zweig_ok = True
    erste = None
    while d > D_B * 1.0000001:
        h = min(1.5 * h_vor, H_MAX, d / 16)
        while True:
            d_b = max(d - h, D_B)
            hb = (d - d_b) / 2
            dm = (d + d_b) / 2
            xm = x_von(dm)
            xrad = rho * arb(hb) * arb("1.0001") + arb(10) ** (-(P32.PREC // 4))
            zc, zr = K.slog_krawczyk(xm, w, xrad)
            if zc is None:
                h /= 2
                continue
            zm, zmr = K.slog_krawczyk(xm, zc, arb(0))
            if zm is None:
                h /= 2
                continue
            if not bool(abs(zm - zc) + zmr < zr):        # (F1, K 21:57): Mittelpunkt-Scheibe muss in der Zellscheibe liegen ⇒ derselbe Zweig
                h /= 2
                continue
            imlo = float((zc.imag - zr).lower())
            imhi = float((zc.imag + zr).upper())
            try:
                if abs(imhi) <= 0.42 and abs(imlo) <= 0.42:
                    art = "T"
                    ball = RD.t_gebiet(zc, zr, 7)
                    mitte = RD.t_gebiet(zm, zmr, 7)
                elif imlo >= 0.39 and float((zm.imag - zmr).lower()) >= 0.39:
                    art = "S"
                    ball, _, _ = RD.struktur(P32.disc(zc, zr), 7)
                    mitte, _, _ = RD.struktur(P32.disc(zm, zmr), 7)
                else:
                    h /= 2
                    continue
            except (AssertionError, ValueError, RuntimeError, ZeroDivisionError):
                h /= 2
                if h < 1e-12 * d:
                    raise
                continue
            c6, s7 = ob(mitte[6]), ob(ball[7])
            c5, s6 = ob(mitte[5]), ob(ball[6])
            if not all(endlich(v) for v in (c6, s7, c5, s6)):
                h /= 2
                if h < 1e-14 * d:
                    raise RuntimeError(f"nicht endlich trotz kleiner Zelle bei d = {dm}")
                continue
            if float(arb(hb) * s7) > max(ACC_ABS, ACC_REL * float(c6)) and h > 1e-12 * d:
                h /= 2
                continue
            break
        if z_vorher is not None:
            zc0, zr0 = z_vorher
            wk, rk = K.slog_krawczyk(x_von(d), zm, arb(0))
            if not (wk is not None and bool(abs(wk - zc0) + rk < zr0) and bool(abs(wk - zc) + rk < zr)):
                zweig_ok = False
                print(f"  WARNUNG Zweigtreue offen bei d = {d:.6g}", flush=True)
        else:
            erste = (zc.str(40), zr.str(10))
        integral += c6 * arb(2 * hb) + arb(hb) ** 2 * s7
        sup5 = P32.amax(sup5, c5 + arb(hb) * s6)
        zellen.append([d, d_b, art, float(c6), float(s7)])
        z_vorher = (zc, zr)
        w = zm
        d = d_b
        h_vor = 2 * hb
        if len(zellen) % 200 == 0:
            print(f"  [{SEITE} {INDEX}] {len(zellen)} Zellen, d = {d:.4e}, Art {art}, h = {2*hb:.2e}, |∂⁶f| = {float(c6):.4g}, "
                  f"sup|∂⁷f| = {float(s7):.4g}, ∫ = {integral.str(6)}  ({time.time() - t0:.0f}s)", flush=True)
    aus = {"seite": SEITE, "index": INDEX, "d_a": D_A, "d_b": d, "bis_L": bis_L, "acc": [ACC_ABS, ACC_REL, H_MAX],
           "integral_obere": float(integral.upper()), "sup5": float(sup5.upper()), "n_zellen": len(zellen), "zweig_ok": zweig_ok,
           "erste_z": erste, "letzte_z": (z_vorher[0].str(40), z_vorher[1].str(10)) if z_vorher else None,
           "zeit": round(time.time() - t0, 1)}
    if bis_L:
        try:
            rest, info = P32.rest_schranke(RD, w, DMIN)
            aus["rest_obere"] = float(rest.upper())
            aus["rest_info"] = info
        except Exception as e:                      # Ergebnis der Zellen trotzdem sichern
            aus["rest_fehler"] = f"{type(e).__name__}: {e}"
            aus["rest_w"] = w.str(40)
    (HIER / f"k2_{SEITE}_{INDEX}.json").write_text(json.dumps(aus, indent=0, ensure_ascii=False), encoding="utf-8")
    (HIER / f"k2_{SEITE}_{INDEX}_zellen.json").write_text(json.dumps(zellen), encoding="utf-8")
    print(f"[{SEITE} {INDEX}] FERTIG: {len(zellen)} Zellen, ∫ ≤ {float(integral.upper()):.6g}, Rest {aus.get('rest_obere')}, "
          f"Zweigtreue {zweig_ok}  ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
