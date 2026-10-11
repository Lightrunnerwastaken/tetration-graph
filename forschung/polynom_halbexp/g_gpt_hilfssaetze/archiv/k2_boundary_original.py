"""Interval bound for Kneser's sixth angular derivative.

Reads the existing P-V certificates; writes only in this workspace.
The statement is conditional on the mathematical validity of those certificates.
"""
import sys, json, math, time
from pathlib import Path
sys.dont_write_bytecode = True
sys.stdout.reconfigure(encoding="utf-8")
BASE = Path(r"D:\Users\janis\Documents\Matura\tetration-graph\forschung\polynom_halbexp\p_pruefung")
sys.path.insert(0, str(BASE))
import p32_rand_ableitungen as old
from flint import acb, acb_series, arb, ctx

old.PREC = 700
OUT = Path(__file__).resolve().parents[1] / "outputs"

def finite(v):
    if not v.is_finite():
        raise ArithmeticError("Non-finite interval encountered")
    return v

def direct_P_jet(self, y, m):
    """Taylor coefficients of P at every point in y, by scaling and exp.

    |P-T_n| <= C |w|^(n+1) on |w| <= r0. Cauchy's estimate on
    the disk of radius rc about each center gives the jet remainder.
    Scaling P(Lw)=exp(P(w)) then gives all derivatives at the original
    centers. No numerical differentiation is used.
    """
    R = self.R
    kk = 0
    scale = acb(1)
    yy = y
    while not old.ob(yy) <= arb("0.1"):
        yy /= R.L
        scale *= R.L
        kk += 1
        if kk > 100:
            raise ArithmeticError("Scaling did not terminate")
    rc = arb("0.12")
    assert old.ob(yy) + rc <= R.r0
    coeff = self.shift(R.p, yy, m)
    rem = R.C * (old.ob(yy) + rc) ** (R.n0 + 1)
    ctx.cap = m + 1
    s = acb_series([(coeff[j] + old.disc(acb(0), rem / rc ** j)) / scale ** j
                    for j in range(m + 1)])
    for _ in range(kk):
        s = s.exp()
    vals = s.coeffs()
    vals += [acb(0)] * (m + 1 - len(vals))
    for v in vals:
        finite(v)
    return vals

def sharp_P_jet(self, y, m):
    if old.ob(y) <= arb("0.12"):
        return self._P_jet_direkt(y, m)
    yc = acb(y.real.mid(), y.imag.mid())
    ry = (y.real.rad() ** 2 + y.imag.rad() ** 2).sqrt() * arb("1.00001") + arb(10) ** -200
    M2 = m + 24
    ph = self._P_jet_punkt(yc, M2)
    radius = old.amax(16 * ry, abs(yc).lower() / 40)
    mr = finite(self.R._grob_sup_kreis(self.R.P_ball, yc, radius, nbox=16))
    ratio = ry / radius
    assert ratio < 1
    delta = old.disc(acb(0), ry)
    out = []
    for j in range(m + 1):
        value = acb(0)
        for i in range(M2, j - 1, -1):
            value = value * delta + ph[i] * math.comb(i, j)
        # Cauchy: |p_i| <= mr/radius^i. Bound the derivative tail
        # using binom(n+l,j) <= binom(n,j)*binom(l+j,j).
        error = (mr / radius ** j * math.comb(M2 + 1, j)
                 * ratio ** (M2 + 1 - j) / (1 - ratio) ** (j + 1))
        out.append(finite(value + old.disc(acb(0), error)))
    return out

old.Rand.P_jet = sharp_P_jet

def angular_by_reversion(self, X, F, m=6):
    """Solve X(v(s))=X(0)*exp(i*s), then evaluate F(v(s))."""
    ctx.cap = m + 1
    xc = X.coeffs() + [acb(0)] * (m + 1 - len(X.coeffs()))
    target = [xc[0] * acb(0,1)**j / math.factorial(j) for j in range(m+1)]
    inverse = [acb(0)] * (m + 1)
    def compose_polynomial(series, argument):
        value = acb_series([0])
        for coefficient in reversed(series.coeffs()[:m+1]):
            value = value * argument + coefficient
        return value
    for j in range(1, m+1):
        composed = compose_polynomial(X, acb_series(inverse)).coeffs()
        value = composed[j] if j < len(composed) else acb(0)
        inverse[j] = (target[j] - value) / xc[1]
    result = compose_polynomial(F, acb_series(inverse)).coeffs()
    result += [acb(0)] * (m + 1 - len(result))
    return [finite(result[j] * math.factorial(j)) for j in range(m+1)]

# Retain the original angular differential operator for point jets.
# The optional reversion variant above is an unsuccessful experiment.

def cauchy_sup6(rd, xm, zm, zmr, half, distance):
    """Bound the sixth derivative on a real angular cell via point jets.

    Cauchy bounds the derivative tail of an angular Taylor series. The
    implicit-function disk is checked to stay in a certified sexp chart.
    """
    radius = min(arb("0.06"), distance / 3)
    Mtheta = None
    for _ in range(12):
        if not half < radius / 3:
            raise ArithmeticError("Angular cell too wide")
        xrad = rd.R.absL * (radius.exp() - 1)
        zc, zr = rd.K.slog_krawczyk(xm, zm, xrad)
        if zc is not None:
            im = abs(zc.imag)
            if im <= rd.K.H:
                chart_ok = im + zr < arb("0.45")
            else:
                chart_ok = zc.imag - zr > arb("0.36")
            if chart_ok:
                val, _ = rd.K.sexp_ball((zc + arb("0.5"), zr))
                Mtheta = finite(old.ob(val))
                break
        radius /= 2
    if Mtheta is None:
        raise ArithmeticError("No certified angular disk")
    Mjet = 26
    if abs(zm.imag) + zmr <= arb("0.42"):
        point = rd.t_gebiet(zm, zmr, Mjet)
    else:
        point, _, _ = rd.struktur(old.disc(zm, zmr), Mjet)
    center = finite(old.ob(point[6]))
    polynomial = sum((old.ob(point[j]) * half ** (j-6) / math.factorial(j-6)
                      for j in range(6,Mjet+1)),arb(0))
    ratio = half / radius
    tail = (math.factorial(6) * Mtheta / radius ** 6 * math.comb(Mjet+1,6)
            * ratio ** (Mjet+1-6) / (1-ratio) ** 7)
    result = finite(polynomial + tail)
    if not result <= arb("1.5") * center + 10:
        raise ArithmeticError("Derivative enclosure too wide")
    return result, center, tail

def run_segment(rd, side, cutoff):
    K, rho = rd.K, rd.R.absL
    # Finite-precision midpoint is a rational number, not the exact corner.
    # Its distance to arg L is carried explicitly into the remainder radius.
    corner = arb(rd.L.imag.mid())
    corner_error = old.ob(acb(rd.L.imag - corner))
    assert corner_error < arb("1e-180")
    sign = 1 if side == "right" else -1
    theta = arb(0) if sign == 1 else arb.pi()
    stop = corner - sign * cutoff
    w = acb("0.31") if sign == 1 else acb("-1.75")
    total = arb(0)
    cells = []
    tstart = time.time()
    while sign * (stop - theta) > 0:
        distance = abs(corner - theta)
        h = min(0.012, float(distance.lower()) / 12)
        remaining = sign * (stop - theta)
        if arb(h) > remaining:
            h = float(remaining.lower())
        if h < 1e-190:
            break
        accepted = False
        for tries in range(40):
            theta2 = theta + sign * arb(h)
            mid = (theta + theta2) / 2
            half = abs(theta2 - theta) / 2
            xm = rho * (acb(0, 1) * mid).exp()
            xr = rho * half * arb("1.00001")
            zc, zr = K.slog_krawczyk(xm, w, xr)
            if zc is None:
                h /= 2
                continue
            zm, zmr = K.slog_krawczyk(xm, zc, arb(0))
            if zm is None:
                h /= 2
                continue
            imlo = (zc.imag - zr).lower()
            imhi = (zc.imag + zr).upper()
            try:
                if distance > arb("0.12"):
                    sup6, c6, cauchy_tail = cauchy_sup6(rd,xm,zm,zmr,half,distance)
                    kind = "C"
                    s7 = arb(0)
                    term = finite(2 * half * sup6)
                    total += term
                    cells.append({"a":str(theta.mid()),"b":str(theta2.mid()),
                                  "method":kind,"c6_upper":str(c6.upper()),
                                  "sup6_upper":str(sup6.upper()),
                                  "cauchy_tail_upper":str(cauchy_tail.upper()),
                                  "integral_upper":str(term.upper())})
                    w, theta, accepted = zm, theta2, True
                    break
                if abs(imlo) <= arb("0.42") and abs(imhi) <= arb("0.42"):
                    kind = "T"
                    ball = rd.t_gebiet(zc, zr, 7)
                    point = rd.t_gebiet(zm, zmr, 7)
                elif imlo >= arb("0.39") and (zm.imag - zmr).lower() >= arb("0.39"):
                    kind = "S"
                    ball, _, _ = rd.struktur(old.disc(zc, zr), 7)
                    point, _, _ = rd.struktur(old.disc(zm, zmr), 7)
                else:
                    h /= 2
                    continue
                c6 = finite(old.ob(point[6]))
                s7 = finite(old.ob(ball[7]))
                if not half * s7 <= c6 + arb(1):
                    h /= 2
                    continue
                term = finite(2 * half * c6 + half ** 2 * s7)
                total += term
                cells.append({"a":str(theta.mid()),"b":str(theta2.mid()),
                              "method":kind,"c6_upper":str(c6.upper()),
                              "s7_upper":str(s7.upper()),"integral_upper":str(term.upper())})
                w = zm
                theta = theta2
                accepted = True
                break
            except (ArithmeticError, ValueError) as exc:
                if not cells and tries == 0:
                    print(f"First cell rejected: {type(exc).__name__}: {exc}",flush=True)
                h /= 2
        if not accepted:
            raise RuntimeError(f"Cell failed at {theta}, h={h}")
        if len(cells) % 100 == 0:
            print(f"{side}: {len(cells)} cells, theta={float(theta):.10g}, "
                  f"distance={float(abs(corner-theta)):.3g}, {kind}, integral<={total.str(8)}, "
                  f"elapsed={time.time()-tstart:.1f}s", flush=True)
    gap = abs(corner-theta) + corner_error
    # Enlarge to include the rational endpoint uncertainty.
    rem, info = old.rest_schranke(rd, w, gap * arb("1.00001"))
    finite(rem)
    bound = total + rem
    record = {"side":side,"precision_bits":ctx.prec,"corner_error":str(corner_error.upper()),
              "cutoff_upper":str(gap.upper()),"cell_sum_upper":str(total.upper()),
              "remainder_upper":str(rem.upper()),"total_upper":str(bound.upper()),
              "rest_info":info,"cells":cells,
              "status":"interval bound conditional on upstream P-V/structure certificates"}
    (OUT / f"k2_{side}.json").write_text(json.dumps(record,ensure_ascii=False,indent=1),encoding="utf-8")
    print(f"FINAL {side}: integral <= {bound.str(12)}, rest <= {rem.str(8)}",flush=True)
    return bound

if __name__ == "__main__":
    rd = old.Rand()
    # The upstream JSON rounds its bound to five significant digits.
    rd.dth = arb("4.4e-167")
    side = sys.argv[1] if len(sys.argv)>1 else "right"
    run_segment(rd, side, arb("1e-8"))
