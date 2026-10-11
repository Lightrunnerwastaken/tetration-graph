"""η₆ und e₆: grösste Basis b, für die die (Kneser-)Hexation H_b einen reellen Fixpunkt hat; e₆ = Berührstelle (H_b(e₆) = e₆, H_b'(e₆) = 1).

 Für jede Basis: hexinit (Fork, exp-079), grobes Gitter für H_b(x) − x (Überläufe per iferr abgefangen), dann Newton auf H_b'(x) = 1
 (zentrale Differenzen), m(b) = H_b(x*) − x*. Sekanten-Verfahren in b auf m(b) = 0. Aufruf: eta6.py b0 b1 [dps] [schritte]
"""
import sys, time, json
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe

b0, b1 = sys.argv[1], sys.argv[2]
dps = int(sys.argv[3]) if len(sys.argv) > 3 else 30
schritte = int(sys.argv[4]) if len(sys.argv) > 4 else 8
mp.mp.dps = dps + 10
t0 = time.time()
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=dps)


def log(*a):
    print(f"[{time.time() - t0:6.1f}s]", *a, flush=True)


def m_von(b):
    bs = mp.nstr(b, dps + 5)
    kopf = gp.eval_batch(bs, ["hexinit(sexp(1),0,1)", "hexreg", "hexerr"])
    if int(kopf[1].real) != 2:
        raise RuntimeError(f"b = {bs}: keine Kneser-Hexation (hexreg {kopf[1]})")
    H = lambda xs: [w.real for w in gp.eval_batch(bs, [f"iferr(hex({mp.nstr(x, dps + 5)}),E,10^100)" for x in xs])]
    gitter = [mp.mpf(2) + mp.mpf(k) / 8 for k in range(33)]          # 2 … 6
    g = [h - x for h, x in zip(H(gitter), gitter)]
    i = min(range(len(g)), key=lambda j: g[j])
    x = gitter[i]
    h = mp.mpf(10) ** -(dps // 2 - 2)
    for _ in range(40):
        hm, h0, hp = H([x - h, x, x + h])
        d1 = (hp - hm) / (2 * h)
        d2 = (hp - 2 * h0 + hm) / h ** 2
        dx = (d1 - 1) / d2
        dx = max(min(dx, mp.mpf("0.25")), mp.mpf("-0.25"))
        x -= dx
        if abs(dx) < mp.mpf(10) ** -(dps // 2 - 4):
            break
    m = H([x])[0] - x
    return m, x, kopf[2].real, d1


protokoll = []
try:
    bs_, ms = [mp.mpf(b0), mp.mpf(b1)], []
    for b in bs_:
        m, x, err, d1 = m_von(b)
        ms.append(m)
        log(f"b = {mp.nstr(b, 20)}: Stellen {mp.nstr(err, 4)}, x* = {mp.nstr(x, 18)}, H'(x*) = {mp.nstr(d1, 8)}, m = {mp.nstr(m, 10)}")
        protokoll.append([str(b), str(x), str(m), float(err)])
    for s in range(schritte):
        b = bs_[-1] - ms[-1] * (bs_[-1] - bs_[-2]) / (ms[-1] - ms[-2])
        m, x, err, d1 = m_von(b)
        bs_.append(b); ms.append(m)
        log(f"b = {mp.nstr(b, dps)}: Stellen {mp.nstr(err, 4)}, x* = {mp.nstr(x, 20)}, m = {mp.nstr(m, 6)}")
        protokoll.append([str(b), str(x), str(m), float(err)])
        if abs(bs_[-1] - bs_[-2]) < mp.mpf(10) ** -(dps - 8):
            break
finally:
    gp.close()
    json.dump(protokoll, open(f"eta6_{dps}.json", "w"), indent=1)
if len(ms) > 2:
    log("η₆ ≈", mp.nstr(bs_[-1], dps - 8), "  e₆ ≈", mp.nstr(x, dps // 2 - 4))
