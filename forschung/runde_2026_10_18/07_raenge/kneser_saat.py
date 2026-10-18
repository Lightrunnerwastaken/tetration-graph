"""Rangleiter fuer e^(1/e) < b < eta_4: Saat H_4 = Kneser-Tetration (gp-Motor, eigener Zwischenspeicher),
ab H_5 kanonische Koenigs-Stufen (koenigs_leiter.leiter).
python kneser_saat.py b R dps grad
"""
import sys, json, os, time
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import koenigs_leiter as KL

GP_SRC = Path(r"D:\Users\janis\Documents\Tetration\gp-tetration\src")


def kneser_werte(basis, xs, dps):
    cf = HERE / 'cache' / f'sexp_{basis}_{dps}.json'
    cf.parent.mkdir(exist_ok=True)
    c = json.loads(cf.read_text()) if cf.is_file() else {}
    keys = [mp.nstr(x, dps + 8, min_fixed=-mp.inf, max_fixed=mp.inf) for x in xs]
    fehlt = [k for k in keys if k not in c]
    if fehlt:
        sys.path.insert(0, str(GP_SRC))
        from fatou_backend import FatouGP
        m = FatouGP(dps=dps, fatou_gp="fork", eval_timeout=500.0, init_timeout=500.0)
        t = time.time()
        s = m.session(basis)
        res = s.eval_batch([f"sexp({k})" for k in fehlt])
        print('gp', time.time() - t, 's', flush=True)
        for k, w in zip(fehlt, res):
            c[k] = mp.nstr(mp.re(w), dps + 8)
        cf.write_text(json.dumps(c, indent=0))
    return [mp.mpf(c[k]) for k in keys]


if __name__ == '__main__':
    basis = sys.argv[1]
    R = int(sys.argv[2])
    mp.mp.dps = int(sys.argv[3])
    grad = int(sys.argv[4])
    b = mp.mpf(basis)
    X = 2
    n = grad
    ts = [mp.cos(mp.pi * (k + mp.mpf(1) / 2) / n) for k in range(n)]
    xs = [X * (t + 1) / 2 for t in ts]
    vals = kneser_werte(basis, xs, mp.mp.dps)
    table = dict(zip([mp.nstr(x, 25) for x in xs], vals))
    seed = KL.Cheb(lambda x: table[mp.nstr(x, 25)], 0, X, n)
    # Kontrolle: Normierung und Funktionalgleichung der Saat
    print('H4(0)=', mp.nstr(seed(0), 15), ' H4(1)=', mp.nstr(seed(1), 15), ' H4(1.5)-b^H4(0.5)=',
          mp.nstr(seed(mp.mpf(1.5)) - b ** seed(mp.mpf(0.5)), 3))
    out, _ = KL.leiter(b, R, X=X, grad=grad, seed4=seed)
    fn = HERE / 'daten' / f'koenigs_b{basis}_kneser_R{R}_d{mp.mp.dps}_g{grad}.json'
    json.dump(dict(b=basis, dps=mp.mp.dps, grad=grad, seed='kneser', rows=out), open(fn, 'w'), indent=1)
    print('->', fn)
