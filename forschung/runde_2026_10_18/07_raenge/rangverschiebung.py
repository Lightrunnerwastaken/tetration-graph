"""Rangverschiebung zwischen zwei Koenigs-Leitern (z. B. ungestoert / theta-gestoert).
Fuer jede Groesse Q (delta_r=p_r-b, lam_r, H_{r+1}(0.25), H_{r+1}(0.5)) wird c_r so bestimmt, dass
Q'(r) = Q(r + c_r), mit kubischer Interpolation von log|Q - Q_inf| in r.
python rangverschiebung.py basis.json gestoert.json
"""
import sys, json
import mpmath as mp
mp.mp.dps = 30

A = json.load(open(sys.argv[1]))
B = json.load(open(sys.argv[2]))
b = mp.mpf(eval(A['b'], {'sqrt': mp.sqrt}))


def series(d, key):
    out = {}
    for r in d['rows']:
        v = mp.mpf(r[key])
        if key == 'p':
            v = v - b
        elif key.startswith('next'):
            v = b - v   # Abstand vom Plateau (Grenzwert b fuer x > b-1)
        out[r['r']] = mp.log(abs(v))
    return out


def shift(qa, qb, r):
    # finde t mit qa(t) = qb(r), Lagrange-Interpolation ueber 4 Nachbarn
    rs = sorted(qa)
    target = qb[r]
    def interp(t):
        # nimm die 4 naechsten Stuetzstellen
        near = sorted(rs, key=lambda s: abs(s - t))[:4]
        tot = 0
        for i in near:
            w = 1
            for j in near:
                if j != i:
                    w *= (t - j) / (i - j)
            tot += w * qa[i]
        return tot
    t = mp.findroot(lambda t: interp(t) - target, mp.mpf(r))
    return t - r


for key in ['p', 'lam', 'next_quarter', 'next_half']:
    qa, qb = series(A, key), series(B, key)
    cs = []
    for r in sorted(qb)[2:-1]:
        try:
            cs.append((r, shift(qa, qb, r)))
        except Exception:
            pass
    print(key.ljust(13), ' '.join(f'{r}:{mp.nstr(c,8)}' for r, c in cs))
