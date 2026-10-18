"""Leiter auf [-1, 2] (negative Hoehen). python neg_lauf.py b R dps grad"""
import sys, json, os
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import koenigs_leiter as KL
b, R = sys.argv[1], int(sys.argv[2]); mp.mp.dps = int(sys.argv[3]); g = int(sys.argv[4])
out, _ = KL.leiter(mp.mpf(eval(b, {'sqrt': mp.sqrt})), R, grad=g, X0=-1)
for r in out:
    print(r['r'] + 1, *(mp.nstr(mp.mpf(r[f'neg_{t}']), 10) for t in ('-0.75', '-0.5', '-0.25')), 'res_neg', r['resid_neg'])
fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'daten', f'koenigs_b{b.replace("(","").replace(")","")}_neg_R{R}_d{mp.mp.dps}_g{g}.json')
json.dump(dict(b=b, dps=mp.mp.dps, grad=g, X0=-1, rows=out), open(fn, 'w'), indent=1)
