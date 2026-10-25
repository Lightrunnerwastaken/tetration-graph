# e4: Pruefung der Jet-Identitaeten aus T'(x+1) = T'(x) T(x+1) bei x = -1 (Basis e, T(0)=1):
#   slog'(0) = 1/a1,  slog''(0) = (a1^2 - 2 a2)/a1^3,  T''(0) = 2 a2;  Basis 2: slog_2'(0) = ln2 / a1(2) (T'(0) = ln b T'(-1)).
import json, mpmath as mp
mp.mp.dps = 200
d = json.load(open('cache/konstanten_200.json'))
g = lambda k: mp.mpf(str(d[k]))
a1, a2 = g('e.taylor0_a1'), g('e.taylor0_a2')
for name, lhs, rhs in [("slog'(0)*a1 - 1", g("e.slog'(0)") * a1, 1),
                       ("slog''(0) - (a1^2-2a2)/a1^3", g("e.slog''(0)"), (a1**2 - 2*a2) / a1**3),
                       ("T''(0) - 2a2", g("e.T''(0)"), 2*a2)]:
    print('%-32s log10|Diff| = %.1f' % (name, float(mp.log10(abs(lhs - rhs) + mp.mpf(10)**-199))))
