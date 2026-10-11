"""b40 (Runde 4): Prüfung der Ableitungsidentität (exakt, rationale Arithmetik, Reihen mod x^M):
  D²Φ(f)[a,b] = −J_∞[(ab)'/(2f')] + (1/(2f'))·[b r_a' − b' r_a + a r_b' − a' r_b − (g'/g)(a r_b + b r_a)] + (r_a r_b)'/(2 g f'),
  r_a := J_∞a = a∘f + g a,  g := f'∘f,  D²Φ(f)[a,b] = (a'∘f) b + (b'∘f) a + (f''∘f) a b.
Für beliebige Polynome f (f(0) ≠ 0, f'(0) ≠ 0), a, b — die Identität ist rein algebraisch."""
import sys, random
import sympy as sp
sys.stdout.reconfigure(encoding="utf-8")
x, y = sp.symbols("x y")
M = 14
random.seed(3)


def rnd_poly(deg, c0=None):
    co = [sp.Rational(random.randint(-9, 9), random.randint(1, 5)) for _ in range(deg + 1)]
    if c0 is not None:
        co[0] = c0
    return sum(c * x**k for k, c in enumerate(co))


def ser(e):
    return sp.series(e, x, 0, M).removeO()


def comp(p, q):          # p∘q für Polynome p, q
    return sp.expand(p.subs(x, q))


ok = True
for versuch in range(3):
    f = rnd_poly(4, sp.Rational(1, 2)) + x          # f(0) = 1/2, f'(0) ≠ 0 (fast sicher)
    if sp.diff(f, x).subs(x, 0) == 0:
        continue
    a = rnd_poly(5); b = rnd_poly(5)
    fp = sp.diff(f, x); fpp = sp.diff(fp, x)
    g = comp(fp, f)
    J = lambda u: sp.expand(comp(u, f) + g * u)
    ra, rb = J(a), J(b)
    links = sp.expand(comp(sp.diff(a, x), f) * b + comp(sp.diff(b, x), f) * a + comp(fpp, f) * a * b)
    v = sp.diff(a * b, x) / (2 * fp)
    Jv = ser(v.subs(x, f)) if False else None
    # J_∞ v als Reihe: v∘f braucht v als Reihe um f(0) — v ist rational, also exakt einsetzbar
    Jv = sp.together(v.subs(x, f) + g * v)
    gp = sp.diff(g, x)
    rechts = -Jv + (b * sp.diff(ra, x) - sp.diff(b, x) * ra + a * sp.diff(rb, x) - sp.diff(a, x) * rb
                    - gp / g * (a * rb + b * ra)) / (2 * fp) + sp.diff(ra * rb, x) / (2 * g * fp)
    diff = sp.simplify(sp.together(links - rechts))
    print(f"Versuch {versuch}: links − rechts = {diff}")
    ok = ok and diff == 0
print("Identität exakt (als rationale Funktionen):", ok)
