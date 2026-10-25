# Pruefung der Keime der linken Hyperoperationen in Log-Koordinaten (exakt, Bruchzahlen, schnell).
# T_k(z, s) := ln( e^z [k'] s ),  linke Raenge: y[k+1'](n+1) = (y[k+1']n)[k'] y,  y[k+1']0 = y.
# Behauptung: T_k(z, s) = z + s z^2 + (A_k s^2 + B_k s) z^3 + O(z^4) mit
#   A_4 = 1/2, B_4 = 0;  A_k = 1, B_k = k - 9/2 fuer k >= 5.
# Formale gebrochene Iterierte von F(m) = m e^m: F^[s](m) = m + s m^2 + (s^2 - s/2) m^3 + O(m^4)
# (Gruppeneigenschaft; fuer s = 1, 2 direkt nachgerechnet).
# Geprueft wird fuer ganze s = 1..4 durch explizites Zusammensetzen; die z^3-Koeffizienten sind
# Polynome vom Grad 2 in s, drei Werte legen sie fest, der vierte ist Gegenprobe.
from fractions import Fraction as Q

N = 4  # Reihen bis z^3 (Index 0..3)

def mul(a, b):
    c = [Q(0)]*N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i):
                c[i+j] += x*b[j]
    return c

def add(a, b):
    return [x + y for x, y in zip(a, b)]

def scal(k, a):
    return [k*x for x in a]

def const(k):
    return [Q(k)] + [Q(0)]*(N - 1)

def exp_s(a):
    # exp einer Reihe mit a[0] = 0
    assert a[0] == 0
    r = const(1); t = const(1)
    for n in range(1, N):
        t = scal(Q(1, n), mul(t, a)); r = add(r, t)
    return r

Z = [Q(0), Q(1), Q(0), Q(0)]

def Fs(s, w):
    # F^[s](w) = w + s w^2 + (s^2 - s/2) w^3 ; s darf eine Reihe sein
    w2 = mul(w, w); w3 = mul(w2, w)
    return add(add(w, mul(s, w2)), mul(add(mul(s, s), scal(Q(-1, 2), s)), w3))

def T4(w, s):
    return mul(w, exp_s(mul(s, w)))

def T5(w, s):
    # ln(e^w [5'] s) = F^[s](w e^w) e^{-w} = F^[s+1](w) e^{-w}
    return mul(Fs(add(s, const(1)), w), exp_s(scal(-1, w)))

def iterate(T, n, z):
    # ln(e^z [k+1'] n): n-fach w -> T(w, y) mit y = e^z fest, Start w = z
    y = exp_s(z)
    w = z
    for _ in range(n):
        w = T(w, y)
    return w

def show(name, ser, n):
    print(name, 'n=%d' % n, [str(x) for x in ser])

# Gegenprobe der Formel F^[s] fuer s=2: F(F(z)) direkt
Fz = mul(Z, exp_s(Z))
print('F(F(z)) =', [str(x) for x in mul(Fz, exp_s(Fz))], ' Formel s=2:', [str(x) for x in Fs(const(2), Z)])

print('Rang 4\': T4(z,s) z^3-Koeff. = s^2/2 (direkt).')
T6 = lambda w, s: None
for n in range(1, 5):
    t5 = iterate(T4, n, Z)          # ln(e^z [5'] n)
    print("T5(z,%d)" % n, [str(x) for x in t5], ' Vorhersage z^3:', Q(n*n) + Q(n, 2))
for n in range(1, 5):
    t6 = iterate(T5, n, Z)          # ln(e^z [6'] n)
    print("T6(z,%d)" % n, [str(x) for x in t6], ' Vorhersage z^3:', Q(n*n) + Q(3*n, 2))

# Rang 7: T6(w, s) fuer Reihen-s braucht die gebrochene Iterierte von T5(., y); bis z^3 genuegt
# die Formel fuer einen Keim w + a w^2 + b w^3:  s-fach = w + s a w^2 + (s b + s(s-1) a^2) w^3.
def iter_germ(a, b, s, w):
    w2 = mul(w, w); w3 = mul(w2, w)
    return add(add(w, mul(mul(s, a), w2)), mul(add(mul(s, b), mul(mul(s, add(s, const(-1))), mul(a, a))), w3))

def T6s(w, s):
    # T6(w, s) = s-fach Iteration von u -> T5(u, e^w) ab u = w; Keim T5(u,y) = u + y u^2 + (y^2 + y/2) u^3
    y = exp_s(w)
    return iter_germ(y, add(mul(y, y), scal(Q(1, 2), y)), s, w)

for n in range(1, 5):
    t7 = iterate(T6s, n, Z)
    print("T7(z,%d)" % n, [str(x) for x in t7], ' Vorhersage z^3:', Q(n*n) + Q(5*n, 2))
