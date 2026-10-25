# Fatou-Koordinate der Normalform F(m) = m*e^m am parabolischen Fixpunkt 0.
# Ansatz: Phi(m) = -1/m + beta*ln m + sum_{k>=1} c_k m^k,  Phi(F(m)) = Phi(m) + 1.
# Wegen ln F(m) = ln m + m (exakt) und 1/F(m) = e^{-m}/m wird die Gleichung rein algebraisch:
#   (1 - e^{-m})/m + beta*m + sum_k c_k (F(m)^k - m^k) = 1.
# Ausgabe: beta und c_k exakt (rational), dazu Wachstum |c_k| fuer die Divergenzfrage.
import sympy as sp, sys

K = int(sys.argv[1]) if len(sys.argv) > 1 else 24
m = sp.symbols('m')
N = K + 2
Fm = sp.series(m*sp.exp(m), m, 0, N + 1).removeO()
lhs0 = sp.series((1 - sp.exp(-m))/m, m, 0, N + 1).removeO()
beta = sp.Symbol('beta')
cs = sp.symbols('c1:%d' % (K + 1))
expr = lhs0 + beta*m - 1
powF = [sp.Integer(1)]
for k in range(1, K + 1):
    powF.append(sp.expand(sp.series(powF[-1]*Fm, m, 0, N + 1).removeO()))
for k in range(1, K + 1):
    expr += cs[k-1]*(powF[k] - m**k)
expr = sp.expand(expr)
sol = {}
# Ordnung m^1 bestimmt beta, Ordnung m^{k+1} bestimmt c_k
e1 = expr.coeff(m, 1)
sol[beta] = sp.solve(e1, beta)[0]
for k in range(1, K + 1):
    ek = sp.expand(expr.coeff(m, k + 1).subs(sol))
    sol[cs[k-1]] = sp.solve(ek, cs[k-1])[0]
print('beta =', sol[beta])
for k in range(1, K + 1):
    c = sol[cs[k-1]]
    print(k, c, float(c))
# Wachstum: |c_k| * (2 pi)^k / (k-1)!  sollte bei Gevrey-1 mit Radius 2pi beschraenkt sein
print('--- Wachstum')
for k in range(2, K + 1):
    c = sol[cs[k-1]]
    if c != 0:
        print(k, float(abs(c)*(2*sp.pi)**k/sp.factorial(k-1)))
