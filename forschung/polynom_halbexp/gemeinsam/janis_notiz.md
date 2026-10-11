# Functional Square Root of the Exponential and Kneser Tetration

## Numerical construction via truncated power series

### 1. Problem

We seek a real analytic function \(f\) satisfying

\[
f(f(x)) = e^x.
\]

Such a function is a **functional square root** (or half-iterate) of the exponential map.

If \(f\) coincides with the Kneser half-iterate of the exponential, then

\[
f(x)=\exp^{\circ 1/2}(x),
\]

and therefore

\[
f(1)=T_e\!\left(\frac12\right),
\]

where \(T_e\) denotes Kneser's real tetration to base \(e\).

The numerical value used for comparison is

\[
T_e\!\left(\frac12\right)
=
1.6463542337511945809719240315921146159\ldots
\]

The purpose of this experiment is to investigate whether a purely local power-series construction of a solution of

\[
f\circ f = \exp
\]

selects the Kneser half-iterate.

---

## 2. Truncated power-series ansatz

For degree \(N\), define

\[
f_N(x)
=
a_0+a_1x+a_2x^2+\cdots+a_Nx^N.
\]

Substitute \(f_N(x)\) into \(f_N\) itself:

\[
f_N(f_N(x)).
\]

Expand this composition as a Taylor series around \(x=0\), and compare it with

\[
e^x
=
1+x+\frac{x^2}{2!}
+\frac{x^3}{3!}
+\cdots.
\]

The coefficients are chosen so that

\[
f_N(f_N(x))
=
\sum_{k=0}^{N}\frac{x^k}{k!}
+
O(x^{N+1}).
\]

Equivalently, for every

\[
k=0,1,\ldots,N,
\]

the coefficient of \(x^k\) in \(f_N(f_N(x))\) is required to equal

\[
\frac1{k!}.
\]

This gives \(N+1\) nonlinear equations for the \(N+1\) unknown coefficients

\[
a_0,a_1,\ldots,a_N.
\]

Because the nonlinear system can have several solution branches, the branch used here is the real branch that can be followed continuously from the low-degree solutions.

---

## 3. Degree \(N=1\)

Let

\[
f_1(x)=a_0+a_1x.
\]

Then

\[
f_1(f_1(x))
=
a_0+a_1a_0+a_1^2x.
\]

Matching

\[
1+x
\]

gives

\[
a_1^2=1,
\]

and on the relevant real branch

\[
a_1=1.
\]

Then

\[
a_0+a_0=1,
\]

so

\[
a_0=\frac12.
\]

Thus

\[
f_1(x)=\frac12+x,
\]

and

\[
f_1(1)=1.5.
\]

---

## 4. Degree \(N=2\)

Take

\[
f_2(x)=a_0+a_1x+a_2x^2.
\]

Matching the coefficients of

\[
f_2(f_2(x))
\]

with those of

\[
1+x+\frac{x^2}{2}
\]

gives the system

\[
a_0+a_0a_1+a_0^2a_2=1,
\]

\[
a_1^2+2a_0a_1a_2=1,
\]

\[
a_1a_2+a_1^2a_2+2a_0a_2^2=\frac12.
\]

The relevant real solution branch gives approximately

\[
a_0=0.497894079065,
\]

\[
a_1=0.878112905194,
\]

\[
a_2=0.261795456736.
\]

Hence

\[
f_2(1)
=
a_0+a_1+a_2
=
1.637802440995079.
\]

Compared with Kneser's value,

\[
T_e\!\left(\frac12\right)
=
1.64635423375119458097\ldots,
\]

the error is approximately

\[
-8.55\times10^{-3}.
\]

---

## 5. Degree \(N=3\)

For

\[
f_3(x)
=
a_0+a_1x+a_2x^2+a_3x^3,
\]

the relevant numerical solution is approximately

\[
f_3(x)
\approx
0.498573862009
+0.876383165884x
+0.247391361453x^2
+0.0241167272112x^3.
\]

Therefore

\[
f_3(1)
=
1.646465116557290.
\]

The difference from Kneser's value is

\[
f_3(1)-T_e\!\left(\frac12\right)
\approx
1.11\times10^{-4}.
\]

Already at cubic order, the approximation is very close.

---

## 6. Degree \(N=5\)

For degree five, the computed polynomial is approximately

\[
\begin{aligned}
f_5(x)\approx{}&
0.498563604143\\
&+0.876334764185x\\
&+0.247546649625x^2\\
&+0.024575890794x^3\\
&-0.000932010410x^4\\
&+0.000265441394x^5.
\end{aligned}
\]

By construction,

\[
f_5(f_5(x))
=
1+x+\frac{x^2}{2}
+\frac{x^3}{6}
+\frac{x^4}{24}
+\frac{x^5}{120}
+O(x^6).
\]

Evaluating at \(x=1\),

\[
f_5(1)
=
1.646354339729278.
\]

The difference from Kneser's value is only

\[
f_5(1)-T_e\!\left(\frac12\right)
\approx
1.06\times10^{-7}.
\]

Thus a polynomial determined only from the first six Taylor coefficients of

\[
f(f(x))=e^x
\]

already reproduces the Kneser half-tetration value at \(x=1\) to roughly six to seven decimal places.

---

## 7. Numerical convergence

The same solution branch was continued to higher degree.

| Degree \(N\) | \(f_N(1)\) | Error \(f_N(1)-T_e(1/2)\) |
|---:|---:|---:|
| 1 | 1.500000000000000 | \(-1.46\times10^{-1}\) |
| 2 | 1.637802440995079 | \(-8.55\times10^{-3}\) |
| 3 | 1.646465116557290 | \(+1.11\times10^{-4}\) |
| 4 | 1.646331568406921 | \(-2.27\times10^{-5}\) |
| 5 | 1.646354339729278 | \(+1.06\times10^{-7}\) |
| 8 | 1.646354176216273 | \(-5.75\times10^{-8}\) |
| 10 | 1.646354226766234 | \(-6.98\times10^{-9}\) |
| 12 | 1.646354233382163 | \(-3.69\times10^{-10}\) |
| 14 | 1.646354233959211 | \(+2.08\times10^{-10}\) |
| 15 | 1.646354233718347 | \(-3.28\times10^{-11}\) |

For comparison,

\[
T_e\!\left(\frac12\right)
=
1.6463542337511945809719240315921146159\ldots
\]

The degree-15 approximation gives

\[
f_{15}(1)
=
1.646354233718347,
\]

which agrees with the Kneser value to about ten decimal digits.

The small oscillation in the final digits at higher degrees is plausibly caused by numerical conditioning of the nonlinear coefficient system rather than by failure of convergence. Higher-precision arithmetic would be needed to study the asymptotic behavior reliably.

---

# 8. Conjecture

## Conjecture: power-series square-root branch equals the Kneser half-iterate

Let

\[
f_N(x)
=
\sum_{k=0}^{N}a_k^{(N)}x^k
\]

be the real solution branch obtained by solving

\[
f_N(f_N(x))
=
e^x+O(x^{N+1})
\]

around \(x=0\), with the branch chosen continuously from the low-order solution

\[
f_1(x)=x+\frac12.
\]

Then the sequence \(f_N\) converges, in an appropriate analytic sense, to a function \(f\) satisfying

\[
f(f(x))=e^x,
\]

and this limiting function is the Kneser half-iterate of the exponential:

\[
\boxed{
f(x)=\exp^{\circ 1/2}(x)
}
\]

on the real domain where Kneser's iteration is defined.

In particular,

\[
\boxed{
\lim_{N\to\infty}f_N(1)
=
T_e\!\left(\frac12\right)
}
\]

and numerically

\[
\boxed{
\lim_{N\to\infty}f_N(1)
=
1.6463542337511945809719240315921146159\ldots
}
\]

---

## 9. Stronger possible form of the conjecture

A stronger statement would be that the coefficient vectors themselves converge:

\[
a_k^{(N)}
\longrightarrow
a_k
\qquad
(N\to\infty,\; k\text{ fixed}),
\]

so that

\[
f(x)
=
\sum_{k=0}^{\infty}a_kx^k
\]

has a nonzero radius of convergence and satisfies

\[
f(f(x))=e^x
\]

analytically.

One would then need to prove that this analytic functional square root satisfies the extra regularity, reality, monotonicity, normalization, and continuation properties that uniquely identify the Kneser fractional iterate.

---

## 10. Why this is nontrivial

The equation

\[
f(f(x))=e^x
\]

alone does **not** obviously determine a unique function.

Functional iteration equations can admit multiple branches, and a general functional square root of the exponential need not automatically be the Kneser half-iterate.

Therefore the main mathematical problem is not merely to prove convergence of the polynomial approximations, but to explain why this particular real analytic branch is selected and why its analytic continuation agrees with Kneser's canonical tetration.

Possible proof directions include:

1. Prove existence and convergence of the coefficient sequence \(a_k^{(N)}\).
2. Show that the limiting \(f\) is real analytic and strictly increasing on a suitable real interval.
3. Construct from \(f\) a continuous iteration group \(f_t\) satisfying
   \[
   f_s\circ f_t=f_{s+t}.
   \]
4. Relate this construction to an Abel function \(A\) satisfying
   \[
   A(e^x)=A(x)+1.
   \]
5. Prove that the resulting Abel function has the normalization and complex-analytic properties that characterize the Kneser solution.
6. Alternatively, establish uniqueness of the real analytic square root under suitable side conditions.

---

## 11. Interpretation

The numerical experiment is striking because no tetration formula, Abel function, Schröder function, or known Kneser construction is inserted into the calculation.

The only imposed condition is the local coefficient identity

\[
f(f(x))=e^x.
\]

Nevertheless, following one natural real solution branch produces

\[
f_N(1)
\]

that appears to converge rapidly to

\[
T_e\!\left(\frac12\right).
\]

This suggests that the canonical Kneser half-iterate may be recoverable directly as a distinguished analytic solution of the functional square-root equation

\[
\boxed{
f\circ f=\exp.
}
\]

At present, the numerical evidence supports this statement, but a proof would require a rigorous treatment of convergence, branch selection, analyticity, and uniqueness.
