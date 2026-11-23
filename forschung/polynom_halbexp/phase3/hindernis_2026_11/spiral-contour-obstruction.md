# A harmonic-measure obstruction for the adjoint heads

For the standard \(170^\circ\) contour, the two proposed estimates (A-lok) and (FA) force the normalized adjoint heads to grow faster than every power of \(n\) somewhere on each circle of the proposed corner scale. This gives a compatibility test for the proof of convergence of the polynomial square roots. The convergence of the nonlinear branch remains open.

## 1. The estimate on the standard contour

Let \(\mathcal C\) be the fixed Jordan contour surrounding \(0\) in the problem sheet. Write \(U\) for its exterior on the Riemann sphere. Suppose that near \(L\) the two sides of \(U\) are logarithmic spirals. In terms of \(Z=(t-L)/\rho\), their equations on a continuous angular lift are

\[
\arg Z-T\log|Z|=\theta_-,\ \theta_+,\qquad
\theta_+-\theta_-=\Theta.
\]

The interval between these values describes \(U\). Thus \(\Theta=170^\circ\) for the standard window from \(-90^\circ\) to \(80^\circ\). The approximately \(132^\circ\) angle belongs to the enlarged approximation set \(E\). Here the domain is the exterior of \(\mathcal C\) itself. [Contour and hypotheses, §§1–2](https://github.com/Lightrunnerwastaken/tetration-graph/blob/main/forschung/polynom_halbexp/phase3/HAUPTSATZ_KR0.md)

Let \(h_n(t)=\Lambda_n(t)/\Lambda_n(L)\), whenever this normalization exists. Each \(h_n\) is a polynomial in \(1/t\) of degree at most \(n\), is holomorphic at infinity, and satisfies \(h_n(L)=1\). Use the following consequences of (A-lok) and (FA) for all sufficiently large \(n\).

* On \(\mathcal C\), \(\log |h_n|\le b_n\), where \(0\le b_n=O(\log(n+1))\).
* On a fixed closed boundary arc \(F\) of positive length away from both corners, \(\log |h_n|\le-\eta n\) for some \(\eta>0\).

The estimates on the two conjugate corners give the same normalization in modulus because the heads have real coefficients. Any nondegenerate smooth subarc of the far part with the stated bound can serve as \(F\).

**Proposition.** Put \(p=\pi(1+T^2)/\Theta\), with \(\Theta\) in radians. Fix \(A>0\) and \(R>0\). If \(Ap<1\), then for some \(c_R>0\) and all sufficiently large \(n\),

\[
\log\max_{|t-L|=R\rho n^{-A}} |h_n(t)|
\ \ge\ c_R n^{1-Ap}.
\tag{1}
\]

**Proof.** Let \(\omega(t)\) be the harmonic measure of \(F\) in \(U\). Subharmonic comparison in \(U\), including the regular point at infinity, gives

\[
\log |h_n(t)|\le b_n-(b_n+\eta n)\omega(t).
\tag{2}
\]

The exponent of \(\omega\) at the corner can be computed exactly. Take the branch of \(\log Z\) in the spiral sector and set

\[
w=e^{-i\pi\theta_-/\Theta}
\exp\!\left(\frac{\pi}{\Theta}(1-iT)\log Z\right).
\]

Its argument is \(\pi(\arg Z-T\log|Z|-\theta_-)/\Theta\). Its modulus is \(|Z|^p\) times a factor bounded above and below by positive constants, since \(\arg Z-T\log|Z|\) ranges over a bounded interval. This sends a neighborhood of the tip in \(U\) onto a neighborhood of \(0\) in the upper half-plane. There \(\omega\) is positive harmonic and vanishes on the adjacent real interval. Reflection and the boundary point lemma give \(\omega\asymp\operatorname{Im}w\) near \(0\). In particular, on the middle half of the angular window,

\[
\omega(t)\ge c|Z|^p.
\tag{3}
\]

On the circle of radius \(R\rho n^{-A}\), this middle window occupies a fixed fraction \(q=\Theta/(4\pi)>0\) of the entire circle. Equations (2) and (3) bound \(\log |h_n|\) there by \(b_n-c'_R n^{1-Ap}\). The small disk about \(L\) excludes the only possible pole, \(0\). Writing \(M_n\) for the maximum of \(|h_n|\) on the circle, the submean inequality therefore gives

\[
0=\log|h_n(L)|
\le q(b_n-c'_R n^{1-Ap})+(1-q)\log M_n .
\]

Since \(1-Ap>0\), this proves (1). The proof uses no properties of the Gamma model.

## 2. Consequence for the proposed corner profile

For \(a=1+iL/(4\pi)\), \(A=\operatorname{Re}a=0.8935861640\ldots\) and \(T=\operatorname{Im}a/A=0.0283309006\ldots\). The critical opening is

\[
\Theta_*=\pi A(1+T^2)=160.97461055\ldots^\circ .
\]

At \(170^\circ\) one has \(1-Ap=0.053090526\ldots\). Equation (1) rules out even a polynomial upper bound on the normalized heads on a single circle \(|t-L|=R\rho n^{-A}\) if both contour estimates hold. In particular the functions \(u\mapsto h_n(L+\rho n^{-a}u)\) cannot converge locally uniformly near \(u=0\) to the proposed Gamma profile, or to any finite holomorphic function. The rotation \(n^{-i\operatorname{Im}a}\) leaves the circular maximum unchanged. The detailed derivation constructs the Gamma model and distinguishes it from the matching estimate for the actual finite heads. [Corner analysis, §§19, 27 and 28](https://github.com/Lightrunnerwastaken/tetration-graph/blob/main/forschung/polynom_halbexp/phase3/k_kern/BERICHT.md)

There is also a quantitative converse. Assume (A-lok) and \(M_n\le C n^B\) on one such circle. For an index at which the whole far part has supremum \(m_n<1\), apply the same proof with \(-\eta n\) replaced by \(\log m_n\). It gives

\[
-\log m_n \ \le\ C_R n^{Ap}\log(n+1).
\]

At \(170^\circ\), \(Ap=0.946909474\ldots<1\), so (FA) would fail. This is a conditional exclusion of (FA) for the actual heads. The two hypotheses alone do not identify which estimate, or which matching claim, eventually fails.

For the narrower window from \(-90^\circ\) to \(60^\circ\), \(\Theta=150^\circ\) and \(1-Ap<0\). The argument gives no contradiction there.

## 3. An unconditional exclusion for windows wider than a half-plane

If \(\Theta>\pi\), the two contour estimates are incompatible for **every** sequence of rational functions of the displayed form with \(h_n(L)=1\). No bound on rescaled heads is needed.

Indeed, let \(I\) be the bounded interior of \(\mathcal C\) and \(G_I(t,0)\) its positive Green function with pole at \(0\). The pole order of \(h_n\) gives \(\log|h_n|\le b_n+nG_I(t,0)\) in \(I\) by the maximum principle. In the interior spiral sector the same straightening argument gives

\[
G_I(t,0)\ \le\ C|Z|^{p_i},\qquad p_i=\frac{\pi(1+T^2)}{2\pi-\Theta}.
\]

Outside, (2) bounds \(\log|h_n|\) by \(b_n\) everywhere and by \(b_n-cn r^p\) on the middle angular window of \(|Z|=r\). Thus the full circle mean is at most \(b_n+Cnr^{p_i}-c'nr^p\). For \(\Theta>\pi\), \(p_i>p\). Choose \(r>0\) small and fixed so that \(Cr^{p_i}<c'r^p/2\), and then let \(n\) grow. The mean becomes negative, contradicting \(h_n(L)=1\). In particular the exploratory window from \(-110^\circ\) to \(95^\circ\) in §35.1 of the corner analysis cannot satisfy both estimates.

## 4. What the recorded data distinguish

The cubic-model samples of the far lens and outer circular arc (the column labelled **R_A fern**, with the connecting arcs recorded separately) give the following values. [Verification report, §3.22](https://github.com/Lightrunnerwastaken/tetration-graph/blob/main/forschung/polynom_halbexp/phase3/v_pruefstand/BERICHT.md)

| \(n\) | Sampled far maximum \(m\) | \(n^{5.734123}m\) |
|---:|---:|---:|
| 1400 | \(1.8\times10^{-7}\) | \(1.975\times10^{11}\) |
| 1800 | \(4.2\times10^{-8}\) | \(1.947\times10^{11}\) |
| 2000 | \(2.3\times10^{-8}\) | \(1.951\times10^{11}\) |

The exponent here is \(\operatorname{Re}(a(2+2\pi i/L))\), computed from the first adjoint Mellin pole independently of these three values. The agreement is a finite-degree observation based on rounded samples. It neither supplies the circle bound in §2 nor proves a power law for a contour supremum.

For the linear convergence theorem, polynomial control on the entire contour would already suffice with an appropriate exponent. In fact the rank-one estimate bounds the kernel increment by \(|L|G_n\sup_{\mathcal C}|h_n|\). If \(G_n=O(n^{\kappa_G})\) and \(\sup_{\mathcal C}|h_n|=O(n^b)\) with \(\kappa_G+b<-1\), the increments are summable in that estimate. Thus losing strict far-field exponential decay would still allow the proposed linear approximation rate. The polynomial bounds and the nonlinear family estimate remain to be proved.
