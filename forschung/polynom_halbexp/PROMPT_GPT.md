# Research task: polynomial half-iterates of exp converge to Kneser's half-exponential

You are asked to work on an open research problem in complex analysis and approximation theory. A team has spent a long session on it and
reduced it to one hard missing piece. Your job is to **prove the conjecture**, or **prove the missing piece**, or find a **different route**
that avoids it. Take as long as you need, compute where it helps, and be completely honest about what you have actually proved.

---

## 1. The conjecture

For each degree N, look for a real polynomial

  p_N(x) = a₀ + a₁x + … + a_Nx^N  with  p_N(p_N(x)) = eˣ + O(x^{N+1})  at x = 0,

i.e. the Taylor coefficients of p_N∘p_N of order 0…N equal 1/k!. These are N+1 polynomial equations in N+1 unknowns. The relevant
solution branch starts at p₁(x) = x + ½ and is continued degree by degree: the padded p_N is the starting point for p_{N+1}, joined by a
homotopy (Newton from the padded p_N follows it).

Let L = 0.31813150520476413 + 1.3372357014306895i be the fixed point of exp in the upper half-plane (e^L = L), ρ := |L| = 1.3745570107,
and D := {|x| < ρ}. Let f_K := sexp_K(slog_K(x) + ½) be **Kneser's real half-iterate of exp** (f_K∘f_K = exp, holomorphic on D, singular
at L and L̄). In particular f_K(1) = T_e(½) = 1.6463542337511945809719240315921…, and f_K(0) = 0.4985632879411144346796…

**Conjecture.** p_N → f_K uniformly on D̄ (in particular every coefficient converges and p_N(1) → T_e(½)).

**Numerical status (very strong).**
- At N = 194, p_N(1) agrees with T_e(½) to 19 digits. The branch is followed up to N = 1500 with high precision.
- Precise error law: p_N − f_K ≈ c_N·δ₋₁ + conj, with δ∓1 := T′(A + ½)·e^{∓2πiA}, A = slog_K(x), and |c_N| ≍ N^{1−2Re γ} = N^{−7.894}.
  Here γ := 2πi/L = 4.44695 + 1.05794i.
- The boundary error on |x| = ρ decays like N^{−4.84}.
- A blind prediction from this law (sign change of p_N(1) − T_e(½) at N ≈ 215–218) was confirmed.
- max_{|x|=ρ}|p_N| ≈ 2.2327 for all N tested.

---

## 2. Known facts you may use (proved unless marked)

**(a) Kneser structure at L (proved, Kneser / Paulsen–Cowgill + our derivation).**
f_K − L = √L·t·(1 + Σc_m t^m) + Σ_{n≥1} d_n·t^{1+nγ}(1 + Σ_j e_{n,j}t^j), with t = x − L. The analytic part is the regular (Schröder)
half-iterate at L; the singular exponents are 1 + nγ + j, and there are no logarithms. f_K(L) = L and f_K′(L) = μ := √L, with |μ| = 1.1717.
Taylor coefficients at 0: a_k = 2Re[D₁L^{−k}Γ(k−1−γ)/(Γ(−1−γ)k!)]·(1 + O(1/k)), with D₁ = −0.0320033 − 0.0112441i, so that
|a_k| ≲ 71.57·k^{−6.447}·ρ^{−k}. f_K is C^{5.447} up to the boundary at L.

**(b) Uniqueness via a corridor (proved, modulo the numerically checked shape of slog(D)).**
Let V₀ := {z ∈ D : f_K(z) ∈ D}. If f is holomorphic on D, real on ℝ ∩ D, satisfies f∘f = exp near 0, and maps V₀ into D, then f = f_K.
Any other real half-iterate violates this, because its perturbation grows like |z − L|^{−4.447n} near L.
Consequence: if sup_D|p_N| ≤ B for all N, |p_N| ≤ ρ on each compact subset of V₀ for large N, and p_N(0) ≤ ρ₀ < ρ, then p_N → f_K
locally uniformly on D (Montel + uniqueness). So uniform **boundedness plus the corridor condition** would already prove the
locally uniform version.

**(c) Finite certification (rigorous, interval arithmetic).**
- The branch is certified by Krawczyk up to N = 193, with boxes 2^{−250}·ρ^{−k}. A naive box fails from N ≈ 106 on.
- Kneser's sexp for base e is certified to ±1.3e−168, and a₀…a₆₄₀ of f_K are certified with error·ρⁿ ≤ 2.9e−127.
- For N = 140…193: ‖p_N − Π_N f_K‖_{ℓ¹_ρ} ≤ 2.2e−21, where Π_N is defined in (d).

**(d) Linearization.** Define
- J_∞u := u∘f + (f′∘f)·u, with f = f_K;
- J_N := T_N J_∞|_{P_N}, where P_N = polynomials of degree ≤ N and T_N = truncation after degree N;
- Π_N u := J_N^{−1}T_N J_∞u.

In the weighted space ℓ²_ρ (coefficient k weighted by ρ^k):
- **Bulk coercivity, essentially rigorous:** σ_min(J_∞|_{P_N}) ≥ 0.6553 for all N, via a positive-real-part (Gårding-type) estimate
  Re⟨𝒯x, x⟩ ≥ 0.9359‖x‖² for a suitable boundary operator 𝒯. The true value is 0.7319 = min over the arc Γin of |f′∘f|. The only
  remaining input is a rigorous tail bound |a_kρ^k| ≤ C·k^{−6} beyond k = 640 (numerically C ≈ 124).
- Rectangular sections T_{N+d}J_∞|_{P_N} have σ_min bounded below only for d ≳ 0.02N.
- **Kantorovich-type ball arguments fail in every norm:** random directions are amplified by ≳ 1.17^N. Only coherent directions can work.

**(e) Conditional main theorem (proved, conditional on (H)).**
(Kneser structure) + (H) + (A) ⟹ p_N → f_K uniformly on D̄, with rate N^{κ−4.69}. Here:
- **(A):** p_{N₀} lies in the ball B_{N₀} = {sup_{D̄}|p − f_K| ≤ R_{N₀}}. This is a finite check: with the constants provable today it is
  needed for N₀ ≈ 318, and the certificate can be extended to N = 400.
- **(H):** a linear statement about the pole responses of the linearizations at all near-solutions P in the ball:
  - (H-i) invertibility;
  - (H-ii) |Π^P_N k_ζ(z) − k_ζ(z)| ≤ C N^κ e^{−εN g(ζ)}/dist(ζ, K) for z ∈ D̄ and ζ on the slits L(1+u), L̄(1+u) and on a level line Γ_d.
    Here k_ζ(z) = 1/(ζ − z), g is the Green function of the complement of a compact K ⊃ D̄, and Π^P is the secant analogue of Π_N;
  - (H-iii) head stability.
- **(H₀):** (H) at P = f_K alone. The extension (H₀) ⟹ (H) to coherent perturbations is also open, but numerically harmless.

The proof structure:
- f_K − Π_N f_K is controlled by a Cauchy integral over Γ_d and the slit banks. The slit contribution is
  ∫u^{5.447}N^κe^{−εNcu^{1/α}}u^{−1}du = O(N^{κ−α·5.447}), with outer angle απ of K at L.
- A homotopy λ ↦ (1−λ)f_K + Q with a covering argument gives existence and uniqueness of p*_N in the ball.
- Induction in N identifies p*_N with the branch.

Useful compact sets: W̄ = D̄ ∪ D̄(a₀, |L − a₀|), outer angle 158.934° at L, α = 0.882964; and K′ = D̄ ∪ f_K(D̄), outer angle 141.691° = π − arg μ.

---

## 3. The obstacle: corner modes at L (this is where the proof is stuck)

J_N^{−1} is **not** uniformly bounded. It has r(N) ≈ N/250 exponentially large singular values ("corner modes"):
- σ₁ ≍ N^{−0.82}·e^{cN}, with c = 0.0118–0.0120.
- The rest of the spectrum is bounded: σ_{r+1} → 1.3663 = 1/0.7319.
- They come in complex-conjugate pairs. The left singular vectors live at arg L, in a boundary layer of width ~1/N and height ~√N.
- The adjoint vectors sit on coefficient indices k ≈ 0.95N.
- They depend essentially only on μ = f′(L): a real cubic f with the same fixed-point data gives the same numbers to 2–3 digits.

The exact solution never excites them dangerously (the error law above is clean), but no one has proved why. In the reduced form,
everything comes down to:

**Hypothesis (H₀).** Write the adjoint corner modes y_i := σ_i v_i as Laurent polynomials Y_i(t) = Σ_{k≤N} y_{i,k}t^{−k} (t = x/ρ).
There is a fixed closed curve γ* through L/ρ and L̄/ρ such that

  |Y_i(t)| ≤ C·N^κ  for t ∈ γ*, all corner modes i, all N ≥ N₀.

The curve γ* runs on one side inside D (through the "bay" h(A-arc), h = branch of f^{−1} with h(L) = L, and the A-arc is the part of the
circle |x − a₀| = |L − a₀| outside D), and on the other side just outside D along Γin = {|x| = ρ, |f(x)| < ρ}.
- It is proved that (H₀) implies the needed pole-response bound.
- Measured up to N = 1000: C = e^{−4.06}, κ = 0.88.
- The **polynomial** form is essential. A bound e^{o(N)} is not enough, because the singularity of f_K at L only decays polynomially.

**Exact identity (proved).** With F̃(t) = f(ρt)/ρ, m̃ = (f′∘f)(ρ·), 𝒰(t) = Σu_kt^{−k−1}, and P a polynomial:

  T(t) := m̃(t)Y(t)/t + Σ_{F̃(z)=t} Y(z)/(zF̃′(z)) − 𝒰(t) + P(t) = O(t^{−N−2}) at ∞.

T is rational with poles only at 0 and F̃(0). This is equivalent to J_Nᵀy = u.

**Limit problem (numerical).** (1/N)·log|Y| and (1/N)·log|T| converge to a function H that is harmonic on a glued surface ℳ:
- ℳ is the plane minus the region between a "bay" curve B ⊂ D and the zero curve Z_T of T outside D, with B glued to Z_T by f (H = H∘h).
- ℳ has genus 0; H ~ −log|t| at 0 and at ∞; H = 0 on the seam; B is a free boundary.
- At L, ℳ is a cone: H ≈ Re(c₁ζ^{κ₁}) with κ₁ = 2πi/(2πi − log μ) = 1.118 − 0.032i.
- The free-boundary problem predicts c = 0.0118 ± 0.0004, the boundary profile to 3e−4, and ∂c/∂μ = 0.044 − 0.040i (measured 0.042 − 0.035i).
- The zeros of Y_N follow the harmonic measure of the domain bounded by B and a thin lens outside Γin; at N = 1000 they agree to KS
  distance 0.005.
- Near L the profiles collapse on the scale |t − L| ~ N^{−0.894} = N^{−1/Re κ₁}; prefactor N^{1/2}.

**Approaches that failed (do not repeat them unless you see what was missed):**
1. *Local arguments.* Via the identity, the bay bound is equivalent to a bound for T on the outer boundary, so the circularity is
   intrinsic. The corner zone needs angle-resolved, polynomially accurate data from outside.
2. *Variational / extremality in ℓ²_ρ.* ‖J_N‖ = e^{0.47N}. A test vector built from the model overlaps v₁ by 0.997, but its Rayleigh
   quotient is 10^{90}.
3. *Gonchar–Rakhmanov / Stahl in the plane.* The exact interpolants Q_n (Q_n∘f + (f′∘f)Q_n = O(xⁿ)) satisfy ∮Q_nψ_k ds = 0, but each
   condition couples a moment at f(0), transported by f, with a moment at 0 of the same index. That is not a system {s^jW(s)}. On ℳ the
   problem would have GR form, but ℳ exists only in the limit (gluing errors e^{−Nδ}, with δ → 0 at L).
4. *n-th root asymptotics.* It is proved that rate limits exist along subsequences and satisfy a "tropical balance": the maximum term
   rate in the identity is attained twice. Uniqueness of the limit configuration is open, and in any case gives only e^{o(N)}.
5. *Changing the weight radius ρ′ ≠ ρ.* Every ρ′ still shows exponential growth.
6. *Local double-scale / Wiener–Hopf model at L.* It is not well-posed on its own; the corner is a matching problem between the inner
   layer and the outer cone.

---

## 4. What we want from you

Primary goals (any one is a major result):
1. **Prove (H₀)**, or any substitute that suffices for the conditional theorem in §2(e). One candidate: strong asymptotics of the corner
   modes via a Riemann–Hilbert / Deift–Zhou-type analysis of a "welding" problem with a non-local perturbation and a local parametrix at
   the cone point L.
2. **Find a route that avoids (H₀).** Ideas we have not pursued, for inspiration only:
   - a nonlinear argument that shows directly that the truncation residual T_{>N}(p_N∘p_N) has no component along the corner modes;
   - a norm or conformal coordinate adapted to the corner (for example ζ = log((x−L)/(x−L̄)), in which D becomes an exact strip and f
     acts like a shift ζ ↦ ζ + log μ near L);
   - working at the level of Abel/slog functions instead of half-iterates;
   - proving uniform boundedness plus the corridor condition (§2(b)) directly, for example by a maximum-principle or subordination
     argument on the polynomial equation.
3. **A weaker but genuine theorem.** For example: for all sufficiently large N there is a solution of the truncated system near f_K, and
   these converge to f_K (without identifying them with the branch from p₁); or convergence locally uniformly on D.

If you come to believe the conjecture is false in some form (for example uniform convergence on D̄ versus only on D), say so and give the
evidence.

**Rules for your answer.**
- Label every claim **proved** / **proof sketch (with gaps named)** / **numerical** / **conjecture**. Do not present a sketch as a proof.
- State theorems precisely, with all hypotheses. When you cite literature, say whether you are sure of the exact statement.
- If you compute, use rigorous or high-precision arithmetic (for example arb / python-flint, mpmath). Known pitfalls:
  - Double-precision evaluation of p_N inside D is pure noise.
  - Coefficient data need about 0.07·N extra digits, because errors grow like 1.17^N.
  - python-flint's acb_series/arb_series silently truncate at ctx.cap = 10 terms unless you raise it.
- End with:
  1. a short summary of what is new and proved;
  2. the precise remaining gap, if any;
  3. the single most promising next step.
