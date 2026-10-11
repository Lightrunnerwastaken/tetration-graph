# Open problem: corner modes of finite sections of a composition operator (Hypothesis (H₀))

Status 11 Oct 2026 (English version of OFFENES_PROBLEM.md, same content). Labels: **proved** / **numerical** / **conjecture**.

## 1. Setting
Let L ≈ 0.3181 + 1.3372i be the fixed point e^L = L, ρ := |L|, D := {|x| < ρ}. We look for polynomials p_N of degree N with T_N[p_N∘p_N − exp] = 0
(T_N = truncation after degree N) that converge to Kneser's half-iterate f_K (f_K∘f_K = exp, singular at L, L̄). Newton's method needs the finite sections
J_N := T_N J|_{P_N} of the linearised operator **J δ := δ∘f + (f'∘f)·δ**, as matrices in ℓ²_ρ (coefficient k weighted by ρ^k). Test model: a real cubic f
with the same fixed-point data (f(L) = L, f'(L) = μ := √L); all numbers below refer to this model, Kneser agrees to 2–3 digits.
Geometry: Γin := {|x| = ρ, |f(x)| < ρ} (the arc from L through −ρ to L̄), Γout the rest; h := the branch of f⁻¹ with h(L) = L;
"A-arc" := the arc of the circle |x − f(0)| = |L − f(0)| outside D (it runs from L to L̄ past Γout).

## 2. Observation (numerical, exact arithmetic up to N = 1000)
J_N⁻¹ has r(N) ≈ N/250 "corner modes" with σ_i ~ e^{cN}, **c ≈ 0.0120** (cubic and Kneser alike). The left singular vectors are concentrated at arg L
(width ~1/N, height ~√N), the right (adjoint) ones on the "bar" (index k ≈ 0.95N, angle near π). Write the adjoint modes y_i := σ_iv_i as Laurent
polynomials Y_i(t) = Σ_{k≤N} y_{i,k}t^{−k}, t = x/ρ.

## 3. Hypothesis (H₀)
There is a fixed closed curve γ* through L, L̄ (on the Γout side inside D through h(A-arc), on the Γin side just outside D) such that
**|Y_i(t)| ≤ C·N^κ for t ∈ γ*, all corner modes i, all N ≥ N₀.**
**Proved (Theorem W1'')**: (H₀) implies, for the corner components of the pole responses T_i(ζ) = ⟨u_i, J_N⁻¹T_NJ(k_ζ − T_Nk_ζ)⟩, the bound
|T_i(ζ)| ≤ C₁N^{κ'}e^{−(N+2)g_K(ζ)} outside a compact set K (g_K its Green function), also on the slit near L. Since the singularity of f_K sits at L and
its data decay only polynomially, the convergence proof needs the **polynomial** form; e^{o(N)} instead of C·N^κ is not enough.
**Measured up to N = 1000**: C = e^{−4.06}, κ = 0.88 on γ*; corner constants C_E = 1, κ_E = 0.6; reserve ε ≥ 0.60.

## 4. Exact starting identity (proved)
With F̃(t) := f(ρt)/ρ, m̃ := (f'∘f)(ρ·), 𝒰(t) := Σ u_k t^{−k−1}, P a polynomial:
**T(t) := m̃(t)Y(t)/t + Σ_{F̃(z)=t} Y(z)/(z F̃'(z)) − 𝒰(t) + P(t) = O(t^{−N−2}) at ∞**,
T rational with poles only at 0 and F̃(0) — equivalent to J_Nᵀy = u. The sum runs over the three preimages.

## 5. Limit problem (numerical, no proof)
(1/N)ln|Y| and (1/N)ln|T| converge to a function H that is harmonic on the **glued surface ℳ**: cut out of the plane the region between the "bay" curve
B (inside D, Γout side) and the zero curve Z_T of T (outside D), and glue B to Z_T by f (H = H∘h). ℳ has genus 0; at 0 and ∞, H ~ −ln|t|; on the seam
H = 0; B is a free boundary. At L, ℳ is a cone with deck translation 2πi − ln μ: H ≈ Re(c₁ζ^{κ₁}), **κ₁ = 2πi/(2πi − ln μ) = 1.118 − 0.032i**
(data: 1.114 − 0.018i). The free boundary problem, with no free parameters, gives **c = 0.0118 ± 0.0004** (measured 0.0117–0.0120), the boundary
profile beyond 8° to within 3e-4, and ∂c/∂μ = 0.044 − 0.040i (measured 0.042 − 0.035i).
Inside, H is the Green function g_{Ω_y}(·, 0) of the domain Ω_y bounded by B and the zero line of a thin "lens" outside Γin, and the N zeros of Y_N
distribute according to the harmonic measure of Ω_y (N = 1000: fraction on B 0.470 vs 0.473; angular KS distance 0.005; distances ~1.4 ln N/N).

## 6. What does not work (findings)
- **Local arguments:** via the identity, the bay bound is EQUIVALENT to a bound for T on ∂K (the circularity is intrinsic). The corner zone around L
  needs angle-resolved boundary data of polynomial accuracy at its outer edge; the inner disk |t − L| ≤ N^{−1/Re κ₁} then follows by the maximum principle.
- **Extremality in ℓ²_ρ:** ‖J_N‖ = e^{0.47N}; a test vector built purely from the model has overlap 0.997 with v₁ (N = 800), but Rayleigh quotient 10^{90}.
- **Gonchar–Rakhmanov in the plane:** the exact interpolants Q_n (Q_n∘f + (f'∘f)Q_n = O(x^n)) satisfy ∮Q_n(s)ψ_k(s)ds = 0 with
  ψ_k = R_k + (f'∘f)(s)s^{−k−1}, R_k(s) = Σ_{f(x)=s}x^{−k−1}/f'(x). Each condition couples a moment at f(0), transported by f, with a moment at 0 of the
  same index — not a system {s^jW(s)}. On ℳ the problem would have GR form, but ℳ exists only in the limit (the gluing T ≈ Y∘h is a dominant balance
  with relative errors e^{−Nδ}, δ > 0 except at L, L̄).
- **Other scaling circles:** in ℓ²_{ρ'} with ρ' ≠ ρ there is no ρ' without exponential growth.

## 7. Questions for experts
1. **n-th root level:** limits of the rates satisfy (proved) a "tropical balance": in the exact identity the maximum of the term rates is always
   attained twice. Is the limit configuration (with normalisation and the singularities at 0, ∞) unique, and is c its principal value (no configuration
   with a larger rate)? Structurally: an obstacle problem (contact set = bay) with the non-local coupling H = H∘h.
2. **Strong asymptotics:** is there a Deift–Zhou-type framework for a "welding RH problem" (two domains glued along a curve by a piecewise analytic
   boundary map) with a small NON-LOCAL perturbation, with a local parametrix at the cone point L (inner layer of width 1/N, prefactor N^{1/2})?
3. **Literature:** finite sections of weighted composition operators with a fixed point on the truncation circle; polynomial solutions of approximate
   functional equations with composition; relation to Muttalib–Borodin biorthogonality (composition x ↦ x^θ) — which, however, lacks the summand
   (f'∘f)·δ that creates the corner at L here.

Files and evidence: w_grenzschicht/BERICHT.md §15–§26 (scripts w38–w53), in German.
