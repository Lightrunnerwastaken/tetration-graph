# Open problem 2: a polynomial square root of exp → Kneser — conditional theorem on K_{r₀} and the remaining coherence statement

Status 12 Oct 2026, strand K (phase 3), numbers from K and V. English version of OFFENES_PROBLEM_2.md, same content.
Labels: **proved** (V: checked by a second strand) / **numerical** / **conjecture**.
The older problem sheet (w_grenzschicht/OPEN_PROBLEM_en.md, hypothesis (H₀)) remains valid. This sheet records the state after phase 3.

## 1. Setting
- L ≈ 0.31813 + 1.33724i is the fixed point e^L = L, ρ := |L|, D := {|x| < ρ}, μ := √L.
- f_K is Kneser's half-iterate (f_K∘f_K = exp, singular at L, L̄). p_N are polynomials of degree N with T_N[p_N∘p_N − exp] = 0 (T_N truncates after degree N),
  connected by a σ-homotopy from p₁. The branch up to N = 400 is rigorously enclosed; sup_{D̄}|p₄₀₀ − f_K| ≤ 1.75e-12 **(proved, V)**.
- Linearisation: **J δ := δ∘f + m·δ, m := f′∘f** (composition + multiplication); sections J_N := T_NJ|_{P_N}, rescaled J̃ = DJD⁻¹, D = diag(ρ^k).
- Test model: a real cubic f with f(L) = L, f′(L) = μ. Locally at L it agrees with Kneser to 4e-4; the numbers below agree to 2–3 digits.
- Heads (rank-one bordering): c_n = J̃_n⁻¹e_n, ℓ_n = J̃_n^{−T}e_n, δ_n(z) = Σc_{n,k}(z/ρ)^k, Λ_n(t) = Σℓ_{n,k}(ρ/t)^k.
- Kernel K_n(z,t) = ε_zᵀJ̃_n⁻¹ẽ_t = t·(J_n⁻¹T_n[1/(t − ·)])(z), with ΔK_n = δ_n(z)Λ_n(t)/ℓ_nn **(proved, V)**.
- Corner pole response u_n := J_n⁻¹T_n[1/(L − ·)], so ΔK_n(z,L) = L·Δu_n(z). Corner scale d_n = ρn^{−a}, a = 1 + iL/(4π).
- Geometry: Γin = {|x| = ρ, |f(x)| < ρ} (the arc through −ρ, the "bar" at 180°), Γout the rest. K_{r₀} := D̄ ∖ (U_{r₀}(L) ∪ U_{r₀}(L̄)).
- Contour 𝒞: logarithmic spirals Σ₁ (outside) and Σ₂ (inside) at L (pitch 1.623° per e-fold), a lens |t| = ρ(1 + S) over Γin, the Γout arc, and the conjugates at L̄.
  E := closure of int𝒞 ∪ f_K(int𝒞) is a quasidisk with a spiral corner (opening ≈ 132°) **(numerical)**.

## 2. Conditional main theorem (K_{r₀} version)
**Theorem A (proved up to the hypotheses; logic checked by V).** Under (EP), (FA), (S_E): sup_{K_{r₀}}|Π_Nf_K − f_K| ≤ C·N^{γ − s_E}, γ = max(0, 1 + κ_G + κ_A).
Here Π_N is the Galerkin projection. Measured: κ_G ≈ −3.2, κ_A = 0, hence γ = 0 and the rate is N^{−3.96}.
**Main theorem (conditional).** If these hypotheses hold uniformly for the family of Newton neighbours (F1–F4), then for N ≥ 400 the polynomial p_N is the unique solution in a
ball around f_K, and sup_{K_{r₀}}|p_N − f_K| → 0 polynomially. Consequently a_k^{(N)} → a_k for every k, and p_N(1) → sexp_K(½).
**Uniform convergence on D̄ is not claimed.** In a model with a fast spiral (same structure), Π_NU − U stops converging in a fixed Γout zone near L from N ≈ 1000 on
(plateau up to N = 2800; V, test T1). Reason: the primal head grows there in a sector that rotates by 1.45° per e-fold in n. Hence r₀ > r_* (Kneser: r_* ~ 1e-6, estimate).

## 3. The open statements (precise)
- **(G-P)** G_n := sup_{z∈K_{r₀}}|Δu_n(z)| ≤ C·n^{κ_G}. Numerically κ_G ≈ −3.2 up to n = 1800 (cubic = Kneser to 1 %).
- **(A-lok)** sup_{t∈Σ₁∪Σ₂, |t−L|≤r₁}|Λ_n(t)/Λ_n(L)| ≤ C·n^{κ_A}. Numerically 1.008 up to n = 2000.
- **(FA)** lim sup_n max_{𝒞∖U_{r₁}(L,L̄)} (1/n)ln|Λ_n/Λ_n(L)| < 0. Numerically: lens (S = 0.03) −0.011 … −0.013, Γout arc −0.037 … −0.053.
  It follows from a limit problem of strand W (W12′).
- **(S_E)** dist_E(f_K, P_N) ≤ C·N^{−5.417α_E}, α_E = exterior angle/π. **Proved** via Andrievskii–Pritsker–Varga, Theorem 2, up to the numerically checked geometry of E.
- **Family (F3′/(BU))**, shown with a bootstrap skeleton (logic proved):
  - Along the homotopy, suppose CF(P) := min_g sup_γ|P − f_K + x^{N+2}g| ≤ δ₀ on the lens contour γ.
  - Then CF ≤ δ₀/2 for all N, provided **(BU)** holds: ‖T_{N+2}(Σ_{n≤N}[w^m](δ_n∘k)·π_n/ℓ_nn)‖ ≤ C·N^{κ_B − s_E}.
  - Here k: 𝔻 → int γ is conformal and π_n = ⟨Λ_n, data⟩ are the head pairings.
  - Termwise (using only |π_n|), (BU) closes at best with margin 0, and only on one critical lens. **(BU) with a reserve is a statement about the phases of the π_n.**
- **Common core:** at the bar, primal head × adjoint are exactly balanced at the n-th-root level (the sum of the rates tends to 0).
  - Every absolute-value estimate therefore ends at places with margin 0. What decides is the coherence of the data and the corner growth zone near L outside D.
  - For the family: the best representative is more coherent than a single head by 0.0024 per degree (V).

## 4. Exact tools
1. **(CD1) Christoffel–Darboux identity (proved, V; cubic and Kneser to 1e-176):**
   (z − t)·K_n(z,t) = q_{n+1}(z)·λ_n(t) − t·ν_n(z) + Ξ_n(z,t), Ξ_n(z,t) = (1/2πi)∮K_n(z,s)·(f(s) − s)·K_n(f(s),t)ds/s.
   - q_{n+1} is the monic primal head of order n+1, λ_n the adjoint head, ν_n = J_n⁻¹·1.
   - Derivation: J S − S J = T_nM_{f−x}C_f − (T_nf^{n+1})e_nᵀ (first-order displacement of the composition channel).
   - Without composition this is the Toeplitz CD formula. So: kernel = rank-2 integrable kernel + a Dyson term carrying the factor f(s) − s.
2. **Explicit corner profiles (proved, V):** heads and corner pole responses near L follow F_β(u) = Γ(β)⁻¹∫₀^∞t^{β−1}e^{−t−ut^{−iα}}dt, with α = L/(4π), β_P = ½ + iα, β_A = ½ − 2iα.
   Corner product P_n = ΔK_n(L,L) = aCL/(n d_n)·(1 + O(1/n)), with C explicit (derivation sketched, numerically 1e-5).
3. **K19 (proved, V):** sup_{K_{r₀}}|Π_Nf − f| ≤ (1 + c(1+M)Σ_{n≤N}Q_n)·dist_E(f, P_N). The data side is classical approximation on E; there are no truncated poles.
4. **Rank-one bordering and the LU decomposition of the heads (proved, V).** Lemma K20: the primal n-th-root rate relative to the corner is ≡ 0 on D̄
   (proved, V; pointwise up to polar sets).
5. **Bar formula (numerical):** ΔK_n(t,t)·m(t) → κ(t) (the local radial exponent of the head, 0.87 … 0.96).
   Moreover K_n(t,t)·m(t) = t·(log δ_n)′(t) + Δ(t) + o(1) with Δ independent of n. The circle mean of ΔK_n(t,t) is exactly tr J̃_n⁻¹ − tr J̃_{n−1}⁻¹ (→ 0.6215).

## 5. Measured values (cubic unless stated otherwise)
| Quantity | Value |
|---|---|
| ‖J̃_N‖, ‖J̃_N⁻¹‖ | e^{0.47N} (composition overflow), e^{0.0118N} (corner mode) |
| \|P_n\| | 0.535·n^{−0.106} from both sides, up to n = 2000, Kneser the same |
| R_P(r₀ = 0.05) = sup_{K_{r₀}}\|δ_n/δ_n(L)\| | ~ n^{−3.1} (1000 … 1800) |
| primal on Σ₁ at r | rate ≈ 0.15·r per degree (growth in the corner zone outside D) |
| bar (lens S): primal / adjoint / sum | ≈ 0.957S − 0.012 / ≈ −(0.957S − 0.012) / +0.0002, decreasing (n ≤ 600) |
| truncation remainder (I − T_n)Ju_n on Γout | rate ≈ ½·ln(\|f\|/ρ), up to +0.252 per degree |
| p_N − f_K on the lens S = 0.03 / at the bar / on Γout, Σ₂ | 1.75 at N = 2000 (e^{0.015N}) / 1e-4 / converges (V) |
| CF representative on the lens S | rate 0.90·(S − 0.0147); CF + adjoint bar rate ≈ −0.0024, independent of S (V, N ≤ 2000) |
| s_E | measured 3.9–4.1, Theorem K22: 3.96 |
| Dyson fixed point from (CD1), linearised | spectral radius > 1 from n ≈ 58 on, ≈ e^{0.004n} (norm-independent) |
| Ξ_n / generator on the corner scale | 0.3 … 3.7, slowly increasing |

## 6. Failed approaches (and why)
- Parametrix from the bulk inverse, perturbation of solved models, explicit parametrices (Taylor, approximant on E):
  loss of norm through composition transport; primal quantities at f(Γout) grow.
- Absolute estimates on a contour, also with mixed data representations: the adjoint has a reserve only outside D, and exactly there the primal grows.
  At the bar the two balance exactly (margin 0).
- Recursions in n (Volterra, generator Schur), absolute near-orthogonality, positivity: cancellations, corner coherence, logarithmic phase rotation.
- IIKS from displacement rank 2 (G40): for the inverse the composition term is not of small rank (the Dyson term in CD1).
- Wronskian pairing on the Abel torus at the corner: the primal and adjoint heads are homogeneous in complementary sectors.
- (FA) via a superharmonic comparison function: impossible (Laplace mass of the heads).
- Identity via the truncation remainder: tautological and exponential.
- Dyson contraction from CD1: the spectral radius exceeds 1. Taking both channels together only gives X = J⁻¹.
- Locally integrable corner kernel: on the corner scale the Dyson term is homogeneous of degree 0. Every local parametrix needs both channels.

## 7. Questions for experts
1. **Riemann–Hilbert with composition.** CD1 = rank-2 integrable kernel + the Dyson term ∮K(z,s)(f(s) − s)K(f(s),t)ds/s.
   - The head conditions [x^k](δ∘f + mδ) = 0 (k < n) form a two-point Hermite–Padé problem (points 0 and f(0), coupled through f).
   - Is there an RH framework for such kernels, for instance 3×3, or with a Carleman shift or conformal welding? It would have to contain a corner where both channels are of
     the same order: a spiral cone point with Mellin exponent a = 1 + iL/(4π).
2. **Coherent Carathéodory–Fejér.** Are there theorems saying that the germs of Galerkin- or Padé-type errors have representatives which, on a domain beyond the disc of
   convergence, are better than the moduli of the individual contributions allow? Concretely (BU): the phases of the pairings π_n give a reserve of 0.0024 per degree.
3. **Christoffel functions for composition operators.** Outside D, ΔK_n(t,t)·m(t) converges to the local exponent of the heads
   (as a Christoffel function converges to an equilibrium density). Is there a potential theory (a vector equilibrium problem for the zeros of the primal and adjoint heads)
   that explains the exact cancellation h_P + h_A = 0 at the bar?
4. **Sections of weighted composition operators** u ↦ u∘f + m·u with a repelling boundary fixed point on the truncation circle.
   - Known results on uniform invertibility on subsets away from the fixed point (instead of on the whole disc)?
   - Known results on the failure of uniform convergence up to the fixed point?

Evidence (in German): phase3/HAUPTSATZ_KR0.md (version 3), phase3/KERNPROBLEM.md §4, phase3/k_kern/BERICHT.md §26–§43, phase3/v_pruefstand/BERICHT.md.
