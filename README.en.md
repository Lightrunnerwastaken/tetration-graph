> English translation of the German README (state of 2026-10-12); the German text is authoritative.

# Tetration Knowledge Graph

Verifiable, statement-oriented knowledge graph on tetration and its connections
(functional equations, iteration theory, analysis, arithmetic, logic, dynamics).

**View:** open `ausgabe/tetration-graph.html` directly in a browser (a single file, offline).
**Audit report:** `ausgabe/pruefbericht.md` · **Map as text:** `ausgabe/landkarte.md`.

**As of 2026-10-12:** 367 statements · 637 edges · 54 experiments (all current) · 292 sources · 9 areas · 8 main ideas with 37 sub-ideas.
Checker: 0 errors, 0 warnings, 26 hints. Numbers change with every build; authoritative are the head of `ausgabe/pruefbericht.md`
and `ausgabe/graph.json` (`statistik`).

```
python baue.py                 # check + build
python baue.py --nur-pruefen   # audit report only
python baue.py --experimente   # first recompute all experiments
python -m pytest -q tests      # rules of the checker
```
Requires Python 3 with `pyyaml` for checking and building (`pip install pyyaml pytest`). The experiments additionally need
`mpmath`, `numpy`, `scipy`, `matplotlib`; some (Kneser engine, hexation) also PARI/GP. If you only want to check and build without experiments,
`pyyaml` is enough.

## Structure

| Location | Contents |
|---|---|
| `daten/aussagen/*.yaml` | one file per area (`00_`…`10_`) and supplementary files `lit_<gebiet>_*.yaml` of the literature map; each statement with objects, domain, hypotheses (each with `benutzt: ja/nein/unklar`, i.e. used: yes/no/unclear), conclusion, proof + gaps, sources, edges, status, assessment (correctness/novelty/usefulness kept separate), five check questions, versions |
| `daten/quellen.yaml`, `daten/quellen/*.yaml` | literature with audit note `original` / `sekundaer` / `nein` (checked against the original / only via a secondary source / not checked); supplements of the literature map are in `quellen/lit_L<n>.yaml`; the loader merges them, a key in two files is an error |
| `daten/experimente.yaml`, `daten/experimente/*.yaml` | experiments → script, result, image, supported statements |
| `daten/objekte.yaml`, `daten/objekte/*.yaml` | vocabulary of the objects |
| `wissen/` | schema (vocabulary), loading, checking/deriving, output |
| `experimente/e*.py` | reproducible computations with self-checks; the result JSON carries the script hash |
| `daten/straenge.yaml`, `daten/straenge/<gebiet>.yaml` | research strands per area with a one-line “what happens there” (middle zoom level of the map, reading paths, textbook) |
| `daten/ideen.yaml` | idea map: 8 main ideas, 37 sub-ideas; each statement belongs to exactly one sub-idea via its strand (or a field `idee:`) |
| `wissen/landkarte.py` | computes coupling, position and fill level of the ideas from the edges; writes `ausgabe/landkarte.md` |
| `forschung/` | working folders of the research strands (`polynom_halbexp`, `literaturkarte`, `eta_hyperoperationen`); the graph only adopts what is established from them |
| `docs/` | research reports, drafts and plans (`docs/superpowers/specs`, `plans`) |
| `ansicht/` | template, style, modules `js/01_basis … 06_steuerung` (merged into one file at build time; KaTeX embedded) |

## Idea map (start view, since 2026-10-12)

The map opens with the **mathematical structure** rather than with places: eight main ideas (The equation · Dynamics of b^z and the fixed point L ·
Infinite tower and base plane · Fractional continuation · Base change · Higher levels · Numbers and arithmetic · Bridges), whose position
is computed from the coupling (edges between their statements, shared objects). A click on a main idea shows its sub-ideas,
a click on a sub-idea its statements (main statements large, auxiliary statements small), a click on a statement the content card. Path bar at the top,
Esc goes back.

- **Fill-level ring:** shares by status (proved/conditional/numerical/conjecture/open/discarded). **Provenance ring:** literature / forum / own work only.
- **Pieces of a conjecture:** For a conjecture, a conditional theorem or a sketch, the content card shows the tree of hypotheses and marks where something is missing.
- **Text version** for research agents: `ausgabe/landkarte.md`.
- Rules: a statement without an idea → hint “ohne-idee” (no idea); unknown idea or a strand in two sub-ideas → error (tests: `tests/test_ideen.py`, `tests/test_landkarte.py`).
- Design and plan: `docs/superpowers/specs/2026-10-12-ideen-landkarte-design.md`, `docs/superpowers/plans/2026-10-12-ideen-landkarte.md`.

## Map (second view)

The earlier zoomable map (areas → progressions → statements → content) is kept as a second view. A zoomable map on a dark blueprint; what you see depends on the zoom level (bar on the left):

1. **Areas** — islands with a status ring, bands with flowing light show the cross-connections.
2. **Progressions** — research strands as regions with a name and a one-line progression; light flows along the proofs from hypothesis to conclusion.
3. **Statements** — titles, important ones first, overlapping ones stay hidden.
4. **Content** — cards with the typeset statement next to the nodes.

Click on a statement: the camera flies there, neighbors light up, second-degree neighbors half, the rest recedes; the detail
panel opens only then. Esc or a click into empty space cancels the focus. With “reduced motion”, animations are omitted.

Layers: **A** every solution of the functional equation · **B** class of solutions · **C** concrete construction · **X** auxiliary lemma from another area.

Level (for the public library): **Einstieg** (entry level: school mathematics plus derivative/integral; area “Foundations” and marked
basic statements) · **Vertiefung** (in-depth: university study) · **Forschung** (research: open questions, own work). Without an entry the rule is: Foundations = Einstieg,
base change and conjectures/observations = Forschung, otherwise Vertiefung. Entry-level statements carry `einfach:` (explanation in everyday language),
identities `formen:` (equivalent notations).

**English version (tools since 2026-10-12):** Translations are stored separately from the German texts in `daten/uebersetzung/<gebiet>.yaml`
(per statement `en` and optionally `en_theorie`, plus a hash of the translated German texts, provenance and review status). The checker enforces that
formulas, `[[…]]` references and identifiers remain unchanged, reports outdated translations (the German text has changed), deviations in numbers and
word choice, and the progress per area. Commands (`stand`, `export`, `einlesen`), glossary and pilot: `docs/uebersetzung/`.
**State 2026-10-12:** all 367 statements (core, proof, check questions, theory), areas, strands, ideas, textbook chapters, objects and experiment titles
are translated (`daten/uebersetzung_struktur/<kind>.yaml`, commands `struktur-stand|-export|-einlesen`); the README is available as `README.en.md`.
Only part of it has been proof-read so far (status: `python -m wissen.uebersetzung stand`). Open: DE/EN switch of the view, interface texts,
result summaries of the experiments, `ausgabe/landkarte.md` and `docs/`. Inconsistencies in the German text noticed along the way:
`docs/uebersetzung/DEUTSCH_ZU_KLAEREN.md`.

## What the checker enforces

- References, vocabulary, no cycles in proof chains.
- “proved” (`bewiesen`) on a weaker foundation automatically becomes “conditional” (`bedingt`); for every conditional statement the
  **root conditions** (unproved ends of the chain) are computed and must be named.
- An assumption counts as discharged as soon as a proved theorem discharges it.
- If something depends (transitively) on a refuted/superseded statement → “re-check”, until a `nachgeprueft` entry lifts it.
- Unused hypotheses → hint “abschwächbar” (can be weakened). Conjectures need a documented counterexample search,
  sketches named gaps, “presumably new” a research note.
- Experiments: a result counts only if the script hash matches; failed self-checks are reported.

## Entering a new statement

```yaml
- id: kn.satz.beispiel            # prefix = area
  titel: Short title
  art: satz                        # definition|satz|lemma|korollar|konstruktion|annahme|vermutung|beobachtung|gegenbeispiel|gescheitert|korrektur|frage
  ebene: B                         # A|B|C|X
  status: bewiesen                 # see wissen/schema.py
  objekte: [kneser_T]
  bereich: '$b>\eta$'
  voraussetzungen:
    - {text: 'T strictly increasing', benutzt: ja}
  aussage: |
    Text with $\TeX$.
  beweis: {art: vollstaendig, text: '…', luecken: []}
  stuetzt_sich_auf: [kn.konstr.kneser]
  kanten: [{typ: spezialfall_von, ziel: fg.satz.allgemeine_loesung}]
  quellen: [{ref: kneser1950, stelle: 'Theorem 1', geprueft: nein}]
  bewertung: {richtigkeit: '…', neuheit: bekannt, nutzen: '…', nutzen_stufe: mittel}
  niveau: einstieg                 # einstieg|vertiefung|forschung (optional)
  einfach: >-                      # required for entry level: in plain words
    …
  formen:                          # equivalent notations
    - {tex: '…', notiz: '…'}
  pruefung: {f1_folgt_unmittelbar: '…', f2_benutzte_voraussetzungen: '…', f3_spezialfall: '…', f4_uebertragbar: '…', f5_neue_frage: '…'}
```

## Textbook (Theory tab, since 2026-10-08)

The “Theory” tab is a textbook (the toggle “Dependencies” shows the old diagram). Chapter order: `daten/lehrbuch.yaml`
(chapter 0 “What is tetration?” as running text, then strands per chapter). Detailed explanations per statement: `daten/theorie/<gebiet>.yaml`
with the parts `worum`, `genau`, `warum`, `beispiel`, `bilder`, `einordnung` (what it is about, in detail, why, example, images, placement; Markdown with `$…$`, `[[id]]`, `![[bild:name]]`).
Figures: `python lehrbuch/bilder.py` → `lehrbuch/bilder/*.svg`. The checker reports dead references, missing parts and the
progress per area. **As of 2026-10-12: complete** — all 367 statements have a theory block (as of 2026-10-08: 218 statements in chapters 0–10, at that time per chapter:
Foundations 32, fixed points/tower 7, functional equation 11, iteration 8, below η 19, Kneser 19, base change 65, base plane 25,
numbers 16, bridges 16). New statements need a block in `daten/theorie/<gebiet>.yaml`, otherwise the check reports
“theorie-fortschritt” (theory progress). The phase picture k7_phase
recomputes slog_e(T₂(n+θ)) itself (Chebyshev interpolants from Kneser values, from n = 5 via the logarithm ladder) and matches
phi_modes.json to 4·10⁻²⁴.
Display formulas of the form “A,\qquad B” or with “,\quad” are set one beneath the other by the view when they do not fit side by side (`formelnEinpassen`).

## Comparison with wikis and the Tetration Forum (2026-10-08)

Report: `docs/recherche_wikis_forum.md` (Wikipedia en/de/ru, Googology, MathWorld, Mizugadro archive, forum: FAQ, TPID list,
1260 topic titles, ~60 topics read). **Round 1 done:** prior work from the forum entered — base-change constant μ
(Fox 2007), values m(2), m(3), m(10) and the “wobble” (Levenstein 2009), two Abel functions for √2 (Asimov–Hickerson–Schroeppel around
1991, “Bummer!” 2007), expansion at η (Trappmann–Kouznetsov 2012); counter-conjecture “base change nowhere analytic”
(Levenstein 2011) as `bw.verm.nirgends_analytisch`. The loader now reports duplicate YAML keys (previously `justus_p2` appeared twice).
**Round 2 done:** 13 core statements — super-roots (strand `gr.umkehrungen`, arithmetic), infinite super-root, TPID 6 with a
new rate $n^2(x_n-\eta)\to2\pi^2\eta/e$, upper super-exponential function and cheta, q-binomial series (Reshetnikov) with rate $q^{\operatorname{Re}z+2}$,
conjecture “T′ completely monotone ⇔ Koenigs coefficients negative”, Möbius iteration, Szekeres' monotonicity, Shishikura/Lavaurs;
experiments e28–e31. **Round 3 done:** strand `kn.hoehere_stufen` — pentation as the regular superfunction of T_e at the fixed point
−1.8503545290… (λ = 6.4607; Fox 2007 confirmed to 22 digits), η₄ = 1.63532449671527639934… and e₄ = 3.08853227180671765448… (Levenstein 2010
confirmed to 32 digits) with three characterizations; experiments e32–e33, helper module `experimente/stufen.py`. **Round 4 done:** open
conjectures with sample checks (e34–e37) — TPID 7 holds for Kneser only up to b₇ = 10.2210634158 (counterexample above it); negative Koenigs coefficients
would answer TPID 7 and Nixon's half-plane question at the same time; η₅ = 1.7373559252, e₅ = 3.3014072785 (TPID 10); Carleman eigenvalues → λ^k;
signs of the derivatives; Shell–Thron area 8.0509332809 (Robbins 2009 confirmed); history in chapter 0. With this, all four rounds of the report are done.

## Hexation (2026-10-09)

Regular hexation below η₅ at the attracting fixed point of pentation (`experimente/hexation.py`, `e38_hexation.py`, `exp.hexation` 10/10):
1.7↑↑↑↑½ = 1.3974472613, 1.7↑↑↑↑∞ = 2.3876372817; every hexation has H(−3) = −2 and a singularity in (−4, −3); the infinite hexation
converges exactly up to η₅; at an attracting complex fixed point there is no Kneser construction. Kneser-like hexation above η₅
via the engine knhexinit in the gp-tetration fork (exp-079): 2↑↑↑↑1.5 = 2.6728939593426194891907964…, and thus η₆ = 1.79912123612881268672356700703,
e₆ = 3.44223464771689323 (`e39_hexation_kneser.py`, `exp.hexation_kneser` 11/11; fourth term for TPID 10, 2 − η_n ≈ 2.8·n^−1.47).
Sheet tracking (pentation continued along paths) built and tested; found with it: the Kneser-like hexation is branched in the upper
half-plane (at −0.7408 + m + 2.7458i; logarithmic winding 2πi/ln 2, power exponent ln λ_T/ln 2 = 2.3707) — Kouznetsov's rectangle method is thereby
ruled out. Prototypes in `experimente/prototypen_hexation/LIESMICH.md` (they write their own
cache `kneser_cache_proto.json`; the shared `experimente/daten/kneser_cache.json` has grown to ~190 MB through earlier prototype runs).

## Literature map (since 2026-10-12)

Task: map the work of others as precisely as our own research (who worked on what, with which method and which
result, where several approached the same question differently). Before: 321 statements, of which 158 had only own sources. Three strands
(L1 iteration/functional equations/fixed points, L2 Kneser/base change, L3 foundations/arithmetic/bridges) put their supplements into
separate files (`daten/aussagen/lit_*.yaml`, `daten/theorie/lit_*.yaml`, `daten/quellen/lit_L*.yaml`, experiment number ranges e60–e89 reserved);
block 3 (freely available full texts: Euler E489, Cho–Kim, Kouznetsov–Trappmann, TPID 16/17) is completed, as is the candidate list.
Maps, reports and the overall index are in `forschung/literaturkarte/` (`GESAMTVERZEICHNIS.md`, `L<n>_*/KARTE_<gebiet>.md`).
Result: 46 new statements (367), 60 instead of 35 statements “bewiesen_literatur” (proved in the literature), **103 citations not yet checked against the original** (previously 68).

## Polynomial root of exp and higher hyperoperations (2026-10-10 to 2026-10-12)

- **Polynomial root** (strand `kn.halbe_exp`, working folder `forschung/polynom_halbexp/`): phases 2–3 have been adopted into the graph
  (experiment `exp.polynom_phase3`). The main theorem stands as a **sketch** with named gaps (assumption (H₀) on the corner modes at L is open;
  tenable version `kn.satz.polynom_hauptsatz_kr0`; Kneser's sexp to base e is rigorously certified, `kn.satz.kneser_e_zertifiziert`).
  Phase 3 was ended at Janis' request on 2026-10-12; resumption only on commission (entry point: `forschung/polynom_halbexp/README.md`).
- **η of higher hyperoperations** (`forschung/eta_hyperoperationen/`, standalone since 2026-10-12): η table for higher ranks with
  accuracy statements, explicitly without a general Kneser theorem for all ranks; report and mathematics in `BERICHT.md` / `MATHEMATIK.md`.

## Recovered investigations (July/August 2026, incorporated 2026-10-08)

Sources: `Downloads/Base_Change_Shared_Fixed_Point.zip` (2026-07-30), `Imaginary_Direction.zip`, `Periode_2_Region.zip` (2026-08-07).
Key numbers independently recomputed in `exp.funde` (14/14). The reconstruction `kl.satz.basiswechsel_gleiches_L` agrees with the original.

- **Shared fixed point:** A_b = σ·A_d + Ψ(κ_d) exactly; sample pair −19/25 ↔ d = 1.425186736 − 0.046203369i, σ = −0.200461 + 0.003486i.
  Attracting partners only in the wedge Re b < 0 (principal branch: only m = ±1); 2, e, 3, 4 have none (`fp.beob.partner_keil`).
- **The wave is damped:** pseudo-period 2π/|arg λ|, damping |λ|^{2π/|arg λ|}; the classical 1-periodic wave is the limiting case |λ| → 1 (`kl.satz.invariante_daempfung`, interpretation).
- **Imaginary direction:** A_b(T_d(iy)) → A_b(L_d), e.g. A_e(L_2) = 0.326448651629 + 0.523513968910i (`bw.satz.imaginaere_richtung`).
- **Period-2 region:** closed form μ = (s/(2 sinh(s/2)))², no circle (covering instead of bijection), b_left = −0.08052815762548356,
  tip √6/e·n^{−1/2}, below e^{−e} the phase becomes a square wave with |c_k| = |z_A − z_B|/(πk) (`fp.periode2`).

## Open points addressed (2026-10-08)

- **Lambert → digamma** (`kl.satz.moebius_qdigamma`, `exp.lambert_zeta` 8/8): In the Möbius model the Lambert fractional sum is exactly
  a q-digamma function with q = λ; λ → 1 gives ψ. The true fractional sums S_b(1/2) tend to S_η(1/2), error O(η − b).
- **Continuity at η** (`kl.beob.stetig_eta`): T_b(1/2) → T_η(1/2) from below (regular) and from above (Kneser), both linear in |b − η|
  and with the same slope 0.438579 / 0.438581 — numerically C¹ across η.
- **Fractional sum of the slog** (`kn.satz.abel_bruchsumme`, `exp.abel_bruchsumme` 9/9): converges because A″ = O(1/y²); S(1/2) = −0.1150976090828,
  γ_A = 0.4951728249. The relation to μ remains open.

## July synthesis incorporated (2026-10-08)

Source: `Downloads/tetration_shell_thron_fractional_sums_base_classes.md` (own synthesis, July 2026). 13 new statements, `exp.koenigs_kalkuel` (12/12), `exp.eta_digamma` (12/12).

- Koenigs spectrum below η: differentiating, integrating, shift sums as factors per mode; fractional sum as a Lambert series (`kl.satz.koenigs_spektral`).
- Value-side base change T_b∘A_d = power χ_d^ρ in Koenigs coordinates, branching (a_d − y)^ρ, Jacobi formula via fractional sums;
  quasiconformal base change and the obstruction: holomorphic + single-valued + time-preserving only for equal multiplier (`kl.satz.qk_basiswechsel`).
- Base η: Fatou series with rational coefficients (1/18, 1/135, −1/972, …), C_T = 2.79824815423138766245… (integer heights only);
  fractional sum S_η(x) = e x − 2e log x − 2eγ + Σ[T(k) − e + 2e/k] + O(log x/x), Σ[T(k) − e + 2e/k] = 10.7975622749393143…
- Reciprocal tower sums C_b = Σ 1/(b↑↑k) (C_2 = 0.8125152587890625…), normal form e^w + κ, atlas of the base classes.

## Reconstruction “linear base change at the same L” (2026-10-08)

The original investigation (claude.ai chat, early August 2026) is not available locally; derived anew in `exp.gleiches_L` (21/21).

- A fixed point L has a ladder of bases u_k = (Log L + 2πik)/L. Real pairs: e^{π/2+2πk} share L = i
  (4.8105 diverges, 0.00898 has period 2), √2 and −√2 share L = 2 and 4, plus 8 mixed positive/negative pairs (`fp.satz.gleiches_L`).
- At the common L, A_b(T_d(x)) = C₁ + C₂·x + O(λ_d^x), C₂ = log λ_d / log λ_b: the deviation falls like |λ_d|^x, no
  wave remains. For different L, A_b(T_d(x)) becomes constant. Exactly linear only for b = d (`kl.satz.basiswechsel_gleiches_L`).

## Findings of the round “regions, fractional sums, tower integrals” (2026-10-08)

Basis: the user's images 1 and 2 and his/her note *Complex Tetration Dynamics* (2026-01-05). 25 new statements, experiments e16–e20.

- **Regions of the base plane** (`fp.beob.regionen`, `exp.parameterraum`): period 1 precisely for |t(b)| < 1 with t = −W₀(−Log b),
  no counterexample in either direction; buds p/q rotate with rotation number p/q (17/17); at rational boundary points convergence like n^{−1/q}.
- **Note confirmed:** Theorem 1 (holds even for every radius), Theorem 2, corollaries 1–2, all numbers for b = i. The disk-mean identity (5)
  is now a conditional theorem (Weyl + mean value theorem, assumption: the Siegel boundary is a Jordan curve). The conjecture “circle mean = 1 also for
  non-integer heights” is **refuted** for the regular family (branching w^t).
- **Below η:** half-step via the limit (scaling by λ^x at the fixed point), fractional sums following Müller–Schleicher, and
  the derivative formula T′(x+h)/T′(x) = (ln b)^h ∏ T(x+k) with fractional product holds **precisely** for regular tetration (`kl.satz.bruchprodukt`).
- **Linear base change** (reconstruction, `kl.satz.koenigs_basiswechsel`): two bases with an attracting fixed point are conjugate via their
  Koenigs functions; the height becomes exactly k·x with k = ln λ_d / ln λ_b — without a periodic part.
- **Tower integrals:** lim I₂ₙ = 0.67671427377181279519…, lim I₂ₙ₊₁ = 0.63130937080198819539…, gap 0.0454049029698245998…;
  image 2 has the gap right, but both limits too small by 2.428·10⁻¹¹. Convergence like n^{−3/2} (predicted, measured 1.502).

## Findings of the foundations round (2026-10-08)

- **Inflection points:** Every measured Kneser tetration has exactly one inflection point (b = 1.5: 3.987, 2: −0.030, e: −0.5218, 3: −0.604,
  5: −0.847, 10: −0.997, 100: −1.182); regular ones below η have none. Half proved: convexity is inherited upward (`gr.verm.ein_wendepunkt`).
- **Numbers:** T′ₑ(0) = 1.0917673512583…, ∫₀¹Tₑ = 1.7145219588906204881…, ∫₁ᵉ slogₑ = 1.0037598695684…; slog′ₑ(0) = 1/T′ₑ(0) = 0.91594605649953…
  is at the same time the orbit value −V₁ from Paper V.
- **Own errors found and fixed:** the first formula for T″(x+1) had one factor T(x+1) too many; a coarse PSLQ search produced
  a spurious relation for ∫₀¹Tₑ (both documented as a version and a check, respectively).

## Findings from the build-up (2026-10-08)

- **New conjecture, strongly supported numerically** (`kl.verm.streifenbreite`, `kl.verm.faktor`): For $1<b<\eta$ the
  Fourier modes of the difference of the regular Abel functions at the two fixed points decay like $C(b)\,k^{\approx-5/4}e^{-2\pi^2k/|\ln\lambda_a|}$;
  rate confirmed to $\le2.6\cdot10^{-4}$ (5 bases), remaining factor equal across the bases to $4\cdot10^{-4}$. Model case for the strip/Gevrey question of the base-change phases.
- **Structure:** 8 of the 13 conditional theorems (Papers I, II, IV) hang on a single root, Hyp. 2.2 (real logarithm channel);
  discharging it via the stopped chain (Paper III) is only sketched (`bw.satz.phasengesetz_gestoppt`).
- **Erratum Paper IV:** $\pi\sqrt{2\eta/e}=3.238930\ldots$, not 3.23907 (`bw.korr.konstante`).
- **Data quality:** `mu_hub.json` carries $m(2)$ only to ~31 of the 40 displayed digits (`bw.beob.hub_genauigkeit`).
- **Weakened hypotheses:** Paper I Prop. 7.1 needs only differentiability at the point −1; A1 theorems 5, 6 and Cor. 1 hold
  for whole solution classes, not just for Kneser; Paper II Thm. B is a special case of Young's area identity.
- **Carleman slog:** converges to Kneser, but only algebraically ($\approx n^{-3.8}$, 8 digits at $n=128$).
- **Own error found:** the first version of `fp.satz.lambert` had the sign of $\operatorname{Im}L$ wrong (`exp.fixpunkte`).

## License

Code (`baue.py`, Python scripts, JavaScript/CSS of the view) under the MIT License ([`LICENSE`](LICENSE)); statements, texts,
reports, data, figures and the generated view under CC BY 4.0. Details and exceptions (KaTeX, quotations from third-party
sources) in [`LICENSE-INHALT.md`](LICENSE-INHALT.md).
