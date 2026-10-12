# Translation guide: German → English (Tetration knowledge graph)

Audience: English-speaking mathematicians and advanced amateurs of the tetration community (tetration forum, MathOverflow, googology).
Register: plain, precise mathematical English. Not literary, not chatty. Keep sentence structure close to the German where the
logic depends on it; do not smooth over hedges.

## Hard rules (a script checks these; violations are errors)

1. **Every `$…$` and `$$…$$` formula, every `\tex` string, every number, every identifier like `kn.konstr.kneser` or `exp.funde`
   stays character-for-character identical.** Do not translate, "improve", re-typeset or re-space formulas. German words *inside*
   a formula (`\text{mit}`, `\mathrm{reg}`) are translated only if they are `\text{…}`/`\quad\text{…}` prose (e.g. `\text{mit}` → `\text{with}`),
   never `\mathrm{…}` labels (K, reg, …).
2. Markdown (`**bold**`, `[[id]]`, line breaks, `(i)`, `(ii)`) keeps its structure.
3. Do not add or drop content. Do not add explanations. Do not merge or split fields.
4. Keep names and citations: Kneser, Écalle, Voronin, Baker–Rippon, Boshernitzan, Hardy, van der Hoeven, Fibonacci Quart. …
5. Decimal separator is already a point in the source. German "Tausender"-style spacing: keep as is.
6. Keys and enumerated values (`art`, `status`, `vollstaendig`, `literatur`, …) are **not** translated – only the text fields listed in the task.

## Epistemic words: must map exactly (never upgrade or downgrade)

| German | English |
|---|---|
| bewiesen | proved (adj.: proven) |
| bewiesen_literatur | proved in the literature |
| bedingt | conditional |
| numerisch | numerical / numerically supported |
| Vermutung | conjecture |
| offen | open |
| Annahme | assumption |
| Skizze | sketch |
| gescheitert | failed |
| widerlegt | refuted |
| überholt | superseded |
| Beobachtung | observation |
| Gegenbeispiel | counterexample |
| Satz / Lemma / Korollar / Hilfssatz | theorem / lemma / corollary / auxiliary lemma |
| genau dann, wenn | if and only if |
| hinreichend / notwendig | sufficient / necessary |
| vermutlich, wohl, offenbar | presumably, probably, apparently (keep the strength) |
| "nur numerisch", "nicht bewiesen" | "only numerical", "not proved" |

Preserve quantifiers and scope exactly ("jede" = every, "eine" = a/one, "genau ein" = exactly one, "für alle … " = for all …).

## Terms

| German | English |
|---|---|
| Tetration, Pentation, Hexation | tetration, pentation, hexation |
| Superfunktion | superfunction |
| Superlogarithmus (slog) | super-logarithm (slog) |
| Superwurzel | super-root |
| Abel-Funktion / Abel-Gleichung | Abel function / Abel equation |
| Schröder-Funktion, Böttcher-Funktion | Schröder function, Böttcher function |
| Koenigs-Funktion | Koenigs function |
| reguläre Iteration / reguläre Tetration | regular iteration / regular tetration |
| Kneser-Tetration, Knesers Lösung | Kneser's tetration, Kneser's solution |
| Fixpunkt (anziehend / abstossend / parabolisch / indifferent) | fixed point (attracting / repelling / parabolic / indifferent) |
| Multiplikator | multiplier |
| Keim | germ |
| Basis, Basiswechsel | base, base change |
| Basen-Ebene | base plane (parameter plane of the base) |
| Türme, Turm, unendlicher Turm | towers, tower, infinite tower |
| Strang | strand (line of research) |
| Leiter, Logarithmenleiter | ladder, logarithm ladder |
| Becken | basin |
| Julia-Menge | Julia set |
| Rechteckwelle | square wave |
| Phase, Phasengesetz | phase, phase law |
| Mode, Moden (Fourier) | mode, modes |
| Streifen | strip |
| Halbe Exponentialfunktion, Halbexponentiale | half-exponential (function) |
| Exponentialität | exponentiality |
| Transreihe, Hardy-Feld | transseries, Hardy field |
| Zeit-1-Abbildung | time-1 map (digits stay digits) |
| Hornabbildungen | horn maps |
| Blätter (parabolische Dynamik, Écalle–Voronin) | petals (Fatou petals) |
| Stufe / Stufenleiter | level / ladder of levels |
| Gebiet (komplexe Analysis) | domain |
| Bereich (Feld `bereich`) | domain of validity / scope (follow the surrounding sentence) |
| Folge | sequence (of numbers) **or** consequence (of a statement) – decide by context; "Messbare Folge" = measurable consequence |
| Naht | seam |
| Kontinuumssumme | continuum sum |
| Rückwärtsdifferenz | backward difference |
| gebrochene Exponentialität | fractional exponentiality |
| Rechenfehler | rounding / computational error (by context) |
| Engine, Rechenwerk | engine |
| Ausgleichsgerade | best-fit line |
| Gegenprobe | control test |
| Barrows Satz | Barrow's theorem |
| Gegenbeispielsuche | counterexample search |
| Prüffrage | check question |
| Voraussetzung | hypothesis (use "assumption" only for Annahme) |
| Bereich | domain (of validity) |
| Bruchsumme | fractional sum |
| Bruchprodukt | fractional product |
| Ableitung / Steigung | derivative / slope |
| Wendepunkt | inflection point |
| Grenzfall | borderline case / limiting case (by context) |
| Reihenentwicklung | series expansion |
| gleichmässig | uniformly |
| ganzzahlig | integer-valued |
| Hauptzweig | principal branch |
| holomorph | holomorphic |

Spelling: American English (neighbor, behavior, normalize, analyze – a script warns on British forms). "gleichberechtigt" = "on an equal footing". Swiss "ss" in source (gleichmässig) is just spelling.
Titles: sentence case, no trailing period, keep any formula or term in the title.
`einfach` texts: simple everyday English for a reader who has had school mathematics plus derivative/integral; keep it short and warm, but never loosen the claim.
