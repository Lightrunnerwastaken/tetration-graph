# Ist die Rang-5-Basiswechsel-Welle nirgends differenzierbar? — Bericht 10, Runde 25.10.2026

## Kurz in einfacher Sprache

Die Vermutung aus der letzten Runde war: Die Rang-5-Welle Φ₅ ist wie eine Weierstrass-Funktion, also überall rau und nirgends differenzierbar.
Die Rechnung sagt etwas anderes, und zwar etwas Merkwürdigeres. Die Welle ist auf jeder Skala fast glatt: Ihr Hölder-Exponent ist überall 1, also so gut wie bei einer glatten Funktion.
Trotzdem ist die übersetzte Höhe Ĝ = id + Φ₅ eine **singuläre Funktion**, so wie die Cantor-Funktion: Sie steigt streng, aber ihre Ableitung ist fast überall 0, also Φ₅′ = −1 fast überall.
Einen endlichen positiven Ableitungswert hat sie an keiner Stelle. Grund: Auf jeder Turm-Ebene multipliziert die Rang-4-Welle die Steigung mit einem Faktor 1 + P₄′. Dieser Faktor ist im Mittel 1, sein Logarithmus aber im Mittel negativ.
Das wird erst nach etwa 5 Millionen Turm-Ebenen sichtbar, also auf keiner Skala, die man je ausrechnen kann. Es folgt trotzdem aus einem Hellinger-/Kakutani-Argument, wie bei den Riesz-Produkten.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| 1 | **Wackel-Lemma (Produktformel über die Skalen):** für jedes y und h = s/U_n′(y): (H(y+h)−H(y))/h = Π_n(y)·(1 + W(U_n,s))·(1+O(s/U_n)+O(1/T_d′(U_n))), W(u,s) = (P₄(u+s)−P₄(u))/s, Π_n unabhängig von s | **hergeleitet mit Lücke** (Lücke L2: Verzerrungsschranken nur skizziert; sonst bedingt wie bw.satz.phasengesetz_rang5 und R18-B7) und **numerisch** | (e,2), Ebene 3, θ=0.02/0.05/0.1: Rest nach Vorhersage 1.1·10⁻⁶ / 1.3·10⁻⁹ / 4.4·10⁻²⁵ ≈ (0.4…0.9)/U₃. Gegenprobe (e,3), Ebene 2, θ=0.02/0.04: 3.2·10⁻²¹ / 1.9·10⁻³³ ≈ 0.1/U₂. Wackel-Amplitude selbst 8·10⁻⁴ bzw. 1.3·10⁻⁴ | neu |
| 2 | **Hölder-Exponent überall 1:** c·h·(1−δ₋)^{n(h)} ≤ H(y+h)−H(y) ≤ C·h·(1+δ₊)^{n(h)}, n(h) = slog(1/h)+O(1), δ± = max(±P₄′); also H bi-Hölder zu jedem Exponenten 1±ε, Stetigkeitsmodul h^{1−o(1)} (genauer h·e^{O(slog 1/h)}) | **hergeleitet mit Lücke** (wie #1) | (e,2): P₄′ ∈ [−8.874·10⁻⁴, 9.196·10⁻⁴] | neu; korrigiert R18-B11 („Weierstrass-artig“): kein Exponent < 1 |
| 3 | **Satz A: Nirgends ein endlicher positiver Ableitungswert.** An jeder Stelle y ∈ (L_d,∞) gilt: H′(y) ist 0, +∞ oder existiert nicht. Grund: Auf jeder Ebene schwankt D(y,·)/Π_n um mindestens c₀ = min_u max_{s∈[1,2]} \|W(u,s)\| > 0 | **hergeleitet mit Lücke** (wie #1; braucht nur: P₄ stetig, nicht konstant) | c₀(e,2) = 9.815·10⁻⁵ = 0.342·Spannweite; c₀/Spannweite = 0.342 für alle 9 Paare (P₄ fast reine Sinuswelle) | neu |
| 4 | **Satz B: H ist singulär** (streng wachsend, stetig, H′ = 0 Lebesgue-fast-überall). Darum ist Ĝ = spog_b∘H∘P_d singulär und **Φ₅′ = −1 fast überall**. Beweis: Martingal Z_n = dH/dλ auf den Ebenen-Zellen, E_λ√Z_n ≤ ∏(A_K + δ_j) → 0 | **hergeleitet mit Lücke** (L1 geerbt: Bedingungen von bw.satz.phasengesetz_rang5 und R18-B7; L2: Verzerrungsschranken δ_j → 0 nicht ausgeschrieben) | Hellinger-Abstand je Ebene 1−A_K (e,2): 1.0·10⁻⁸ (K=2), 4.8·10⁻⁸ (K=8), → 5.0945·10⁻⁸ = Var(P₄′)/8; M=64 und M=128 gleich auf alle Stellen | neu; **widerlegt R18-B11 im Wortlaut** (Φ₅ ist fast überall differenzierbar), ersetzt durch „singulär“ |
| 5 | Folgerung aus #4: dH-fast-überall gilt H′ = +∞; das Bildmass dH hat überall die lokale Dimension 1 (singulär, aber volldimensional); H bildet eine Lebesgue-Nullmenge auf ein Intervall voller Länge ab | **hergeleitet** aus #2 und #4 (Standard-Masstheorie) | — | neu |
| 6 | Zerfallsrate: für λ-fast alle y gilt (1/n)·log Π_n(y) → −χ, χ = −∫₀¹log(1+P₄′) = 2.0377869·10⁻⁷ für (e,2), ≈ Var(P₄′)/2; unter dH → +KL = 2.0378·10⁻⁷. Sichtbar erst ab n* = 1/χ ≈ 4.9·10⁶ Ebenen (Skala h ~ 1/(2↑↑(5·10⁶))) | χ: **numerisch** (Gitter 4000/8000: rel. 10⁻¹³; M=64/128 gleich); Grenzwert: **hergeleitet mit Lücke** (starkes Gesetz der grossen Zahlen für die fast unabhängigen Ziffern {U_j}) | Tabelle unten, 9 Paare: n* von 1.0·10⁶ (e,10) bis 1.5·10¹⁰ (e,2.7) | neu |
| 7 | Dichte Mengen (konstruktiv, Schachtelung der Ebenen): Stellen mit H′ = 0, mit H′ = +∞ und ohne Ableitung liegen je dicht | **hergeleitet mit Lücke** (wie #1; braucht P₄ ∈ C¹) | — | neu |
| 8 | Dichotomie über die Ränge: Der Rang-r-Konjugator ist glatt, wenn der Schälschritt eine konstante Phase hat (Rang 4: log_b∘exp_d = α·id), und singulär, wenn dessen Phase nicht konstant ist. Also singulär für alle r ≥ 5 (Rang 6: c₆(V) = Φ₅({V}) aus R18-B15 ist selbst singulär) | **Vermutung** für r ≥ 6 (Beweis wie #4 eine Stufe höher, Lücke für gerade Ränge wie R18-B3); r = 5 siehe #4 | R18-B15: c₆ = Φ₅ ab V ≈ 3.3 | neu |
| 9 | **Negativ:** Die numerisch gemessene „Ableitung“ Ĝ′ ≈ 1.0657 (R18-B10) ist kein Ableitungswert, sondern ein Plateau Π_n. Keine ausrechenbare Skala zeigt den Unterschied zu einer Lipschitz-Funktion. Der Fourier-Schwanz \|a_k\| ~ k^{−1.3} (k ≤ 64) gehört zum Vorlauf. Sicher ist nur a_k = O(1/k), weil dĜ ein endliches Mass ist | **negativ** | — | — |

Kennzahlen der Rang-4-Welle P₄ = Φ₄ − μ₄ (`s1_konstanten.py`, `daten/s1_konstanten.json`; Moden aus Runde 18, nur gelesen):

| (b,d) | Spannweite P₄ | P₄′ min / max | χ = −∫log(1+P₄′) | Var(P₄′)/2 | n* = 1/χ |
|---|---|---|---|---|---|
| (e,2) | 2.873·10⁻⁴ | −8.874·10⁻⁴ / 9.196·10⁻⁴ | 2.0377869·10⁻⁷ | 2.0378192·10⁻⁷ | 4.9·10⁶ |
| (2,e) | 2.873·10⁻⁴ | −9.188·10⁻⁴ / 8.881·10⁻⁴ | 2.0378029·10⁻⁷ | 2.0377711·10⁻⁷ | 4.9·10⁶ |
| (e,3) | 7.292·10⁻⁵ | −2.333·10⁻⁴ / 2.253·10⁻⁴ | 1.3126468·10⁻⁸ | 1.3126415·10⁻⁸ | 7.6·10⁷ |
| (2,3) | 3.602·10⁻⁴ | −1.152·10⁻³ / 1.114·10⁻³ | 3.2033025·10⁻⁷ | 3.2032397·10⁻⁷ | 3.1·10⁶ |
| (e,10) | 6.344·10⁻⁴ | −2.029·10⁻³ / 1.961·10⁻³ | 9.9356391·10⁻⁷ | 9.9352887·10⁻⁷ | 1.0·10⁶ |
| (e,4) | 2.518·10⁻⁴ | −8.054·10⁻⁴ / 7.780·10⁻⁴ | 1.5649939·10⁻⁷ | 1.5649721·10⁻⁷ | 6.4·10⁶ |
| (e,1.8) | 4.181·10⁻⁴ | −1.291·10⁻³ / 1.339·10⁻³ | 4.3156955·10⁻⁷ | 4.3157948·10⁻⁷ | 2.3·10⁶ |
| (e,2.8) | 2.267·10⁻⁵ | −7.253·10⁻⁵ / 7.003·10⁻⁵ | 1.2684647·10⁻⁹ | 1.2684631·10⁻⁹ | 7.9·10⁸ |
| (e,2.7) | 5.263·10⁻⁶ | −1.626·10⁻⁵ / 1.684·10⁻⁵ | 6.8377352·10⁻¹¹ | 6.8377372·10⁻¹¹ | 1.5·10¹⁰ |

(Die Dateinamen `9_5`, `14_5`, `27_10` der Runde 18 stehen für die Basen 9/5 = 1.8, 14/5 = 2.8, 27/10 = 2.7.)
Für alle Paare gilt χ ≈ 4π²|a₁|² (dominante Mode, P₄′ ≈ 4π|a₁|cos): (e,2) |a₁| = 7.18·10⁻⁵ gibt 2.04·10⁻⁷.

## Beweise und Details

**Bezeichnungen** wie R18-01: F = slog_b∘T_d, H = lim slog_b^n∘T_d^n (Konjugator, T_b∘H = H∘T_d), U_n = T_d^n(y), c(U) = H(U) − U.
P₄ = Φ₄ − μ₄ ist die mittelwertfreie Rang-4-Welle. Aus R18-B7 (bedingt bewiesen): c(U) = μ₄ + P₄(U) + r(U) mit sup_{U≥V}|r(U)| → 0.
Genauer gilt r(U) = ε(U) + [slog_b(T_dU + c(T_dU)) − slog_b(T_dU)]. Dabei ist ε(U) = F(U) − U − μ₄ − P₄(U) → 0 (Rang-4-Phasengesetz), und die Klammer ist ≤ C·slog_b′(T_dU − C) (Mittelwertsatz, |c| ≤ C).

**#1 (Wackel-Lemma).** Exakt gilt H(y) = slog_b^n(H(U_n)) und H(U_n) = U_n + c(U_n). Für h = s/U_n′:
H(y+h) − H(y) = (slog_b^n)′(ξ)·[ΔU + P₄(U_n+ΔU) − P₄(U_n) + Δr].
Dabei ist ΔU = U_n(y+h) − U_n(y) = s(1 + O(s·U_n″/U_n′²)) und U_n″/U_n′² = O(1/U_n) (Lücke L2: aus T″/T′² ≈ 1/T für grosse Argumente, nicht ausgeschrieben).
Weiter ist |Δr| ≤ 2 sup|r|, und (slog_b^n)′(ξ)/(slog_b^n)′(H(U_n)) = 1 + O(s/U_n) (Lücke L2).
Mit Π_n(y) := (slog_b^n)′(H(U_n))·U_n′(y) folgt die Formel. Zwei Grenzfälle: Für s → 0 gibt sie Π_n(1+P₄′(U_n)) = Π_{n+1}, das ist R18-B9. Für s ≫ 1 gibt sie Π_n.
Numerisch (`s2_wackeln.py`) ist Q(s) = D(s)/(1+W(U_n,s)) konstant in s ∈ [10⁻³, 10³] bis auf (0.1…0.9)/U_n. Das ist genau die vorhergesagte Fehlerordnung.
Auf Ebene 2 von (e,2) (U₂ = 6…290) überdeckt der Fehler 1/U₂ die Welle (10⁻³). Dort ist der Test nicht aussagekräftig, und er wird nicht gewertet.

**#2.** Auf dem Plateau zwischen zwei Ebenen ist D ≈ Π_n. Es gilt Π_{n+1}/Π_n = (1+P₄′(U_n))(1+O(1/U_n)). Dazu kommt das Wackeln 1+W mit \|W\| ≤ max\|P₄′\|. Daraus folgen beide Schranken.

**#3 (Satz A).** Angenommen, H′(y) = L ∈ (0,∞) existiert. Dann gilt D(y, s/U_n′)/D(y, S/U_n′) → 1 für n → ∞, gleichmässig für s ∈ [1,2] und festes S.
Nach #1 ist dieser Quotient (1+W(U_n,s))/(1+W(U_n,S)) + o(1). Es gilt \|W(u,S)\| ≤ Spannweite/S.
Ausserdem gibt es zu jedem u ein s ∈ [1,2] mit \|W(u,s)\| ≥ c₀ ≥ Spannweite/4, weil u+s für s ∈ [1,2] eine volle Periode überstreicht.
Für S > 8 ist das ein Widerspruch. Gebraucht wird nur, dass P₄ stetig und nicht konstant ist (a₁(e,2) = 7.2·10⁻⁵ auf 30 Stellen).

**#4 (Satz B).** Zellen der Ebene n sind die y-Intervalle, auf denen ⌊K·U_j⌋ für j ≤ n fest ist. Sie schachteln sich, und ihre Länge geht gegen 0, also erzeugen sie die Borelmengen.
Setze μ = dH und Z_n = μ(C_n(y))/λ(C_n(y)). Unter λ ist das ein positives Martingal mit E_λ[√Z_{n+1} | F_n] = √Z_n·Σ_i √(p_i q_i).
Fasst man die Kinder einer Zelle nach ihrem Ziffernfach i = ⌊K{U_{n+1}}⌋ zusammen, wird die Affinität grösser (Vergröberung), also ist die Fassung über die Fächer eine obere Schranke.
In der U_{n+1}-Koordinate überdeckt eine Zelle sehr viele Einheitsintervalle. Die Gewichte 1/U_{n+1}′ (für λ) und (slog_b^{n+1})′(H(U)) (für μ) ändern sich auf einem Einheitsintervall nur um O(1/U_{n+1}) (L2).
Daraus folgt q_i = 1/K + δ und, mit #1 integriert, p_i = 1/K + P₄((i+1)/K) − P₄(i/K) + δ, wobei δ_j = O(1/U_{j+1}) + O(K/T_d′(U_j)) + O(sup|r|) → 0.
Also gilt E_λ√Z_n ≤ ∏_j (A_K + δ_j) → 0, denn A_K < 1, sobald P₄ nicht 1/K-periodisch ist (schon K = 2 reicht für alle Paare: 1−A₂ = 1.0·10⁻⁸ bei (e,2)).
Daraus folgt Z_∞ = 0 λ-fast-überall, also μ ⊥ λ, und H′ = dμ_ac/dλ = 0 fast überall. Ĝ und H unterscheiden sich nur durch Verkettung mit glatten Diffeomorphismen (spog_b, P_d), darum überträgt sich alles.
**Verwandtschaft:** Das ist das Kakutani-Dichotomie-Argument. Es entspricht auch Zygmunds Satz über lakunäre Riesz-Produkte ∏(1+a_k cos n_k x): Diese sind singulär, wenn Σa_k² = ∞.
Hier ist a_k ≡ Amplitude von P₄′ konstant, die „Frequenzen“ U_j(y) wachsen turmartig und sind nicht linear in y.

**#5.** Für einen singulären Anteil gilt Dμ = +∞ μ-fast-überall (Standard). Die lokale Dimension folgt aus #2: Es gilt log μ(B(y,h))/log h = 1 + O(n(h))/log h → 1, weil log h turmartig gegen −∞ geht und n(h) = slog(1/h).

**#6.** Es gilt log Π_n = Σ_{j<n} log(1+P₄′({U_j})) + (konvergente Reihe), und die Ziffern {U_j} sind unter λ fast unabhängig und gleichverteilt (Fehler δ_j wie in #4).
Mit dem starken Gesetz der grossen Zahlen folgt die Rate. Nach Jensen ist χ > 0, weil ∫P₄′ = 0. Lücke: das Gesetz der grossen Zahlen für die schwach abhängigen Ziffern ist nicht ausgeschrieben (Standard für ψ-mischende Folgen).

**#7.** Liegt U_j(I) für ein Intervall I in einem Fenster der Länge η, dann überdeckt U_{j+1}(I) eine Länge ≥ η·T_d′(U_j) ≫ 1. Darum kann man {U_{j+1}} in jedes Fenster lenken (Schachtelung).
Lenkt man in die Nähe der Minimalstelle von P₄′, fällt Π_n geometrisch, und es folgt H′ = 0. Bei der Maximalstelle folgt H′ = ∞. Lenkt man abwechselnd, bleibt Π_n beschränkt, und nach #3 gibt es keine Ableitung.

**#8.** Auf Rang 4 lautet der Schälschritt log_b∘exp_d = α·id. Seine „Phase“ ist konstant, die Ziffern-Verteilung des Konjugators ist gleichverteilt, und die Affinität ist 1. Damit passt der Rang 4 zur glatten (analytischen) Rang-4-Welle.
Ab Rang 5 ist die Phase des Schälschritts nicht konstant, und es folgt Singularität. Auf Rang 6 ist der Schälschritt G₅ = id + Φ₅, und dΦ₅ + dθ = dĜ ist bereits singulär. Schon eine einzige Ebene trägt dann die Affinität A_K → 0 für K → ∞.
Lücke: der Satz B auf Rang 6 setzt den Rang-6-Ebenensatz c₆ = Φ₅ + o(1) voraus. Den gibt es nur numerisch (R18-B15), und dazu kommt die Lücke für gerade Ränge.

## Skripte (jede Zahl reproduzierbar)

| Skript | liefert |
|---|---|
| `lib.py`, `lib6.py` | Kopien aus Runde 18 (Kneser-T aus Tschebyschow-Daten in `cache/`, kopiert, nicht neu gerechnet) |
| `s1_konstanten.py` → `daten/s1_konstanten.json` | Kennzahlen von P₄ für 9 Paare (#3, #6, Tabelle) |
| `s2_wackeln.py b d n dps θ…` → `daten/s2_*.json` | Wackel-Lemma (#1). Läufe: `exp(1) 2 3 80 0.05`, `exp(1) 2 3 100 0.02 0.1`, `exp(1) 3 2 80 0.02 0.04` (je ≤ 10 s) |
| `s3_hellinger.py` | Affinitäten A_K (#4), Gegenprobe mit M=64/128 |

## Vorschläge für Grapheinträge

* `bw.satz.wackel_lemma` — Differenzenquotient von H auf Skala s/U_n′ = Π_n·(1+W(U_n,s)) — satz — skizze (Lücke L2), numerisch auf O(1/U_n)
* `bw.satz.hoelder_eins` — Rang-5-Phase hat überall den Hölder-Exponenten 1, Modul h·e^{O(slog 1/h)} — satz — skizze
* `bw.satz.kein_endlicher_ableitungswert` — H′ ∈ {0, ∞} oder existiert nicht, an jeder Stelle — satz — bedingt/skizze
* `bw.satz.rang5_singulaer` — Rang-5-Konjugator und Ĝ singulär, Φ₅′ = −1 fast überall (Hellinger-Martingal) — satz — bedingt (L1, L2)
* `bw.beob.singulaer_unsichtbar` — Drift χ = −∫log(1+P₄′) ≈ Var(P₄′)/2 = 2.04·10⁻⁷ für (e,2); sichtbar erst nach ~5·10⁶ Turm-Ebenen — beobachtung — numerisch
* `bw.verm.rang5_nirgends_diffbar` — **Status ändern:** im Wortlaut widerlegt (bedingt), ersetzt durch bw.satz.rang5_singulaer — vermutung → negativ
* `bw.verm.dichotomie_raenge` — Konjugator glatt ⇔ Phase des Schälschritts konstant; singulär für alle r ≥ 5 — vermutung
* Kante: `bw.satz.rang5_singulaer` analog_zu Kakutani-Dichotomie / Zygmund (Riesz-Produkte); stuetzt_sich_auf bw.satz.phasengesetz_rang5 und R18-B7

## Offene Fragen

1. L2 sauber ausschreiben: gleichmässige Schranken für U_n″/U_n′² und (slog_b^n)″/(slog_b^n)′ an Turm-Höhen. Das sollte mit T_b′ wachsend und T″/T′² ≈ 1/T Routine sein.
2. Ist dĜ ein Rajchman-Mass, gehen also die Fourier-Stieltjes-Koeffizienten gegen 0? Bei linearen lakunären Riesz-Produkten nicht, hier sind die „Frequenzen“ U_j(y) nicht linear.
3. Exakte Hausdorff-Dimension der Trägermenge von dĜ: Sie ist 1, aber welche Mass-Funktion (Hausdorff-Funktion h·e^{−χ slog(1/h)}?) ist die richtige?
4. Gibt es eine natürliche, *glatte* Rang-5-Übersetzung, zum Beispiel über einen anderen Konjugator oder komplexe Basen mit gemeinsamem Fixpunkt (R18 Frage 5)? Nach #8 bräuchte sie einen wellenfreien Schälschritt.
5. Hängt χ(b,d) für d → b quadratisch von ln log_b d ab, wie Var(P₄′) ∝ \|a₁\|² nahelegt (lineare Antwort, bw.satz.antwortkern)? Wegen χ ≈ 4π²|a₁|² folgt dann χ ≈ 4π²·K₄²·(ln log_b d)², mit K₄ ≈ 1.94·10⁻⁴ … 1.98·10⁻⁴.
