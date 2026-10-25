# Runde 25.10.2026 — 06 Riemannfläche von sexp: Nebenblätter und Verzweigungspunkte (Kneser, Basis e)

## Kurz in einfacher Sprache

- Die offene Frage aus Runde 18 ist beantwortet: **Ja, die Verzweigungspunkte liegen dicht in ℂ.** Neben jedem Punkt der Ebene, auch neben
  einem ganz harmlosen Punkt wie 0.3, sitzen auf irgendeinem Nebenblatt beliebig nahe Verzweigungspunkte. Es sind aber nur abzählbar viele.
- Es gibt eine einfache Regel dafür, wo genau sie liegen: z ist ein Verzweigungsort genau dann, wenn die Bahn sexp(z+M) irgendwann
  **genau** auf die Bahn einer Zahl −2, −3, −4, … fällt (sexp(z+M) = sexp(M−n), n ≥ 2). n−1 ist die „Tiefe“: Tiefe 1 ist ein reiner Logarithmus, Tiefe 2 ein log log usw.
- Für jede feste Tiefe ist die Menge dagegen dünn (nirgends dicht). Dicht wird es erst, wenn man alle Tiefen zusammennimmt.
  Und jedes einzelne Blatt sieht in jedem beschränkten Gebiet nur endlich viele Verzweigungspunkte.
- Das Gittermuster gibt es: Oben in der Ebene liegen die Punkte in Reihen mit Abstand genau γ = 2πi/L = 4.44695 + 1.05794i.
  Die Abweichung vom exakten Gitter schrumpft von Glied zu Glied um den Faktor e^{2πiγ} ≈ −0.0012262 + 0.0004246i — vorhergesagt und gemessen.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| R1 | Blattstruktur: Die Riemannfläche R von sexp ist die **volle** Logarithmus-Überlagerung von (R−1) ohne Nullstellen. Die Werte aller regulären Zweige über z sind V(z) = exp⁻¹(V(z+1)∖{0}) = ⋃_K {alle Ketten log∘…∘log von sexp(z+K)} (Hauptblatt, keine Zwischenstufe 0) | **bewiesen** | Beweis §1 | neu |
| R2 | Charakterisierung: z ist Projektion eines Verzweigungspunkts ⇔ ∃M, n ≥ 2: sexp(z+M) = sexp(M−n) (Hauptblatt). Nullstellenmenge aller Zweige Z = {ζ : sexp(ζ+M) = sexp(M−1)}; alle Nullstellen einfach; Verzweigungsorte P = ⋃_{m≥1}(Z−m), über Z−m vom Typ log^[m] | **bewiesen** | §1; schliesst die Lücke von Runde 18 (B8: „Turm nur für n ≤ 2 ausgeführt“) | neu |
| R3 | Rechtsschranke: ζ ∈ Z, ζ+N im Hauptgebiet ⇒ Re sexp(ζ+N+k) ≥ sexp(N+k−1) für alle k ≥ 0. Insbesondere Re sexp(w) ≥ 0 für jede Nullstelle-im-Nebenblatt w (die Werte 2πik liegen genau auf dem Rand) | **bewiesen** | §2, Rückwärtsinduktion mit Re log u = log|u| ≥ log Re u | neu |
| R4 | Für jede feste Tiefe m ist Z−m nirgends dicht (Abschluss nirgends dicht) | **bewiesen** | §2 (aus R3 + Fatoumenge von exp leer) | neu |
| R5 | **P ist dicht in ℂ**; P ist abzählbar, also mager und Lebesgue-Nullmenge | **bewiesen** | §3 (Rückwärtsbahn von 1 unter exp dicht, Julia-Menge von exp = ℂ) | beantwortet r18 f5 / Offene Frage 1 |
| R6 | Jeder einzelne Zweig hat in jedem beschränkten Gebiet nur endlich viele Singularitäten | **bewiesen** | §3 | neu |
| R7 | Gitter: Lösungen von sexp(w) = 2πik oben bilden Reihen w_j mit w_{j+1}−w_j → γ = 2πi/L; Fehler δ_j = w_{j+1}−w_j−γ erfüllt δ_{j+1}/δ_j → e^{2πiγ} = −0.00122619+0.00042455i | **hergeleitet** (aus Kneser-Darstellung, Lücke: Restglied nicht abgeschätzt) + **numerisch** | gemessen −0.0012284+0.0004250i (k=1), −0.0012326+0.0004181i (k=2); w auf 22 Stellen gleich bei dps 40 und 50 (k=1) | neu (Gitter γ war aus B7 bekannt, das Fehlergesetz nicht) |
| R8 | Beispiele naher Verzweigungspunkte: Tiefe 1 bei 0.300141218682681 − 0.006243993284628i (Abstand 0.0062 von 0.3; sexp(z+4) = sexp(2)); Tiefe 2 bei 1.302765581658078 + 0.103084277039998i (Abstand 0.0041 von 1.3+0.1i); Tiefe 3 bei −0.518732284642257 + 0.215501282693660i (0.039 von −0.5+0.25i) | **numerisch** (eine Einstellung dps 40; Kettenprobe 10⁻¹¹⁹) | `s2_dicht.py`; Abstände für 0.3 nach Stufe s = 2, 3, 4: 1.24, 0.29, 0.0062 (überexponentiell kleiner) | neu |
| R9 | Annäherung über feste Tiefe ist nicht gleichmässig: nahe 0.3 liegt ein Punkt der **Tiefe 1** im Abstand 0.006, obwohl Tiefe-1-Punkte nirgends dicht sind | **numerisch** (Beobachtung) | R8 | — |
| R10 | Negativ: Eine „Gitterperiode“ von sexp selbst gibt es nicht. Verschieben um γ erhält die Klasse nur asymptotisch (Im → ∞), Verschieben um 1 wechselt die Klasse (Z+1 ist die Klasse von 0, nicht von −1) | **bewiesen** | §4 | — |

## 1. Blattstruktur und Charakterisierung (R1, R2)

Bezeichnungen: F_N(z) := sexp(z+N) auf dem Hauptblatt, holomorph auf dom_N = {z : z+N ∉ (−∞,−2]}; E_k := exp^[k](0) = sexp(k−1), also E_0 = 0, E_1 = 1, E_2 = e, …

**(a) Volle Logarithmus-Überlagerung.** Sei R₁ die Fläche von z ↦ sexp(z+1) (= R um −1 verschoben) und S die Funktion darauf. Fortsetzung von sexp längs γ existiert
genau dann, wenn S längs γ fortsetzbar ist und nicht 0 wird; dann ist sexp = log S (aus sexp(z+1) = exp sexp(z) und Eindeutigkeit). Also ist R die Zusammenhangskomponente
des Hauptkeims in {(p, λ) : p ∈ R₁, e^λ = S(p)}. Diese Menge ist zusammenhängend: R₁ ohne die (diskreten) Nullstellen von S ist zusammenhängend, und ein Umlauf um die einfache
Nullstelle z = −2 von S (Hauptblatt, sexp′(−1) ≠ 0) ändert λ um 2πi. Also liegen über jedem p alle Logarithmen von S(p) in R.

**(b) Endliche Ketten.** Jeder Punkt von R über z wird über einen beschränkten Weg γ erreicht. Für N mit γ+N ⊂ {Re > −2} bleibt sexp(·+N) längs γ auf dem Hauptblatt (geschlitzte
Ebene ist einfach zusammenhängend). Also ist jeder Zweigwert über z eine Kette h_0 = log h_1, …, h_{N−1} = log h_N, h_N = F_N(z), und nach (a) kommt jede solche Kette ohne Zwischenwert 0 vor.

**(c) Nullstellen.** Ein Zweig ist bei ζ null ⇔ h_0 = 0 ⇔ h_i = E_i für alle i (exp ist eindeutig) ⇔ F_K(ζ) = E_K für ein K. Die Zwischenwerte E_i ≥ 1 sind nie 0, die Kette ist also zulässig.
Ableitung des Zweigs: h_0′ = F_K′(ζ)/(E_1E_2⋯E_K) ≠ 0 (sexp′ nullstellenfrei, kn.satz.ableitung_nullstellenfrei) — einfache Nullstelle. Daraus Z = {ζ : ∃M sexp(ζ+M) = sexp(M−1)}
(Bedingung ist stabil in M). Über ζ−1 hat der Zweig log(h_0) einen reinen logarithmischen Verzweigungspunkt; über ζ−m den Typ log^[m] (jeder weitere Logarithmus eines
Werts → ∞ bleibt singulär). Umgekehrt ist eine Singularität eines Zweigs über z eine Stelle, an der ein Zwischenwert der Kette 0 wird, also z+m ∈ Z für ein m ≥ 1. Ergebnis:
P = ⋃_{m≥1}(Z−m) = {z : ∃M, n ≥ 2 : sexp(z+M) = sexp(M−n)}. ∎ (Runde 18 hatte „⊆“ und die Erreichbarkeit nur für n ≤ 2; (a) liefert sie für alle n.)

**Monodromie.** Um einen Punkt von Z−1 wirkt sie als λ ↦ λ+2πi auf dem betroffenen Logarithmus-Turm (unendlich zyklisch); um einen Punkt von Z−m wie beim log^[m]:
der innerste Index wird erhöht, dabei entstehen weitere Verzweigungen über den Stellen, wo ein Zwischenwert 1, e, … wird. Die Faser über einem regulären Punkt ist abzählbar
unendlich (Ketten mit endlich vielen Zweignummern ≠ 0).

## 2. Rechtsschranke und Dünnheit fester Tiefe (R3, R4)

Sei exp^[k](c) = E_{N+k} für ein k. Mit a_i = exp^[i](c): Re a_{k−1} = log E_{N+k} = E_{N+k−1}. Ist Re a_{i+1} ≥ E_{N+i+1} ≥ 1, so Re a_i = log|a_{i+1}| ≥ log Re a_{i+1} ≥ E_{N+i}.
Also liegt G_N := {c : ∃k exp^[k](c) = E_{N+k}} in A_N := {c : Re exp^[i](c) ≥ E_{N+i} ∀i ≥ 0}. A_N ist abgeschlossen und hat kein Inneres (auf einer offenen Menge in A_N
gingen alle Iterierten gleichmässig gegen ∞, die Familie wäre normal, aber die Fatoumenge von exp ist leer, Misiurewicz 1981). Weil Z ∩ dom_N = F_N⁻¹(G_N) ⊂ F_N⁻¹(A_N) und F_N offen ist,
ist der Abschluss von Z in jedem dom_N nirgends dicht; Verschiebungen erhalten das. Für N = 0: Re sexp(w) ≥ 0 für alle w ∈ Z im Hauptgebiet; die Runde-18-Punkte
(sexp(w) = 2πik) liegen genau auf dem Rand Re = 0. ∎

## 3. Dichtheit (R5) und Endlichkeit je Zweig (R6)

Ist c ∈ exp^{−[s]}(1), dann F_{N+s}(z) = E_1 für jedes z mit F_N(z) = c, also z ~ −(N+s) und z ∈ P für N+s ≥ 2. Die Rückwärtsbahn ⋃_s exp^{−[s]}(1) ist dicht in J(exp) = ℂ
(1 ist kein Ausnahmewert; Satz über Rückwärtsbahnen ganzer Funktionen, Bergweiler 1993, Bull. AMS, Thm. 4). Ein Urbild einer dichten Menge unter der offenen Abbildung F_N ist dicht in dom_N;
⋃ dom_N = ℂ. Also ist P dicht. Abzählbar: jede Bedingung sexp(z+M) = sexp(M−n) hat auf jedem kompakten Stück von dom_M endlich viele Lösungen. ∎

Endlichkeit je Zweig: Über einer beschränkten Scheibe D ⊂ dom_N ist ein Zweig eine Kette von Logarithmen von F_{N'} (N' ≥ N fest für diesen Zweig); singulär nur, wo F_{N'} = E_k
mit k < N' — endlich viele Punkte in D̄. Erst die Vereinigung über alle Zweige ist dicht. ∎

## 4. Gittermuster (R7, R10)

Oben gilt sexp = Ψ(χ), χ(z) = κ e^{L(z+θ(z)−a₀)}, θ 1-periodisch, θ(z) = a₀ + a₁e^{2πiz} + … (Kneser). sexp(w) = c mit festem Urbild t ∈ Ψ⁻¹(c) heisst
w + θ(w) = const + jγ, γ = 2πi/L. Also w_{j+1} − w_j − γ = −(θ(w_{j+1}) − θ(w_j)) ≈ a₁ e^{2πiw_j}, und das Verhältnis aufeinanderfolgender Fehler ist e^{2πiγ}
(|e^{2πiγ}| = e^{−2π Im γ} = 0.0012976). Gemessen (`s1_familie.py`, Motor-Newton):

| k | j | w_j (sexp(w_j) = 2πik) | sexp′(w_j) |
|---|---|---|---|
| 1 | 0 | 2.004312482594688003944 + 0.709413540869403626679i | −11.600149 + 6.494471i |
| 1 | 1 | 6.451336487382556769423 + 1.766298794622866278511i | −11.520432 + 6.456817i |
| 1 | 2 | 10.89828756564794741946 + 2.824240112601807370779i | −11.520514 + 6.456897i |
| 2 | 0 | 2.155871087462956998251 + 0.490641727266906844928i | −30.486886 + 25.869768i |
| 2 | 1 | 6.60644184445492460808 + 1.546446382998737726742i | −29.490303 + 26.208759i |
| 2 | 2 | 11.05338899528792567913 + 2.60439051951416481697i | −29.491659 + 26.208780i |

δ₀ = 7.328·10⁻⁵ − 1.0547·10⁻³i, δ₁ = 3.582·10⁻⁷ + 1.3268·10⁻⁶i (k=1); δ₁/δ₀ = −0.0012284 + 0.0004250i gegen e^{2πiγ} = −0.0012262 + 0.0004246i.
Bild: Nahe i∞ ist Z ≈ (1/L)·log(Ψ⁻¹(G₀)/κ) + γℤ, und P = Z − {1,2,…}. Das Gitter Λ = ℤ + γℤ ordnet also die Punkte oben; γ erhält die Klasse (asymptotisch),
1 verschiebt die Klasse um eins (Z+1 ist die Klasse von 0). sexp′ entlang einer Reihe strebt gegen L·t·Ψ′(t).

## Skripte (alle in diesem Ordner)
`motor.py` (Kopie aus 03_basis_ebene, eigener Speicher `daten/motor_cache.json`, Speichern mit Wiederholung bei Windows-Dateisperre),
`s1_familie.py DPS` (Reihen, Ergebnis `daten/familie_40.json`; Gegenprobe dps 50 in `daten/s1_dps50.log`, k=1 auf alle 22 gedruckten Stellen gleich, Lauf danach an einer
Dateisperre abgebrochen), `s2_dicht.py 40` (nahe Verzweigungspunkte, `daten/dicht_40.json`).

## Vorschläge für Grapheinträge
- kn.satz.riemannflaeche_logueberlagerung — „Fläche von sexp ist volle Log-Überlagerung von (R−1); Zweigwerte = alle Log-Ketten von sexp(z+K)“ — satz — bewiesen
- kn.satz.verzweigungsorte_charakterisierung — „z verzweigt ⇔ sexp(z+M) = sexp(M−n), n ≥ 2; Typ log^[n−1]“ — satz — bewiesen
- kn.satz.nebenblatt_nullstellen_rechts — „Nullstellen der Nebenzweige: Re sexp(ζ+N+k) ≥ sexp(N+k−1); feste Tiefe nirgends dicht“ — satz — bewiesen
- kn.satz.verzweigung_dicht — „Verzweigungsorte aller Blätter dicht in ℂ, abzählbar; je Zweig lokal endlich“ — satz — bewiesen
- kn.beob.verzweigung_gitter — „Reihen mit Abstand γ = 2πi/L, Fehlerverhältnis e^{2πiγ}“ — beobachtung — hergeleitet + numerisch
- Kante: kn.satz.verzweigung_dicht beantwortet r18-Frage f5 (kn.satz.nebenblatt_zweigpunkte).

## Offene Fragen
1. Abstand d_m(z) vom Punkt z zum nächsten Verzweigungsort der Tiefe ≤ m: vermutlich d_m ≍ 1/|F_m′(z)| (überexponentiell); Schranke von unten fehlt.
2. Hausdorff-Dimension des Abschlusses von Z (Tiefe 1): Er liegt in sexp⁻¹(A₀), A₀ ist ein „Haarbündel“ von exp; Dimension von A₀?
3. Gibt es Punkte von Z im Streifen −2 < Re ≤ −1? (R3 sagt dort Re sexp ≥ 0, sexp(ζ+1) ∈ 1 + … ; nicht entschieden.)
4. Andere Basen 1 < b < e^{1/e}: dort hat exp_b eine Fatoumenge, also könnte P **nicht** dicht sein — Vermutung, nicht geprüft.
