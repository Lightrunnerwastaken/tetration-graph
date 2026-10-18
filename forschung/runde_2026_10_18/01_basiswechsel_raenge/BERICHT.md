# Basiswechsel auf höheren Rängen (Pentation, Hexation, Rang 7) — Bericht 01, Runde 18.10.2026

## Kurz in einfacher Sprache

Bei der Tetration weiss man: Übersetzt man die Höhe einer Basis d in die Höhe einer Basis b (slog_b(T_d(x))), entsteht
„x + Konstante + eine winzige Welle“. Ich habe dasselbe eine, zwei und drei Stufen höher gerechnet (Pentation, Hexation, Rang 7).

1. **Das Phasengesetz gilt auch dort.** spog_b(P_d(x)) − x läuft gegen μ + Welle; die Welle ist 1-periodisch.
   Der Beweis ist sogar einfacher als bei der Tetration: Man schält statt mit Logarithmen mit dem Rang-4-Basiswechsel ab.
   Aus dem Phasengesetz auf Rang 4 folgt es so für jeden höheren Rang.
2. **Die Welle ist auf Rang 5 etwa 60-mal grösser als auf Rang 4.** Für (e, 2) misst die Spannweite 0.0193 statt 0.00029.
   Auf Rang 6 und 7 wächst sie nicht weiter (0.017 und 0.014), aber die Welle wird von Rang zu Rang rauer.
3. **Die Rang-5-Welle ist nicht glatt.** Die Rang-4-Welle steckt auf jeder tieferen Stufe wieder in der Rang-5-Welle,
   mit immer schnellerer Schwingung (Frequenz wie ein Turm). Die „Ableitung“ der Rang-5-Phase ist darum ein unendliches Produkt
   ∏(1 + Φ₄′(P_d(θ+j))). Die Rechnung bestätigt das auf 3·10⁻⁹. Vermutlich ist die Phase ab Rang 5 nirgends differenzierbar,
   ähnlich wie eine Weierstrass-Funktion. Ihre Fourier-Moden fallen nur noch wie eine Potenz ab, auf Rang 4 dagegen exponentiell.
   Bei der Tetration streitet man noch, ob die Welle analytisch oder nur unendlich oft differenzierbar ist. Ab Rang 5 ist sie
   nach diesen Befunden nicht einmal einmal differenzierbar.
4. **Fixpunkte:** Der Konjugator H, der T_d in T_b überführt, schickt den Fixpunkt L_d immer auf L_b. Die Welle misst genau,
   wie weit H dort von einer reinen Potenz abweicht (Schröder-Koordinaten). Zwei verschiedene reelle Basen haben nie einen
   gemeinsamen Fixpunkt von T. Eine wellenfreie reelle Basis gibt es daher nicht. Gemessen ist die Welle immer ≈ 0.010–0.015·|ln log_b d|.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| 1 | Abschäl-Identität: spog_b(P_d(x+1)) = 1 + spog_b(F_{b,d}(P_d(x))), F = slog_b∘T_d; allgemein A_{r,b}∘H_{r,d}(x+1) = 1 + A_{r,b}(G_{r−1}(H_{r,d}(x))) | **bewiesen** (Abel-Gleichung) | — | Schreibweise |
| 2 | Rang-5-Phasengesetz: H := lim slog_b^n∘T_d^n existiert, T_b∘H = H∘T_d, Ĝ := spog_b∘H∘P_d erfüllt Ĝ(x+1)=Ĝ(x)+1 **für alle x**, G(x) − Ĝ(x) → 0 superexponentiell; Φ₅ = Ĝ − id stetig, 1-periodisch | **bewiesen unter** bw.satz.phasengesetz (nur Beschränktheit von F(u)−u gebraucht) und den gemessenen Pentationsvoraussetzungen (λ>1, L einziger Fixpunkt) | Konvergenz n=0,1,2,∞ für (e,2): −0.055, −0.415, −0.839, −0.93476; unabhängig vom Abbruch X (40/60/75) auf 10⁻⁴⁰; Motor-`ipentx` stimmt mit eigener Kette auf 30 Stellen | vermutlich neu (Forum tid 970 erwähnt „base-change pentation“ nur spekulativ) |
| 3 | Induktion über den Rang: Phasengesetz auf Rang r−1 (beschränkt) ⇒ Phasengesetz auf Rang r (r ≥ 5) | **hergeleitet mit Lücke** (für gerade r: Kneser-Rang, Reellwertigkeit/Monotonie der Umkehrung nur gemessen) | numerisch für r=6 (Kneser-Hexation) und r=7 (regulär am linken Fixpunkt der Hexation): periodisch, X-unabhängig auf 10⁻³⁰ | vermutlich neu |
| 4 | Kozykel und inverse Gleichung auf Rang 5: Φ_{b,d}(θ) = Φ_{c,d}(θ) + Φ_{b,c}(θ+Φ_{c,d}(θ)); μ₅(b,d) = −μ₅(d,b) | **bewiesen** (bedingt wie #2; Antisymmetrie wie bw.satz.antisymmetrie) | Kozykel-Rest (2,e,3) < 10⁻¹¹⁷; μ₅(e,2)+μ₅(2,e) = 2·10⁻⁸ bei M=128 (= Quadraturfehler) | Übertragung |
| 5 | Fixpunkt-Deutung: H(L_d) = L_b, und σ_b(H(y)) = C·σ_d(y)^κ·exp(ln λ_b·P₅(spog_d y)), κ = ln λ_b/ln λ_d. Welle = 0 ⇔ H ist in Schröder-Koordinaten eine reine Potenz | **bewiesen** (Umformung von #2) | Ĝ(x) = x+Φ₅({x}) geprüft bis x = −15.7 auf 10⁻¹¹⁰; κ(e,2) = 1.13541; G(x) selbst → spog_b(L_d) = −3.41245 (verschiedene Fixpunkte ⇒ Konstante, wie Rang 4) | vermutlich neu |
| 6 | Kein gemeinsamer Fixpunkt reeller Basen: L_b fällt streng mit b (1.8: −1.6938, 2: −1.7439, e: −1.8504, 3: −1.8745, 4: −1.9252, 5: −1.9501, 10: −1.9859); also keine wellenfreien reellen Paare | **numerisch** (7 + 4 Basen; Monotonie nicht bewiesen) | Fixpunkte gegen Motor-`pentx` auf 10⁻⁷⁰ | — |
| 7 | Ebenen-Zerlegung: Φ₅(θ) = spog_b(slog_b^m(U_m + c(U_m))) − θ exakt, U_m = P_d(θ+m), c = H − id, und c(U) − Φ₄({U}) → 0 superexponentiell | **bewiesen** (bedingt wie #2) | c(U)−μ₄ = P₄(U) auf alle Stellen ab U ≈ 5.4 (d=2); bei U=2: 0.18 (Anlauf) | vermutlich neu |
| 8 | Einspeisegesetz: Φ₅ − Φ₅^{(m)} = (1+Φ₅′(θ))·P₄(P_d(θ+m)) / P_d′(θ+m) + O(P₄²) (Φ₅^{(m)}: tiefere Welle durch ihr Mittel ersetzt) | **hergeleitet** (Kettenregel, erste Ordnung), **numerisch** bestätigt | Quotient gemessen/vorhergesagt = 1 ± 10⁻³ an 12 Stellen θ ∈ [0.25, 0.94] | neu |
| 9 | Produktformel: H_{n+1}′(y)/H_n′(y) = F′(U_n)·(1+o(1)) = 1 + Φ₄′(U_n) + o(1); formal H′(y) = ∏_j F′(T_d^j y) | **bewiesen** (exakte Kettenregel; o(1) superexponentiell) | θ=0.05, U₃ = 5.3·10⁸: 8.1454952·10⁻⁴ gegen Φ₄′(U₃) = 8.1454722·10⁻⁴ | neu |
| 10 | Skalensprung der Ableitung: Differenzenquotient von Φ₅ bei θ₀=0.05 ist bis h≈2⁻¹⁸ konstant (Ĝ′ = 1.0657272), springt bei h≈1/U₃′ und ist ab h=2⁻⁵⁰ wieder konstant (1.0665953); Verhältnis 1.000815 = 1+Φ₄′(U₃) | **numerisch** (dps 120) | auf 10⁻⁶ | neu |
| 11 | **Vermutung:** Φ₅ (und alle Φ_r, r ≥ 5) ist stetig, aber nirgends (mindestens fast nirgends) differenzierbar; Stetigkeitsmodul ≈ h·(Anzahl Ebenen mit U_j′h<1) ~ h·slog(1/h) | **Vermutung** | #9, #10; Fourier-Schwanz algebraisch ~k^{−1.3} (Rang 5, k ≤ 64) statt exponentiell (Rang 4: |a₆₄| = 10⁻⁴⁹) | neu |
| 12 | Wellengrösse: |a₁|(Φ₅)/|a₁|(Φ₄) ≈ 54–77 für 13 Paare; Grenzwert d→b: r(e) = 62.1, r(2) ≈ 75 | **numerisch** | (e,2.7): 62.31, (e,2.8): 61.95; (2,1.9): 76.9, (2,2.1): 72.4 | neu |
| 13 | Lineare Antwort: |a₁(Φ₅)| ≈ K₅(b)·|ln log_b d|, K₅(e) = 0.01207, K₅(2) ≈ 0.0148 (basisabhängig); Rang 4: K₄ ≈ 1.94·10⁻⁴ … 1.98·10⁻⁴ (fast universell, vgl. bw.satz.antwortkern) | **numerisch** | Paare mit |ε| = 0.007 … 0.08 | neu (Rang 5) |
| 14 | Die Rang-5-Welle stammt fast ganz (~99 %) aus dem **Anlauf** von F_{b,d} bei kleinen Höhen (u ∈ [d, ᵈd]), nicht aus der Rang-4-Welle (Beitrag ≲ 10⁻⁴ von 10⁻²) | **numerisch** (über #7/#8) | Ebene 1: S₁·(c−μ₄) = 0.09 … 0.002; P₄-Anteil ≤ 7·10⁻⁵ | neu |
| 15 | Rang 6/7: Wellen gleich gross wie Rang 5 (Spannweite (e,2): 0.0193 → 0.0174 → 0.0144), aber rauer: |a₂/a₁| = 0.0096 (r=4), 0.090 (5), 0.33 (6), 0.51 (7); c₆(V) = H₆(V)−V = Φ₅({V}) exakt ab V ≈ 3.3 (Rang-5-Welle wird eingespeist) | **numerisch** | Antisymmetrie-Rest 4·10⁻⁶ (r=6, M=64), 3·10⁻⁵ (r=7, M=32) = Quadraturfehler | neu |
| 16 | Mittel über die Ränge (e,2): μ₄…μ₇ = −1.1284038, −0.9432744(2), −0.920180(4), −0.91917(3); (e,3): 0.1924223, 0.1258095(1), 0.102464(1), 0.08882(2) | **numerisch** (Fehler in Klammern aus M-Vergleich und Antisymmetrie) | — | neu |
| 17 | Hub-Additivität μ(2,3) ≈ μ(e,3) − μ(e,2): Rest 3·10⁻¹⁰ (r=4), 1·10⁻⁵ (r=5), 2·10⁻⁵ (r=6), Ordnung Welle² wie vom Kozykel verlangt | **numerisch** | — | Übertragung |

## Methode (alles reproduzierbar)

* **Tetration:** Kneser-T_b aus dem Motor (fatou-Fork, nur gelesen, eigener Zwischenspeicher `cache/`, auch `FATOU_CACHE_DIR`) an 150 Tschebyschow-Knoten
  auf [−1,0], 80 Stellen (`holen.py`). T, T′, slog über exp/log-Verschiebung und Newton (`lib.py`). Prüfung: μ_{e,2} = −1.128403776628924009879217902036 (alle 30 Stellen wie bw.beob.m2).
* **Pentation:** selbst gebaut, regulär am Fixpunkt L_b ∈ (−2,−1): Taylor von T an L+1 aus dem Tschebyschow-Polynom, log-Reihe, Schröder-Reihe (48 Glieder).
  Gegenprobe gegen Motor-`pentx` bei 7 Basen: 10⁻⁷⁰ (`t5_basen.py`); e↑↑↑½ = 1.6323247404360631 wie kn.beob.pentation_werte.
* **Rang-4-Grundbaustein** R(u) = slog_b(T_d(u)) für beliebig grosse u: exaktes Abschälen log_b(c·exp_d w) = log_b c + α w mit log1p-Korrekturen;
  Glieder unter 10^(−10⁶) entfallen. Abbruch bei E_{j₀}, j₀ ≤ 2, damit kein Logarithmus einer negativen Zahl entsteht.
* **Rang 5:** H(y) = slog_b^{n−1}(R(U_{n−1})), Kette U_j = T_d^j(y) bis T_d(U) > 10^X (Rest ~10^−X), dps 120–150.
* **Rang 6:** Kneser-Hexation (`hexx`) und **Rang 7** `hyper(7,·)` (regulär am linken Fixpunkt der Hexation) aus dem Motor an 120 Knoten, 60 Stellen; Umkehrungen per Newton;
  H₆, H₇ wie H mit G₅ bzw. G₆ als Abschälschritt (`lib6.py`, `lib7.py`). Gegenproben: hex(x+1) = P(hex(x)) auf 10⁻⁶⁹, Periodizität, X-Unabhängigkeit.
* Mittel per Trapezregel; wegen der rauen Wellen konvergiert sie auf Rang ≥ 5 nur langsam (M=32/64/128: −0.94327417/−0.94327429/−0.94327436).
  Die Fehlerangaben in #16 kommen daher.

| Skript | liefert |
|---|---|
| `t0_pruefen.py`, `t2_gegenprobe.py`, `t5_basen.py` | Grundprüfungen, Motor-Gegenproben |
| `t1_rang5.py`, `t3_fein.py` | erste Rang-5-Werte, Feinabtastung (Wellen in U₂) |
| `e1_paare.py M b d …` → `daten/e1_*.json` | Φ₄, Φ₅, Moden (#12, #13, #16) |
| `e2_rang6.py`, `e3_rang7.py` → `daten/e2_*`, `e3_*` | Rang 6, 7 (#15, #16) |
| `t4_ebenen.py`, `t11_gewicht.py` | #7, #8 |
| `t8_produkt.py` | #9, #10 |
| `t9_kozykel.py` | #4 |
| `t13_fixpunkt.py` | #5 |
| `t10_rang6_skalen.py`, `t14_rang6_ebenen.py` | Rang-6-Skalen, c₆ = Φ₅ |
| `zusammenfassung.py` | Tabellen aus `daten/`, Bild `wellen_raenge.png` |

## Beweise und Details

**Bezeichnungen.** b,d > η₅ = 1.73736 (sonst bleibt die Pentation beschränkt). T = Kneser-Tetration. P_b = Pentation regulär an L_b (kn.konstr.pentation),
spog_b = P_b⁻¹ auf (L_b,∞), spog_b∘T_b = spog_b + 1. F := F_{b,d} = slog_b∘T_d. Rang-4-Phasengesetz (bw.satz.phasengesetz): F(u) = u + Φ₄({u}) + ε(u), ε → 0.

**#1.** spog_b(T_d z) = spog_b(slog_b(T_d z)) + 1 = 1 + spog_b(F(z)); mit z = P_d(x) folgt die Identität. Iteriert: G(x+n) − n = spog_b(H_n(P_d x)), H_n = slog_b^n∘T_d^n.
Das ist das „Abschälen über die Tetration“. Auf Rang 4 ist der Schälschritt log_b∘exp_d = α·id **linear**.
Ab Rang 5 ist er F = id + μ₄ + P₄ + Anlauf: eine **Verschiebung mit Welle**. Daraus folgen alle Unterschiede.

**#2 (Beweis).** H_{n+1}(y) − H_n(y) = slog_b^n(F(U_n)) − slog_b^n(U_n) mit U_n = T_d^n y. Nach dem Mittelwertsatz ist das (slog_b^n)′(ξ)·(F(U_n) − U_n).
|F(u) − u| ≤ C für u ≥ u₀ (Phasengesetz). Per Induktion bleiben die Zwischenpunkte slog_b^i(ξ) ≥ U_{n−i} − C′.
T_b′ wächst auf [0,∞), also fällt slog_b′ dort, und es folgt (slog_b^n)′(ξ) ≤ slog_b′(U_n − C′).
U_n → ∞ turmartig (L_d einziger Fixpunkt, T_d(y) > y für y > L_d). Darum ist Σ|H_{n+1} − H_n| < ∞, gleichmässig auf Kompakta von (L_d,∞).
H_n(T_d y) = T_b(H_{n+1}(y)) gibt im Limes H∘T_d = T_b∘H. Daraus folgt Ĝ(x+1) = spog_b(T_b(H(P_d x))) = Ĝ(x) + 1,
und G(n+θ) − n = spog_b(H_n(P_dθ)) → Ĝ(θ). Benutzt wurde vom Rang 4 nur, dass F(u) − u beschränkt ist.
Strenge Monotonie von H kommt aus dem Gegenpaar (wie bw.satz.homoeo).
**#3:** Derselbe Beweis eine Stufe höher. Er braucht nur, dass G_{r−1}(u) − u beschränkt ist (aus dem Rang-(r−1)-Gesetz) und dass H_{r−1,b}′ → ∞ wächst.
**Lücke für gerade r ≥ 6:** Reellwertigkeit und Monotonie der Kneser-Umkehrung auf dem benutzten Bereich und Eindeutigkeit des Fixpunkts sind nur gemessen.

**#5.** Ĝ ist exakt (nicht nur asymptotisch) x + Φ₅({x}). Für x → −∞ ist P_d(x) = Ψ_d(σ_d(1)λ_d^x) → L_d und H(P_d x) → L_b.
Einsetzen von spog = log(σ/σ(1))/log λ in Ĝ(x) = x + Φ₅({x}) gibt die Schröder-Form. Sie ist eine reine Potenz genau dann, wenn Φ₅ konstant ist.
Die Welle misst also, wie „nicht-regulär“ der von +∞ kommende Konjugator H am Fixpunkt ist. Das ist das Gegenstück zu kl.satz.basiswechsel_gleiches_L:
Dort sorgt ein gemeinsamer Fixpunkt für eine wellenfreie Übersetzung. Hier verbindet H die beiden Fixpunkte, aber nicht regulär.

**#7–#9.** H = slog_b^m∘H∘T_d^m gibt die exakte Zerlegung.
c(U) − F(U) + U = slog_b(T_dU + c(T_dU)) − slog_b(T_dU) = O(slog_b′(T_dU)). Mit F(U) − U − Φ₄({U}) = ε(U) folgt #7.
Die Kettenregel liefert exakt H_{n+1}′/H_n′ = F′(U_n)·(slog_b^n)′(F(U_n))/(slog_b^n)′(U_n).
Der zweite Faktor ist 1 + O(1/U_n), und F′ = 1 + Φ₄′ + o(1) (C¹-Phasengesetz, gemessen). Das gibt #9.
Ableiten nach c_m mit spog_b′(H)H′ = Ĝ′/P_d′ ergibt das Gewicht S_m = Ĝ′(θ)/(P_d′(θ+m)(1+c′(U_m))) in #8.

**Warum #11 plausibel ist.** Ebene j trägt die Schwingung P₄(U_j)/U_j′ mit der Frequenz U_j′. Ihr Beitrag zur Steigung ist P₄′(U_j) und nimmt mit j **nicht** ab.
Das entspricht dem Grenzfall a·b = 1 der Weierstrass-Funktion, nur mit turmartigen statt geometrischen Lücken.
Hardy: dort nirgends differenzierbar. Ob H′(y) = ∏(1+Φ₄′(T_d^j y)+o(1)) existiert, hängt an Σ_j Φ₄′(T_d^j y).
Für fast alle y sollten die Brüche {T_d^j y} gleichverteilt und „fast unabhängig“ sein (Koksma, Salem–Zygmund für sehr lückenhafte Folgen).
Dann divergiert die Reihe wie eine Irrfahrt. **Lücken des Beweises:** (a) Koksmas metrischer Satz für die Folge T_d^j; (b) ein Lemma vom Hardy-Typ:
Die Differenzenquotienten folgen den Teilprodukten, wie in #10 gesehen. Messbar ist nur eine Ebene mehr (U₄ ist schon ~2^(5·10⁸)).
Rang 4 hat diese Struktur nicht: Der Schälschritt log_b∘exp_d hat die konstante Ableitung α, das Produkt konvergiert absolut.
Das passt zu bw.verm.streifen bzw. bw.verm.nirgends_analytisch (dort geht es um C^∞ gegen analytisch).

**#12/#14.** Die Rang-5-Welle (erste Mode ~5·10⁻³ bei (e,2)) entsteht in den Ebenen 0–1, also bei y = P_d(θ) ∈ [1,d] und U₁ ∈ [d, ᵈd].
Dort ist F noch weit von seiner Asymptotik entfernt: c(U₁) − μ₄ = 0.18 bei U₁ = 2. Die Rang-4-Welle (1.4·10⁻⁴) trägt höchstens S₁·|P₄| ≈ 7·10⁻⁵ bei.
Antwort auf „Komposition der Wellen?“: **Nicht in führender Ordnung.** Die Rang-4-Welle sitzt in allen tieferen Ebenen, aber klein (Gewicht 1/P_d′(θ+m)).
Sie bestimmt die **Glattheit**, nicht die **Grösse**.

## Fehlschläge / Ehrliches

* Mittelwerte μ₅…μ₇ sind nur auf 10⁻⁷ … 10⁻⁵ genau. Grund ist die Rauheit #11 (Trapezregel konvergiert algebraisch). Eine Zerlegungsquadratur (Ebene 2 im U-Gitter) wäre der nächste Schritt.
* Eine geschlossene Form für das Verhältnis r(b) (62 bei e, 75 bei 2) oder für K₅(b) habe ich nicht gefunden. Die Nähe zu λ_e/λ₂ = 1.249 gegenüber K₅(2)/K₅(e) = 1.22 ist nicht belastbar.
* Mein erster Rang-4-Baustein nahm den Logarithmus negativer Zahlen (b=e, d=1.8) und rechnete einen Turm 10↑↑3.x direkt aus (Absturz). Beides ist behoben (Abbruch bei E_{j₀}).
  Ergebnisse vor der Korrektur wurden für (e,2) und (e,3) bitgleich bestätigt.
* `timeout` in Git-Bash beendet Windows-Prozesse nicht. Ein hängender Lauf wurde von Hand per PID beendet; danach habe ich jeden Lauf geprüft.

## Vorschläge für Grapheinträge

* `bw.satz.phasengesetz_rang5` — Phasengesetz für die Pentation über den Konjugator H = lim slog_b^n∘T_d^n — satz — bedingt (bw.satz.phasengesetz)
* `bw.satz.phasengesetz_induktion` — Phasengesetz auf Rang r folgt aus Rang r−1 — satz — skizze (Lücke: gerade Ränge)
* `bw.prop.abschaelen_rang` — Abschäl-Identität A_{r,b}∘H_{r,d}(x+1) = 1 + A_{r,b}∘G_{r−1}∘H_{r,d}(x) — satz — bewiesen
* `bw.satz.fixpunkt_konjugator` — Welle = Nicht-Regularität von H am Fixpunkt (Schröder-Form) — satz — bedingt
* `bw.satz.produktformel_rang5` — H_{n+1}′/H_n′ = 1 + Φ₄′(T_d^n y) + o(1) — satz — bedingt; dazu `bw.beob.skalensprung` — numerisch
* `bw.verm.rang5_nirgends_diffbar` — Rang-5-Phase stetig, nirgends differenzierbar (Weierstrass-Typ mit Turm-Lücken) — vermutung
* `bw.beob.wellen_raenge` — Wellengrösse und Moden Rang 4–7, r(e) = 62.1, Rauheit wächst mit dem Rang — beobachtung — numerisch
* `bw.beob.anlauf_ursache` — Rang-5-Welle stammt aus dem Anlauf von F_{b,d}, nicht aus der Rang-4-Welle — beobachtung — numerisch
* `bw.beob.mu_raenge` — μ_r(e,2), μ_r(e,3) für r = 4…7 — beobachtung — numerisch
* `kn.beob.fixpunkt_monoton` — L_b fällt streng in b (keine gemeinsamen Fixpunkte reeller Basen) — beobachtung — numerisch

## Offene Fragen

1. Beweis von #11: Gleichverteilung von {T_d^j(y)} für fast alle y, plus ein Hardy-Lemma. Ist Φ₅ in der Zygmund-Klasse?
2. Was ist lim_{r→∞} μ_r(b,d)? Für (e,2) scheint die Folge schnell zu konvergieren (Abstände 0.185, 0.023, 0.001), für (e,3) nicht so klar.
3. Geschlossene Formel für den Rang-5-Antwortkern K₅ (Gegenstück zu bw.satz.antwortkern), mit Ṗ_b/P_b′ und dem Anlauf von F.
4. Basen unter η₅ (beschränkte Pentation): Konstante statt Phasengesetz wie bei Rang 4 unter η? Und Hexations-Basiswechsel mit gemeinsamen Pentations-Fixpunkten.
5. Komplexe Basen mit gemeinsamem Fixpunkt von T_b und T_d: Ist dann eine wellenfreie Rang-5-Übersetzung möglich (Gegenstück zu kl.frage.gleiches_L_kneser)?

Quelle zur Neuheitsprüfung (eine Suche): [Tetration Forum tid 970](https://tetrationforum.org/printthread.php?tid=970) („base-change pentation“, nur Spekulation),
[tid 236](https://www.tetrationforum.org/showthread.php?tid=236), [tid 937](https://tetrationforum.org/printthread.php?tid=937).
