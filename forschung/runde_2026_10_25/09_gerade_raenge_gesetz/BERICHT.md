# Runde 25.10.2026 · 09 · Gesetz für die θ-Singularität bei geraden Rängen

## Kurz gesagt (einfache Sprache)

Bei den geraden Rechenarten (Rang 6, 8, 10, …) hat die Kneser-Welle θ pro Periode genau eine scharfe Stelle. Vermutet war: Sie liegt dort, wo H_r den Fixpunkt zwei Stufen tiefer trifft, und ihre Schärfe ist α_r = ln λ_{r−4} / ln λ_{r−2}.
Ich habe das bei Rang 8 und neu bei Rang 10 nachgerechnet. Der Ort stimmt. Der Exponent stimmt auf 4 Stellen, sobald man die Korrekturglieder mitnimmt, die die Herleitung vorhersagt: Rang 8 gemessen 0.7302 (vorhergesagt 0.73014), Rang 10 gemessen 0.8360 (vorhergesagt 0.83607).
Den Kern der Herleitung habe ich als eigenes Lemma herausgezogen und mit dem Motor auf 6–7 Stellen bestätigt, samt der exakten Konstante. Das Lemma sagt: Die Umkehrung zwei Stufen tiefer läuft für V → −∞ mit der Rate λ_{r−4}^V gegen einen festen Punkt. Damit ist die Vermutung für alle geraden Ränge ab 6 hergeleitet. Offen bleibt dieselbe Zweigwahl-Lücke wie schon bei Rang 4 und 6.

Rechenbasis überall b = 2. Motor: gp-Fork, dps 20, Rang 8/9/10 mit 29.9 Stellen laut `hyperdigits`.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| G1 | **Rang 10, Ort:** θ₁₀ (Kneser-Dekation des Motors am abstossenden Fixpunkt kp₁₀ = −6.446551033773 + 1.135522857734i von H₉, μ = 6.846572362362 + 3.728976625797i, \|μ\| = 7.7962) ist genau bei x₁₀* + ℤ singulär. Dabei ist x₁₀* = H₁₀⁻¹(p₈) = −6.693605066320422, und p₈ = −5.683035382152679 ist der reelle Fixpunkt von H₈. | numerisch | arg(a_n/a_{n−1}) = −1.99724 (n=2), −1.93027 (8), −1.92666 (16), −1.925615 (31); Ziel −2π·frac(x₁₀*) = −1.925136; der Abstand fällt etwa wie n^{−1.7}. Negative Moden ≤ 8.4·10⁻¹¹. H₁₀(x₁₀*) − p₈ = 4·10⁻²⁵ | neu (Rang 10 war ungerechnet) |
| G2 | **Kern-Lemma, Teil Zahlen:** Für s = 4, 6, 8 gilt H_s⁻¹(V) − z_s ≈ κ_s·λ_{s−2}^V für Re V → −∞. Dabei ist z_s = x_s* − 1 (mit x₄* = −1, also z₄ = −2), und die Konstanten sind κ₄ = 1/sexp′(−1), κ₆ = λ₄^{pcent}/A₆, κ₈ = hyc₇/A₈. | numerisch | Quotient aufeinanderfolgender Werte (V = −k + it): s=4: 2.000130 (k=12, Ziel 2); s=6: 5.1717970 (k=12, Ziel 5.1717973); s=8: 9.493036 (k=10, Ziel 9.493087). Konstante: 0.779424 gegen 0.779373; 29.3046705 gegen 29.3046709; 11330.710 gegen 11330.717. Die Abweichungen fallen um den Faktor λ_{s−2} pro Schritt, wie das O(λ^{2V})-Glied verlangt | neu als eigenständige Aussage |
| G3 | **Kern-Lemma, Beweis:** Ist H_s(x_s*) = p_{s−2} mit A_s = H_s′(x_s*) ≠ 0 und gilt nahe p_{s−2}: H_{s−1}⁻¹(w) = log(σ_{s−2}(w − p_{s−2})/c_{s−1})/ln λ_{s−2} mit σ′(0) = 1, dann gibt es ein Φ_s, analytisch bei 0, mit H_s⁻¹(V) = z_s + Φ_s(λ_{s−2}^V) und Φ_s′(0) = c_{s−1}/A_s. | bewiesen (lokal, auf dem Zweig, auf dem die Funktionalgleichung gilt) | Abschnitt 1 | neu |
| G4 | **Struktursatz für gerade r ≥ 6:** Nahe x_r* gilt θ_r(z) = −z + Θ_r(u) mit u = λ_{r−4}^{H_r(z−1)} = (z − x_r*)^{α_r}·g_r(z). Dabei sind Θ_r und g_r analytisch, g_r(x_r*) = (A_r/c_{r−1})^{α_r} ≠ 0 und α_r = ln λ_{r−4}/ln λ_{r−2}. Insbesondere ist θ* = 2 − x_r* + h_r(ψ_r) mit ψ_r = 1 + H_{r−1}⁻¹(z_{r−2}), und die singulären Exponenten sind genau {kα_r + m : k ≥ 1, m ≥ 0}. | hergeleitet mit Lücke (Zweigwahl: h_r(w) = 1 + h_r(H_{r−1}⁻¹ w) längs des Wegs; wie bei kn.satz.theta_singularitaet) | Abschnitt 2; r = 6 stimmt mit θ*, K aus r18/02 B3 überein | neu (verallgemeinert r18/02 B3) |
| G5 | **Koeffizienten mit Korrekturen:** a_n·e^{2πinx*}·n^{1+α} = c₁ + c₂n^{−α} + c₃n^{−1} + c₄n^{−2α} + …, c₁ = KΓ-Faktor. Wegen n^{−α} konvergiert der rohe lokale Exponent sehr langsam: Rang 10 bei n = 30 nur 1.8275 statt 1.8361. | hergeleitet (aus G4 mit dem Übertragungssatz) / numerisch | siehe G6 | neu |
| G6 | **Exponent Rang 8 und 10:** Der Fit mit den Gliedern aus G5 gibt α₈ = 0.7302 (Ziel 0.730137) und α₁₀ = 0.8360 (Ziel 0.836070). Gegenprobe Rang 10: Mit dem vorhergesagten α ist der relative Rest 7.8·10⁻⁶ (n ≥ 4), mit α = 0.73 sind es 6.1·10⁻⁴, mit α = 1 sind es 2.4·10⁻². | numerisch | `p5_fit.py`; Rang 8 aus den r18/02-Daten (nur \|a_n\|, n ≤ 24), Rang 10 aus den neuen komplexen a_n (n ≤ 31). Je drei Startwerte nmin (3, 5, 8) und drei Modelltiefen, stabil auf ±0.0002 | neu; r18 hatte bei Rang 8 nur 1.713 |
| G7 | **Folge der Exponenten** (Basis 2): α₄ = 0, α₆ = 0.42182244, α₈ = 0.73013714, α₁₀ = 0.83606974. Dazu ln λ₂, ln λ₄, ln λ₆, ln λ₈ = 0.693147, 1.643220, 2.250564, 2.691837 (λ₈ = 14.75876993730426, neu). | numerisch | `p2_lemma.py` | λ₈, α₁₀ neu |
| G8 | **Annäherungen bei grossem Rang:** x_r* − p_{r−2} = −1.02111, −1.01446, −1.01057 (r = 6, 8, 10) und A_r = H_r′(x_r*) = 1.03250, 1.00523, 0.99838. Der Abstand z_{r−2} − p_{r−2} = x_{r−2}* − 1 − p_{r−2} ist −0.06561 (r = 8) und −0.03084 (r = 10). **Vermutung:** α_r → 1 und z_{r−2} − p_{r−2} → 0. Dann rückt der Verzweigungspunkt von H_{r−1}⁻¹ an ψ_r heran, und der Konvergenzradius von Θ_r schrumpft. | numerisch (3 Werte) / Vermutung | `p2_lemma.py`, `p3_rang10.py` | neu |
| G9 | **Negativ:** Der rohe lokale Exponent (Quotient benachbarter \|a_n\|) taugt bei hohen Rängen nicht zum Messen von α. Bei n ≤ 31 liegt er um 0.01–0.015 zu tief, weil das Korrekturglied n^{−α} mit α nahe 1 kaum abfällt. | negativ | G5/G6 | – |

## 1. Das Kern-Lemma (G2, G3)

Bezeichnungen der Motor-Familie: Ungerade Ränge sind regulär am reellen Fixpunkt p_{r−1} von H_{r−1}, mit λ_{r−1} := H_{r−1}′(p_{r−1}) > 1. Gerade Ränge sind Kneser am abstossenden komplexen Fixpunkt von H_{r−1}. Für ungerade q gilt nahe p_{q−1}

H_q⁻¹(w) = log(σ_{q−1}(w − p_{q−1})/c_q)/ln λ_{q−1}, mit σ′(0) = 1 und c_q = `hyc[q]` (c₃ = 1, H₃⁻¹ = log_b; c₅ = λ₄^{pcent}, weil der Motor ipent = log_{λ₄}σ − pcent schreibt; c₇ = `hyc[7]` = 11389.99447).

**Beweis von G3.** Für gerades s gilt H_s(z) = H_{s−1}⁻¹(H_s(z+1)) und H_s(z_s + 1) = H_s(x_s*) = p_{s−2}. Also ist

λ_{s−2}^{H_s(z)} = σ_{s−2}(H_s(z+1) − p_{s−2})/c_{s−1} =: E(z).

E ist bei z_s analytisch und hat dort eine einfache Nullstelle mit E′(z_s) = A_s/c_{s−1} ≠ 0. Also ist E lokal umkehrbar, z = z_s + E⁻¹(u). Setzt man u = λ_{s−2}^V, folgt H_s⁻¹(V) = z_s + Φ_s(λ_{s−2}^V) mit Φ_s = E⁻¹, Φ_s′(0) = c_{s−1}/A_s. ∎

Gültig ist das auf dem Zweig von H_s⁻¹, der die Funktionalgleichung erfüllt, für V mit beschränktem Imaginärteil und Re V → −∞. Der Motor zeigt genau diese Konstanten (G2). s = 4 ist derselbe Satz mit H₃⁻¹ = log_b, p₂ = 0, c₃ = 1, λ₂ = b und A₄ = sexp′(−1). Für s = 6 schreibt der Motor die Normierung als pent(z) = Ψ(λ₄^{z+pcent}); daher kommt κ₆ = λ₄^{pcent}/A₆. Für s = 8 ist κ₈ = hyc₇/A₈ = 11389.99447/1.00523157 = 11330.7170.

## 2. Herleitung des Gesetzes (G4)

Sei r ≥ 6 gerade, h = h_r die Abel-Funktion von F = H_{r−1} am Kneser-Fixpunkt, also θ_r(z) = h(H_r z) − z. Zweimal die Abel-Gleichung h(w) = 1 + h(F⁻¹ w) und einmal F⁻¹(V) = 1 + F⁻¹(H_{r−2}⁻¹ V) ergeben

θ_r(z) = 2 − z + h(1 + F⁻¹(H_{r−2}⁻¹(V))), mit V = F⁻¹(H_r z) = H_r(z−1).

Für z → x_r* (von oben) gilt H_r z → p_{r−2}, also V = log(σ(H_r z − p_{r−2})/c_{r−1})/ln λ_{r−2} → −∞ mit Im V → π/(2 ln λ_{r−2}). Nach dem Kern-Lemma (s = r−2) ist H_{r−2}⁻¹(V) = z_{r−2} + Φ(u), mit

u = λ_{r−4}^V = [σ(H_r z − p_{r−2})/c_{r−1}]^{ln λ_{r−4}/ln λ_{r−2}} = (z − x_r*)^{α_r}·g_r(z).

Also ist θ_r(z) = −z + Θ_r(u) mit Θ_r(u) = 2 + h(1 + F⁻¹(z_{r−2} + Φ(u))). Θ_r ist bei u = 0 analytisch, wenn (i) F⁻¹ bei z_{r−2} analytisch ist und (ii) h bei ψ_r = 1 + F⁻¹(z_{r−2}) analytisch ist. (i) gilt, weil z_{r−2} ≠ p_{r−2} (Abstand −0.066 bzw. −0.031). Da z_{r−2} < p_{r−2}, ist F⁻¹(z_{r−2}) komplex mit Im = π/ln λ_{r−2}; für r = 6 ist das psing des Motors, Im = 1.91185 = π/ln λ₄. (ii) folgt aus h′ = 1/F′(h) ≠ 0 wie in r18/02. Der lineare Koeffizient ist

K = h′(ψ_r)·(F⁻¹)′(z_{r−2})·(c_{r−3}/A_{r−2})·(A_r/c_{r−1})^{α_r}.

Für r = 6 ist das genau das K aus r18/02 B3 (c₃ = 1, A₄ = sexp′(−1), λ₄^{−pcent·α} = b^{−pcent}). Daraus folgen die singulären Exponenten kα_r + m. Mit dem Übertragungssatz (Flajolet–Odlyzko, angewandt auf (1 − Q e^{…})^{kα+m}) folgt die Reihe in G5: Das Glied u² liefert n^{−1−2α}, das Glied (z − x*)·u liefert n^{−2−α}. Bei Rang 4 hat H₂⁻¹(V) = V/b keinen endlichen Grenzwert, die Kette bricht nie ab, und es entsteht der superlogarithmische Fall (α₄ = 0).

**Lücke (wie bei kn.satz.theta_singularitaet und r18/02 B3).** Dass die Abel-Gleichung und die Hauptzweige von F⁻¹ und H_{r−2}⁻¹ längs H_r(x* + iε) die richtige analytische Fortsetzung sind, ist nicht bewiesen. Es wird nur dadurch gestützt, dass Ort, Exponent und Konstanten numerisch passen. Ebenfalls nur numerisch ist, dass es auf ℝ pro Periode **keine weitere** Singularität gibt; das stützt die Phasenkonvergenz (G1, r18/02 B2/B6).

## 3. Zahlen zu Rang 8 und 10 (G1, G6)

| | Rang 8 | Rang 10 |
|---|---|---|
| p_{r−2} | −3.69940781692481576 | −5.68303538215267880 |
| λ_{r−2} | 9.49308711431589 | 14.7587699373042601 |
| x_r* | −4.71387560712689689 | −6.69360506632042193 |
| A_r = H_r′(x_r*) | 1.00523157491 | 0.998384186281 |
| kp_r | −4.294771722864 + 1.283451257072i | −6.446551033773 + 1.135522857734i |
| α_r vorhergesagt | 0.73013714 | 0.83606974 |
| α_r Fit (G5-Modell) | 0.7302 (\|a_n\|, n ≤ 24) | 0.8360 (komplex), 0.8360 (\|a_n\|) |
| roher lokaler Exponent − 1 | 0.7164 (n = 23) | 0.8275 (n = 30) |
| c₁ (Rang 10, komplex) | – | 0.00214 − 0.02680i |

Rang 10, θ-Lauf: dps 20, y = 0.03, N = 64, Taylor-Kreis um kp₁₀ mit ρ = 0.03. Mit ρ = 0.08 scheitert er: Der Kreis trifft eine Stelle, an der der Rang-8-Newton aus dem Bogen läuft («knhex: Newton bleibt nicht im Bogen»). Grösster Schritt von θ auf der Linie 0.0047, also eindeutige Zweigverfolgung. Die Gegenprobe mit einer zweiten Einstellung (andere Linie y oder dps 30) fehlt für Rang 10 aus Zeitgründen. Die Phasenkonvergenz und der Fit-Rest von 8·10⁻⁶ sprechen aber gegen Rechenfehler auf dem Niveau von 10⁻⁴.

## Vorschläge für Grapheinträge

- `kn.satz.kernlemma_umkehr_minus_unendlich`: H_s⁻¹(V) = z_s + Φ(λ_{s−2}^V), Φ′(0) = c_{s−1}/A_s, für gerade s ≥ 4 | satz | bewiesen (lokal), numerisch für s = 4, 6, 8
- `kn.satz.theta_gerade_raenge_struktur`: θ_r = −z + Θ_r((z−x_r*)^{α_r} g_r(z)), singuläre Exponenten kα_r + m | satz | hergeleitet mit Lücke (Zweigwahl)
- Hochstufung `kn.verm.theta_gerade_raenge` → hergeleitet mit Lücke; Nachprüfung Rang 8 (α = 0.7302) und Rang 10 (Ort, α = 0.8360) eintragen | ergänzung | numerisch
- `kn.beob.theta10_koeffizienten`: Rang-10-Daten (x₁₀*, kp₁₀, μ, c₁, λ₈) | beobachtung | numerisch
- `kn.verm.alpha_r_gegen_1`: α_r → 1 und z_{r−2} − p_{r−2} → 0 für r → ∞ | vermutung | 3 Werte

## Offene Fragen

1. Zweigwahl-Lücke schliessen (gemeinsam mit Rang 4 und 6). Vielleicht genügt eine Intervallprüfung an einem Punkt plus Eindeutigkeit der Fortsetzung.
2. Rang 10 mit einer zweiten Einstellung nachrechnen (dps 30 oder y = 0.05) und K bzw. c₁ aus der Formel in Abschnitt 2 vorhersagen und vergleichen.
3. Rang 12 (braucht Rang 11/12, kostet wohl 3–5 min): α₁₂ = ln λ₈/ln λ₁₀ voraussagen und prüfen.
4. Lässt sich G8 beweisen, also H_r ≈ z + 1 auf einem wachsenden Bereich links und damit α_r → 1? Was bedeutet das Zusammenrücken von z_{r−2} und p_{r−2} für die Grösse von θ_r?
5. Andere Basen: Das Gesetz sollte über den Schwellen η_r unverändert gelten, mit λ₂ = b.

## Reproduktion (aus diesem Ordner, `C:/Users/janis/miniconda3/python.exe`, `PYTHONIOENCODING=utf-8`)

```
python p2_lemma.py 20 rang9        # G2, G7, Konstanten Rang 8/9 (~1.5 min) -> daten/lemma_20.json
python p3_rang10.py 20             # x10*, kp10, A10 (~2 min) -> daten/rang10_20.json
python p4_theta_r.py 20 0.03 64 0.03 10   # theta_10-Koeffizienten (~2 min) -> daten/theta10_20_0.03_64.json
python p5_fit.py                   # G5, G6 (Sekunden; liest auch r18/02/daten/theta8_20_0.03_64.json)
```

`motor.py` ist die Schnittstelle aus r18/02 mit `state_cache=False`. Der Zwischenspeicher liegt nur in `cache/`. `lib_r.py` ist eine unveränderte Kopie aus r18/02, `p4_theta_r.py` eine auf Rang R umgestellte Kopie von r18/02 `s6_theta8.py`. `p1_konstanten.py` ist ein früher Versuch: Er überschritt die Zeitgrenze und lieferte keine Daten, `p2`/`p3` ersetzen ihn. Es laufen keine Hintergrundprozesse, es gibt keine Commits.
