# 01 Halbebenen-Eindeutigkeit (Q*), Runde 25.10.2026

**Kurz und einfach.** Die Frage Q* ist weiter offen: Ist jede reelle Tetration, die auf Re z > −2 holomorph ist, schon Kneser?
Ich habe genau herausgearbeitet, welche Bedingung bei i∞ fehlt. Es fehlt nur, dass die periodische Abweichung p_T von Kneser keinen
Fourier-Anteil e^{−2πinz} (n ≥ 1) hat, also langsamer wächst als e^{2π Im z}. Diese Schwelle ist scharf. Ausserdem gilt: Eine Verbiegung
von Kneser, die die obere Halbebene in sich abbildet, ist immer trivial (Schwarz-Lemma). Ein Gegenbeispiel müsste also auf Nebenblätter von
sexp_K ausweichen. Für alle 48 geprüften Verbiegungen mit trigonometrischem Polynom liegt die Singularität schon auf dem Hauptblatt.

Zeit: 13:22–13:44 (`date`). Basis e, H = {Re z > −2}, U = {Im z > 0}, F = sexp_K (Kneser), q = e^{2πiz}, L = 0.31813… + 1.33724…i,
a = Re L, b = Im L = Im log L. Jede reelle Tetration T auf (−2, ∞) schreibt sich als T = F∘σ, σ = id + p, p = slog_K∘T − id reell-analytisch
und 1-periodisch (slog_K ist auf ganz ℝ analytisch, also braucht es keine Monotonie). Hol(H) = Menge der reellen Tetrationen, die auf H holomorph sind.

## Befunde

| # | Aussage | Status | Beleg/Genauigkeit | Neuheit-Vermutung |
|---|---|---|---|---|
| 1 | **Satz A (scharfes Liouville-Kriterium):** T ∈ Hol(H) ist Kneser, sobald p_T sich eindeutig holomorph auf U fortsetzt und einseitig −Im p_T(z) ≤ o(e^{2π Im z}) gilt, gleichmässig in Re z (statt −Im geht auch ±Re oder +Im) | bewiesen | Laurent in q und Borel–Carathéodory, unten | Formulierung neu, Kern klassisch |
| 2 | Die Schwelle e^{2π Im z} ist scharf: p = ε sin 2πz erreicht genau \|p\| ~ (ε/2)e^{2π Im z} | bewiesen | direkte Rechnung | – |
| 3 | **Satz A′ (Schröder-Form, Abschwächung von Paulsen–Cowgill):** Liegt T(z) für Re z ∈ [c, c+1+δ] (c ∈ (−2,−1]), Im z ≥ Y in einer kleinen Kreisscheibe B um L (Ψ_L dort umkehrbar), dann gilt T → L, also T = Kneser | bewiesen (mit PC) | E(q) = Ψ_L(T(z))·L^{−z} wächst höchstens wie \|q\|^{−b/2π}, b/2π = 0.21282 < 1 | vermutlich neu |
| 4 | **Klassifikation der Enden:** Ist E_T eindeutig nahe q = 0, so gilt: E_T hebbar ⇔ T → L. Eine Nullstelle der Ordnung k ↔ ζ ~ e^{z(log L+2πik)}, ein Pol der Ordnung k ↔ ζ ~ e^{z(log L−2πik)} (T hat keinen Grenzwert) | hergeleitet mit Lücke | Lücke: E_T eindeutig verlangt, dass T nahe i∞ die singulären Werte {0, 1, e, e^e, …} von Ψ_L umgeht | neu |
| 5 | Es gibt kein T ∈ Hol(H) mit Ende „→ L mit k ≥ 1 Zusatzwindungen“ (ζ ~ e^{z(log L+2πik)}) | bewiesen (aus PC) | T → L, Spiegelung gibt → L̄, PC ⇒ Kneser, das aber k = 0 hat | neu |
| 6 | Hätte T ∈ Hol(H) den Grenzwert L̄ bei +i∞ („vertauschtes Kneser“), so müsste der Schröder-Exponent bei L̄ ein Nebenzweig log L̄ + 2πik mit k ≥ 1 sein | bewiesen (bei eindeutiger Schröder-Koordinate) | \|L̄^z\| = e^{ax+by} wächst, also E(0) = 0 | neu; PC schliesst diesen Fall **nicht** aus |
| 7 | **Satz B (Schwarz-Starrheit):** Bildet σ = id + p (p reell, 1-periodisch, holomorph auf U) die Halbebene U in U ab, so ist p ≡ 0 | bewiesen | g(q) = q·e^{2πip} ist eine Selbstabbildung der Kreisscheibe mit \|g\| = 1 am Rand | neu (für diese Frage) |
| 8 | **Lemma C:** Ist p nicht konstant, so reicht jede Zusammenhangskomponente von {Im σ < 0} ∩ U bis i∞, und σ ist auf ihrem Rand nach unten und oben unbeschränkt. Die an ℝ anliegende Komponente von {Im σ > 0} ist nach oben unbeschränkt | bewiesen | Maximumprinzip (ein Punkt ∞ ist polar) bzw. eigentliche Abbildung auf die Scheibe | neu |
| 9 | Für p = 0.1 sin 2πz, 0.01 sin 2πz, 0.1(cos 2πz−1), 0.03 sin 4πz und eine Mischung hat F∘σ eine **Log-Singularität auf dem Hauptblatt** in H, z.B. z* = −1.34447550088905 + 0.43879214577273i für ε = 0.1 | numerisch | \|σ(z*)+2\| < 10^{−40}; senkrechter Weg mit 4000 bzw. 12000 Schritten, dps 30 bzw. 40, kein Kreuzen von (−∞,−2] | macht die Beobachtung aus Runde 18 zu einem Nachweis auf dem Hauptblatt |
| 10 | Vermutung: Jede Verbiegung F∘(id+p) mit nicht konstantem reellem trigonometrischem Polynom p (σ′ > 0 auf ℝ) hat eine Hauptblatt-Singularität in H | Vermutung | 43 von 43 zufälligen monotonen p bis Grad 4 (seed 7), dazu 5 Beispiele; Gegenbeispiel-Suche ohne Fund | neu |
| 11 | Auf dem Blatt von F, das man nach einmaligem Kreuzen von (−4,−3) erreicht, ist −2 **regulär**: F = Log(Log F(w+2) ± 2πi) → log 2π ± iπ/2 | bewiesen (lokale Rechnung) | – | erklärt, warum das Zählen der Lösungen von σ = −2 allein nichts beweist |
| 12 | Q* selbst: weder Beweis noch Gegenbeispiel | negativ | Die Lücke sitzt genau in Befund 4 (Eindeutigkeit von E_T) und bei Enden ohne Grenzwert oder mit Grenzwert L̄ bzw. anderen Fixpunkten | – |

## Beweise

**Satz A.** Weil p_T auf U eindeutig und 1-periodisch ist, ist P(q) := p_T(z) holomorph auf 0 < |q| < 1. Wegen der reellen Analytizität ist P auch
auf einem Ring um |q| = 1 holomorph, und die Spiegelung P(1/q̄) = conj P(q) setzt P auf ganz ℂ* fort: P = Σ c_n q^n mit c_{−n} = conj c_n.
Sei f = iP mit Koeffizienten a_n = i c_n und A(r) = max(0, sup_{|q|=r} Re f). Für n ≠ 0 gilt
(1/π)∫ Re f(re^{it}) e^{−int} dt = a_n r^n + conj(a_{−n}) r^{−n}, und |Re f| ≤ 2A − Re f. Daraus folgt |a_{−m}| ≤ r^m(4A(r) + 2|a_0|) + |a_m| r^{2m}.
Mit A(r) = o(1/r) (das ist −Im p ≤ o(e^{2π Im z}), denn |q| = e^{−2π Im z}) geht die rechte Seite für r → 0 gegen 0. Also ist c_{−m} = 0 für m ≥ 1,
wegen der Symmetrie auch c_m = 0, und p_T ist konstant. Aus T(0) = 1 folgt p_T = 0. ∎
*Was das für Q* heisst:* Ein Gegenbeispiel T ∈ Hol(H) hat entweder ein p_T, das auf U mehrdeutig ist (T trifft dort Singularitäten der Riemannschen
Fläche von slog_K), oder es gilt |p_T| ≥ c·e^{2π Im z} auf einer Folge. Die fehlende Bedingung bei i∞ ist also genau: „p_T ist eindeutig und o(e^{2π Im z})“.
Das ist viel schwächer als der Grenzwert L bei Paulsen–Cowgill.

**Satz A′.** Sei Ψ die Schröder-Funktion bei L (Ψ(e^w) = LΨ(w), Ψ(L) = 0, univalent auf B). Für eine kleine Scheibe B ist {w ∈ B : e^w ∈ B} zusammenhängend,
also gilt dort die Schröder-Gleichung. Dann ist ζ(z) = Ψ(T(z)) auf dem Halbstreifen V holomorph mit ζ(z+1) = Lζ(z), und E(q) = ζ(z)L^{−z} ist eindeutig
auf einer punktierten Scheibe um q = 0. Es gilt |E| ≤ C·e^{−ax+by} ≤ C′|q|^{−b/2π} mit b/2π = 0.21282 < 1, also ist q = 0 hebbar.
Damit gilt ζ → 0, also T → L gleichmässig auf V. Mit T(z+1) = e^{T(z)} folgt der Grenzwert L für jedes x ≥ c. Für x ∈ [c−1, c] gilt e^{T(x+iy)} → L
gleichmässig, also liegt T(x+iy) schliesslich nahe der diskreten Menge L + 2πiℤ. Wegen des Zusammenhangs ist der Index konstant, und bei x = c ist er 0.
Wegen c − 1 ≤ −2 deckt das jedes x > −2 ab (deshalb c ≤ −1).
Spiegelung gibt L̄ bei −i∞, und Paulsen–Cowgill (kn.satz.paulsen_cowgill) liefert T = Kneser. ∎
(Für c > −1 bekäme man links nur T(x+iy) → L + 2πim. Deshalb muss der Streifen am linken Rand beginnen.)

**Befund 4/6 (Enden).** Ist E eindeutig und hebbar mit Nullstelle der Ordnung k ≥ 0, so gilt ζ = e^{z(log L+2πik)}E_1(q) mit E_1(0) ≠ 0.
Ein Pol der Ordnung k gibt e^{z(log L−2πik)}, und wegen 2πk > b geht dann |ζ| → ∞. Beim Ende an L̄ mit Ψ̄ ist |L̄^z| = e^{ax+by}.
Beschränktes ζ → 0 erzwingt dann E(0) = 0, also den Exponenten log L̄ + 2πik mit k ≥ 1. Kandidaten für ein Gegenbeispiel zu Q* sind deshalb nur
(i) Enden ohne Grenzwert (E mit Pol oder wesentlicher Singularität), (ii) Enden an L̄ mit Nebenzweig k ≥ 1, (iii) Enden an anderen Fixpunkten L_j von exp,
(iv) Fälle, in denen T nahe i∞ die singulären Werte {0, 1, e, …} von Ψ umläuft. Die Werte von Ψ_L (= P^{−1}, P ganz, P(Lw) = e^{P(w)}) sind singulär genau dort:
P lässt 0 aus, hat keine kritischen Punkte (P′(Lw) = e^{P(w)}P′(w)/L, also folgt P′ ≠ 0 aus P′(0) = 1), und seine asymptotischen Werte enthalten die Bahn 0, 1, e, e^e, ….
Dass es *nur* diese sind, habe ich nicht bewiesen.

**Satz B.** g(q) = e^{2πiσ(z)} = q·e^{2πip} ist holomorph auf D*. Wegen σ(U) ⊂ U gilt |g| < 1, also ist 0 hebbar. Die Windungszahl von g auf kleinen Kreisen ist 1,
weil p eindeutig ist. Daher ist g(0) = 0 eine einfache Nullstelle, und g hat in D keine weitere. Da σ auf ℝ reell ist, gilt |g| = 1 auf |q| = 1 (g ist dort holomorph).
Also ist g/q nullstellenfrei mit Betrag 1 am Rand. Maximum- und Minimumprinzip geben g = e^{iα}q, also p konstant, also p = 0. ∎
*Folgerung:* Jede nicht triviale Verbiegung F∘σ verlässt mit σ die obere Halbebene. Solange σ dabei (−2, ∞) kreuzt, bleibt man auf dem Hauptblatt.

**Lemma C.** (a) Sei W_0 eine Komponente von {Im σ < 0} ∩ U, und angenommen σ ≥ m auf ∂W_0. Dann bildet u = σ − m + 1 die Menge W_0 in die untere Halbebene ab
und den Rand nach [1, ∞). Die Funktion Im log u liegt in (−π, 0), ist harmonisch und beschränkt und hat am endlichen Rand den Wert 0. Der Randpunkt ∞ ist polar,
also folgt Im log u ≡ 0, ein Widerspruch. Ebenso nach oben. Eine beschränkte Komponente (ohne i∞) widerspricht dem Maximumprinzip für |g|.
(b) Die an ℝ anliegende Komponente C_0 von {Im σ > 0}: Wäre sie nach oben beschränkt, so wäre g: Ĉ_0 → D eigentlich, also surjektiv, und träfe 0.
Aber g ≠ 0 auf D*. ∎

## Rechnungen

`hauptblatt.py` (nur mpmath, kein Kneser-Wert nötig): Newton für σ(z) = −2 aus Zufallsstarts in (−2, 0) × (0.05, 2.5). Für jede Wurzel prüft das Skript den senkrechten
Weg von Re z* nach z*: Wechselt Im σ das Vorzeichen, muss dort Re σ > −2 sein. Dann bleibt man auf dem Hauptblatt, und F∘σ wird bei z* unbeschränkt (σ′(z*) ≠ 0, Log-Punkt von F bei −2).
- `python hauptblatt.py 0` gibt die 5 Beispiele, `log_hauptblatt.txt`, etwa 2 min.
- Gegenprobe mit dps 40 und 12000 Wegschritten: `log_gegenprobe.txt`, alle 5 bestätigt, |σ(z*)+2| ≤ 10^{−40}.
- `python hauptblatt.py 120 7` zieht 120 zufällige p bis Grad 4. Davon sind 43 monoton, und alle 43 haben einen Hauptblatt-Treffer. `log_zufall.txt`, etwa 2 min.

## Vorschläge für Grapheinträge

- br.satz.halbebene_liouville_scharf: „Q* gilt für jedes T, dessen p_T auf U eindeutig ist und o(e^{2π Im z}) wächst; die Schwelle ist scharf“ (satz, bewiesen)
- br.satz.schroeder_kugel_genuegt: „Bleibt T im oberen Halbstreifen nahe L, so ist T Kneser (Abschwächung von Paulsen–Cowgill)“ (satz, bewiesen, stützt sich auf kn.satz.paulsen_cowgill)
- br.satz.enden_klassifikation: „Mögliche Enden einer Tetration in Hol(H) bei i∞: Schröder-Exponent log L + 2πik bzw. log L̄ + 2πik; nur L, k = 0 ist realisiert“ (satz, hergeleitet mit Lücke)
- br.satz.verbiegung_schwarz: „Eine Verbiegung sexp_K∘(id+p), die U in U abbildet, ist trivial“ (satz, bewiesen)
- br.lemma.zungen_bis_unendlich: „Komponenten von {Im σ < 0} reichen bis i∞, σ ist auf ihrem Rand beidseitig unbeschränkt“ (lemma, bewiesen)
- br.beob.hauptblatt_singularitaet_trig: „Verbiegungen mit trigonometrischem Polynom haben eine Log-Singularität auf dem Hauptblatt in H“ (beobachtung bzw. vermutung, numerisch an 48 Fällen)
- br.satz.nebenblatt_minus2_regulaer: „Nach Kreuzen von (−4,−3) ist −2 ein regulärer Punkt von sexp_K“ (satz, bewiesen)

## Offene Fragen

1. Q* in der Form von Befund 4: Kann E_T bei q = 0 einen Pol oder eine wesentliche Singularität haben? Und kann T nahe i∞ die Werte 0, 1, e, … umlaufen?
2. Gibt es ein „vertauschtes Kneser“ (Grenzwert L̄ bei +i∞, Exponent log L̄ + 2πi) in Hol(H)? Paulsen–Cowgill schliesst das nicht aus.
   Das ist der natürlichste Kandidat für ein Gegenbeispiel: Kneser-Konstruktion mit einem anderen Zweig des Multiplikator-Logarithmus.
3. Beweis von Vermutung 10 für alle nicht konstanten ganzen p. Fehlender Schritt: Eine Randkurve von σ^{−1}(ℂ∖(−∞,−2]) mit σ ∈ (−n−1, −n) muss bis zu einem
   ganzzahligen Wert laufen, bevor sie H verlässt oder nach i∞ geht.
4. Sind die singulären Werte von Ψ_L genau {0, 1, e, e^e, …}?
