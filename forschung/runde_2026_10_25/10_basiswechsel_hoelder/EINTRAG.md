# Eintrag Bericht 10 (Basiswechsel Rang 5: Hölder, Singularität) — 25.10.2026

Dateien: `daten/aussagen/r25_10_bw.yaml`, `daten/theorie/r25_10_bw.yaml`. `python baue.py --nur-pruefen`: 0 Fehler; 2 Warnungen „widerspruch“ zu `bw.verm.rang5_nirgends_diffbar` (gewollt: die alte Vermutung steht noch auf `vermutung`, siehe Änderungen am Bestand).

## Angelegte ids

| id | Titel kurz | Status |
|---|---|---|
| bw.lem.wackel | Wackel-Lemma: Differenzenquotient auf Skala s/U_n′ = Π_n·(1 + W(U_n,s)) | skizze |
| bw.satz.hoelder_eins | Hölder-Exponent überall 1, Modul h·e^{O(slog 1/h)} | skizze |
| bw.satz.kein_endlicher_ableitungswert | An jeder Stelle H′ ∈ {0, ∞} oder keine Ableitung; alle drei Mengen dicht (B3 + B7) | skizze |
| bw.kor.rang5_fast_ueberall_diffbar | Ĝ monoton ⇒ Φ₅ fast überall differenzierbar (Lebesgue) — **aus der Gegenprobe**, folgt aus B4 | bedingt |
| bw.satz.rang5_singulaer | Ĝ singulär (Cantor-artig), Φ₅′ = −1 fast überall; zwei Beweiswege | skizze |
| bw.kor.singulaer_mass | Ĝ′ = +∞ dĜ-f. ü., lokale Dimension überall 1 (B5) | skizze |
| bw.beob.singulaer_unsichtbar | Drift χ ≈ Var(P₄′)/2 = 2.04·10⁻⁷, sichtbar erst nach ~5·10⁶ Ebenen; 1.0657 ist ein Plateau (B6 + B9) | numerisch |
| bw.verm.dichotomie_raenge | glatt ⇔ Schälschritt mit konstanter Phase; singulär für alle r ≥ 5 (B8) | vermutung |

## Entscheidung: Runde 18 („Weierstrass-artig, nirgends differenzierbar“) gegen Bericht 10 („singulär, Cantor-artig“)

**Bericht 10 hat recht; die Vermutung der Runde 18 ist im Wortlaut widerlegt.** Die Entscheidung hängt nicht an den Lücken des Berichts:

1. **Lebesgue-Argument (neu in der Gegenprobe, lückenlos).** Ĝ(x) = lim (G(x+n) − n) mit G = spog_b∘P_d streng wachsend, also ist Ĝ monoton. Eine monotone Funktion ist nach Lebesgue
   fast überall differenzierbar. Also ist Φ₅ = Ĝ − id fast überall differenzierbar — „nirgends (mindestens fast nirgends) differenzierbar“ ist unmöglich. Gebraucht wird nur die Existenz von Φ₅,
   also genau die Voraussetzung, unter der die Vermutung überhaupt formuliert ist (bw.satz.phasengesetz_rang5). Eingetragen als `bw.kor.rang5_fast_ueberall_diffbar`.
2. **Wert der Ableitung.** Satz A des Berichts (an keiner Stelle eine Ableitung in (0, ∞)) braucht nur das Wackel-Lemma. Zusammen mit 1. folgt Ĝ′ = 0 fast überall — ein zweiter, kürzerer
   Beweis von Satz B ohne Martingal. Der Hellinger/Kakutani-Weg des Berichts ist ein unabhängiger dritter Weg (und liefert die Grössen je Ebene). Ein vierter, heuristischer: log Π_n ist eine Irrfahrt mit Drift −χ < 0 (Jensen), also Π_n → 0 fast überall, und nach dem Wackel-Lemma folgt dann H′ = 0.
3. **Wo die Runde 18 fehlging.** Ihr Bild war additiv: Φ₅′ ≈ Σ_j P₄′(U_j), eine Irrfahrt ohne Drift ⇒ „keine Ableitung“ (Weierstrass/Hardy, ab = 1). Exakt ist die Struktur
   multiplikativ: Ĝ′ = ∏(1 + P₄′(U_j)) mit lauter positiven Faktoren. Wegen log(1+x) = x − x²/2 + … hat log Π_n die Drift −χ ≈ −Var(P₄′)/2; sie überwiegt die Schwankung ab n ≈ 2/χ ≈ 10⁷ Ebenen,
   und das Produkt geht gegen 0, nicht in eine Oszillation. Der Skizzentext der Runde 18 („divergiert die Reihe wie eine Irrfahrt“) hatte also die richtige Zutat, aber den falschen Schluss: Divergenz von Σ log(1+P₄′) nach −∞ heisst Ableitung 0.
4. **Was von der Runde 18 bleibt:** Φ₅ ist nicht C¹ (sogar: nirgends ein Ableitungswert in (−1, ∞)), und der Modul h·slog(1/h) ist auf allen rechenbaren Skalen richtig.
   Befund #2 des Berichts („korrigiert R18-B11: kein Exponent < 1“) ist **keine** Korrektur: die Runde 18 hatte nie einen Exponenten < 1 behauptet; beide Fassungen sind Hölder-1.

Restunsicherheit: Satz A und damit der Wert Ĝ′ = 0 hängen am Wackel-Lemma (Lücke L2, Verzerrungsschranken, Routine) und an den bedingten Sätzen der Runde 18 (L1). Die **Widerlegung** von „nirgends differenzierbar“ hängt an keiner dieser Lücken.

## Gegenprobe (Pflicht)

| Befund | geprüft wie | Ergebnis |
|---|---|---|
| #1 Wackel-Lemma (numerisch) | `s2_wackeln.py exp(1) 2 3 90 0.03 0.08` — neue Stellen θ und neue Präzision (dps 90), Ausgabe `daten/s2_exp1_2_n3_dps90.json` | Rest nach Vorhersage 2.5·10⁻⁷ (U₃ = 4.5·10⁶) und 4.9·10⁻¹⁶ (U₃ = 7.3·10¹⁴), also 1.1/U₃ und 0.36/U₃ — vorhergesagte Ordnung O(1/U_n); Wackel-Amplitude 8.3·10⁻⁴ / 6.4·10⁻⁴ = W(U₃,s) auf 3·10⁻⁷. **Bestanden.** (Ein Versuch bei θ = 0.3 scheitert erwartungsgemäss: U₃ ≈ 10^(2·10³⁰⁰), h unter der Rechengenauigkeit.) |
| #4 Affinitäten | `s3_hellinger.py` erneut | 1 − A₂ = 1.0005·10⁻⁸, 1 − A₈ = 4.8377·10⁻⁸, Grenzwert 5.094487·10⁻⁸, M = 64 und 128 gleich — wie im Bericht. **Bestanden.** |
| #6 Drift χ | drittes Gitter (Mittelpunktregel, 12345 Punkte) und Parseval aus den Moden, eigene Rechnung | χ(e,2) = 2.0377869066·10⁻⁷ (Bericht 2.0377869·10⁻⁷); Var(P₄′)/2 = 2.0378192·10⁻⁷; P₄′ ∈ [−8.8736·10⁻⁴, 9.1963·10⁻⁴]; |a₁| = 7.1833·10⁻⁵ passt zu K₄·|ln log_e 2| = 1.96·10⁻⁴·0.3665 aus bw.beob.wellen_raenge. **Bestanden.** |
| #3 Satz A | Beweis gelesen | Schluss richtig: c₀ ≥ Spannweite/4, weil u+s über s ∈ [1,2] eine volle Periode überstreicht; Fehlerglieder O(s/U_n) und sup\|r\|/s → 0. Lücke L2 bleibt → `skizze`. |
| #4 Satz B | Beweis gelesen | Martingal- und Vergröberungsschritt richtig (Cauchy–Schwarz: Zusammenfassen erhöht die Affinität). Zweiter Beweis über Lebesgue + Satz A gefunden (oben). Status statt „bedingt“ `skizze` (L2 offen). |
| #5 | Beweis gelesen | Standard (de la Vallée Poussin für den singulären Anteil; Dimension aus #2); erbt die Lücken → `skizze` statt „hergeleitet“. |
| #7, #8 | gelesen | #7 in die Aussage zu Satz A eingearbeitet (Skizze); #8 als Vermutung. |
| #9 | gelesen | Fourier-Schranke geprüft: \|a_k\| ≤ 1/(2πk) folgt durch partielle Integration, weil dĜ ein Wahrscheinlichkeitsmass ist; in `bw.beob.singulaer_unsichtbar` eingearbeitet. |

Kein Befund fiel durch. Herabgestuft: #4 (bedingt → skizze), #5 (hergeleitet → skizze). Prozesse: nur eigene Läufe im Vordergrund (je < 10 s), keine Reste.

## Änderungen am Bestand (Vorschläge, nicht selbst ausgeführt)

1. **`bw.verm.rang5_nirgends_diffbar`** (`daten/aussagen/r18_01_bw.yaml`)
   - `status`: `vermutung` → `widerlegt`.
   - `bewertung.richtigkeit`: 'Vermutung mit zwei numerischen Belegen' → 'Widerlegt 25.10.2026: Ĝ = id + Φ₅ ist monoton, also nach Lebesgue fast überall differenzierbar (bw.kor.rang5_fast_ueberall_diffbar); richtig ist: singulär, Φ₅′ = −1 fast überall, nirgends ein endlicher positiver Ableitungswert von Ĝ (bw.satz.rang5_singulaer, bw.satz.kein_endlicher_ableitungswert).'
   - neues Feld `versionen`: `- {datum: '2026-10-25', aenderung: 'Status vermutung → widerlegt.', grund: 'Monotonie von Ĝ + Satz von Lebesgue; das additive Weierstrass-Bild Σ P₄′(U_j) übersieht die Drift −χ von Σ log(1 + P₄′(U_j)); Bericht forschung/runde_2026_10_25/10_basiswechsel_hoelder.'}`
   - `naechste_frage` → 'Beantwortet: Hardy-Lemma = bw.lem.wackel; Ergebnis bw.satz.rang5_singulaer.'
   - Kanten ergänzen: `{typ: ersetzt, ziel: …}` ist nicht nötig, die neuen Aussagen zeigen schon mit `widerlegt`/`korrigiert` auf sie. Titel und Aussage als Dokument der Runde 18 stehen lassen.
   - Theorie `daten/theorie/r18_01_bw.yaml`, Block `bw.verm.rang5_nirgends_diffbar`: an `genau` anhängen „Widerlegt (25.10.2026): fast überall differenzierbar, sogar singulär — [[bw.satz.rang5_singulaer]].“
2. **`bw.beob.skalensprung`**: in `aussage` „($\hat G'=1.0657272$)“ → „(Plateau $\Pi_n=1.0657272$ des Differenzenquotienten — kein Ableitungswert, vgl. bw.beob.singulaer_unsichtbar)“ und „($1.0665953$)“ entsprechend als nächstes Plateau $\Pi_{n+1}$; Kante `{typ: stuetzt, ziel: bw.verm.rang5_nirgends_diffbar}` → `{typ: stuetzt, ziel: bw.lem.wackel}`.
3. **`bw.satz.produktformel_rang5`**: `pruefung.f5_neue_frage` 'Existiert ∏(1 + Φ₄′(T_d^j y)) für fast alle y?' → 'Beantwortet: das Produkt geht für fast alle y gegen 0 (Drift −χ, bw.beob.singulaer_unsichtbar), also H′ = 0 f. ü. (bw.satz.rang5_singulaer).'; `bewertung.nutzen` „… der die Rang-5-Welle rau macht“ → „… der die Rang-5-Welle singulär macht“.
4. **`bw.satz.abschaelen_rang`** und **`bw.satz.phasengesetz_rang5`**: in `f5_neue_frage` den Verweis `bw.verm.rang5_nirgends_diffbar` → `bw.satz.rang5_singulaer`. Theorie `r18_01_bw.yaml`: die `einordnung`-Links „Glattheit: [[bw.verm.rang5_nirgends_diffbar]]“ (Blöcke phasengesetz_rang5, produktformel_rang5, skalensprung, anlauf_ursache) → [[bw.satz.rang5_singulaer]].
5. Übersetzungen (`daten/uebersetzung/bw.yaml`, `uebersetzung_struktur/idee.yaml`) der geänderten Felder nachziehen.

## Strang und Ideenkarte

- **Strang `bw.hoehere_raenge`** (bestehend): alle acht neuen ids als Mitglieder. Neuer `verlauf`: „Abschälen mit dem Rang-4-Wechsel, Phasengesetz für Rang 5 und (skizziert) alle Ränge; die tiefere Welle wird auf jeder Ebene multiplikativ eingespeist. Die Rang-5-Phase ist auf jeder Skala fast glatt (Hölder 1), aber singulär wie die Cantor-Funktion: Φ₅′ = −1 fast überall, nirgends ein endlicher positiver Ableitungswert von Ĝ — „nirgends differenzierbar“ ist widerlegt.“
- **Ideenkarte**: Unteridee `stufen.basiswechsel` (bestehend), `kurz` ändern: „Ab Rang 5 wird der Schälschritt selbst wellig; die tiefere Welle wird auf jeder Ebene eingespeist — keine Weierstrass-Rauheit, sondern eine singuläre, Cantor-artige Phase (Ableitung fast überall 0), die man auf keiner rechenbaren Skala sieht.“
