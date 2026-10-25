# Eintrag Bericht 14 „Andere Klammerungen“ (25.10.2026)

Dateien: `daten/aussagen/r25_14_br.yaml`, `daten/aussagen/r25_14_it.yaml`, `daten/theorie/r25_14_br.yaml`, `daten/theorie/r25_14_it.yaml`.
`python baue.py --nur-pruefen`: 0 Fehler, keine Warnung aus diesen Dateien.

**Wichtig beim Eintragen bemerkt:** $F(m)=me^m$ ist genau der Schritt der Produkt-Tetration zur Basis $e$ (br.def.produkt_tetration). Die Fatou-Reihe von Befund 4
steht mit denselben Koeffizienten schon in **br.beob.produkt_iteration_dreieck**, $F^{[1/2]}(1)=\mathrm{prodtet}_e(\tfrac12)=1.5134281085$ schon in br.beob.produkt_tetration_werte,
und das Residuum ½ ≠ ⅓ (Befund 9 für Rang 5′) schon in br.satz.produkt_parabolisch. Befund 4 ist darum **nicht** als eigene Aussage eingetragen; die Einträge verweisen auf den Bestand.

## Angelegte ids

| id | Titel kurz |
|---|---|
| br.satz.klammerung_log_reduktion | ausgewogene Tetration und linke Pentation = Iteration von m·eᵐ (= Produkt-Tetration) (bewiesen) |
| br.kor.linke_pentation_ausgewogen | x[5′]t = (x^x [4] 2^t)^{1/x} (bewiesen) |
| br.satz.klammerung_fortsetzung | gebrochene Höhen über die Fatou-Koordinate, ohne Kneser (Skizze) |
| br.beob.klammerung_werte | 2[4]3, 2[4]√2, 2[5′]½, 2[6′]½ … auf 22–30 Stellen (numerisch) |
| br.satz.linke_raenge_keim | Bauart der linken Ränge, Residuum ½ bzw. (9−2k)/(2x) (bedingt) |
| br.kor.linke_varianten_nicht_eta | keine linke Variante formal konjugiert zur η-Tetration oder eᶻ−1 (bedingt) |
| br.beob.ausgewogene_raenge_entartung | ausgewogene Pentation: dreifacher Fixpunkt bei 1 (numerisch) |
| it.verm.mexp_halbiterierte_divergent | formale Halbiterierte von m·eᵐ divergiert; C^∞, nicht analytisch bei 1 (Vermutung) |

## Gegenprobe

| Befund | geprüft wie | Ergebnis |
|---|---|---|
| 1, 2, 3 | Rechnungen nachvollzogen; Konvention $a_5b$ für $b=1,2$ eingesetzt | richtig, **bewiesen**. Befund 13 (Fixpunkte $2\pi ik$) elementar nachgerechnet und in br.satz.klammerung_log_reduktion eingearbeitet. |
| 4 Fatou-Reihe | eigener Bruchzahl-Code (unabhängig vom sympy-Skript) | β = ½, c₁…c₈ identisch. Aber schon im Bestand (s. oben) → nicht neu eingetragen. |
| 5 Fortsetzung | Argument gelesen | Konstante fällt in $\Phi^{-1}(\Phi+t)$ heraus; Eindeutigkeit nur zitiert → **Skizze** (wie Bericht „mit Lücke“). Auf $m>0$ ist es die Produkt-Tetration. |
| 6, 14 Werte | `kanonisch.py 30 300 12` neu gerechnet (`gegenprobe_kanonisch_A.txt`, identisch mit `kanonisch_A.txt`); **zweites Verfahren** (formale Halbiterierte bis Ordnung 14 nahe 0, $n=2000/4000$ Schritte mit W, dps 40; für $2[6']\tfrac12$ die Halbiterierte von $G$, $n=800/1500$) | $F^{[1/2]}(1)=1.51342810850016187455237845958$, $2[4]\sqrt2=2.55195159006647331599487105508$, $2[5']\tfrac12=3.30490800075364753129041404525$ (30 Stellen gleich), $2[6']\tfrac12=8.476232183507540378741322$ (25 Stellen gleich). Bestanden. |
| 7, 8 Keime | Induktion von Hand nachgerechnet (Iterationsformel für $n=2$, Rekursion $B_{k+1}=A_k+B_k$, $B_k=k-\tfrac92$); `rang_keime.py` neu gerechnet (identisch) | richtig für ganze Iterationszahlen; für Rang ≥ 6′ mit nicht ganzer Basis die benannte Lücke → **bedingt**. |
| 9 | Invarianz von β, $\eta^{e+w}$-Entwicklung ($b/a^2=\tfrac23$) nachgerechnet | richtig; Status wie 7 (bedingt), für Rang 5′ unbedingt. |
| 10, 11 | $g_2\dots g_8$ mit eigenem Code bestätigt; Folgerung geprüft (Eindeutigkeit formaler Halbiterierter, $x\ln x$ bei 1 umkehrbar) | Divergenz nur numerisch → **Vermutung** mit bewiesener Folgerung. |
| 12 | Zahlen in `kanonisch_A.txt`/`kanonisch_B.txt` nachgelesen | **Korrektur:** der Bericht nennt $(\ln P-\ell)/\ell^3=1.43932$ bei $\ell=0.005$ und das Mittel 1.44284; die Ausgabe zeigt **1.439155**, das Mittel ist **1.442762** (näher an $1/\ln2=1.442695$), Steigung −0.7214. Mit den richtigen Zahlen eingetragen (numerisch für Rang 5, formal für höhere). |
| 15 | — | nicht bearbeitet, nicht eingetragen. |

## Änderungen am Bestand (Vorschläge, nicht ausgeführt)

1. **br.beob.produkt_iteration_dreieck**, Feld `aussage`, nach der Fatou-Reihe ergänzen: „Alle Koeffizienten sind rational, weil $\ln F=\ln x+x$ exakt gilt: $\Phi(F)-\Phi=1$ wird zu $(1-e^{-x})/x+\beta x+\sum c_k(F^k-x^k)=1$ (22 Glieder in forschung/runde_2026_10_25/14_klammerungen/fatou_reihe_ausgabe.txt); $|c_k|(2\pi)^k/(k-1)!$ bleibt für $k\le22$ zwischen 0.03 und 5.4.“ Kante `verwandt` → it.verm.mexp_halbiterierte_divergent. Begründung: Befund 4, beim Eintragen unabhängig nachgerechnet.
2. **br.satz.produkt_parabolisch**, Feld `aussage`: „Die Reihe divergiert im Allgemeinen (it.satz.baker)“ → „Die formale Halbiterierte divergiert numerisch wie $(k-2)!/(2\pi)^k$ (it.verm.mexp_halbiterierte_divergent)“; ausserdem ergänzen: „Dieselbe Abbildung ist im Logarithmus die ausgewogene Tetration und die linke Pentation (br.satz.klammerung_log_reduktion).“ Kante `verwandt` → br.satz.klammerung_log_reduktion.
3. **br.satz.hyperop_varianten**, Feld `aussage` (ii), ergänzen: „Im Logarithmus ist das die Iteration von $me^m$, also der Produkt-Tetration; damit gibt es eine natürliche gebrochene Höhe ohne Kneser (br.satz.klammerung_fortsetzung).“ Feld `pruefung.f5_neue_frage` ergänzen um den Verweis auf br.satz.linke_raenge_keim.

## Stränge und Ideenkarte

- **Neuer Strang** `br.klammerungen` (gebiet br), Name „Andere Klammerungen — linke und ausgewogene Hyperoperationen“, Mitglieder: die sieben br-ids oben.
  `verlauf`: „Im Logarithmus sind ausgewogene Tetration und linke Pentation die Iteration von m·eᵐ (wie die Produkt-Tetration); der parabolische Fixpunkt gibt gebrochene Höhen ohne Kneser, und die linken Ränge unterscheiden sich nur im iterativen Residuum.“
- it.verm.mexp_halbiterierte_divergent in den bestehenden Strang **it.parabolisch** (gleiches Gebiet); `verlauf` ergänzen: „… auch die Halbiterierte von m·eᵐ divergiert (numerisch).“
- Ideenkarte: Unteridee **stufen.hyperop** um `straenge: [br.klammerungen]` ergänzen (bisher nur `aussagen`).
