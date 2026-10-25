# Eintrag Bericht 29 (Transzendenz und Irrationalitätsmasse) in den Graphen

Stand 25.10.2026 17:55. Dateien: `daten/aussagen/r25_29_ar.yaml`; Theorie: `daten/theorie/r25_29_ar.yaml`.
`python baue.py --nur-pruefen`: 0 Fehler in diesen Dateien (erwartete Hinweise; dazu der gewollte Hinweis „gr.satz.kehrwert_turmsumme wird korrigiert“, siehe unten).

## Vollständigkeit des Berichts (Lauf um 14:00 gestoppt)
BERICHT.md ist um 14:00:12 zuletzt geschrieben und formal vollständig (Befunde, Beweise, Vorschläge, Fragen, Dateien); die Zeile „Ende etwa 14:15“ ist eine Schätzung.
Alle vier Skripte liegen vor; nur `e1` hat eine Ergebnisdatei, `e2`–`e4` geben auf den Bildschirm aus. Alle vier habe ich neu laufen lassen (unten) — sie bestätigen die Zahlen.
Nicht fertig und darum **nicht** als eigene Aussage eingetragen: B8′ (T̃(½) und T̃′(0) gleichzeitig frei; Positivität der Kerne nicht ausgeschrieben) und B10′ (fast alle Basen; Analytizität in b ungeprüft) — beide nur als offene Punkte in den Prüffragen erwähnt.

## Angelegte ids
- ar.satz.kehrwertsumme_liouville | Σ 1/(b↑↑k) für ganze b ≥ 2 Liouville, expliziter Kettenbruch durch Faltung (B1) — bewiesen
- ar.satz.rang_dichotomie_kehrwerte | Kehrwertsummen: divergent (Rang ≤ 2), 1/(b−1) (Rang 3), Liouville (≥ 4) (B2) — bewiesen
- ar.satz.ableitungsorbit | T′(1) = e·T′(0) nicht beide algebraisch; algebraische / ganze Basis (B3–B5) — bewiesen
- ar.satz.schanuel_ableitungen | unter Schanuel höchstens ein algebraisches T′(n), trdeg ≥ n (B6) — bedingt
- ar.satz.jet_identitaeten | slog′(0) = 1/a₁, slog″(0) = (a₁² − 2a₂)/a₁³ (B7) — bewiesen
- ar.satz.keine_arith_info_ableitung | sexp′(0) beliebig für exotische Lösungen; Logik „SC ⇒ transzendent für alle“ ⇔ ¬SC (B8, B9) — bewiesen
- ar.satz.metrisch_fast_alle | fast alle Stellen: 2n unabhängige Zahlen, μ = 2, Chintschin/Lévy (B10) — bewiesen
- ar.satz.fixpunkt_transzendent | L, λ = L transzendent; L_b bei algebraischem b; λ_b unter Schanuel (B11) — bewiesen
- ar.beob.kf_population | 47 Kneser-Werte, 21 485 Glieder wie Zufallszahlen; Ableitungs-Konstanten (B12, B13) — numerisch
- ar.frage.kehrwertsumme_e | Ist C_e = e⁻¹ + e⁻ᵉ + … transzendent? (B15, Frage 1) — offen

Nicht eigens eingetragen: B14 (Literaturstand, steckt in ar.satz.keine_arith_info_ableitung und ar.beob.ausschluss), Fragen 2–4 (in den Prüffragen von ar.satz.fixpunkt_transzendent bzw. ar.satz.ableitungsorbit).

## Gegenprobe
| Befund | geprüft wie | Ergebnis |
|---|---|---|
| B1 | Beweis gelesen (gekürzter Nenner b^{e_N} wegen ≡ 1 mod b; Rest < 2b^{−e_{N+1}}); `e2_kehrwertsumme_liouville.py` neu gerechnet (exakte Brüche) | Faltungsregel stimmt für b = 2 bis N = 5, b = 3, 4, 5, 7 bis N = 3, b = 10 bis N = 2; C₂-Kettenbruch wie angegeben → bewiesen |
| B2 | H_r(b,k) = b↑↑M_k für r ≥ 4 nachgeprüft (auch b = 2, wo H_r(2,2) = 4 für alle r) | bestanden → bewiesen |
| B3–B5 | Quotienten nachgerechnet; Lücke gefunden: braucht T′(0) ≠ 0 (sonst T′(0) = T′(1) = 0 beide algebraisch) | als Voraussetzung ergänzt → bewiesen |
| B6 | Schanuel-Orbit mit α = 1, Monom-Argument, Körpergleichheit nachgeprüft | bestanden → bedingt (Schanuel) |
| B7 | Rechnung nachvollzogen; `e4_jet_identitaeten.py` neu gerechnet | alle Differenzen < 10⁻¹⁹⁹ → bewiesen, Neuheit folklore (slog′(0) = 1/T′(0) stand schon in gr.beob.kneser_ableitung) |
| B8 | Poisson-Kern-Konstruktion nachgeprüft; für γ < 1 muss r so gewählt werden, dass (1−r)/(1+r) < γ (im Bericht stillschweigend) | ergänzt → bewiesen |
| B9 | Logik: ∀T (SC → P(T)) ⇔ (SC → ∀T P(T)), und ∀T P(T) ist falsch wegen ar.satz.keine_arith_info | bestanden → bewiesen (mit B8 in einer Aussage) |
| B10 | Nullstellen-Argument und Rückzug von Nullmengen geprüft (braucht T′ > 0) | Voraussetzung T′ > 0 ergänzt → bewiesen |
| B11 | Hermite–Lindemann / Gelfond–Schneider; Schanuel-Fall selbst nachgerechnet (abhängiger Fall gibt L_b = q rational, λ = ln q) | bestanden; (iii) präziser formuliert |
| B12 | `e1_kf_statistik.py` neu gerechnet (11 s) | Ergebnisdatei bitgleich: χ² = 12.03, Chintschin 2.6643, Lévy 1.1800 (sd 0.037), μ_max 2.2308, grösstes Glied 154 488, Kontrollen e (780), √2 (2995) |
| B13 | `e3_ableitungen_kf.py` neu gerechnet | gleiche Anfänge, Gegenprobe gegen 160er-Stufe bestanden, μ_max ≤ 2.1483 |

## Änderungen am Bestand
Die „offen“-Stelle steht nicht in 07_arithmetik.yaml oder lit_ar_1.yaml, sondern in `gr.satz.kehrwert_turmsumme` (`daten/aussagen/00_grundlagen.yaml`). In 07_arithmetik.yaml sagt `ar.satz.ef_liouville` (f3) schon das Richtige, aber ohne „ganze b“ — der Graph widerspricht sich also derzeit.
1. `gr.satz.kehrwert_turmsumme`, Feld `pruefung.f4_uebertragbar`: alt „Auf jede rasch wachsende Folge (Engel-artige Reihen); ob C_b transzendent ist, ist offen.“ → neu
   „Auf jede rasch wachsende Folge (Engel-artige Reihen). Für ganze b ≥ 2 ist C_b eine Liouville-Zahl, also transzendent (ar.satz.kehrwertsumme_liouville); für b = e und nicht ganze b offen (ar.frage.kehrwertsumme_e).“
   Begründung: B1, Liouville-Beweis vollständig und nachgerechnet.
2. `gr.satz.kehrwert_turmsumme`, Feld `pruefung.f5_neue_frage`: alt „Ist $C_2$ irrational? (Glieder $2^{-2\uparrow\uparrow k}$: liouvilleartig schnell — vermutlich transzendent.)“ → neu
   „Ist $C_e$ transzendent, ist $C_b$ für rationales, nicht ganzes $b$ irrational? ($C_2$ ist eine Liouville-Zahl: ar.satz.kehrwertsumme_liouville.)“
3. `gr.satz.kehrwert_turmsumme`, neue Felder: `kanten: + {typ: verwandt, ziel: ar.satz.kehrwertsumme_liouville}` und
   `nachgeprueft: [{gegen: ar.satz.kehrwertsumme_liouville, ergebnis: 'Aussage und Beweis unberührt; nur f4/f5 angepasst.'}]` (beseitigt den Hinweis „fassung … prüfen“; Status bewiesen bleibt richtig).
   Die englische Fassung in `daten/uebersetzung/gr.yaml` (f4/f5) muss danach nachgezogen werden.
4. `ar.satz.ef_liouville` (07_arithmetik.yaml), Feld `pruefung.f3_spezialfall`: alt „Ja: Σ 1/(b↑↑n) ist aus demselben Grund eine Liouville-Zahl.“ → neu
   „Ja: Σ 1/(b↑↑n) ist für ganze b aus demselben Grund eine Liouville-Zahl (ar.satz.kehrwertsumme_liouville, mit explizitem Kettenbruch); für b = e offen (ar.frage.kehrwertsumme_e).“
   Begründung: ohne „ganze b“ ist die Aussage für b = e unbelegt.
5. `ar.beob.keine_relationen` (id_ar_1.yaml), Feld `aussage`: Satz anhängen „Nicht unabhängig sind $\mathrm{slog}'(0)=1/a_1$ und $\mathrm{slog}''(0)=(a_1^2-2a_2)/a_1^3$ (ar.satz.jet_identitaeten) — rationale Beziehungen, die eine lineare Relationssuche nicht anzeigt.“
   Begründung: B7; die Liste der „50 Konstanten“ enthält zwei abhängige.
6. `ar.satz.keine_arith_info` (07_arithmetik.yaml), Feld `pruefung.f5_neue_frage`: alt „Jeder Beweis über $T_e^{\mathrm K}(1/2)$ muss die komplexe Kneser-Bedingung benutzen.“ → neu
   „Jeder Beweis über $T_e^{\mathrm K}(1/2)$ muss die komplexe Kneser-Bedingung benutzen; als Äquivalenz geschärft und auf $T'(0)$ ausgedehnt in ar.satz.keine_arith_info_ableitung.“

## Quellen, die im Bestand fehlen (nicht angelegt)
- A. J. van der Poorten, J. Shallit, Folded continued fractions, J. Number Theory 40 (1992) — Faltungslemma für B1. Angaben aus dem Gedächtnis.

## Vorschlag Stränge
- Bestehender Strang `ar.klassisch`: aufnehmen ar.satz.kehrwertsumme_liouville, ar.satz.rang_dichotomie_kehrwerte, ar.satz.fixpunkt_transzendent, ar.frage.kehrwertsumme_e.
  `verlauf` ergänzen: „… Die Kehrwertsumme eines Turms ganzer Zahlen ist eine Liouville-Zahl, die zur Basis e ist offen.“
- Bestehender Strang `ar.gleichung`: aufnehmen ar.satz.ableitungsorbit, ar.satz.schanuel_ableitungen, ar.satz.jet_identitaeten, ar.satz.keine_arith_info_ableitung, ar.satz.metrisch_fast_alle.
  Neue `verlauf`-Zeile: „Von zwei Nachbarwerten ist einer transzendent, ebenso von zwei Nachbarsteigungen (T′(1) = e·T′(0)); die reelle Gleichung legt weder T(½) noch T′(0) fest, ‚Schanuel ⇒ transzendent für alle Lösungen‘ ist gleichwertig mit ¬Schanuel, und fast überall ist alles transzendent.“
- Bestehender Strang `ar.ks` (dort steht ar.beob.ausschluss): aufnehmen ar.beob.kf_population.

## Vorschlag Ideenkarte
- Unteridee `zahlen.transzendenz` (Stränge ar.klassisch, ar.gleichung, ar.zaehlen): `kurz` ergänzen zu „Gelfond–Schneider, Schanuel und was die Gleichung über Zahlwerte und Steigungen erzwingt — und warum die Kehrwertsumme eines Turms transzendent ist.“
- ar.beob.kf_population läuft über `zahlen.ks` mit.
