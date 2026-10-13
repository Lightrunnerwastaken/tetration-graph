> **English:** [README.md](README.md). Die Ansicht startet auf Englisch; der Knopf **DE/EN** schaltet die Sprache um (die Wahl wird gemerkt).

# Tetration-Wissensgraph

Überprüfbarer, aussagenorientierter Wissensgraph zur Tetration und ihren Verbindungen
(Funktionalgleichungen, Iterationstheorie, Analysis, Arithmetik, Logik, Dynamik).

**Ansehen:** `ausgabe/tetration-graph.html` direkt im Browser öffnen (eine Datei, offline).
**Prüfbericht:** `ausgabe/pruefbericht.md` · **Landkarte als Text:** `ausgabe/landkarte.md`.

**Stand 12.10.2026:** 367 Aussagen · 637 Kanten · 54 Experimente (alle aktuell) · 292 Quellen · 9 Gebiete · 8 Hauptideen mit 37 Unterideen.
Prüfprogramm: 0 Fehler, 0 Warnungen, 26 Hinweise. Zahlen ändern sich mit jedem Bau; massgebend sind der Kopf von `ausgabe/pruefbericht.md`
und `ausgabe/graph.json` (`statistik`).

```
python baue.py                 # prüfen + bauen
python baue.py --nur-pruefen   # nur Prüfbericht
python baue.py --experimente   # vorher alle Experimente neu rechnen
python -m pytest -q tests      # Regeln des Prüfprogramms
```
Benötigt Python 3 mit `pyyaml` für Prüfen und Bauen (`pip install pyyaml pytest`). Die Experimente brauchen zusätzlich
`mpmath`, `numpy`, `scipy`, `matplotlib`; manche (Kneser-Motor, Hexation) auch PARI/GP. Wer ohne Experimente nur prüfen und bauen will,
braucht nur `pyyaml`.

## Aufbau

| Ort | Inhalt |
|---|---|
| `daten/aussagen/*.yaml` | eine Datei pro Gebiet (`00_`…`10_`) und Ergänzungsdateien `lit_<gebiet>_*.yaml` der Literaturkarte; jede Aussage mit Objekten, Bereich, Voraussetzungen (je „benutzt: ja/nein/unklar“), Schluss, Beweis + Lücken, Quellen, Kanten, Status, Bewertung (Richtigkeit/Neuheit/Nutzen getrennt), fünf Prüffragen, Fassungen |
| `daten/quellen.yaml`, `daten/quellen/*.yaml` | Literatur mit Prüfvermerk `original` / `sekundaer` / `nein` (Ergänzungen der Literaturkarte in `quellen/lit_L<n>.yaml`; der Lader führt zusammen, ein Schlüssel in zwei Dateien ist ein Fehler) |
| `daten/experimente.yaml`, `daten/experimente/*.yaml` | Experimente → Skript, Ergebnis, Bild, gestützte Aussagen |
| `daten/objekte.yaml`, `daten/objekte/*.yaml` | Vokabular der Objekte |
| `wissen/` | Schema (Vokabular), Laden, Prüfen/Ableiten, Ausgabe |
| `experimente/e*.py` | reproduzierbare Rechnungen mit Selbstprüfungen; Ergebnis-JSON trägt den Skript-Hash |
| `daten/straenge.yaml`, `daten/straenge/<gebiet>.yaml` | Forschungsstränge je Gebiet mit einer Zeile „was dort passiert“ (mittlere Zoomstufe der Karte, Lesewege, Lehrbuch) |
| `daten/ideen.yaml` | Ideen-Landkarte: 8 Hauptideen, 37 Unterideen; jede Aussage gehört über ihren Strang (oder ein Feld `idee:`) zu genau einer Unteridee |
| `wissen/landkarte.py` | berechnet aus den Kanten Kopplung, Lage und Füllstand der Ideen; schreibt `ausgabe/landkarte.md` |
| `forschung/` | Arbeitsordner der Forschungsstränge (`polynom_halbexp`, `literaturkarte`, `eta_hyperoperationen`); der Graph übernimmt nur Gesichertes daraus |
| `docs/` | Recherche-Berichte, Entwürfe und Pläne (`docs/superpowers/specs`, `plans`) |
| `ansicht/` | Vorlage, Stil, Module `js/01_basis … 06_steuerung` (werden beim Bau zu einer Datei zusammengefügt; KaTeX eingebettet) |

## Ideen-Landkarte (Startansicht, ab 12.10.2026)

Die Karte öffnet mit der **mathematischen Struktur** statt mit Orten: acht Hauptideen (Die Gleichung · Dynamik von b^z und der Fixpunkt L ·
Unendlicher Turm und Basen-Ebene · Fraktionale Fortsetzung · Basiswechsel · Höhere Stufen · Zahlen und Arithmetik · Brücken), deren Lage
sich aus der Kopplung (Kanten zwischen ihren Aussagen, gemeinsame Objekte) berechnet. Klick auf eine Hauptidee zeigt ihre Unterideen,
Klick auf eine Unteridee ihre Aussagen (Hauptaussagen gross, Hilfsaussagen klein), Klick auf eine Aussage die Inhaltskarte. Pfadleiste oben,
Zurück mit Esc.

- **Füllstand-Ring:** Anteile nach Status (bewiesen/bedingt/numerisch/Vermutung/offen/verworfen). **Herkunfts-Ring:** Literatur / Forum / nur eigene Arbeit.
- **Stücke einer Vermutung:** Bei Vermutung, bedingtem Satz oder Skizze zeigt die Inhaltskarte den Baum der Voraussetzungen und markiert, wo es fehlt.
- **Textfassung** für Forschungs-Agenten: `ausgabe/landkarte.md`.
- Regeln: Aussage ohne Idee → Hinweis „ohne-idee“; unbekannte Idee oder Strang in zwei Unterideen → Fehler (Tests: `tests/test_ideen.py`, `tests/test_landkarte.py`).
- Entwurf und Plan: `docs/superpowers/specs/2026-10-12-ideen-landkarte-design.md`, `docs/superpowers/plans/2026-10-12-ideen-landkarte.md`.

## Karte (zweite Ansicht)

Die frühere zoombare Karte (Gebiete → Verläufe → Aussagen → Inhalt) bleibt als zweite Ansicht erhalten. Eine zoombare Karte auf dunkler Blaupause; was man sieht, hängt von der Zoomstufe ab (Leiste links):

1. **Gebiete** — Inseln mit Statusring, Bänder mit fliessendem Licht zeigen die Querverbindungen.
2. **Verläufe** — Forschungsstränge als Flächen mit Name und einer Zeile Verlauf; Licht fliesst entlang der Beweise von der Voraussetzung zur Folgerung.
3. **Aussagen** — Titel, wichtige zuerst, Überlappendes bleibt verborgen.
4. **Inhalt** — Karten mit der gesetzten Aussage neben den Knoten.

Klick auf eine Aussage: Kamera fliegt hin, Nachbarn leuchten, Nachbarn zweiten Grades halb, Rest tritt zurück; die Detailtafel
öffnet sich nur dann. Esc oder Klick ins Leere hebt den Fokus auf. Bei „reduzierter Bewegung“ entfallen Animationen.

Ebenen: **A** jede Lösung der Funktionalgleichung · **B** Klasse von Lösungen · **C** konkrete Konstruktion · **X** Hilfssatz aus anderem Gebiet.

Niveau (für die öffentliche Bibliothek): **Einstieg** (Schulmathematik plus Ableitung/Integral; Gebiet „Grundlagen“ und markierte
Basisaussagen) · **Vertiefung** (Studium) · **Forschung** (offene Fragen, eigene Arbeit). Ohne Eintrag gilt: Grundlagen = Einstieg,
Basiswechsel und Vermutungen/Beobachtungen = Forschung, sonst Vertiefung. Einstiegsaussagen tragen `einfach:` (Erklärung in Alltagssprache),
Identitäten `formen:` (gleichwertige Schreibweisen).

**Englische Fassung (Werkzeuge ab 12.10.2026):** Übersetzungen stehen getrennt von den deutschen Texten in `daten/uebersetzung/<gebiet>.yaml`
(je Aussage `en` und optional `en_theorie`, dazu Hash der übersetzten deutschen Texte, Herkunft und Prüfstand). Das Prüfprogramm erzwingt, dass
Formeln, `[[Verweise]]` und Kennungen unverändert bleiben, meldet veraltete Übersetzungen (deutscher Text hat sich geändert), Zahl- und
Wortwahl-Abweichungen und den Fortschritt je Gebiet. Befehle (`stand`, `export`, `einlesen`), Glossar und Pilot: `docs/uebersetzung/`.
**Stand 12.10.2026:** alle 367 Aussagen (Kern, Beweis, Prüffragen, Theorie), Gebiete, Stränge, Ideen, Lehrbuch-Kapitel, Objekte und Experiment-Titel
sind übersetzt (`daten/uebersetzung_struktur/<art>.yaml`, Befehle `struktur-stand|-export|-einlesen`); die englische README ist `README.md` (Hauptfassung auf GitHub).
Gegengelesen ist erst ein Teil (Prüfstand: `python -m wissen.uebersetzung stand`). Offen: DE/EN-Umschalter der Ansicht, Oberflächentexte,
Ergebnis-Zusammenfassungen der Experimente, `ausgabe/landkarte.md` und `docs/`. Unstimmigkeiten im deutschen Text, die dabei aufgefallen sind:
`docs/uebersetzung/DEUTSCH_ZU_KLAEREN.md`.

## Was das Prüfprogramm erzwingt

- Verweise, Vokabular, keine Kreise in Beweisketten.
- „bewiesen“ auf schwächerem Fundament wird automatisch „bedingt“; für jede bedingte Aussage werden die
  **Wurzelbedingungen** (unbewiesene Enden der Kette) berechnet und müssen genannt sein.
- Eine Annahme gilt als entlastet, sobald ein bewiesener Satz sie entlastet.
- Hängt etwas (transitiv) an einer widerlegten/überholten Aussage → „neu prüfen“, bis ein `nachgeprueft`-Eintrag es aufhebt.
- Unbenutzte Voraussetzungen → Hinweis „abschwächbar“. Vermutungen brauchen eine dokumentierte Gegenbeispielsuche,
  Skizzen benannte Lücken, „vermutlich neu“ eine Recherche-Angabe.
- Experimente: Ergebnis gilt nur, wenn der Skript-Hash passt; fehlgeschlagene Selbstprüfungen werden gemeldet.

## Neue Aussage eintragen

```yaml
- id: kn.satz.beispiel            # Präfix = Gebiet
  titel: Kurzer Titel
  art: satz                        # definition|satz|lemma|korollar|konstruktion|annahme|vermutung|beobachtung|gegenbeispiel|gescheitert|korrektur|frage
  ebene: B                         # A|B|C|X
  status: bewiesen                 # siehe wissen/schema.py
  objekte: [kneser_T]
  bereich: '$b>\eta$'
  voraussetzungen:
    - {text: 'T streng wachsend', benutzt: ja}
  aussage: |
    Text mit $\TeX$.
  beweis: {art: vollstaendig, text: '…', luecken: []}
  stuetzt_sich_auf: [kn.konstr.kneser]
  kanten: [{typ: spezialfall_von, ziel: fg.satz.allgemeine_loesung}]
  quellen: [{ref: kneser1950, stelle: 'Satz 1', geprueft: nein}]
  bewertung: {richtigkeit: '…', neuheit: bekannt, nutzen: '…', nutzen_stufe: mittel}
  niveau: einstieg                 # einstieg|vertiefung|forschung (optional)
  einfach: >-                      # Pflicht bei Einstieg: in einfachen Worten
    …
  formen:                          # gleichwertige Schreibweisen
    - {tex: '…', notiz: '…'}
  pruefung: {f1_folgt_unmittelbar: '…', f2_benutzte_voraussetzungen: '…', f3_spezialfall: '…', f4_uebertragbar: '…', f5_neue_frage: '…'}
```

## Lehrbuch (Theorie-Reiter, ab 08.10.2026)

Der Reiter „Theorie“ ist ein Lehrbuch (Umschalter „Abhängigkeiten“ zeigt das alte Diagramm). Kapitelfolge: `daten/lehrbuch.yaml`
(Kapitel 0 „Was ist Tetration?“ als Fliesstext, dann je Kapitel Stränge). Ausführliche Erklärungen je Aussage: `daten/theorie/<gebiet>.yaml`
mit den Teilen `worum`, `genau`, `warum`, `beispiel`, `bilder`, `einordnung` (Markdown mit `$…$`, `[[id]]`, `![[bild:name]]`).
Abbildungen: `python lehrbuch/bilder.py` → `lehrbuch/bilder/*.svg`. Das Prüfprogramm meldet tote Verweise, fehlende Teile und den
Fortschritt je Gebiet. **Stand 12.10.2026: vollständig** — alle 367 Aussagen haben einen Theorie-Block (Stand 08.10.: 218 Aussagen in Kapitel 0–10, damals je Kapitel:
Grundlagen 32, Fixpunkte/Turm 7, Funktionalgleichung 11, Iteration 8, unter η 19, Kneser 19, Basiswechsel 65, Basen-Ebene 25,
Zahlen 16, Brücken 16). Neue Aussagen brauchen einen Block in `daten/theorie/<gebiet>.yaml`, sonst meldet die Prüfung
„theorie-fortschritt“. Das Phasenbild k7_phase
rechnet slog_e(T₂(n+θ)) selbst nach (Tschebyschow-Interpolanten aus Kneser-Werten, ab n = 5 über die Logarithmenleiter) und trifft
phi_modes.json auf 4·10⁻²⁴.
Anzeigeformeln der Form „A,\qquad B“ oder mit „,\quad“ setzt die Ansicht untereinander, wenn sie nebeneinander nicht passen (`formelnEinpassen`).

## Abgleich mit Wikis und Tetration-Forum (08.10.2026)

Bericht: `docs/recherche_wikis_forum.md` (Wikipedia en/de/ru, Googology, MathWorld, Mizugadro-Archiv, Forum: FAQ, TPID-Liste,
1260 Thementitel, ~60 Themen gelesen). **Runde 1 erledigt:** Vorarbeiten aus dem Forum eingetragen — Basiswechsel-Konstante μ
(Fox 2007), Werte m(2), m(3), m(10) und das „Wackeln“ (Levenstein 2009), zwei Abel-Funktionen zu √2 (Asimov–Hickerson–Schroeppel um
1991, „Bummer!“ 2007), Entwicklung bei η (Trappmann–Kouznetsov 2012); Gegen-Vermutung „Basiswechsel nirgends analytisch“
(Levenstein 2011) als `bw.verm.nirgends_analytisch`. Der Lader meldet jetzt doppelte YAML-Schlüssel (vorher stand `justus_p2` doppelt).
**Runde 2 erledigt:** 13 Kern-Aussagen — Superwurzeln (Strang `gr.umkehrungen`, Arithmetik), unendliche Superwurzel, TPID 6 mit
neuer Rate $n^2(x_n-\eta)\to2\pi^2\eta/e$, obere Superexponentialfunktion und cheta, q-Binomial-Reihe (Reshetnikov) mit Rate $q^{\operatorname{Re}z+2}$,
Vermutung „T′ vollständig monoton ⇔ Koenigs-Koeffizienten negativ“, Möbius-Iteration, Szekeres' Monotonie, Shishikura/Lavaurs;
Experimente e28–e31. **Runde 3 erledigt:** Strang `kn.hoehere_stufen` — Pentation als reguläre Superfunktion von T_e am Fixpunkt
−1.8503545290… (λ = 6.4607; Fox 2007 auf 22 Stellen bestätigt), η₄ = 1.63532449671527639934… und e₄ = 3.08853227180671765448… (Levenstein 2010
auf 32 Stellen bestätigt) mit drei Kennzeichnungen; Experimente e32–e33, Hilfsmodul `experimente/stufen.py`. **Runde 4 erledigt:** offene
Vermutungen mit Stichproben (e34–e37) — TPID 7 gilt für Kneser nur bis b₇ = 10.2210634158 (Gegenbeispiel darüber); negative Koenigs-Koeffizienten
würden TPID 7 und Nixons Halbebenen-Frage zugleich beantworten; η₅ = 1.7373559252, e₅ = 3.3014072785 (TPID 10); Carleman-Eigenwerte → λ^k;
Ableitungs-Vorzeichen; Shell–Thron-Fläche 8.0509332809 (Robbins 2009 bestätigt); Geschichte in Kapitel 0. Damit sind alle vier Runden des Berichts erledigt.

## Hexation (09.10.2026)

Reguläre Hexation unter η₅ am anziehenden Fixpunkt der Pentation (`experimente/hexation.py`, `e38_hexation.py`, `exp.hexation` 10/10):
1.7↑↑↑↑½ = 1.3974472613, 1.7↑↑↑↑∞ = 2.3876372817; jede Hexation hat H(−3) = −2 und eine Singularität in (−4, −3); die unendliche Hexation
konvergiert genau bis η₅; an einem anziehenden komplexen Fixpunkt gibt es keine Kneser-Konstruktion. Kneser-artige Hexation über η₅
über das Rechenwerk knhexinit im gp-tetration-Fork (exp-079): 2↑↑↑↑1.5 = 2.6728939593426194891907964…, und damit η₆ = 1.79912123612881268672356700703,
e₆ = 3.44223464771689323 (`e39_hexation_kneser.py`, `exp.hexation_kneser` 11/11; viertes Glied für TPID 10, 2 − η_n ≈ 2.8·n^−1.47).
Blattverfolgung (Pentation entlang von Wegen fortgesetzt) gebaut und geprüft; damit gefunden: die Kneser-artige Hexation ist in der oberen
Halbebene verzweigt (bei −0.7408 + m + 2.7458i; Log-Umlauf 2πi/ln 2, Potenz-Exponent ln λ_T/ln 2 = 2.3707) — Kouznetsovs Rechteck-Methode fällt
damit weg. Prototypen in `experimente/prototypen_hexation/LIESMICH.md` (schreiben einen eigenen
Zwischenspeicher `kneser_cache_proto.json`; der gemeinsame `experimente/daten/kneser_cache.json` ist durch frühere Prototyp-Läufe auf ~190 MB gewachsen).

## Literaturkarte (ab 12.10.2026)

Auftrag: die Arbeiten anderer so genau kartieren wie die eigene Forschung (wer woran gearbeitet hat, mit welcher Methode und welchem
Ergebnis, wo mehrere dieselbe Frage verschieden angegangen sind). Vorher: 321 Aussagen, davon 158 nur mit eigenen Quellen. Drei Stränge
(L1 Iteration/Funktionalgleichungen/Fixpunkte, L2 Kneser/Basiswechsel, L3 Grundlagen/Arithmetik/Brücken) legten ihre Ergänzungen in
eigene Dateien (`daten/aussagen/lit_*.yaml`, `daten/theorie/lit_*.yaml`, `daten/quellen/lit_L*.yaml`, Experiment-Nummernbereiche e60–e89 vorgesehen);
Block 3 (freie Volltexte: Euler E489, Cho–Kim, Kouznetsov–Trappmann, TPID 16/17) ist abgeschlossen, ebenso die Kandidatenliste.
Karten, Berichte und das Gesamtverzeichnis stehen in `forschung/literaturkarte/` (`GESAMTVERZEICHNIS.md`, `L<n>_*/KARTE_<gebiet>.md`).
Folge: 46 neue Aussagen (367), 60 statt 35 Aussagen „bewiesen_literatur“, **103 Zitate noch nicht am Original geprüft** (vorher 68).

## Polynom-Wurzel von exp und höhere Hyperoperationen (10.–12.10.2026)

- **Polynom-Wurzel** (Strang `kn.halbe_exp`, Arbeitsordner `forschung/polynom_halbexp/`): Phasen 2–3 sind in den Graphen übernommen
  (Experiment `exp.polynom_phase3`). Der Hauptsatz steht als **Skizze** mit benannten Lücken (Annahme (H₀) über die Eckmoden bei L ist offen;
  tragfähige Fassung `kn.satz.polynom_hauptsatz_kr0`; Knesers sexp zur Basis e ist streng zertifiziert, `kn.satz.kneser_e_zertifiziert`).
  Phase 3 wurde auf Janis' Wunsch am 12.10. beendet; Wiederaufnahme nur auf Auftrag (Einstieg: `forschung/polynom_halbexp/README.md`).
- **η höherer Hyperoperationen** (`forschung/eta_hyperoperationen/`, eigenständig ab 12.10.): η-Tabelle für höhere Ränge mit
  Genauigkeitsangaben, ausdrücklich ohne allgemeinen Kneser-Satz für alle Ränge; Bericht und Mathematik in `BERICHT.md` / `MATHEMATIK.md`.

## Wiedergefundene Untersuchungen (Juli/August 2026, eingearbeitet 08.10.2026)

Quellen: `Downloads/Base_Change_Shared_Fixed_Point.zip` (30.07.), `Imaginary_Direction.zip`, `Periode_2_Region.zip` (07.08.).
Kernzahlen unabhängig nachgerechnet in `exp.funde` (14/14). Die Rekonstruktion `kl.satz.basiswechsel_gleiches_L` stimmt mit dem Original überein.

- **Geteilter Fixpunkt:** A_b = σ·A_d + Ψ(κ_d) exakt; Musterpaar −19/25 ↔ d = 1.425186736 − 0.046203369i, σ = −0.200461 + 0.003486i.
  Anziehende Partner nur im Keil Re b < 0 (Hauptzweig: nur m = ±1); 2, e, 3, 4 haben keinen (`fp.beob.partner_keil`).
- **Die Welle ist gedämpft:** Pseudoperiode 2π/|arg λ|, Dämpfung |λ|^{2π/|arg λ|}; die klassische 1-periodische Welle ist der Grenzfall |λ| → 1 (`kl.satz.invariante_daempfung`, Deutung).
- **Imaginäre Richtung:** A_b(T_d(iy)) → A_b(L_d), z. B. A_e(L_2) = 0.326448651629 + 0.523513968910i (`bw.satz.imaginaere_richtung`).
- **Periode-2-Region:** geschlossene Form μ = (s/(2 sinh(s/2)))², kein Kreis (Überlagerung statt Bijektion), b_links = −0.08052815762548356,
  Spitze √6/e·n^{−1/2}, unter e^{−e} wird die Phase zur Rechteckwelle mit |c_k| = |z_A − z_B|/(πk) (`fp.periode2`).

## Offene Punkte bearbeitet (08.10.2026)

- **Lambert → Digamma** (`kl.satz.moebius_qdigamma`, `exp.lambert_zeta` 8/8): Im Möbius-Modell ist die Lambert-Bruchsumme exakt eine
  q-Digamma-Funktion mit q = λ; λ → 1 gibt ψ. Die echten Bruchsummen S_b(1/2) laufen gegen S_η(1/2), Fehler O(η − b).
- **Stetigkeit bei η** (`kl.beob.stetig_eta`): T_b(1/2) → T_η(1/2) von unten (regulär) und von oben (Kneser), beide linear in |b − η|
  und mit gleicher Steigung 0.438579 / 0.438581 — numerisch C¹ über η hinweg.
- **Bruchsumme des slog** (`kn.satz.abel_bruchsumme`, `exp.abel_bruchsumme` 9/9): konvergiert, weil A″ = O(1/y²); S(1/2) = −0.1150976090828,
  γ_A = 0.4951728249. Offen bleibt der Bezug zu μ.

## Juli-Synthese eingearbeitet (08.10.2026)

Quelle: `Downloads/tetration_shell_thron_fractional_sums_base_classes.md` (eigene Synthese, Juli 2026). 13 neue Aussagen, `exp.koenigs_kalkuel` (12/12), `exp.eta_digamma` (12/12).

- Koenigs-Spektrum unter η: Ableiten, Integrieren, Shift-Summen als Faktoren je Mode; Bruchsumme als Lambert-Reihe (`kl.satz.koenigs_spektral`).
- Wertseitiger Basiswechsel T_b∘A_d = Potenz χ_d^ρ in Koenigs-Koordinaten, Verzweigung (a_d − y)^ρ, Jacobi-Formel über Bruchsummen;
  quasikonformer Basiswechsel und das Hindernis: holomorph + eindeutig + zeiterhaltend nur bei gleichem Multiplikator (`kl.satz.qk_basiswechsel`).
- Basis η: Fatou-Reihe mit rationalen Koeffizienten (1/18, 1/135, −1/972, …), C_T = 2.79824815423138766245… (nur ganze Höhen);
  Bruchsumme S_η(x) = e x − 2e log x − 2eγ + Σ[T(k) − e + 2e/k] + O(log x/x), Σ[T(k) − e + 2e/k] = 10.7975622749393143…
- Kehrwert-Turmsummen C_b = Σ 1/(b↑↑k) (C_2 = 0.8125152587890625…), Normalform e^w + κ, Atlas der Basenklassen.

## Rekonstruktion „linearer Basiswechsel bei gleichem L“ (08.10.2026)

Die Original-Untersuchung (claude.ai-Chat, Anfang August 2026) ist lokal nicht vorhanden; neu hergeleitet in `exp.gleiches_L` (21/21).

- Zu einem Fixpunkt L gehört eine Leiter von Basen u_k = (Log L + 2πik)/L. Reelle Paare: e^{π/2+2πk} teilen L = i
  (4.8105 divergiert, 0.00898 hat Periode 2), √2 und −√2 teilen L = 2 und 4, dazu 8 gemischte Paare positiv/negativ (`fp.satz.gleiches_L`).
- Am gemeinsamen L ist A_b(T_d(x)) = C₁ + C₂·x + O(λ_d^x), C₂ = log λ_d / log λ_b: die Abweichung fällt wie |λ_d|^x, es bleibt
  keine Welle. Bei verschiedenen L wird A_b(T_d(x)) konstant. Exakt linear nur für b = d (`kl.satz.basiswechsel_gleiches_L`).

## Befunde der Runde „Regionen, Bruchsummen, Turm-Integrale“ (08.10.2026)

Grundlage: Bild 1 und 2 des Nutzers und seine Notiz *Complex Tetration Dynamics* (5. Jan. 2026). 25 neue Aussagen, Experimente e16–e20.

- **Regionen der Basen-Ebene** (`fp.beob.regionen`, `exp.parameterraum`): Periode 1 genau für |t(b)| < 1 mit t = −W₀(−Log b),
  ohne Gegenbeispiel in beide Richtungen; Knospen p/q drehen mit Drehzahl p/q (17/17); an rationalen Randpunkten Konvergenz wie n^{−1/q}.
- **Notiz bestätigt:** Satz 1 (gilt sogar für jeden Radius), Satz 2, Korollare 1–2, alle Zahlen zu b = i. Die Scheibenmittel-Identität (5)
  ist jetzt ein bedingter Satz (Weyl + Mittelwertsatz, Annahme: Siegel-Rand ist Jordankurve). Die Vermutung „Kreismittel = 1 auch für
  nichtganze Höhen“ ist für die reguläre Familie **widerlegt** (Verzweigung w^t).
- **Unter η:** Halbschritt über den Grenzwert (Streckung um λ^x am Fixpunkt), Bruchsummen nach Müller–Schleicher, und
  die Ableitungsformel T′(x+h)/T′(x) = (ln b)^h ∏ T(x+k) mit Bruchprodukt gilt **genau** für die reguläre Tetration (`kl.satz.bruchprodukt`).
- **Linearer Basiswechsel** (Rekonstruktion, `kl.satz.koenigs_basiswechsel`): zwei Basen mit anziehendem Fixpunkt sind über ihre
  Koenigs-Funktionen konjugiert; die Höhe wird exakt k·x mit k = ln λ_d / ln λ_b — ohne periodischen Teil.
- **Turm-Integrale:** lim I₂ₙ = 0.67671427377181279519…, lim I₂ₙ₊₁ = 0.63130937080198819539…, Lücke 0.0454049029698245998…;
  Bild 2 hat die Lücke richtig, beide Grenzwerte aber um 2.428·10⁻¹¹ zu klein. Konvergenz wie n^{−3/2} (vorhergesagt, gemessen 1.502).

## Befunde der Grundlagen-Runde (08.10.2026)

- **Wendepunkte:** Jede gemessene Kneser-Tetration hat genau einen Wendepunkt (b = 1.5: 3.987, 2: −0.030, e: −0.5218, 3: −0.604,
  5: −0.847, 10: −0.997, 100: −1.182); reguläre unter η keinen. Halb bewiesen: Konvexität vererbt sich nach oben (`gr.verm.ein_wendepunkt`).
- **Zahlen:** T′ₑ(0) = 1.0917673512583…, ∫₀¹Tₑ = 1.7145219588906204881…, ∫₁ᵉ slogₑ = 1.0037598695684…; slog′ₑ(0) = 1/T′ₑ(0) = 0.91594605649953…
  ist zugleich der Orbitwert −V₁ aus Paper V.
- **Eigene Fehler gefunden und behoben:** erste Formel für T″(x+1) hatte einen Faktor T(x+1) zu viel; eine grobe PSLQ-Suche lieferte
  eine Scheinrelation für ∫₀¹Tₑ (beides als Fassung bzw. Prüfung dokumentiert).

## Befunde aus dem Aufbau (08.10.2026)

- **Neue Vermutung, numerisch stark gestützt** (`kl.verm.streifenbreite`, `kl.verm.faktor`): Für $1<b<\eta$ fallen die
  Fourier-Moden der Differenz der regulären Abel-Funktionen an beiden Fixpunkten wie $C(b)\,k^{\approx-5/4}e^{-2\pi^2k/|\ln\lambda_a|}$;
  Rate auf $\le2.6\cdot10^{-4}$ bestätigt (5 Basen), Restfaktor über die Basen gleich bis $4\cdot10^{-4}$. Modellfall für die Streifen-/Gevrey-Frage der Basiswechsel-Phasen.
- **Struktur:** 8 der 13 bedingten Sätze (Papers I, II, IV) hängen an einer einzigen Wurzel, Hyp. 2.2 (reeller Logarithmenkanal);
  deren Ablösung über die gestoppte Kette (Paper III) ist nur skizziert (`bw.satz.phasengesetz_gestoppt`).
- **Erratum Paper IV:** $\pi\sqrt{2\eta/e}=3.238930\ldots$, nicht 3.23907 (`bw.korr.konstante`).
- **Datenqualität:** `mu_hub.json` trägt $m(2)$ nur auf ~31 von 40 angezeigten Stellen (`bw.beob.hub_genauigkeit`).
- **Abgeschwächte Voraussetzungen:** Paper I Prop. 7.1 braucht nur Differenzierbarkeit im Punkt −1; A1 Sätze 5, 6 und Kor. 1 gelten
  für ganze Lösungsklassen, nicht nur für Kneser; Paper II Thm. B ist Spezialfall der Youngschen Flächenidentität.
- **Carleman-slog:** konvergiert gegen Kneser, aber nur algebraisch ($\approx n^{-3.8}$, 8 Stellen bei $n=128$).
- **Eigener Fehler gefunden:** erste Fassung von `fp.satz.lambert` hatte das Vorzeichen von $\operatorname{Im}L$ falsch (`exp.fixpunkte`).

## Lizenz

Code (`baue.py`, Python-Skripte, JavaScript/CSS der Ansicht) unter der MIT-Lizenz ([`LICENSE`](LICENSE)); Aussagen, Texte,
Berichte, Daten, Abbildungen und die erzeugte Ansicht unter CC BY 4.0. Einzelheiten und Ausnahmen (KaTeX, Zitate aus fremden
Quellen) in [`LICENSE-INHALT.md`](LICENSE-INHALT.md).
