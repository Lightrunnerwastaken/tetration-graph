# Eintrag Bericht 05 (Knesers Welle für komplexe Basen) — 25.10.2026

Dateien: `daten/aussagen/r25_05_kn.yaml`, `daten/theorie/r25_05_kn.yaml`. `python baue.py --nur-pruefen`: 0 Fehler (nur erwartete Hinweise).

## Angelegte ids

| id | Titel kurz | Status |
|---|---|---|
| kn.satz.komplex_radius | Komplexe Basen: Radius ≤ 1, Singularität bei q = 1 bleibt, beide Wellen; keine wandernden Nullstellen (K1, K2) | bedingt |
| kn.beob.komplex_koeffizientengesetz | Koeffizientengesetz mit c = log(±iA/(2π ln b)) oben und unten, Exponent genau 1 (K3, K4, K5) | numerisch |
| kn.beob.komplex_wellen_unsymmetrisch | Obere und untere Welle nicht konjugiert (K6) | numerisch |
| kn.beob.rand_zusatzsingularitaet | Zusatzsingularität von θ_U in der oberen Halbebene bei ln b = 1.1ie^{−i}, nicht in den Vergleichsrichtungen (K7, K8) | numerisch |
| kn.verm.zusatzsingularitaet_bahnwert | z* liegt über dem Bahnwert sexp(4) (Knospe 1/4, fast geschlossene Bahn) (K10, beim Eintragen geschärft) | vermutung |
| kn.verm.a0_rand | a₀ ≈ (7…7.6)/ln\|λ₊\| zum Rand (K9) | vermutung |

K11 (negativ: h₊ bricht bei r = 1.05, h₋ nahe der Achse) ist in die Felder `f2`/`f5`/`gegenbeispielsuche` eingearbeitet.

## Gegenprobe

| Befund | geprüft wie | Ergebnis |
|---|---|---|
| K1 | Beweis gelesen, auch der Fall R = 1 (Fortsetzung in eine Umgebung von q = 1 gibt denselben Widerspruch) | richtig, aber nur gegeben die Kneser-Darstellung für komplexe Basen (Paulsen 2019, Volltext nicht geprüft) → **bedingt** statt „bewiesen“. |
| K2 | gelesen | richtig in der vorsichtigen Form „wo sexp(z) = b^{sexp(z−1)} gilt“; der Satz „die einzige Nullstelle ist z = −1“ ist für (−2, −1) nicht begründet, darum so nicht übernommen. |
| K3, K4 | `e2_welle.py "exp(2*(sin(1)+I*cos(1)))" 30 0.07 512 0 U` (neue Linie y = 0.07) | n·rel. Rest 0.126 (n = 64), 0.120 (96), 0.116 (128); n ln n \|aₙ\| = 0.1385 (n = 128); aₙ gegen y = 0.05 relativ 0 (n ≤ 16), 1·10⁻²¹ (64), 2·10⁻¹⁵ (96), 3·10⁻⁹ (128) — bestätigt. Anmerkung: Die im Bericht genannte Gegenprobe y = 0.08 für r = 2 hat keine Datei in `daten/` (die Motorwerte stehen im Zwischenspeicher, der Lauf ist vermutlich bei der unteren Welle abgebrochen, vgl. K11). |
| K5 | nicht nachgerechnet | bleibt „nur eine Einstellung“, so eingetragen. |
| K6 | Fixpunkte aus dem Lauf y = 0.07 | L₊ = 0.18477 + 0.65020i (\|λ₊\| = 1.35189), L₋ = −0.54770 − 0.90505i (\|λ₋\| = 2.11575), a₀^U = 22.9110 + 3.5238i — wie im Bericht. |
| K7 | `e2_welle.py "exp(1.1*(sin(1)+I*cos(1)))" 30 0.06 512 0 U` und `e3_prony.py … 100 180` | Quotienten r_{n+1}/r_n bei n = 40, 80, 120, 160 auf 8 Stellen wie bei y = 0.05; z* = 0.03388 + 0.01022i (1 Term), 0.03393 + 0.01070i (2 Terme) — bestätigt. |
| K10 | eigene Rechnung `g1_zstern_gegenprobe.py` (Ausgabe `daten/log_g1_zstern.txt`) | Bahn von 0: sexp(3) = 0.03722 + 0.00249i, **sexp(4) = 1.03321 + 0.02524i** — das ist der im Bericht gemessene Wert sexp(z*) ≈ 1.034 + 0.024i, und sexp(7) = 0.04962 − 0.00987i ist der Wert „nach 3 Schritten 0.049 − 0.010i“. Newton für sexp(z) = sexp(4) nahe 0: z = 0.0338189 + 0.0113749i (Probe 9·10⁻³⁹). Realteil wie Prony auf 3·10⁻⁵; der Imaginärteil liegt über den Prony-Werten, passend zur Unterschätzung des Quotiententests bei log-artiger Singularität (\|ρ\|(1 − (1 + 1/ln n)/n) sagt 1.0677 bei n = 199 und 1.0661 bei n = 160 voraus, gemessen 1.0673 und 1.0659). b^z hat einen anziehenden 4-Zyklus (\|Multiplikator\| 0.524). → K10 als geschärfte **Vermutung** eingetragen. |
| K8 | dasselbe Skript für die Vergleichsbasen | Bei ln b = 1.2ie^{−i} (auch 4-Zyklus) und ln b = 1.1e^{i−e^{i}} liegen die Urbilder der Bahnwerte nahe 0 **unter** der Achse — passt zu „keine Zusatzsingularität“. **Korrektur der Deutung:** ln b = 1.1e^{i−e^{i}} ist nicht „irrational ohne Knospe“; b^z hat dort einen anziehenden 6-Zyklus. Entscheidend scheint die Lage des Urbilds, nicht die Rationalität. |
| K9 | a₀ und \|λ₊\| aus dem Lauf y = 0.06 | a₀ = 150.901 + 4.147i, ln\|λ₊\| = 0.04657, Produkt 7.03 — wie im Bericht. |

Neue Dateien der Gegenprobe (nur im Berichtsordner): `g1_zstern_gegenprobe.py`, `daten/log_g1_zstern.txt`, `daten/welle_exp_1.1__sin_1__I_cos_1____30_0.06_512_m0_U.json`, `daten/welle_exp_2__sin_1__I_cos_1____30_0.07_512_m0_U.json`, zusätzliche Einträge in `cache/`.

## Änderungen am Bestand (Vorschläge, nicht ausgeführt)

1. **kn.konstr.kneser**, `pruefung.f4_uebertragbar`: „Auf komplexe Basen (Paulsen 2019) …“ → ergänzen: „; dort muss die Welle θ_U nicht in der ganzen oberen Halbebene holomorph sein (kn.beob.rand_zusatzsingularitaet).“ Begründung: K7 mit Gegenprobe.
2. **kn.satz.paulsen_komplex**, `kanten` ergänzen: `{typ: verwandt, ziel: kn.beob.rand_zusatzsingularitaet, notiz: 'nahe dem Shell–Thron-Rand hat θ_U Singularitäten in der oberen Halbebene'}` und `{typ: verwandt, ziel: kn.satz.komplex_radius}`.
3. **kn.satz.theta_singularitaet**, `pruefung.f4_uebertragbar`: ergänzen „(i)–(ii) gelten wörtlich auch für komplexe Basen ausserhalb Shell–Thron (kn.satz.komplex_radius).“
4. **kn.satz.theta_koeffizienten**, `pruefung.f4_uebertragbar`: ergänzen „Für komplexe Basen numerisch bestätigt, mit c = log(±iA/(2π ln b)) (kn.beob.komplex_koeffizientengesetz); nicht in Knospen nahe dem Rand (kn.beob.rand_zusatzsingularitaet).“
5. **fp.verm.knospen**, `pruefung.f5_neue_frage` („Welche Tetration gehört zu Basen mit anziehendem q-Zyklus?“): Hinweis auf kn.verm.zusatzsingularitaet_bahnwert (Kneser-Welle bekommt dort Singularitäten über den Bahnwerten).

## Stränge

Vorschlag: neuer Strang **`kn.theta_komplex`** (Gebiet kn), Name „Knesers Welle für komplexe Basen“, Mitglieder: alle sechs ids oben.
`verlauf: Für komplexe Basen bleiben die Singularität bei den ganzen Zahlen und das Koeffizientengesetz, oben und unten mit eigener Konstante; nahe dem Shell–Thron-Rand, in einer Knospe, erscheint knapp über der Achse eine neue Singularität über dem Bahnwert sexp(4).`
(Alternative: an `kn.theta` anhängen.)

## Ideenkarte

Unteridee **`fortsetzung.kneser`** („Kneser über η“) — `straenge` um `kn.theta_komplex` ergänzen. Falls man die Welle lieber bei den Formeln führt: `zahlen.identitaeten` (dort steht `kn.theta`).

## Quellen

Keine neuen angelegt. Für kn.verm.zusatzsingularitaet_bahnwert wäre eine Quelle für den Satz von Iversen (ausgelassene Werte ganzer Funktionen sind asymptotische Werte) gut, z. B. F. Iversen, *Recherches sur les fonctions inverses des fonctions méromorphes*, Thèse, Helsingfors 1914 (nicht geprüft).
