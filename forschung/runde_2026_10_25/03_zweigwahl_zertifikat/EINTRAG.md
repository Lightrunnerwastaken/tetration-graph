# Eintrag Bericht 03 (Zweigwahl-Zertifikat in Intervallarithmetik) — 25.10.2026

Dateien: `daten/aussagen/r25_03_kn.yaml`, `daten/theorie/r25_03_kn.yaml`; neu im Berichtsordner `gegenprobe_motor_gegen_zert.py` (Gegenprobe, nur Lesen).
`baue.py --nur-pruefen`: 0 Fehler (nur die erwarteten Hinweise).

## Angelegte Aussagen

| id | Titel kurz |
|---|---|
| kn.lem.zweig_ein_ball | S4 an z₀ ⟸ e^{ℓ(z₀+θ(z₀))} ∈ σ(ℍ); geprüft mit Koebe-¼ um einen freien Testpunkt, ohne Ψ; a₀-Zweig egal (bewiesen) |
| kn.lem.sigma_de_branges | rigoroses σ, σ′ über λⁿσ(log_bⁿ w) und de-Branges-Schranken nahe L (bewiesen) |
| kn.satz.zweigwahl_zertifikat | Zertifikat für Basis e und 2 an je 4 Punkten, Abstand ≈ 10⁻⁶¹ ≪ ρ; Toleranz τ = 0.027 (bedingt auf Motor-Annahme E) |
| kn.satz.sexp_schlicht_nahe_null | Teil A rigoros: sexp schlicht auf \|z\| < 0.1, ε₀, ε₁, sexp′/A, Basis e und 2; Korrektur ε₀ ≤ 0.0277 (bedingt auf E2) |

Befunde 7, 8, 9, 10 sind in kn.lem.zweig_ein_ball bzw. kn.satz.zweigwahl_zertifikat eingearbeitet, 11–13 in kn.satz.sexp_schlicht_nahe_null.

## Gegenprobe

| Befund | geprüft wie | Ergebnis |
|---|---|---|
| B1, B7 (Logik) | nachgelesen: P = σ(w′) ⇒ Ψ(P) = w′ = sexp z₀ ⇒ σ(sexp z₀) = P; Koebe-¼ für die normierte Funktion um w₁; Toleranzformel; P invariant unter θ → θ + 2πik/ℓ | richtig → `bewiesen` |
| B2–B6, B9 | `zertifikat.py` erneut gelaufen (26 s) | alle 8 Prüfungen bestanden, `daten/zertifikat.json` inhaltlich identisch: e, 0.05i: \|P − σ(w₁)\| ≤ 2.16·10⁻⁶¹, ρ ≥ 0.1212, τ_max ≥ 0.027266; Basis 2: 1.56·10⁻⁶¹, 0.0980, 0.027269 |
| Code zertifikat.py | gelesen | L per Intervall-Newton (Einschluss + Eindeutigkeit); jede Log-Iterierte auf Im > 0 geprüft; de-Branges-Schranken richtig (Σ k tᵏ, Σ k² t^{k−1}); ρ mit unterer Schranke von \|σ′\|; w₁ als exakte Zahl; τ_max mit oberen Schranken im Nenner; θ_num summiert alle 4097 Datei-Koeffizienten (a₀ … a₄₀₉₆) in Bällen. **Kein Fehler gefunden.** Die Zweigzahl-Bestimmung (n_h) ist Gleitkomma, aber nur Konvention (Befund 10) |
| B11–B13 | `teil_a.py` erneut gelaufen (6 s), Code gelesen | identisch (e: NW 0.060971 < 1.09177, ε₀ ≤ 0.0276164, ε₁ ≤ 0.0338161, \|sexp′/A − 1\| ≤ 0.0558586, Toleranz 4.39; 2: 0.0046221, 0.00294655, 0.0480237, 0.00514667, 3.77). Noshiro–Warschawski, Restschranken und Bogenüberdeckung richtig. Kleinigkeiten: der Docstring von teil_a.py nennt k ≤ 256, die Datei hat 128 Koeffizienten (k ≤ 127, wie im Bericht); normiert wird mit A aus den Daten (Abweichung ≤ 10⁻³⁰, in den Rundungen enthalten) |
| (E) selbst | **zusätzlich**: Motor-θ gegen das zertifizierte θ̃ aus kn.satz.kneser_e_zertifiziert (`polynom_halbexp/p_pruefung/p26_iter6.json`, \|θ_K − θ̃\| ≤ 4.4·10⁻¹⁶⁷ auf Im z ≥ 0.36), `gegenprobe_motor_gegen_zert.py`, Gleitkomma dps 80 | θ_motor − θ̃ = 3·2πi/L bis auf ≤ 3.2·10⁻⁷⁹ an 0.36i, 0.4i, 0.3+0.5i, i; die ersten 60 Fourier-Koeffizienten gleich auf 5.3·10⁻⁸². Starker Beleg, dass der Motor Knesers θ rechnet — aber kein Beweis von (E) bei z₀ = 0.05i (dafür braucht man θ_K nahe der Achse) |

Fazit: **Das Zertifikat hält.** Die frühere Lücke (Gleitkomma, nur Basis e) ist geschlossen. Es bleibt die Motor-Annahme (E): die Motor-θ-Reihe ist an einem Punkt
auf 0.027 Knesers Welle (mod 2πi/ℓ), dazu für die Hilfsschranken (E2).

## Änderungen am Bestand (Vorschläge, nicht ausgeführt)

**Zum Status von kn.satz.theta_singularitaet und kn.satz.theta_koeffizienten (bedingt → bewiesen?).**
Empfehlung: **Status `bedingt` behalten, aber den `bedingungen_text` durch die viel schwächere Bedingung ersetzen.** Grund: Es bleibt eine unbewiesene Annahme über die
Motordaten. Genau dieselbe Art Annahme führt im Graphen bei gr.satz.ein_wendepunkt_e2 („Annahme H: Motorwerte … auf 10⁻⁵ genau“) zum Status `bedingt`; ein Wechsel auf
`bewiesen` wäre damit uneinheitlich. Wer Motor-Annahmen mit dieser Reserve (Toleranz 0.027, gemessen 10⁻⁵⁹, Abgleich mit dem zertifizierten Kneser 10⁻⁷⁹) als erledigt
zählen will, sollte beide Einträge und gr.satz.ein_wendepunkt_e2 gemeinsam auf `bewiesen` setzen.
**Weg zu `bewiesen` (Basis e), ohne neue Motordaten:** S4 gilt automatisch weit oben (σ(sexp) und e^{L(z+θ)} streben beide gegen 0, Ψ ist dort injektiv) bzw. an einem Punkt mit
Im z ≥ 0.36 per kn.lem.zweig_ein_ball mit dem zertifizierten θ̃; dann überträgt das Hochhebungs-Lemma S4 längs des Wegs iy, 0.05 ≤ y ≤ 0.4, sobald Im sexp_K(iy) > 0 dort
zertifiziert ist (aus dem Polynom T̃ von kn.satz.kneser_e_zertifiziert, Fehler 1.3·10⁻¹⁶⁸). Teil A ebenso mit T̃ statt der Taylor-Daten (streicht E2). Geschätzt eine kurze
Rechnung mit den vorhandenen Werkzeugen (kc_rig, zertifikat.py).

1. **kn.satz.theta_singularitaet**
   - `bedingungen_text`: alt „Zweigwahl an einem Punkt (z₀ = 0.05i, Abstand 4·10⁻⁴⁸ bei nötigen < 0.12) nur in Gleitkomma geprüft, nicht mit Intervallen; nur Basis e“ →
     neu „Motor-Annahme (E): die Motor-θ-Reihe trifft Knesers θ bei z₀ = 0.05i bis auf 0.027 (mod 2πi/ℓ); die Zweigwahl selbst ist mit Intervallarithmetik bewiesen, Basis e und 2
     (kn.satz.zweigwahl_zertifikat); Hilfsschranken nahe 0 unter (E2) (kn.satz.sexp_schlicht_nahe_null). Belege für (E): zwei Motoreinstellungen auf 2·10⁻⁵⁹ gleich; Basis e stimmt
     auf Im z ≥ 0.36 mit dem zertifizierten θ_K (kn.satz.kneser_e_zertifiziert) auf 3·10⁻⁷⁹“.
   - `beweis.luecken`: alt „Zweigwahl: bewiesen ist, … nur in Gleitkomma gerechnet, nicht mit Intervallen, und nur für Basis e.“ →
     neu „Zweigwahl: mit Intervallen zertifiziert für Basis e und 2 (kn.satz.zweigwahl_zertifikat); offen ist nur die Motor-Annahme (E) mit Toleranz 0.027.“
     (Falls der Status auf `bewiesen` geht: `beweis.art` von `skizze` auf `vollstaendig` und die Lücke streichen.)
   - `aussage` (iii): nach der Formel einfügen „(modulo $2\pi i/L$; der Zweig von $h$ ist der, der an $z_0=0.05i$ den Datenwert von θ liefert — die Zweigzahl 0 von `lib3.Basis.h` ist eine Konvention)“. Grund: Befund 10.
   - `pruefung.f4_uebertragbar`: anhängen „Zweigwahl für Basis 2 jetzt ebenfalls zertifiziert.“
   - `bewertung.richtigkeit`: alt „Kern bewiesen, Lokalform numerisch auf 10⁻²⁷“ → neu „Kern bewiesen; Zweigwahl mit Intervallen zertifiziert (Basis e und 2) unter Motor-Annahme (E); Lokalform numerisch auf 10⁻²⁷“.
   - `kanten`/`stuetzt_sich_auf`: `kn.satz.zweigwahl_zertifikat` und `kn.satz.sexp_schlicht_nahe_null` ergänzen.
2. **kn.satz.theta_koeffizienten**
   - `bedingungen_text`: alt „dieselbe Zweigprüfung wie bei kn.satz.theta_singularitaet“ → neu „dieselbe Bedingung wie kn.satz.theta_singularitaet: Motor-Annahme (E); die Zweigwahl ist mit Intervallen zertifiziert (kn.satz.zweigwahl_zertifikat), Basis e und 2“.
   - `bewertung.richtigkeit`: alt „bewiesen, bedingt auf eine in Gleitkomma geprüfte Zweigwahl; numerisch bis n = 4000 auf 5·10⁻⁹“ → neu „bewiesen unter Motor-Annahme (E) (Zweigwahl mit Intervallen zertifiziert); numerisch bis n = 4000 auf 5·10⁻⁹“.
   - `pruefung.f2_benutzte_voraussetzungen`: „(bedingt auf die Zweigprüfung in kn.satz.theta_singularitaet)“ → „(bedingt auf Motor-Annahme E, Zweigwahl zertifiziert)“.
3. Kein Grapheintrag, aber zu korrigieren: `forschung/identitaeten_2026_10_14/luecken/BERICHT_LUECKEN_ANALYTISCH.md`, Satz Z: |ε₀| ≤ 0.027 → 0.0277 (zertifiziert), in (M3) 1.028 statt 1.027.
4. Hinweis für die parallel eingetragene Aussage br.satz.momente_cauchy (r25_04_br.yaml): deren `bedingungen_text` nennt noch die „nur in Gleitkomma geprüfte Zweigwahl“; nach diesem Eintrag gilt dort ebenfalls „Motor-Annahme (E)“.

## Stränge

Bestehender Strang **kn.theta**: die vier ids aufnehmen; `verlauf` ergänzen um „…; die Zweigwahl dahinter ist mit Intervallarithmetik zertifiziert (Basis e und 2), es bleibt eine Motor-Annahme mit Toleranz 0.027.“
(Alternative: eigener kleiner Strang kn.theta_zertifikat „Zweigwahl streng“.)

## Ideenkarte

Unteridee **zahlen.identitaeten** (enthält kn.theta); keine Änderung nötig.
