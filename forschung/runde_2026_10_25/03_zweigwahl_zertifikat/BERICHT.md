# 03 Zweigwahl-Zertifikat in Intervallarithmetik (25.10.2026, 13:22–13:44)

**Kurz gesagt.** Die beiden Sätze kn.satz.theta_singularitaet und kn.satz.theta_koeffizienten hingen an einer einzigen Prüfung:
Stimmt an einem Punkt z₀ = 0.05i die Kneser-Welle mit der Schröder-Funktion überein? Das war bisher nur mit gewöhnlichen
Kommazahlen gerechnet. Jetzt ist es mit Intervallen (python-flint, arb/acb) bewiesen, für Basis e **und** Basis 2, an vier Punkten.
Dabei zeigt sich: Man braucht die Poincaré-Funktion Ψ gar nicht, und die Motordaten müssen nur auf etwa **2 Hundertstel** stimmen
(gemessen stimmen zwei Motoreinstellungen auf 10⁻⁵⁹ überein). Auch die Hilfsschranken für sexp nahe 0 (Teil A) sind jetzt
rigoros, zum ersten Mal auch für Basis 2.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit |
|---|---|---|---|---|
| 1 | Zurückführung: S4 an z₀ folgt schon aus P := e^{ℓ(z₀+θ_K(z₀))} ∈ σ(ℍ); Ψ muss nicht ausgewertet werden | bewiesen | Abschnitt 1 (drei Zeilen aus Lemma U(c) und V1) | neu (vereinfacht Satz Z) |
| 2 | Basis e, z₀ = 0.05i: \|P − σ(w₁)\| ≤ 2.2·10⁻⁶¹ < ρ ≥ 0.1212 (Koebe-¼-Radius) | bewiesen (Intervall), bedingt nur auf (E) | `zertifikat.py`, prec 800 Bit | ersetzt das Gleitkomma-Zertifikat (4·10⁻⁴⁸) |
| 3 | Basis 2, z₀ = 0.05i: \|P − σ(w₁)\| ≤ 1.6·10⁻⁶¹ < ρ ≥ 0.0980 | bewiesen (Intervall), bedingt nur auf (E) | `zertifikat.py` | neu (Basis 2 war nie gerechnet) |
| 4 | Dieselbe Prüfung an 0.07e^{iπ/3}, 0.08e^{3iπ/4}, 0.5+0.05i, beide Basen: alle bestanden (ρ ≥ 0.031) | bewiesen (Intervall), bedingt auf (E) | Tabelle unten | Gegenprobe |
| 5 | Toleranz der Motor-Annahme (E): es genügt \|θ_K(z₀) − θ_num(z₀)\| ≤ τ_max = 0.0273 (e) bzw. 0.0273 (2), modulo 2πi/ℓ | bewiesen | `zertifikat.py`, rigorose untere Schranke | neu |
| 6 | (E) selbst: zwei Motoreinstellungen (dps 80, y = 0.003 / dps 60, y = 0.002) geben θ(0.05i) bis auf 2.1·10⁻⁵⁹ (e) bzw. 1.6·10⁻⁵⁹ (2) gleich | numerisch | 57 Grössenordnungen unter τ_max | – |
| 7 | Die Konstante a₀ (und damit ihr Zweig) spielt für S4 keine Rolle: P ändert sich nicht bei θ → θ + 2πik/ℓ | bewiesen | trivial | Klarstellung |
| 8 | Nebenprodukt: sexp(0.05i) liegt im Ball um w₁ mit Radius 2.4·10⁻⁶² (e), 1.8·10⁻⁶² (2), aus θ allein (Koebe-Wachstum); stimmt mit den unabhängigen Taylor-Daten auf 2.5·10⁻⁶³ (e) bzw. 1.3·10⁻⁶³ (2) | bewiesen (bedingt auf θ_num exakt) / Abgleich numerisch | `zertifikat.py` | Gegenprobe zweier Datenarten |
| 9 | „Zweigzahl 0“ für θ = 3 − z + h(u₃) mit dem h aus lib3/kern: an allen 8 Punkten k = 0, Abstand zur Sprungstelle ≥ 0.07 (in Einheiten von 1) | bewiesen (Intervall), bedingt auf (E) | Im Y/2π ∈ {0.31, 0.33, 0.29, 0.43} (e), {−0.17, −0.16, −0.19, −0.08} (2) | – |
| 10 | Diese Zweigzahl ist eine **Konvention**: h aus lib3 hängt von der Abbruchregel n_h ab (n_h = 23 für e, 37 für 2); ändert man n_h, springt h um 2πi/ℓ. Inhaltlich zählt nur S4 (ohne Logarithmus) | negativ (Klarstellung) | Abschnitt 3 | – |
| 11 | Teil A rigoros, Basis e: Σ_{k≥2} k\|c_k\| r^{k−1} ≤ 0.06098 < A = 1.0918 (r = 0.1) ⇒ sexp schlicht auf \|z\| < 0.1; \|ε₀\| ≤ 0.0277, \|ε₁\| ≤ 0.0339, \|sexp′/A − 1\| ≤ 0.0559 | bewiesen (Intervall), bedingt auf (E2) | `teil_a.py`, 16384 Bälle | bestätigt 0.034 / 0.056; ε₀ nur ≤ 0.0277 statt 0.027 zertifiziert |
| 12 | Teil A rigoros, Basis 2: Summe ≤ 0.00463 < A = 0.88936; \|ε₀\| ≤ 0.0030, \|ε₁\| ≤ 0.0481, \|sexp′/A − 1\| ≤ 0.0052 | bewiesen (Intervall), bedingt auf (E2) | `teil_a.py` | neu |
| 13 | Toleranz von (E2): die Schlichtheit hält, solange die Taylor-Koeffizienten bis auf 4.39 (e) bzw. 3.77 (2) stimmen | bewiesen | `teil_a.py` | neu |

## 1. Logik des Zertifikats (warum ein Ball reicht und Ψ entfällt)

Bezeichnungen wie in luecken/BERICHT_LUECKEN_ANALYTISCH.md: λ = L ln b, ℓ = Log λ (Hauptzweig; für e ist ℓ = L),
F(s) = Ψ(e^{ℓs}), σ = Schröder-Funktion mit σ(w) = λσ(log_b w), σ′(L) = 1. Lemma U: σ ist schlicht auf ℍ und Ψ∘σ = id auf ℍ.

**Behauptung (\*):** P := e^{ℓ(z₀+θ_K(z₀))} liegt in σ(ℍ). Dann gilt S4 an z₀:
P = σ(w′) mit w′ ∈ ℍ ⇒ Ψ(P) = w′ (U(c)); nach V1 ist Ψ(P) = F(z₀+θ_K(z₀)) = sexp z₀ ⇒ σ(sexp z₀) = σ(w′) = P. ∎
Mit Lemma H (Ψ lokal umkehrbar) und Teil A (sexp(D) ⊂ ℍ, damit σ∘sexp auf D stetig) folgt S4 auf ganz D.

**Beweis von (\*):** Wähle einen beliebigen Testpunkt w₁ ∈ ℍ (hier per Gleitkomma-Newton mit σ(w₁) ≈ P; keine Annahme, w₁ wird
als exakte Zahl weiterverwendet). σ ist auf B(w₁, Im w₁) ⊂ ℍ schlicht, also enthält σ(B(w₁, Im w₁)) nach Koebe-¼ die Scheibe
B(σ(w₁), ρ), ρ = |σ′(w₁)| Im w₁/4. Es genügt also |P − σ(w₁)| < ρ.

**Rigoroses σ.** σ(w₁) = λⁿσ(v_n), v_{k+1} = log_b v_k (acb-Bälle, alle in ℍ geprüft). Nahe L: f(t) = σ(L+Rt)/R ist auf der
Einheitsscheibe schlicht (R = 1 < Im L), f(0) = 0, f′(0) = 1. Mit de Branges (|a_k| ≤ k):
|σ(L+u) − u| ≤ R(t/(1−t)² − t), |σ′(L+u) − 1| ≤ (1+t)/(1−t)³ − 1, t = |u|/R.
n = 450 (e) bzw. 700 (2) Schritte; der Rest ist dann ≈ 10⁻⁶¹. L selbst: Intervall-Newton für b^w = w (Einschluss und Eindeutigkeit
im Ball vom Radius 10⁻⁴⁰).

**Toleranz.** Mit θ_K(z₀) = θ_num(z₀) + δ, |δ| ≤ τ ist |P_K − P_num| ≤ |P_num|(e^{|ℓ|τ} − 1). Also genügt
|P_num − σ(w₁)| + |P_num|(e^{|ℓ|τ} − 1) < ρ; daraus τ_max = log1p((ρ − dist)/|P|)/|ℓ| (rigoros nach unten gerundet).

**Nebenprodukt sexp(z₀).** g(ζ) = (σ(w₁+r₁ζ) − σ(w₁))/(r₁σ′(w₁)) ist normiert schlicht; aus |g(ζ*)| ≤ ε und Koebe-Wachstum
|g(ζ)| ≥ |ζ|/(1+|ζ|)² folgt |sexp z₀ − w₁| ≤ r₁t*, t*/(1+t*)² = ε.

## 2. Zahlen (`daten/zertifikat.json`)

| Basis | z₀ | \|P − σ(w₁)\| ≤ | ρ ≥ | τ_max ≥ | \|θ₈₀ − θ₆₀\| | \|w₁ − sexp_Taylor\| |
|---|---|---|---|---|---|---|
| e | 0.05i | 2.2·10⁻⁶¹ | 0.1212 | 0.0273 | 2.1·10⁻⁵⁹ | 2.5·10⁻⁶³ |
| e | 0.07e^{iπ/3} | 2.1·10⁻⁶¹ | 0.1164 | 0.0268 | 2.1·10⁻⁵⁹ | 3.1·10⁻⁶³ |
| e | 0.08e^{3iπ/4} | 1.8·10⁻⁶¹ | 0.0958 | 0.0234 | 2.1·10⁻⁵⁹ | 3.0·10⁻⁶³ |
| e | 0.5+0.05i | 1.9·10⁻⁶¹ | 0.0382 | 0.0092 | 2.1·10⁻⁵⁹ | 3.5·10⁻¹⁷ (Taylor um 0.4) |
| 2 | 0.05i | 1.6·10⁻⁶¹ | 0.0980 | 0.0273 | 1.6·10⁻⁵⁹ | 1.3·10⁻⁶³ |
| 2 | 0.07e^{iπ/3} | 1.5·10⁻⁶¹ | 0.0945 | 0.0269 | 1.6·10⁻⁵⁹ | 1.5·10⁻⁶³ |
| 2 | 0.08e^{3iπ/4} | 1.4·10⁻⁶¹ | 0.0789 | 0.0234 | 1.6·10⁻⁵⁹ | 3.1·10⁻⁶³ |
| 2 | 0.5+0.05i | 1.4·10⁻⁶¹ | 0.0306 | 0.0092 | 1.6·10⁻⁵⁹ | 2.2·10⁻¹⁷ (Taylor um 0.4) |

Die Distanz 10⁻⁶¹ ist durch die Zahl der Log-Schritte begrenzt, nicht durch die Daten (θ_num wird aus den Datei-a_n exakt summiert).
Gegenprobe der Einstellung: mit prec 320 und n = 300 (e) bzw. 420 (2) kommt 1.1·10⁻³³ bzw. 1.2·10⁻³² heraus, ebenfalls weit unter ρ;
zu viele Schritte bei zu wenig Bits sprengen den Ball (n = 900 bei 320 Bit) — dann schlägt die Prüfung fehl, statt Falsches zu behaupten.

## 3. Was bleibt vorausgesetzt (ehrlich)

- **(V1) Kneser-Darstellung**: sexp holomorph auf ℂ∖(−∞,−2], reell, sexp = F(z+θ) mit θ periodisch und holomorph in ℍ (Kneser 1950,
  Paulsen–Cowgill 2017). Unverändert.
- **(E) Motor, Fourier-Seite**: |θ_K(z₀) − Σ_{n≤4096} a_n q₀ⁿ| ≤ 0.027 (mod 2πi/ℓ), a_n aus etappe2/daten/a_{e,2}_80_0.003_16384.json.
  Belegt durch zwei Einstellungen auf 10⁻⁵⁹; die Annahme ist also rund 57 Grössenordnungen schwächer als das, was der Motor liefert.
  Nicht rigoros ist: dass der Motor (fatou.gp) wirklich Knesers Lösung berechnet und die FFT auf der Linie Im z = 0.003 treu ist.
- **(E2) Motor, Taylor-Seite**: |c_k(Kneser) − c_k(Datei)| ≤ 10⁻³⁰ für k ≤ 127 und |c_k| ≤ 1 für k ≥ 128 (Dateien etappe3/daten/taylor_130_256.json,
  taylor_b2_80_256.json). Für die Schlichtheit allein reicht Fehler ≤ 3.7.
- **Zitierte Sätze**: de Branges (Bieberbach-Vermutung, 1985), Koebe-¼ und Koebe-Wachstum, Noshiro–Warschawski, Hurwitz, Schwarz–Pick;
  für die Koeffizienten zusätzlich Flajolet–Odlyzko (wie bisher).
- **Kleine Korrektur**: In luecken/BERICHT_LUECKEN_ANALYTISCH.md, Satz Z, steht |ε₀| ≤ 0.027 (Gittermaximum 0.02688). Rigoros zertifiziert ist
  hier nur |ε₀| ≤ 0.0277. Wahrscheinlich gilt 0.027, aber für die Konstanten in (M3) sollte man 1.028 statt 1.027 schreiben. Alle Aussagen
  bleiben gleich, nur eine Konstante c ändert sich minimal.
- **Zweigzahl (Befund 10)**: „Hauptzweig h = lib3.Basis.h“ ist keine innere Eigenschaft. Für die Koeffizientenformel ist das egal (ein
  anderer Zweig ändert nur b₀, wie im Transferbeweis schon gesagt). Für die Lokalform (iii) sollte man sagen: „für den Zweig von h, der
  an z₀ = 0.05i den Datenwert θ(z₀) liefert“, oder einfach „mod 2πi/ℓ“.

## 4. Vorschläge für Grapheinträge

- kn.satz.theta_singularitaet — Statusänderung: bedingt → bewiesen (rechnergestützt, Intervallarithmetik; Motordaten nur auf 0.027 nötig); bedingungen_text ersetzen durch „Motor-Annahme (E) mit τ = 0.027 und (E2); Basis e und 2“.
- kn.satz.theta_koeffizienten — Statusänderung: bedingt → bewiesen (rechnergestützt), gleiche Bedingungen.
- kn.lemma.zweig_ein_ball — Lemma, bewiesen: „S4 an z₀ ⟸ e^{ℓ(z₀+θ(z₀))} ∈ σ(ℍ); geprüft über Koebe-¼ um einen freien Testpunkt, ohne Ψ“.
- kn.methode.sigma_de_branges — Methode, bewiesen: rigoroses σ und σ′ über λⁿσ(log_bⁿ w) und de Branges nahe L (Fehler ~ |λ|⁻ⁿ).
- kn.beob.teil_a_basis2 — Beobachtung/rechnergestützt bewiesen: sexp₂ schlicht auf |z| < 0.1, ε₀ ≤ 0.0030, ε₁ ≤ 0.0481, |sexp′/A − 1| ≤ 0.0052.
- Korrektur in luecken-Bericht (kein Grapheintrag): |ε₀| ≤ 0.0277 statt 0.027.

## 5. Offene Fragen

- Lässt sich (E) selbst rigoros machen? Eine Möglichkeit: sexp(iy) für y ≥ 0.05 über die θ-Reihe mit rigorosem Ψ einschliessen und zeigen,
  dass der Weg in ℍ bleibt; dann folgt die Übereinstimmung der Hochhebungen bei y → ∞ (beide → 0, Ψ dort injektiv) ganz ohne Punktprüfung.
  Das ersetzt aber nur eine Motor-Annahme durch eine andere, solange Kneser selbst nicht rigoros berechnet ist.
- Die Schranke ε₀ ≤ 0.027 mit zentrierter Form (Taylor-Verschiebung je Bogen) statt reiner Ballauswertung zertifizieren.
- Dasselbe Ein-Ball-Zertifikat für die Kneser-artige Pentation (kn.beob.einheitliche_raenge): dort fehlt Lemma U (Schlichtheit), das ist die eigentliche Hürde.

## Reproduzieren

```
C:/Users/janis/miniconda3/python.exe zertifikat.py   # Befunde 1–10, < 2 min, schreibt daten/zertifikat.json
C:/Users/janis/miniconda3/python.exe teil_a.py       # Befunde 11–13, ~30 s, schreibt daten/teil_a.json
```
Kein Motorlauf; nur Lesezugriff auf forschung/identitaeten_2026_10_14/etappe2/daten und etappe3/daten.
