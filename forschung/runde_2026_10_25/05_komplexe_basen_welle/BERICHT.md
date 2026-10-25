# Runde 25.10.2026 — 05 Knesers Welle für komplexe Basen (Paulsen–Cowgill)

Zeit: 13:23–13:52 (aus `date`). Motor: fatou.gp-Fork über eigenes `kern.py` (eigener Zwischenspeicher `cache/`).

**Kurz in einfacher Sprache.** Für komplexe Basen ausserhalb des Shell–Thron-Gebiets gibt es zwei Wellen: eine obere (zum oberen Fixpunkt L₊)
und eine untere (zum unteren Fixpunkt L₋). Sie sind nicht mehr Spiegelbilder voneinander. Die scharfe Stelle der Welle bei den ganzen Zahlen
bleibt genau dort, weil sexp(−1) = 0 für jede Basis gilt. Das Koeffizientengesetz aₙ = [qⁿ] h((c + log(1−q))/ln b) gilt weiter, für beide Wellen, mit
c = log(±iA/(2π ln b)). Neu und unerwartet: geht man zum Shell–Thron-Rand in der Richtung, in der der Multiplikator gegen i läuft (Drehung um 1/4),
dann wandert bei |λ₊| ≈ 1.05 eine **zusätzliche Singularität in die obere Halbebene** (bei z* ≈ 0.034 + 0.0104i mod 1). Dann wachsen die Koeffizienten
geometrisch und das Gesetz bricht. In einer anderen Richtung (λ → e^i) mit gleichem |λ₊| passiert das nicht.

## Befunde

| Nr | Aussage | Status | Beleg / Genauigkeit | Neuheit |
|---|---|---|---|---|
| K1 | Für jede Basis b (ausserhalb Shell–Thron) mit sexp = F₊(z + θ_U(z)) im oberen Bereich hat die q-Reihe von θ_U den Radius **R ≤ 1**, und z = −1 (q = 1) ist singulär oder liegt ausserhalb. Grund: sexp(−1) = 0 exakt, F₊(s) = Ψ₊(λ₊^s) ist ganz und lässt 0 aus (Ψ(λu) = b^{Ψ(u)} ≠ 0), auch für komplexes b. Dasselbe für θ_L mit F₋. | **bewiesen** (gegeben die Darstellung sexp = F₊(z+θ_U) für grosses Im z) | Beweis unten | Übertragung von kn.satz.theta_singularitaet auf komplexe b; neu |
| K2 | sexp hat für komplexe Basen keine weiteren Nullstellen dort, wo sexp(z) = b^{sexp(z−1)} gilt (b^w ≠ 0). Die Singularität q = 1 **wandert nicht**. Newton-„Nullstellen“ nahe 3 bzw. 5 (|sexp| ≈ 10⁻¹⁴ bei dps 30) sind Fehlkonvergenzen. | **bewiesen** (ein Satz) / **negativ** für die Newton-Funde | e-Schritt in diesem Bericht | — |
| K3 | **Koeffizientengesetz der oberen Welle bleibt:** aₙ = [qⁿ] h₊((c_U + log(1−q))/ln b) · (1 + O(1/n)), c_U = log(iA/(2π ln b)), A = sexp′(0) komplex. Relativer Rest n·\|aₙ−bₙ\|/\|aₙ\| ≈ 0.11 (r=2), 0.12 (r=1.5), 0.16 (r=1.2) bei n = 200; ebenso φ = 1: 0.17 (r=1.1), 0.13 (r=1.2). | **numerisch** | dps 30, N = 512, Linie y = 0.05; Gegenprobe y = 0.08: aₙ auf 10 Stellen gleich bis n = 96 (r = 2) | neu |
| K4 | **Exponent genau 1 bleibt:** n·ln n·\|aₙ\| ≈ 0.139 (r=2), 0.116 (r=1.5), 0.080 (r=1.2), 0.095 (φ=1, r=1.1) bei n = 200, langsam veränderlich; lokaler Exponent 1.18–1.25 sinkt wie im reellen Fall. | **numerisch** | wie K3 | neu |
| K5 | **Untere Welle:** θ_L = Σ dₙ pⁿ, p = e^{−2πiz}, erfüllt das Gesetz mit h₋ und c_L = log(−iA/(2π ln b)); n·rel ≈ 0.128 bei n = 64 (r = 2). | **numerisch** (nur eine Einstellung: dps 60, y = 0.25, N = 256; Gegenprobe fehlt) | `welle_..._60_0.25_256_m0_L.json` | neu |
| K6 | Die Wellen sind nicht mehr konjugiert: r = 2: a₀^U = 22.911 + 3.524i, a₀^L = 10.370 − 4.063i; \|a₁^U\| = 0.08979, \|a₁^L\| = 0.09317; arg aₙ^U ≈ −1.47…−1.27, arg dₙ ≈ +1.52…+1.43. | numerisch | wie K3/K5 | — |
| K7 | **Neue Singularität in der oberen Halbebene nahe am Rand:** ln b = 1.1·i·e^{−i} (b = 2.0907 + 1.4130i, λ₊ = −0.0487 + 1.0465i): die Reste rₙ = aₙ − bₙ wachsen wie ρⁿ, ρ → 1.0432 − 0.2253i (n = 199), also Singularität von θ_U bei **z* ≈ 0.0339 + 0.0104i (mod 1)** (Prony 1–3 Terme: Im 0.0103–0.0108; Quotient steigt noch langsam ⇒ kein Pol, eher Verzweigung). Radius der q-Reihe R ≈ 0.936 < 1. | **numerisch** | zwei Linien y = 0.05 und 0.08: aₙ auf 10 Stellen gleich bis n = 128, also kein Rechenfehler; Rauschen ≤ 10⁻³ bei n = 200 gegen \|rₙ\| = 70 | neu |
| K8 | Bei gleichem \|λ₊\| ≈ 1.047 in der Richtung λ → e^{i} (φ = 1 rad, irrationale Drehung) tritt **keine** solche Singularität auf (Gesetz hält bis n = 200). Bei φ = π/2, r = 1.2 (\|λ₊\| = 1.091) auch nicht. | numerisch | K3-Läufe | neu |
| K9 | a₀ (der Mittelwert von θ_U) wächst zum Rand hin: 22.9 (r=2), 38.7 (1.5), 82.7 (1.2), 150.9 (1.1) bei φ = π/2; 90.9 (1.2), 166.4 (1.1) bei φ = 1. Etwa ∝ 1/log\|λ₊\|: a₀·ln\|λ₊\| ≈ 6.9, 7.2, 7.2, 7.0 (φ = π/2) und 7.6, 7.6 (φ = 1). | numerisch / Vermutung (Skalierung) | 4 + 2 Basen | neu |
| K10 | Zusammenhang zu K7 (Vermutung): für diese Basis kehrt die Bahn 0 → 1 → b → b^b → sexp(3) fast zu 0 zurück (\|sexp(3)\| = 0.0373, 4 Schritte = Nenner der Drehzahl 1/4); die Bahn von sexp(z*) trifft nach 3 Schritten 0.049 − 0.010i. sexp(z*) ≈ 1.034 + 0.024i ist weder 0, 1, b noch L±. | **Vermutung** | eine Basis | — |
| K11 | Bei r = 1.05 (φ = π/2) bricht mein h₊ (Hauptzweig-Iteration von log_b) auf der Linie ab: Bahnen erreichen L₊ nicht. Ebenso h₋ auf Im z = −0.05 für r = 2 (beide Fixpunkte haben Zweig k = 0, die Iteration landet bei L₊). Verschieben nach links (sexp(z−m)) hilft nicht (Motorwerte bei Re z < −2 springen zwischen den Becken). | **negativ** | Fehlermeldungen in den Läufen | — |

## Beweise und Details

**K1.** Sei R > 1. Dann ist θ_U auf einem Streifen um die reelle Achse holomorph. Beide Seiten von sexp(z) = F₊(z + θ_U(z)) sind auf der
zusammenhängenden Menge (Streifen ∪ oberer Bereich) holomorph und stimmen oben überein, also auch bei z = −1. Daraus folgt F₊(−1 + θ_U(−1)) = 0.
Aber F₊ = Ψ₊(λ₊^s) mit Ψ₊ ganz und Ψ₊(λu) = b^{Ψ₊(u)} ≠ 0, also lässt F₊ die 0 aus. Widerspruch, also R ≤ 1. Für θ_L geht es genauso mit
F₋ und p = e^{−2πiz}. Benutzt wird nur ln b ≠ 0 und |λ±| > 1. Ψ′ ≠ 0 gilt auch, mit demselben Häufungsargument wie bei kn.satz.theta_singularitaet,
und dort spielt komplexes b keine Rolle.

**K2.** Wo sexp(z) = b^{sexp(z−1)} = exp(ln b · sexp(z−1)) gilt, ist sexp(z) ≠ 0. Die einzige Nullstelle ist z = −1, wo sexp(z−1) = sexp(−2) singulär ist.
Das Bild „die Singularitäten wandern“ stimmt also für die q = 1-Stelle nicht. Neue Singularitäten (K7) kommen nicht von Nullstellen.

**K3/K5 Herleitung.** Wie im reellen Fall: θ_U(z) = 2 + h₊(log_b log_b sexp z) − z, und log_b(sexp z) ≈ Az/ln b.
Mit z ≈ i(1−q)/(2π) für q → 1 aus |q| < 1 ergibt sich log_b log_b sexp z ≈ (log(iA/(2π ln b)) + log(1−q))/ln b, mit Hauptzweigen.
Unten gilt z ≈ −i(1−p)/(2π), daher das Minuszeichen in c_L. Der Rest ist relativ O(1/n), so wie die δ₁-Korrektur im reellen Fall.

**K7 Ort.** Quotienten r_{n+1}/r_n: |·| = 1.0596 (n=80), 1.0642 (120), 1.0659 (160), 1.0673 (199); Argument fest bei −0.2132 rad ⇒ Re z* = 0.0339.
Da |ρ| > 1, liegt z* in der **oberen** Halbebene, also im angenommenen Holomorphiegebiet von θ_U. Wegen F₊′ ≠ 0 muss sexp(z*) über einem
asymptotischen Wert von F₊ liegen (Iversen). Die bekannten asymptotischen Werte sind 0, 1, b, …, L₊. sexp(z*) passt zu keinem davon. Das ist offen.

## Skripte (jede Zahl reproduzierbar)
- `kern.py`: Motor-Zugang (eigener Speicher), Fixpunkt mit Newton, Schröder-Reihe, h, FFT.
- `e1_basis.py "<basis>" <dps>`: Grunddaten (zeigte: b = 2+i liegt **innerhalb** von Shell–Thron, |λ₊| = 0.9877).
- `e2_welle.py "<basis>" <dps> <y> <N> [m] [U|L|UL]`: Koeffizienten und Modellvergleich → `daten/welle_*.json`.
  Läufe: `exp(r*(sin(1)+I*cos(1)))` mit r = 2, 1.5, 1.2, 1.1 (dps 30, y = 0.05, N = 512, U), r = 1.1 auch y = 0.08, r = 2 L mit dps 60, y = 0.25, N = 256;
  `exp(r*exp(I)*exp(-exp(I)))` mit r = 1.2, 1.1.
- `e3_prony.py <welle-json>`: Ort der Zusatzsingularität.

## Vorschläge für Grapheinträge
- kn.satz.komplex_radius — „Für komplexe Basen hat die q-Reihe der Kneser-Welle Radius ≤ 1, die Stelle q = 1 bleibt“ — satz — bewiesen
- kn.beob.komplex_koeffizientengesetz — „Koeffizientengesetz aₙ = [qⁿ]h₊((c+log(1−q))/ln b) gilt für komplexe Basen, oben und unten (c = log(±iA/2π ln b))“ — beobachtung — numerisch
- kn.beob.komplex_wellen_unsymmetrisch — „Obere und untere Welle sind für komplexe b verschieden (a₀, |a₁|, arg)“ — beobachtung — numerisch
- kn.beob.rand_zusatzsingularitaet — „Nahe dem Shell–Thron-Rand bei λ → i wandert eine Singularität von θ_U in die obere Halbebene (z* ≈ 0.034+0.0104i bei r = 1.1)“ — beobachtung — numerisch
- kn.vermutung.a0_rand — „a₀ der oberen Welle wächst wie ≈ (7…7.6)/ln|λ₊| zum Rand“ — vermutung — numerisch, 6 Basen

## Offene Fragen
1. Woher kommt z* (K7)? Ist sexp(z*) ein asymptotischer Wert von F₊ aus der Nähe-Parabolik (λ₊ ≈ i), etwa ein Punkt eines nahen 4-Zyklus oder sein Urbild?
   Mein Newton fand innerhalb 0.3 um L₊ keinen 4-Zyklus. Prüfen: z* für r = 1.15, 1.12, 1.08 verfolgen. Wann überquert es Im z = 0?
2. Gibt es solche Zusatzsingularitäten bei jeder rationalen Randrichtung p/q, und ab welchem Abstand? Bei φ = 1 rad gab es keine (K8).
3. Für r → 1 braucht man ein h₊, das nicht über die Hauptzweig-Iteration läuft (K11), zum Beispiel Fatou-Koordinaten.
4. Für die untere Welle fehlt die zweite Einstellung (K5). Ausserdem fehlt ein h₋, das nahe der reellen Achse funktioniert.
