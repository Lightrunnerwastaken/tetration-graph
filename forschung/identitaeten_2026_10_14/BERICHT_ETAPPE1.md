# Identitätssuche, Etappe 1 (14.10.2026)

## Kurz gesagt

- **Konstanten:** 50 Zahlen aus der Tetration wurden auf rund 160 sichere Stellen berechnet: Kneser-Werte zu den Basen e und 2, Taylor-Koeffizienten, Integrale, die exponentielle Fakultät, Fatou-Konstanten und Fourier-Koeffizienten von Knesers θ. Keine dieser Zahlen hat eine einfache Formel mit π, e, γ, ln 2, ζ(3), Catalan, Γ(¼), Γ(⅓), Ω, dem Fixpunkt L und ähnlichen Konstanten. Jede solche Formel bräuchte riesige ganze Koeffizienten (Schranken siehe Tabelle 2). Die Methode selbst funktioniert: drei eingebaute Kontrollformeln wurden alle wiedergefunden.
- **Beinahe-Treffer:** Mehrere Zahlen sehen auf 4–5 Stellen bekannt aus, zum Beispiel slog′(0) = 0.915946 und Catalan G = 0.915966. Alle fallen bei voller Genauigkeit weg. Dass drei Grössen nahe bei ½ liegen, kommt daher, dass sexp auf [−1, 0] fast eine Gerade ist.
- **Kombinatorik:** Hier gibt es echte, bewiesene Zusammenhänge.
  1. Die Störungsreihe g_n(y) aus e81 zählt markierte Bäume einer bestimmten Höhe. Ihr höchster Anteil ist genau die Tangensfunktion (Tangentenzahlen). Für den zweithöchsten Anteil gibt es eine geschlossene Formel mit sec² und ln cos.
  2. Die gebrochenen Iterierten von eˣ−1 sind Werte des Zeta-Polynoms des Partitionsverbands. Das ist klassisch.
  3. Für x·eˣ entsteht das entsprechende Dreieck. Es steht nicht in der OEIS.
- **Neu und offen:** Knesers θ hat Fourier-Koeffizienten, die nur **algebraisch** abfallen (ungefähr n^(−1.2), mit langsam fallendem Exponenten). θ ist also an der reellen Achse singulär. Das ist der natürliche Einstieg für Etappe 2.

Alle Zahlen sind reproduzierbar mit den Skripten in diesem Ordner. Interpreter: `C:/Users/janis/miniconda3/python.exe`, Ausgaben liegen in `daten/`.

| Schritt | Aufruf | Dauer |
|---|---|---|
| Kneser-Konstanten | `a_konstanten.py 120 / 160 / 200` | je < 5 min |
| Fatou-Konstanten | `a_fatou.py` | < 1 min |
| θ-Fourier | `a_theta.py 120 / 160 / 200` | je < 1 min, Zwischenspeicher `daten/motor_cache.json` |
| Relationssuche | `a_relationen.py` | 1.5 min |
| Kombinatorik | `b_kombinatorik.py` (≈ 8 min), `b_gewichte.py`, `b_baeume.py`, `b_diagonalen.py` | |

Motor: fatou.gp-Fork über `motor.py`. Der Code nutzt bewusst nicht den 190-MB-Zwischenspeicher von `experimente/gemeinsam.py`.

**Genauigkeit:** Die Stufen dps 160 und 200 stimmen bei allen Kneser-Grössen auf ≤ 4·10⁻¹⁵⁹ überein (bei a₁₅: 4·10⁻¹⁵⁹). Die Stufen 120 und 160 stimmen auf ≈ 10⁻¹³⁰ überein. Verwendet wird **P = 160 Stellen**. Die Fatou-Konstanten sind auf 185 Stellen genau (Abstiegstiefe 4000 gegen 8000 und Reihenlänge 60 gegen 70). Die θ-Koeffizienten sind auf ≥ 162 Stellen genau (Stufe 160 gegen 200), und a_n hängt bis 10⁻²⁰³ nicht von der Abtasthöhe y ab.

**Nebenbefund zum Motor:** Der Fork bricht bei dps ≥ 60 für |Im z| ≳ 0.3 mit `exp: overflow` bzw. `log: domain error` in `invabel` ab. Bei dps 40/50 läuft er. Für |Im z| ≤ 0.2 geht alles. Deshalb laufen alle Cauchy-Kreise mit Radius 0.2.

## Tabelle 1: Befunde (positiv)

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit |
|---|---|---|---|---|
| K1 | **Höchster Anteil der Störungsreihe = Tangens.** Für G(y) = e^{G(y−1)} − 1 + u, G = Σ g_n(y) uⁿ gilt Σ_n LC(g_n)·uⁿ y^{2n−1} = √(2u)·tan(y√(u/2)). Also ist der Leitkoeffizient LC(g_n) = 2^{1−n} T_n/(2n−1)! mit den Tangentenzahlen T_n = 1, 2, 16, 272, … | **bewiesen** (siehe Beweis B1); exakt geprüft für n ≤ 12 | `b_kombinatorik.py`, `b_gewichte.py` | im Geist klassisch (parabolischer Grenzfall x + x²/2 + ε, „parabolic implosion“); in dieser Form für die Tetrations-Störungsreihe wohl nicht notiert |
| K2 | **Zweithöchster Anteil, geschlossen:** Mit H(y) = G(y−1) und θ = y√(u/2) gilt H₁ = (u/3)(sec²θ·ln cos θ + 1 − 4 sec²θ). Daraus folgen die Koeffizienten von y^{2n−2} in g_n(y−1): −1, −3/4, −13/48, −13/160, −99/4480, … | **bewiesen** (B1, zweite Schicht); exakt geprüft für n ≤ 12 | `b_gewichte.py` | wahrscheinlich neu in dieser Form |
| K3 | **Baumdeutung:** n!·g_n(k) zählt Bäume mit n markierten Blättern und Höhe ≤ k. Innere Knoten sind nichtleere Mengen von Kindern, Einzelkinder sind erlaubt. Grund: G(k) = f^{∘k}(0) mit f = u + eˣ − 1 ist eine Spezies-Komposition. Folgerungen: n!g_n(2) = 1, n!g_n(3) = Bell(n+1), g_n(−1) = −1/n. Das Newton-Dreieck N(n,j) ist ganzzahlig, seine Diagonale N(n,2n−1) = n!·T_n/2^{n−1} (1, 2, 24, 816, 59520) folgt aus K1. | bewiesen (EGF-Argument); Werte exakt geprüft | `b_kombinatorik.py`, `b_baeume.py` | Bell-Fall trivial. Newton-Dreieck, Zeilensummen 1, 3, 70, 4640, 668316, Diagonale sowie n!g_n(4) = 4, 14, 87, 658, 5801 und n!g_n(5): **nicht in der OEIS** |
| K4 | **Gebrochene Iterierte von eˣ−1 = Zeta-Polynom des Partitionsverbands:** n![xⁿ](eˣ−1)^{∘t} = Σ_j B(n,j)·C(t,j) = Z(Π_n, t). Dabei ist B(n,j) die Zahl der echten Ketten der Länge j. Diagonale n!(n−1)!/2^{n−1} (maximale Ketten). Bei t = −1 kommt die Möbius-Zahl (−1)^{n−1}(n−1)! heraus. Halbiterierte = Z(Π_n, ½). | bewiesen (Spezies + Binomialentwicklung), exakt geprüft für n ≤ 16 | `b_diagonalen.py` | **klassisch**: Dreieck A008826, Halbiterierte A052104, Iterierte A000110/A000258/A000307 |
| K5 | **Gleiches für x·eˣ (Produkt-Tetration):** Binomialdreieck 2; 3, 12; 4, 84, 144; 5, 510, 2760, 2880; … mit Diagonale n!(n−1)!. Bei t = −1 kommt (−n)^{n−1} heraus, die Koeffizienten von Lambert-W, wie bei der Möbius-Zahl der „pointed partitions“. | Identitäten bewiesen bzw. für n ≤ 16 geprüft; die Deutung als Ketten punktierter Partitionen (Chapoton–Vallette) ist nicht nachgeprüft | `b_diagonalen.py` | Iterierte ×2 = A185298 (bekannt). Dreieck, Zeilensummen 2, 15, 232, 6155, ×3/×4-Iterierte, Julia-Koeffizienten (2, −3, 10, −50, 321, −2422) und Halbiterierte: **nicht in der OEIS** |
| K6 | **Fatou-Koordinaten exakt:** α(x·eˣ) = α(x)+1 mit α = −1/x + ½ ln x − x/6 + x²/16 − 19x³/540 + x⁴/48 − 41x⁵/4200 + 37x⁶/103680 + … Für eˣ−1: α = −2/x + ⅓ ln x − x/36 + x²/540 + x³/7776 − 71x⁴/435456 + … | bewiesen (Koeffizientenrekursion) | `a_fatou.py` → `daten/fatou.json` (70 Koeffizienten) | eˣ−1: aus der Literatur bekannt (Écalle). Beide Brüche nicht in der OEIS |
| K7 | **0!₂ ist eine Ableitung:** lim_{x↓0} x!₂ = F′(1) mit F(x) = x!₂ (Begründung: 1!₂ = 1^{0!₂} und L'Hôpital), F′(1) = 0.5730489722129730653505964537816700999925… Ausserdem D = slog_e(262144·ln 5 + ln ln 6) − 4 bis auf 10^{−180000}, D = −1.0552140087126523890449559784636570718988… | Beweis elementar; 160 Stellen | `a_konstanten.py` | Präzisiert e82 (dort 0.573049 durch Extrapolation) |
| K8 | **Kneser-θ zur Basis e:** sexp(z) = F(z+θ(z)) mit F = σ⁻¹(e^{Lz}), σ′(L) = 1. Fourier-Koeffizienten: a₁ = 0.0038515787354112044204 − 0.0905647779541168204455 i, a₂ = 0.0021009167869028148730 − 0.0379562481524750855419 i, a₃ = 0.0014475477584217851203 − 0.0228172971300237550290 i, a₀ ≡ −0.48841508844379660331 − 0.92230686292605693600 i (mod 2πi/L). arg a_n → nahe −π/2. \|a_n\| fällt **algebraisch**: lokaler Exponent 1.254 (n=2), 1.241 (8), 1.208 (32), 1.192 (64), und n·\|a_n\| fällt langsam. | numerisch, 162 Stellen (a₁–a₃), a₆₄ auf 10⁻⁹⁷ | `a_theta.py`, `daten/theta_*.json` | nach kurzer Suche keine veröffentlichten hochgenauen θ-Koeffizienten bekannt |

### Beweis B1 (K1, K2)

Wir setzen H(y) = G(y−1). Dann gilt H(y+1) − H(y) = Σ_{m≥2} H^m/m! + u. Wir geben y das Gewicht 1 und u das Gewicht −2. Ein Term uⁿy^k hat dann das Gewicht k − 2n, also ist g_n vom Grad ≤ 2n−1 genau dann, wenn H nur Gewichte ≤ −1 hat (Induktion über n). Wir zerlegen H = H₀ + H₁ + … nach Gewichten −1, −2, … Die Differenz ist H(y+1) − H(y) = H′ + ½H″ + … Jede weitere Ableitung senkt das Gewicht um 1.

- **Gewicht −2:** H₀′ = H₀²/2 + u mit H₀(0) = 0. Also ist H₀ = √(2u) tan θ, und eine Reihenentwicklung ergibt K1.
- **Gewicht −3:** H₁′ + ½H₀″ = H₀H₁ + H₀³/6. Mit H₀″ = H₀H₀′ wird daraus H₁′ − H₀H₁ = −H₀³/12 − uH₀/2. Die Anfangsbedingung kommt aus H(0) = G(−1) = log(1−u): H₁(0) = −u. Wir lösen mit dem integrierenden Faktor cos²θ. Dabei ist ∫(sin³/cos)dθ = −ln cos θ + cos²θ/2. Das ergibt K2.

Beide Schichten stimmen mit den exakten Polynomen überein (`b_gewichte.py`: True/True bis n = 12). ∎

## Tabelle 2: Negative Ergebnisse (Relationssuche, `daten/relationen.json`)

**Methode:** LLL (python-flint) auf dem Gitter [I | round(10^{P_find}·x_i)] mit P_find = P − 20. Ein Treffer zählt nur, wenn der Rest bei voller Genauigkeit P < 10^{−(P_find+15)}·‖m‖ ist. Die Schranke ist streng für das LLL-Gitter: ‖m‖₂ ≥ min‖b_i*‖/√(1+n/4). „10^5.6“ heisst also: **Jede ganzzahlige Relation hat einen Koeffizientenvektor mit euklidischer Norm ≥ 4·10⁵.**

**Basen:**

- **lin (23 Elemente):** 1, π, π², e, eᵉ, γ, ln 2, ln 3, ln π, ζ(3), G, Γ(¼), Γ(⅓), Ω = W(1), Re L, Im L, \|L\| sowie ln γ, ln ζ(3), ln G, ln Γ(¼), ln Γ(⅓). Bei Basis 2 kommen L₂ (Re, Im, Betrag), ln ln 2 und 1/ln 2 dazu. Dabei gilt ln\|L\| = Re L, arg L = Im L und ln Ω = −Ω, diese Grössen sind also schon enthalten.
- **klein:** vier Fünferbasen: {1, π, e, γ, ln 2}, {1, π, Re L, Im L, \|L\|}, {1, Ω, e, π, γ} und {1, π, π², ζ(3), G}.
- **log:** ln\|x\| gegen 14 Logarithmen (Produkte von Potenzen): ln 2, ln 3, ln 5, ln π, 1, ln γ, ln ζ(3), ln G, ln Γ(¼), ln Γ(⅓), −Ω, Re L, ln Im L, e.
- **poly:** algebraisch vom Grad 2, 4 und 6.

**Positivkontrollen** (alle gefunden): 3π − 2γ + 5 ln 2/7 + Im L/11 (Norm 290), Γ(¼)²e³/π^{3/2} (log), Wurzel von x⁵ − x − 1 (poly6).

| Zielgrösse (Wert, 40 Stellen in `daten/relationen.json`) | P | lin | klein (min) | log | poly2 | poly6 |
|---|---|---|---|---|---|---|
| Kneser e: T(½) = 1.6463542337511945809…, T(−½), T(⅓), slog(½), slog(2), slog(e/2) | 160 | 10^5.6 | 10^22.9 | 10^8.9 | 10^46.1 | 10^19.7 |
| Taylor-Koeffizienten a₁…a₁₅ von sexp_e bei 0 (a₁ = T′(0) = 1.0917673512583209918…), T″(0), slog′(0), slog″(0) | 158–160 | 10^5.5 | 10^22.4 | 10^8.9 | 10^45.6 | 10^19.2 |
| ∫₀¹T_e = 1.7145219588906204882…, ∫₋₁⁰T_e = 0.4992181828299913397… | 160 | 10^5.5 | 10^23.0 | 10^9.0 | 10^46.5 | 10^19.6 |
| D = −1.0552140087126523890…, 0!₂ = F′(1) = 0.5730489722129730654… | 160 | 10^5.6 | 10^23.0 | 10^8.9 | 10^46.2 | 10^19.5 |
| Basis 2: T₂(±½), ∫₀¹T₂, ∫₋₁⁰T₂, a₁…a₇ | 160 | 10^4.4 | 10^22.9 | 10^8.3 | 10^46.2 | 10^19.6 |
| prodtet_e(½) = 1.5134281085001618746…, α_{x·eˣ}(1) = −1.1259817765744955784…, α_{eˣ−1}(1) = −2.0259120986858638825…, Halbiterierte (eˣ−1)^{∘½}(1) = 1.2710274138899515214… | 175 | 10^6.2 | 10^25.4 | 10^9.9 | 10^51.3 | 10^21.9 |
| θ: Re/Im a₀ (mod 2πi/L), Re/Im/\|·\|/arg von a₁, a₂, a₃ | 160 | 10^5.5 | 10^22.9 | 10^8.9 | 10^46.0 | 10^19.6 |
| Gemeinsame Relation unter 18 Kneser-e-Zahlen + 1, π, Re L, Im L | 158 | 10^5.8 | – | – | – | – |

**Was die Schranken bedeuten:** Bei den kleinen Basen gibt es keine Relation mit Koeffizienten unter ≈ 10²² (Norm). Insbesondere ist **∫₀¹T_e** („Ist ∫₀¹ T_e eine bekannte Konstante?“) **keine** rationale Kombination von 1, π, e, γ, ln 2 mit Nennern bis etwa 10²². Sie ist auch nicht algebraisch vom Grad ≤ 6 mit Höhe ≲ 10¹⁹.

**Wörtliche Beinahe-Treffer** (x ≈ (p/q)·b mit \|p\| ≤ 36, q ≤ 12; relativer Abstand < 2·10⁻⁴). Bei 160 Stellen hält keiner davon:

| Zielgrösse | Beinahe-Wert | relativer Abstand |
|---|---|---|
| slog′(0) | Catalan G | 2.1·10⁻⁵ |
| T_e(½) | (19/8)·ln 2 | 7.9·10⁻⁵ |
| T₂(½) | (12/11)·Im L | 1.4·10⁻⁵ |
| ∫₀¹T₂ | 8·ln ζ(3) | 1.5·10⁻⁵ |
| 0!₂ | (3/7)·Im L | 9·10⁻⁵ |
| arg θ-a₁ | −(8/7)·Im L | 1.6·10⁻⁵ |
| ∫₀¹T_e | 12/7 | 1.4·10⁻⁴ |

Der alte Kandidat 0!₂ ≈ γ liegt um 4.2·10⁻³ daneben. Die Häufung T(−½) = 0.49856, slog″(0) = 0.49871, ∫₋₁⁰T = 0.49922 bei ½ kommt daher, dass T_e(x) auf [−1, 0] fast gleich x + 1 ist; eine Identität steckt nicht dahinter. `mpmath.identify` liefert bei 15 Stellen für jede Zahl eine Fünfterm-Formel. Alle diese Formeln weichen bei voller Genauigkeit um 10⁻¹³ bis 10⁻¹⁵ ab und sind damit reine Anpassungsartefakte (Liste in `relationen.json`, Schlüssel `identify`).

## Einordnung: „echte Identität“ im Sinn der Aufgabe?

- **Echt, aber klassisch:** K4 sowie K3 im Bell-Fall.
- **Echt und vermutlich neu formuliert:** K1 und K2. Sie verbinden die Tetrations-Störungsreihe mit Tangens- und Sekansfunktionen und den Tangentenzahlen. Dabei geht es nur um die Ränder der Polynome, nicht um den Kneser-Wert selbst. K5 als Zählfolge.
- **Keine Identität:** Für keinen Kneser- oder θ-Wert gibt es eine Verbindung zu klassischen Konstanten (Tabelle 2).

## Vorschläge für Etappe 2

1. **Singularitätsanalyse von θ:** Den algebraischen Abfall von a_n mit 256–1024 Koeffizienten bei dps 120 genau vermessen (die Daten im Motor-Speicher reichen fast). Ist es n^{−1}(ln n)^{−p} oder n^{−1−c}? Dann den Exponenten gegen Kandidaten aus L prüfen, etwa Re(1/L) = 0.169 oder 1 + Re(2πi/L)/2π. Ein exakt bestimmbarer Exponent wäre die erste echte Konstante, die Kneser kennzeichnet.
2. **q-Reihen-Weg:** Reguläre Tetration an L als Reihe in q = L^z schreiben (Reshetnikovs q-Binomial-Reihe, Aussage in 04_loesungsklassen) und θ = Σ a_n e^{2πinz} als Korrektur dazu. Dann prüfen, ob die a_n die Gestalt Σ_k c_k (L^k)ⁿ / (q;q)-Pochhammer-Summen haben. Dazu die erzeugende Funktion Σ a_n w^n auf [−1, 1] numerisch analytisch fortsetzen und ihre Singularitäten bei w = e^{2πi·(…)} suchen.
3. **η-Grenzfall → Digamma:** Die Störungsreihe K1/K2 ist genau die Entwicklung um den parabolischen Fall. Mit derselben Gewichtsmethode für b = e^{1/e}·(1+ε) die zweite Schicht ausrechnen. Sie sollte die Digamma-Grenzwerte aus 10_bruchsummen (Lambert-Reihen → ψ) als ln-cos/sec²-Glieder erklären.
4. **Borel-Summe der Störungsreihe:** Die Borel-Transformierte von Σ g_n(y) uⁿ bei u = 1 mit den Singularitäten aus K1 vergleichen. Die tan-Pole bei u = (π/y)²·2(k+½)² sind Kandidaten für die Stokes-Daten, und zu klären ist, ob ihre Lage die Kneser-Lösung auswählt (Aussage „Welche Tetration summiert die Störungsreihe?“).
5. **OEIS-Einträge:** K3 (Newton-Dreieck der Störungsreihe, n!g_n(4)) und K5 (Iterations-Dreieck von x·eˣ, Julia-Koeffizienten) als neue Folgen vorbereiten, falls die Literaturprüfung (Chapoton–Vallette) nichts findet.
