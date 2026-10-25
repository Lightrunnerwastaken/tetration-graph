# Gebiet 02 — Die Schwankung von arg a_n der θ-Welle, explizit (Runde 25.10.2026)

**In einfacher Sprache.** Knesers Tetration enthält eine kleine Welle θ mit Koeffizienten a_n. Etappe 2 hatte gezeigt, dass der Winkel
(das Argument) von a_n nie zur Ruhe kommt. Jetzt ist klar, *wovon* er abhängt: Der Winkel von a_n ist bis auf π und eine kleine Korrektur
π/(2 ln n) gleich dem Winkel von 1 + θ′ an der **reellen** Stelle s ≈ 1/(2πn). Er wird also von der Welle selbst an einem Punkt der reellen
Achse bestimmt. Über eine Periode betrachtet ist dieser Winkel im Mittel genau 0. An den Rändern der Periode wiederholt sich das ganze Bild in
iterierten Logarithmen, jeweils um π verschoben. Deshalb wandert arg a_n ohne Grenzwert. Alles ist mit den gemessenen a_n bis n = 4096 und mit
einer zweiten Rechnung (direkt über h) auf 10⁻⁶ bis 10⁻⁸ geprüft.

Bezeichnungen wie in Etappe 2 (`forschung/identitaeten_2026_10_14/etappe2/BERICHT_ETAPPE2.md`): h ist die reguläre Abelfunktion am Fixpunkt L
(Hauptzweig, holomorph in der oberen Halbebene), A = sexp′(0), c = ln(A/2π) + iπ/2, θ(z) = Σ a_n qⁿ mit q = e^{2πiz}. Führend gilt
a_n = −h′(c − ln n)/n · (1 + O(1/ln n)) (kn.satz.theta_koeffizienten). Neu eingeführt werden

- **Φ(s) := arg(1 + θ′(s + i0))** für reelles s ∉ ℤ, die Phase der Welle auf der reellen Achse (1-periodisch),
- **P(Y) := arg h′(Y + iπ)**, die Phase von h′ auf der Geraden Im w = π (dort ist h holomorph),
- **σ_n := 1 + slog(A/(2πn))**, das ist ≈ 1/(2πn) (σ_n → 0⁺).

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit-Vermutung |
|---|---|---|---|---|
| 1 | Φ(s) = arg h′(sexp s + i0) = arg h′(sexp(s−2) + i0) = P(ln(−sexp(s−2))) − π (mod 2π) für s ∈ (0,1). Eine Periode der Welle entspricht genau der ganzen Geraden Im w = π. | hergeleitet mit Lücke (nur die Zweigwahl S4 aus Etappe 2) | Teil B: Newton für s + θ(s) = h(x+i0) liefert reelles s (\|Im s\| ≤ 1.3·10⁻⁷ für x = 1.2…2.3). sexp′(s) = (1+θ′)/h′ ist reell (Im/Re ≤ 2·10⁻⁶), und die Argumente stimmen auf 2·10⁻⁶ überein. p3: P(Y) − π gegen Teil B stimmt auf 2·10⁻⁶ | neu |
| 2 | **arg a_n = π + Φ(σ_n) + arctan(π/(2(ln n − ln(A/2π)))) + O(1/(ln n·ln ln n))**. Für den Betrag gilt \|a_n\| ≈ \|1+θ′(σ_n)\|/(n·sexp′(σ_n−2)) ≈ \|1+θ′(σ_n)\|/(2πn²), kompakt a_n ≈ −(1+θ′(σ_n))/(2πn²)·e^{iπ/(2 ln n)} | hergeleitet mit Lücke (Fehlerordnung nicht streng) + numerisch | Abstand zum führenden Term F0: −5.2·10⁻² (n = 4096), −4.0·10⁻⁶ (ln ln n = 10), −3.5·10⁻⁸ (ln ln n = 15.154, also auch an einem Turmwert), −9.8·10⁻¹⁵ (30). (F1−F0)·ln n·ln ln n bleibt in [−3.1, −0.9]. Basis 2 verhält sich gleich (−4.9·10⁻²). Modulquotient Achse/führend: 0.895, 0.937, 0.948, 0.956 (n = 16…4096) → 1 | neu |
| 3 | Selbstähnliche Form: arg a_n = Φ(slog n) − π/ln ln n + O(1/(ln ln n·ln ln ln n)), gleichwertig a_n/h′(n+i0) → 1. Gilt **nur**, wenn slog n mod 1 von 0 wegbleibt; nicht gleichmässig | hergeleitet + numerisch | F3−F0: +5.1·10⁻⁴ (ln ln n = 1000), +2·10⁻⁵ (10⁴). An Turmstellen versagt die Form: −1.17 (ln ln ln n = e), −1.43 (= eᵉ). Bei n = 4096 ist sie unbrauchbar (\|a_n/h′(n)\| = 0.39). Bei Basis 2 sind die Turmwerte 2, 4, 16 selbst ganze Zahlen n. Dort ist h′(n+i0) singulär, und das Verhältnis \|a_n/h′(n)\| fällt auf 0. | neu |
| 4 | Φ ist auf [0.02, 0.98] streng fallend, von +1.0970 auf −0.9930. Es gilt Φ(½) = −0.0165, und das Minimum von \|1+θ′\| liegt bei ≈ 0.684 nahe s = ½. Basis 2 ist fast gleich: +1.0844 und −1.0020, Φ(½) = −0.0118, \|1+θ′(½)\| = 0.686. | numerisch | zwei Datensätze (dps 80/y = .003 und dps 60/y = .002) stimmen auf 10⁻⁸ überein, η → 0 quadratisch extrapoliert. Teil B für Basis 2: \|Im s\| ≤ 1.1·10⁻⁷ für x = 1.2, 1.5, 1.8, 2.3, 2.6. Bei x = 2, einem Turmwert von Basis 2, versagt das Newton-Verfahren erwartungsgemäss. | neu |
| 5 | Mittelwert über eine Periode: ∫₀¹ Φ = 0 und ∫₀¹ ln\|1+θ′\| = 0. Das ist das Mittelwertprinzip kn.satz.theta_mittelwerte, angewandt auf log(1+θ′). | bewiesen, unter der Voraussetzung, dass 1+θ′ in der oberen Halbebene keine Nullstelle hat; diese Voraussetzung ist numerisch geprüft | Umlaufzahl von 1+θ′ auf \|q\| = e^{−2πη}: 0.000000 für η = 0.003 und 0.002. Unterhalb davon gilt sexp′ > 0 auf ℝ, und Nullstellen von 1+θ′ sind genau Nullstellen von sexp′, weil h′ ≠ 0 ist | Folgerung aus Bekanntem |
| 6 | Randgesetz: Für s → 0⁺ gilt Φ(s) ≈ Φ(slog ln Y) **+ π** mit Y = ln(−sexp(s−2)), also eine Kopie, die um π verschoben ist. Für s → 1⁻ entsteht eine Kopie ohne Verschiebung. Über die aufgelösten Ebenen reicht der stetige Zweig von −1.17 bis +4.79 (Spannweite 5.96 > π). | numerisch (mod 2π hergeleitet) | p3, adaptive Spur mit 3025 Punkten: Φ(Y = e^{3.5}) = 3.448, Φ(e^{10}) = 2.526, Φ(e^{16}) = 4.251, Φ(e^{19}) = 3.994. Der Übergang bei ln Y = eᵉ ist eine stetige Drehung um ≈ +2.1 | neu |
| 7 | Der stetige Zweig von Φ ist bei s → 0⁺ nach oben unbeschränkt, weil jede tiefere Ebene +π hinzufügt. Damit wächst auch der stetige Zweig von arg a_n (in n) unbeschränkt, aber unvorstellbar langsam. | Vermutung | Selbstkonsistenz: Ist R der Wertebereich, so müsste R + π ⊆ R̄ gelten, was für ein beschränktes R unmöglich ist. Die Verschiebung +π hält am ersten aufgelösten Übergang. Die zweite Ebene liegt bei ln Y ≈ e^{eᵉ} und ist nicht gerechnet | neu |
| 8 | Eine reine Taylorentwicklung an der Achse auf der Ebene von n (Achse + iε mit ε ≈ π/ln ln n) konvergiert für n ≤ 4096 nicht, weil ε ≈ 1.5 grösser ist als der Abstand zum Turmpunkt 1. | negativ | ln ln 4096 = 2.12, ε ≈ 1.48, x = 0.75 | – |

## Herleitung

**(1) Achsenformel.** Auf der oberen Halbebene gilt h(sexp z) = z + θ(z), also 1 + θ′(z) = h′(sexp z)·sexp′(z). Für reelles s ∉ ℤ ist
sexp′(s) > 0, also Φ(s) = arg h′(sexp s + i0). Die Kettenregel h′(w) = h′(log w)/w ändert das Argument nicht, solange w > 0 ist. Zweimal
angewandt (sexp s > 0 und sexp(s−1) > 0 für s ∈ (0,1)) folgt Φ(s) = arg h′(sexp(s−2) + i0) mit sexp(s−2) = −X < 0. Ein dritter Schritt
log(−X + i0) = ln X + iπ liefert arg h′(−X+i0) = P(ln X) − π. Durchläuft s die Periode (0,1), so läuft ln X über ganz ℝ. Lücke: Dass der
Randwert des Hauptzweigs die richtige Fortsetzung ist (S4 aus Etappe 2), ist nur numerisch gesichert. Teil B bestätigt es hier erneut.

**(2) Formel bei σ_n.** Für n a_n ≈ −h′(v) mit v = Y_n + iπ/2 und Y_n = ln(A/2πn) < 0 ist h nahe der negativen Achse holomorph, denn log w
bleibt dort in der oberen Halbebene. Zerlegt man h′(Y+iπ/2) = h′(log(Y+iπ/2))/(Y+iπ/2), so unterscheidet sich der äussere Faktor von dem bei
Y + i0 genau um arctan(π/2\|Y\|) im Argument. Der innere Punkt verschiebt sich um −i·arctan(π/2\|Y\|). Das ergibt einen Beitrag
(π/2\|Y\|)·Re(h″/h′)(ln\|Y\| + iπ) = O(1/(\|Y\| ln\|Y\|)), weil h″/h′ ≈ −1/w bei grossem Re w. Schliesslich ist arg h′(Y_n + i0) = Φ(σ_n),
denn sexp(σ_n − 2) = ln(A/2πn) nach Definition von σ_n. Für den Betrag folgt dasselbe mit sexp′(σ_n−2) = sexp′(σ_n)/(sexp σ_n · sexp(σ_n−1)) ≈ 1/σ_n ≈ 2πn.

**(3) Ebene von n.** Zwei Logarithmen bringen v auf w₂ = ln ln ln n + iπ/ln ln n + …. Dann ist arg h′(w₂) ≈ arg h′(ln ln ln n + i0) + (π/ln ln n)·Re(h″/h′).
Dabei gilt arg h′(ln ln ln n + i0) = arg h′(n + i0) = Φ(slog n), weil n, ln n und ln ln n positiv sind. Die Korrektur −arg log v = −π/ln ln n
überwiegt. Die Form ist nicht gleichmässig: Liegt ln ln ln n nahe einem Turmwert eᵉ…, so ist h′ dort an der Achse singulär.

**(6/7) Ränder.** Für Y → +∞ gilt P(Y) = arg h′(ln Y + iπ/Y) − arctan(π/Y). Der Punkt nähert sich der Achse bei ln Y, also ist
Φ(s) ≡ Φ(slog ln Y) − π (mod 2π). Die stetige Spur in p3 legt die Verschiebung als +π fest. Für Y → −∞ landet man nach zwei Logarithmen
wieder nahe der Achse, ohne Verschiebung. Weil arg a_n ≈ π + Φ(σ_n) mit σ_n → 0⁺ gilt, ist die Langzeitschwankung von arg a_n genau
dieses Randgesetz von Φ bei 0⁺. Das erklärt die Werte aus Etappe 2, Abschnitt 3.3 (z. B. −0.208 bei ln ln n = 15.154 und +0.322 bei 30).

## Dateien (alle Zahlen reproduzierbar, `C:/Users/janis/miniconda3/python.exe`, `PYTHONIOENCODING=utf-8`)

- `kern.py`: Kopie des Etappe-2-Kerns (nur die Klasse `Fix` wird benutzt; kein Motoraufruf, kein Zwischenspeicher).
- `p1_phase.py [e|2]`: Teil A (Φ aus den a_n, Umlaufzahl, Mittel), Teil B (reelle Achse, Newton), Teil C. Ausgabe `daten/p1_e.json`, ~1 min.
- `p2_formeln.py [e|2]`: F0/F1/F3 gegen die Daten und für astronomische n. Ausgabe `daten/p2_*.json`, < 1 min.
- `p3_gerade.py`: Φ über eine ganze Periode als P(Y) − π, mit adaptiver stetiger Spur. Ausgabe `daten/p3_e.json`, `daten/p3_phi.npy`, ~1 min.
  Hinweis: In der Ausgabetabelle meint die Zeile „0.000“ die Stelle Y = 0, alle anderen Zeilen meinen \|Y\| = e^{Spalte}.
- `p4_modul.py [e|2]`: Modulquotient Achse/führend.
- Eingabedaten (nur gelesen): `forschung/identitaeten_2026_10_14/etappe2/daten/a_*_16384.json`, `korrektur_*.json`.

## Vorschläge für Grapheinträge

- „Phase der θ-Koeffizienten = Phase von 1+θ′ bei 1/(2πn)“: arg a_n = π + Φ(σ_n) + arctan(π/2(ln n − ln(A/2π))) + O(1/(ln n ln ln n)), mit a_n ≈ −(1+θ′(σ_n))/(2πn²) — Satz, hergeleitet mit Lücke + numerisch
- „Eine Periode von θ′ ist die Gerade Im w = π von h′“: Φ(s) = arg h′(Y+iπ) − π, Y = ln(−sexp(s−2)) — Satz, hergeleitet mit Lücke (S4)
- „Selbstähnliche Form arg a_n ≈ arg h′(n) − π/ln ln n, nicht gleichmässig an Turmstellen“ — Beobachtung/Satz, hergeleitet + numerisch
- „Mittlere Phase von 1+θ′ über eine Periode ist 0“ — Folgerung aus kn.satz.theta_mittelwerte, bewiesen unter der Voraussetzung, dass 1+θ′ in der oberen Halbebene keine Nullstelle hat (numerisch geprüft)
- „Randkopien von Φ mit Verschiebung +π; stetiger Zweig von arg a_n unbeschränkt?“ — Frage/Vermutung
- Antwort-Ergänzung zu kn.frage.theta_exponent: der offene Punkt (ii) „Schwankung von arg a_n explizit durch θ′“ ist mit Befund 2 beantwortet.

## Offene Fragen

1. Befund 7: Bleibt die Verschiebung +π auf allen tieferen Ebenen erhalten? Dann wäre der stetige Zweig unbeschränkt. Rechnen müsste man auf der
   zweiten Ebene, also bei ln Y ≈ e^{eᵉ} mit einer Auflösung von etwa e^{−eᵉ} in Im. Eine Möglichkeit ist, P direkt über eine verschachtelte
   Kettenregel zu verfolgen.
2. Fehlerordnung in Befund 2 streng machen. Sie hängt an einer Schranke für h″/h′ auf der Geraden Im = π und an Teil S8 aus Etappe 2.
3. Gibt es eine Vollentwicklung, also h′(Y+iπ/2) als Reihe in (iπ/2)^k h^{(k+1)}(Y+i0), übersetzt in Ableitungen von θ an σ_n? Die Reihe
   konvergiert für \|Y\| > π/2, also für alle n ≥ 2.
