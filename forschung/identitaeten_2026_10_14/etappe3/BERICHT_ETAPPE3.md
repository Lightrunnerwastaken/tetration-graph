# Identitätssuche, Etappe 3 (14.10.2026): Integralformeln für Knesers sexp

## Kurz gesagt

- **Die Ableitung von Knesers sexp wird nie null** (ausser am Schnitt (−∞, −2]). Der Beweis braucht nur das Argumentprinzip, die Turmgleichung und Knesers Grenzwerte L, L̄ oben und unten. Folge: z ↦ z + θ(z) ist in der oberen Halbebene lokal umkehrbar, und log sexp′ ist dort eine eindeutige Funktion.
- **Mittelwerte über eine Periode sind multiplikativ.** Mit θ(x) = h(sexp x) − x (h ist die reguläre Abelfunktion am Fixpunkt L) gilt auf der reellen Achse zum Beispiel
  - ∫₀¹ θ² dx = (∫₀¹ θ dx)²,
  - ∫₀¹ σ(sexp x)·e^{−Lx} dx · ∫₀¹ e^{Lx}/σ(sexp x) dx = 1 (σ ist die Schröder-Funktion an L),
  - **∫₁ᵉ slog′(w)² / h′(w) dw = 1**,
  - ∫₀¹ ln sexp′(x) dx = −∫₀¹ ln|h′(sexp x)| dx.
  
  Für jede andere reelle Tetration sind diese Formeln falsch. Die Gegenprobe T(x) = sexp(x + ε sin 2πx) weicht schon bei ε = 10⁻⁴ um etwa 10⁻⁵ ab. Die Formeln sind also **echte Kennzeichen von Kneser**. Geprüft auf **10⁻⁹⁸** (Basis e) und 10⁻⁵⁸ (Basis 2).
- **Senkrechte Integrale haben geschlossene Werte:** Auf jeder senkrechten Geraden gilt ∫_{−∞}^{∞} [e^{sexp(x₀+it)} − sexp(x₀+it)] dt = **2 Im L**. Allgemeiner ergibt jede ganze Funktion Φ den Wert (Φ(L) − Φ(L̄))/i. Dazu kommt eine Raabe-artige Formel: ∫ₓ₀^{x₀+1} sexp = Re L + ∫₀^∞ Im[sexp(x₀+1+it) − sexp(x₀+it)] dt. Für Basis 2 heisst die rechte Seite 2 Im L₂.
- **Eine Binet-artige Formel für Knesers Asymptotik-Konstante:** sexp(z) ≈ L + κ·e^{Lz} für Im z → ∞, mit
  κ = exp(∫₀¹ ln sexp′ − 3L/2 + i∫₀^∞(sexp(it) − L) dt) = ∫₀¹ σ(sexp x) e^{−Lx} dx = 1.7176401959431580370… − 2.3844430117304858788… i.
- **Nebenergebnis zur offenen Vermutung „genau ein Wendepunkt“ (gr.verm.ein_wendepunkt), Basis e:** Das Argumentprinzip zählt die Nullstellen von sexp″ in Streifen. Im ganzen Gebiet [−1.8, 1.4] × [−4, 4] gibt es genau eine, den reellen Wendepunkt −0.5218. Zusammen mit zwei elementaren Fortpflanzungsargumenten folgt: **Kneser-sexp zur Basis e hat genau einen Wendepunkt.** Das ist rechnergestützt, aber nicht mit Intervallarithmetik zertifiziert.
- **Was nicht geklappt hat:** Abel–Plana für Σ f(sexp n) scheitert. sexp strebt abseits der Achse gegen L und L̄ statt gegen 0, deshalb fehlt die nötige Abklingbedingung. Eine Mellin-Transformierte mit eigener Funktionalgleichung habe ich nicht gefunden. Die Cauchy-Darstellungen von sexp′(0) und Ähnlichem sind richtig, aber nur die Cauchy-Formel selbst, keine neue Identität.

**Einordnung:** Die Formeln sind exakt. Das Muster dahinter ist das von Raabe und Binet bei Γ: Funktionalgleichung plus Verhalten im Unendlichen ergibt geschlossene Integralwerte. Hier kommt alles aus zwei Kneser-Eigenschaften. Erstens: θ ist in der oberen Halbebene holomorph und periodisch mit Grenzwert a₀. Zweitens: sexp → L, L̄. Keine der Formeln gibt eine Kneser-Zahl in klassischen Konstanten. Sie verbinden Kneser-Grössen untereinander und mit L.

## 1. Bezeichnungen (Basis e; für allgemeines b siehe `s5_basis_b.py`)

- L: Fixpunkt e^L = L mit Im L > 0. λ = L ist der Multiplikator, log λ = L (Hauptzweig).
- σ: Schröder-Funktion an L mit σ(e^w) = Lσ(w) und σ′(L) = 1. Fortgesetzt wird über Hauptzweig-Logarithmen.
- Ψ = σ⁻¹ ist ganz. F(s) = Ψ(e^{Ls}) ist die reguläre Superfunktion, h = F⁻¹, also e^{L h(w)} = σ(w).
- sexp(z) = F(z + θ(z)) mit θ(z) = Σ aₙ e^{2πinz}. θ ist holomorph und 1-periodisch in ℍ und strebt gegen a₀ für Im z → ∞ (Knesers Konstruktion; Paulsen–Cowgill).
- Auf (0, 1) bezeichnet θ(x) den Randwert von oben, θ(x) = h(sexp x) − x, auf dem stetigen Zweig. Er ist an den ganzen Zahlen singulär (Etappe 2).
- M[G] := ∫₀¹ G(x) dx ist der Periodenmittelwert.

## 2. Ergebnisse

| # | Aussage | Status | Genauigkeit | trivial / Kneser-spezifisch | Neuheit (Vermutung) |
|---|---|---|---|---|---|
| A | sexp′(z) ≠ 0 für alle z ∈ ℂ∖(−∞,−2]; gleichwertig 1+θ′ ≠ 0 in ℍ, und log sexp′(z+1) = sexp(z) + log sexp′(z) gilt in ganz ℍ mit dem stetigen Zweig | **bewiesen** (Beweis B-A) | Zählung numerisch: Windungszahl 0 auf 10⁻³⁰ in 12 Streifen `s4` | Kneser-spezifisch (braucht sexp → L oben) | möglicherweise bekannt, weil es aus der Univalenz in Knesers Konstruktion folgt; als eigenständiger Beweis per Argumentprinzip nicht gesehen |
| Z0 | **Mittelwertsatz:** Für jede in ℍ holomorphe 1-periodische G mit G → G_∞ oben und integrierbarer Majorante nahe ℝ ist M[G(x+i0)] = G_∞. Also ist M auf solchen Funktionen multiplikativ: M[G₁G₂] = M[G₁]M[G₂]. | **bewiesen**, mit Lücke bei der Majorante (B-Z) | – | Kneser-spezifisch (θ ∈ ℍ) | Prinzip klassisch (Hardy-Raum); Anwendung auf θ neu formuliert |
| Z1 | M[θ²] = M[θ]², M[θ³] = M[θ]³, M[e^{cθ}] = e^{c M[θ]} für jedes c ∈ ℂ; M[θ] = a₀ | hergeleitet (Z0) | 10⁻⁹⁸ (θ²), 10⁻⁹⁶ (θ³), 10⁻⁸⁸ (e^θ, Werte ~10¹¹), 10⁻¹⁰² (e^{iθ}) | Kneser-spezifisch; Gegenprobe ε = 10⁻⁴: Abweichung 9·10⁻⁶ | vermutlich neu |
| Z2 | **∫₀¹ σ(sexp x) e^{−Lx} dx · ∫₀¹ e^{Lx}/σ(sexp x) dx = 1**; allgemeiner M[(σ(sexp x)e^{−Lx})^k] = κ^k | hergeleitet (Z0 mit G = e^{kLθ}, ohne Zweigwahl) | 3·10⁻¹⁰⁰ | Kneser-spezifisch (ε = 10⁻⁴: 1.7·10⁻⁵) | vermutlich neu |
| Z3 | **∫₁ᵉ slog′(w)²/h′(w) dw = 1** (= M[1/(1+θ′)]), allgemein M[(sexp′·h′∘sexp)^s] = 1 für Re s < 1 | hergeleitet (Z0 + A) | 3·10⁻⁹⁹ (s = −1), 6·10⁻⁹⁹ (s = −2+i), 8·10⁻⁵³ (s = ½, Quadratur langsamer) | Kneser-spezifisch (ε = 10⁻⁴: 3.6·10⁻⁴) | vermutlich neu |
| Z4 | M[log(sexp′·h′∘sexp)] = 0, also **∫₀¹ ln sexp′ = −∫₀¹ ln\|h′(sexp x)\| dx** und ∫₀¹ arg h′(sexp x) dx = 0 (stetiger Zweig); ∫₀¹ ln sexp′ = 0.4984395595245179939040904770938110… | hergeleitet (Z0 + A) | 2.4·10⁻⁹⁹ | Kneser-spezifisch (ε = 10⁻⁴: 1.8·10⁻⁴) | vermutlich neu |
| Z5 | Randwerte von θ und 1/(1+θ′) haben keine negativen Fourier-Moden: ∫₀¹ θ(x) e^{2πimx} dx = 0 (m ≥ 1). Gleichwertig: Im θ und Re θ sind ein periodisches Hilbert-Paar (Kramers–Kronig für h∘sexp). Ausserdem aₘ = ∫₀¹ θ(x) e^{−2πimx} dx aus rein reellen Werten. | hergeleitet (Z0) | Moden m = 1, 2, 3: ≤ 1.2·10⁻⁹⁹; a₁…a₃ gegen Etappe 1 (alter Motor, y = 0.2, 160 Stellen): **5·10⁻¹⁰⁰** | Kneser-Kennzeichen (Umformulierung von „θ holomorph in ℍ“) | Kennzeichnung bekannt; konkrete Prüfung neu |
| V1 | **∫_ℝ [Φ(e^{s}) − Φ(s)] dt = (Φ(L) − Φ(L̄))/i** mit s = sexp(x₀+it), für jedes x₀ > −2 und jede ganze Funktion Φ; Φ(w) = w ergibt 2 Im L = 2.67447140286… | **bewiesen** (Cauchy am Halbstreifen, B-V) | 10⁻³⁹ (dps 40; Φ = w, w², 1/w; x₀ = 0, −0.4, 0.2); Basis 2: 2 Im L₂ = 3.13486…, 2·10⁻⁵⁹ | Kneser-spezifisch (Grenzwerte L, L̄) | vermutlich neu als Formel, Beweis eine Zeile |
| V2 | Raabe-Form: ∫ₓ₀^{x₀+1} sexp = Re L + ∫₀^∞ Im[sexp(x₀+1+it) − sexp(x₀+it)] dt; für x₀ = 0 ergibt das 1.71452195889062048816597476… (wie gr.beob.integrale_kneser) | **bewiesen** (B-V) | 4·10⁻⁴⁰ | Kneser-spezifisch | vermutlich neu |
| V3 | **Binet-Form der Asymptotik-Konstante:** sexp(z) = L + κe^{Lz} + o(e^{−y Im L}) mit κ = e^{La₀} = M[σ(sexp x)e^{−Lx}] = exp(∫₀¹ ln sexp′ − 3L/2 + i∫₀^∞ (sexp(it)−L) dt). Allgemein b: κ_b = λ^{a₀} = exp(∫₀¹ ln sexp′ + i ln b·K₀)/(√λ·log λ) mit K₀ = ∫₀^∞(sexp(it)−L)dt | **bewiesen** bis auf Z0-Lücke (B-V3) | e: Ganzzahligkeit der Zweigzahl 6 auf 10⁻⁴¹, κ auf 10⁻³⁵ (Druckgenauigkeit); Basis 2: Zweigzahl −7 auf 5·10⁻⁶⁰ | Kneser-spezifisch | vermutlich neu |
| V4 | Re ∫₀^∞ (sexp(x₀+it) − L) dt = Im L·x₀ + const | bewiesen (Ableitung nach x₀) | 2·10⁻⁴¹ | gilt für jede in ℍ holomorphe, auf ℝ reelle Funktion mit Grenzwert L | trivial |
| W | **Basis e: sexp″ hat in [−1.8, 1.4] × [−4, 4] genau eine Nullstelle, den reellen Wendepunkt.** Daraus folgt: genau ein Wendepunkt auf (−2, ∞). | rechnergestützt (Windungszahlen ganzzahlig auf 10⁻³⁰, Schritt ≤ 0.26 rad); **nicht** mit Intervallen zertifiziert | Zählung: N = 1 für x₀ ∈ {−1.4, …, −0.6}, sonst 0 | Kneser-spezifisch | löst gr.verm.ein_wendepunkt für b = e (bis auf Zertifizierung) |
| N | Abel–Plana für Σ f(sexp n), Mellin-Funktionalgleichung, Cauchy-Formeln für sexp^{(k)}(0) | negativ bzw. nur Umformulierung | – | – | – |

Konstanten (Basis e, 40 Stellen; Dateien in `daten/`):
- κ = 1.717640195943158037022581479756274116225 − 2.384443011730485878838953547660210549918 i, |κ| = 2.938682752341599794679794273585155142936, arg κ = −0.9465409639478231197117247491705600960889
- K₀ = ∫₀^∞(sexp(it) − L) dt = 1.0593125881982109936400184656200058 − 1.056719135810549650082724099018668 i
- M[θ] = a₀ = 26.193289231958250359366265551049680594540 + 5.425333084015578166577642116636883044617 i (der Zweig von Etappe 2; a₀ ist nur modulo 2πi/L bestimmt, κ nicht)
- κ₂ (Basis 2) = −0.0912950128173059192138757198914 − 3.00920981975562728783175111341 i

## 3. Beweise

**B-A (sexp′ ≠ 0).** Sei φ = sexp″/sexp′. Das ist meromorph auf ℂ∖(−∞,−2]. Aus sexp′(z+1) = e^{sexp z}·sexp′(z) folgt exakt

  φ(z+1) = φ(z) + sexp′(z).

Wähle x₀ > −2 so, dass sexp′ auf der Geraden x₀ + iℝ keine Nullstelle hat. Mit z also auch z+1, weil e^{sexp} ≠ 0; dafür reicht, abzählbar viele x₀ zu meiden. Sei R = [x₀, x₀+1] × [0, Y], und Y sei so gewählt, dass auf der Oberkante keine Nullstelle liegt. Die Unterkante liegt auf ℝ. Dort ist sexp′ > 0, und ihr Beitrag ist ln sexp′(x₀+1) − ln sexp′(x₀) = sexp(x₀). Die beiden Seiten zusammen ergeben i∫₀^Y sexp′(x₀+it) dt = sexp(x₀+iY) − sexp(x₀). Also gilt

  2πi·N(Y) = ∮_R φ = sexp(x₀+iY) − ∫_{x₀}^{x₀+1} φ(x+iY) dx.

Wegen sexp = F(z+θ) gilt φ = θ″/(1+θ′) + (1+θ′)·L·(1 + uΨ″(u)/Ψ′(u)) mit u = e^{L(z+θ)}. Für Y → ∞ gehen θ′ und θ″ gegen 0 (Fourier-Reihe) und u ebenfalls, denn |u| ≈ e^{−Y Im L}. Also strebt φ gleichmässig gegen L, und die rechte Seite strebt gegen L − L = 0. N(Y) ist ganzzahlig, nicht negativ und wachsend in Y, also ist N ≡ 0. Hätte sexp′ eine Nullstelle z₀ ∈ ℍ, so auch z₀ + n für alle n ≥ 1. Eine davon läge in einem solchen Halbstreifen, Widerspruch. ℍ⁻ folgt durch Spiegelung, ℝ ist bekannt. ∎
Nebenprodukt: Für N = 0 ist die Gleichung exakt. Also ist log sexp′(x₀+1+iY) − log sexp′(x₀+iY) = sexp(x₀+iY), ohne 2πik.

**B-Z (Mittelwertsatz).** G sei holomorph und 1-periodisch in ℍ. Dann gilt d/dy ∫₀¹ G(x+iy) dx = i[G(1+iy) − G(iy)] = 0. Für y → ∞ ist der Mittelwert G_∞, für y → 0 der Randmittelwert, sofern eine integrierbare Majorante existiert (dominierte Konvergenz). Anwendung:
- G = Φ(θ) mit ganzem Φ (Z1, Z2). Hier θ(z) = 2 − z + h(log(Az) + O(z)) (Etappe 2, S6), und h wächst über jede iterierte Logarithmenstufe hinweg nur wie slog. Deshalb ist |θ| ≲ log log log(1/|z|), und e^{c|θ|} ist integrierbar.
- G = (1+θ′)^s = exp(s log(1+θ′)) und G = log(1+θ′) (Z3, Z4). Nach A ist 1+θ′ ≠ 0, also ist log(1+θ′) eindeutig, und es gilt G_∞ = 1 bzw. 0. Nahe 0 ist |1+θ′| ≍ |h′(log Az)|/|z| ≲ 1/(|z| log(1/|z|)). Damit ist |G| ≲ |z|^{−Re s}, integrierbar für Re s < 1.
- Fourier-Moden G·e^{2πimz}, m ≥ 1, mit G_∞ = 0 (Z5).

**Lücke:** Die Majoranten sind aus der lokalen Form abgelesen, gleichmässig in einem Sektor über 0. Eine ausgeschriebene Abschätzung von h und h′ auf dem Bild log(Az), 0 < arg z < π, fehlt.

Umrechnungen:
- Z3 entsteht durch w = sexp x, dx = slog′(w) dw und 1/(sexp′(x) h′(sexp x)) = slog′(w)/h′(w).
- Z2 benutzt e^{Lθ(x)} = e^{L h(sexp x)} e^{−Lx} = σ(sexp x) e^{−Lx}.
- Z4: Der Realteil von log(1+θ′) = ln sexp′ + log h′(sexp) ergibt die erste Gleichung, der Imaginärteil die zweite, weil sexp′ > 0.

**B-V (V1, V2).** Cauchy auf dem Halbstreifen [x₀, x₀+1] × [0, ∞) für Φ∘sexp: Unterkante + i∫₀^∞[Φ(sexp(x₀+1+it)) − Φ(sexp(x₀+it))] dt − Φ(L) = 0. Die Oberkante geht gegen Φ(L), und die Konvergenz ist exponentiell, ~e^{−t Im L}. Spiegelung ergibt die untere Hälfte mit Φ(L̄). Zusammen hebt sich die Unterkante weg (V1); mit Φ = id ist der Imaginär- bzw. Realteil V2. ∎

**B-V3.** E(z) = log sexp′(z) − L·h(sexp z) = log(1+θ′) + log L + log Ψ′(u). Der Zweig ist durch A eindeutig. Es gilt E(z+1) − E(z) = sexp(z) − L und E → log L = L oben. Das Mittelwertargument liefert für den Mittelwert m(y) = ∫₀¹ E(x+iy) dx die Ableitung m′(y) = i(sexp(iy) − L). Also ist ∫₀¹ E(x) dx = L − i∫₀^∞(sexp(it) − L) dt. Mit ∫₀¹ h(sexp x) dx = a₀ + ½ (Z1) folgt V3 modulo 2πi. Die Zweigzahl misst nur, welchen Repräsentanten von a₀ man nimmt; κ = e^{La₀} ist eindeutig. ∎

**W (Wendepunkt, Basis e).** Nach A sind die Nullstellen von sexp″ genau die von φ. Die Windungszahl von φ auf ∂([x₀, x₀+1] × [−4, 4]) ist gleich Δarg/π längs des oberen Halbwegs (Spiegelsymmetrie). Sie ist 0 für x₀ = −1.8, −1.6, −0.4, …, 0.4 und 1 für x₀ = −1.4, …, −0.6. Oberhalb Im = 4 gilt |φ − L| ≪ |L|, dort liegt keine Nullstelle. Daraus folgt:
- sexp″ < 0 auf [−1.6, −0.6], weil keine Nullstelle und sexp″(−1) < 0.
- Genau eine (reelle, einfache) Nullstelle in [−0.6, 0.4].
- sexp″ > 0 auf [−0.4, 0.6] und nach gr.kor.konvex_nach_oben auf ganz [−0.4, ∞).
- Auf (−2, −1.6) gilt sexp″(x) = [sexp″(x+1)·sexp(x+1) − sexp′(x+1)²]/sexp(x+1)² < 0, weil sexp″(x+1) < 0.

Also gibt es genau einen Wendepunkt (−0.5218…). **Lücke:** Die Windungszahlen sind in Gleitkomma (dps 30) gerechnet, nicht in Intervallarithmetik.

## 4. Numerik und Motor

- **Reelle Achse:** Die Taylor-Reihen von sexp um −0.4, −0.2, 0, 0.2, 0.4 kommen aus Cauchy-Kreisen mit r = 0.2 (|Im z| ≤ 0.2, sicherer Motorbereich), mit dps 130 und M = 256. Sie stimmen mit dps 100/M = 128 überein, und sexp′(0) = A wie in Etappe 1. Ausgewertet wird nur für |x − c| ≤ 0.1. Damit ist sexp(x) − 1 auch bei x ~ 10⁻¹⁰⁰ ohne Auslöschung bekannt.
- **Quadratur:** tanh-sinh auf (0, 1) mit Schrittweite 2^{−k}, k = 5, 6, 7, und Knoten bis zum Randabstand 10⁻¹⁰¹. θ wird stetig gemacht (modulo 2πi/L), log(1+θ′) modulo 2πi.
- **Obere Halbebene:** sexp = F(z + θ(z)) mit θ aus den Fourier-Koeffizienten von Etappe 2 (n ≤ 4096). Gegen den Motor bei Im z = 0.05…0.2 stimmt das auf 4·10⁻⁶⁰ (`t_lib.py`).
- **Motor-Hinweis:** Ein Test des aktuellen Forks mit Im z = 0.8 … 6 bei dps 60 lief länger als 10 min und wurde abgebrochen (eigener Prozess). Für grosse Im wurde deshalb nur die θ-Reihe benutzt. Die Grenze „|Im z| > ¼“ war damit nicht nötig.
- **Querprüfung alter Aufbau:** a₁, a₂, a₃ aus den reellen Achsenintegralen (dps 100) stimmen mit Etappe 1 (alter Motor, Kreis y = 0.2, 160 Stellen) auf 5·10⁻¹⁰⁰ überein.

## 5. Dateien und Aufrufe

Interpreter `C:/Users/janis/miniconda3/python.exe`, `PYTHONIOENCODING=utf-8`; jeder Lauf < 1 min ausser wo vermerkt.

| Skript | Inhalt | Ausgabe |
|---|---|---|
| `k3.py`, `lib3.py` | Motorzugang (Speicher `cache/`), L, σ, Ψ, F, h, h′, θ-Reihe | – |
| `s1_taylor.py 130 256`, `s1b_taylor_links.py 80 256` | Taylor-Daten (je ~20 s) | `daten/taylor_*.json` |
| `s2_nullmoden.py 100` | Z1–Z5 (~25 s) | `daten/nullmoden_100.json`, `daten/log_s2_100.txt` |
| `s2b_kontrolle.py 30` | Gegenprobe T = sexp(x + ε sin 2πx) | Konsole |
| `s3_senkrecht.py 40` | V1–V4 (~3 min) | Konsole |
| `s4_zaehlung.py` | Nullstellenzählung sexp′, sexp″ (~2 min) | `daten/zaehlung.json` |
| `s5_basis_b.py 2 50` / `e 40` | Z, V1, V3 für allgemeine Basis | `daten/taylor_b*.json` |

## 6. Offen / Vorschläge

1. Z0-Majorante ausschreiben, also h, h′ auf {log(Az) : 0 < arg z < π} abschätzen. Dann sind Z1–Z5 und V3 voll bewiesen.
2. W zertifizieren: Windungszahl mit Intervallarithmetik (python-flint arb) auf den Taylor- und θ-Daten. Dasselbe für Basis 2, wo der Wendepunkt bei −0.0302 liegt, und als Familie in b.
3. Literatur: Prüfen, ob sexp′ ≠ 0 (A) bei Kneser, Kouznetsov oder Paulsen steht.
4. Grapheinträge: A (Satz), Z0 bis Z5 (Satz/Korollare, „Kneser-Kennzeichen“), V1/V2 (Satz), V3 (Korollar mit κ), W (Beobachtung, rechnergestützt) zu gr.verm.ein_wendepunkt. Für V1 und Z3 gibt es einfache Fassungen: „die Fläche zwischen e^T und T auf einer senkrechten Geraden ist 2 Im L“ und „∫₁ᵉ slog′²/h′ = 1“.
