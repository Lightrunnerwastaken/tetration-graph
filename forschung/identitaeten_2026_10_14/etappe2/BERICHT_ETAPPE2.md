# Identitätssuche, Etappe 2 (14.10.2026): Warum Knesers θ an der reellen Achse singulär ist

## Kurz gesagt

- **Wo:** θ ist genau an den **ganzen Zahlen** singulär (x = …, −1, 0, 1, …) und sonst überall auf der reellen Achse analytisch. Grund: sexp(x) trifft dort die Werte 0, 1, e, eᵉ, … Die reguläre Lösung F am Fixpunkt L nimmt den Wert 0 **nie** an. Deshalb muss z + θ(z) bei z → −1 ins Unendliche laufen.
- **Welche Art:** Es ist weder eine Potenz noch ein Logarithmus. Es ist eine **„Superlogarithmus“-Singularität**. Nahe z = 0 gilt

  θ(z) = 2 − z + h(log(A·z) + κ·z + O(z²)),

  mit h = F⁻¹ (die reguläre Abelfunktion an L, holomorph in der ganzen oberen Halbebene), A = sexp′(0) und κ = sexp″(0)/(2A) − A/2. θ wächst also unbeschränkt, aber unvorstellbar langsam: von y = 10⁻² bis y = 10⁻³⁰ ändert sich θ(iy) nur um 0.14.
- **Folge für die Koeffizienten (die gesuchte Identität):**

  **a_n = [qⁿ] h(c + log(1−q)) + O(n⁻²),  c = ln(A/2π) + iπ/2.**

  Mit dem Korrekturglied δ₁ = ½ + iκ/2π wird der Rest O(n⁻³). Gemessen bis n = 4000: Der Rest fällt von 2·10⁻⁵ (erste Ordnung) auf 5·10⁻⁹ (zweite Ordnung, relativ). Für Basis 2 gilt dasselbe mit ihrem L, A und κ.
- **Exponent:** Der wahre Exponent ist **genau 1**. Führend gilt a_n ≈ −h′(c − ln n)/n, also |a_n| ≈ C/(n · ln n · ln ln n · …), ein Produkt iterierter Logarithmen. Der lokale Exponent sinkt wie 1 + 1/ln n gegen 1: Er beträgt 1.254 (n=2), 1.184 (64), 1.137 (1024) und 1.121 (4096). Damit ist die Drift aus Etappe 1 vollständig erklärt. Die Kandidaten aus Etappe 1, n^{−1−c} und n^{−1}(ln n)^{−p}, sind beide falsch.
- **arg a_n strebt nicht gegen −π/2.** Das Argument wandert langsam: −1.53 (n=1), −1.44 (64), −1.37 (4096). Für astronomische n schwankt es weiter, ohne Grenzwert. Der Grund ist, dass die Singularität sich selbst ähnlich ist (siehe 3.3).
- **Basis 2:** Die Struktur ist dieselbe, und das Modell hält ebenso gut. Die Zahlen liegen bei beiden Basen auffallend nahe beieinander (|a_4096|: 4.603·10⁻⁶ für e und 4.579·10⁻⁶ für 2). Das zeigt: Der Abfall ist universell, und L geht nur über h und die Konstante c ein.
- **q-Reihe:** θ setzt sich über jeden Bogen des Einheitskreises ausser q = 1 analytisch fort. |q| = 1 ist also **keine** natürliche Grenze. Eine Struktur wie bei (q;q)-Pochhammer- oder Lambert-Reihen, die typischerweise eine natürliche Grenze bei |q| = 1 haben, ist damit unwahrscheinlich. Eine genauere Prüfung habe ich nicht gemacht.

## 1. Bezeichnungen

Basis b (hier e oder 2). L ist der Fixpunkt b^L = L mit Im L > 0, und λ = L·ln b ist der Multiplikator. Ψ ist die Poincaré-Funktion mit Ψ(λu) = b^{Ψ(u)}, Ψ(0) = L und Ψ′(0) = 1. Sie ist ganz. F(s) = Ψ(λ^s) ist die reguläre Superfunktion, und h = F⁻¹ ist die reguläre Abelfunktion. Numerisch ist h(w) = m + log σ(u_m)/log λ mit Hauptzweig-Logarithmen u_{k+1} = log_b u_k (das ist `kern.Fix.h`). Kneser: sexp(z) = F(z + θ(z)), also θ(z) = h(sexp z) − z (modulo 2πi/log λ), und θ(z) = Σ a_n qⁿ mit q = e^{2πiz}.

## 2. Herleitung

| Schritt | Aussage | Status |
|---|---|---|
| S1 | Ψ′ hat keine Nullstelle. Daher hat F′ keine Nullstelle (F′(s) = log λ·λ^s·Ψ′(λ^s)). | **bewiesen**: Aus λΨ′(λu) = ln b·b^{Ψ(u)}Ψ′(u) und b^{Ψ} ≠ 0 folgt: Ψ′(u₀) = 0 ⇒ Ψ′(u₀λ^{−k}) = 0 für alle k. Diese Nullstellen häufen sich bei 0, also wäre Ψ′ ≡ 0. Widerspruch zu Ψ′(0) = 1. |
| S2 | Ψ und F lassen den Wert 0 aus. | **bewiesen**: Ψ(λu) = b^{Ψ(u)} ≠ 0. |
| S3 | Die Fortsetzung von h ist nur über den Werten 0, 1, b, b^b, … (Bahn von 0) behindert. | Standard für Poincaré-Funktionen: Die singulären Werte sind die Bahn der singulären Werte von b^w, und 0 ist der einzige solche Wert. Hier zitiert, nicht neu bewiesen. Auf der oberen Halbebene ist der Hauptzweig von h holomorph. Grund: Hauptzweig-log bildet die obere Halbebene in den Streifen 0 < Im < π ab, und die Iteration konvergiert gegen L. |
| S4 | θ ist auf ℝ∖ℤ analytisch fortsetzbar. | **hergeleitet**: sexp ist auf (−2, ∞) streng monoton und trifft die Bahn von 0 nur bei x ∈ {−1, 0, 1, …}. An allen anderen x ist h∘sexp nach S1/S3 lokal holomorph. Lücke: Dass der Hauptzweig die richtige Fortsetzung ist, ist numerisch geprüft (stetiger Verlauf ohne Sprünge in allen Abtastungen), aber nicht allgemein bewiesen. |
| S5 | θ ist an jeder ganzen Zahl singulär (unbeschränkt). | **bewiesen** (mit S2): sexp(z) → 0 für z → −1. Wäre z + θ(z) dort beschränkt, so gäbe es einen Häufungspunkt s* mit F(s*) = 0. Das widerspricht S2. Wegen der Periode 1 gilt dasselbe an allen ganzen Zahlen. |
| S6 | Lokale Form: θ(z) = 2 − z + h(ℓ₂(z)), ℓ₂ = log_b log_b sexp z. Für Basis e gilt ℓ₂ = log(Az) + κz + O(z²). | Identität trivial (h(w) = 1 + h(log_b w)); Taylor-Schritt elementar. **Numerisch** bestätigt für z = iy, siehe 3.1. |
| S7 | In q mit u = 1−q: z = i(u + u²/2 + …)/2π. Daraus folgt ℓ₂·ln b = c + log u + δ₁u + O(u²) mit c = ln(A/(2π ln b)) + iπ/2 und δ₁ = ½ + iκ/2π. | elementar |
| S8 | Übertragung auf Koeffizienten: a_n = b_n + O(n⁻²·h′) mit b_n = [qⁿ] h((c + log(1−q))/ln b). Der nächste Schritt ist b2_n = b_n − δ₁[(n+1)b_{n+1} − 2n b_n + (n−1)b_{n−1}]. Asymptotisch gilt a_n ≈ (1/n) Σ_{k≥1} g_k (−1/ln b)^k h^{(k)}((c − ln n)/ln b), mit 1/Γ(s) = Σ g_k s^k. | **hergeleitet mit Lücke**: Die Δ-Analytizität bei q = 1 folgt aus S4. Die Transfersätze von Flajolet–Odlyzko sind für „langsam variierende“ Funktionen der Form h(c + log u) aber nicht wörtlich bewiesen, weil h keine endliche Log-Potenz-Entwicklung hat. Die Zahlen bestätigen es sehr deutlich (3.2). |
| S9 | Führend gilt a_n ≈ −h′(c − ln n)/n (Basis e). Mit h′(v) = h′(log v)/v folgt \|a_n\| ≍ 1/(n·ln n·ln ln n ⋯). Der Exponent ist genau 1, der lokale Exponent ≈ 1 + 1/ln n + …. | folgt aus S8 |

**Selbstähnlichkeit:** Für v → −∞ + iπ/2 ist h(v) = 2 + h(log log v), und log log v liegt knapp über der reellen Achse bei grossem Realteil. Dort ist h(x) = slog(x) + θ(slog x). Die Singularität von θ bei 0 enthält also θ selbst, ausgewertet nahe der reellen Achse an der Stelle slog(ln ln ln(1/|z|)). Die feinere Asymptotik von a_n (Argument, Vorfaktor) ist deshalb keine Konstante, sondern schwankt periodisch in slog(ln ln n). Wenn ln ln ln n nahe bei sexp(k) liegt, gibt es Ausschläge (3.3). Das ist hergeleitet, aber nur qualitativ geprüft.

## 3. Numerik

**Motor:** fatou.gp-Fork (`kern.py`, eigener Zwischenspeicher `etappe2/cache/`). Alle Abtastungen liegen bei |Im z| ≤ 0.2, nahe der Achse (y = 0.002 / 0.003, dazu Kreise mit Radius 0.1 für A und B). Damit ist der Bereich der Motorgrenze für |Im z| > ¼ und des neuen ctedge-Fixes nicht berührt. Gegenprobe gegen Etappe 1 (alter Motor, dps 120, y = 0.2): a₁…a₆₄ stimmen relativ auf **3.5·10⁻⁸⁰** überein (`daten/abgleich_koeff.txt`).

### 3.0 Koeffizienten bis n = 4096 (`p2_koeff.py`)

Zwei Läufe je Basis: (dps 80, y = 0.003, N = 16384) und (dps 60, y = 0.002, N = 16384). Grösster relativer Unterschied: 5·10⁻⁶⁰ (n ≤ 512), 6·10⁻⁵⁶ (≤ 2048), 1.4·10⁻⁴⁴ (≤ 4096, begrenzt durch dps 60·e^{2πny}); bei Basis 2 2.3·10⁻⁴⁴. Die Koeffizienten sind also auf ≥ 44 Stellen sicher. Laufzeit je Lauf 1.5–3 min.

### 3.1 Lokale Form (`p5_lokal.py`, Basis e, dps 50)

Modelle: M1 = 2 + h(log(iAy)) und M2 = 2 − iy + h(log(iAy) + κiy).

| y | \|θ−M1\| | \|θ−M2\| |
|---|---|---|
| 10⁻² | 1.0·10⁻² | 1.7·10⁻⁶ |
| 10⁻³ | 1.0·10⁻³ | 9.8·10⁻⁹ |
| 10⁻⁵ | 1.0·10⁻⁵ | 4.9·10⁻¹³ |
| 10⁻⁸ | 1.0·10⁻⁸ | 2.6·10⁻¹⁹ |
| 10⁻¹² | 1.0·10⁻¹² | 1.5·10⁻²⁷ |

Der Fehler fällt wie y (M1) bzw. y² (M2), genau wie von S6 vorhergesagt. Ab y ≈ 10⁻²⁰ begrenzt die Genauigkeit von w − 1 bei dps 50 das Ergebnis. θ(iy) wandert dabei nur von 26.205 + 5.227i (y = 10⁻²) nach 26.330 + 4.763i (y = 10⁻³⁰).

### 3.2 Koeffizienten-Modell (`p3_modell.py`, `p4_korrektur.py`)

b_n kommt aus dem Cauchy-Integral von h(c + log(1−q)) auf |q| = e^{−2π·0.003} mit M = 16384 Punkten. Dafür braucht es keinen Kneser-Motor, nur h und A. T8 ist die Transfer-Entwicklung mit 8 Termen.

**Basis e** (A = 1.0917673512583209918, B = sexp″(0)/2 = 0.2714832129016945953, beide gleich wie in Etappe 1; κ = −0.29721969736603008485):

| n | \|a_n\| | lok. Exp. | arg a_n | \|a−b\|/\|a\| | \|a−b2\|/\|a\| | n²·\|a−b2\|/\|a\| | T8/a_n |
|---|---|---|---|---|---|---|---|
| 2 | 3.8014·10⁻² | 1.2538 | −1.5155 | 2.1·10⁻¹ | – | – | 0.9945 − 0.032i |
| 16 | 2.8899·10⁻³ | 1.2163 | −1.4720 | 9.3·10⁻³ | 7.3·10⁻⁴ | 0.186 | 1.00005 − 0.00087i |
| 64 | 5.4756·10⁻⁴ | 1.1845 | −1.4426 | 1.9·10⁻³ | 3.4·10⁻⁵ | 0.140 | 1.00001 − 0.00018i |
| 256 | 1.0802·10⁻⁴ | 1.1582 | −1.4148 | 4.0·10⁻⁴ | 1.7·10⁻⁶ | 0.111 | 1.000004 − 3.7·10⁻⁵i |
| 1024 | 2.2013·10⁻⁵ | 1.1375 | −1.3891 | 8.4·10⁻⁵ | 8.6·10⁻⁸ | 0.091 | 1.0000014 − 8.1·10⁻⁶i |
| 4096 | 4.6026·10⁻⁶ | 1.1210 | −1.3656 | 1.8·10⁻⁵ | – | – | 1.0000006 − 1.8·10⁻⁶i |
| 4000 | | | | 1.9·10⁻⁵ | 4.8·10⁻⁹ | 0.076 | |

**Basis 2** (A = 0.8893649546209763728, B = 0.0086765489653699340, κ = −0.43492658472636039403): Bei n = 4096 ist der lokale Exponent 1.1205 und arg a_n = −1.3843. Die relativen Reste sind \|a−b\|/\|a\| = 1.8·10⁻⁵ und \|a−b2\|/\|a\| = 4.8·10⁻⁹ (n = 4000), dazu n²·Rest 0.0765 und T8/a = 1.0000005. Die vollständigen Tabellen stehen in `daten/log_modell_2.txt` und `daten/korrektur_2.json`.

**Deutung:** n·(a−b)/a und n²·(a−b2)/a fallen langsam wie 1/ln n, genau so, wie S8 es für den nächsten Term (∝ h″/h′) vorhersagt. Der Rest schwankt nicht. Das spricht gegen weitere Singularitäten auf |q| = 1 mit vergleichbarer Stärke und passt zu S4.

### 3.3 Astronomische n (`p6_arg.py`, führender Term)

| ln ln n | n ln n \|a_n\| | arg a_n |
|---|---|---|
| 2.12 (n ≈ 4000) | 0.173 | −1.381 |
| 5 | 0.093 | −0.959 |
| 15.15 (= ln eᵉᵉ) | 0.040 | −0.209 |
| 30 | 0.0073 | +0.322 |
| 1000 | 5.7·10⁻⁵ | −0.297 |
| 3.8·10⁶ (ln ln ln n = eᵉ) | 3.1·10⁻³ | +0.294 |

n ln n |a_n| → 0 (Faktor 1/ln ln n und weitere). Das Argument hat keinen Grenzwert. Ausschläge kommen, wenn ln ln ln n nahe einem Turmwert eᵉ, … liegt; das ist die Selbstähnlichkeit aus Abschnitt 2. Die Werte sind eine Auswertung der Formel −h′(c − ln n) und nicht gemessen. Für so grosse n gibt es keine direkten Koeffizienten.

## 4. Was ist neu, was bleibt offen

- **Neu (nach meinem Wissen):** die exakte Beschreibung der Singularität von Knesers θ an den ganzen Zahlen (S5/S6), die geschlossene Asymptotik a_n = [qⁿ] h(c + log(1−q)) + O(n⁻²) mit c = ln(A/(2π ln b)) + iπ/2 und die Erklärung des Exponenten (genau 1, Drift 1/ln n). Literatur geprüft habe ich nicht. Paulsen–Cowgill (2017) behandeln θ, aber nach meiner Erinnerung nicht diese Asymptotik.
- **Lücken:** S3 zitiert, S4 Zweigwahl nur numerisch, S8 Transfer nicht streng.
- **Keine neue Konstante:** Der Abfall wird durch A = sexp′(0) und die Funktion h bestimmt. Eine einzelne „Kneser-Exponentenzahl“ gibt es nicht. Die Frage kn.frage.theta_exponent ist damit beantwortet: Der Exponent ist 1, mit iterierten log-Korrekturen.
- **Offen:** (i) S8 streng beweisen; (ii) die Schwankung des Arguments als Funktion von slog(ln ln n) explizit durch θ′ auf der reellen Achse ausdrücken; (iii) q-Reihen-Frage (Ziel 4) nur grob behandelt.

## Dateien

`kern.py` (Motor, L, Schröder, h), `p1_sonde.py` (θ(iy)), `p2_koeff.py` (a_n), `p3_modell.py` (b_n, Transfer), `p4_korrektur.py` (zweite Ordnung), `p5_lokal.py` (lokale Form), `p6_arg.py` (grosse n). Daten liegen in `daten/` dieses Ordners, Motor-Zwischenspeicher in `cache/` (~14 MB). Aufrufe:

```
python p2_koeff.py e 80 0.003 16384 4096 ; python p2_koeff.py e 60 0.002 16384 4096   (je ~2 min; dasselbe mit 2)
python p3_modell.py e 80 0.003 16384 daten/a_e_80_0.003_16384.json                       (~1 min)
python p4_korrektur.py e daten/a_e_80_0.003_16384.json daten/modell_e_80_0.003_16384.json
python p5_lokal.py ; python p6_arg.py
```

Interpreter: `C:/Users/janis/miniconda3/python.exe`, mit `PYTHONIOENCODING=utf-8`.
