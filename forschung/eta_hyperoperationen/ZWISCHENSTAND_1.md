# η höherer Hyperoperationen — Zwischenstand 1

Stand: 12.10.2026. Das bisherige H₀-Projekt ist pausiert; diese Untersuchung ist unabhängig.

## Einheitliche Indizierung

Wir schreiben H₃(b,x)=bˣ, H₄=Tetration, H₅=Pentation, H₆=Hexation, H₇=Heptation und H₈=Octation. Für r≥3 gilt

Hᵣ₊₁(b,0)=1 und Hᵣ₊₁(b,x+1)=Hᵣ(b,Hᵣ₊₁(b,x)).

ηᵣ bezeichnet die positive Berührbasis **der Abbildung Hᵣ**: Hᵣ(ηᵣ,eᵣ)=eᵣ und ∂ₓHᵣ(ηᵣ,eᵣ)=1. Sie ist unter den üblichen Monotonie- und Fixpunkthypothesen die Grenze beschränkter ganzzahliger Iteration für **Hᵣ₊₁**. Die vom Autor erinnerte Reihe beginnt somit mit η₃=e^(1/e).

## Bisherige Ergebnisse

| Berührende Abbildung | Grenze für | ηᵣ | eᵣ | Herkunft |
|---|---|---:|---:|---|
| H₃, Exponentialfunktion | Tetration | 1.444667861010 | 2.718281828459 | exakt |
| H₄, Tetration | Pentation | 1.635324496715 | 3.088532271807 | vorhandene Projektdaten |
| H₅, Pentation | Hexation | 1.737355925152 | 3.301407278523 | vorhandene Projektdaten |
| H₆, Hexation | Heptation | 1.799121236129 | 3.442234647717 | vorhandene Projektdaten |
| H₇, Heptation | Octation | **1.83975866072333057691675** | **3.5434897855274943163** | neue Berechnung, zwei Kontrollläufe |

η₇: Sekante in der Basis, Newton für ∂ₓH₇−1, Ableitungen durch Cauchy-Abtastung. Läufe mit angeforderten 25/45 Dezimalstellen (GP tatsächlich 38/57), 48/64 Punkten, Radius 0.04/0.025 stimmen in η₇ bis etwa 25 Stellen überein. Das ist **numerische Evidenz**, kein Intervallbeweis. Laufzeiten 401/320 Sekunden.

H₇ bei Basis 2 ist durch reguläre Iteration am negativen reellen Fixpunkt der Hexation konstruiert. H₇(2,1.5)≈2.65161799652261. Eine beschränkte H₈ bei Basis 1.83 wurde am positiven anziehenden Fixpunkt von H₇ konstruiert; ihr ganzzahliger Grenzwert ist ungefähr 2.89902557268021.

Die unbeschränkte H₈ bei Basis 2 ist inzwischen experimentell durch einen eigenen, getrennten Zwei-Bogen-Abel-Aufbau konstruiert: H₈(2,1.5)≈2.63586067475121. Geschätzte Genauigkeit aus zwei Bögen rund 16 Stellen; zusätzliche Kontrollen stehen noch aus. Der erste Versuch scheiterte an zu engen, von der Hexation übernommenen Zeitgrenzen des Bogenplaners. Das war ein Rechenhindernis und kein Nichtexistenzresultat.

## Auswahl und bisherige Literatur

Die Zahlen gelten für die hier verwendete Auswahl: Kneser-Tetration, reguläre Pentation, ein Zwei-Bogen-Abel-Aufbau für Hexation, reguläre Heptation und ein entsprechender experimenteller Abel-Aufbau für unbeschränkte Octation. Für die höheren geraden Stufen wird kein allgemeiner Kneser-Existenz- oder Eindeutigkeitssatz behauptet.

Bereits [Sheldons eigener Bericht von 2017](https://www.tetrationforum.org/showthread.php?pid=9061&tid=1177) beschreibt eine Heptationskonstruktion und nennt für die nächste Berührbasis ungefähr 1.839. Der zusätzliche Wert wird hier genauer nachgerechnet; er ist keine erstmalige Entdeckung der nächsten Stufe.

Die Vermutungen ηᵣ→2 und eᵣ→4 sind weiterhin offen. Ohne zusätzliche Auswahlbedingungen sind sie nicht universell: Eine Verschiebung des Höhenparameters durch x↦x+ε sin(2πx) erhält Rekursion und alle ganzzahligen Werte, verändert aber die reellen Zwischenwerte.

## Als Nächstes

Octation genauer gegenprüfen, ihre eigene positive Berührbasis suchen, weitere reguläre Stufen versuchsweise aufbauen und die η-Folge mit Vorhersagen aus den früheren Gliedern vergleichen. Jede fehlgeschlagene Konstruktion wird mit Ursache und Reichweite protokolliert.
