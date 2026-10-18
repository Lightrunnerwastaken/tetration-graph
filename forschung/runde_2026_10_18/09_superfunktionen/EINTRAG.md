# Eintrag Bericht 09 (Superfunktionen-Zoo) in den Graphen

Dateien: `daten/aussagen/r18_09_kn.yaml`, `r18_09_it.yaml`, `r18_09_br.yaml`; Theorie `daten/theorie/r18_09_{kn,it,br}.yaml`.
Prüfung `baue.py --nur-pruefen`: 0 Fehler.

## Angelegte ids
- kn.satz.zoo_singulaerer_ort | θ singulär genau über den singulären Werten von f (bedingt: Zweigwahl)
- kn.satz.zoo_puiseux | aₙ ~ C₁ n^(−1−1/d) bei kritischem Wert vom Grad d (bedingt: (V) + Zweigwahl)
- kn.beob.zoo_mittelwerte | Knesers Mittelwertgesetze universell (x²+1, x⁴+1, cosh; numerisch)
- kn.beob.zoo_ableitung_nullstellen | S′ ≠ 0 nur bei eˣ; Nullstelle bei x²+1; universell 1+θ′ ≠ 0 (numerisch)
- kn.verm.zoo_lokal_konvex | Re(1+θ′) > 0 bei d < ∞, falsch bei eˣ (Vermutung)
- kn.beob.zoo_wendepunkt | genau ein Wendepunkt für x²+1, cosh, x⁴+1 (rechnergestützt)
- it.satz.poincare_polynom_surjektiv | Poincaré-Funktion von x^d + c lässt keinen Wert aus (bewiesen, folklore)
- br.satz.boettcher_nicht_kneser | Böttcher-Superfunktion ≠ Kneser (bewiesen)
- br.beob.boettcher_kneser_abstand | Abstand 4.5·10⁻⁷ bei x²+1, Defekt E = 3.9·10⁻¹⁴; schwaches Breitengesetz (numerisch)

Nicht eingetragen: K1/K2 (trivial), Z0 (Löser, steht in den Quellenangaben). Fehlende Quelle: Flajolet–Odlyzko 1990 (Transfersatz) ist nicht in `daten/quellen.yaml` — bei Bedarf anlegen.

## Vorschlag Stränge
- Neuer Strang `kn.zoo` (Gebiet kn), Name „Kneser für andere Abbildungen“, verlauf: „Knesers Bauweise für x²+1, x⁴+1, cosh: θ singulär über den singulären Werten, Abfall n^(−1−1/d), Mittelwertgesetze universell, S′ ≠ 0 nicht.“ Mitglieder: die sechs kn.*zoo*-ids.
- it.satz.poincare_polynom_surjektiv → bestehender Strang `it.hyperbolisch`.
- br.satz.boettcher_nicht_kneser, br.beob.boettcher_kneser_abstand → neuer Strang `br.boettcher` „Zwei natürliche Superfunktionen“, verlauf: „Polynome haben neben Kneser eine Böttcher-Superfunktion von ∞ her; verschieden, aber nur um 10⁻⁷.“ (alternativ in den kn.zoo-Strang).

## Vorschlag Ideenkarte
- Unteridee `fortsetzung.andere` (andere Konstruktionen) für die Böttcher-Einträge; `verwandte.elementar` oder eine neue Unteridee `verwandte.superfunktionen` („Kneser-Superfunktionen anderer Abbildungen: was an eˣ besonders ist und was universell“) für die kn.*zoo*-Einträge.
