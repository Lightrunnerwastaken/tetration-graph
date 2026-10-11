# Versuchsgeschichte und offene Prüfungen

Arbeitsstand vom 11.10.2026. Es wird kein vollständiger Beweis der ursprünglichen Grad-N-Vermutung beansprucht.

## 1. Hermite-Konstruktion

**Herleitung abgeschlossen, unabhängige P-Prüfung offen.** Statt nur N+1 Koeffizienten werden bis Grad 5N verwendet.
Eine Korrektur vom Grad höchstens 2N+1 erhält den N-Jet in 0 und stellt den N-Jet in α wieder her.
Die explizite Cauchy-Schranke hat geometrische Basis 0.675622435625 < 1.
Der Beweis und die endliche Konstruktion stehen in `BERICHT.md`, §1.

**Numerisch:** Rechnungen bei N = 5, 10, 20, 40 mit 160 Dezimalstellen und zertifizierten Koeffizientenmittelpunkten.
Die Ergebnisse stehen vollständig in `ergebnisse/lemma_checks.json`. Keine Zertifizierung des gesamten numerischen Fehlers;
der bekannte Taylorschwanz endet bei 640. Die ursprüngliche Folge p_N wird durch diese Konstruktion nicht identifiziert.

## 2. Adjungierte Gleichung

**Herleitung abgeschlossen, unabhängige P-Prüfung offen.** Die rationale Darstellung gilt für beliebige holomorphe Keime.
Bei nichtpolynomialem m muss der holomorphe Term B_y erhalten bleiben; er kann im Allgemeinen nicht durch ein Polynom ersetzt werden.
Für die lokale Darstellung mit inversen Ästen kommt zudem der Konturterm I_F vor.
Eine exakte rationale Gegenprüfung bei N = 5 bestätigt die Gradschranke, ersetzt aber nicht den allgemeinen Beweis.

## 3. K2-Randintegral: erfolglose und unvollständige Versuche

**Status: kein neues Zertifikat.** Aus alten Projektlogs wurden nichtendliche Integralwerte festgestellt.
Diese zählen nicht als Schranken. Der eigene Code weist nichtendliche Intervalle ausdrücklich zurück.

Die während des Durchgangs untersuchten Varianten waren:

1. Direkte skalierte Poincaré-Jets mit Cauchy-Rest: Einschließungen wurden zu breit.
2. Punkt-Jets plus Cauchy-Schranke für die Verschiebung innerhalb einer Zelle: schärfer, aber kein vollständiger Randlauf.
3. Winkelreihen durch formale Inversion: ein Kompositionsproblem bei der Nullreihe wurde durch Horner-Auswertung umgangen;
   die Variante führte dennoch zu keinem vollständigen Zertifikat und ist nicht aktiv geschaltet.
4. Letzter Hybridstand: außerhalb der Ecknähe Punkt-Jets bis Ordnung 26 mit Cauchy-Tail; nahe der Ecke
   ursprüngliche Winkelableitungen mit Schranke der siebten Ableitung. Der letzte Lauf wurde vor Abschluss unterbrochen.

`experimente/k2_boundary.py` enthält den letzten Code, einschließlich der inaktiven Funktionen der früheren Varianten.
Die unveränderte Ursprungsversion steht in `archiv/k2_boundary_original.py`.
Frühere Zwischenstände sind keine separat gespeicherten Versionen; vollständige Laufprotokolle wurden nicht als Dateien gesichert.
Es gibt deshalb keine `k2_right.json` oder `k2_left.json` aus einem abgeschlossenen eigenen Lauf.

Noch zu prüfen:

- Gültigkeit der reellen und imaginären Chartgrenzen für sämtliche verwendeten Kreisscheiben und verschobenen Argumente.
- Identifikation und lückenlose Verklebung des inversen Zweigs entlang des ganzen Bogens.
- Gültigkeit aller Cauchy-Disk-Einschließungen und der Restbehandlung am exakten Eckwinkel.
- Vollständige Überdeckung beider oberen Teilbögen, danach die für den ganzen Rand nötige Symmetrie und Normierung.
- Voraussetzungen für sechs partielle Integrationen, einschließlich absoluter Stetigkeit und Randtermen.

Die im Skript vorhandene Prüfung imaginärer Chartgrenzen schließt diese offenen Punkte nicht automatisch.
Insbesondere sind weder C₆ ≤ 300 noch die vollständige Bulk-Schranke hier neu bewiesen.

## 4. Literatur

Eine Carleman-Verschiebungsarbeit wurde nur anhand ihres Verlagsabstracts gesichtet (Verweis im Bericht).
Es wurde kein übertragbarer externer Satz gefunden oder als Fundament verwendet; der Literaturhinweis ist kein Beweisbaustein.
