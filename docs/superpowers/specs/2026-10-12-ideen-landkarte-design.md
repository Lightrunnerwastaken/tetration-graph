# Ideen-Landkarte — Entwurf (12.10.2026)

## Ziel

Die Karte soll die **mathematische Struktur der Tetration** zeigen, nicht nur Aussagen an Orten mit Verzweigungen.
Zwei Zwecke, Überblick zuerst:

1. **Überblick:** Welche grossen Ideen gibt es, wie hängen sie zusammen, wo ist die Tetration gut erforscht und wo liegt die Grenze?
2. **Einem Ergebnis nachgehen:** Woraus besteht ein Satz oder eine Vermutung, welche Stücke sind sicher, wo fehlt es?

Leser: Janis, spätere Leser der öffentlichen Bibliothek und Forschungs-Agenten (diese über eine Textfassung).

Erfolg heisst: Beim Öffnen sieht man acht benannte Hauptideen, deren Lage die Kopplung wiedergibt; zwei Klicks führen zu einzelnen
Aussagen; Haupt- und Hilfsaussagen sind unterscheidbar; Erforschtheit und Herkunft (eigene Arbeit / Literatur) sind ohne Lesen erkennbar.

## 1. Daten: Ideen und Zuordnung

Neue Datei `daten/ideen.yaml`, zwei Stufen:

```yaml
- id: fortsetzung                 # Hauptidee
  name: Fraktionale Fortsetzung
  kurz: Wie man T(x) zwischen den ganzen Höhen festlegt.
- id: fortsetzung.kneser          # Unteridee
  eltern: fortsetzung
  name: Kneser über η
  kurz: Holomorphe Lösung über die Fixpunkte L, L̄.
  straenge: [kn.konstruktion, kn.rand]    # Vorgabe für alle Aussagen dieser Stränge
```

- Jede Aussage gehört zu genau einer Unteridee: über ihren Strang (Feld `straenge` der Unteridee) oder ausnahmsweise über ein Feld
  `idee: <unteridee>` in der Aussage selbst (überschreibt den Strang).
- Die Stränge bleiben unverändert als Lesewege und für das Lehrbuch.
- Prüfung (`wissen/graph.py`): Aussage ohne Idee → Hinweis „ohne-idee“; unbekannte Idee, Unteridee ohne gültige Eltern, Strang in zwei
  Unterideen → Fehler.

Hauptideen (Ausgangspunkt, von Janis bestätigt):

| Hauptidee | Unterideen (Ausgangspunkt) |
|---|---|
| Die Gleichung T(x+1) = a^T(x) | Begriffe; Freiheit (periodisches θ); fest und frei; Werte, Ableitungen, Integrale; slog; Superwurzeln |
| Dynamik von b^z und der Fixpunkt L | Fixpunkte; anziehend/abstossend (Koenigs); parabolisch (Fatou, Grenzfall η); indifferent (Siegel) |
| Unendlicher Turm und Basen-Ebene | Konvergenz; Periode 2 unter e^{-e}; Regionen und Fraktale; Spiralen und Mittelwerte |
| Fraktionale Fortsetzung | regulär unter η; Grenzfall η; Kneser über η; andere Konstruktionen; Eindeutigkeit; reelle Kennzeichnung; halbe Exponentialfunktion |
| Basiswechsel T_b ∘ slog_d | gleicher Fixpunkt; Kreisabbildung; Phasengesetz und Welle; Log-Peeling und slog(C); Moden; Antwortkern (Riemann–Hilbert); Rekonstruktion |
| Höhere Stufen | Pentation, Hexation; η_n → 2, e_n → 4 (TPID 10); Hyperoperationen |
| Zahlen und Arithmetik | Transzendenz; Kneser–Schanuel; Kongruenzen, Türme modulo m |
| Brücken | Logik und Wachstum; Fourier und Kreis |

Wo ein Strang zwei Richtungen hat, wird er über `idee:`-Felder einzelner Aussagen geteilt.

## 2. Lage (beim Bauen berechnet)

Neues Modul `wissen/landkarte.py`, aufgerufen von `baue.py`; Ergebnis deterministisch (fester Startwert), Koordinaten gehen mit den
Daten in die HTML-Datei.

- **Kopplung** zweier Unterideen U, V:
  w(U,V) = Σ Kanten zwischen ihren Aussagen (stuetzt_sich_auf, verallgemeinert, ersetzt, korrigiert, widerlegt: 1.0; stuetzt, motiviert: 0.5;
  verwandt, analog_zu und übrige: 0.25) + 0.5 · Jaccard(Objekte(U), Objekte(V)) · min(|U|,|V|).
- **Hauptideen:** Kraftlayout auf den summierten Kopplungen (Anziehung ∝ w, Abstossung zwischen allen, Grösse ∝ √Anzahl Aussagen,
  keine Überlappung).
- **Unterideen:** Kraftlayout in der Blase der Hauptidee; zusätzliche schwache Anziehung zu gekoppelten Unterideen anderer Hauptideen
  (dadurch rücken z. B. „Phasengesetz und Welle“ an den Rand zur fraktionalen Fortsetzung). Bleiben innerhalb der Blase.
- **Aussagen:** in der Unteridee nach Tiefe im Abhängigkeitsgraphen: tragende (wenig eigene Voraussetzungen, viele Abhängige) innen,
  aufbauende aussen; Winkel nach Kopplung.

## 3. Ansicht und Bedienung

Die Ideen-Karte wird Startansicht; die heutige Zoom-Karte (Gebiete → Stränge → Aussagen) bleibt als zweite Ansicht wählbar.

1. Start: acht Hauptideen als Kreise, Verbindungslinien mit Dicke ∝ Kopplung.
2. Klick auf Hauptidee: öffnet sich mit Kamerafahrt, zeigt ihre Unterideen; andere Hauptideen klein am Rand, Linien zu gekoppelten
   Unterideen anderswo bleiben.
3. Klick auf Unteridee: zeigt ihre Aussagen. **Hauptaussagen** gross und hell (Art satz/vermutung/frage/konstruktion und hohe
   `_tragweite` oder `nutzen_stufe: hoch`), **Hilfsaussagen** klein und gedämpft (lemma, korollar, beobachtung, definition ohne
   Abhängige). Die Schwelle wird so gewählt, dass je Unteridee höchstens etwa ein Drittel Hauptaussagen sind.
4. Klick auf Aussage: Inhaltskarte wie heute.

Pfadleiste oben („Fraktionale Fortsetzung › Kneser über η“), Zurück mit Esc oder Klick auf den Pfad. Mausrad zoomt nur innerhalb der
offenen Stufe. Suche und Fokusmodus der heutigen Karte springen in die Ideen-Karte an die richtige Stelle.

## 4. Grenzen sichtbar machen

- **Füllstand-Ring** um jede Idee: Anteile nach Status in den bestehenden Statusfarben (bewiesen/Literatur/Definition; bedingt/Skizze;
  numerisch; Vermutung/Annahme; offen hohl; widerlegt/gescheitert rot).
- **Herkunfts-Ring** (dünn, innen): Anteil Aussagen mit Literatur-, Forums- und nur eigenen Quellen.
- **Stücke einer Vermutung:** Beim Öffnen einer Aussage mit Status vermutung/bedingt/skizze zeigt die Inhaltskarte den Baum ihrer
  Voraussetzungen (stuetzt_sich_auf, rekursiv, Tiefe ≤ 3): bewiesen voll, bedingt halb durchsichtig, offen/annahme als hohler Umriss mit
  „hier fehlt es“. Rein aus Kanten und Status berechnet.
- **Textfassung für Agenten:** `baue.py` schreibt `ausgabe/landkarte.md` mit Hauptideen, Unterideen, Füllständen, Herkunft, den
  stärksten Kopplungen und je Unteridee den offenen Aussagen (Vermutungen, Fragen, Annahmen).

## Nicht Teil dieses Entwurfs

Englische Fassung, mobile Bedienung über die heutige hinaus, neue Inhalte im Graphen.

## Prüfung

- `baue.py --nur-pruefen`: 0 Fehler, kein „ohne-idee“-Hinweis.
- Tests für `wissen/landkarte.py`: Kopplung an einem Mini-Graphen von Hand nachgerechnet; Layout deterministisch; keine Überlappung von
  Hauptideen; jede Unteridee in ihrer Blase.
- Im Browser (Vorschau Port 5196): drei Klicks bis zur Inhaltskarte, Pfadleiste, Rückweg, alte Karte weiter erreichbar, keine
  Konsolenfehler; Bildschirmfoto als Beleg.
