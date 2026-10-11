"""Erlaubte Werte und Bedeutungen. Einzige Stelle, an der das Vokabular steht."""

ARTEN = {
    "definition": "Definition",
    "satz": "Satz",
    "lemma": "Lemma",
    "korollar": "Korollar",
    "konstruktion": "Konstruktion",
    "annahme": "Annahme / Hypothese",
    "vermutung": "Vermutung",
    "beobachtung": "Beobachtung",
    "gegenbeispiel": "Gegenbeispiel",
    "gescheitert": "Gescheiterter Ansatz",
    "korrektur": "Korrektur",
    "frage": "Offene Frage",
}

# Arten, deren Inhalt eine mathematische Behauptung mit Beweisanspruch ist
BEHAUPTUNGEN = {"satz", "lemma", "korollar", "konstruktion", "gegenbeispiel", "korrektur"}
# Arten, für die die fünf Prüffragen Pflicht sind
PRUEFPFLICHT = BEHAUPTUNGEN | {"vermutung", "beobachtung"}

EBENEN = {
    "A": "allgemeine Funktionalgleichung (jede Lösung)",
    "B": "Klasse von Lösungen (durch Bedingungen beschrieben)",
    "C": "konkrete Konstruktion (z. B. Kneser, regulär an einem Fixpunkt)",
    "X": "Hilfssatz aus einem anderen Gebiet",
}

# rang: Beweiskraft für die Weitergabe. 5 = trägt einen Beweis.
STATUS = {
    "definition":         {"label": "Definition",                   "rang": 5},
    "bewiesen_literatur": {"label": "bewiesen (Literatur)",         "rang": 5},
    "bewiesen":           {"label": "bewiesen (hier/eigene Arbeit)", "rang": 5},
    "entlastet":          {"label": "Annahme entlastet",            "rang": 5},
    "bedingt":            {"label": "bedingt bewiesen",             "rang": 3},
    "skizze":             {"label": "Beweisskizze, Lücken offen",   "rang": 2},
    "numerisch":          {"label": "numerisch beobachtet",         "rang": 1},
    "vermutung":          {"label": "Vermutung",                    "rang": 1},
    "annahme":            {"label": "Annahme (unbewiesen)",         "rang": 1},
    "offen":              {"label": "offen",                        "rang": 0},
    "gescheitert":        {"label": "gescheitert",                  "rang": -1},
    "widerlegt":          {"label": "widerlegt",                    "rang": -2},
    "ueberholt":          {"label": "überholt (ersetzt)",           "rang": -2},
}
BEWIESEN = {"bewiesen", "bewiesen_literatur", "definition", "entlastet"}
TOT = {"widerlegt", "ueberholt"}

NEUHEIT = {
    "klassisch": "klassisch, lange bekannt",
    "bekannt": "bekannt aus der Literatur",
    "folklore": "vermutlich bekannt, ohne gefundene Quelle",
    "vermutlich_neu": "nach dokumentierter Recherche vermutlich neu",
    "neu_ungeprueft": "als neu bezeichnet, Recherche nicht dokumentiert",
    "nur_schreibweise": "neue Schreibweise, kein neuer Inhalt",
    "unklar": "unklar",
}
NUTZEN = {"hoch", "mittel", "gering", "unklar"}

QUELLEN_PRUEFUNG = {
    "original": "am Original geprüft",
    "sekundaer": "über Zusammenfassung/Rezension geprüft",
    "nein": "nicht geprüft (aus dem Gedächtnis)",
}
QUELLEN_ARTEN = {"literatur", "eigene", "software", "forum", "notiz"}

BENUTZT = {"ja", "nein", "unklar"}

# Kantentypen: Gewicht für die Verwandtschaft, ob sie Beweiskraft trägt
KANTEN = {
    "stuetzt_sich_auf":      {"label": "stützt sich auf (Beweis)", "gewicht": 3.0, "beweis": True},
    "folgt_unmittelbar_aus": {"label": "folgt unmittelbar aus",    "gewicht": 3.0, "beweis": True},
    "spezialfall_von":       {"label": "Spezialfall von",          "gewicht": 3.0, "beweis": False},
    "verallgemeinert":       {"label": "verallgemeinert",          "gewicht": 2.5, "beweis": False},
    "ersetzt":               {"label": "ersetzt (neuere Fassung)", "gewicht": 3.0, "beweis": False},
    "korrigiert":            {"label": "korrigiert",               "gewicht": 3.0, "beweis": False},
    "widerlegt":             {"label": "widerlegt",                "gewicht": 3.0, "beweis": False},
    "entlastet":             {"label": "entlastet (beweist Annahme)", "gewicht": 3.0, "beweis": False},
    "stuetzt":               {"label": "stützt (Evidenz)",         "gewicht": 2.0, "beweis": False},
    "uebertragbar_auf":      {"label": "übertragbar auf",          "gewicht": 2.0, "beweis": False},
    "analog_zu":             {"label": "analog zu",                "gewicht": 1.5, "beweis": False},
    "motiviert":             {"label": "motiviert",                "gewicht": 1.0, "beweis": False},
    "verwandt":              {"label": "verwandt",                 "gewicht": 1.0, "beweis": False},
}

PRUEFFRAGEN = {
    "f1_folgt_unmittelbar": "1. Folgt es unmittelbar aus bekannten Aussagen?",
    "f2_benutzte_voraussetzungen": "2. Welche Voraussetzungen benutzt der Beweis wirklich?",
    "f3_spezialfall": "3. Spezialfall eines allgemeineren Satzes?",
    "f4_uebertragbar": "4. Übertragbar auf andere Konstruktionen?",
    "f5_neue_frage": "5. Welche neue Frage, Abschätzung oder Berechnung ermöglicht es?",
}

PFLICHTFELDER = ["id", "titel", "art", "ebene", "status", "aussage"]

# Für wen eine Aussage gedacht ist (Bibliothek für die Öffentlichkeit).
NIVEAU = {
    "einstieg": "Einstieg — Schulmathematik plus Ableitung und Integral",
    "vertiefung": "Vertiefung — Studium (Funktionentheorie, Iteration)",
    "forschung": "Forschung — offene Fragen, eigene Arbeit, Spezialliteratur",
}
# Felder, die eine spätere englische Fassung unter `en:` tragen darf.
UEBERSETZBAR = {"titel", "aussage", "bereich", "einfach", "beweis", "formen", "naechste_frage", "gegenbeispielsuche"}

# Lehrbuch: Teile eines Theorie-Blocks (daten/theorie/<gebiet>.yaml), in Lesereihenfolge.
THEORIE_TEILE = {
    "worum": "Worum es geht",
    "genau": "Genau",
    "warum": "Warum es gilt",
    "beispiel": "Beispiel",
    "bilder": "Bilder",
    "einordnung": "Einordnung",
}
THEORIE_PFLICHT = ["worum", "genau", "warum", "einordnung"]
