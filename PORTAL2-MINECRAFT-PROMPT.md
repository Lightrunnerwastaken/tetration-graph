# Auftrag: Portal 2 in Minecraft (Fabric-Mod)

## Ziel
Ich möchte Minecraft spielen, aber mit den Mechaniken von Portal 2: Source-Engine-Movement, Portal-Gun und die Portal-2-Spielelemente. Alles, was Minecraft ausmacht, soll erhalten bleiben: Blöcke abbauen und platzieren, Inventar, Crafting, Überleben und den Enderdrachen besiegen.

Es geht nicht darum, ein paar vorhandene Mods zusammenzustecken. Die Mechaniken sollen möglichst genau dem Original entsprechen, so wie andere Leute Mechaniken aus einem Spiel (z. B. Mario) in ein anderes Spiel übertragen haben. Die Mod soll dafür die echten Formeln und Konstanten verwenden und nicht nur ungefähr so wirken.

## Vorgehen

### 1. Recherche (kurz)
Prüfe, was es schon gibt, und nutze Vorhandenes als Referenz oder Ausgangspunkt, soweit die Lizenz es erlaubt:
- Source-Movement für Minecraft: Bunnyhop, Minehop/MinehopLite, MomentumReforged, BunnyHop+
- Portal-Gun: Portal Gun (Fabric) mit Immersive Portals, PortalGunMod (Fabric 1.21.11), Portal Gun (Reforged)
- Gels: Portal Gels (NeoForge), Modular Gel (Forge)
- Nicht gefunden: Excursion Funnels, Thermal Discouragement Beams, Lichtbrücken, Turrets, Würfel, Knöpfe

### 2. Werte aus dem Original ermitteln
- **Source-Movement:** Quelle ist `gamemovement.cpp` aus dem öffentlichen Source SDK 2013. Lies Bodenbeschleunigung, Air-Accelerate, Reibung, Stopp-Geschwindigkeit, Sprunggeschwindigkeit, Stufenhöhe, Gravitation und Duck-Verhalten direkt aus dem Quellcode.
- **Portal-spezifisches** (Portal-Transformation, Impulserhalt, Platzierungsregeln, Gels, Funnels, Lichtbrücken): Analysiere die Portal-2-Binärdateien (`server.dll`, `client.dll`) mit REA (https://github.com/morluto/rea, Installation mit `npx rea-agents setup`). Dafür muss Portal 2 lokal installiert sein. Leite jeden Wert aus der Analyse ab und belege ihn.

### 3. Umsetzung
- Fabric-Mod für Minecraft 1.21.11 (Java). Wenn eine andere Version sinnvoller ist, begründe es kurz.
- Reihenfolge: (a) Source-Movement, (b) Portal-Gun mit Portalen und Impulserhalt, (c) Gels, (d) Funnels und Lichtbrücken, (e) Rest wie Würfel, Knöpfe und Turrets.
- Die Mod enthält nur eigenen Code und Zahlenwerte. Keine Texturen, Modelle oder Sounds von Valve ins Repo. Verwende eigene oder frei lizenzierte Assets.
- Das Ergebnis muss in einem normalen Survival-Spiel laufen. Singleplayer reicht zunächst, Multiplayer später.

### 4. Überprüfung
Vergleiche die Mod mit dem Original anhand von Messwerten:
- Geschwindigkeit nach dem Durchgehen eines Portals, bei verschiedenen Portalwinkeln
- Sprunghöhe und -weite, Air-Strafe-Beschleunigung
- Sprung auf Repulsion-Gel, Beschleunigung auf Propulsion-Gel

Schreibe die Messwerte und Abweichungen in eine Datei im Repo.

## Regeln für die Arbeit
- Starte erst mit dem Movement und zeig mir, dass es sich wie in Source anfühlt, bevor du Portale baust.
- Behaupte nichts, was du nicht getestet hast. Sag klar, was ungetestet ist oder nicht funktioniert.
- Wenn dir für eine Analyse etwas fehlt (Spieldateien, Zugriff), frag nach, statt zu raten.
- Nutze ein eigenes neues Repo. Das Repo `tetration-graph` hat damit nichts zu tun.

## Hinweis zu den Beispielen
Ich habe gehört, dass Leute Mario-Mechaniken in Minecraft übertragen haben, aber keine Quelle dafür. Wenn du konkrete Projekte findest, sind sie als Vorbild hilfreich.
