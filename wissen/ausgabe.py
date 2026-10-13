"""Prüfbericht (Markdown) und die eigenständige HTML-Ansicht."""
import base64
import json
import re
from collections import defaultdict
from pathlib import Path


def bericht_md(g, stat, befunde):
    z = []
    z.append("# Prüfbericht Tetrations-Wissensgraph\n")
    z.append(f"{len(g.knoten)} Aussagen · {len(g.kanten)} Kanten · {len(g.experimente)} Experimente\n")
    st = {k.split(":", 1)[1]: v for k, v in stat.items() if k.startswith("status:")}
    z.append("**Effektiver Status:** " + ", ".join(f"{k} {v}" for k, v in sorted(st.items(), key=lambda x: -x[1])) + "\n")
    z.append(f"**Zitate nicht am Original geprüft:** {stat.get('zitate_ungeprueft', 0)}\n")
    nach_stufe = defaultdict(list)
    for b in befunde:
        nach_stufe[b["stufe"]].append(b)
    for stufe in ("FEHLER", "WARNUNG", "HINWEIS"):
        z.append(f"\n## {stufe} ({len(nach_stufe[stufe])})\n")
        regeln = defaultdict(list)
        for b in nach_stufe[stufe]:
            regeln[b["regel"]].append(b)
        for r, bs in sorted(regeln.items()):
            z.append(f"\n**{r}** ({len(bs)})\n")
            for b in bs:
                z.append(f"- `{b['id']}` — {b['text']}")
    z.append("\n## Beweislücken\n")
    for i, a in g.knoten.items():
        for l in (a.get("beweis") or {}).get("luecken") or []:
            z.append(f"- `{i}` — {l}")
    z.append("\n## Bedingte Aussagen (Wurzelbedingungen = unbewiesene Enden der Beweiskette)\n")
    for i, a in g.knoten.items():
        if a["_status_eff"] == "bedingt":
            w = a.get("_wurzelbedingungen") or []
            extra = f" · zusätzlich: {a['bedingungen_text']}" if a.get("bedingungen_text") else ""
            z.append(f"- `{i}` — {', '.join(w) or '—'}{extra}")
    z.append("\n## Experimente\n")
    for e in g.experimente.values():
        erg = e.get("_ergebnis") or {}
        n = len(erg.get("pruefungen") or [])
        f = len(erg.get("fehlgeschlagen") or [])
        z.append(f"- `{e['id']}` — {e.get('_zustand')}, {n - f}/{n} Prüfungen. {erg.get('zusammenfassung', '')}")
    # GitHub stellt \operatorname nicht dar (KaTeX in der Ansicht schon) — im Bericht \mathrm verwenden
    return re.sub(r"\\operatorname\*?\{([^}]*)\}", r"\\mathrm{\1}", "\n".join(z) + "\n")


def _katex_css_inline(vendor: Path) -> str:
    css = (vendor / "katex.min.css").read_text(encoding="utf-8")

    def ersetze(m):
        name = m.group(1)
        pfad = vendor / "fonts" / name
        if pfad.suffix == ".woff2" and pfad.is_file():
            data = base64.b64encode(pfad.read_bytes()).decode()
            return f"url(data:font/woff2;base64,{data}) format(\"woff2\")"
        return ""

    css = re.sub(r"url\(fonts/([^)]+)\)\s*format\(\"[a-z0-9]+\"\)", ersetze, css)
    css = re.sub(r",\s*(?=[;}])", "", css)         # leere Fallbacks (woff/ttf) entfernen
    return css


def _app_js(ansicht: Path) -> str:
    """Alle Module aus ansicht/js in Dateireihenfolge, in einer gemeinsamen Funktion."""
    teile = [f.read_text(encoding="utf-8") for f in sorted((ansicht / "js").glob("*.js"))]
    return "(() => {\n\"use strict\";\n" + "\n".join(teile) + "\n})();"


def html_bauen(wurzel: Path, daten: dict) -> str:
    ansicht = wurzel / "ansicht"
    vendor = ansicht / "vendor" / "katex"
    vorlage = (ansicht / "vorlage.html").read_text(encoding="utf-8")
    teile = {
        "/*__KATEX_CSS__*/": _katex_css_inline(vendor),
        "/*__STIL__*/": (ansicht / "stil.css").read_text(encoding="utf-8"),
        "/*__KATEX_JS__*/": (vendor / "katex.min.js").read_text(encoding="utf-8"),
        "/*__APP_JS__*/": _app_js(ansicht),
        "/*__DATEN__*/": "window.GRAPH = " + json.dumps(daten, ensure_ascii=False).replace("</", "<\\/") + ";",
    }
    for platz, inhalt in teile.items():
        vorlage = vorlage.replace(platz, inhalt)
    return vorlage
