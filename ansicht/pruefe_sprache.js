// Prüft die englische Ansicht auf Restdeutsch: klickt sich durch alle Reiter und sammelt Texte, die noch deutsch aussehen.
//   NODE_PATH=/opt/node22/lib/node_modules node ansicht/pruefe_sprache.js ausgabe/tetration-graph.html [rest.json] [--de]
// Ohne Chromium-Pfad nimmt playwright seinen eigenen (PLAYWRIGHT_BROWSERS_PATH). Rückgabewert 1, wenn Reste oder Fehler gefunden wurden.
const { chromium } = require("playwright");
const fs = require("fs");
const pfad = require("path").resolve(process.argv[2] || "ausgabe/tetration-graph.html");
const ausgabe = process.argv[3] && !process.argv[3].startsWith("--") ? process.argv[3] : null;
const SPR = process.argv.includes("--de") ? "de" : "en";

(async () => {
  const b = await chromium.launch();
  const p = await b.newPage({ viewport: { width: 1500, height: 950 } });
  const fehler = [];
  p.on("console", (m) => { if (m.type() === "error") fehler.push(m.text()); });
  p.on("pageerror", (e) => fehler.push("PAGEERROR " + e.message));
  await p.goto("file://" + pfad + "?lang=" + SPR);
  await p.waitForTimeout(600);
  const rest = new Map(); // Text -> {n, wo}
  const sammle = async (wo) => {
    const funde = await p.evaluate(() => {
      const S = window.__sprache, aus = [];
      const w = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
      while (w.nextNode()) {
        const n = w.currentNode;
        if (/^(SCRIPT|STYLE)$/.test(n.parentNode.nodeName)) continue;
        const s = n.nodeValue.trim();
        if (s && S.sieht_deutsch_aus(s)) aus.push(s);
      }
      for (const el of document.body.querySelectorAll("[title],[aria-label],[placeholder],[alt]"))
        for (const a of ["title", "aria-label", "placeholder", "alt"]) { const v = el.getAttribute(a); if (v && S.sieht_deutsch_aus(v)) aus.push(v); }
      if (S.sieht_deutsch_aus(document.title)) aus.push(document.title);
      return aus.concat(S.canvas.splice(0).map((c) => "[canvas] " + c));
    });
    for (const s of funde) { const k = s.length > 400 ? s.slice(0, 400) + "…" : s; const e = rest.get(k) || { n: 0, wo }; e.n++; rest.set(k, e); }
  };
  const tabs = await p.$$eval("nav.ansichten button", (bs) => bs.map((x) => x.dataset.ansicht));
  for (const tab of tabs) {
    await p.click(`nav.ansichten button[data-ansicht="${tab}"]`);
    await p.waitForTimeout(500);
    await sammle(tab);
    const geklickt = new Set();
    for (let runde = 0; runde < 70; runde++) {
      const h = await p.evaluateHandle((geklickt) => {
        const sicht = (el) => { const r = el.getBoundingClientRect(); const st = getComputedStyle(el); return r.width > 2 && r.height > 2 && st.visibility !== "hidden" && st.display !== "none" && !el.closest("[hidden]"); };
        const sig = (el) => `${el.tagName}|${el.className}|${(el.dataset.id || el.dataset.exp || el.dataset.ansicht || "")}|${(el.textContent || "").trim().slice(0, 30)}`;
        const liste = [...document.querySelectorAll("button, summary, [role=tab], .verweis, .eintrag, li[data-id], [data-id], a[href]")].filter((el) => el.id !== "sprache-knopf" && sicht(el) && !geklickt.includes(sig(el)));
        const el = liste[0];
        if (!el) return null;
        el.setAttribute("data-geklickt-sig", sig(el));
        return el;
      }, [...geklickt]);
      const knoten = h.asElement();
      if (!knoten) break;
      const sig = await knoten.getAttribute("data-geklickt-sig");
      geklickt.add(sig);
      await knoten.click({ timeout: 1500, force: true }).catch(() => {});
      await p.waitForTimeout(120);
      await sammle(tab);
    }
    // Karten: Mausklicks und Zoom auf der Leinwand
    for (const sel of ["#karte", "#ideen-flaeche"]) {
      const sichtbar = await p.$eval(sel, (c) => { const r = c.getBoundingClientRect(); return r.width > 10 && !c.hidden && getComputedStyle(c).display !== "none"; }).catch(() => false);
      if (!sichtbar) continue;
      const box = await (await p.$(sel)).boundingBox();
      for (let i = 0; i < 24; i++) {
        const x = box.x + box.width * (0.15 + 0.7 * ((i * 37) % 100) / 100), y = box.y + box.height * (0.15 + 0.7 * ((i * 61) % 100) / 100);
        await p.mouse.move(x, y); await p.waitForTimeout(40);
        if (i % 3 === 0) { await p.mouse.click(x, y).catch(() => {}); await p.waitForTimeout(150); }
        if (i % 8 === 7) { await p.mouse.wheel(0, -500).catch(() => {}); await p.waitForTimeout(250); }
        await sammle(tab + ":" + sel);
      }
    }
    await p.keyboard.press("Escape").catch(() => {});
  }
  // Suche
  await p.fill("#suche", "Kneser"); await p.waitForTimeout(250); await sammle("suche");
  await p.fill("#suche", "qqqq"); await p.waitForTimeout(250); await sammle("suche-leer");
  await p.fill("#suche", ""); 
  await p.click("#legende-knopf").catch(() => {}); await p.waitForTimeout(200); await sammle("legende");
  await p.click("#tafel-knopf").catch(() => {}); await p.waitForTimeout(200); await sammle("tafel");
  const liste = [...rest.entries()].sort((a, b) => b[1].n - a[1].n);
  if (ausgabe) fs.writeFileSync(ausgabe, JSON.stringify(liste.map(([t, e]) => ({ text: t, n: e.n, wo: e.wo })), null, 1));
  console.log(`${liste.length} deutsch aussehende Texte, ${fehler.length} Konsolenfehler`);
  for (const [t, e] of liste.slice(0, 60)) console.log(`  [${e.wo}] ×${e.n}  ${JSON.stringify(t.slice(0, 120))}`);
  for (const f of fehler.slice(0, 10)) console.log("  FEHLER", f.slice(0, 200));
  await b.close();
  process.exit(liste.length || fehler.length ? 1 : 0);
})();
