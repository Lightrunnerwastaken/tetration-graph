// Rauchtest des DE/EN-Umschalters: beide Sprachen laden, Reiter und Kapitelüberschrift prüfen, Knopf klicken.
//   NODE_PATH=/opt/node22/lib/node_modules node ansicht/pruefe_umschalter.js ausgabe/tetration-graph.html
const { chromium } = require('playwright');
(async () => {
  const b = await chromium.launch(); const errs = [];
  for (const lang of ['de', 'en']) {
    const p = await b.newPage({ viewport: { width: 1400, height: 900 } });
    p.on('console', m => { if (m.type() === 'error') errs.push(lang + ': ' + m.text()); });
    p.on('pageerror', e => errs.push(lang + ' PAGEERROR ' + e.message));
    await p.goto('file://' + require('path').resolve(process.argv[2] || 'ausgabe/tetration-graph.html') + '?lang=' + lang);
    await p.waitForTimeout(700);
    console.log(lang, 'lang attr', await p.getAttribute('html', 'lang'), '| tabs', await p.$$eval('nav.ansichten button', bs => bs.map(x => x.textContent.trim()).join(', ')), '| Knopf', await p.textContent('#sprache-knopf'));
    await p.click('nav.ansichten button[data-ansicht="theorie"]'); await p.waitForTimeout(500);
    console.log(lang, 'Theorie-Überschrift:', (await p.$$eval('#buch h1, #buch h2, #buch .lb-kapitel-kopf', els => els.slice(0,3).map(e => e.textContent.trim().slice(0, 60)))).join(' | '));
    if (lang === 'de') { await p.click('#sprache-knopf'); await p.waitForTimeout(800); console.log('nach Klick auf EN-Knopf, Tab 1:', await p.$eval('nav.ansichten button', b => b.textContent.trim()), '| url', p.url().split('?')[1]); }
    await p.close();
  }
  console.log('Fehler:', errs); await b.close();
})();
