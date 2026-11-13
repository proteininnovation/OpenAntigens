// Browser walk-through of a built OpenAntigens site, both portals.
//
//   NODE_PATH=<dir with playwright> node scripts/check_regression_site.cjs <public_site> [out_dir]
//
// Serves <public_site> over local HTTP (as the web server would), blocks every outside request, and
// visits the static pages and every report page listed in each portal's index. On each report:
// no page or console errors; no "Traceback"/"undefined"/"NaN"/"[object Object]" in the text; the
// construct cards match the index's construct_count and every construct tab shows cards; the 3D
// viewer reaches "viewer ready" where the index says a structure exists; the boundary editor
// applies; every relative link resolves to a file. Writes a screenshot per report and summary.json
// to out_dir (default <public_site>/../browser_check) and exits non-zero on any failure.
const { chromium } = require('playwright');
const fs = require('node:fs');
const http = require('node:http');
const path = require('node:path');
const zlib = require('node:zlib');

const site = path.resolve(process.argv[2] || '');
const outDir = path.resolve(process.argv[3] || path.join(site, '..', 'browser_check'));
if (!fs.existsSync(path.join(site, 'index.html'))) {
  console.error(`usage: check_regression_site.cjs <public_site> [out_dir]; no index.html in ${site}`);
  process.exit(2);
}
const STATIC_PAGES = ['index.html', 'builder.html', 'constructs.html', 'help.html', 'methods.html', 'downloads.html', 'calculator.html'];
const TYPES = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.json': 'application/json',
  '.svg': 'image/svg+xml', '.png': 'image/png', '.woff2': 'font/woff2', '.pdb': 'chemical/x-pdb', '.cif': 'chemical/x-cif',
  '.txt': 'text/plain', '.tsv': 'text/tab-separated-values', '.gz': 'application/gzip' };

function serve(root) {
  const server = http.createServer((req, res) => {
    const url = new URL(req.url, 'http://localhost');
    let file = path.join(root, decodeURIComponent(url.pathname));
    if (!file.startsWith(root)) { res.writeHead(403); return res.end(); }
    if (fs.existsSync(file) && fs.statSync(file).isDirectory()) file = path.join(file, 'index.html');
    if (!fs.existsSync(file)) { res.writeHead(404); return res.end(); }
    res.writeHead(200, { 'Content-Type': TYPES[path.extname(file)] || 'application/octet-stream' });
    fs.createReadStream(file).pipe(res);
  });
  return new Promise(resolve => server.listen(0, '127.0.0.1', () => resolve(server)));
}

function loadIndex(portalDir) {
  const file = path.join(portalDir, 'downloads', 'agdesign2_portal_index.json');
  if (!fs.existsSync(file)) return null;
  let data = fs.readFileSync(file);
  if (data[0] === 0x1f && data[1] === 0x8b) data = zlib.gunzipSync(data);
  return JSON.parse(data.toString('utf8'));
}

(async () => {
  fs.mkdirSync(outDir, { recursive: true });
  const server = await serve(site);
  const base = `http://127.0.0.1:${server.address().port}/`;
  const browser = await chromium.launch({
    headless: true,
    executablePath: process.env.CHROME_EXECUTABLE || undefined,
    args: ['--enable-unsafe-swiftshader'],
  });
  const failures = [];
  const visited = { human: 0, mouse: 0 };
  const fail = (portal, page, message) => failures.push({ portal, page, message });
  try {
    const context = await browser.newContext({ viewport: { width: 1440, height: 1100 } });
    await context.route('**/*', route => (route.request().url().startsWith(base) ? route.continue() : route.abort()));
    const page = await context.newPage();
    let errors = [];
    page.on('pageerror', error => errors.push(`pageerror: ${error.message}`));
    page.on('console', message => {
      if (message.type() === 'error' && !/net::ERR_FAILED|Failed to load resource/.test(message.text())) {
        errors.push(`console: ${message.text()}`);
      }
    });

    for (const [portal, prefix] of [['human', ''], ['mouse', 'mouse/']]) {
      const portalDir = path.join(site, prefix);
      const index = loadIndex(portalDir);
      if (!index) { fail(portal, 'downloads/agdesign2_portal_index.json', 'portal index missing'); continue; }
      for (const name of STATIC_PAGES) {
        if (!fs.existsSync(path.join(portalDir, name))) { fail(portal, name, 'page missing'); continue; }
        errors = [];
        await page.goto(base + prefix + name, { waitUntil: 'load' });
        errors.forEach(error => fail(portal, name, error));
      }
      await page.goto(base + prefix + 'index.html', { waitUntil: 'load' });
      await page.locator('#searchBox').fill(index[0].gene_symbol || index[0].entry_name);
      await page.waitForTimeout(300);
      if (await page.locator('#resultsTable tbody tr').count() < 1) fail(portal, 'index.html', 'search finds nothing for the first target');

      fs.mkdirSync(path.join(outDir, portal), { recursive: true });
      for (const row of index) {
        const name = row.detail_page;
        if (!name) { fail(portal, row.entry_name, 'index row has no detail page'); continue; }
        errors = [];
        visited[portal]++;
        await page.goto(base + prefix + name, { waitUntil: 'load' });
        const text = await page.locator('body').innerText();
        for (const bad of ['Traceback', '[object Object]']) {
          if (text.includes(bad)) fail(portal, name, `page text contains "${bad}"`);
        }
        // A leaked missing value fills a whole cell or label; methods copy may use the words in prose.
        const leaked = await page.$$eval('body *', elements => elements
          .filter(el => !el.children.length && el.offsetParent !== null)
          .map(el => el.textContent.trim())
          .filter(value => ['undefined', 'NaN', 'null', 'None', 'NaN%', 'undefined aa'].includes(value)));
        for (const value of new Set(leaked)) fail(portal, name, `an element shows "${value}"`);
        for (const id of ['overview', 'notes', 'homology']) {
          if (await page.locator(`#${id}`).count() === 0) fail(portal, name, `section #${id} missing`);
        }
        const cards = await page.locator('#construct-details .construct-card').count();
        if (Number(row.construct_count || 0) !== cards) fail(portal, name, `index says ${row.construct_count} constructs, page shows ${cards}`);
        const tabs = page.locator('#construct-details [role="tab"]');
        for (let i = 0; i < await tabs.count(); i++) {
          await tabs.nth(i).click();
          if (await page.locator('#construct-details .construct-card:visible').count() < 1) {
            fail(portal, name, `construct tab "${(await tabs.nth(i).innerText()).trim()}" shows no cards`);
          }
        }
        if (Number(row.has_alphafold_structure || 0) && await page.locator('#viewerInfo').count()) {
          try {
            await page.waitForFunction(() => /viewer ready/i.test(document.getElementById('viewerInfo').textContent), null, { timeout: 20000 });
          } catch { fail(portal, name, `3D viewer not ready: ${(await page.locator('#viewerInfo').innerText()).slice(0, 120)}`); }
        }
        const editor = page.locator('.construct-workbench-boundary-editor input');
        if (await editor.count() >= 2) {
          await editor.nth(0).fill('2');
          await editor.nth(1).fill('5');
          await page.getByRole('button', { name: 'Apply boundaries' }).first().click();
          if (!/2-5\t/.test(await page.locator('#selectionTsv').inputValue())) fail(portal, name, 'boundary editor did not apply 2-5');
        }
        const hrefs = await page.$$eval('a[href]', anchors => anchors.map(a => a.getAttribute('href')));
        for (const href of new Set(hrefs)) {
          if (!href || /^(https?:|mailto:|tel:|#|javascript:|data:)/.test(href)) continue;
          const target = path.normalize(path.join(portalDir, path.dirname(name), href.split(/[?#]/)[0]));
          if (!fs.existsSync(target)) fail(portal, name, `broken link ${href}`);
        }
        await page.screenshot({ path: path.join(outDir, portal, path.basename(name, '.html') + '.png') });
        errors.forEach(error => fail(portal, name, error));
      }
    }
  } finally {
    await browser.close();
    server.close();
  }
  const summary = { site, visited, failures };
  fs.writeFileSync(path.join(outDir, 'summary.json'), JSON.stringify(summary, null, 2));
  console.log(`visited ${visited.human} human and ${visited.mouse} mouse report pages; ${failures.length} failures`);
  for (const item of failures.slice(0, 60)) console.log(`  [${item.portal}] ${item.page}: ${item.message}`);
  if (failures.length > 60) console.log(`  ... ${failures.length - 60} more in ${path.join(outDir, 'summary.json')}`);
  process.exit(failures.length ? 1 : 0);
})().catch(error => { console.error(error); process.exit(2); });
