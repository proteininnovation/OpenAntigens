// Run after: python -m tests.build_browser_fixture /path/to/fixture
const { chromium } = require('playwright');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const { pathToFileURL } = require('node:url');

(async () => {
  const root = path.resolve(process.argv[2]);
  const browser = await chromium.launch({
    headless: true,
    executablePath: process.env.CHROME_EXECUTABLE || undefined,
    args: ['--enable-unsafe-swiftshader'],
  });
  try {
    const page = await browser.newPage({ viewport: { width: 1440, height: 1100 } });
    await page.addInitScript(() => {
      window.paeAllocations = 0;
      const create = CanvasRenderingContext2D.prototype.createImageData;
      CanvasRenderingContext2D.prototype.createImageData = function (...args) {
        if (this.canvas.id === 'paePlot') window.paeAllocations++;
        return create.apply(this, args);
      };
    });
    const errors = [];
    page.on('pageerror', error => errors.push(error.message));
    await page.route('https://**/*', route => route.abort());
    await page.route('**/3Dmol-min.js', route => route.abort());
    const visit = name => page.goto(pathToFileURL(path.join(root, name)).href);
    await visit('site/reports/test_human.html');
    await page.locator('.construct-mutation-checkbox').check();
    assert.match(await page.locator('.construct-export-tsv').inputValue(), /MMSMMMMMMM/);
    assert.match(await page.locator('.construct-export-fasta').inputValue(), /MMSMMMMMMM/);
    const initialPaeAllocations = await page.evaluate(() => window.paeAllocations);
    const inputs = page.locator('.construct-workbench-boundary-editor input');
    await inputs.nth(0).fill('2');
    await inputs.nth(1).fill('5');
    await page.getByRole('button', { name: 'Apply boundaries' }).click();
    assert.match(await page.locator('#selectionTsv').inputValue(), /2-5\tMCMM/);
    await page.evaluate(() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))));
    assert.equal(await page.evaluate(() => window.paeAllocations), initialPaeAllocations);
    assert.equal(await page.locator('#injected-annotation').count(), 0);
    assert.equal(await page.evaluate(() => Boolean(window.annotationExecuted)), false);
    await page.locator('.selection-mutation-checkbox').check();
    assert.match(await page.locator('#selectionFasta').inputValue(), /MSMM/);
    assert.equal(await page.locator('.construct-workbench-species-row').count(), 1);
    assert.equal(await page.locator('#constructWorkbenchPlddtPlot').count(), 0);

    // The same UI also operates with the real bundled 3D library.
    await page.unroute('**/3Dmol-min.js');
    await visit('site/reports/test_human.html');
    await page.waitForFunction(() => document.getElementById('viewerInfo').textContent.includes('viewer ready'));
    for (const width of [1440, 1920, 320, 390]) {
      await page.setViewportSize({ width, height: 1100 });
      assert.equal(await page.locator('.construct-workbench-current').evaluate(panel => {
        const bounds = panel.getBoundingClientRect();
        const controls = [...panel.querySelectorAll('.selection-live-header button, .selection-export-help')];
        return controls.length === 3 && controls.every(control => {
          const rect = control.getBoundingClientRect();
          return rect.left >= bounds.left && rect.right <= bounds.right && rect.height < 60;
        });
      }), true, `builder export controls fit at ${width}px`);
    }
    await page.setViewportSize({ width: 1440, height: 1100 });
    await page.getByRole('button', { name: 'Apply boundaries' }).click();
    await page.locator('#interactive-construct-builder').screenshot({ path: path.join(root, 'builder.png') });

    await visit('tabs-site/reports/test_human.html');
    await page.waitForFunction(() => document.getElementById('viewerInfo').textContent.includes('viewer ready'));
    const details = page.locator('#construct-details');
    const tabs = details.getByRole('tab');
    assert.deepEqual(await tabs.allTextContents(), [
      'Full ectodomain (1)', 'PDB (2)', 'Annotated domains / repeats (3)', 'Strict (2)',
      'Lenient (1)', 'Full-length multipass (1)', 'Trimmed multipass (4)',
    ]);
    const categoryCounts = [1, 2, 3, 2, 1, 1, 4];
    for (let index = 0; index < categoryCounts.length; index++) {
      await tabs.nth(index).click();
      assert.equal(await details.getByRole('tabpanel').count(), 1);
      assert.equal(await details.locator('.construct-card:visible').count(), categoryCounts[index]);
    }
    assert.equal(await details.locator('.construct-card').count(), 14);
    await tabs.nth(1).click();
    assert.equal(await details.getByRole('heading', { name: 'Soluble PDB constructs' }).isVisible(), true);
    assert.equal(await details.getByRole('heading', { name: 'Membrane PDB constructs' }).isVisible(), true);
    await tabs.nth(0).click();
    const baseline = details.getByRole('tabpanel');
    await baseline.locator('.construct-mutation-checkbox').check();
    const editedTsv = await baseline.locator('.construct-export-tsv').inputValue();
    const editedFasta = await baseline.locator('.construct-export-fasta').inputValue();
    assert.match(editedFasta, /MMSMMMMMMM/);
    await tabs.nth(3).click();
    await details.getByRole('button', { name: 'View in 3D structure' }).first().click();
    assert.match(await page.locator('#selectionTsv').inputValue(), /MMCM/);
    const viewerInfo = await page.locator('#viewerInfo').textContent();
    await tabs.nth(0).click();
    assert.equal(await baseline.locator('.construct-mutation-checkbox').isChecked(), true);
    assert.equal(await baseline.locator('.construct-export-tsv').inputValue(), editedTsv);
    assert.equal(await baseline.locator('.construct-export-fasta').inputValue(), editedFasta);
    assert.equal(await page.locator('#viewerInfo').textContent(), viewerInfo);
    await tabs.nth(0).focus();
    for (const [key, selected] of [['ArrowLeft', 6], ['ArrowRight', 0], ['End', 6], ['Home', 0], ['ArrowRight', 1]]) {
      await page.keyboard.press(key);
      assert.equal(await tabs.nth(selected).getAttribute('aria-selected'), 'true');
      assert.equal(await tabs.nth(selected).evaluate(tab => tab === document.activeElement), true);
      assert.equal(await details.locator('[role="tab"][tabindex="0"]').count(), 1);
    }
    await details.screenshot({ path: path.join(root, 'construct-tabs-desktop.png') });
    await page.setViewportSize({ width: 390, height: 844 });
    await tabs.nth(6).click();
    assert.equal(await details.evaluate(section => section.getBoundingClientRect().right <= innerWidth), true);
    assert.equal(await details.evaluate(section => section.scrollWidth <= section.clientWidth), true);
    assert.equal(await details.locator('.construct-tabs').evaluate(bar => bar.scrollWidth > bar.clientWidth), true);
    assert.equal(await details.getByRole('tabpanel').locator('.selection-table-wrap').first().evaluate(table => table.scrollWidth > table.clientWidth), true);
    assert.equal(await tabs.nth(6).isVisible(), true);
    await details.locator('h2').scrollIntoViewIfNeeded();
    await page.screenshot({ path: path.join(root, 'construct-tabs-mobile.png') });
    await page.setViewportSize({ width: 1440, height: 1100 });

    await page.emulateMedia({ reducedMotion: 'reduce' });
    const sectionButton = page.getByRole('button', { name: 'Sections', exact: true });
    const sectionMenu = page.locator('#report-sections-menu');
    assert.equal(await sectionMenu.locator('a').count(), 14);
    assert.equal(await page.evaluate(() => getComputedStyle(document.documentElement).scrollBehavior), 'auto');
    const selectedCategory = await details.locator('[role="tab"][aria-selected="true"]').textContent();
    const selectedSequence = await page.locator('#selectionTsv').inputValue();
    for (const id of ['interactive-construct-builder', 'construct-details']) {
      await sectionButton.click();
      await sectionMenu.locator(`a[href="#${id}"]`).click();
      await page.waitForFunction(id => document.activeElement === document.getElementById(id).querySelector('h2'), id);
      assert.equal(await sectionMenu.isVisible(), false);
    }
    assert.equal(await details.locator('[role="tab"][aria-selected="true"]').textContent(), selectedCategory);
    assert.equal(await page.locator('#selectionTsv').inputValue(), selectedSequence);
    assert.equal(await page.locator('#construct-panel-0 .construct-export-tsv').inputValue(), editedTsv);

    const helpLink = page.getByRole('link', { name: 'Don’t know how to pick your construct?' });
    await helpLink.focus();
    const popupPromise = page.waitForEvent('popup');
    await page.keyboard.press('Enter');
    const guidePage = await popupPromise;
    await guidePage.waitForLoadState();
    assert.equal(new URL(guidePage.url()).hash, '#choose-constructs');
    assert.equal(await guidePage.locator('main section').first().getAttribute('id'), 'choose-constructs');
    const guide = guidePage.locator('#choose-constructs');
    await guidePage.waitForFunction(() => Math.abs(document.getElementById('choose-constructs').getBoundingClientRect().top - 16) < 2);
    assert.equal(await guidePage.locator('.construct-choice-workflow li').count(), 4);
    const desktopSteps = await guidePage.locator('.construct-choice-workflow li').evaluateAll(steps => steps.map(step => step.getBoundingClientRect().top));
    assert.equal(new Set(desktopSteps).size, 1);
    await guide.screenshot({ path: path.join(root, 'construct-guide-desktop.png') });
    for (const width of [320, 390]) {
      await guidePage.setViewportSize({ width, height: 844 });
      assert.equal(await guidePage.evaluate(() => document.documentElement.scrollWidth <= innerWidth), true);
      assert.equal(await guide.evaluate(section => section.scrollWidth <= section.clientWidth), true);
      assert.equal(await guidePage.locator('.construct-choice-workflow li').evaluateAll(steps => steps.every((step, i) => !i || step.getBoundingClientRect().top > steps[i - 1].getBoundingClientRect().bottom)), true);
      await guide.screenshot({ path: path.join(root, `construct-guide-${width}.png`) });
    }
    for (const hash of ['#construct-classes', '#construct-feature-review']) {
      await guide.locator(`a[href="${hash}"]`).click();
      assert.equal(new URL(guidePage.url()).hash, hash);
      assert.equal(await guidePage.locator(hash).isVisible(), true);
    }
    await guidePage.close();
    assert.equal(await details.locator('[role="tab"][aria-selected="true"]').textContent(), selectedCategory);
    assert.equal(await page.locator('#selectionTsv').inputValue(), selectedSequence);
    assert.equal(await page.locator('#construct-panel-0 .construct-export-tsv').inputValue(), editedTsv);
    assert.equal(await page.locator('#construct-panel-0 .construct-export-fasta').inputValue(), editedFasta);
    assert.equal(await page.locator('#construct-panel-0 .construct-mutation-checkbox').isChecked(), true);
    assert.equal(await page.locator('#viewerInfo').textContent(), viewerInfo);

    await visit('navigation-site/reports/test_human.html');
    await page.waitForFunction(() => document.querySelectorAll('#report-sections-menu a').length === 18);
    const navigationIds = await sectionMenu.locator('a').evaluateAll(links => links.map(link => link.hash.slice(1)));
    assert.equal(new Set(navigationIds).size, 18);
    const layout = await page.locator('main.page').boundingBox();
    for (const id of navigationIds) {
      await sectionButton.click();
      await sectionMenu.locator(`a[href="#${id}"]`).click();
      await page.waitForFunction(id => document.activeElement === document.getElementById(id).querySelector('h2'), id);
      assert.equal(new URL(page.url()).hash, '#' + id);
      assert.equal(await sectionMenu.isVisible(), false);
      await page.waitForFunction(id => document.querySelector('#report-sections-menu [aria-current="location"]')?.hash === '#' + id, id);
    }
    await page.goBack();
    assert.equal(new URL(page.url()).hash, '#construct-details');
    await page.waitForFunction(() => document.activeElement === document.querySelector('#construct-details h2'));
    await sectionButton.focus();
    await page.keyboard.press('Enter');
    assert.equal(await sectionMenu.isVisible(), true);
    assert.equal((await page.locator('main.page').boundingBox()).width, layout.width);
    await page.keyboard.press('Escape');
    assert.equal(await sectionMenu.isVisible(), false);
    assert.equal(await sectionButton.evaluate(button => button === document.activeElement), true);
    await sectionButton.click();
    await page.mouse.click(2, 2);
    assert.equal(await sectionMenu.isVisible(), false);
    await page.goto(pathToFileURL(path.join(root, 'navigation-site/reports/test_human.html')).href + '#homology');
    await page.waitForFunction(() => Math.abs(document.getElementById('homology').getBoundingClientRect().top - 16) < 2);
    await page.waitForFunction(() => document.activeElement === document.querySelector('#homology h2'));
    await sectionButton.click();
    await sectionMenu.evaluate(menu => { menu.scrollTop = 0; });
    await page.screenshot({ path: path.join(root, 'section-menu-desktop.png') });
    await page.keyboard.press('Escape');
    for (const width of [320, 390]) {
      await page.setViewportSize({ width, height: 844 });
      await sectionButton.click();
      assert.equal(await sectionMenu.evaluate(menu => {
        const rect = menu.getBoundingClientRect();
        const button = document.querySelector('.report-sections-button').getBoundingClientRect();
        return rect.left >= 0 && rect.right <= innerWidth && rect.top >= 0 && rect.bottom <= button.top - 8 && menu.scrollHeight > menu.clientHeight && button.height >= 44;
      }), true);
      await sectionMenu.evaluate(menu => { menu.scrollTop = 0; });
      await page.screenshot({ path: path.join(root, `section-menu-${width}.png`) });
      await page.keyboard.press('Escape');
    }
    await page.emulateMedia({ reducedMotion: 'no-preference' });
    assert.equal(await page.evaluate(() => getComputedStyle(document.documentElement).scrollBehavior), 'smooth');
    const historyLength = await page.evaluate(() => history.length);
    await page.evaluate(() => window.scrollTo({ top: 0, behavior: 'instant' }));
    await page.waitForFunction(() => document.querySelector('#report-sections-menu [aria-current="location"]')?.hash === '#overview');
    assert.equal(await page.evaluate(() => history.length), historyLength);
    await page.setViewportSize({ width: 1440, height: 1100 });

    await visit('ptprc-site/reports/ptprc_human.html');
    await page.locator('.construct-workbench-option').filter({ hasText: 'strict_domain_2' }).click();
    const expected = JSON.parse(fs.readFileSync(path.join(__dirname, '../tests/fixtures/ptprc_mapping.json'))).construct;
    const mouse = page.locator('.construct-workbench-species-row').filter({ hasText: 'MOUSE' });
    assert.match(await mouse.innerText(), /375-470/);
    assert.equal((await mouse.locator('code').innerText()).trim(), expected.sequence);
    assert.match(await page.locator('#selectionTsv').inputValue(), new RegExp(expected.sequence));
    assert.match(await page.locator('#selectionFasta').inputValue(), new RegExp(expected.sequence.slice(0, 60)));
    await page.locator('#interactive-construct-builder').screenshot({ path: path.join(root, 'ptprc.png') });

    await visit('no-structure-site/reports/test_human.html');
    await page.locator('.construct-mutation-checkbox').check();
    assert.match(await page.locator('.construct-export-tsv').inputValue(), /MMSMMMMMMM/);

    for (const [file, mw, status] of [
      [path.join(root, 'site/calculator.html'), '50000', '#calcStatus'],
      [path.join(__dirname, '../docs/protein_concentration_calculator.html'), '50000', '#status'],
    ]) {
      await page.goto(pathToFileURL(file).href);
      await page.locator('#molecularWeight').fill(mw);
      await page.locator('#volumeUnit').selectOption('ml');
      await page.locator('#massConcentration').fill('1');
      await page.locator('#volume').fill('1');
      await page.locator('#massQuantityUnit').selectOption('mg');
      await page.locator('#massQuantity').fill('100');
      assert.equal(await page.locator(status).innerText(), 'Check inputs');
      await page.locator('#massQuantity').fill('1');
      assert.equal(await page.locator(status).innerText(), 'Calculated');
      for (const invalid of ['-1', '0']) {
        await page.locator('#massQuantity').fill(invalid);
        assert.equal(await page.locator(status).innerText(), 'Check inputs');
        assert.equal(await page.locator('#massQuantity').inputValue(), invalid);
        assert.equal(await page.locator('#factMass').innerText(), '-');
      }
      await page.locator('#massQuantity').fill('1');
      assert.equal(await page.locator(status).innerText(), 'Calculated');
    }
    await page.setViewportSize({ width: 390, height: 844 });
    for (const name of ['help', 'builder', 'constructs', 'methods', 'downloads', 'calculator', 'terms']) {
      await visit(`site/${name}.html`);
      assert.equal(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), true, `${name}: mobile page overflow`);
    }
    assert.deepEqual(errors, []);
    console.log('PASS: construct selection guide, help link and preserved builder state, section menu, anchors, focus, dismissal, browser Back, reduced motion, mobile layout, construct tabs, persistent exports, PTPRC mapping, real 3D viewer, annotation escaping, and calculators.');
  } finally {
    await browser.close();
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
