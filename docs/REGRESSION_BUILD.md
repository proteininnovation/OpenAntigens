# Regression panel build

Run this after any change to the pipeline, the portal or the build, before a full rebuild or a deploy. It builds 123 real targets through the production build and checks both portals three ways. A run takes about 30 minutes on a laptop with warm caches.

```bash
scripts/regression.sh build snapshots/openantigen-2026-11-01          # -> snapshots/regression-<date>-<sha>
scripts/regression.sh check snapshots/regression-<date>-<sha> snapshots/<previous regression build>
```

## One-time setup

The project's declared dependencies must be installed. A venv that lacks `mistune` or `jsonschema`, or that has an older Biopython, fails about 170 unit tests and can break the build. Create one from `pyproject.toml`:

```bash
python3 -m venv .venv && .venv/bin/pip install -e ".[test]" pytest
```

The browser walk-through needs Playwright 1.62.1, the version CI uses:

```bash
npm install --prefix ~/.cache/oa-playwright playwright@1.62.1
~/.cache/oa-playwright/node_modules/.bin/playwright install chromium
export PLAYWRIGHT_NODE_PATH=~/.cache/oa-playwright/node_modules
```

`regression.sh` reads `PYTHON` (default `.venv/bin/python`), `SNAPSHOT_ROOT` (default `snapshots`), `PANEL` (default `data/regression/panel.tsv`) and `PLAYWRIGHT_NODE_PATH`.

## Checklist for every change

This applies to people and coding agents alike; `AGENTS.md` makes it a rule.

1. Run the unit tests: `.venv/bin/python -m pytest -q`.
2. If the change adds a feature or fixes a case, add a row to `data/regression/panel.tsv` that would fail without it. Give it a category and any extra expectation; the tokens are listed at the top of `src/agdesign2/regression.py`.
3. Build: `scripts/regression.sh build snapshots/<latest production snapshot>`. The build must finish with `validate-release-data: ok`.
4. Check: `scripts/regression.sh check snapshots/regression-<date>-<sha> <baseline>`. Use the previous regression build as the baseline, or the production snapshot when there is none. The exit status must be 0.
5. Read `<snapshot>/regression_check.md`. Every change listed against the baseline must be one the work intended. Explain anything else, or fix it.
6. Open a few screenshots in `<snapshot>/browser_check/` for the pages the change touches, human and mouse. To browse the whole build, serve it with `python3 -m http.server --directory <snapshot>/public_site`.
7. Report the commands, the snapshot name, the three results and the change list. A failed expectation is either a regression or a wrong expectation: read the report JSON or the page before editing the panel, and say which it was.

The panel build lives under `snapshots/`, which git ignores. Keep the latest one as the next baseline and delete older ones; cloned caches share disk blocks with the reference snapshot on APFS.

## The build is the production build

`regression.sh build` runs the production warm-rebuild command with one more flag:

```bash
agdesign2 build-fresh-snapshot --target-panel data/regression/panel.tsv \
  --resume --reanalyze-reports --ortholog-jobs 6
```

`--target-panel` narrows the target TSV after the normal source-CSV and bucket selection, so a panel target that the production build would not include stops the build instead of being tested. Every later step is the production code: UniProt prefetch, the ortholog, family and paralog tables, analysis, assets, Open Targets, cross-reactivity, both portals, and the publish gate `validate-release-data`. The gate keeps all its checks; only its catalog minimums scale, from 5,000 human / 4,500 mouse pages to the panel size and the same mouse-to-human ratio.

The build first clones `cache/` and `data/` from the reference production snapshot, as a production warm rebuild does. That isolates code changes from upstream data drift. The same choice means a defect that only a fresh UniProt or AlphaFold record would expose cannot show up here.

Three outputs depend on which targets are in the build, so the panel differs from production there by design:

- The family-alignment and paralog tables are built from the target list, so family context and cross-reactivity families contain only panel proteins. Judge those in a full rebuild.
- When several human targets share one mouse ortholog, the mouse page goes to the winner among the targets in the build. Compared with a full production build, the panel can therefore move a co-ortholog's mouse page. On 2026-11-13 that moved OR4M2, SIGLEC7 and SIGLEC12.

`check-regression-panel --baseline` also accepts a production snapshot, which shows what a change would alter relative to the live build, plus the co-ortholog shifts above.

## What the checks cover

`regression.sh check` runs three checks and fails if any one does.

1. **Panel expectations** (`agdesign2 check-regression-panel`). Each target in `data/regression/panel.tsv` has a category, and each category sets default expectations (`CATEGORY_DEFAULTS` in `src/agdesign2/regression.py`): topology class, whether a design region exists, a minimum construct count, the GPCR and assembly flags. A row can add its own, such as VEGFA's `region_start=207 region_end=395`. Given a baseline panel build, the command also lists every target whose topology class, region, construct count, flags or mouse page changed. A change is expected when the change under test aims at it; anything else needs explaining before release. The report lands in `<snapshot>/regression_check.md`.
2. **Production pre-publish checks** (`scripts/release_checks.py`): the 15 HTML and JSON checks run before every deploy, on both portals.
3. **Browser walk-through** (`scripts/check_regression_site.cjs`, Playwright). It serves the site over local HTTP with all outside traffic blocked, and visits the static pages and every report page in both portals. On each report it checks for:
   - page and console errors;
   - a leaked missing value (an element reading `undefined`, `NaN`, `null` or `None`);
   - the required sections;
   - construct cards matching the index's `construct_count`;
   - every construct tab showing cards;
   - the 3D viewer reaching "viewer ready" where a structure exists;
   - the boundary editor applying;
   - every relative link resolving.

   It saves a screenshot of each page under `<snapshot>/browser_check/` for a visual pass.

## Use cases in the panel

| Category | Targets | What it guards |
|---|---|---|
| gpcr | 26 | class A, B1, C, F, adhesion, taste and olfactory receptors; 7TM advice; paralog identity (LGR4/GPR63); mouse co-orthologs (TAS2Rs) |
| single_pass | 18 | type I, II and III receptors; shed forms (CSF1, TNF); furin processing; mouse ortholog admission (IL6R, EGFR); co-ortholog labelling (SIGLECs) |
| single_pass_inferred | 3 | orientation inferred without UniProt topology (SDK1, KCNE1, ITM2C) |
| very_large | 6 | construct size warning; processed chains on NOTCH1 |
| multipass_compact | 6 | full-length membrane constructs only |
| multipass_mixed | 10 | multipass with a major extracellular region (PIEZO1 cap) |
| gpi, gpi_with_tm, gpi_untrimmed | 15 | omega-site boundaries; GPI plus TM; untrimmed warning |
| secreted, secreted_propeptide, secreted_short | 13 | signal-peptide handling; propeptides kept; a defensin below the 50-residue construct minimum gets a region but no construct; mouse winner rule (CCL2) |
| secreted_fallback | 5 | leaderless secreted fallbacks: chain, full length, propeptide |
| secreted_processed_chain, curated_override | 3 | IL16 and PDZD2 processed chains; VEGFA curated region |
| ig_constant | 1 | constant-region genes kept, germline V/D/J excluded |
| assembly, complex_partner | 8 | integrins; TCR-CD3 and BCR complexes |
| beta_stranded | 4 | MACPF, gasdermin and VDAC beta strands not read as helices |
| no_design_region | 2 | intracellular proteins: the page says no region was found |
| no_alphafold | 3 | pages without a structure |

To add a case, add a row to `panel.tsv` with the category whose defaults fit and any extra expectation, using the tokens documented at the top of `regression.py`. The target must be in the production universe.
