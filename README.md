<p align="center">
  <img src="assets/branding/openantigens-banner.svg" alt="OpenAntigens: structure-aware antigen construct design for human cell-surface and secreted proteins" width="860">
</p>

<p align="center">
  <a href="https://openantigens.org"><img alt="Portal: openantigens.org" src="https://img.shields.io/badge/portal-openantigens.org-4ec3dc?style=flat-square&labelColor=0d1f3c"></a>
  <a href="LICENSE"><img alt="Code license: Apache-2.0" src="https://img.shields.io/badge/code-Apache--2.0-007fa3?style=flat-square&labelColor=0d1f3c"></a>
  <a href="DATA_LICENSE.md"><img alt="Data license: CC BY 4.0" src="https://img.shields.io/badge/data-CC_BY_4.0-4ec3dc?style=flat-square&labelColor=0d1f3c"></a>
  <img alt="Python 3.11+" src="https://img.shields.io/badge/python-3.11%2B-82a0b9?style=flat-square&labelColor=0d1f3c">
</p>

**OpenAntigens** is built around one question in antibody, reagent, and binder work: which region of a protein should you actually express? For human cell-surface and secreted proteins, completed reports gather UniProt topology, PDB precedent, InterPro and Pfam domains, mouse and cynomolgus orthologs, paralog and family context, Open Targets disease links, and local BLAST sequence similarity for cross-reactivity review. When a compatible AlphaFold DB model is available, reports also include pLDDT, PAE, structure-aware construct boundaries, and browser-based boundary adjustment.

<p align="center">
  <b>5,328</b> human targets in the canonical portal universe
</p>

The portal is free, needs no login, and is live at **[openantigens.org](https://openantigens.org)**.

The implementation is the `agdesign2` Python package and CLI; this repository builds the portal and serves it.

## What a build contains

The portal is the primary output. Each build produces:

- a searchable target index with status, topology class, ectodomain/design-region boundaries, construct count, structure availability, cross-reactivity hit count, PubTator count, and top Open Targets disease association;
- disease-aware search and autocomplete when Open Targets data are present;
- per-target report pages with summary cards, construct tables, ortholog and paralog context, PDB evidence, cross-reactivity results, and disease-association tables;
- an integrated construct workbench linking suggested constructs, live species sequences, boundary editing, AlphaFold structure, pLDDT, PAE, cysteine and furin warnings, PTM annotations, and copy-ready TSV/FASTA export;
- downloadable index, manifest, and release metadata, plus browser TSV/FASTA copy actions for selected regions;
- support pages for construct methodology, methods, downloads, protein concentration calculation, help, terms, and privacy.

Static builds are written to `outputs/<batch>/portal/`. Dynamic serving reads the same JSON and report files live and exposes the same pages plus JSON APIs.

## Website analytics and privacy

The public portal uses Plausible Analytics to count page visits, clicks on static TSV/JSON downloads, and five browser-generated export actions. The integration uses no analytics cookies, session replay, advertising identifiers, cross-site tracking, or custom event properties. Search terms, protein sequences, residue selections, clipboard contents, names, and email addresses are not attached to analytics events.

Every portal page links to `privacy.html`, which explains what is measured and how aggregate results may be reported to funders. Dashboard access remains private. Funder reports suppress breakdowns with fewer than five unique visitors. See [`docs/ANALYTICS.md`](docs/ANALYTICS.md) for account setup, the event dictionary, testing, and quarterly reporting rules.

## Quick start

Use the local Python requested for this project:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/pip install -e . pytest notebook nbformat nbconvert
```

Prepare BLAST databases used for ortholog fallback and cross-reactivity:

```bash
.venv/bin/agdesign2 setup-blast
```

Build a small development portal:

```bash
.venv/bin/agdesign2 analyze-tsv data/inputs/accessible_human_topology_refined_accessible_only.csv \
  --output-dir outputs/surfy_batch \
  --limit 25 \
  --verbose

.venv/bin/agdesign2 build-portal outputs/surfy_batch/batch_summary.json \
  --bundle-vendor-assets
```

Open:

```text
outputs/surfy_batch/portal/index.html
```

Serve a live dynamic portal during development:

```bash
.venv/bin/agdesign2 serve-portal outputs/surfy_batch/batch_summary.json \
  --host 127.0.0.1 \
  --port 8000
```

Then visit:

```text
http://127.0.0.1:8000/
```

Build and serve the SQLite-backed portal artifact:

```bash
.venv/bin/agdesign2 build-portal-db outputs/surfy_batch/batch_summary.json \
  --db outputs/surfy_batch/openantigen.sqlite

.venv/bin/agdesign2 serve-portal-db outputs/surfy_batch/openantigen.sqlite \
  --host 127.0.0.1 \
  --port 8000
```

`openantigen.sqlite` is the canonical serving artifact for the SQLite-backed portal. Existing JSON report files remain a transition/import format while parity with the static JSON-backed portal is verified.

## Recommended build paths

### Fresh public snapshot

Use this when you want a clean, reproducible portal build from current source databases:

```bash
.venv/bin/agdesign2 build-fresh-snapshot \
  --snapshot-root snapshots \
  --snapshot-name openantigen-YYYY-MM-DD \
  --verbose
```

This creates a snapshot-local target TSV, cache, data directory, reports, static portal, `public_site/`, `snapshot_benchmark.json`, and `snapshot_manifest.json`. For a test run:
Static PyMOL construct images are excluded by default. pLDDT/PAE quality plots remain enabled and constructs with identical boundaries share one plot.

```bash
.venv/bin/agdesign2 build-fresh-snapshot \
  --snapshot-name benchmark-20 \
  --limit 20 \
  --force \
  --verbose
```

### Existing batch refresh

Use this when `outputs/surfy_batch/batch_summary.json` already exists and you need to rebuild its precomputed data:

```bash
.venv/bin/agdesign2 prepare-deployment-data \
  --target-tsv accessible_secreted_gpi_singlepass.tsv \
  --summary-path outputs/surfy_batch/batch_summary.json \
  --output-dir outputs/surfy_batch \
  --verbose
```

This runs the precompute path: UniProt cache prefetch, BLAST setup, ortholog table, family alignments, paralog references, Open Targets disease scores, batched cross-reactivity refresh, and static portal rebuild. Use `--dry-run` to inspect the planned steps and `--skip-*` flags to rerun only selected pieces.

### Resumable local refresh job

Build the active accessible target TSV:

```bash
.venv/bin/agdesign2 build-accessible-targets \
  --source-csv data/inputs/accessible_human_topology_refined_accessible_only.csv \
  --output-tsv accessible_secreted_gpi_singlepass.tsv
```

Run the restart-safe accessible portal refresh:

```bash
.venv/bin/agdesign2 refresh-accessible-portal \
  --source-csv data/inputs/accessible_human_topology_refined_accessible_only.csv \
  --output-tsv accessible_secreted_gpi_singlepass.tsv \
  --output-dir outputs/surfy_batch \
  --state-path outputs/surfy_batch/accessible_portal_refresh_state.json \
  --verbose
```

The refresh job rewrites the filtered TSV, persists state, resumes after interruption, periodically rebuilds `batch_summary.json`, and periodically rebuilds the static portal.

## Pipeline pieces

Single target analysis:

```bash
.venv/bin/agdesign2 analyze EGFR --output-dir outputs --verbose
```

Batch analysis:

```bash
.venv/bin/agdesign2 analyze-tsv accessible_secreted_gpi_singlepass.tsv \
  --output-dir outputs/surfy_batch \
  --jobs 0 \
  --verbose
```

Retry failed or incomplete rows:

```bash
.venv/bin/agdesign2 retry-batch accessible_secreted_gpi_singlepass.tsv \
  --output-dir outputs/surfy_batch \
  --summary-path outputs/surfy_batch/batch_summary.json \
  --verbose
```

Rebuild a batch summary from existing reports:

```bash
.venv/bin/agdesign2 rebuild-batch-summary accessible_secreted_gpi_singlepass.tsv \
  --output-dir outputs/surfy_batch
```

Precompute ortholog/RefSeq references:

```bash
.venv/bin/agdesign2 build-ortholog-table accessible_secreted_gpi_singlepass.tsv \
  --output-path outputs/surfy_batch/ortholog_reference_table.tsv \
  --verbose
```

Precompute InterPro-family ectodomain alignments:

```bash
.venv/bin/agdesign2 build-family-alignments accessible_secreted_gpi_singlepass.tsv \
  --output-dir outputs/surfy_batch/family_alignments \
  --verbose
```

Precompute broader paralog references:

```bash
.venv/bin/agdesign2 build-paralog-reference accessible_secreted_gpi_singlepass.tsv \
  --output-dir outputs/surfy_batch/paralogs \
  --verbose
```

Build Open Targets disease associations:

```bash
.venv/bin/agdesign2 build-open-targets outputs/surfy_batch/batch_summary.json \
  --verbose
```

Refresh report modules in place:

```bash
.venv/bin/agdesign2 refresh-report-module outputs/surfy_batch/batch_summary.json \
  --module complex-portal \
  --module cross-reactivity \
  --module family-context \
  --verbose
```

Generate or refresh static construct assets:

```bash
.venv/bin/agdesign2 generate-assets outputs/surfy_batch/batch_summary.json \
  --data-dir .agdesign2/data \
  --jobs 0 \
  --verbose
```

Build the static portal:

```bash
.venv/bin/agdesign2 build-portal outputs/surfy_batch/batch_summary.json \
  --bundle-vendor-assets \
  --verbose
```

Build the static portal from the SQLite serving artifact:

```bash
.venv/bin/agdesign2 build-portal-from-db outputs/surfy_batch/openantigen.sqlite \
  --output-dir outputs/surfy_batch/public_site \
  --bundle-vendor-assets
```

Force-refresh reports before publishing:

```bash
.venv/bin/agdesign2 build-portal outputs/surfy_batch/batch_summary.json \
  --refresh-reports \
  --bundle-vendor-assets \
  --verbose
```

Package the rendered portal for public static hosting:

```bash
python3 scripts/build_public_release.py \
  --portal-dir outputs/surfy_batch/portal \
  --output-root releases \
  --checksums
```

Release validation rejects index data with no populated PubTator links, so an offline render must use cached literature metadata rather than silently publishing `n/a` for every target.

## Key outputs

Per-target analysis:

- `outputs/<entry_name>_report.json`
- `outputs/<entry_name>_report.md`
- `outputs/<entry_name>_report_assets/`

Batch analysis:

- `outputs/surfy_batch/batch_summary.json`
- `outputs/surfy_batch/batch_summary.md`
- `outputs/surfy_batch/openantigen.sqlite`
- `outputs/surfy_batch/assets/`

Precomputed analysis data:

- `outputs/surfy_batch/ortholog_reference_table.tsv`
- `outputs/surfy_batch/family_alignments/`
- `outputs/surfy_batch/paralogs/`
- `outputs/surfy_batch/open_targets_disease_associations.tsv`
- `outputs/surfy_batch/open_targets_disease_associations.json`

Mouse portal builds copy the ortholog reference table into `mouse_data/` so reciprocal human construct sequences do not depend on incidental cache coverage.

Static portal:

- `outputs/surfy_batch/portal/index.html`
- `outputs/surfy_batch/portal/reports/*.html`
- `outputs/surfy_batch/portal/downloads/`
- `outputs/surfy_batch/portal/portal_metadata.json`

Public release bundle:

- `releases/openantigen-YYYY-MM-DD/index.html`
- `releases/openantigen-YYYY-MM-DD/release_manifest.json`
- `releases/openantigen-YYYY-MM-DD/README_RELEASE.md`

Runtime caches and generated artifacts live under `outputs/`, `.agdesign2/`, `snapshots/`, and `releases/`; they are intentionally ignored by Git.

## Documentation

Detailed documentation lives in [docs/DETAILED_GUIDE.md](docs/DETAILED_GUIDE.md). It covers:

- end-to-end portal builds
- individual pipeline commands, inputs, and outputs
- report module refreshes
- static versus dynamic portal serving
- fresh snapshots and public release packaging
- construct boundary logic
- ortholog, paralog, family, disease, BLAST, PDB, and AlphaFold interpretation
- portal pages and downloadable data

Repository organization is summarized in [docs/REPOSITORY_STRUCTURE.md](docs/REPOSITORY_STRUCTURE.md).

## Project layout

- `src/agdesign2/`: package source for data retrieval, analysis, module refreshes, and portal rendering/serving
- `data/inputs/`: canonical and legacy target input tables
- `assets/branding/`: source branding assets copied into portal builds
- `scripts/`: release packaging, figure generation, and measurement utilities
- `tests/`: unit and regression tests
- `docs/`: detailed guides and small support pages

## Tests

Run:

```bash
.venv/bin/python -m unittest discover -s tests -v
```

or:

```bash
.venv/bin/pytest
```

## Notes

- Fresh human snapshots query cached Complex Portal human ComplexTab data for every target. Integrin assembly warnings remain limited to canonical ITGA and ITGB1-ITGB8 chains.
- `build-portal` refreshes stale per-gene reports by default when they are older than precomputed ortholog or family files. Use `--no-refresh-stale-reports` for a pure HTML rebuild.
- BLAST sequence similarity is a ranking aid for cross-reactivity review, not a hard pass/fail classifier.
- Cynomolgus monkey support uses taxon `9541` and favors canonical RefSeq protein records.

## Citation

A database paper describing OpenAntigens is in preparation for the *Nucleic Acids Research* Database Issue. Until it appears, cite the portal and the snapshot build date, following the attribution in [DATA_LICENSE.md](DATA_LICENSE.md):

> OpenAntigens, Institute for Protein Innovation, created by Andre A. R. Teixeira. https://openantigens.org. Include the portal build date and software version.

## License

The OpenAntigens software is released under the [Apache License 2.0](LICENSE). OpenAntigens-generated annotations are released under CC BY 4.0. Third-party source data keep their original licenses, terms of use, and citation requirements. See [DATA_LICENSE.md](DATA_LICENSE.md) for the full terms and recommended attribution.

Report security issues through the process in [SECURITY.md](SECURITY.md).
