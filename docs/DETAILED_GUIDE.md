# OpenAntigens Detailed Guide

## Overview

This document explains three things:

1. how to run each piece of the OpenAntigens pipeline
2. how to build, serve, and package the portal
3. how each major result in a per-gene report is generated and how to interpret it

All commands below assume you are running from the repository root.

## Portal Build Workflow

The normal deliverable is a portal, not an isolated report. A complete portal build has these stages:

1. create the Python environment
2. build or choose a target TSV
3. prepare local caches and BLAST databases
4. precompute ortholog, family, paralog, and disease context
5. analyze target reports
6. generate optional static assets
7. refresh selected report modules when needed
8. build the static portal or serve the dynamic portal
9. package the rendered portal for public hosting

The CLI can run these stages one by one, or through combined commands such as `prepare-deployment-data` and `build-fresh-snapshot`.

## Pipeline Command Reference

### 1. Create the environment

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/pip install -e . pytest notebook nbformat nbconvert
```

### 2. Build the active portal target TSV

For the current cell-surface/secreted portal target universe:

```bash
.venv/bin/agdesign2 build-accessible-targets \
  --source-csv data/inputs/accessible_human_topology_refined_accessible_only.csv \
  --output-tsv accessible_secreted_gpi_singlepass.tsv
```

This filters the refined accessible-topology CSV to the default `Secreted`, `GPI`, `Single-pass`, and `Multipass` buckets. The output TSV is the common input for precomputes and batch analysis.

The canonical CSV can also be passed directly to batch and precompute commands. Use `--limit` for smaller test builds.

### 3. Prepare BLAST databases

```bash
.venv/bin/agdesign2 setup-blast
```

This builds:

- `.agdesign2/data/blastdb`
  - reviewed UniProt FASTA sets used for cross-reactivity ranking
- `.agdesign2/data/ortholog_blastdb`
  - broader organism proteomes used as ortholog fallback search space

### 4. Prefetch UniProt entries

Target-list prefetching is used by `build-fresh-snapshot` internally. For a general local cache of reviewed taxa:

```bash
.venv/bin/agdesign2 prefetch-uniprot \
  --taxa 9606 10090 \
  --cache-dir .agdesign2/cache \
  --page-size 500
```

This reduces repeated UniProt traffic during batch runs. It is optional for small development batches, but useful for large portal refreshes.

### 5. Precompute the ortholog reference table

```bash
.venv/bin/agdesign2 build-ortholog-table accessible_secreted_gpi_singlepass.tsv \
  --output-path outputs/surfy_batch/ortholog_reference_table.tsv \
  --verbose
```

This writes:

- `outputs/surfy_batch/ortholog_reference_table.tsv`
- `outputs/surfy_batch/ortholog_reference_table.json`

### 6. Precompute family alignments and identity matrices

```bash
.venv/bin/agdesign2 build-family-alignments accessible_secreted_gpi_singlepass.tsv \
  --output-dir outputs/surfy_batch/family_alignments \
  --verbose
```

This writes:

- `outputs/surfy_batch/family_alignments/family_alignment_index.tsv`
- `outputs/surfy_batch/family_alignments/family_alignment_index.json`
- `outputs/surfy_batch/family_alignments/families/*.json`
- `outputs/surfy_batch/family_alignments/families/*.fasta`

### 7. Precompute broader paralog reference sets

```bash
.venv/bin/agdesign2 build-paralog-reference accessible_secreted_gpi_singlepass.tsv \
  --output-dir outputs/surfy_batch/paralogs \
  --verbose
```

This writes:

- `outputs/surfy_batch/paralogs/paralog_index.tsv`
- `outputs/surfy_batch/paralogs/paralog_index.json`
- `outputs/surfy_batch/paralogs/targets/*.json`

These per-target records combine:

- Ensembl Compara human paralog calls
- precomputed InterPro-family membership when available
- HGNC family membership as a supporting source

### 8. Run per-gene analyses

Single gene:

```bash
.venv/bin/agdesign2 analyze EGFR --output-dir outputs --verbose
```

Whole portal TSV:

```bash
.venv/bin/agdesign2 analyze-tsv accessible_secreted_gpi_singlepass.tsv \
  --output-dir outputs/surfy_batch \
  --jobs 0 \
  --verbose
```

Development subset:

```bash
.venv/bin/agdesign2 analyze-tsv accessible_secreted_gpi_singlepass.tsv \
  --output-dir outputs/surfy_batch \
  --limit 25 \
  --verbose
```

Retry errored or incomplete rows:

```bash
.venv/bin/agdesign2 retry-batch accessible_secreted_gpi_singlepass.tsv \
  --output-dir outputs/surfy_batch \
  --summary-path outputs/surfy_batch/batch_summary.json \
  --verbose
```

Rebuild `batch_summary.json` from an existing TSV plus report files:

```bash
.venv/bin/agdesign2 rebuild-batch-summary accessible_secreted_gpi_singlepass.tsv \
  --output-dir outputs/surfy_batch
```

Per-target analysis writes report JSON, report Markdown, and report assets. Batch analysis maintains `batch_summary.json` and `batch_summary.md`.

### 9. Generate or refresh static report assets

```bash
.venv/bin/agdesign2 generate-assets outputs/surfy_batch/batch_summary.json \
  --data-dir .agdesign2/data \
  --jobs 0 \
  --verbose
```

This can render construct structure images and quality plots after reports exist. Use `--no-structure-images` or `--no-quality-plots` to skip either class of asset.

### 10. Precompute Open Targets disease associations

```bash
.venv/bin/agdesign2 build-open-targets outputs/surfy_batch/batch_summary.json --verbose
```

This writes:

- `outputs/surfy_batch/open_targets_disease_associations.tsv`
- `outputs/surfy_batch/open_targets_disease_associations.json`

The default is the scalable bulk-data path: OpenAntigens resolves the latest Open Targets bulk release at runtime, downloads `target`, `disease`, `association_overall_direct`, and `association_by_datasource_indirect`, then joins them locally by Ensembl target ID and disease ID. The resolved release is recorded in `open_targets_disease_associations.json`. The portal index uses the broader indirect disease score for the top disease column and disease search, while target pages show both indirect and direct overall scores. Use `--release VERSION` to pin a specific Open Targets release for reproducible rebuilds, or `--source api` only for small debugging runs because it resolves targets and pages disease associations through GraphQL one target at a time.

### 11. Refresh report modules in place

Use module refreshes for targeted updates without rerunning every report stage:

```bash
.venv/bin/agdesign2 refresh-report-module outputs/surfy_batch/batch_summary.json \
  --module complex-portal \
  --module cross-reactivity \
  --module family-context \
  --verbose
```

Limit to specific targets with repeated `--query` values:

```bash
.venv/bin/agdesign2 refresh-report-module outputs/surfy_batch/batch_summary.json \
  --module cross-reactivity \
  --query UNC5D_HUMAN \
  --verbose
```

Supported modules:

- `complex-portal`: records a cached human Complex Portal lookup for every report. Empty successful lookups are marked `no_hits`; errors remain explicit.
- `cross-reactivity`
- `family-context`
- `target-metadata`

Whole-batch cross-reactivity refreshes use a batched local BLAST path when the standard `BlastClient` is active.

### 12. Build the static portal

```bash
.venv/bin/agdesign2 build-portal outputs/surfy_batch/batch_summary.json \
  --bundle-vendor-assets \
  --verbose
```

Default behavior:

- stale reports older than the precomputed ortholog or family files are refreshed before publishing
- error rows are retried before publishing when possible
- browser libraries are copied/downloaded into the portal directory when `--bundle-vendor-assets` is used

Force-refresh every report:

```bash
.venv/bin/agdesign2 build-portal outputs/surfy_batch/batch_summary.json \
  --refresh-reports \
  --bundle-vendor-assets \
  --verbose
```

Disable stale-report refresh for a pure HTML rebuild:

```bash
.venv/bin/agdesign2 build-portal outputs/surfy_batch/batch_summary.json \
  --no-refresh-stale-reports
```

Main portal page:

- `outputs/surfy_batch/portal/index.html`

Per-gene HTML pages:

- `outputs/surfy_batch/portal/reports/*.html`

Portal support pages:

- `outputs/surfy_batch/portal/help.html`
- `outputs/surfy_batch/portal/builder.html`
- `outputs/surfy_batch/portal/constructs.html`
- `outputs/surfy_batch/portal/methods.html`
- `outputs/surfy_batch/portal/downloads.html`
- `outputs/surfy_batch/portal/calculator.html`
- `outputs/surfy_batch/portal/terms.html`
- `outputs/surfy_batch/portal/privacy.html`

### 13. Serve a dynamic portal

```bash
.venv/bin/agdesign2 serve-portal outputs/surfy_batch/batch_summary.json \
  --host 127.0.0.1 \
  --port 8000
```

This mode reads the current batch summary, per-gene JSON reports, report asset folders, and local AlphaFold files directly at request time.

The dynamic server keeps the same overall page layout and report rendering logic as the static portal, but serves:

- `/index.html`
- `/reports/<entry>.html`
- `/api/genes`
- `/api/gene/<entry_name>`

### 14. Package a public static release

After the static portal has been built and checked locally, package only the publishable website files:

```bash
python3 scripts/build_public_release.py \
  --portal-dir outputs/surfy_batch/portal \
  --output-root releases \
  --checksums
```

The default release name is `openantigen-YYYY-MM-DD`. The release directory is
ready to serve directly with nginx, Apache, institutional hosting, or object
storage plus a CDN:

```text
releases/openantigen-YYYY-MM-DD/index.html
```

The release package includes:

- rendered portal pages
- per-target report HTML files
- copied portal structures/assets/download files
- `portal_metadata.json` with the visible release snapshot date, software version, target counts, status counts, and topology counts
- `release_manifest.json`
- `README_RELEASE.md`
- `LICENSE` and `DATA_LICENSE.md`, when present

The footer of every portal page displays the OpenAntigens software version and
release snapshot date so screenshots, citations, and downloaded files can be
traced back to a specific build.

The release package intentionally excludes:

- `.agdesign2` caches and BLAST databases
- refresh state and logs
- intermediate JSON/Markdown files outside the rendered portal
- hidden desktop metadata such as `.DS_Store`

Use `--checksums` when creating a formal archived release. This computes a
SHA-256 digest for every published file and records it in
`release_manifest.json`; it is slower for multi-GB releases.

### 15. Optional helpers

Rebuild only the Markdown batch dashboard:

```bash
.venv/bin/agdesign2 summarize-batch outputs/surfy_batch/batch_summary.json
```

Prefetch reviewed UniProt entries into cache:

```bash
.venv/bin/agdesign2 prefetch-uniprot --taxa 9606 10090 --cache-dir .agdesign2/cache --page-size 500
```

## What Gets Written Per Gene

For a target like `EGFR_HUMAN`, the workflow writes:

- `outputs/egfr_human_report.json`
- `outputs/egfr_human_report.md`
- `outputs/egfr_human_report_assets/`

The JSON is the full machine-readable record. The Markdown is the human-readable report. The assets directory contains the rendered structure images and the pLDDT/PAE inspection plots.

## Input Resolution

The main entrypoint accepts:

- UniProt accession, for example `P00533`
- UniProt entry name, for example `EGFR_HUMAN`
- official gene symbol, for example `EGFR`

Resolution is done through UniProt. The pipeline prefers the reviewed human record for gene-symbol input. If the result is ambiguous, it raises an explicit error rather than silently picking a questionable target.

Important distinction:

- `analyze` resolves live from the query you give it
- `analyze-tsv` prefers the `gene` column, then `uniprot_name`, then `uniprot_id`
- `build-ortholog-table` and `build-family-alignments` use `uniprot_name` first, then `gene`, then `uniprot_id`

That difference is intentional because the precompute jobs are meant to track exact entries from the membrane reference file.

## Topology And The Design Region

The design region is derived from UniProt features.

### Single-pass membrane proteins

The pipeline looks for:

- `TRANSMEM`
- `TOPO_DOM`
- `SIGNAL`

Rules:

- exactly one transmembrane region is expected for the single-pass path
- if the extracellular side is clearly on the N-terminus, the ectodomain is `signal-end+1` to `tm-start-1`
- if the extracellular side is clearly on the C-terminus, the ectodomain is `tm-end+1` to sequence end
- if topology is incomplete but a signal peptide exists, the code can still infer an N-terminal ectodomain with a warning
- if topology is contradictory or missing in a way that prevents confident assignment, the ectodomain is left unresolved and downstream structure-led design is skipped

### Secreted proteins

If there is no transmembrane region, the workflow checks whether the mature secreted chain can be defined.

Rules:

- if a signal peptide is present, the design region starts just after it
- if topology annotations indicate extracellular context at the mature-chain start, those are used as supporting evidence
- if the curated target universe classifies a no-transmembrane target as secreted but UniProt lacks signal/topology boundaries, the workflow uses UniProt chain or peptide boundaries when available, then propeptide-trimmed sequence when available, and otherwise the full canonical sequence with a warning note

## Precomputed References Versus Live Lookups

Per-gene analysis now prefers local precomputed reference files when they exist:

- ortholog reference table:
  - `outputs/surfy_batch/ortholog_reference_table.tsv`
- paralog reference index:
  - `outputs/surfy_batch/paralogs/paralog_index.tsv`
- family alignment index:
  - `outputs/surfy_batch/family_alignments/family_alignment_index.tsv`

The logic is:

1. resolve the target live in UniProt
2. look for a matching precomputed ortholog row
3. use that row for canonical family and ortholog data when available
4. look for a matching precomputed paralog record
5. if none exists, look for the corresponding precomputed family alignment set
6. use live services only as fallback when the target is missing from the precomputed records

This reduces repeated HGNC/InterPro/NCBI traffic and makes batch reporting more stable and reproducible.

### Deployment precompute command

For public deployment, prefer the combined bulk/precompute entry point instead
of manually running each long step one at a time:

```bash
.venv/bin/agdesign2 prepare-deployment-data \
  --target-tsv accessible_secreted_gpi_singlepass.tsv \
  --summary-path outputs/surfy_batch/batch_summary.json \
  --output-dir outputs/surfy_batch \
  --verbose
```

The command performs the deployment-oriented steps in this order:

1. bulk-prefetch reviewed UniProt entries for human and mouse into the local HTTP cache
2. prepare local BLAST databases
3. build the ortholog/RefSeq reference table
4. build InterPro-family alignment matrices
5. build broader paralog reference matrices
6. build Open Targets disease-association tables
7. refresh cross-reactivity with batched local BLAST
8. rebuild the static portal

Use `--dry-run` to print the exact planned steps without starting network,
BLAST, or report-writing work:

```bash
.venv/bin/agdesign2 prepare-deployment-data --dry-run
```

Use `--skip-*` flags to rerun only selected parts. For example, after changing
cross-reactivity parameters:

```bash
.venv/bin/agdesign2 prepare-deployment-data \
  --skip-uniprot \
  --skip-blast-setup \
  --skip-orthologs \
  --skip-family-alignments \
  --skip-paralogs \
  --skip-open-targets \
  --summary-path outputs/surfy_batch/batch_summary.json \
  --output-dir outputs/surfy_batch \
  --verbose
```

That keeps only the batched cross-reactivity refresh and portal rebuild.

### Fresh snapshot from current source databases

For a public release, use the fresh-snapshot command when the goal is to avoid
reusing local caches and regenerate the website from current source databases:

```bash
.venv/bin/agdesign2 build-fresh-snapshot \
  --snapshot-root snapshots \
  --snapshot-name openantigen-YYYY-MM-DD \
  --verbose
```

Fresh public snapshots exclude static PyMOL construct images by default. Use
`--structure-images` only for a private build that explicitly needs them.

Fresh-snapshot behavior:

1. use the local surface/secreted protein universe CSV as the only local biological input
2. create `snapshots/<name>/cache` and `snapshots/<name>/data`
3. build a target TSV inside the snapshot
4. bulk-prefetch current UniProt records for those targets
5. download/build fresh BLAST databases
6. build ortholog, family, and paralog precomputed references
7. analyze every target into snapshot-local report JSON/Markdown/assets
8. build Open Targets disease tables
9. refresh cross-reactivity with batched local BLAST
10. build the portal and copy a publishable `public_site/`
11. write `snapshot_benchmark.json` and `snapshot_manifest.json`

For benchmarking or development:

```bash
.venv/bin/agdesign2 build-fresh-snapshot \
  --snapshot-name benchmark-20 \
  --limit 20 \
  --force \
  --verbose
```

By default, UniProt prefetching uses target-list bulk queries. This avoids
mirroring the entire reviewed human and mouse proteomes when building a small
subset or when only the portal target universe is needed. Use
`--uniprot-prefetch-mode taxa` only when a full reviewed-proteome UniProt cache
is deliberately required.

## InterPro And Pfam Annotations

The workflow pulls both InterPro and Pfam annotations from the InterPro API.

### What is kept for design

For construct design, the pipeline filters annotations aggressively:

- only annotations fully contained inside the ectodomain are kept
- only these types are considered:
  - `domain`
  - `repeat`
  - `homologous_superfamily`
- redundant overlaps are pruned
- when InterPro and Pfam describe the same span, InterPro is preferred
- Pfam entries with the same integrated InterPro assignment are treated as secondary evidence rather than separate design objects

This is why the design-oriented InterPro/Pfam section should be much cleaner than the raw upstream annotation stream.

## Canonical Family Assignment

Each protein may receive one canonical family assignment.

### Where it comes from

The family is selected from annotations of type `family`, with strong preference for InterPro families.

### How it is ranked

Candidate families are scored using:

- source preference
  - InterPro is favored over other databases
- coverage of the full protein sequence
  - merged family fragments covering more of the protein score higher
- fragment count
  - multi-fragment families get a small bonus
- accession preference
  - `IPR*` gets a small boost
- name similarity to the target protein name and gene symbol
  - very broad full-length labels are down-ranked when they look like """this exact protein""" rather than a meaningful family

The selected family is stored in the ortholog reference table as:

- `canonical_family_accession`
- `canonical_family_name`

Example:

- `PTPRC` should resolve to `Protein-Tyrosine Phosphatase (IPR050348)`

## Precomputed Family Alignment Sets

These are generated by:

```bash
.venv/bin/agdesign2 build-family-alignments data/inputs/accessible_human_topology_refined_accessible_only.csv --output-dir outputs/surfy_batch/family_alignments --verbose
```

### What goes into them

For each input row:

1. the target is resolved in UniProt
2. the ectodomain is derived from topology
3. the canonical family is selected
4. the row is kept only if that canonical family is an InterPro family with an accession starting with `IPR`

Pfam-only families are intentionally ignored here.

### Sequence used

The stored sequence is the human ectodomain sequence only, not the full-length protein.

### Identity matrix meaning

The family identity matrix is directional.

Cell `X -> Y` means:

- number of exact matches in the global alignment between `X` and `Y`
- divided by the length of `X`
- multiplied by 100

That means the upper-right and lower-left cells are intentionally different when one ectodomain is shorter or largely contained within another.

This solves the earlier symmetric-matrix problem where a short protein could appear to have """100% identity""" to multiple longer relatives simply because it aligned perfectly over its own full length.

When a precomputed paralog record is available, the report also carries a directional coverage matrix.

Cell `X -> Y` in that matrix means:

- how much of ectodomain `X` is aligned to ectodomain `Y`
- reported as a percent of the full length of `X`

This lets you distinguish:

- high identity over nearly the whole ectodomain
- high identity over only a smaller contained region

### How broader paralog sets are assembled

The precomputed paralog builder combines three evidence streams:

1. Ensembl Compara human paralog calls
2. the precomputed InterPro family-alignment set for the canonical family, when one exists
3. HGNC family membership as supporting fallback evidence

Candidates are deduplicated, resolved to reviewed human UniProt entries when possible, reduced to ectodomain sequences, and then ranked by:

- evidence source strength
- target-to-candidate directional ectodomain identity
- target-to-candidate directional ectodomain coverage

The stored paralog context is intentionally capped so the report remains readable.

### What is written per family

Each family JSON includes:

- family accession and name
- member list
- directional identity matrix
- directional coverage matrix when available
- pairwise alignments

Each family FASTA includes the human ectodomain sequences of the members.

### When family context is omitted from reports

If the family is too large, the report suppresses family-context display. The current limit is 20 members.

## Ortholog Reference Table

The ortholog reference table is generated by:

```bash
.venv/bin/agdesign2 build-ortholog-table data/inputs/accessible_human_topology_refined_accessible_only.csv --output-path outputs/surfy_batch/ortholog_reference_table.tsv --verbose
```

### One-row-per-input semantics

This table is meant to be audit-friendly.

Rules:

- one output row per input row
- input order is preserved
- duplicates are preserved
- failures are written as `status=error` rows rather than being dropped
- TSV and JSON are checkpointed after every row
- resume is keyed by `input_index`

### Main columns

Every row preserves the input identifiers:

- `input_index`
- `input_uniprot_id`
- `input_uniprot_name`
- `input_gene`
- `input_prot_family`
- `query`
- `source_column`
- `status`
- `error`

Biology columns include:

- `human_gene_symbol`
- `mouse_gene_symbol`
- `macaca_fascicularis_gene_symbol`
- `canonical_family_accession`
- `canonical_family_name`
- `human_refseq_accession`
- `mouse_refseq_accession`
- `macaca_fascicularis_refseq_accession`
- `human_uniprot_accession`
- `human_refseq_sequence`
- `mouse_refseq_sequence`
- `macaca_fascicularis_refseq_sequence`
- `mouse_identity_to_human`
- `macaca_fascicularis_identity_to_human`
- species-specific notes columns

### How ortholog symbols are obtained

Human, mouse, and cynomolgus monkey ortholog relationships are seeded from HGNC HCOP.

### How RefSeq proteins are selected

For each species, the code queries NCBI RefSeq protein records and ranks candidates.

Preference order:

- `MANE Select` when available
- `RefSeq Select` when available
- curated `NP_` accessions ahead of predicted `XP_`
- lower isoform number ahead of higher isoform number

For cynomolgus monkey, the workflow still uses RefSeq, but there is often no `MANE Select` equivalent, so it picks the best-ranked available canonical-like protein candidate.

### Pairwise identity in the table

The table-level species identity values are full-sequence global-alignment identities, not ectodomain-only identities.

The per-report ortholog section can then transfer the human ectodomain onto the species sequence and calculate ectodomain-focused values.

## Orthologs In Per-Gene Reports

### Preferred source

If the target exists in the precomputed ortholog table, the report uses that first.

### How ectodomain boundaries are transferred

When precomputed ortholog data is available:

1. align the full human target sequence to the precomputed species sequence
2. transfer the human ectodomain boundaries through the full-sequence alignment
3. define the species ectodomain on that transferred span
4. align human ectodomain to species ectodomain to calculate identity and coverage

### Construct-level ortholog mapping

For each construct:

1. calculate its position relative to the human ectodomain
2. project that span into the species ectodomain alignment
3. extract the species construct sequence
4. compute identity of that construct to the human construct sequence

That is why each construct can show:

- human boundary and sequence
- mouse equivalent boundary and sequence
- cynomolgus monkey equivalent boundary and sequence
- identity to the human construct

### Fallback behavior

If a target is not in the precomputed ortholog table, the older live path is used:

- same-name UniProt entry where possible
- gene-symbol search in the correct species
- proteome BLAST fallback if needed

## AlphaFold Data And Structural Region Calling

If an ectodomain is available, the workflow tries to fetch AlphaFold artifacts:

- PDB
- PAE JSON

### pLDDT structured-region calling

The workflow first identifies residues above a pLDDT threshold and merges nearby gaps.

There are two passes:

- lenient pass
  - pLDDT threshold: `60`
  - allowed disordered gap: `12`
- strict seed pass
  - pLDDT threshold: `70`
  - allowed disordered gap: `8`

Other defaults:

- minimum structured segment: `25 aa`

Interpretation:

- lenient regions are intended to preserve larger folded extracellular units
- strict seeds are intended to isolate more domain-like cores

### PAE-based domain splitting

Strict seed regions are then recursively split using the PAE matrix.

For each candidate cut, the code computes:

- `inter_block_pae`
  - mean PAE across the left-versus-right blocks
- `intra_block_pae`
  - average of the mean intra-block PAEs of the left and right blocks
- `separation`
  - `inter_block_pae - intra_block_pae`
- `boundary_pae`
  - local inter-block PAE around the proposed boundary
- `linker_plddt`
  - mean pLDDT in a local boundary window

The split is considered only if:

- `inter_block_pae >= pae_domain_threshold`
- `separation >= pae_separation_threshold`

Default thresholds:

- `pae_domain_threshold = 12`
- `pae_separation_threshold = 4`
- `domain_linker_window = 5`
- `min_domain_size = 40`
- `max_domain_recursion_depth = 4`

The candidate split score is:

- `separation`
- plus `0.35 * max(0, boundary_pae - intra_block_pae)`
- plus a linker bonus when the local linker pLDDT is low

In plain language:

- a good split has low internal error within each side
- high error between the two sides
- and ideally a softer, less confident linker around the boundary

### Structured-region confidence

Domain confidence is derived from mean intra-domain PAE.

Current transform:

- higher mean intra-domain PAE gives lower confidence
- clipped into the range `0.4` to `0.95`

This confidence is reported for structural regions when available.

## Construct Generation

Constructs are generated as boundary proposals from multiple evidence streams.

### Seed sources

The pipeline creates construct seeds from:

- the full ectodomain
- filtered InterPro/Pfam annotations
- UniProt `DOMAIN` and `REGION` annotations inside the ectodomain
- AlphaFold lenient regions
- AlphaFold strict domains
- ectodomain-filtered PDB construct chains

### Base scores

Current seed scores are:

- full ectodomain: `90`
- PDB construct chain: `85`
- InterPro representative domain: `82`
- Pfam domain seed: `80`
- UniProt domain/region: `80`
- AlphaFold strict domain: `76`
- AlphaFold generic structural region: `75`
- AlphaFold lenient region: `74`

The full ectodomain is penalized by `5` if the ectodomain-level cysteine analysis already contains surface-exposed unpaired-cysteine warnings.

### Construct filtering

After seeds are created:

- duplicate boundaries are deduplicated by exact `(start, end)` span, keeping the highest-scoring version
- all candidates are passed through annotation-aware classification
- candidates outside the ectodomain are removed
- constructs shorter than the configured minimum are removed

Current minimum construct length:

- `50 aa`

### Construct classification

Each surviving construct is classified against the filtered ectodomain design annotations.

Possible classes include:

- `complete_domain`
- `multi_domain_unit`
- `repeat_module`
- `family_module`
- `partial_domain`
- `structured_region`

Interpretation:

- `complete_domain`
  - contains exactly one full annotated domain
- `multi_domain_unit`
  - contains multiple full annotated domains
- `repeat_module`
  - captures a repeat-only unit
- `family_module`
  - maps to a broader family/superfamily region without a cleaner domain call
- `partial_domain`
  - clips at least one annotation boundary
- `structured_region`
  - structure-led region without a complete annotation match

### Score adjustments after classification

After classification:

- `complete_domain`: `+5`
- `multi_domain_unit` and `repeat_module`: `+3`
- `partial_domain`: `-8`

Partial-domain constructs also get an explicit warning telling you which annotation was clipped.

### What the construct score means

The construct score is a ranking heuristic, not a probability and not an experimental success predictor.

It is best read as:

- higher score = stronger support from curated annotation or experimental precedent
- lower score = more speculative or more structurally derived
- penalties usually mean boundary clipping or other design risk

Use it to prioritize review, not as an absolute pass/fail threshold.

## Construct Detail Metrics

The detailed construct sections compute per-construct metrics.

### Sequence and boundaries

Each construct stores:

- `start`
- `end`
- `length`
- exact amino-acid sequence

### Structural metrics

When AlphaFold data is available:

- `mean_plddt`
- `min_plddt`
- `max_plddt`
- `structured_fraction`
  - fraction of residues with pLDDT above the structured threshold
- `mean_intra_pae`
- `max_intra_pae`

Interpretation:

- higher `mean_plddt` is better
- lower `mean_intra_pae` is better
- a high `structured_fraction` means most of the construct looks confidently folded

### PAE split diagnostics

Each construct can also show candidate internal split diagnostics:

- `score`
- `inter_block_pae`
- `intra_block_pae`
- `separation`
- `boundary_pae`
- `linker_plddt`

These are not """construct quality""" scores by themselves. They are diagnostics showing whether a construct contains an internal PAE-supported boundary that might justify splitting it further.

## PDB-Based Construct Evidence

PDB construct evidence is extracted from UniProt cross-references and chain mappings.

Rules:

- only chain spans fully contained in the ectodomain are retained
- PDB entries fully outside the ectodomain are ignored
- if multiple PDB entries use the exact same boundaries, the report shows only one representative example

These are treated as evidence-backed construct suggestions and are also evaluated with the same construct-level metrics used for the other construct types.

## Cysteine Analysis

### How cysteines are paired

Within a construct, cysteines are called from the AlphaFold structure and paired if the sulfur-sulfur distance is `<= 2.4 A`.

### Surface exposure

For an unpaired cysteine, the pipeline estimates whether it is surface exposed using a simple local-neighbor count around the residue.

### Warnings

Warnings are added when:

- a construct contains unpaired cysteines
- an unpaired cysteine also appears surface exposed

These warnings are construct-specific as well as ectodomain-level.

## Candidate furin-like motifs

Furin motifs are scanned only within the ectodomain or secreted design region.

Current motif:

- `R.[KR]R`

For each match, the report includes:

- absolute residue boundaries
- matched motif
- simple mutation suggestions at the first and last motif residues

These are annotations only. The pipeline does not redesign the sequence automatically.

## Ligand Interaction Regions

The report includes UniProt-derived ligand-interaction context when features overlap the design region.

The workflow accepts:

- `BINDING`
- `BINDING_SITE`
- feature records whose text or metadata contains ligand/binding language

These are shown both at the protein level and within any construct that overlaps the annotated interaction span.

Interpretation:

- this section helps flag constructs that may include or exclude ligand-contacting regions
- it is descriptive, not a functional binding prediction

## Assembly / Partner Requirements

This section is intended to separate true partner requirements from generic interaction noise.

### Sources used

- Complex Portal human complex records
- Complex Portal mouse ortholog complex records
- UniProt `SUBUNIT` comments as interaction context

### Current classifications

- `obligatory_partner_requirement`
- `conditional_assembly`
- `interaction_context`
- `heteromeric_assembly`
- `stable_complex_context`
- `mouse_supported_complex_context`

### UniProt text interpretation

UniProt sentences are scanned for assembly language.

UniProt `SUBUNIT` comments can support integrin assembly classification, but
they do not independently create obligatory-partner warnings. Non-obligatory
free-text interaction rows are not shown in the public portal.

### Family rules

The current explicit family elevation is for integrins, but it is applied through human Complex Portal evidence:

- `ITGA*` chains are flagged as requiring an integrin beta partner
- `ITGB*` chains are flagged as requiring an integrin alpha partner

This is used when Complex Portal reports a curated human alpha/beta integrin complex.

### Complex Portal

The public human portal records the Complex Portal lookup status for every
target and shows returned records under **Interactions / Assembly**. Human
Complex Portal records are the source of truth for assembly warnings. Mouse
Complex Portal records are fetched for canonical integrin orthologs as
supporting context; mouse-only evidence does not trigger a human obligatory
warning.

Enable it with:

```bash
.venv/bin/agdesign2 analyze ITGAV --output-dir outputs --complex-portal --verbose
```

When used, it can add:

- curated stable complex membership
- heterodimer evidence
- partner symbols
- mouse ortholog complex support

Human Complex Portal membership that does not establish an obligatory partner is classified as stable complex context.

## Cross-Reactivity BLAST

The cross-reactivity section uses local `blastp`.

### Databases

By default, the workflow prepares reviewed UniProt sets for:

- human
- mouse
- cynomolgus monkey

These are used as the main cross-reactivity search space.

### Query sequence

The query is the target ectodomain or secreted design-region sequence, not the full-length protein.

### Self-hit filtering

Hits are excluded if the BLAST subject identifier contains:

- the target accession
- the target entry name

### Coverage definition

Coverage is:

- aligned query span
- divided by query length
- multiplied by 100
- capped at `100`

This is why coverage should never exceed 100%.

### Hit ranking

Hits are deduplicated by subject and ranked primarily by:

- higher bitscore
- then lower e-value

### Macaque canonicalization

For cynomolgus monkey, the workflow can remap raw hit descriptions to canonical RefSeq proteins by gene symbol so the report is more interpretable than a raw `XP_*` list.

### Interpretation

This section is a ranking aid. It does not apply an automatic pass/fail threshold for """acceptable""" cross-reactivity risk.

## Portal Contents

The portal is built from `batch_summary.json` plus the per-gene report JSON files.

Main page features include:

- processed-gene table
- status
- topology class and design-region boundaries
- ectodomain boundaries
- AlphaFold presence
- PDB presence
- construct count
- cross-reactivity hit count
- PubTator literature count
- Open Targets top disease association
- disease-aware search and autocomplete when Open Targets associations are available
- sorting and filtering
- links to help, construct methodology, methods, downloads, calculator, terms, and the dedicated Interactive Construct Builder support page

Per-gene pages include:

- summary cards
- structure-aware interactive viewer
- interactive sequence selection
- pLDDT and PAE panels
- live selected-region export in TSV and FASTA
- construct cards with copy-ready sequence exports
- PTM, unpaired-cysteine, and furin-site review cues
- family context
- PDB evidence
- ortholog construct equivalents
- cross-reactivity section
- Open Targets disease-association table

Download pages include:

- portal index TSV/JSON
- browser TSV/FASTA copy actions for selected regions
- release/build metadata
- data manifest files that support reproducible downstream analysis

### License and terms

The portal includes `terms.html`, which separates:

- OpenAntigens software: Apache License 2.0
- OpenAntigens-generated annotations: CC BY 4.0
- third-party source data: original source terms apply

OpenAntigens downloads may include source-derived fields such as accessions, sequences, AlphaFold-derived confidence values, RefSeq-derived ortholog sequences, BLAST hit descriptions, and Open Targets disease associations. These fields remain subject to the terms and attribution expectations of the original source databases. Licensing, corrections, and takedown questions should be sent to Andre A. R. Teixeira at `andre.teixeira@proteininnovation.org`.

### Website analytics and privacy

Every rendered page includes the Plausible script assigned to `openantigens.org`. The portal records aggregate page traffic and clicks on static `.tsv` and `.json` files. Five fixed custom events cover browser-generated exports: full PDB, selected-region PDB, PNG, TSV copy, and FASTA copy. The custom events carry no properties.

The integration disables automatic outbound-link and form-submission measurement. It does not add cookies, session replay, advertising identifiers, or an ad-blocker bypass. Localhost, loopback, and `file://` previews are ignored by the Plausible script, so analytics verification must use the deployed domain.

The generated `privacy.html` page is the public disclosure and must ship with every release. The operational specification, Plausible account steps, event dictionary, and quarterly reporting rules are in [`docs/ANALYTICS.md`](ANALYTICS.md).

### Interactive Construct Builder support page

The portal includes `builder.html`, a dedicated support page for the Interactive Construct Builder.

This page explains:

- how sequence, structure, pLDDT, PAE, and Live Selected Region exports are linked
- how to adjust construct boundaries from sequence and plot selections
- how to interpret pLDDT as local residue-level confidence
- how to interpret PAE as relative-domain placement confidence
- why high pLDDT alone is not enough to define a construct boundary
- how low internal PAE supports a coherent construct region
- how high inter-block PAE can suggest flexible linkers or domain split points
- how cysteine and furin warnings should influence final construct choice
- when biological evidence should override automated structural heuristics

For construct design, use pLDDT and PAE together:

- pLDDT is best for identifying locally confident structure and low-confidence/disordered tails.
- PAE is best for identifying whether selected residues behave as one coherent unit or as multiple flexibly coupled units.
- A strong soluble construct usually has acceptable pLDDT across most residues and low PAE within the selected region.
- A large construct with high pLDDT but high PAE between internal blocks may still express, but it may behave as multiple independently moving domains.
- A boundary near low pLDDT and high local/inter-block PAE is often a good candidate split point, provided it does not disrupt required biology.

## Troubleshooting

### The portal still shows old data

The portal is static HTML. If reports changed, rebuild it:

```bash
.venv/bin/agdesign2 build-portal outputs/surfy_batch/batch_summary.json
```

### A report is still using old family context

Make sure both of these exist and are current:

- `outputs/surfy_batch/ortholog_reference_table.tsv`
- `outputs/surfy_batch/family_alignments/family_alignment_index.tsv`

Then rerun the target and rebuild the portal.

### Batch retries

To rerun rows that failed or are missing a local AlphaFold artifact:

```bash
.venv/bin/agdesign2 retry-batch data/inputs/accessible_human_topology_refined_accessible_only.csv --output-dir outputs/surfy_batch --summary-path outputs/surfy_batch/batch_summary.json --verbose
```

## Current Important Defaults

From `AnalysisConfig`:

- `plddt_structured_threshold = 70`
- `lenient_plddt_structured_threshold = 60`
- `pae_domain_threshold = 12`
- `pae_separation_threshold = 4`
- `max_disordered_gap = 8`
- `lenient_max_disordered_gap = 12`
- `min_structured_segment = 25`
- `min_domain_size = 40`
- `min_construct_length = 50`
- `blast_evalue = 1e-5`
- `blast_max_hits = 25`
- `max_family_context_members = 20`

## Code Pointers

If you want to inspect the exact implementation:

- pipeline orchestration: [pipeline.py](../src/agdesign2/pipeline.py)
- structure segmentation: [structure_utils.py](../src/agdesign2/structure_utils.py)
- family precompute: [family_alignments.py](../src/agdesign2/family_alignments.py)
- ortholog precompute: [ortholog_table.py](../src/agdesign2/ortholog_table.py)
- precomputed-reference loader: [precomputed.py](../src/agdesign2/precomputed.py)
- BLAST logic: [blast.py](../src/agdesign2/blast.py)
- report rendering: [reporting.py](../src/agdesign2/reporting.py)
- portal generation: [portal.py](../src/agdesign2/portal.py)
