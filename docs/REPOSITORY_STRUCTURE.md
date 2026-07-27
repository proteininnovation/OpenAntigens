# Repository structure

This repository contains the source code, reference inputs, tests, and documentation used to build OpenAntigens. Generated reports, downloaded biological databases, local caches, and rendered portal artifacts are runtime outputs and should not be committed.

## Source tree

- `src/agdesign2/`: Python package for data retrieval, construct analysis, module refreshes, static portal rendering, and dynamic portal serving.
- `tests/`: unit and regression tests for the package.
- `notebooks/`: exploratory scientist-facing notebooks. These should call the package rather than duplicate pipeline logic.
- `scripts/`: release packaging, figure generation, and measurement utilities.

## Reference inputs

- `data/inputs/accessible_human_topology_refined_accessible_only.csv`: canonical cell-surface and secreted target universe for the public portal.

## Branding and static assets

- `assets/branding/`: source branding files used when rendering the portal.
- Generated portal builds copy required assets into `outputs/<batch>/portal/`.

## Runtime outputs

The following are generated locally and ignored by Git:

- `outputs/`: JSON reports, Markdown reports, construct images, portal HTML, batch summaries, precomputed ortholog/paralog/family tables, logs, and refresh state.
- `.agdesign2/`: downloaded AlphaFold files, BLAST databases, HTTP caches, Matplotlib cache, and other local working data.
- `accessible_secreted_gpi_singlepass.tsv`: generated filtered working TSV for the active cell-surface/secreted portal refresh job.

## Recommended workflow

1. Edit code under `src/agdesign2/`.
2. Keep canonical target/source tables under `data/inputs/`.
3. Run local refresh/build commands into `outputs/`.
4. Publish the rendered portal and downloadable database artifacts from `outputs/`, but do not commit those generated files to the source repository.
5. For GitHub review, commit source, tests, docs, scripts, and small reference input files only.
