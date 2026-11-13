# Canonical repository

- This checkout belongs exclusively to `https://github.com/proteininnovation/OpenAntigens`.
- Before any fetch, push, pull request, merge, release, or deployment, verify that `origin` resolves to `proteininnovation/OpenAntigens` and that the target branch is in that repository.
- Never push, open pull requests, merge, release, or deploy from a personal or private fork. If the remote differs from the canonical repository, stop and correct it before doing any work.
- Changes intended for `main` must go through a pull request in `proteininnovation/OpenAntigens` unless the user explicitly requests a different workflow.
- Never create or use a branch whose name begins with `codex`. Use an informative prefix such as `feat/`, `fix/`, `perf/`, `docs/`, `test/`, `refactor/`, or `chore/`.
- Keep `core.hooksPath` set to `.githooks`; the tracked pre-push hook rejects every non-canonical destination.

# Verify every change with the regression panel

Any change to the pipeline (`src/agdesign2/`), the portal renderers, the build, or the target inputs MUST pass the regression panel before you report it done, open a PR for it, run a full rebuild, or deploy. Unit tests alone do not count: the panel builds 123 real targets through the production build and checks both portals. Full guide: `docs/REGRESSION_BUILD.md`.

```bash
scripts/regression.sh build snapshots/<latest production snapshot>       # ~20 min with warm caches
scripts/regression.sh check snapshots/regression-<date>-<sha> <baseline>   # ~10 min
```

- The panel build is the production command, `build-fresh-snapshot --resume --reanalyze-reports --ortholog-jobs 6`, plus `--target-panel data/regression/panel.tsv`. Never test a change through a separate build path, a custom script, or a hand-rendered page instead; that is not what production runs.
- Use the previous regression build as the baseline, or the latest production snapshot when there is none. Every change the check lists must be one the work intended. Explain any other change before reporting.
- All three layers must pass: panel expectations (`check-regression-panel`), `scripts/release_checks.py`, and the Playwright walk-through of both portals (`scripts/check_regression_site.cjs`).
- A failed expectation is either a regression or a wrong expectation. Read the report JSON or the page before changing `data/regression/panel.tsv`, and say which it was.
- A new feature or a fixed case gets a panel row (and, where needed, a category with default expectations in `src/agdesign2/regression.py`), so the next change cannot silently undo it.
- In the report, give the commands, the snapshot name, each layer's result, and the change list against the baseline.

# Human and mouse portal parity

- Every change to shared page layout, navigation, copy, citations, agent guidance, documentation, or interactive behavior MUST apply to both human and mouse portals. A change is incomplete until both portals are checked. Do not update only the portal named in a request when the change also applies to the other.
- Differences MUST have a species-specific reason: target identity and provenance, available data and homolog mappings, species labels and theme, relative links, or supported features (for example, human disease context). Preserve those differences; do not copy human data or claims into mouse pages to make them look alike.
- Reuse shared renderers and content for common behavior. Inspect both build paths and their overrides before editing; remove stale overrides that prevent shared changes from reaching one portal.
- Verify the actual human and mouse build paths, not only the shared renderer or a custom preview. Refresh both local previews and inspect affected pages, copy, navigation, and links in each portal. Run relevant existing checks and add a regression check when it protects against divergence.
- Before reporting completion, publishing, or deploying, confirm that both portal artifacts contain the changes and state any intentional species-specific differences. A successful human build or preview alone does not establish mouse parity.
