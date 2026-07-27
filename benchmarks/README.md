# Benchmarks

Committed baselines for measuring portal generation and front-end changes.

## Files

- `baseline_full_2026-05-30.json` - per-step generation timings captured from
  the 2026-05-30 full snapshot (`snapshot_benchmark.json`, 5,326 human targets).
  The dominant step is `build-ortholog-table` at ~10.8 h.
- `frontend_baseline.json` - page-weight metrics for that snapshot's
  `public_site`, produced by `scripts/measure_portal_frontend.py`. Headline:
  `portal-index-data.js` 19.7 MB raw / 2.5 MB gz; largest report 9.37 MB;
  site total 25.3 GB.

## Regenerating

### Front-end page weight

The measurement script requires Python 3.10 or newer and no project
dependencies. Run it against a locally built portal:

```bash
python3 scripts/measure_portal_frontend.py \
  outputs/surfy_batch/public_site --out benchmarks/frontend_<label>.json
```

Compare a post-fix run to the baseline by diffing the two JSON files (the
per-page `bytes` / `gzip_bytes` / `inline_payload_bytes` fields are the metrics).

### Generation timings

`build-fresh-snapshot` writes `snapshot_benchmark.json` next to the snapshot.
For fast iteration use the 20-target slice; for sign-off use the full build:

```bash
# fast loop (minutes); also the input set for the golden test
.venv/bin/agdesign2 build-fresh-snapshot --snapshot-name benchmark-20 \
  --limit 20 --force --verbose

# full reference build (hours); use a cold cache for reference numbers
.venv/bin/agdesign2 build-fresh-snapshot --snapshot-name openantigen-YYYY-MM-DD \
  --verbose
```

Copy the resulting `snapshot_benchmark.json` into `benchmarks/` with a dated
name and compare step durations against `baseline_full_2026-05-30.json`.

## Golden portal-render regression

The render safety net lives in `tests/test_golden_portal.py` with fixtures and a
SHA-256 manifest under `tests/golden/`. It is hermetic (no network) and runs in
under a second.

```bash
# check render output is unchanged
PYTHONPATH=src python3 -m unittest tests.test_golden_portal

# after an intended render change, regenerate the manifest and review its diff
UPDATE_GOLDEN=1 PYTHONPATH=src python3 -m unittest tests.test_golden_portal
```

The fixtures are a 6-target faithful slice of the 2026-05-30 snapshot. To extend
coverage to the structure-viewer path, drop a structured report JSON plus its
AlphaFold `.pdb` into `tests/golden/fixtures/` and regenerate; no test code
change needed.
