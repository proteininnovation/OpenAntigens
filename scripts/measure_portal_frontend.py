#!/usr/bin/env python3
"""Measure end-user portal page weight for a built ``public_site``/``portal`` directory.

Records, per representative page, the raw and gzipped transfer size, how many inline
``<script>`` blocks the page carries, and how many bytes of that page are
inline JSON/JS payload (the data the browser must parse before first paint).
These are the metrics the front-end fixes (index blob, detail-page structure
inlining, flat files) are expected to move.

It is intentionally dependency-free: it only reads files and runs gzip, so it
works on any host that has Python 3.10+ (including the remote snapshot host).
Browser timing metrics (FCP/TTI/DOM nodes) require Playwright and a served
site; this tool deliberately stays at the static-size layer so the baseline is
reproducible without a browser. Capture timings separately once the harness
serves a site.

Usage:
    python3 scripts/measure_portal_frontend.py /path/to/public_site \
        --out benchmarks/frontend_baseline.json

The chosen pages are the ones called out in the review: the index (19 MB blob
in the 2026-05-30 snapshot), the largest report (9.4 MB), a typical report, and
the static support pages. Large flat downloads are reported separately because
they are shipped/linked but not parsed on page load.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from pathlib import Path
from typing import Any

# Inline <script>...</script> blocks that have NO src attribute carry payload
# the browser must parse inline. Capture their combined byte length.
_INLINE_SCRIPT = re.compile(
    rb"<script(?![^>]*\bsrc=)[^>]*>(.*?)</script>", re.IGNORECASE | re.DOTALL
)
_ANY_SCRIPT_OPEN = re.compile(rb"<script\b", re.IGNORECASE)
_IMG_TAG = re.compile(rb"<img\b[^>]*>", re.IGNORECASE)
_LAZY_IMG = re.compile(rb"<img\b[^>]*\bloading=[\"']lazy[\"']", re.IGNORECASE)


def _gzip_size(data: bytes) -> int:
    return len(gzip.compress(data, 9))


def measure_file(path: Path) -> dict[str, Any]:
    data = path.read_bytes()
    inline_blocks = _INLINE_SCRIPT.findall(data)
    inline_payload_bytes = sum(len(b) for b in inline_blocks)
    imgs = _IMG_TAG.findall(data)
    lazy = _LAZY_IMG.findall(data)
    return {
        "path": path.name,
        "bytes": len(data),
        "gzip_bytes": _gzip_size(data),
        "script_tags": len(_ANY_SCRIPT_OPEN.findall(data)),
        "inline_script_blocks": len(inline_blocks),
        "inline_payload_bytes": inline_payload_bytes,
        "img_tags": len(imgs),
        "img_lazy": len(lazy),
        "img_eager": len(imgs) - len(lazy),
    }


def _largest(glob_dir: Path, pattern: str) -> Path | None:
    matches = sorted(glob_dir.glob(pattern), key=lambda p: p.stat().st_size, reverse=True)
    return matches[0] if matches else None


def _median(glob_dir: Path, pattern: str) -> Path | None:
    matches = sorted(glob_dir.glob(pattern), key=lambda p: p.stat().st_size)
    return matches[len(matches) // 2] if matches else None


def build_report(site: Path) -> dict[str, Any]:
    pages: list[dict[str, Any]] = []

    def add(label: str, path: Path | None) -> None:
        if path is not None and path.is_file():
            entry = {"label": label}
            entry.update(measure_file(path))
            pages.append(entry)

    # Top-level pages and the index data blob.
    add("index", site / "index.html")
    add("index_data_js", site / "portal-index-data.js")
    add("index_js", site / "portal-index.js")
    for name in ("constructs.html", "calculator.html", "methods.html",
                 "help.html", "downloads.html", "builder.html"):
        add(name.split(".")[0], site / name)

    reports = site / "reports"
    if reports.is_dir():
        add("report_largest", _largest(reports, "*.html"))
        add("report_median", _median(reports, "*.html"))
    scripts = site / "report_scripts"
    if scripts.is_dir():
        add("report_script_largest", _largest(scripts, "*.js"))
        add("report_script_median", _median(scripts, "*.js"))

    # Flat files shipped/linked but not parsed on load — report size only.
    flat: list[dict[str, Any]] = []
    for name in ("open_targets_disease_associations.json",
                 "open_targets_disease_associations.tsv"):
        p = site / name
        if p.is_file():
            flat.append({"path": name, "bytes": p.stat().st_size})

    # Whole-site footprint for archive and static-host sizing.
    total = 0
    for p in site.rglob("*"):
        if p.is_file():
            total += p.stat().st_size

    return {
        "site": str(site),
        "site_total_bytes": total,
        "pages": pages,
        "flat_files": flat,
        "notes": "Browser FCP/TTI/DOM-node metrics require a separate browser benchmark.",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("site", type=Path, help="public_site or portal directory")
    parser.add_argument("--out", type=Path, default=None,
                        help="write JSON report here (default: stdout)")
    args = parser.parse_args(argv)

    if not args.site.is_dir():
        print(f"error: {args.site} is not a directory", file=sys.stderr)
        return 2

    report = build_report(args.site)
    text = json.dumps(report, indent=2)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text + "\n", encoding="utf-8")
        # Human-readable digest to stderr.
        print(f"wrote {args.out}", file=sys.stderr)
        for page in report["pages"]:
            print(f"  {page['label']:24} {page['bytes']:>12,} B  "
                  f"gz {page['gzip_bytes']:>11,} B  "
                  f"inline {page['inline_payload_bytes']:>11,} B",
                  file=sys.stderr)
        for flat in report["flat_files"]:
            print(f"  [flat] {flat['path']:32} {flat['bytes']:>14,} B", file=sys.stderr)
        print(f"  site_total {report['site_total_bytes']:>26,} B", file=sys.stderr)
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
