"""Golden-output regression for portal generation.

This renders a fixed, committed set of real
report-JSON fixtures through ``build_portal`` with **no network access**
(``refresh_reports=False``, ``refresh_stale_reports=False``,
``include_disease_context=False``) and compares every produced file against a
committed SHA-256 manifest.

A refactor that is meant to be behaviour-preserving must leave this test green.
A refactor that intentionally changes rendered output regenerates the manifest
(``UPDATE_GOLDEN=1``) and the diff in the manifest is reviewed line-by-line in
the PR — the point is that output changes become explicit and auditable.

Fixtures live in ``tests/golden/fixtures/`` (a 6-target faithful slice of the
2026-05-30 snapshot). They currently exercise the non-structure render paths
(constructs, family context, cross-reactivity, PTMs, cysteines, downloads,
index). To extend coverage to the structure-viewer path, drop a structured
report JSON plus its AlphaFold ``.pdb`` into the fixtures dir and regenerate;
no test code change is needed.

Run:
    PYTHONPATH=src python3 -m unittest tests.test_golden_portal
Regenerate the manifest after an intended change:
    UPDATE_GOLDEN=1 PYTHONPATH=src python3 -m unittest tests.test_golden_portal
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

FIXTURES = ROOT / "tests" / "golden" / "fixtures"
MANIFEST = ROOT / "tests" / "golden" / "manifest.json"

from agdesign2.portal import build_portal  # noqa: E402

# Files whose bytes legitimately differ between machines/runs (timestamps,
# absolute paths, tool versions). They are still built and checked for
# existence, but excluded from the byte-for-byte hash comparison.
_VOLATILE = {
    "portal_metadata.json",   # embeds build date / version
    "download_manifest.json",  # embeds generated_at timestamp
}


def _render_to(tmp: Path) -> Path:
    """Build the portal from the committed fixtures into ``tmp``; return portal dir."""
    batch_dir = tmp / "batch"
    batch_dir.mkdir(parents=True, exist_ok=True)
    for f in FIXTURES.iterdir():
        if f.is_file():
            shutil.copy(f, batch_dir / f.name)
    portal_dir = tmp / "portal"
    # Pin the footer "Release snapshot" date via the build-date hook so rendered
    # bytes are deterministic across days; otherwise every page's hash changes
    # daily and the manifest comparison is unrunnable after the day it was made.
    prior_build_date = os.environ.get("OPENANTIGEN_BUILD_DATE")
    os.environ["OPENANTIGEN_BUILD_DATE"] = "2026-01-01"
    try:
        # build_portal returns the rendered index.html path, not the directory.
        build_portal(
            batch_dir / "batch_summary.json",
            output_dir=portal_dir,
            refresh_reports=False,
            refresh_stale_reports=False,
            include_disease_context=False,
            bundle_vendor_assets=False,
            # No live PubTator fetch: the count comes from the committed fixtures,
            # so the render is deterministic across machines and over time.
            fetch_literature=False,
            # Mouse sibling link points at files that do not exist in this slice;
            # disabling keeps the rendered nav deterministic and self-contained.
            sibling_link=None,
            detail_sibling_link=None,
        )
    finally:
        if prior_build_date is None:
            os.environ.pop("OPENANTIGEN_BUILD_DATE", None)
        else:
            os.environ["OPENANTIGEN_BUILD_DATE"] = prior_build_date
    return portal_dir


def _manifest_for(portal_dir: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for path in sorted(portal_dir.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(portal_dir).as_posix()
        if path.name in _VOLATILE:
            out[rel] = "<volatile>"
            continue
        out[rel] = hashlib.sha256(path.read_bytes()).hexdigest()
    return out


class GoldenPortalTests(unittest.TestCase):
    def test_portal_render_matches_golden(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            portal_dir = _render_to(Path(raw))
            produced = _manifest_for(portal_dir)

            if os.environ.get("UPDATE_GOLDEN"):
                MANIFEST.write_text(json.dumps(produced, indent=2, sort_keys=True) + "\n",
                                    encoding="utf-8")
                self.skipTest(f"UPDATE_GOLDEN set; wrote {len(produced)} entries to {MANIFEST}")

            self.assertTrue(
                MANIFEST.exists(),
                "tests/golden/manifest.json missing — run with UPDATE_GOLDEN=1 once to create it.",
            )
            golden = json.loads(MANIFEST.read_text(encoding="utf-8"))

            produced_keys = set(produced)
            golden_keys = set(golden)
            self.assertEqual(
                produced_keys, golden_keys,
                f"file set changed.\n  added: {sorted(produced_keys - golden_keys)}\n"
                f"  removed: {sorted(golden_keys - produced_keys)}",
            )
            mismatches = [
                k for k in golden_keys
                if golden[k] != "<volatile>" and produced[k] != golden[k]
            ]
            self.assertFalse(
                mismatches,
                "rendered bytes changed for:\n  " + "\n  ".join(sorted(mismatches))
                + "\nIf intended, rerun with UPDATE_GOLDEN=1 and review the manifest diff.",
            )

    def test_fixtures_present(self) -> None:
        reports = list(FIXTURES.glob("*_report.json"))
        self.assertGreaterEqual(len(reports), 6, "expected the committed fixture slice")
        self.assertTrue((FIXTURES / "batch_summary.json").exists())


if __name__ == "__main__":
    unittest.main()
