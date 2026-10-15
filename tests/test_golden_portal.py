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

Fixtures live in ``tests/golden/fixtures/`` (human) and
``tests/golden/fixtures/mouse/`` (mouse), copied verbatim from the snapshots
behind the 2026-09-30 build: ``openantigen-2026-09-15-cah4-fix-local`` for
human, ``openantigen-2026-07-23`` for mouse. Between them they exercise
constructs, cysteines, InterPro annotations, cross-reactivity, family context,
GPCR annotation and engineering variants, multipass and GPI topology, homology
rows (including blank ones), PTMs and the empty-state paths — in **both**
portal themes, so a render change that reaches only one portal fails here.
To extend coverage to the structure-viewer path, drop a structured report JSON
plus its AlphaFold ``.pdb`` into the fixtures dir and regenerate; no test code
change is needed.

Run:
    PYTHONPATH=src python3 -m unittest tests.test_golden_portal
Regenerate the manifest after an intended change:
    UPDATE_GOLDEN=1 PYTHONPATH=src python3 -m unittest tests.test_golden_portal
"""

from __future__ import annotations

import hashlib
import json
import os
import re
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
MOUSE_FIXTURES = FIXTURES / "mouse"
MANIFEST = ROOT / "tests" / "golden" / "manifest.json"

from agdesign2.portal import build_portal  # noqa: E402

# Files whose bytes legitimately differ between machines/runs (timestamps,
# absolute paths, tool versions). They are still built and checked for
# existence, but excluded from the byte-for-byte hash comparison.
_VOLATILE = {
    "portal_metadata.json",   # embeds build date / version
    "download_manifest.json",  # embeds generated_at timestamp
}

# ``_PORTAL_CSS_VERSION`` (portal.py:40) is sha256 of portal.py's own bytes and
# rides on every page as ``?v=<16 hex>``. Left in, a one-character edit to
# portal.py "changes" every rendered file and buries the handful that really
# did change — which trains the reader to regenerate the manifest unread. The
# stylesheet's own content is still covered: portal.css is hashed as a file.
_CSS_VERSION_QUERY = re.compile(rb"\?v=[0-9a-f]{16}")


def _render_to(tmp: Path) -> tuple[Path, Path]:
    """Build both portals from the committed fixtures; return (human, mouse) dirs."""
    # Pin the footer "Release snapshot" date via the build-date hook so rendered
    # bytes are deterministic across days; otherwise every page's hash changes
    # daily and the manifest comparison is unrunnable after the day it was made.
    prior_build_date = os.environ.get("OPENANTIGEN_BUILD_DATE")
    prior_alignment_backend = os.environ.get("AGDESIGN2_ALIGNMENT_BACKEND")
    os.environ["OPENANTIGEN_BUILD_DATE"] = "2026-01-01"
    os.environ["AGDESIGN2_ALIGNMENT_BACKEND"] = "python"
    try:
        human_dir = _build_one(
            tmp / "human",
            FIXTURES,
            portal_title="OpenAntigens",
            theme="human",
            sibling_link=("OpenAntigens Mouse", "mouse/index.html"),
            detail_sibling_link=("OpenAntigens Mouse", "../mouse/index.html"),
        )
        mouse_dir = _build_one(
            tmp / "mouse",
            MOUSE_FIXTURES,
            portal_title="OpenAntigens Mouse",
            portal_subtitle="Antigen constructs for mouse orthologs of human targets.",
            theme="mouse",
            sibling_link=("OpenAntigens Human", "../index.html"),
            detail_sibling_link=("OpenAntigens Human", "../../index.html"),
        )
    finally:
        if prior_build_date is None:
            os.environ.pop("OPENANTIGEN_BUILD_DATE", None)
        else:
            os.environ["OPENANTIGEN_BUILD_DATE"] = prior_build_date
        if prior_alignment_backend is None:
            os.environ.pop("AGDESIGN2_ALIGNMENT_BACKEND", None)
        else:
            os.environ["AGDESIGN2_ALIGNMENT_BACKEND"] = prior_alignment_backend
    return human_dir, mouse_dir


def _build_one(
    tmp: Path,
    fixtures: Path,
    *,
    portal_title: str,
    theme: str,
    sibling_link: tuple[str, str],
    detail_sibling_link: tuple[str, str],
    portal_subtitle: str | None = None,
) -> Path:
    batch_dir = tmp / "batch"
    batch_dir.mkdir(parents=True, exist_ok=True)
    for f in fixtures.iterdir():
        if f.is_file():
            shutil.copy(f, batch_dir / f.name)
    portal_dir = tmp / "portal"
    # build_portal returns the rendered index.html path, not the directory.
    build_portal(
        batch_dir / "batch_summary.json",
        output_dir=portal_dir,
        refresh_reports=False,
        refresh_stale_reports=False,
        include_disease_context=False,
        bundle_vendor_assets=False,
        portal_title=portal_title,
        portal_subtitle=portal_subtitle,
        theme=theme,
        # No live PubTator fetch: the count comes from the committed fixtures,
        # so the render is deterministic across machines and over time.
        fetch_literature=False,
        # The real cross-portal links, as mouse_portal.py:218 and
        # deployment.py:364 pass them: the human<->mouse nav is exactly the
        # surface AGENTS.md parity is about, so it has to be in the manifest.
        # The targets do not exist in this slice; only the markup is compared.
        sibling_link=sibling_link,
        detail_sibling_link=detail_sibling_link,
    )
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
        out[rel] = hashlib.sha256(
            _CSS_VERSION_QUERY.sub(b"?v=<css-version>", path.read_bytes())
        ).hexdigest()
    return out


class GoldenPortalTests(unittest.TestCase):
    def test_portal_render_matches_golden(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            human_dir, mouse_dir = _render_to(Path(raw))
            produced = {f"human/{k}": v for k, v in _manifest_for(human_dir).items()}
            produced.update({f"mouse/{k}": v for k, v in _manifest_for(mouse_dir).items()})

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
        self.assertGreaterEqual(len(reports), 11, "expected the committed fixture slice")
        self.assertTrue((FIXTURES / "batch_summary.json").exists())
        mouse_reports = list(MOUSE_FIXTURES.glob("*_report.json"))
        self.assertGreaterEqual(len(mouse_reports), 4, "expected the mouse fixture slice")
        self.assertTrue((MOUSE_FIXTURES / "batch_summary.json").exists())

    def test_fixtures_cover_the_render_paths(self) -> None:
        """The slice must keep exercising the paths the remediation PRs touch.

        Without this the manifest can stay green while a fixture is swapped for
        one that never reaches constructs, cysteines or the GPCR sections.
        """
        wanted = {
            "construct_details",
            "cysteine_analysis",
            "interpro_annotations",
            "cross_reactivity_hits",
            "family_context",
            "gpcr_annotation",
            "gpcr_engineering_variants",
            "ectodomain_homology",
            "ptms",
            "assembly_requirements",
        }
        for label, fixtures in (("human", FIXTURES), ("mouse", MOUSE_FIXTURES)):
            seen: set[str] = set()
            multipass = False
            for path in fixtures.glob("*_report.json"):
                report = json.loads(path.read_text(encoding="utf-8"))
                seen.update(key for key in wanted if report.get(key))
                topology_class = str((report.get("topology") or {}).get("topology_class") or "")
                multipass = multipass or "multi" in topology_class
            self.assertEqual(wanted - seen, set(), f"{label} fixtures stopped covering these")
            self.assertTrue(multipass, f"{label} fixtures lost the multipass topology case")


if __name__ == "__main__":
    unittest.main()
