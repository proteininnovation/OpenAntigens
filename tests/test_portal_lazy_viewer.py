"""The detail-page 3Dmol structure viewer initializes lazily (on scroll /
first interaction) rather than eagerly on page load. Every viewer.* call site
guards on a null viewer, so deferral cannot crash an early interaction."""

from __future__ import annotations

import re
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.portal import _render_structure_widget

_PDB = (
    "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 95.00           C\n"
    "ATOM      2  CA  ALA A   2       1.000   0.000   0.000  1.00 95.00           C\nEND\n"
)


def _render_widget(root: Path) -> str:
    structures = root / "portal" / "structures"
    structures.mkdir(parents=True)
    pdb = structures / "test_human.pdb"
    pdb.write_text(_PDB, encoding="utf-8")
    report = {
        "target": {"entry_name": "TEST_HUMAN", "accession": "P00001", "sequence": "AA"},
        "ectodomain": {"start": 1, "end": 1, "label": "e", "source": "t"},
        "topology": {"topology_class": "multipass_mixed"},
    }
    return _render_structure_widget({"portal_structure_path": str(pdb)}, report, batch_dir=root)


class LazyStructureViewerTests(unittest.TestCase):
    def test_viewer_init_is_deferred_not_eager(self) -> None:
        with TemporaryDirectory() as tmp:
            html = _render_widget(Path(tmp))

        self.assertIn("function initStructureViewer()", html)
        self.assertIn("new IntersectionObserver(", html)
        # The WebGL viewer creation must live INSIDE initStructureViewer, not run
        # at top level on load.
        init_pos = html.index("function initStructureViewer()")
        create_pos = html.index("createViewer(viewerEl")
        self.assertGreater(create_pos, init_pos, "createViewer should be inside initStructureViewer")
        # First-interaction (capture-phase) init as a fallback to the observer.
        self.assertIn("capture: true", html)
        # Immediate-init fallback where IntersectionObserver is unavailable.
        self.assertRegex(html, r"\}\s*else\s*\{\s*initStructureViewer\(\);")

    def test_all_viewer_calls_remain_guarded(self) -> None:
        # Defensive invariant: deferral is only safe because viewer.* is never
        # called without a null check. Every JS function that uses `viewer.`
        # must contain an `if (!viewer` guard.
        with TemporaryDirectory() as tmp:
            html = _render_widget(Path(tmp))
        # Split on function declarations and check each body that uses viewer.*
        parts = re.split(r"\n        function (\w+)\(", html)
        checked = 0
        for i in range(1, len(parts), 2):
            name, body = parts[i], parts[i + 1]
            # Only consider the body up to the next top-level function (already
            # split) — guard appears within the first lines if present.
            if "viewer." in body and name not in {"initStructureViewer", "attemptStructureViewerInit"}:
                self.assertIn("if (!viewer", body, f"{name} uses viewer.* without a null guard")
                checked += 1
        self.assertGreater(checked, 0, "expected to find guarded viewer functions")


if __name__ == "__main__":
    unittest.main()
