"""refresh_report_modules must only rewrite a report when its content
actually changed (handlers previously always reported 'changed')."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import agdesign2.module_refresh as mr


def _setup(batch_dir: Path) -> Path:
    # Report written in a deliberately non-canonical format so an unnecessary
    # rewrite (to indent=2, sort_keys) would change the bytes and be detected.
    report = {"target": {"accession": "P1", "entry_name": "AAA_HUMAN", "gene_symbol": "AAA"}, "notes": []}
    report_path = batch_dir / "aaa_human_report.json"
    report_path.write_text(json.dumps(report, indent=4), encoding="utf-8")
    summary = [{"query": "AAA", "resolved_entry_name": "AAA_HUMAN",
                "json_report": "aaa_human_report.json", "status": "ok"}]
    (batch_dir / "batch_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return report_path


class ModuleRefreshDirtyCheckTests(unittest.TestCase):
    def test_noop_handler_does_not_rewrite_report(self) -> None:
        with TemporaryDirectory() as tmp:
            batch_dir = Path(tmp)
            report_path = _setup(batch_dir)
            before = report_path.read_bytes()
            # Returns True (like the real cross-reactivity/family-context
            # handlers) but changes nothing — the old code rewrote/reformatted
            # the report anyway; the fix must not.
            with mock.patch.dict(mr._MODULE_HANDLERS, {"noop": lambda tool, payload: True}):
                _, _, updated = mr.refresh_report_modules(
                    batch_dir / "batch_summary.json", modules=["noop"], analyzer=object()
                )
            self.assertEqual(updated, 0)
            self.assertEqual(report_path.read_bytes(), before)  # not rewritten/reformatted

    def test_mutating_handler_rewrites_report(self) -> None:
        def mutate(tool, payload):
            payload["notes"] = (payload.get("notes") or []) + [{"added": True}]

        with TemporaryDirectory() as tmp:
            batch_dir = Path(tmp)
            report_path = _setup(batch_dir)
            with mock.patch.dict(mr._MODULE_HANDLERS, {"mut": mutate}):
                _, _, updated = mr.refresh_report_modules(
                    batch_dir / "batch_summary.json", modules=["mut"], analyzer=object()
                )
            self.assertEqual(updated, 1)
            persisted = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(persisted["notes"], [{"added": True}])


if __name__ == "__main__":
    unittest.main()
