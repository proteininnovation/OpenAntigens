"""F09: a note must be written once and rendered once.

Two producers, both covered here. ``_append_note`` stacked a second copy of
every mouse note when ``_finalize_mouse_batch_results`` re-entered on a resumed
build (``skipped_existing`` rows), which put duplicate notes on 4,060 of the
4,702 mouse pages in the 2026-09-30 build. The renderer also collapses
duplicates, so pages already built that way read correctly before the rebuild.
"""

from __future__ import annotations

import unittest

from agdesign2.mouse_portal import _append_note
from agdesign2.portal import _dedupe_notes


class AppendNoteIdempotencyTests(unittest.TestCase):
    def test_repeated_append_writes_one_note(self) -> None:
        report: dict = {}
        for _ in range(3):
            _append_note(report, "info", "Generated from the mouse ortholog.", "OpenAntigens Mouse")
        self.assertEqual(len(report["notes"]), 1)

    def test_distinct_notes_are_all_kept(self) -> None:
        report: dict = {}
        _append_note(report, "info", "first", "OpenAntigens Mouse")
        _append_note(report, "info", "second", "OpenAntigens Mouse")
        _append_note(report, "warning", "first", "OpenAntigens Mouse")
        _append_note(report, "info", "first", "UniProt")
        self.assertEqual(len(report["notes"]), 4)

    def test_preexisting_duplicates_are_not_disturbed(self) -> None:
        # Already-built reports carry the duplicate; appending must not add a third.
        note = {"severity": "info", "message": "m", "source": "s"}
        report = {"notes": [dict(note), dict(note)]}
        _append_note(report, "info", "m", "s")
        self.assertEqual(len(report["notes"]), 2)


class RenderNoteCollapseTests(unittest.TestCase):
    def test_duplicate_notes_collapse_in_order(self) -> None:
        notes = [
            {"severity": "info", "message": "a", "source": "UniProt"},
            {"severity": "info", "message": "b", "source": "topology"},
            {"severity": "info", "message": "a", "source": "UniProt"},
        ]
        self.assertEqual(_dedupe_notes(notes), notes[:2])

    def test_non_dict_notes_pass_through(self) -> None:
        self.assertEqual(_dedupe_notes(["a", "a"]), ["a", "a"])

    def test_missing_notes_render_empty(self) -> None:
        self.assertEqual(_dedupe_notes(None), [])


if __name__ == "__main__":
    unittest.main()
