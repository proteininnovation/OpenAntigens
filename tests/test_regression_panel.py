import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from agdesign2.deployment import build_fresh_snapshot
from agdesign2.regression import (
    PanelRow,
    check_target,
    diff_against_baseline,
    load_panel,
    scaled_mouse_minimum,
    select_panel_rows,
)

REPO = Path(__file__).resolve().parents[1]


def _write(root: str, text: str) -> Path:
    path = Path(root) / "panel.tsv"
    path.write_text("gene_symbol\tcategory\texpect\twhy\n" + text, encoding="utf-8")
    return path


def _summary(**overrides):
    report = {
        "ectodomain": {"start": 207, "end": 395},
        "topology": {"topology_class": "secreted"},
        "construct_recommendations": [{}] * 3,
        "notes": [{"message": "Curated design region 207-395"}],
    }
    summary = {"topology_class": "secreted", "region": "207-395", "constructs": 3, "gpcr": False,
               "assembly": 0, "mouse": True, "notes": ["Curated design region 207-395"], "_report": report}
    summary.update(overrides)
    return summary


class PanelFileTests(unittest.TestCase):
    def test_shipped_panel_loads_and_covers_the_use_cases(self) -> None:
        panel = load_panel(REPO / "data" / "regression" / "panel.tsv")
        self.assertGreaterEqual(len(panel), 100)
        categories = {row.category for row in panel}
        for required in ("gpcr", "single_pass", "multipass_compact", "multipass_mixed", "secreted", "gpi",
                         "secreted_fallback", "single_pass_inferred", "curated_override", "no_design_region"):
            self.assertIn(required, categories)

    def test_panel_errors_surface(self) -> None:
        with TemporaryDirectory() as root:
            with self.assertRaisesRegex(ValueError, "unknown category"):
                load_panel(_write(root, "EGFR\tnot_a_category\t\t\n"))
            with self.assertRaisesRegex(ValueError, "listed twice"):
                load_panel(_write(root, "EGFR\tsingle_pass\t\t\negfr\tsingle_pass\t\t\n"))
            with self.assertRaisesRegex(ValueError, "unknown expectation token"):
                load_panel(_write(root, "EGFR\tsingle_pass\tregion_starts=3\t\n"))

    def test_selection_keeps_panel_order_and_refuses_missing_targets(self) -> None:
        rows = [{"gene_symbol": "EGFR"}, {"gene_symbol": "IL6"}, {"gene_symbol": "VEGFA"}]
        panel = [PanelRow("VEGFA", "curated_override", (), ""), PanelRow("EGFR", "single_pass", (), "")]
        self.assertEqual([row["gene_symbol"] for row in select_panel_rows(rows, panel)], ["VEGFA", "EGFR"])
        with self.assertRaisesRegex(ValueError, "TRBC1"):
            select_panel_rows(rows, panel + [PanelRow("TRBC1", "ig_constant", (), "")])

    def test_mouse_minimum_scales_with_the_production_ratio(self) -> None:
        self.assertEqual(scaled_mouse_minimum(123, human_minimum=5000, mouse_minimum=4500), 110)


class ExpectationTests(unittest.TestCase):
    def test_tokens_pass_and_fail_with_the_observed_value(self) -> None:
        row = PanelRow("VEGFA", "curated_override",
                       ("class=secreted", "region", "region_start=207", "constructs>=1", "mouse", "note~curated_design"), "")
        self.assertEqual(check_target(row, _summary()).failures, [])
        failing = PanelRow("VEGFA", "curated_override", ("class=single_pass", "region_start=1", "constructs=0", "gpcr"), "")
        failures = check_target(failing, _summary()).failures
        self.assertEqual(len(failures), 4)
        self.assertIn("observed secreted", failures[0])
        self.assertIn("observed 207", failures[1])

    def test_missing_report_is_a_failure(self) -> None:
        self.assertEqual(check_target(PanelRow("X", "secreted", (), ""), None).failures, ["no human report in this build"])

    def test_baseline_diff_lists_changed_fields_only(self) -> None:
        before = {"VEGFA": _summary(region="", constructs=0, topology_class="non_membrane")}
        now = {"VEGFA": _summary()}
        changes = diff_against_baseline(now, before, ["VEGFA"])
        self.assertEqual({key for _, key, _, _ in changes}, {"region", "constructs", "topology_class"})


class BuildIntegrationTests(unittest.TestCase):
    def test_panel_build_plans_the_production_steps_including_the_gate(self) -> None:
        with TemporaryDirectory() as root:
            result = build_fresh_snapshot(
                target_panel=REPO / "data" / "regression" / "panel.tsv", snapshot_root=root,
                snapshot_name="panel", dry_run=True,
            )
            production = build_fresh_snapshot(snapshot_root=root, snapshot_name="full", dry_run=True)
        self.assertEqual([t["step"] for t in result.timings], [t["step"] for t in production.timings])
        self.assertIn("validate-release-data", [t["step"] for t in result.timings])

    def test_panel_cannot_be_combined_with_other_selections(self) -> None:
        with TemporaryDirectory() as root, self.assertRaisesRegex(ValueError, "--target-panel"):
            build_fresh_snapshot(target_panel=REPO / "data" / "regression" / "panel.tsv", limit=5,
                                 snapshot_root=root, dry_run=True)


if __name__ == "__main__":
    unittest.main()
