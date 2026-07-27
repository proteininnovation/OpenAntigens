from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from scripts.generate_alphafold3_webserver_inputs import load_reports


class AlphaFoldServerInputTests(unittest.TestCase):
    def test_report_manifest_paths_are_relative(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            reports_dir = Path(tmpdir) / "reports"
            reports_dir.mkdir()
            report_path = reports_dir / "test_human_report.json"
            report_path.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "P00001",
                            "gene_symbol": "TEST",
                            "sequence": "ACDE",
                        }
                    }
                ),
                encoding="utf-8",
            )

            reports = load_reports(reports_dir)

        self.assertEqual(reports[0]["report_path"], "test_human_report.json")


if __name__ == "__main__":
    unittest.main()
