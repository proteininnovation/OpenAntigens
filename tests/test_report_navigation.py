"""Report navigation targets include only sections actually rendered."""
import re
import tempfile
import unittest
from pathlib import Path

from agdesign2.portal import render_detail_page
from tests.test_portal_db import _write_fixture


SECTION_LABELS = [
    "Overview", "Obligatory complex context", "Construct Builder", "Construct Summary",
    "GPCR Engineering", "Disease Associations", "Homology", "Cross-reactivity", "Cysteines",
    "Furin Cleavage Sites", "Post-translational Modifications", "Family Context",
    "Interactions / Assembly", "Membrane Engineering", "GPCR Annotation", "InterPro / Pfam",
    "Construct Details", "Notes",
]


def with_optional_sections(report):
    return {
        **report,
        "assembly_requirements": [{"obligatory": True, "partners": ["TEST_PARTNER"], "summary": "Synthetic test context."}],
        "advanced_membrane_suggestions": [{"title": "Test membrane suggestion"}],
        "gpcr_engineering_variants": [{"name": "test_gpcr_variant", "start": 1, "end": 10}],
        "gpcr_annotation": {"entry_name": "test_human"},
    }


def section_targets(html):
    return re.findall(r'<(?:section|article)[^>]*id="([^"]+)" data-report-section="([^"]+)"', html)


class ReportNavigationTests(unittest.TestCase):
    def test_sections_have_unique_stable_targets_and_keep_empty_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, report = _write_fixture(root)
            for optional in (False, True):
                for disease in (False, True):
                    with self.subTest(optional=optional, disease=disease):
                        html = render_detail_page(
                            {}, with_optional_sections(report) if optional else report,
                            batch_dir=root, include_disease_context=disease,
                        )
                        omitted = set() if optional else {"Obligatory complex context", "GPCR Engineering", "Membrane Engineering", "GPCR Annotation"}
                        if not disease:
                            omitted.add("Disease Associations")
                        targets = section_targets(html)
                        self.assertEqual([label for _, label in targets], [label for label in SECTION_LABELS if label not in omitted])
                        self.assertEqual(len({anchor for anchor, _ in targets}), len(targets))
                        for anchor, _ in targets:
                            self.assertEqual(html.count(f'id="{anchor}"'), 1)
                        self.assertIn(("interactive-construct-builder", "Construct Builder"), targets)
                        self.assertIn(("construct-details", "Construct Details"), targets)
                        self.assertIn("No constructs available.", html)
                        self.assertIn('popovertarget="report-sections-menu"', html)
                        self.assertIn('<div id="report-sections-menu" popover>', html)


if __name__ == "__main__":
    unittest.main()
