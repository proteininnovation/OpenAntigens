from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.portal import _ANALYTICS_EVENT_NAMES, _portal_analytics_head, render_privacy_page


class PortalAnalyticsTests(unittest.TestCase):
    def test_analytics_head_uses_fixed_privacy_conscious_configuration(self) -> None:
        html = _portal_analytics_head()

        self.assertEqual(
            html.count('src="https://plausible.io/js/pa-Mda6DF7h4_8b17ZUlD_9E.js"'),
            1,
        )
        self.assertIn('fileExtensions: ["tsv", "json"]', html)
        self.assertIn("formSubmissions: false", html)
        self.assertIn("outboundLinks: false", html)
        self.assertNotIn("props:", html)
        self.assertIn(
            json.dumps(list(_ANALYTICS_EVENT_NAMES), separators=(",", ":")),
            html,
        )
        self.assertIn("throw new Error(`Unknown OpenAntigens analytics event:", html)

    def test_privacy_page_discloses_measurement_limits_and_reporting_policy(self) -> None:
        html = render_privacy_page()

        for expected in (
            "Plausible Analytics",
            "does not use analytics cookies",
            "do not send site-search terms",
            "Counts are approximate",
            "fewer than five unique visitors",
            "andre.teixeira@proteininnovation.org",
        ):
            self.assertIn(expected, html)
        self.assertEqual(html.count("data-openantigens-analytics"), 1)
        self.assertIn('href="privacy.html"', html)


if __name__ == "__main__":
    unittest.main()
