"""PubTator counts: retry transient failures, never hide a failure, follow renamed genes."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from agdesign2 import portal
from agdesign2.cache import FileCache
from agdesign2.exceptions import ExternalServiceError
from agdesign2.release_validation import _validate_catalog


class FlakyHttp:
    def __init__(self, responses):
        self.responses = list(responses)
        self.calls = 0
        self.urls = []

    def fetch_json(self, url, **_):
        self.calls += 1
        self.urls.append(url)
        response = self.responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response


class PubTatorSearchCountTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.cache = FileCache(Path(self.tmp.name))

    def test_transient_failures_are_retried_with_backoff(self):
        http = FlakyHttp([ExternalServiceError("HTTP error: 429"), ValueError("bad json"), {"count": 5676}])
        with patch.object(portal, "sleep") as sleep:
            self.assertEqual(portal._pubtator_search_count("APLN", cache=self.cache, http=http), 5676)
        self.assertEqual([call.args[0] for call in sleep.call_args_list], [2, 5])

    def test_hyphenated_symbol_queries_the_underscore_gene_id(self):
        http = FlakyHttp([{"count": 542}])
        self.assertEqual(portal._pubtator_search_count("ERVFRD-1", cache=self.cache, http=http), 542)
        self.assertTrue(http.urls[0].endswith("text=%40GENE_ERVFRD_1"), http.urls[0])

    def test_raises_once_retries_are_spent(self):
        http = FlakyHttp([ExternalServiceError("HTTP error: 429")] * 4)
        with patch.object(portal, "sleep"), self.assertRaises(ExternalServiceError):
            portal._pubtator_search_count("APLN", cache=self.cache, http=http)
        self.assertEqual(http.calls, 4)


class PubTatorLiteratureInfoTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        portal._PUBTATOR_COUNT_CACHE.clear()
        self.addCleanup(portal._PUBTATOR_COUNT_CACHE.clear)

    def _info(self, counts, previous, *, human=True):
        def count(symbol, **_):
            value = counts[symbol]
            if isinstance(value, Exception):
                raise value
            return value

        with patch.object(portal, "_pubtator_search_count", side_effect=count), patch.object(
            portal, "_hgnc_previous_symbols", return_value=previous
        ):
            return portal._get_pubtator_literature_info("ACKR5", cache_dir=Path(self.tmp.name), human=human)

    def test_zero_falls_back_to_hgnc_previous_symbol(self):
        info = self._info({"ACKR5": 0, "ADMR": 0, "GPR182": 760}, ["ADMR", "GPR182"])
        self.assertEqual(info["count"], 760)
        self.assertIn("GENE_GPR182", info["query_url"])

    def test_report_link_uses_the_underscore_gene_id(self):
        portal._PUBTATOR_COUNT_CACHE.clear()
        with patch.object(portal, "_pubtator_search_count", return_value=542) as count:
            info = portal._get_pubtator_literature_info("ERVFRD-1", cache_dir=Path(self.tmp.name), human=True)
        count.assert_called_once()
        self.assertIn("GENE_ERVFRD_1", info["query_url"])

    def test_no_fallback_outside_human(self):
        info = self._info({"ACKR5": 0, "GPR182": 760}, ["GPR182"], human=False)
        self.assertEqual(info["count"], 0)
        self.assertIn("GENE_ACKR5", info["query_url"])

    def test_failure_is_none_not_zero_and_is_reported(self):
        with patch("sys.stderr") as stderr:
            info = self._info({"ACKR5": ExternalServiceError("HTTP error: 429")}, [])
        self.assertIsNone(info["count"])
        self.assertIn("PubTator count unavailable for ACKR5", "".join(c.args[0] for c in stderr.write.call_args_list))


class PubTatorReleaseGateTests(unittest.TestCase):
    def test_missing_count_fails_the_release(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "reports").mkdir()
            (root / "reports/a.html").write_text("", encoding="utf-8")
            rows = [{"status": "ok", "detail_page": "reports/a.html", "pubtator_count": None}]
            with self.assertRaisesRegex(ValueError, "1 targets without a PubTator count"):
                _validate_catalog(
                    rows, label="human", catalog_root=root, metadata_path=root / "m.json",
                    minimum_targets=0, strict=True,
                )


if __name__ == "__main__":
    unittest.main()
