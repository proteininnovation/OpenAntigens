from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path
from unittest import mock
from urllib.parse import parse_qs, urlparse

from agdesign2.uniprot_prefetch import _fetch_page_with_retry, prefetch_targets_from_tsv


class UniProtPrefetchTests(unittest.TestCase):
    def test_target_prefetch_respects_live_uniprot_query_limits(self) -> None:
        urls: list[str] = []

        class FakeHttp:
            def __init__(self, **kwargs) -> None:
                pass

            def fetch_json(self, url, **kwargs):
                urls.append(url)
                return {"results": []}

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            target_tsv = root / "targets.tsv"
            with target_tsv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["uniprot_id"], delimiter="\t")
                writer.writeheader()
                writer.writerows({"uniprot_id": f"P{index:05d}"} for index in range(101))
            with mock.patch("agdesign2.uniprot_prefetch.HttpClient", FakeHttp):
                summary = prefetch_targets_from_tsv(
                    tsv_path=target_tsv,
                    cache_dir=root / "cache",
                    chunk_size=500,
                )

        self.assertEqual(summary["chunks"], 2)
        self.assertEqual(summary["requested_identifiers"], 101)
        for url in urls:
            params = parse_qs(urlparse(url).query)
            self.assertLessEqual(params["query"][0].count(" OR ") + 1, 100)
            self.assertLessEqual(int(params["size"][0]), 100)

    def test_exhausted_retries_are_not_reported_as_success(self) -> None:
        class FailingHttp:
            def fetch_json(self, *args, **kwargs):
                raise OSError("network down")

        with self.assertRaisesRegex(RuntimeError, "failed after 2 attempts"):
            _fetch_page_with_retry(
                http=FailingHttp(),
                url="https://example.test/uniprot",
                cache_namespace="test",
                cache_key="test",
                attempts=2,
            )


if __name__ == "__main__":
    unittest.main()
