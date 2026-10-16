"""F04b: HTTP 204 from the InterPro API is "no entries", not a failed lookup.

``InterProClient.fetch_annotations`` queries two databases in one loop. The API
answers 204 with an empty body when a protein has no entries in one of them;
``json.loads("")`` raised, the exception escaped, and ``pipeline.py:213-222``
recorded a lookup failure and discarded the entries already read from the other
database. On the 2026-09-30 build that is 157 human + 89 mouse reports, of which
a live re-probe found 179 do have real InterPro annotations.
"""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.cache import FileCache  # noqa: E402
from agdesign2.exceptions import ExternalServiceError  # noqa: E402
from agdesign2.http import HttpClient  # noqa: E402
from agdesign2.interpro import InterProClient  # noqa: E402

_INTERPRO_PAYLOAD = {
    "results": [
        {
            "metadata": {
                "accession": "IPR005461",
                "name": "A domain",
                "type": "domain",
                "source_database": "interpro",
            },
            "proteins": [{"entry_protein_locations": [{"fragments": [{"start": 10, "end": 90}]}]}],
        }
    ],
    "next": None,
}


def _responses(interpro: bytes, pfam: bytes):
    def fake(self, url, *, headers=None):  # noqa: ANN001
        return (pfam if "/pfam/" in url else interpro), {}

    return fake


class InterProEmptyBodyTests(unittest.TestCase):
    def test_empty_pfam_body_keeps_the_interpro_entries(self) -> None:
        client = InterProClient(HttpClient(cache=None))
        with patch.object(
            HttpClient,
            "_request_response",
            _responses(json.dumps(_INTERPRO_PAYLOAD).encode(), b""),
        ):
            annotations = client.fetch_annotations("Q9HBG7")
        self.assertEqual([a.accession for a in annotations], ["IPR005461"])

    def test_empty_body_is_cached_as_an_object_not_as_empty_bytes(self) -> None:
        # A cached "" would raise on every later read (poisoned cache).
        with TemporaryDirectory() as tmp:
            http = HttpClient(cache=FileCache(Path(tmp)))
            client = InterProClient(http)
            with patch.object(
                HttpClient,
                "_request_response",
                _responses(json.dumps(_INTERPRO_PAYLOAD).encode(), b""),
            ):
                client.fetch_annotations("Q9HBG7")
            # Second call must be served from cache without touching the network.
            with patch.object(
                HttpClient,
                "_request_response",
                lambda *a, **k: (_ for _ in ()).throw(AssertionError("network used")),
            ):
                annotations = client.fetch_annotations("Q9HBG7")
        self.assertEqual([a.accession for a in annotations], ["IPR005461"])

    def test_connection_errors_still_raise(self) -> None:
        client = InterProClient(HttpClient(cache=None))

        def boom(self, url, *, headers=None):  # noqa: ANN001
            raise ExternalServiceError("Connection error for %s: offline" % url)

        with patch.object(HttpClient, "_request_response", boom):
            with self.assertRaises(ExternalServiceError):
                client.fetch_annotations("Q9HBG7")

    def test_non_empty_non_json_body_still_raises(self) -> None:
        client = InterProClient(HttpClient(cache=None))
        with patch.object(
            HttpClient, "_request_response", _responses(b"<html>maintenance</html>", b"")
        ):
            with self.assertRaises(json.JSONDecodeError):
                client.fetch_annotations("Q9HBG7")


if __name__ == "__main__":
    unittest.main()
