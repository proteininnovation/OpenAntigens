"""Hardening regressions for cache poisoning, retry scope, and path containment."""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.cache import FileCache
from agdesign2.exceptions import ExternalServiceError
from agdesign2.http import HttpClient
from agdesign2.dynamic_portal import _is_within_directory
from agdesign2.refseq import _is_retryable_external_error


class HttpCachePoisoningTests(unittest.TestCase):
    def test_non_http_url_is_rejected(self) -> None:
        client = HttpClient()
        with self.assertRaisesRegex(ExternalServiceError, "Unsupported URL"):
            client.fetch_text("file:///etc/passwd")

    def test_non_json_200_body_is_not_cached(self) -> None:
        with TemporaryDirectory() as tmp:
            client = HttpClient(cache=FileCache(Path(tmp)))
            calls = {"n": 0}

            def fake_request(url, *, headers=None):
                calls["n"] += 1
                return b"<html>maintenance</html>"

            client._request = fake_request  # type: ignore[assignment]
            with self.assertRaises(json.JSONDecodeError):
                client.fetch_json("https://example.org/x")
            # The bad body must not have been written to the cache, so a second
            # call re-requests rather than re-raising from a poisoned file.
            with self.assertRaises(json.JSONDecodeError):
                client.fetch_json("https://example.org/x")
            self.assertEqual(calls["n"], 2)

    def test_valid_json_is_cached(self) -> None:
        with TemporaryDirectory() as tmp:
            client = HttpClient(cache=FileCache(Path(tmp)))
            calls = {"n": 0}

            def fake_request(url, *, headers=None):
                calls["n"] += 1
                return b'{"ok": true}'

            client._request = fake_request  # type: ignore[assignment]
            self.assertEqual(client.fetch_json("https://example.org/y"), {"ok": True})
            self.assertEqual(client.fetch_json("https://example.org/y"), {"ok": True})
            self.assertEqual(calls["n"], 1)  # second served from cache


class HttpAtomicCacheConcurrencyTests(unittest.TestCase):
    def test_concurrent_writes_to_same_cache_file_stay_valid(self) -> None:
        # Many threads fetch the same uncached URL at once; atomic writes must
        # leave a single, valid (non-interleaved) cache file. Required for the
        # parallel ortholog-table build which shares one HTTP cache.
        import json as _json
        import threading
        from concurrent.futures import ThreadPoolExecutor

        with TemporaryDirectory() as tmp:
            client = HttpClient(cache=FileCache(Path(tmp)))
            payload = _json.dumps({"value": "x" * 5000}).encode("utf-8")
            barrier = threading.Barrier(16)

            def fake_request(url, *, headers=None):
                barrier.wait()  # maximize overlap of the cache writes
                return payload

            client._request = fake_request  # type: ignore[assignment]
            with ThreadPoolExecutor(max_workers=16) as ex:
                results = list(ex.map(lambda _: client.fetch_json("https://example.org/same"), range(16)))

            self.assertTrue(all(r == {"value": "x" * 5000} for r in results))
            # The cache file must be complete, valid JSON (no torn write).
            cache_file = client._cache_path("json", "https://example.org/same", ".json")
            self.assertEqual(_json.loads(cache_file.read_text(encoding="utf-8")), {"value": "x" * 5000})


class RefseqRetryScopeTests(unittest.TestCase):
    def test_retryable_errors(self) -> None:
        for code in ("429", "500", "502", "503", "504"):
            self.assertTrue(
                _is_retryable_external_error(ExternalServiceError(f"HTTP error for u: {code}")),
                code,
            )
        self.assertTrue(
            _is_retryable_external_error(ExternalServiceError("Connection error for u: timed out"))
        )

    def test_non_retryable_errors(self) -> None:
        for code in ("400", "404", "401"):
            self.assertFalse(
                _is_retryable_external_error(ExternalServiceError(f"HTTP error for u: {code}")),
                code,
            )


class PathContainmentTests(unittest.TestCase):
    def test_sibling_prefix_directory_is_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            served = root / "run"
            sibling = root / "run-secrets"
            served.mkdir()
            sibling.mkdir()
            secret = sibling / "secret.txt"
            secret.write_text("x", encoding="utf-8")
            # String startswith("/.../run") would accept this; is_relative_to must not.
            self.assertFalse(_is_within_directory(secret, served))
            self.assertTrue(_is_within_directory(served / "a" / "b.txt", served))

    def test_parent_traversal_is_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            served = root / "run"
            served.mkdir()
            self.assertFalse(_is_within_directory((served / ".." / "etc").resolve(), served))


if __name__ == "__main__":
    unittest.main()
