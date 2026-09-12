from __future__ import annotations

import json
import hashlib
import re
import os
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

from .cache import FileCache
from .exceptions import ExternalServiceError


def _atomic_write(path: Path, data: bytes) -> None:
    """Write ``data`` to ``path`` atomically (temp file in the same dir + replace).

    os.replace is atomic on POSIX, so a reader (or a concurrent writer in another
    thread) never observes a half-written cache file. This makes the file cache
    safe to share across threads — required by the parallel ortholog-table build.
    The temp name includes the pid+thread id so concurrent writers to the same
    final path do not clobber each other's temp file.
    """
    tmp = path.with_name(f"{path.name}.{os.getpid()}.{_thread_id()}.tmp")
    try:
        tmp.write_bytes(data)
        os.replace(tmp, path)
    finally:
        try:
            tmp.unlink()
        except OSError:
            pass


def _thread_id() -> int:
    import threading

    return threading.get_ident()


class HttpClient:
    def __init__(self, cache: FileCache | None = None, timeout: int = 60) -> None:
        self.cache = cache
        self.timeout = timeout

    def fetch_json(
        self,
        url: str,
        *,
        headers: dict[str, str] | None = None,
        cache_namespace: str = "json",
        cache_key: str | None = None,
    ) -> Any:
        cache_path = self._cache_path(cache_namespace, cache_key or url, ".json")
        if cache_path and cache_path.exists():
            return json.loads(cache_path.read_text(encoding="utf-8"))
        payload = self._request(url, headers=headers)
        text = payload.decode("utf-8")
        # Parse before caching: a 200 response with a non-JSON body (e.g. an
        # HTML maintenance page) must not be written to the cache, or every
        # later read would re-raise JSONDecodeError forever (poisoned cache).
        value = json.loads(text)
        if cache_path:
            _atomic_write(cache_path, text.encode("utf-8"))
        return value

    def fetch_text(
        self,
        url: str,
        *,
        headers: dict[str, str] | None = None,
        cache_namespace: str = "text",
        cache_key: str | None = None,
        suffix: str = ".txt",
    ) -> str:
        cache_path = self._cache_path(cache_namespace, cache_key or url, suffix)
        if cache_path and cache_path.exists():
            return cache_path.read_text(encoding="utf-8")
        payload = self._request(url, headers=headers)
        text = payload.decode("utf-8")
        if cache_path:
            _atomic_write(cache_path, text.encode("utf-8"))
        return text

    def download(
        self,
        url: str,
        destination: Path,
        *,
        headers: dict[str, str] | None = None,
    ) -> Path:
        receipt = destination.with_name(destination.name + ".download.json")
        if destination.exists() and receipt.exists():
            recorded = json.loads(receipt.read_text())
            if recorded.get("url") == url and recorded.get("sha256") == hashlib.sha256(destination.read_bytes()).hexdigest():
                return destination
        destination.parent.mkdir(parents=True, exist_ok=True)
        payload = self._request(url, headers=headers)
        _atomic_write(destination, payload)
        _atomic_write(receipt, json.dumps({"url": url, "sha256": hashlib.sha256(payload).hexdigest()}).encode())
        return destination

    def fetch_json_page(self, url: str) -> tuple[Any, str | None]:
        # Cursor pages are fetched live; caching them can outlive the server cursor.
        payload, headers = self._request_response(url)
        value = json.loads(payload)
        next_link = re.search(r'<([^>]+)>;\s*rel="next"', headers.get("Link", ""))
        return value, next_link.group(1) if next_link else None

    def _cache_path(self, namespace: str, key: str, suffix: str) -> Path | None:
        if not self.cache:
            return None
        return self.cache.path_for(namespace, key, suffix)

    def _request(self, url: str, *, headers: dict[str, str] | None = None) -> bytes:
        return self._request_response(url, headers=headers)[0]

    def _request_response(self, url: str, *, headers: dict[str, str] | None = None) -> tuple[bytes, Any]:
        parsed = urlsplit(url)
        if parsed.scheme not in {"http", "https"} or not parsed.netloc:
            raise ExternalServiceError(f"Unsupported URL: {url}")
        request = Request(
            url,
            headers={
                "User-Agent": "agdesign2/0.1",
                "Accept": "application/json",
                **(headers or {}),
            },
        )
        try:
            with urlopen(request, timeout=self.timeout) as response:  # nosec B310
                return response.read(), response.headers
        except HTTPError as exc:
            raise ExternalServiceError(f"HTTP error for {url}: {exc.code}") from exc
        except URLError as exc:
            raise ExternalServiceError(f"Connection error for {url}: {exc.reason}") from exc
