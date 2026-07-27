from __future__ import annotations

import hashlib
from pathlib import Path


class FileCache:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)

    def path_for(self, namespace: str, key: str, suffix: str) -> Path:
        digest = hashlib.sha1(key.encode("utf-8"), usedforsecurity=False).hexdigest()
        folder = self.root / namespace
        folder.mkdir(parents=True, exist_ok=True)
        return folder / f"{digest}{suffix}"
