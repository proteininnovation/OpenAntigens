from __future__ import annotations

import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HOOK = ROOT / ".githooks" / "pre-push"
CANONICAL = "https://github.com/proteininnovation/OpenAntigens.git"


class PrePushHookTests(unittest.TestCase):
    @staticmethod
    def _run(local_ref: str, remote_ref: str, url: str = CANONICAL) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [str(HOOK), "origin", url],
            input=f"{local_ref} local-oid {remote_ref} remote-oid\n",
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )

    def test_allows_informative_branch_on_canonical_remote(self) -> None:
        result = self._run("refs/heads/chore/test", "refs/heads/chore/test")
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_rejects_any_codex_prefix_on_either_branch(self) -> None:
        for ref in ("refs/heads/codex", "refs/heads/codex/fix", "refs/heads/codex-fix"):
            with self.subTest(ref=ref):
                result = self._run(ref, ref)
                self.assertNotEqual(result.returncode, 0)
        result = self._run("refs/heads/feat/test", "refs/heads/codex-old")
        self.assertNotEqual(result.returncode, 0)

    def test_rejects_noncanonical_remote(self) -> None:
        result = self._run("refs/heads/feat/test", "refs/heads/feat/test", "https://github.com/example/wrong.git")
        self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
