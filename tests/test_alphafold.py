from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.alphafold import AlphaFoldClient
from agdesign2.config import AnalysisConfig
from agdesign2.exceptions import ExternalServiceError


class FakeHttpClient:
    def __init__(self) -> None:
        self.downloaded: list[str] = []

    def fetch_json(self, url: str, *, headers=None, cache_namespace="json", cache_key=None):
        if url == "https://alphafold.ebi.ac.uk/api/prediction/Q9BYF1":
            raise ExternalServiceError("metadata unavailable")
        if url == "https://alphafold.com/api/prediction/Q9BYF1":
            return []
        if url == "https://www.alphafold.com/api/prediction/Q9BYF1":
            return []
        raise AssertionError(f"Unexpected metadata URL: {url}")

    def download(self, url: str, destination: Path, *, headers=None) -> Path:
        self.downloaded.append(url)
        destination.parent.mkdir(parents=True, exist_ok=True)
        if url.endswith("model_v4.pdb") or url.endswith("predicted_aligned_error_v4.json"):
            raise ExternalServiceError(f"HTTP error for {url}: 404")
        if url.endswith("model_v4.pdb.gz"):
            with gzip.open(destination, "wb") as handle:
                handle.write(b"ATOM      1  N   MET A   1      11.000  13.000  15.000  1.00 90.00           N  \n")
            return destination
        if url.endswith("predicted_aligned_error_v4.json.gz"):
            with gzip.open(destination, "wb") as handle:
                handle.write(b'{"predicted_aligned_error":[[0.0]]}')
            return destination
        raise AssertionError(f"Unexpected download URL: {url}")


class CanonicalSelectingHttpClient:
    def __init__(self) -> None:
        self.downloaded: list[str] = []

    def fetch_json(self, url: str, *, headers=None, cache_namespace="json", cache_key=None):
        if url.endswith("/api/prediction/P08195") or url.endswith("/api/prediction/P08195-2"):
            return [
                {
                    "modelEntityId": "AF-P08195-F1",
                    "uniprotAccession": "P08195",
                    "sequence": "ABCDEFGH",
                    "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P08195-F1-model_v6.pdb",
                    "paeDocUrl": "https://alphafold.ebi.ac.uk/files/AF-P08195-F1-predicted_aligned_error_v6.json",
                },
                {
                    "modelEntityId": "AF-P08195-2-F1",
                    "uniprotAccession": "P08195-2",
                    "sequence": "CDEFGH",
                    "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P08195-2-F1-model_v6.pdb",
                    "paeDocUrl": "https://alphafold.ebi.ac.uk/files/AF-P08195-2-F1-predicted_aligned_error_v6.json",
                },
            ]
        return []

    def download(self, url: str, destination: Path, *, headers=None) -> Path:
        self.downloaded.append(url)
        destination.parent.mkdir(parents=True, exist_ok=True)
        if url.endswith("AF-P08195-2-F1-model_v6.pdb"):
            destination.write_text(
                "ATOM      1  CA  CYS A   1       0.000   0.000   0.000  1.00 90.00           C  \n"
                "ATOM      2  CA  ASP A   2       0.000   0.000   0.000  1.00 90.00           C  \n"
                "ATOM      3  CA  GLU A   3       0.000   0.000   0.000  1.00 90.00           C  \n"
                "ATOM      4  CA  PHE A   4       0.000   0.000   0.000  1.00 90.00           C  \n"
                "ATOM      5  CA  GLY A   5       0.000   0.000   0.000  1.00 90.00           C  \n"
                "ATOM      6  CA  HIS A   6       0.000   0.000   0.000  1.00 90.00           C  \n",
                encoding="utf-8",
            )
            return destination
        if url.endswith("AF-P08195-2-F1-predicted_aligned_error_v6.json"):
            destination.write_text('{"predicted_aligned_error":[[0.0]]}', encoding="utf-8")
            return destination
        if url.endswith(".pdb"):
            destination.write_text(
                "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C  \n",
                encoding="utf-8",
            )
            return destination
        if url.endswith(".json"):
            destination.write_text('{"predicted_aligned_error":[[0.0]]}', encoding="utf-8")
            return destination
        raise AssertionError(f"Unexpected download URL: {url}")


class AlphaFoldClientTests(unittest.TestCase):
    def test_disabled_fetch_uses_compatible_local_artifacts_without_remote_lookup(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            local_paths = (root / "local.pdb", root / "local.pae.json")
            config = AnalysisConfig(
                data_dir=root,
                fetch_alphafold=False,
                enable_local_af3_fallback=True,
            )
            client = AlphaFoldClient(FakeHttpClient(), config)
            with (
                patch(
                    "agdesign2.alphafold.find_local_af3_artifacts",
                    return_value=local_paths,
                ),
                patch.object(
                    client,
                    "_ensure_alphafold_db_artifacts",
                ) as remote_lookup,
            ):
                self.assertEqual(
                    client.ensure_artifacts("P00001", canonical_sequence="ACD"),
                    local_paths,
                )
            remote_lookup.assert_not_called()

    def test_disabled_fetch_fails_when_no_compatible_local_artifact_exists(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            config = AnalysisConfig(
                data_dir=Path(tmpdir),
                fetch_alphafold=False,
            )
            client = AlphaFoldClient(FakeHttpClient(), config)
            with (
                patch(
                    "agdesign2.alphafold.find_local_af3_artifacts",
                    return_value=None,
                ),
                patch.object(
                    client,
                    "_ensure_alphafold_db_artifacts",
                ) as remote_lookup,
            ):
                with self.assertRaisesRegex(
                    ExternalServiceError,
                    "remote retrieval is disabled",
                ):
                    client.ensure_artifacts("P00001", canonical_sequence="ACD")
            remote_lookup.assert_not_called()

    def test_tries_multiple_hosts_and_gzip_fallbacks(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            config = AnalysisConfig(data_dir=Path(tmpdir))
            client = AlphaFoldClient(FakeHttpClient(), config)
            pdb_path, pae_path = client.ensure_artifacts("Q9BYF1")
            self.assertTrue(pdb_path.exists())
            self.assertTrue(pae_path.exists())
            self.assertIn("ATOM", pdb_path.read_text())
            self.assertIn("predicted_aligned_error", pae_path.read_text())

    def test_selects_model_matching_canonical_isoform_sequence(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            config = AnalysisConfig(data_dir=Path(tmpdir))
            http = CanonicalSelectingHttpClient()
            client = AlphaFoldClient(http, config)
            pdb_path, pae_path = client.ensure_artifacts(
                "P08195",
                canonical_sequence="CDEFGH",
                canonical_isoform_id="P08195-2",
            )
            self.assertTrue(pdb_path.exists())
            self.assertTrue(pae_path.exists())
            self.assertTrue(any(url.endswith("AF-P08195-2-F1-model_v6.pdb") for url in http.downloaded))
            self.assertFalse(any(url.endswith("AF-P08195-F1-model_v6.pdb") for url in http.downloaded))

    def test_raises_when_no_model_matches_canonical_sequence(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            config = AnalysisConfig(data_dir=Path(tmpdir))
            http = CanonicalSelectingHttpClient()
            client = AlphaFoldClient(http, config)
            with self.assertRaises(ExternalServiceError):
                client.ensure_artifacts(
                    "P08195",
                    canonical_sequence="ZZZZZZ",
                    canonical_isoform_id="P08195-2",
                )


if __name__ == "__main__":
    unittest.main()
