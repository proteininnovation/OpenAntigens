from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.blast import BlastClient
from agdesign2.config import AnalysisConfig


class FakeHttp:
    def download(self, url: str, destination: Path, **_: object) -> Path:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text("", encoding="utf-8")
        return destination


class StubBlastClient(BlastClient):
    def __init__(self, responses: dict[tuple[str, str], list], config: AnalysisConfig | None = None) -> None:
        super().__init__(FakeHttp(), config or AnalysisConfig())
        self.responses = responses

    def ensure_databases(self) -> dict[str, Path]:
        return {species: Path(f"/tmp/{species}") for species in self.config.blast_species}

    def ensure_ortholog_databases(self) -> dict[str, Path]:
        return {species: Path(f"/tmp/{species}_ortholog") for species in self.config.ortholog_blast_species}

    def _search_database(
        self,
        *,
        query_sequence: str,
        query_path: Path,
        db_prefix: Path,
        species: str,
        query_length: int,
        target_accession: str,
        target_entry_name: str,
        output_path: Path,
    ):
        namespace = "ortholog" if db_prefix.name.endswith("_ortholog") else "reviewed"
        return list(self.responses.get((species, namespace), []))


class StubBulkBlastClient(BlastClient):
    def __init__(self, responses: dict[tuple[str, str], dict[str, list]], config: AnalysisConfig | None = None) -> None:
        super().__init__(FakeHttp(), config or AnalysisConfig(blast_species=("human", "mouse")))
        self.responses = responses
        self.batch_calls: list[tuple[str, list[str]]] = []

    def ensure_databases(self) -> dict[str, Path]:
        return {species: Path(f"/tmp/{species}") for species in self.config.blast_species}

    def ensure_ortholog_databases(self) -> dict[str, Path]:
        return {species: Path(f"/tmp/{species}_ortholog") for species in self.config.ortholog_blast_species}

    def _search_database_many(
        self,
        *,
        query_path: Path,
        db_prefix: Path,
        species: str,
        query_lengths: dict[str, int],
        target_accessions: dict[str, str],
        target_entry_names: dict[str, str],
        output_path: Path,
    ):
        namespace = "ortholog" if db_prefix.name.endswith("_ortholog") else "reviewed"
        self.batch_calls.append((species, sorted(query_lengths)))
        return {key: list(value) for key, value in self.responses.get((species, namespace), {}).items()}


class BlastTests(unittest.TestCase):
    def test_cli_compatible_path_uses_private_user_alias_directory(self) -> None:
        client = BlastClient(FakeHttp(), AnalysisConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source_dir = root / "path with spaces"
            source_dir.mkdir()
            source = source_dir / "mouse"
            source.write_text("", encoding="utf-8")
            with patch("agdesign2.blast.tempfile.gettempdir", return_value=tmpdir):
                alias = client._cli_compatible_path(source)

            alias_root = root / f"agdesign2_blast_aliases_{os.getuid()}"
            self.assertTrue(alias.parent.is_symlink())
            self.assertEqual(alias.resolve(), source.resolve())
            self.assertEqual(alias_root.stat().st_mode & 0o777, 0o700)

    def test_parses_tabular_hits_and_filters_self(self) -> None:
        client = BlastClient(FakeHttp(), AnalysisConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "hits.tsv"
            path.write_text(
                "\n".join(
                    [
                        "sp|P00533|EGFR_HUMAN\t99.0\t100\t1\t100\t1\t100\t0.0\t400\tAAAAAAAAAA\tAAAAAAAAAA\tsp|P00533|EGFR_HUMAN Epidermal growth factor receptor OS=Homo sapiens OX=9606",
                        "sp|Q01279|EGFR_MOUSE\t88.0\t90\t2\t91\t5\t94\t1e-30\t250\tACDEFGHIKL\tACD-FGHIKL\tsp|Q01279|EGFR_MOUSE Epidermal growth factor receptor OS=Mus musculus OX=10090",
                    ]
                ),
                encoding="utf-8",
            )
            hits = client._parse_hits(
                path,
                species="mouse",
                query_length=100,
                target_accession="P00533",
                target_entry_name="EGFR_HUMAN",
            )
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].species, "Mus musculus")
        self.assertEqual(hits[0].query_alignment, "ACDEFGHIKL")
        self.assertEqual(hits[0].subject_alignment, "ACD-FGHIKL")
        self.assertEqual(hits[0].alignment_source, "BLAST")

    def test_coverage_uses_query_span_and_does_not_exceed_100_percent(self) -> None:
        client = BlastClient(FakeHttp(), AnalysisConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "hits.tsv"
            path.write_text(
                "sp|Q01279|EGFR_MOUSE\t88.0\t120\t1\t100\t5\t124\t1e-30\t250\tAAAAAAAAAA\tAAAAAAAAAA\tsp|Q01279|EGFR_MOUSE Epidermal growth factor receptor OS=Mus musculus OX=10090",
                encoding="utf-8",
            )
            hits = client._parse_hits(
                path,
                species="mouse",
                query_length=100,
                target_accession="P00533",
                target_entry_name="EGFR_HUMAN",
            )
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].coverage, 100.0)

    def test_cross_reactivity_falls_back_to_macaque_proteome(self) -> None:
        client = StubBlastClient(
            responses={
                ("mouse", "reviewed"): [],
                ("human", "reviewed"): [],
                ("macaca_fascicularis", "reviewed"): [],
                (
                    "macaca_fascicularis",
                    "ortholog",
                ): [
                    client_hit(
                        subject_id="tr|A0A2K5WK39|A0A2K5WK39_MACFA",
                        species="Macaca fascicularis",
                        identity=98.7,
                    )
                ],
            }
        )
        hits = client.search("AAAA", target_accession="P00533", target_entry_name="EGFR_HUMAN")
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].species, "Macaca fascicularis")
        self.assertIn("MACFA", hits[0].subject_id)

    def test_dedupes_hits_by_subject(self) -> None:
        client = BlastClient(FakeHttp(), AnalysisConfig())
        hits = client._dedupe_hits(
            [
                client_hit(subject_id="sp|P08069|IGF1R_HUMAN", species="Homo sapiens", bitscore=200.0),
                client_hit(subject_id="sp|P08069|IGF1R_HUMAN", species="Homo sapiens", bitscore=250.0),
            ]
        )
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0].bitscore, 250.0)

    def test_search_many_batches_queries_by_species(self) -> None:
        client = StubBulkBlastClient(
            responses={
                ("human", "reviewed"): {
                    "q1": [client_hit(subject_id="sp|P08069|IGF1R_HUMAN", species="Homo sapiens", bitscore=180.0)]
                },
                ("mouse", "reviewed"): {
                    "q2": [client_hit(subject_id="sp|Q01279|EGFR_MOUSE", species="Mus musculus", bitscore=220.0)]
                },
            },
        )
        hits = client.search_many(
            [
                {"id": "q1", "sequence": "AAAA", "target_accession": "P00533", "target_entry_name": "EGFR_HUMAN"},
                {"id": "q2", "sequence": "BBBB", "target_accession": "Q99999", "target_entry_name": "TEST_HUMAN"},
            ]
        )
        self.assertEqual(set(hits), {"q1", "q2"})
        self.assertEqual(hits["q1"][0].subject_id, "sp|P08069|IGF1R_HUMAN")
        self.assertEqual(hits["q2"][0].subject_id, "sp|Q01279|EGFR_MOUSE")
        self.assertEqual([call[0] for call in client.batch_calls], ["human", "mouse"])

    def test_blast_tool_uses_explicit_bin_dir_when_path_is_missing(self) -> None:
        client = BlastClient(FakeHttp(), AnalysisConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            tool = Path(tmpdir) / "blastp"
            tool.write_text("#!/bin/sh\n", encoding="utf-8")
            tool.chmod(0o755)
            with patch.dict("os.environ", {"AGDESIGN2_BLAST_BIN_DIR": tmpdir, "PATH": ""}, clear=True):
                self.assertEqual(client._blast_tool("blastp"), str(tool))

    def test_blast_tool_reports_missing_tool(self) -> None:
        from agdesign2.exceptions import BlastDatabaseError

        client = BlastClient(FakeHttp(), AnalysisConfig())
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.dict("os.environ", {"AGDESIGN2_BLAST_BIN_DIR": tmpdir, "PATH": ""}, clear=True):
                with self.assertRaisesRegex(BlastDatabaseError, "definitely-not-blast was not found"):
                    client._blast_tool("definitely-not-blast")


def client_hit(
    *,
    subject_id: str,
    species: str,
    identity: float = 90.0,
    coverage: float = 100.0,
    bitscore: float = 300.0,
):
    from agdesign2.models import BlastHit

    return BlastHit(
        subject_id=subject_id,
        description=subject_id,
        species=species,
        identity=identity,
        coverage=coverage,
        alignment_length=100,
        evalue=1e-20,
        bitscore=bitscore,
        query_start=1,
        query_end=100,
        subject_start=1,
        subject_end=100,
    )


if __name__ == "__main__":
    unittest.main()
