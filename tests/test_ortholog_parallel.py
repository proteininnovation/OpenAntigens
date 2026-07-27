"""The parallel ortholog-table build must produce byte-identical output
to the sequential build, regardless of worker completion order."""

from __future__ import annotations

import sys
import time
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2 import ortholog_table as ot

N = 25


def _fake_build_table_row(*, analyzer, resolver, batch_row, input_index):
    # Deterministic row from the index, plus a sleep that makes higher indices
    # finish FIRST under concurrency — this would corrupt output order if results
    # were appended by completion rather than placed by input index.
    time.sleep((N - input_index) * 0.002)
    data = ot._base_row_data(batch_row=batch_row, input_index=input_index)
    data["status"] = "ok"
    data["human_gene_symbol"] = f"GENE{input_index}"
    data["human_refseq_sequence"] = "M" * ((input_index % 5) + 1)
    data["mouse_identity_to_human"] = round(input_index / N, 2)
    return ot.OrthologTableRow(**data)


class _FakeAnalyzer:
    class _Config:
        ortholog_fasta_dir = Path("/tmp")
        species_taxonomy = {"human": 9606, "mouse": 10090, "macaca_fascicularis": 9541}

    class _RefSeq:
        http = None

    config = _Config()
    refseq_client = _RefSeq()


def _build(tmp: Path, jobs: int) -> tuple[bytes, bytes]:
    tsv = tmp / "in.tsv"
    tsv.write_text("uniprot_name\n" + "\n".join(f"P{i}_HUMAN" for i in range(1, N + 1)) + "\n", encoding="utf-8")
    out = tmp / "ortholog_reference_table.tsv"
    # Stub the proteome resolver so the parity test stays hermetic (no downloads).
    with mock.patch.object(ot, "_build_table_row", side_effect=_fake_build_table_row), mock.patch.object(
        ot, "OrthologProteomeResolver"
    ):
        ot.build_ortholog_table_from_tsv(
            analyzer=_FakeAnalyzer(), tsv_path=tsv, output_path=out, resume=False, jobs=jobs
        )
    return out.read_bytes(), out.with_suffix(".json").read_bytes()


class OrthologParallelParityTests(unittest.TestCase):
    def test_parallel_output_matches_sequential_byte_for_byte(self) -> None:
        with TemporaryDirectory() as a, TemporaryDirectory() as b:
            seq_tsv, seq_json = _build(Path(a), jobs=1)
            par_tsv, par_json = _build(Path(b), jobs=8)
        self.assertEqual(seq_tsv, par_tsv, "parallel TSV differs from sequential")
        self.assertEqual(seq_json, par_json, "parallel JSON differs from sequential")
        # Sanity: the rows are actually in input order in the output.
        self.assertIn(b"GENE1\t", seq_tsv)
        self.assertLess(seq_tsv.index(b"GENE1\t"), seq_tsv.index(b"GENE25\t"))


if __name__ == "__main__":
    unittest.main()
