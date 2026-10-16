"""PR-9: the mouse portal's own derivations.

* **SF3** (207 mouse pages) — several human genes can claim one mouse
  accession, and the first in file order won. Pcdhgb5's page named PCDHGB3 as
  its origin while its own Homology table named PCDHGB5, and 53 human targets
  whose symbol exactly matches the mouse gene were skipped as duplicates.
* **SF7** (188) — a build-time network failure shipped as "No homolog records."
* **F19-2** (1,815) — the cynomolgus gene symbol was guessed as the first
  all-capital token of a RefSeq FASTA header, which is "X1".
* **F19-3** (266) — TrEMBL entries have no recommendedName, so the protein name
  fell back to the entry id.
* **SF4b** (1,563) — with no design region the mouse live path left identity
  and coverage None, where the human path computes a full-length number.
* **SF6** (3,925 accessions) — unversioned RefSeq accessions on mouse pages,
  versioned on human pages, for the same protein.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.mouse_portal import MOUSE_PIPELINE_VERSION, _select_mouse_claim_winners  # noqa: E402


def _claim(index: int, human_symbol: str, mouse_symbol: str, key: str = "Q0001"):
    human_report = {"target": {"gene_symbol": human_symbol, "entry_name": f"{human_symbol}_HUMAN"}}
    mouse_target = {"gene_symbol": mouse_symbol, "accession": key}
    return (index, {"query": human_symbol}, human_report, mouse_target, key, key, key)


class ClaimWinnerTests(unittest.TestCase):
    """SF3."""

    def test_the_same_symbol_wins_regardless_of_file_order(self) -> None:
        claims = [_claim(1, "PCDHGB3", "Pcdhgb5"), _claim(2, "PCDHGB5", "Pcdhgb5")]
        winner = _select_mouse_claim_winners(claims, {})["Q0001"]
        self.assertEqual(winner[0], 2)

    def test_higher_table_identity_wins_when_no_symbol_matches(self) -> None:
        claims = [_claim(1, "ALPHA", "Gamma"), _claim(2, "BETA", "Gamma")]
        rows = {
            "ALPHA_HUMAN": {"mouse_identity_to_human": "51.0"},
            "BETA_HUMAN": {"mouse_identity_to_human": "88.0"},
        }
        winner = _select_mouse_claim_winners(claims, rows)["Q0001"]
        self.assertEqual(winner[0], 2)

    def test_file_order_is_only_the_last_resort(self) -> None:
        claims = [_claim(1, "ALPHA", "Gamma"), _claim(2, "BETA", "Gamma")]
        winner = _select_mouse_claim_winners(claims, {})["Q0001"]
        self.assertEqual(winner[0], 1)

    def test_unclaimed_accessions_are_unaffected(self) -> None:
        claims = [_claim(1, "ALPHA", "Alpha", key="Q1"), _claim(2, "BETA", "Beta", key="Q2")]
        winners = _select_mouse_claim_winners(claims, {})
        self.assertEqual({k: v[0] for k, v in winners.items()}, {"Q1": 1, "Q2": 2})


class RefSeqGeneSymbolTests(unittest.TestCase):
    """F19-2."""

    def test_no_symbol_is_guessed_from_the_fasta_header(self) -> None:
        from agdesign2 import portal

        with TemporaryDirectory() as tmp:
            batch_dir = Path(tmp)
            url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein"
                "&id=XP_005563650&rettype=fasta&retmode=text"
            )
            path = portal._cached_response_path("refseq_fasta", url, ".fasta", batch_dir=batch_dir)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(
                ">XP_005563650.1 interleukin-7 isoform X1 [Macaca fascicularis]\nMKVLAA\n",
                encoding="utf-8",
            )
            result = portal._fetch_refseq_target("XP_005563650", batch_dir=batch_dir)
        self.assertIsNotNone(result)
        self.assertNotEqual(result["gene_symbol"], "X1")
        self.assertIsNone(result["gene_symbol"])


class SubmittedNameTests(unittest.TestCase):
    """F19-3."""

    def test_a_trembl_entry_gets_its_submitted_name(self) -> None:
        from agdesign2.uniprot import UniProtClient

        client = UniProtClient.__new__(UniProtClient)
        entry = {
            "primaryAccession": "Q8VF03",
            "uniProtkbId": "Q8VF03_MOUSE",
            "proteinDescription": {"submissionNames": [{"fullName": {"value": "Olfactory receptor"}}]},
            "sequence": {"value": "M" * 10},
            "organism": {"scientificName": "Mus musculus", "taxonId": 10090},
        }
        target = client._target_from_entry(entry)
        self.assertEqual(target.protein_name, "Olfactory receptor")

    def test_the_recommended_name_still_wins(self) -> None:
        from agdesign2.uniprot import UniProtClient

        client = UniProtClient.__new__(UniProtClient)
        entry = {
            "primaryAccession": "P1",
            "uniProtkbId": "X_MOUSE",
            "proteinDescription": {
                "recommendedName": {"fullName": {"value": "Interleukin-7"}},
                "submissionNames": [{"fullName": {"value": "Something else"}}],
            },
            "sequence": {"value": "M" * 10},
            "organism": {"scientificName": "Mus musculus", "taxonId": 10090},
        }
        self.assertEqual(client._target_from_entry(entry).protein_name, "Interleukin-7")


class HomologLookupFailureTests(unittest.TestCase):
    """SF7."""

    def _render(self, notes: list[dict]) -> str:
        import json

        from agdesign2.portal import _build_index_entry, render_detail_page

        report = json.loads(
            (ROOT / "tests" / "golden" / "fixtures" / "sfta3_human_report.json").read_text(encoding="utf-8")
        )
        report["ectodomain_homology"] = []
        report["notes"] = notes
        entry = _build_index_entry({"json_report": "x.json", "status": "ok"}, report, None)
        with TemporaryDirectory() as tmp:
            return render_detail_page(entry, report, batch_dir=Path(tmp), sibling_link=None,
                                      portal_title="OpenAntigens Mouse")

    def test_a_failed_lookup_is_not_shown_as_no_records(self) -> None:
        html = self._render([
            {"severity": "warning", "source": "UniProt",
             "message": "Species homolog lookup failed: Connection error"},
        ])
        self.assertIn("The species homolog lookup failed when this report was built", html)

    def test_a_genuine_empty_result_still_says_no_records(self) -> None:
        html = self._render([])
        self.assertIn("No homolog records.", html)
        self.assertNotIn("lookup failed when this report was built", html)


class VersionSentinelTests(unittest.TestCase):
    def test_the_mouse_pipeline_version_was_bumped(self) -> None:
        # Without the bump a resumed mouse build reuses every July report and
        # none of PR-1..PR-9 reaches a mouse page.
        self.assertNotEqual(MOUSE_PIPELINE_VERSION, "2026-07-no-default-af3")


if __name__ == "__main__":
    unittest.main()
