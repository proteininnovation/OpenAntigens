from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.complex_portal import (
    COMPLEX_PORTAL_HUMAN_COMPLEXTAB,
    COMPLEX_PORTAL_HUMAN_PREDICTED_COMPLEXTAB,
    COMPLEX_PORTAL_MOUSE_COMPLEXTAB,
    ComplexPortalClient,
)
from agdesign2.exceptions import ExternalServiceError


COMPLEXTAB_HEADER = (
    "#Complex ac\tRecommended name\tAliases for complex\tTaxonomy identifier\t"
    "Identifiers (and stoichiometry) of molecules in complex\tEvidence Code\t"
    "Experimental evidence\tGo Annotations\tCross references\tDescription\t"
    "Complex properties\tComplex assembly\tLigand\tDisease\tAgonist\tAntagonist\t"
    "Comment\tSource\tExpanded participant list\n"
)


class FakeHttp:
    def __init__(self) -> None:
        self.fetch_text_calls: list[str] = []
        self.taxon_id = "9606"

    def fetch_text(self, url: str, *, cache_namespace: str, suffix: str) -> str:
        self.fetch_text_calls.append(url)
        if url == COMPLEX_PORTAL_HUMAN_COMPLEXTAB:
            return COMPLEXTAB_HEADER + (
                f"CPX-TEST\tIntegrin complex\t-\t{self.taxon_id}\tCHEBI:29105(2)|P08648(1)|P05556(1)\t"
                "ECO:0000353(test evidence)\t-\t-\t-\t-\tstable|nuclear\tHeterodimer\t"
                "-\t-\t-\t-\t-\tpsi-mi:MI:0469(IntAct)\tP08648(1)|P05556(1)\n"
            )
        if url == COMPLEX_PORTAL_HUMAN_PREDICTED_COMPLEXTAB:
            return COMPLEXTAB_HEADER + (
                "CPX-PRED\tPredicted complex\t-\t9606\tP08648(0)|O11111(0)\t"
                "ECO:0008004(machine learning evidence)\t-\t-\t-\t-\t-\t-\t"
                "-\t-\t-\t-\t-\tpsi-mi:MI:2424(HuMap)\tP08648(0)|O11111(0)\n"
            )
        if url == COMPLEX_PORTAL_MOUSE_COMPLEXTAB:
            return COMPLEXTAB_HEADER + (
                "CPX-MOUSE\tMouse integrin complex\t-\t10090\tP12345(1)|Q12345(1)\t"
                "ECO:0000353(test evidence)\t-\t-\t-\t-\t-\tHeterodimer\t"
                "-\t-\t-\t-\t-\tpsi-mi:MI:0469(IntAct)\tP12345(1)|Q12345(1)\n"
            )
        raise AssertionError(f"unexpected URL: {url}")


class FailingHttp:
    def __init__(self) -> None:
        self.fetch_text_calls = 0

    def fetch_text(self, url: str, *, cache_namespace: str, suffix: str) -> str:
        self.fetch_text_calls += 1
        raise RuntimeError("catalogue unavailable")


class MalformedHttp:
    def __init__(self) -> None:
        self.fetch_text_calls = 0

    def fetch_text(self, url: str, *, cache_namespace: str, suffix: str) -> str:
        self.fetch_text_calls += 1
        return "<html>maintenance</html>"


class EmptyParticipantsHttp:
    def __init__(self) -> None:
        self.fetch_text_calls = 0

    def fetch_text(self, url: str, *, cache_namespace: str, suffix: str) -> str:
        self.fetch_text_calls += 1
        return COMPLEXTAB_HEADER + "\t".join(["CPX-EMPTY", "-", "-", "9606"] + ["-"] * 15)


class ComplexPortalClientTests(unittest.TestCase):
    def test_human_complextab_is_indexed_once_and_parses_catalogue_fields(self) -> None:
        http = FakeHttp()
        client = ComplexPortalClient(http)
        records = client.fetch_complexes_for_target(
            accession="P08648",
            gene_symbol="ITGA5",
            taxon_id=9606,
        )
        self.assertEqual([record.complex_ac for record in records], ["CPX-TEST"])
        self.assertEqual([participant.identifier for participant in records[0].participants], ["P08648", "P05556"])
        self.assertEqual(records[0].participants[1].identifier, "P05556")
        self.assertEqual(records[0].participants[0].stoichiometry, "1")
        self.assertEqual(records[0].participants[0].interactor_type, "protein")
        self.assertEqual(records[0].species, "Homo sapiens; 9606")
        self.assertEqual(records[0].evidence_code, "ECO:0000353")
        self.assertEqual(records[0].evidence_description, "test evidence")
        self.assertEqual(records[0].properties, ["stable", "nuclear"])
        self.assertEqual(records[0].complex_assemblies, ["Heterodimer"])
        self.assertFalse(records[0].predicted_complex)
        self.assertEqual(
            client.fetch_complexes_for_target(
                accession="P00000",
                gene_symbol="MISSING",
                taxon_id=9606,
            ),
            [],
        )
        self.assertEqual(http.fetch_text_calls, [COMPLEX_PORTAL_HUMAN_COMPLEXTAB])

    def test_predicted_catalogue_is_explicit_opt_in(self) -> None:
        http = FakeHttp()
        client = ComplexPortalClient(http)
        records = client.fetch_human_complexes_for_target(accession="P08648", include_predicted=True)

        self.assertEqual([record.complex_ac for record in records], ["CPX-TEST", "CPX-PRED"])
        self.assertTrue(records[1].predicted_complex)
        self.assertEqual(
            http.fetch_text_calls,
            [COMPLEX_PORTAL_HUMAN_COMPLEXTAB, COMPLEX_PORTAL_HUMAN_PREDICTED_COMPLEXTAB],
        )

    def test_mouse_lookup_uses_the_cached_mouse_catalogue(self) -> None:
        http = FakeHttp()
        records = ComplexPortalClient(http).fetch_complexes_for_target(
            accession="P12345",
            gene_symbol="Itga5",
            taxon_id=10090,
        )

        self.assertEqual([record.complex_ac for record in records], ["CPX-MOUSE"])
        self.assertEqual(records[0].species, "Mus musculus; 10090")
        self.assertEqual(http.fetch_text_calls, [COMPLEX_PORTAL_MOUSE_COMPLEXTAB])

    def test_human_catalogue_failure_is_re_raised_without_retries(self) -> None:
        http = FailingHttp()
        client = ComplexPortalClient(http)

        with self.assertRaisesRegex(RuntimeError, "catalogue unavailable"):
            client.fetch_human_complexes_for_target(accession="P08648")
        with self.assertRaisesRegex(RuntimeError, "catalogue unavailable"):
            client.fetch_human_complexes_for_target(accession="P08648")

        self.assertEqual(http.fetch_text_calls, 1)

    def test_malformed_human_catalogue_is_not_treated_as_no_hits(self) -> None:
        http = MalformedHttp()
        client = ComplexPortalClient(http)

        with self.assertRaisesRegex(ExternalServiceError, "Unexpected ComplexTab header"):
            client.fetch_human_complexes_for_target(accession="P08648")
        with self.assertRaisesRegex(ExternalServiceError, "Unexpected ComplexTab header"):
            client.fetch_human_complexes_for_target(accession="P08648")

        self.assertEqual(http.fetch_text_calls, 1)

    def test_human_catalogue_rejects_nonhuman_rows(self) -> None:
        http = FakeHttp()
        http.taxon_id = "10090"

        with self.assertRaisesRegex(ExternalServiceError, "Unexpected ComplexTab taxonomy"):
            ComplexPortalClient(http).fetch_human_complexes_for_target(accession="P08648")

    def test_human_catalogue_with_no_participants_is_not_treated_as_no_hits(self) -> None:
        http = EmptyParticipantsHttp()
        client = ComplexPortalClient(http)

        with self.assertRaisesRegex(ExternalServiceError, "contains no participants"):
            client.fetch_human_complexes_for_target(accession="P08648")
        with self.assertRaisesRegex(ExternalServiceError, "contains no participants"):
            client.fetch_human_complexes_for_target(accession="P08648")

        self.assertEqual(http.fetch_text_calls, 1)

    def test_participant_match_does_not_accept_requested_gene_for_an_unrelated_participant(self) -> None:
        client = ComplexPortalClient.__new__(ComplexPortalClient)
        self.assertFalse(
            client._participant_matches_target(
                identifier="P99999",
                name="UNRELATED",
                accession="P08648",
                gene_symbol="ITGA5",
            )
        )
        self.assertTrue(
            client._participant_matches_target(
                identifier="P08648",
                name="ITGA5",
                accession="P08648",
                gene_symbol="ITGA5",
            )
        )


if __name__ == "__main__":
    unittest.main()
