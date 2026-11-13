"""Design regions for proteins whose UniProt record lacks a signal peptide or extracellular topology.

Before these fallbacks, 138 human pages (61 without a TM, 77 single-pass) and mouse VEGFA had no
design region, so no constructs were offered (snapshot openantigen-2026-11-01).
"""
import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from agdesign2.config import AnalysisConfig
from agdesign2.models import Feature, TargetRecord
from agdesign2.pipeline import AntigenAnalyzer
from agdesign2.topology import derive_ectodomain
from agdesign2.uniprot import membrane_hints_from_entry, secreted_chain_names


def _tm(start: int, end: int) -> Feature:
    return Feature(type="TRANSMEM", start=start, end=end, description="Helical")


def _topo(start: int, end: int, description: str) -> Feature:
    return Feature(type="TOPO_DOM", start=start, end=end, description=description)


def _location_comment(locations: list[tuple[str, str]], molecule: str | None = None) -> dict:
    comment = {
        "commentType": "SUBCELLULAR LOCATION",
        "subcellularLocations": [
            {"location": {"value": loc}, **({"topology": {"value": topo}} if topo else {})} for loc, topo in locations
        ],
    }
    if molecule:
        comment["molecule"] = molecule
    return comment


class SinglePassOrientationTests(unittest.TestCase):
    config = AnalysisConfig()

    def test_cytoplasmic_c_terminus_makes_the_n_terminus_the_ectodomain(self) -> None:
        # KCNE1: TM 44-66, only the cytoplasmic tail 67-129 is annotated.
        result = derive_ectodomain([_tm(44, 66), _topo(67, 129, "Cytoplasmic")], 129, self.config)
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (1, 43))
        self.assertIn("cytoplasmic", result.ectodomain.metadata["orientation_source"])

    def test_cytoplasmic_c_terminus_on_a_large_receptor(self) -> None:
        # SDK1: TM 2010-2030, cytoplasmic 2031-2213 -> 2,009-residue ectodomain.
        result = derive_ectodomain([_tm(2010, 2030), _topo(2031, 2213, "Cytoplasmic")], 2213, self.config)
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (1, 2009))

    def test_vesicular_side_means_an_organelle_protein_and_no_ectodomain(self) -> None:
        # VAMP8: cytoplasmic 1-75, TM 76-96, vesicular 97-100.
        features = [_topo(1, 75, "Cytoplasmic"), _tm(76, 96), _topo(97, 100, "Vesicular")]
        self.assertIsNone(derive_ectodomain(features, 100, self.config).ectodomain)

    def test_no_plasma_membrane_location_means_no_inference(self) -> None:
        # STOML1-like: cytoplasmic C-terminus annotated, but UniProt places it only in endosomes.
        features = [_tm(58, 78), _topo(79, 398, "Cytoplasmic")]
        self.assertIsNone(derive_ectodomain(features, 398, self.config, at_plasma_membrane=False).ectodomain)

    def test_membrane_type_ii_puts_the_c_terminus_outside(self) -> None:
        result = derive_ectodomain([_tm(56, 76)], 267, self.config, membrane_types=("Single-pass type II membrane protein",))
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (77, 267))

    def test_membrane_type_iii_puts_the_n_terminus_outside(self) -> None:
        result = derive_ectodomain([_tm(100, 120)], 132, self.config, membrane_types=("Single-pass type III membrane protein",))
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (1, 99))

    def test_positive_inside_rule_only_when_allowed(self) -> None:
        # Type II signal anchor near the N-terminus with a K/R-rich N-terminal flank -> C-terminal ectodomain.
        sequence = "S" * 30 + "SKRKRSSSSS" + "L" * 21 + "S" * 160
        features = [_tm(41, 61)]
        self.assertIsNone(derive_ectodomain(features, len(sequence), self.config, sequence=sequence).ectodomain)
        result = derive_ectodomain(
            features, len(sequence), self.config, sequence=sequence, allow_predicted_orientation=True
        )
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (62, len(sequence)))
        self.assertEqual(result.ectodomain.confidence, 0.4)

    def test_positive_inside_rule_refuses_implausible_insertions(self) -> None:
        # A 157-residue N-terminal "ectodomain" without a signal peptide, a type II anchor 1,000 residues
        # in (NCKAP1-like), and a mitochondrial protein are all cytoplasmic-protein artefacts.
        long_n = "S" * 157 + "L" * 21 + "KRKRKRSSSS"
        deep_anchor = "S" * 984 + "KRKRK" + "S" * 5 + "L" * 21 + "S" * 100
        for sequence, tm, extra in (
            (long_n, _tm(158, 178), []),
            (deep_anchor, _tm(995, 1015), []),
            ("SKRKR" + "S" * 35 + "L" * 21 + "S" * 100, _tm(41, 61), [Feature(type="TRANSIT_PEPTIDE", start=1, end=30)]),
        ):
            result = derive_ectodomain(
                [tm, *extra], len(sequence), self.config, sequence=sequence, allow_predicted_orientation=True
            )
            self.assertIsNone(result.ectodomain)

    def test_annotated_extracellular_topology_still_wins(self) -> None:
        features = [_topo(1, 40, "Extracellular"), _tm(41, 61), _topo(62, 100, "Cytoplasmic")]
        result = derive_ectodomain(features, 100, self.config, membrane_types=("Single-pass type II membrane protein",))
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (1, 40))
        self.assertEqual(result.ectodomain.metadata, {})


class UniProtLocationHintTests(unittest.TestCase):
    def test_membrane_hints_ignore_isoform_blocks_and_read_plasma_membrane_types(self) -> None:
        entry = {
            "comments": [
                _location_comment([("Cell membrane", "Single-pass type I membrane protein"), ("Cytoplasm", "")]),
                _location_comment([("Cell membrane", "Single-pass type II membrane protein")], molecule="Isoform 2"),
            ]
        }
        self.assertEqual(membrane_hints_from_entry(entry), (("Single-pass type I membrane protein",), True, True))
        self.assertEqual(membrane_hints_from_entry({"comments": [_location_comment([("Cytoplasm", "")])]}), ((), False, False))
        # Also in an organelle membrane (COQ8B-like): types kept, sequence prediction not allowed.
        dual = {"comments": [_location_comment([("Cell membrane", ""), ("Mitochondrion membrane", "")])]}
        self.assertEqual(membrane_hints_from_entry(dual), ((), True, False))

    def test_secreted_chain_names_skip_isoforms_and_exosomes(self) -> None:
        entry = {
            "comments": [
                _location_comment([("Cytoplasm", ""), ("Nucleus", "")], molecule="N-VEGF"),
                _location_comment([("Secreted", "")], molecule="VEGFA"),
                _location_comment([("Secreted", "")], molecule="Isoform VEGF165"),
                _location_comment([("Secreted, extracellular exosome", "")], molecule="Fragment X"),
            ]
        }
        self.assertEqual(secreted_chain_names(entry), ["VEGFA"])


class PipelineDesignRegionTests(unittest.TestCase):
    def _analyzer(self) -> AntigenAnalyzer:
        return AntigenAnalyzer(config=AnalysisConfig(min_construct_length=20, generate_assets=False))

    def _target(self, accession: str, entry_name: str, gene: str, length: int) -> TargetRecord:
        return TargetRecord(
            accession=accession,
            entry_name=entry_name,
            gene_symbol=gene,
            protein_name=gene,
            organism="Homo sapiens",
            taxon_id=9606,
            sequence="A" * length,
        )

    def test_secreted_processed_chain_becomes_the_design_region(self) -> None:
        # IL16: pro-IL-16 1-1332, no signal peptide; UniProt places chain "Interleukin-16" (1212-1332) as secreted.
        analyzer = self._analyzer()
        target = self._target("Q14005", "IL16_HUMAN", "IL16", 1332)
        features = [
            Feature(type="CHAIN", start=1, end=1332, description="Pro-interleukin-16"),
            Feature(type="CHAIN", start=1212, end=1332, description="Interleukin-16"),
        ]
        entry = {"comments": [_location_comment([("Secreted", "")], molecule="Interleukin-16"),
                              _location_comment([("Cytoplasm", "")], molecule="Isoform 1")]}
        topology = derive_ectodomain(features, 1332, analyzer.config)
        self.assertIsNone(topology.ectodomain)
        note = analyzer._apply_secreted_chain_fallback(target=target, topology=topology, features=features, entry=entry)
        self.assertIsNotNone(note)
        self.assertEqual((topology.ectodomain.start, topology.ectodomain.end), (1212, 1332))
        self.assertEqual(topology.topology.topology_class, "secreted")

    def test_curated_override_sets_the_vegfa_design_region(self) -> None:
        with TemporaryDirectory() as root:
            inputs = Path(root) / "data" / "inputs"
            inputs.mkdir(parents=True)
            (inputs / "design_region_overrides.tsv").write_text(
                "accession\tentry_name\tstart\tend\tlabel\tevidence\n"
                "P15692\tVEGFA_HUMAN\t207\t395\tSecreted VEGF-A\talternative initiation + 26-residue signal peptide\n",
                encoding="utf-8",
            )
            cwd = Path.cwd()
            try:
                os.chdir(root)
                analyzer = self._analyzer()
                target = self._target("P15692", "VEGFA_HUMAN", "VEGFA", 395)
                topology = derive_ectodomain([], 395, analyzer.config)
                note = analyzer._apply_design_region_override(target=target, topology=topology)
            finally:
                os.chdir(cwd)
        self.assertIsNotNone(note)
        self.assertEqual((topology.ectodomain.start, topology.ectodomain.end), (207, 395))
        self.assertEqual(topology.topology.topology_class, "secreted")


if __name__ == "__main__":
    unittest.main()
