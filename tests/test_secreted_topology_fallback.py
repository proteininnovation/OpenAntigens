import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from agdesign2.config import AnalysisConfig
from agdesign2.models import Feature, TargetRecord
from agdesign2.pipeline import AntigenAnalyzer
from agdesign2.topology import derive_ectodomain


class SecretedTopologyFallbackTests(unittest.TestCase):
    def test_curated_secreted_target_uses_chain_boundaries_when_topology_is_missing(self) -> None:
        with TemporaryDirectory() as root:
            input_dir = Path(root) / "data" / "inputs"
            input_dir.mkdir(parents=True)
            (input_dir / "accessible_human_topology_refined_accessible_only.csv").write_text(
                "accession,uniprot_entry_name,gene_symbol,primary_topology_class,simplified_topology_bucket,final_accessibility_bucket\n"
                "P01584,IL1B_HUMAN,IL1B,secreted_soluble,Secreted,high_confidence\n",
                encoding="utf-8",
            )
            analyzer = AntigenAnalyzer(config=AnalysisConfig(min_construct_length=20, generate_assets=False))
            target = TargetRecord(
                accession="P01584",
                entry_name="IL1B_HUMAN",
                gene_symbol="IL1B",
                protein_name="Interleukin-1 beta",
                organism="Homo sapiens",
                taxon_id=9606,
                sequence="M" * 269,
            )
            topology = derive_ectodomain([], len(target.sequence), analyzer.config)
            features = [Feature(type="CHAIN", start=117, end=269, description="Interleukin-1 beta")]

            cwd = Path.cwd()
            try:
                import os

                os.chdir(root)
                note = analyzer._apply_secreted_universe_fallback(target=target, topology=topology, features=features)
            finally:
                os.chdir(cwd)

        self.assertIsNotNone(note)
        self.assertIsNotNone(topology.ectodomain)
        self.assertEqual((topology.ectodomain.start, topology.ectodomain.end), (117, 269))
        self.assertEqual(topology.topology.topology_class, "secreted_external_evidence")

    def test_curated_secreted_target_uses_full_length_when_no_boundaries_exist(self) -> None:
        with TemporaryDirectory() as root:
            input_dir = Path(root) / "data" / "inputs"
            input_dir.mkdir(parents=True)
            (input_dir / "accessible_human_topology_refined_accessible_only.csv").write_text(
                "accession,uniprot_entry_name,gene_symbol,primary_topology_class,simplified_topology_bucket,final_accessibility_bucket\n"
                "P01583,IL1A_HUMAN,IL1A,secreted_soluble,Secreted,high_confidence\n",
                encoding="utf-8",
            )
            analyzer = AntigenAnalyzer(config=AnalysisConfig(min_construct_length=20, generate_assets=False))
            target = TargetRecord(
                accession="P01583",
                entry_name="IL1A_HUMAN",
                gene_symbol="IL1A",
                protein_name="Interleukin-1 alpha",
                organism="Homo sapiens",
                taxon_id=9606,
                sequence="M" * 271,
            )
            topology = derive_ectodomain([], len(target.sequence), analyzer.config)

            cwd = Path.cwd()
            try:
                import os

                os.chdir(root)
                note = analyzer._apply_secreted_universe_fallback(target=target, topology=topology, features=[])
            finally:
                os.chdir(cwd)

        self.assertIsNotNone(note)
        self.assertIsNotNone(topology.ectodomain)
        self.assertEqual((topology.ectodomain.start, topology.ectodomain.end), (1, 271))
        self.assertEqual(topology.ectodomain.metadata["fallback_reason"], "full_length")

    def test_isoform_only_secreted_target_is_not_treated_as_curated_secreted(self) -> None:
        # secreted_isoform_only means secreted evidence exists only for a
        # non-canonical isoform; the canonical sequence (what we design) is not
        # the secreted species, so the secreted-universe fallback must NOT fire
        # and must not fabricate a full-canonical secreted design region.
        with TemporaryDirectory() as root:
            input_dir = Path(root) / "data" / "inputs"
            input_dir.mkdir(parents=True)
            (input_dir / "accessible_human_topology_refined_accessible_only.csv").write_text(
                "accession,uniprot_entry_name,gene_symbol,primary_topology_class,simplified_topology_bucket,final_accessibility_bucket\n"
                "Q9ISO1,ISO1_HUMAN,ISO1,secreted_isoform_only,Secreted,secreted_isoform_only_accessible\n",
                encoding="utf-8",
            )
            analyzer = AntigenAnalyzer(config=AnalysisConfig(min_construct_length=20, generate_assets=False))
            target = TargetRecord(
                accession="Q9ISO1",
                entry_name="ISO1_HUMAN",
                gene_symbol="ISO1",
                protein_name="Isoform-only secreted test target",
                organism="Homo sapiens",
                taxon_id=9606,
                sequence="M" * 250,
            )
            topology = derive_ectodomain([], len(target.sequence), analyzer.config)

            cwd = Path.cwd()
            try:
                import os

                os.chdir(root)
                self.assertFalse(analyzer._target_is_curated_secreted(target))
                note = analyzer._apply_secreted_universe_fallback(target=target, topology=topology, features=[])
            finally:
                os.chdir(cwd)

        self.assertIsNone(note)
        self.assertIsNone(topology.ectodomain)


if __name__ == "__main__":
    unittest.main()
