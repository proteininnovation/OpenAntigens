from __future__ import annotations

import copy
import gzip
import hashlib
import json
import os
import socket
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from agdesign2.public_api_data import (
    PublicDataError, build_public_api_bundle, canonical_bytes,
    project_public_target, sequence_fasta, validate_public_api_bundle,
    semantic_hash,
)
from agdesign2.public_api_schema import PUBLIC_SCHEMA, validate_schema


def source_fixture(root: Path) -> tuple[Path, Path, Path]:
    """Small frozen reports and finalized public artifacts, with distinct species sequences."""
    site = root / "site"
    summaries = []
    for species, taxon, accession in (("human", 9606, "P00001"), ("mouse", 10090, "P00002")):
        batch = root / species
        batch.mkdir(parents=True)
        sequence = ("ACDEFGHIKLMNPQRSTVWY" if species == "human" else "MCDEFGHIKLMNPQRSTVWY") * 3
        entry_name = f"SYNTHETIC_{species.upper()}"
        report = {
            "target": {"accession": accession, "entry_name": entry_name, "gene_symbol": "SYNTHETIC", "protein_name": "Synthetic software fixture", "taxon_id": taxon, "sequence": sequence, "alternative_names": ["Fixture alias"], "gene_synonyms": ["SYNT"]},
            "topology": {"topology_class": "single_pass_type_i", "source": "synthetic", "extracellular_regions": [{"start": 5, "end": 30, "label": "Extracellular", "source": "synthetic"}]},
            "ectodomain": {"start": 5, "end": 30, "source": "synthetic", "label": "Design region"},
            "construct_details": [{"name": "strict_fixture", "start": 5, "end": 25, "length": 21, "sequence": sequence[4:25], "classification": "soluble_ecd", "score": 0.75, "warnings": ["Synthetic warning"], "evidence": ["synthetic"], "rationale": "Software fixture"}],
            "interpro_annotations": [{"accession": "IPR000001", "name": "Synthetic domain", "type": "domain", "source_database": "InterPro", "start": 5, "end": 20, "metadata": {"secret": "PRIVATE_CANARY"}}],
            "family_context": {"gene_symbol": "SYNTHETIC", "source": "synthetic", "family_names": ["Synthetic family"], "members": [{"gene_symbol": "SYNT2", "gene_name": "Synthetic family member", "accession": "P00003", "entry_name": "SYNT2_HUMAN", "notes": []}], "metadata": {"secret": "PRIVATE_CANARY"}, "identity_matrix": [[100]]},
            "ectodomain_homology": [{"species": "mouse" if species == "human" else "human", "accession": "P00002" if species == "human" else "P00001", "available": False, "selection_status": "withheld: mapping uncertain", "notes": ["Mapping withheld in synthetic fixture"]}],
            "assembly_requirements": [{"summary": "Synthetic optional assembly", "obligatory": False, "confidence": "low", "source": "synthetic", "partners": ["SYNT2"]}],
            "experimental_constructs": [{"pdb_id": "1ABC", "method": "Synthetic", "resolution": "2.0", "chains": [{"start": 5, "end": 25, "source": "synthetic", "label": "A"}]}],
            "cross_reactivity_hits": [{"subject_id": "P00003", "identity": 50, "coverage": 80, "alignment_length": 20, "query_alignment": "PRIVATE_CANARY", "subject_alignment": "PRIVATE_CANARY"}],
            "ptms": [], "literature_context": {"count": 0, "query_url": "https://www.ncbi.nlm.nih.gov/research/pubtator3/docsum?text=synthetic"},
            "metadata": {"private": "PRIVATE_CANARY"},
        }
        if species == "mouse":
            report["metadata"]["mouse_sequence_provenance"] = {"identifier_source": "human_report_mouse_homolog", "sequence_source": "ortholog_reference_table.mouse_refseq_sequence", "refseq_accession": "NP_000001.1"}
        report_path = batch / "report.json"
        report_path.write_text(json.dumps(report))
        summary = batch / "batch_summary.json"
        summary.write_text(json.dumps([{"query": entry_name, "resolved_entry_name": entry_name, "status": "ok", "json_report": "report.json"}]))
        summaries.append(summary)
        portal = site / ("mouse" if species == "mouse" else "")
        (portal / "downloads").mkdir(parents=True)
        (portal / "reports").mkdir()
        (portal / "reports" / "synthetic.html").write_text("Synthetic public report")
        (portal / "portal_metadata.json").write_text(json.dumps({"target_count": 1}))
        final_entry = {"entry_name": entry_name, "gene_symbol": "SYNTHETIC", "protein_name": "Synthetic software fixture", "status": "ok", "topology_bucket": "GPI", "construct_count": 1, "has_any_structure": 1, "has_pdb": 1, "has_alphafold_structure": 0, "detail_page": "reports/synthetic.html", "pubtator_count": 0}
        (portal / "downloads/agdesign2_portal_index.json").write_text(json.dumps([final_entry]))
    return summaries[0], summaries[1], site


class PublicApiDataTests(unittest.TestCase):
    def build(self, root: Path, *, previous: Path | None = None, snapshot: str = "synthetic") -> Path:
        human, mouse, site = source_fixture(root)
        build_public_api_bundle(human, mouse, root / "bundle", public_site_dir=site, software_commit="0" * 40, snapshot_id=snapshot, previous_bundle_dir=previous)
        return root / "bundle"

    def release(self, bundle: Path) -> dict:
        root = validate_public_api_bundle(bundle)
        return json.loads((bundle / root["releases"][0]["manifest_path"]).read_bytes())

    def test_species_provenance_frozen_fields_exports_and_no_private_leak(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with mock.patch.object(socket, "create_connection", side_effect=AssertionError("network access")):
                bundle = self.build(root)
            release = self.release(bundle)
            for species, identifier in (("human", "9606:uniprot:P00001"), ("mouse", "10090:uniprot:P00002")):
                info = release["species"][species]
                index = json.loads((bundle / info["index_path"]).read_bytes())["data"]
                self.assertEqual(index[0]["target_id"], identifier)
                self.assertEqual(index[0]["track"], "gpi")
                target = json.loads((bundle / info["target_paths"][identifier]).read_bytes())["data"]
                summary = target["summary"]
                self.assertEqual(index[0]["summary_path"], info["target_paths"][identifier])
                self.assertNotIn("aliases", index[0])
                self.assertEqual(summary["analyzed_topology"], "single_pass_type_i")
                self.assertEqual(summary["section_availability"]["ptms"], "unknown")
                self.assertEqual(summary["sequence_accession"], "NP_000001.1" if species == "mouse" else "P00001")
                sequence = json.loads((bundle / target["reference_sequence_path"]).read_bytes())["data"]
                self.assertEqual(sequence["sequence"][0], "M" if species == "mouse" else "A")
                self.assertIn(f"sequence_accession={sequence['sequence_accession']}", sequence_fasta(sequence, release_id=release["release_id"], provenance_url=release["manifest_url"], records_url=info["records_download"]["url"]))
                constructs = json.loads((bundle / target["constructs"]["pages"][0]).read_bytes())["data"]
                self.assertEqual(constructs[0]["category"], "strict")
                self.assertEqual(constructs[0]["classification"], "soluble_ecd")
                for path in bundle.rglob("*.json"):
                    content = path.read_text()
                    self.assertNotIn("PRIVATE_CANARY", content)
                    self.assertNotIn(str(root), content)
                metadata = json.loads((root / "site" / ("mouse" if species == "mouse" else "") / "portal_metadata.json").read_bytes())
                self.assertEqual(metadata["public_api"]["state"], "packaged")
                self.assertNotIn("service_url", metadata["public_api"])
                self.assertNotIn("mcp_url", metadata["public_api"])

    def test_frozen_bundle_is_deterministic_with_ambient_cache_poison(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            human, mouse, site = source_fixture(root)
            first = root / "first"
            build_public_api_bundle(human, mouse, first, public_site_dir=site, software_commit="0" * 40, snapshot_id="synthetic")
            poison = root / "ambient/.agdesign2/cache"
            poison.mkdir(parents=True)
            (poison / "synthetic.json").write_text('{"sequence":"PRIVATE_CANARY"}')
            second = root / "second"
            previous_cwd = Path.cwd()
            try:
                os.chdir(root / "ambient")
                build_public_api_bundle(human, mouse, second, public_site_dir=site, software_commit="0" * 40, snapshot_id="synthetic")
            finally:
                os.chdir(previous_cwd)
            for relative in (path.relative_to(first) for path in first.rglob("*") if path.is_file()):
                self.assertEqual((first / relative).read_bytes(), (second / relative).read_bytes(), str(relative))

    def test_retention_copies_previous_bytes_and_static_downloads(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            previous = self.build(root / "previous")
            current = self.build(root / "current", previous=previous, snapshot="next")
            old = self.release(previous)
            new_root = validate_public_api_bundle(current)
            self.assertEqual(len(new_root["releases"]), 2)
            for relative in old["files"]:
                self.assertEqual((previous / relative).read_bytes(), (current / relative).read_bytes())
            for info in old["species"].values():
                for download in info["downloads"]:
                    from urllib.parse import urlsplit
                    path = root / "current/site" / urlsplit(download["url"]).path.lstrip("/")
                    self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), download["sha256"])

    def test_coordinate_mismatch_and_species_mismatch_fail_before_publication(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            human, mouse, site = source_fixture(root)
            report_path = human.parent / "report.json"
            report = json.loads(report_path.read_bytes())
            report["construct_details"][0]["end"] = 999
            report_path.write_text(json.dumps(report))
            with self.assertRaisesRegex(PublicDataError, "coordinates"):
                build_public_api_bundle(human, mouse, root / "bad", public_site_dir=site, software_commit="0" * 40, snapshot_id="bad")
            report["construct_details"][0]["end"] = 25
            report["target"]["taxon_id"] = 10090
            report_path.write_text(json.dumps(report))
            with self.assertRaisesRegex(PublicDataError, "species"):
                build_public_api_bundle(human, mouse, root / "bad2", public_site_dir=site, software_commit="0" * 40, snapshot_id="bad")

    def test_recursive_schema_rejects_canary_and_unknown_metric(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bundle = self.build(Path(tmp))
            release = self.release(bundle)
            info = release["species"]["human"]
            target = json.loads((bundle / next(iter(info["target_paths"].values()))).read_bytes())["data"]
            constructs = json.loads((bundle / target["constructs"]["pages"][0]).read_bytes())
            constructs["data"][0]["structural_metrics"] = {"private": "PRIVATE_CANARY"}
            with self.assertRaises(PublicDataError):
                validate_schema("construct_page", constructs)

    def test_bundle_rejects_changed_bytes_unlisted_files_and_symlinks(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bundle = self.build(Path(tmp))
            unexpected = bundle / "unexpected.json"
            unexpected.write_text("{}")
            with self.assertRaisesRegex(PublicDataError, "Unexpected"):
                validate_public_api_bundle(bundle)
            unexpected.unlink()
            unexpected.symlink_to(bundle / "manifest.json")
            with self.assertRaises(PublicDataError):
                validate_public_api_bundle(bundle)
            unexpected.unlink()
            release = self.release(bundle)
            target_file = bundle / next(iter(release["species"]["human"]["target_paths"].values()))
            target_file.write_text("{}")
            with self.assertRaisesRegex(PublicDataError, "checksum"):
                validate_public_api_bundle(bundle)

    def test_full_records_fallback_matches_all_public_shards(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bundle = self.build(Path(tmp))
            release = self.release(bundle)
            for species, info in release["species"].items():
                rows = [json.loads(line) for line in gzip.decompress((bundle / info["records_download"]["path"]).read_bytes()).splitlines()]
                expected = {path for path, metadata in release["files"].items() if metadata["schema"] and (metadata["schema"] in {"index", "target_shard", "construct_page", "sequence_shard"} or metadata["schema"].startswith("section_")) and f"/{species}/" in path}
                self.assertEqual({row["path"] for row in rows}, expected)
                for row in rows:
                    self.assertEqual(row["payload"], json.loads((bundle / row["path"]).read_bytes()))

    def test_native_pdb_window_preserved_variant_export_withheld(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            human, _, site = source_fixture(Path(tmp))
            report = json.loads((human.parent / "report.json").read_bytes())
            entry = json.loads((site / "downloads/agdesign2_portal_index.json").read_bytes())[0]
            report["construct_details"][0].update({"name": "pdb_1abc_A", "pdb_id": "1ABC"})
            projection = project_public_target(report, entry, "human", public_site_dir=site)
            self.assertEqual(projection["constructs"][0]["category"], "pdb")
            self.assertEqual(projection["constructs"][0]["sequence_basis"], "reference_slice")
            self.assertEqual(projection["constructs"][0]["export_state"], "available")
            report["construct_details"][0]["sequence"] = "M" + report["construct_details"][0]["sequence"][1:]
            projection = project_public_target(report, entry, "human", public_site_dir=site)
            self.assertEqual(projection["constructs"][0]["export_state"], "unavailable")
            self.assertIsNone(next(iter(projection["construct_sequences"].values())))

    def test_schema_resource_matches_validator(self) -> None:
        path = Path(__file__).resolve().parents[1] / "src/agdesign2/data/public_api/schema.json"
        self.assertEqual(json.loads(path.read_text()), PUBLIC_SCHEMA)

    def test_processing_hazards_survive_dataclass_projection_and_sequence_export(self) -> None:
        from agdesign2.models import CysteineFinding, DomainAnnotation, Feature, FurinSite
        with tempfile.TemporaryDirectory() as tmp:
            human, _, site = source_fixture(Path(tmp))
            report = json.loads((human.parent / "report.json").read_bytes())
            entry = json.loads((site / "downloads/agdesign2_portal_index.json").read_bytes())[0]
            cysteine = CysteineFinding(position=8, paired_with=None, surface_exposed=True, warning="Recorded disulfide partner lies outside the selected range.", partner_outside_range=40).to_dict()
            furin = FurinSite(start=10, end=13, motif="RKKR", suggestions=["R10A"], curated_annotation="Recorded processing site.").to_dict()
            ptm = Feature(type="CARBOHYD", start=15, end=15, description="Recorded glycosylation site.").to_dict()
            domain = DomainAnnotation(accession="IPR000002", name="Recorded clipped domain", type="domain", source_database="InterPro", start=1, end=20).to_dict()
            report.update({"cysteine_analysis": [cysteine], "furin_sites": [furin], "ptms": [ptm]})
            report["construct_details"][0].update({"cysteine_analysis": [cysteine], "ptms": [ptm], "ligand_interactions": [Feature(type="BINDING", start=20, end=20, description="Recorded ligand contact.").to_dict()], "clipped_domain_annotations": [domain]})
            projection = project_public_target(report, entry, "human", public_site_dir=site)
            construct = projection["constructs"][0]
            self.assertEqual(construct["hazards"]["cysteines"][0]["partner_outside_range"], 40)
            self.assertEqual(construct["hazards"]["furin_sites"][0]["curated_annotation"], "Recorded processing site.")
            self.assertEqual(construct["hazards"]["ptms"][0]["type"], "CARBOHYD")
            self.assertEqual(construct["hazards"]["ligand_interactions"][0]["type"], "BINDING")
            self.assertEqual(construct["hazards"]["clipped_domains"][0]["accession"], "IPR000002")
            warnings = projection["construct_sequences"][construct["construct_id"]]["warnings"]
            self.assertTrue(any(text.endswith(cysteine["warning"]) and "reference interval 5-25" in text for text in warnings))
            self.assertTrue(any("Recorded processing site." in text for text in warnings))
            self.assertTrue(any("glycosylation" in text for text in warnings))
            self.assertEqual({row["kind"] for row in projection["sections"]["ptms"]["items"]}, {"feature", "furin", "cysteine"})

    def test_complex_portal_records_and_lookup_states_are_typed(self) -> None:
        from agdesign2.models import ComplexPortalComplex, ComplexPortalParticipant
        with tempfile.TemporaryDirectory() as tmp:
            human, _, site = source_fixture(Path(tmp))
            report = json.loads((human.parent / "report.json").read_bytes())
            entry = json.loads((site / "downloads/agdesign2_portal_index.json").read_bytes())[0]
            report["complex_portal_complexes"] = [ComplexPortalComplex(complex_ac="CPX-0001", name="Synthetic complex", species="human", predicted_complex=True, participants=[ComplexPortalParticipant(identifier="P00001", name="Synthetic partner", stoichiometry="1")]).to_dict()]
            report["complex_portal_lookup"] = {"human": {"status": "ok"}, "mouse": {"status": "error", "error": "/private/PRIVATE_CANARY"}}
            section = project_public_target(report, entry, "human", public_site_dir=site)["sections"]["assembly"]
            complex_row = next(row for row in section["items"] if row["kind"] == "complex")
            self.assertTrue(complex_row["complex"]["predicted_complex"])
            self.assertEqual(complex_row["complex"]["participants"][0]["identifier"], "P00001")
            self.assertNotIn("PRIVATE_CANARY", json.dumps(section))
            report["complex_portal_complexes"] = []
            report["assembly_requirements"] = []
            for status, expected in (("no_hits", "no_results"), ("error", "error"), ("disabled", "unavailable")):
                report["complex_portal_lookup"] = {"human": {"status": status}}
                section = project_public_target(report, entry, "human", public_site_dir=site)["sections"]["assembly"]
                self.assertEqual(section["state"], expected)

    def test_human_primary_sequence_never_inherits_mouse_ortholog_source(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            human, _, site = source_fixture(Path(tmp))
            report = json.loads((human.parent / "report.json").read_bytes())
            entry = json.loads((site / "downloads/agdesign2_portal_index.json").read_bytes())[0]
            entry.update({"mouse_sequence_source": "ortholog_reference_table.mouse_refseq_sequence", "mouse_refseq_accession": "NP_000099.2"})
            summary = project_public_target(report, entry, "human", public_site_dir=site)["summary"]
            self.assertEqual(summary["sequence_accession"], "P00001")
            self.assertEqual(summary["sequence_source"], "uniprot")

    def test_exact_duplicate_constructs_merge_evidence_without_losing_source_count(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            human, _, site = source_fixture(Path(tmp))
            report = json.loads((human.parent / "report.json").read_bytes())
            entry = json.loads((site / "downloads/agdesign2_portal_index.json").read_bytes())[0]
            duplicate = copy.deepcopy(report["construct_details"][0])
            duplicate["evidence"] = ["second recorded evidence"]
            report["construct_details"].append(duplicate)
            entry["construct_count"] = 2
            projection = project_public_target(report, entry, "human", public_site_dir=site)
            self.assertEqual(projection["summary"]["portal_construct_count"], 2)
            self.assertEqual(projection["summary"]["construct_count"], 1)
            self.assertIn("second recorded evidence", projection["constructs"][0]["evidence"])
            self.assertEqual(next(iter(projection["construct_sequences"].values()))["evidence"], projection["constructs"][0]["evidence"])

    def test_release_hash_binds_shards_even_when_transport_checksums_are_rewritten(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bundle = self.build(Path(tmp))
            root = json.loads((bundle / "manifest.json").read_bytes())
            release_pointer = root["releases"][0]
            release = json.loads((bundle / release_pointer["manifest_path"]).read_bytes())
            # An attacker can rewrite a public count and transport hashes; the immutable core must still reject it.
            release["species"]["human"]["counts"]["constructs"] += 1
            content = canonical_bytes(release) + b"\n"
            (bundle / release_pointer["manifest_path"]).write_bytes(content)
            release_pointer["sha256"] = hashlib.sha256(content).hexdigest()
            root_bytes = canonical_bytes(root) + b"\n"
            (bundle / "manifest.json").write_bytes(root_bytes)
            (bundle / "manifest.sha256").write_text(hashlib.sha256(root_bytes).hexdigest() + "  manifest.json\n")
            with self.assertRaisesRegex(PublicDataError, "aggregate construct"):
                validate_public_api_bundle(bundle)

    def test_unresolved_row_is_retained_without_an_invented_target_or_sequence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            human, mouse, site = source_fixture(root)
            summary_rows = json.loads(human.read_bytes())
            summary_rows.append({"query": "UNRESOLVED_HUMAN", "status": "error", "error": "/private/PRIVATE_CANARY"})
            human.write_text(json.dumps(summary_rows))
            index_path = site / "downloads/agdesign2_portal_index.json"
            index_rows = json.loads(index_path.read_bytes())
            index_rows.append({"entry_name": "UNRESOLVED_HUMAN", "status": "error", "error": "analysis failed"})
            index_path.write_text(json.dumps(index_rows))
            build_public_api_bundle(human, mouse, root / "bundle", public_site_dir=site, software_commit="0" * 40, snapshot_id="unresolved")
            release = self.release(root / "bundle")
            info = release["species"]["human"]
            self.assertEqual(info["counts"]["targets"], 2)
            self.assertEqual(info["counts"]["resolved_targets"], 1)
            rows = json.loads((root / "bundle" / info["index_path"]).read_bytes())["data"]
            unresolved = next(row for row in rows if row["target_id"] is None)
            unresolved_summary = json.loads((root / "bundle" / unresolved["summary_path"]).read_bytes())["data"]["summary"]
            self.assertIsNone(unresolved_summary["sequence_sha256"])
            self.assertEqual(unresolved["status"], "error")

    def test_identity_tuple_digest_and_portable_numeric_digest_are_frozen(self) -> None:
        self.assertEqual(semantic_hash({"integer": 50, "decimal": 50.0, "tiny": 1e-7, "negative_zero": -0.0, "boolean": False}), "6541843bb6ae61c70aa65dfcb7aeb531d4a8736b3918b96aa24c930a4b1f9ede")
        with tempfile.TemporaryDirectory() as tmp:
            human, _, site = source_fixture(Path(tmp))
            report = json.loads((human.parent / "report.json").read_bytes())
            entry = json.loads((site / "downloads/agdesign2_portal_index.json").read_bytes())[0]
            projection = project_public_target(report, entry, "human", public_site_dir=site)
            self.assertEqual(projection["constructs"][0]["construct_id"], "c40b78c5de4eaf8ddf49072fcd5c0d282aefd5f7a1ab6c23c4e05598be587e9c")

    def test_rewritten_shards_and_download_hashes_cannot_reuse_release_id(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bundle = self.build(Path(tmp))
            root = json.loads((bundle / "manifest.json").read_bytes())
            pointer = root["releases"][0]
            release = json.loads((bundle / pointer["manifest_path"]).read_bytes())
            info = release["species"]["human"]

            def replace_file(relative: str, content: bytes) -> None:
                (bundle / relative).write_bytes(content)
                metadata = release["files"][relative]
                metadata["sha256"], metadata["bytes"] = hashlib.sha256(content).hexdigest(), len(content)
                for species_info in release["species"].values():
                    for download in [*species_info["downloads"], species_info["records_download"]]:
                        if download["path"] == relative:
                            download.update({key: metadata[key] for key in ("sha256", "bytes")})

            index = json.loads((bundle / info["index_path"]).read_bytes())
            index["data"][0]["search_text"] += " altered"
            replace_file(info["index_path"], canonical_bytes(index) + b"\n")
            target_path = next(iter(info["target_paths"].values()))
            target = json.loads((bundle / target_path).read_bytes())
            target["data"]["summary"]["search_text"] = index["data"][0]["search_text"]
            replace_file(target_path, canonical_bytes(target) + b"\n")
            fallback_path = info["records_download"]["path"]
            fallback = [json.loads(line) for line in gzip.decompress((bundle / fallback_path).read_bytes()).splitlines()]
            for record in fallback:
                record["payload"] = json.loads((bundle / record["path"]).read_bytes())
            replace_file(fallback_path, gzip.compress(b"".join(canonical_bytes(record) + b"\n" for record in fallback), mtime=0))
            release_bytes = canonical_bytes(release) + b"\n"
            (bundle / pointer["manifest_path"]).write_bytes(release_bytes)
            pointer["sha256"] = hashlib.sha256(release_bytes).hexdigest()
            root_bytes = canonical_bytes(root) + b"\n"
            (bundle / "manifest.json").write_bytes(root_bytes)
            (bundle / "manifest.sha256").write_text(hashlib.sha256(root_bytes).hexdigest() + "  manifest.json\n")
            with self.assertRaisesRegex(PublicDataError, "immutable canonical payload"):
                validate_public_api_bundle(bundle)

    def test_metadata_symlinks_fail_before_any_bundle_or_static_mutation(self) -> None:
        for species in ("human", "mouse"):
            with self.subTest(species=species), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                human, mouse, site = source_fixture(root)
                victim = root / "outside.json"
                victim.write_bytes(b'{"sentinel":"preserve bytes"}\n')
                metadata = site / ("mouse" if species == "mouse" else "") / "portal_metadata.json"
                metadata.unlink()
                metadata.symlink_to(victim)
                before = {str(path.relative_to(site)): path.read_bytes() for path in site.rglob("*") if path.is_file()}
                with self.assertRaisesRegex(PublicDataError, "symbolic"):
                    build_public_api_bundle(human, mouse, root / "bundle", public_site_dir=site, software_commit="0" * 40, snapshot_id="symlink")
                self.assertEqual(victim.read_bytes(), b'{"sentinel":"preserve bytes"}\n')
                self.assertFalse((root / "bundle").exists())
                self.assertEqual(before, {str(path.relative_to(site)): path.read_bytes() for path in site.rglob("*") if path.is_file()})

    def test_second_metadata_failure_rolls_back_first_exact_bytes(self) -> None:
        import stat
        from agdesign2.public_api_data import _atomic_public_replace
        for modes in ((0o644, 0o640), (0o640, 0o644)):
            with self.subTest(modes=modes), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                human, mouse, site = source_fixture(root)
                paths = [site / "portal_metadata.json", site / "mouse/portal_metadata.json"]
                for path, mode in zip(paths, modes):
                    path.chmod(mode)
                before = [(path.read_bytes(), stat.S_IMODE(path.stat().st_mode)) for path in paths]
                victim = root / "outside.txt"
                victim.write_bytes(b"Unrelated outside bytes remain unchanged.")
                victim.chmod(0o640)
                outside = (victim.read_bytes(), stat.S_IMODE(victim.stat().st_mode))
                calls = []

                def fail_second(directory, relative, content):
                    calls.append(relative)
                    if len(calls) == 2:
                        raise OSError("synthetic atomic replace failure")
                    return _atomic_public_replace(directory, relative, content)

                with mock.patch("agdesign2.public_api_data._atomic_public_replace", side_effect=fail_second):
                    with self.assertRaisesRegex(OSError, "synthetic atomic"):
                        build_public_api_bundle(human, mouse, root / "bundle", public_site_dir=site, software_commit="0" * 40, snapshot_id="rollback")
                self.assertEqual(before, [(path.read_bytes(), stat.S_IMODE(path.stat().st_mode)) for path in paths])
                self.assertEqual(outside, (victim.read_bytes(), stat.S_IMODE(victim.stat().st_mode)))
                self.assertEqual(calls, ["portal_metadata.json", "mouse/portal_metadata.json", "portal_metadata.json"])

    def test_metadata_success_preserves_human_mouse_modes_and_external_hardlink(self) -> None:
        import stat
        for modes in ((0o644, 0o640), (0o640, 0o644)):
            with self.subTest(modes=modes), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                human, mouse, site = source_fixture(root)
                paths = [site / "portal_metadata.json", site / "mouse/portal_metadata.json"]
                for path, mode in zip(paths, modes):
                    path.chmod(mode)
                victim = root / "outside.json"
                os.link(paths[0], victim)
                outside = (victim.read_bytes(), stat.S_IMODE(victim.stat().st_mode))
                build_public_api_bundle(human, mouse, root / "bundle", public_site_dir=site, software_commit="0" * 40, snapshot_id="mode-preservation")
                self.assertEqual(tuple(stat.S_IMODE(path.stat().st_mode) for path in paths), modes)
                self.assertTrue(all(json.loads(path.read_bytes())["public_api"]["state"] == "packaged" for path in paths))
                self.assertEqual(outside, (victim.read_bytes(), stat.S_IMODE(victim.stat().st_mode)))
                self.assertNotEqual(paths[0].stat().st_ino, victim.stat().st_ino)

    def test_cysteine_source_interval_is_distinct_from_full_reference_export(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            human, _, site = source_fixture(Path(tmp))
            report = json.loads((human.parent / "report.json").read_bytes())
            entry = json.loads((site / "downloads/agdesign2_portal_index.json").read_bytes())[0]
            report["cysteine_analysis"] = [{"position": 8, "paired_with": None, "surface_exposed": True, "partner_outside_range": 40, "warning": "Partner is outside the analyzed range."}]
            projection = project_public_target(report, entry, "human", public_site_dir=site)
            reference = projection["reference_sequence"]
            self.assertFalse(any("outside exported interval 1-60" in text for text in reference["warnings"]))
            self.assertTrue(any("source analysis for reference interval 5-30" in text for text in reference["warnings"]))
            construct = projection["constructs"][0]
            self.assertEqual(construct["hazards"]["cysteines"][0]["analysis_start"], 5)
            self.assertEqual(construct["hazards"]["cysteines"][0]["analysis_end"], 30)
            sequence = projection["construct_sequences"][construct["construct_id"]]
            self.assertTrue(any("outside exported interval 5-25" in text for text in sequence["warnings"]))
            self.assertTrue(any("Recorded cysteine analysis for reference interval 5-30" in text for text in sequence["warnings"]))

    def test_final_model_sequence_and_provenance_checked_in_plain_and_gzip_assets(self) -> None:
        from Bio.SeqUtils import seq3
        with tempfile.TemporaryDirectory() as tmp:
            human, _, site = source_fixture(Path(tmp))
            report = json.loads((human.parent / "report.json").read_bytes())
            entry = json.loads((site / "downloads/agdesign2_portal_index.json").read_bytes())[0]
            entry.update({"has_alphafold_structure": 1, "structure_source": "AlphaFold"})
            structures = site / "structures"
            structures.mkdir()
            model = structures / "synthetic.pdb"

            def pdb(sequence):
                return "".join(f"ATOM  {pos:5d}  CA  {seq3(aa).upper():3s} A{pos:4d}    {float(pos):8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{90.0:6.2f}           C\n" for pos, aa in enumerate(sequence, 1)).encode("ascii")

            model.write_bytes(pdb("A"))
            with self.assertRaisesRegex(PublicDataError, "model sequence"):
                project_public_target(report, entry, "human", public_site_dir=site)
            matching = pdb(report["target"]["sequence"])
            for content in (matching, gzip.compress(matching, mtime=0)):
                model.write_bytes(content)
                projected = project_public_target(report, entry, "human", public_site_dir=site)
                model_row = next(row for row in projected["sections"]["structures"]["items"] if row["kind"] == "predicted_model")
                self.assertEqual(model_row["sha256"], hashlib.sha256(content).hexdigest())
            model.with_suffix(".meta.json").write_text('{"source":"alphafold_server","model_type":"AF3"}')
            with self.assertRaisesRegex(PublicDataError, "excluded source provenance"):
                project_public_target(report, entry, "human", public_site_dir=site)

    def test_literature_query_is_recovered_from_the_exact_frozen_public_row(self) -> None:
        from agdesign2.portal import _render_index_data_js
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            human, mouse, site = source_fixture(root)
            recorded = "https://www.ncbi.nlm.nih.gov/research/pubtator3/docsum?text=%40GENE_SYNTHETIC"
            for species, summary in (("human", human), ("mouse", mouse)):
                report_path = summary.parent / "report.json"
                report = json.loads(report_path.read_bytes())
                report.pop("literature_context")
                report_path.write_text(json.dumps(report))
                portal = site / ("mouse" if species == "mouse" else "")
                index_path = portal / "downloads/agdesign2_portal_index.json"
                rows = json.loads(index_path.read_bytes())
                rows[0]["pubtator_count"] = 7
                index_path.write_text(json.dumps(rows))
                js = _render_index_data_js([{**rows[0], "pubtator_query_url": recorded}], include_disease_context=species == "human").encode("utf-8")
                (portal / "portal-index-data.js").write_bytes(gzip.compress(js, mtime=0) if species == "mouse" else js)
            with mock.patch.object(socket, "create_connection", side_effect=AssertionError("network access")):
                build_public_api_bundle(human, mouse, root / "bundle", public_site_dir=site, software_commit="0" * 40, snapshot_id="frozen-literature")
            release = self.release(root / "bundle")
            for info in release["species"].values():
                target = json.loads((root / "bundle" / next(iter(info["target_paths"].values()))).read_bytes())["data"]
                page = json.loads((root / "bundle" / target["sections"]["literature"]["pages"][0]).read_bytes())["data"]
                self.assertEqual(page["items"][0]["count"], 7)
                self.assertEqual(page["items"][0]["query_url"], recorded)

    def test_fasta_has_unique_primary_ids_and_lossless_release_scoped_metadata(self) -> None:
        from urllib.parse import unquote
        with tempfile.TemporaryDirectory() as tmp:
            bundle = self.build(Path(tmp))
            release = self.release(bundle)
            primary_ids = []
            for species, info in release["species"].items():
                downloads = {row["name"]: row for row in info["downloads"]}
                for filename in ("reference.fasta", "constructs.fasta"):
                    for line in (bundle / downloads[filename]["path"]).read_text().splitlines():
                        if not line.startswith(">"):
                            continue
                        pieces = line[1:].split(" ")
                        primary_ids.append(pieces[0])
                        values = {key: unquote(value) for key, value in (piece.split("=", 1) for piece in pieces[1:])}
                        self.assertEqual(values["release_id"], release["release_id"])
                        self.assertEqual(values["provenance_url"], release["manifest_url"])
                        self.assertEqual(values["records_url"], info["records_download"]["url"])
                        self.assertIn("sequence_source", values)
                        self.assertIn("sequence_basis", values)
                        self.assertIsInstance(json.loads(values["warnings"]), list)
                        self.assertIsInstance(json.loads(values["evidence"]), list)
                        self.assertEqual(json.loads(values["record_selector"])["target_id"], values["target_id"])
            self.assertEqual(len(primary_ids), len(set(primary_ids)))

    def test_three_builds_keep_known_retired_ids_without_retired_scientific_payloads(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = self.build(root / "one", snapshot="one")
            first_id = validate_public_api_bundle(first)["current_release_id"]
            second = self.build(root / "two", previous=first, snapshot="two")
            third = self.build(root / "three", previous=second, snapshot="three")
            final = validate_public_api_bundle(third)
            self.assertEqual(final["retired_releases"], [{"release_id": first_id, "retired_by_release_id": final["current_release_id"], "archive_url": None}])
            self.assertFalse((third / "releases" / first_id).exists())

    def test_embedded_schema_is_closed_offline_and_compatible_with_older_fields(self) -> None:
        from agdesign2.public_api_data import _approved_release_schema
        old = copy.deepcopy(PUBLIC_SCHEMA)
        old.pop("x-openantigens-fasta-version")
        old.pop("x-openantigens-index-version")
        old.pop("x-openantigens-disease-selection")
        old["$defs"]["summary"]["properties"].pop("disease_coverage")
        old["$defs"]["summary"]["required"].remove("disease_coverage")
        old["$defs"]["section_disease"]["properties"]["data"]["properties"].pop("coverage")
        old["$defs"]["section_disease"]["properties"]["data"]["required"].remove("coverage")
        old["$defs"].pop("disease_coverage")
        old["$defs"]["provenance"]["properties"].pop("projection_methods")
        old["$defs"]["index"]["properties"]["data"]["items"] = {"$ref": "#/$defs/summary"}
        old["$defs"].pop("search_row")
        for field in ("analysis_start", "analysis_end"):
            old["$defs"]["cysteine"]["properties"].pop(field)
            old["$defs"]["cysteine"]["required"].remove(field)
        old["$defs"]["root_manifest"]["properties"].pop("retired_releases")
        old["$defs"].pop("retired_release")
        self.assertEqual(_approved_release_schema(old), 1)
        for unsafe in (
            {"$ref": "https://invalid.example/remote.json"},
            {"$dynamicRef": "file:///private/never-read.json"},
        ):
            bad = copy.deepcopy(PUBLIC_SCHEMA)
            bad["$defs"]["feature"]["properties"]["description"] = unsafe
            with mock.patch.object(socket, "create_connection", side_effect=AssertionError("schema network access")):
                with self.assertRaises(PublicDataError):
                    _approved_release_schema(bad)
        bad = copy.deepcopy(PUBLIC_SCHEMA)
        bad["$defs"]["feature"]["additionalProperties"] = True
        with self.assertRaises(PublicDataError):
            _approved_release_schema(bad)

    def test_fasta_preserves_unknown_sequence_accession_and_exact_warning_json(self) -> None:
        from urllib.parse import unquote
        with tempfile.TemporaryDirectory() as tmp:
            _, mouse, site = source_fixture(Path(tmp))
            report = json.loads((mouse.parent / "report.json").read_bytes())
            report["metadata"]["mouse_sequence_provenance"]["refseq_accession"] = None
            entry = json.loads((site / "mouse/downloads/agdesign2_portal_index.json").read_bytes())[0]
            projected = project_public_target(report, entry, "mouse", public_site_dir=site)
            sequence = projected["reference_sequence"]
            sequence["warnings"] = ["Recorded scope; é, apostrophe ' and newline\nremain exact."]
            rid = "0" * 64
            text = sequence_fasta(sequence, release_id=rid, provenance_url=f"https://openantigens.org/data/releases/{rid}/manifest.json", records_url=f"https://openantigens.org/data/releases/{rid}/mouse/records.jsonl.gz")
            fields = {key: unquote(value) for key, value in (token.split("=", 1) for token in text.splitlines()[0].split(" ")[1:])}
            self.assertEqual(fields["sequence_accession"], "null")
            self.assertEqual(fields["identifier_accession"], "P00002")
            self.assertEqual(fields["sequence_source"], "ortholog_reference_table.mouse_refseq_sequence")
            self.assertEqual(json.loads(fields["warnings"]), sequence["warnings"])

    def test_tombstones_keep_only_the_most_recent_thousand(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = self.build(root / "one", snapshot="one")
            second = self.build(root / "two", previous=first, snapshot="two")
            old = validate_public_api_bundle(second)
            history = [{"release_id": hashlib.sha256(f"synthetic-retired-{n}".encode()).hexdigest(), "retired_by_release_id": old["current_release_id"], "archive_url": None} for n in range(1000)]
            old["retired_releases"] = history
            content = canonical_bytes(old) + b"\n"
            (second / "manifest.json").write_bytes(content)
            (second / "manifest.sha256").write_text(hashlib.sha256(content).hexdigest() + "  manifest.json\n")
            third = self.build(root / "three", previous=second, snapshot="three")
            result = validate_public_api_bundle(third)
            self.assertEqual(len(result["retired_releases"]), 1000)
            self.assertEqual(result["retired_releases"][1:], history[:999])
            self.assertNotIn(history[-1], result["retired_releases"])

    def test_raw_reports_are_read_once_and_released_before_the_next_target(self) -> None:
        import weakref
        from agdesign2.public_api_data import _load_sources, _read_json_hashed, project_public_target as real_project

        class TrackedReport(dict):
            pass

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            human, _, site = source_fixture(root)
            prototype = json.loads((human.parent / "report.json").read_bytes())
            original_entry = json.loads((site / "downloads/agdesign2_portal_index.json").read_bytes())[0]
            summary, index = [], []
            report_paths = set()
            for number in range(4):
                report = copy.deepcopy(prototype)
                entry_name = f"SYNTHETIC_{number}_HUMAN"
                report["target"].update({"entry_name": entry_name, "accession": f"P{number + 10:05d}"})
                report["private_unused_blob"] = "x" * (1024 * 1024)
                path = human.parent / f"report-{number}.json"
                path.write_text(json.dumps(report))
                report_paths.add(path.resolve())
                summary.append({"query": entry_name, "resolved_entry_name": entry_name, "json_report": path.name})
                page = f"reports/synthetic-{number}.html"
                (site / page).write_text("Synthetic public record")
                index.append({**original_entry, "entry_name": entry_name, "detail_page": page})
            human.write_text(json.dumps(summary))
            (site / "downloads/agdesign2_portal_index.json").write_text(json.dumps(index))
            references, projection_references, reads = [], [], []

            def tracked_read(path, **kwargs):
                if path in report_paths:
                    self.assertTrue(all(reference() is None for reference in references), "A preceding raw report remains resident")
                value, checksum = _read_json_hashed(path, **kwargs)
                if path in report_paths:
                    value = TrackedReport(value)
                    references.append(weakref.ref(value))
                    reads.append(path)
                return value, checksum

            def tracked_project(*args, **kwargs):
                self.assertTrue(all(reference() is None for reference in projection_references), "A preceding public projection remains resident")
                value = TrackedReport(real_project(*args, **kwargs))
                projection_references.append(weakref.ref(value))
                return value

            with mock.patch("agdesign2.public_api_data._read_json_hashed", side_effect=tracked_read), mock.patch("agdesign2.public_api_data.project_public_target", new=tracked_project):
                projected, hashes = _load_sources(human, "human", site, "https://openantigens.org", staging_dir=root / "public-projection-stage")
            self.assertEqual(len(projected), 4)
            self.assertEqual(set(reads), report_paths)
            self.assertEqual(len(reads), len(report_paths))
            self.assertTrue(all(reference() is None for reference in references))
            self.assertTrue(all(reference() is None for reference in projection_references))
            self.assertEqual(sum(key.startswith("human:report:") for key in hashes), 4)
            self.assertNotIn("private_unused_blob", canonical_bytes(list(projected)).decode("utf-8"))

    def test_correct_hashes_cannot_authorize_private_downloads_or_forged_citations(self) -> None:
        for role in ("private_download", "citation"):
            with self.subTest(role=role), tempfile.TemporaryDirectory() as tmp:
                bundle = self.build(Path(tmp))
                root = json.loads((bundle / "manifest.json").read_bytes())
                pointer = root["releases"][0]
                release = json.loads((bundle / pointer["manifest_path"]).read_bytes())
                species = release["species"]["human"]
                content = b"PRIVATE_CANARY synthetic confidential input\n"
                if role == "private_download":
                    relative = f"releases/{release['release_id']}/human/downloads/private.sqlite"
                    metadata = {"sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content), "media_type": "application/octet-stream", "schema": None, "encoding": None}
                    release["files"][relative] = metadata
                    species["downloads"].append({"name": "private.sqlite", "url": f"https://openantigens.org/data/releases/{release['release_id']}/human/private.sqlite", "path": relative, **metadata, "row_count": None, "license": "synthetic"})
                else:
                    download = next(row for row in species["downloads"] if row["name"] == "openantigens.bib")
                    relative = download["path"]
                    release["files"][relative].update({"sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)})
                    download.update({"sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)})
                (bundle / relative).write_bytes(content)
                body = canonical_bytes(release) + b"\n"
                (bundle / pointer["manifest_path"]).write_bytes(body)
                pointer["sha256"] = hashlib.sha256(body).hexdigest()
                body = canonical_bytes(root) + b"\n"
                (bundle / "manifest.json").write_bytes(body)
                (bundle / "manifest.sha256").write_text(hashlib.sha256(body).hexdigest() + "  manifest.json\n")
                with self.assertRaises(PublicDataError):
                    validate_public_api_bundle(bundle)

    def test_schema_valid_unreachable_file_is_rejected_even_with_complete_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bundle = self.build(Path(tmp))
            root = json.loads((bundle / "manifest.json").read_bytes())
            pointer = root["releases"][0]
            release = json.loads((bundle / pointer["manifest_path"]).read_bytes())
            info = release["species"]["human"]
            target = json.loads((bundle / next(iter(info["target_paths"].values()))).read_bytes())["data"]
            original = target["reference_sequence_path"]
            extra = f"releases/{release['release_id']}/human/sequences/{'f' * 64}.json"
            content = (bundle / original).read_bytes()
            (bundle / extra).write_bytes(content)
            release["files"][extra] = copy.deepcopy(release["files"][original])
            fallback = info["records_download"]
            records = [json.loads(line) for line in gzip.decompress((bundle / fallback["path"]).read_bytes()).splitlines()]
            records.append({"kind": "sequence", "path": extra, "payload": json.loads(content)})
            content = gzip.compress(b"".join(canonical_bytes(row) + b"\n" for row in records), mtime=0)
            (bundle / fallback["path"]).write_bytes(content)
            metadata = {"sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)}
            release["files"][fallback["path"]].update(metadata)
            for download in [*info["downloads"], fallback]:
                if download["name"] == "records.jsonl.gz":
                    download.update({**metadata, "row_count": len(records)})
            body = canonical_bytes(release) + b"\n"
            (bundle / pointer["manifest_path"]).write_bytes(body)
            pointer["sha256"] = hashlib.sha256(body).hexdigest()
            body = canonical_bytes(root) + b"\n"
            (bundle / "manifest.json").write_bytes(body)
            (bundle / "manifest.sha256").write_text(hashlib.sha256(body).hexdigest() + "  manifest.json\n")
            with self.assertRaisesRegex(PublicDataError, "unreachable"):
                validate_public_api_bundle(bundle)

    def test_disease_source_spans_preserve_quoted_multiline_rows_and_detect_changes(self) -> None:
        import csv
        from agdesign2.public_api_data import _DiseaseRows
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "associations.tsv"
            fields = ["entry_name", "disease_id", "disease_name", "score"]
            rows = [
                {"entry_name": "FIRST_HUMAN", "disease_id": "EFO_000001", "disease_name": "Synthetic condition\nwith a quoted continuation", "score": "0.5"},
                {"entry_name": "FIRST_HUMAN", "disease_id": "EFO_000002", "disease_name": "Synthetic\tcondition", "score": "0.0"},
                {"entry_name": "SECOND_HUMAN", "disease_id": "EFO_000003", "disease_name": "Synthetic condition", "score": ""},
            ]
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
                handle.write("\n")
            indexed = _DiseaseRows(path)
            self.assertEqual(indexed.sha256, hashlib.sha256(path.read_bytes()).hexdigest())
            self.assertEqual(indexed.get("FIRST_HUMAN"), rows[:2])
            self.assertEqual(indexed.get("SECOND_HUMAN"), rows[2:])
            self.assertEqual(indexed.get("ABSENT", []), [])
            before = path.stat()
            content = path.read_bytes().replace(b"EFO_000001", b"EFO_000009")
            path.write_bytes(content)
            os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns))
            with self.assertRaisesRegex(PublicDataError, "rows changed"):
                indexed.get("FIRST_HUMAN")

    def test_disease_source_rejects_noncontiguous_target_groups(self) -> None:
        from agdesign2.public_api_data import _DiseaseRows
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "associations.tsv"
            path.write_text("entry_name\tdisease_id\nFIRST\tEFO_000001\nSECOND\tEFO_000002\nFIRST\tEFO_000003\n")
            with self.assertRaisesRegex(PublicDataError, "grouped by target"):
                _DiseaseRows(path)

    def test_compact_search_keeps_recorded_sequence_and_isoform_accessions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            human, mouse, site = source_fixture(root)
            human_report_path = human.parent / "report.json"
            human_report = json.loads(human_report_path.read_bytes())
            human_report["target"]["canonical_isoform_id"] = "P00001-1"
            human_report_path.write_text(json.dumps(human_report))
            mouse_report_path = mouse.parent / "report.json"
            originals = {"human": human_report_path.read_bytes(), "mouse": mouse_report_path.read_bytes()}
            build_public_api_bundle(human, mouse, root / "bundle", public_site_dir=site, software_commit="0" * 40, snapshot_id="recorded-search-accessions")
            release = self.release(root / "bundle")
            for species, queries in (("human", ("P00001", "P00001-1")), ("mouse", ("P00002", "NP_000001.1"))):
                info = release["species"][species]
                index = json.loads((root / "bundle" / info["index_path"]).read_bytes())["data"]
                self.assertEqual(len(index), 1)
                for query in queries:
                    self.assertIn(query.lower(), index[0]["search_text"])
                target = json.loads((root / "bundle" / index[0]["summary_path"]).read_bytes())["data"]
                summary = target["summary"]
                self.assertEqual(summary["search_text"], index[0]["search_text"])
                self.assertEqual(summary["identifier_accession"], "P00001" if species == "human" else "P00002")
                self.assertEqual(summary["sequence_accession"], "P00001" if species == "human" else "NP_000001.1")
                reference = json.loads((root / "bundle" / target["reference_sequence_path"]).read_bytes())["data"]
                original_sequence = json.loads(originals[species])["target"]["sequence"]
                self.assertEqual(reference["sequence"], original_sequence)
                self.assertEqual(reference["sequence_sha256"], hashlib.sha256(original_sequence.encode("ascii")).hexdigest())
                self.assertEqual(release["provenance"]["input_sha256"][f"{species}:report:0"], hashlib.sha256(originals[species]).hexdigest())
                fallback = [json.loads(line) for line in gzip.decompress((root / "bundle" / info["records_download"]["path"]).read_bytes()).splitlines()]
                stored_reference = next(row["payload"]["data"] for row in fallback if row["path"] == target["reference_sequence_path"])
                self.assertEqual(stored_reference, reference)
            self.assertEqual(human_report_path.read_bytes(), originals["human"])
            self.assertEqual(mouse_report_path.read_bytes(), originals["mouse"])

    def test_disease_projection_matches_first_100_source_rows_with_explicit_coverage(self) -> None:
        import csv
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            human, mouse, site = source_fixture(root)
            sources = {path: path.read_bytes() for path in (human.parent / "report.json", mouse.parent / "report.json")}
            rows = [{"entry_name": "SYNTHETIC_HUMAN", "disease_id": f"EFO_SYNTHETIC_{n:04d}", "disease_name": f"Synthetic association {n}", "score": ((n * 17 + 90) % 101) / 100, "direct_score": "", "indirect_score": ((n * 17 + 90) % 101) / 100, "open_targets_url": f"https://platform.opentargets.org/disease/EFO_SYNTHETIC_{n:04d}"} for n in range(107)]
            source = site / "open_targets_disease_associations.tsv"
            with source.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, tuple(rows[0]), delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            original = source.read_bytes()
            for species in ("human", "mouse"):
                metadata = site / ("mouse" if species == "mouse" else "") / "portal_metadata.json"
                metadata.write_text(json.dumps({"target_count": 1, "software_version": "fixture-version", "portal_build_date": "2026-10-16"}))
            build_public_api_bundle(human, mouse, root / "bundle", public_site_dir=site, software_commit="0" * 40, snapshot_id="portal-top-100")
            release = self.release(root / "bundle")
            info = release["species"]["human"]
            target = json.loads((root / "bundle" / next(iter(info["target_paths"].values()))).read_bytes())["data"]
            coverage = target["summary"]["disease_coverage"]
            self.assertEqual(coverage, {"state": "available", "selection": "portal_top_100", "source_count": 107, "projected_count": 100, "limit": 100, "source_reference": {"state": "live", "url": "https://openantigens.org/open_targets_disease_associations.tsv", "sha256": hashlib.sha256(original).hexdigest()}})
            evidence = []
            for relative in target["sections"]["disease"]["pages"]:
                page = json.loads((root / "bundle" / relative).read_bytes())["data"]
                self.assertEqual(page["coverage"], coverage)
                evidence.extend(page["items"])
            self.assertEqual([row["disease_id"] for row in evidence], [row["disease_id"] for row in rows[:100]])
            self.assertEqual([row["score"] for row in evidence], [row["score"] for row in rows[:100]])
            self.assertEqual(target["summary"]["disease_ids"], sorted(row["disease_id"] for row in rows[:100]))
            self.assertEqual(release["provenance"]["input_sha256"]["human:disease_tsv"], hashlib.sha256(original).hexdigest())
            self.assertEqual(release["provenance"]["projection_methods"]["disease"], {"selection": "portal_top_100", "species": "human", "limit": 100, "order": "published_source_order"})
            self.assertEqual(release["provenance"]["source_versions"]["human_portal_software"], "fixture-version")
            self.assertIsNone(release["provenance"]["source_versions"]["open_targets"])
            self.assertEqual(release["provenance"]["source_dates"]["human_portal_build_date"], "2026-10-16")
            self.assertEqual(source.read_bytes(), original)
            for path, content in sources.items():
                self.assertEqual(path.read_bytes(), content)
            self.assertFalse(any("open_targets_disease_associations.tsv" in path for path in release["files"]))

    def test_unknown_disease_query_is_preserved_and_mouse_coverage_is_not_applicable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            human, mouse, site = source_fixture(root)
            build_public_api_bundle(human, mouse, root / "bundle", public_site_dir=site, software_commit="0" * 40, snapshot_id="unknown-disease-query")
            release = self.release(root / "bundle")
            for species, info in release["species"].items():
                target = json.loads((root / "bundle" / next(iter(info["target_paths"].values()))).read_bytes())["data"]
                coverage = target["summary"]["disease_coverage"]
                self.assertEqual(coverage["state"], "unknown" if species == "human" else "not_applicable")
                self.assertIsNone(coverage["source_count"])
                self.assertEqual(coverage["projected_count"], 0)
                self.assertEqual(coverage["source_reference"], {"state": "unknown" if species == "human" else "not_applicable", "url": None, "sha256": None})
            (site / "open_targets_disease_associations.tsv").write_text("entry_name\tdisease_id\tdisease_name\tscore\n")
            build_public_api_bundle(human, mouse, root / "empty-source-bundle", public_site_dir=site, software_commit="0" * 40, snapshot_id="recorded-empty-source")
            empty = self.release(root / "empty-source-bundle")
            info = empty["species"]["human"]
            target = json.loads((root / "empty-source-bundle" / next(iter(info["target_paths"].values()))).read_bytes())["data"]
            self.assertEqual(target["summary"]["disease_coverage"]["source_count"], 0)
            self.assertEqual(target["summary"]["disease_coverage"]["state"], "unknown")
            self.assertEqual(target["sections"]["disease"]["state"], "unknown")

    def test_static_bulk_limit_reports_role_bytes_and_preserves_source_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            human, mouse, site = source_fixture(root)
            metadata_paths = [site / "portal_metadata.json", site / "mouse/portal_metadata.json"]
            original = [path.read_bytes() for path in metadata_paths]
            with mock.patch("agdesign2.public_api_data.EXTERNAL_DOWNLOAD_LIMIT", 256):
                with self.assertRaisesRegex(PublicDataError, r"Public download role index\.jsonl has [0-9]+ bytes; storage limit is 256 bytes") as error:
                    build_public_api_bundle(human, mouse, root / "bundle", public_site_dir=site, software_commit="0" * 40, snapshot_id="static-bulk-budget")
            self.assertNotIn(str(root), str(error.exception))
            self.assertEqual([path.read_bytes() for path in metadata_paths], original)

    def test_schema_metadata_keeps_one_mib_limit_despite_external_bulk_allowance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bundle = self.build(Path(tmp))
            root = json.loads((bundle / "manifest.json").read_bytes())
            pointer = root["releases"][0]
            release = json.loads((bundle / pointer["manifest_path"]).read_bytes())
            relative = next(row["path"] for row in release["species"]["human"]["downloads"] if row["name"] == "schema.json")
            content = (bundle / relative).read_bytes().ljust(1024 * 1024 + 1, b" ")
            (bundle / relative).write_bytes(content)
            metadata = {"sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)}
            release["files"][relative].update(metadata)
            for download in release["species"]["human"]["downloads"]:
                if download["name"] == "schema.json":
                    download.update(metadata)
            body = canonical_bytes(release) + b"\n"
            (bundle / pointer["manifest_path"]).write_bytes(body)
            pointer["sha256"] = hashlib.sha256(body).hexdigest()
            body = canonical_bytes(root) + b"\n"
            (bundle / "manifest.json").write_bytes(body)
            (bundle / "manifest.sha256").write_text(hashlib.sha256(body).hexdigest() + "  manifest.json\n")
            with self.assertRaisesRegex(PublicDataError, r"schema\.json has 1048577 bytes; storage limit is 1048576 bytes"):
                validate_public_api_bundle(bundle)


if __name__ == "__main__":
    unittest.main()
