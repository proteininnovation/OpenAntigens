from __future__ import annotations

import json
import sqlite3
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.af3 import catalogue_af3_structures, copy_af3_catalog_artifacts, find_local_af3_artifacts, prepare_af3_structures
from agdesign2.alphafold import AlphaFoldClient
from agdesign2.config import AnalysisConfig
from agdesign2.exceptions import ExternalServiceError
from agdesign2.portal import portal_index_tsv
from agdesign2.structure_utils import load_pae_matrix, parse_alphafold_pdb


class AlwaysFailingHttpClient:
    def fetch_json(self, url: str, *, headers=None, cache_namespace="json", cache_key=None):
        raise ExternalServiceError("metadata unavailable")

    def download(self, url: str, destination: Path, *, headers=None) -> Path:
        raise ExternalServiceError("download unavailable")


class AF3PrepareTests(unittest.TestCase):
    def test_af3_fallback_is_disabled_by_default(self) -> None:
        self.assertFalse(AnalysisConfig().enable_local_af3_fallback)

    def test_prepare_af3_structures_selects_best_model_and_writes_compatible_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            snapshot = root / "snapshot"
            data_dir = snapshot / "data"
            data_dir.mkdir(parents=True)
            report_path = data_dir / "test_human_report.json"
            report_path.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "P12345",
                            "entry_name": "TEST_HUMAN",
                            "gene_symbol": "TEST",
                            "sequence": "MG",
                        }
                    }
                ),
                encoding="utf-8",
            )
            (data_dir / "batch_summary.json").write_text(
                json.dumps([{"query": "TEST_HUMAN", "json_report": str(report_path), "status": "ok"}]),
                encoding="utf-8",
            )
            zip_path = root / "af3.zip"
            with zipfile.ZipFile(zip_path, "w") as zf:
                base = "openantigen_af3_test_p12345"
                zf.writestr("terms_of_use.md", "AlphaFold Server terms")
                zf.writestr(
                    f"{base}/fold_openantigen_af3_test_p12345_job_request.json",
                    json.dumps(
                        [
                            {
                                "name": "OpenAntigen_AF3_TEST_P12345",
                                "sequences": [{"proteinChain": {"sequence": "MG", "count": 1}}],
                            }
                        ]
                    ),
                )
                for index, score in ((0, 0.2), (1, 0.9)):
                    zf.writestr(f"{base}/fold_openantigen_af3_test_p12345_model_{index}.cif", _mini_cif())
                    zf.writestr(
                        f"{base}/fold_openantigen_af3_test_p12345_summary_confidences_{index}.json",
                        json.dumps({"ranking_score": score, "ptm": 0.1 + index, "fraction_disordered": 0.0}),
                    )
                    zf.writestr(
                        f"{base}/fold_openantigen_af3_test_p12345_full_data_{index}.json",
                        json.dumps({"pae": [[0.0, 1.0], [1.0, 0.0]]}),
                    )

            result = prepare_af3_structures(zip_path, snapshot_dir=snapshot)
            self.assertEqual(len(result.rows), 1)
            self.assertEqual(result.rows[0]["status"], "ok")
            self.assertEqual(result.rows[0]["selected_model_index"], 1)
            pdb_path = result.output_dir / "P12345.pdb"
            pae_path = result.output_dir / "P12345.pae.json"
            meta_path = result.output_dir / "P12345.meta.json"
            self.assertTrue(pdb_path.exists())
            self.assertTrue(pae_path.exists())
            self.assertTrue(meta_path.exists())
            self.assertEqual(_sequence_from_pdb(pdb_path), "MG")
            self.assertEqual(load_pae_matrix(pae_path), [[0.0, 1.0], [1.0, 0.0]])
            self.assertIsNotNone(find_local_af3_artifacts("P12345", data_dir=data_dir, canonical_sequence="MG"))

    def test_alphafold_client_uses_local_af3_after_af2_failure(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            data_dir = Path(tmpdir)
            af3_dir = data_dir / "alphafold3"
            af3_dir.mkdir()
            pdb_path = af3_dir / "P12345.pdb"
            pdb_path.write_text(_mini_pdb(), encoding="utf-8")
            (af3_dir / "P12345.pae.json").write_text('{"pae":[[0.0,1.0],[1.0,0.0]]}', encoding="utf-8")
            (af3_dir / "P12345.meta.json").write_text(
                json.dumps(
                    {
                        "source": "alphafold_server",
                        "model_type": "AF3",
                        "sequence_sha256": "6f496e9738ca13d116978af2d9a3b15922fd01087585ba27f0aae705fa77bc65",
                    }
                ),
                encoding="utf-8",
            )
            client = AlphaFoldClient(
                AlwaysFailingHttpClient(),
                AnalysisConfig(data_dir=data_dir, enable_local_af3_fallback=True),
            )
            pdb, pae = client.ensure_artifacts("P12345", canonical_sequence="MG")
            self.assertEqual(pdb, af3_dir / "P12345.pdb")
            self.assertEqual(pae, af3_dir / "P12345.pae.json")

    def test_catalogue_af3_structures_writes_sqlite_and_selects_best_duplicate(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            zip_a = root / "af3_a.zip"
            zip_b = root / "af3_b.zip"
            _write_catalogue_zip(zip_a, accession="P12345", score=0.4)
            _write_catalogue_zip(zip_b, accession="P12345", score=0.9)

            result = catalogue_af3_structures(root, output_dir=root / "catalog")

            self.assertTrue(result.db_path.exists())
            self.assertTrue(result.manifest_tsv.exists())
            self.assertTrue(result.manifest_json.exists())
            self.assertTrue((result.output_dir / "terms_of_use" / "terms_of_use.md").exists())
            self.assertTrue((result.output_dir / "terms_of_use" / "manifest.json").exists())
            selected = [row for row in result.rows if row.get("selected_for_accession")]
            self.assertEqual(len(selected), 1)
            self.assertEqual(selected[0]["source_zip"], str(zip_b))
            artifact_dir = result.output_dir / "alphafold3"
            self.assertTrue((artifact_dir / "P12345.pdb").exists())
            self.assertTrue((artifact_dir / "P12345.pae.json.gz").exists())
            self.assertTrue((artifact_dir / "P12345.meta.json").exists())
            with sqlite3.connect(result.db_path) as conn:
                target_count = conn.execute("SELECT COUNT(*) FROM targets").fetchone()[0]
                selected_count = conn.execute("SELECT COUNT(*) FROM targets WHERE selected_for_accession = 1").fetchone()[0]
                member_count = conn.execute("SELECT COUNT(*) FROM archive_members").fetchone()[0]
            self.assertEqual(target_count, 2)
            self.assertEqual(selected_count, 1)
            self.assertGreater(member_count, 0)

    def test_copy_af3_catalog_artifacts_filters_to_target_tsv(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            _write_catalogue_zip(root / "af3_a.zip", accession="P12345", score=0.9)
            _write_catalogue_zip(root / "af3_b.zip", accession="Q99999", score=0.8)
            result = catalogue_af3_structures(root, output_dir=root / "catalog")
            target_tsv = root / "targets.tsv"
            target_tsv.write_text("uniprot_id\tuniprot_name\nP12345\tTEST_HUMAN\n", encoding="utf-8")

            manifest = copy_af3_catalog_artifacts(result.db_path, data_dir=root / "snapshot_data", target_tsv=target_tsv)

            af3_dir = root / "snapshot_data" / "alphafold3"
            self.assertEqual(manifest["copied_accessions"], ["P12345"])
            self.assertTrue((af3_dir / "P12345.pdb").exists())
            self.assertTrue((af3_dir / "P12345.pae.json.gz").exists())
            self.assertFalse((af3_dir / "Q99999.pdb").exists())

    def test_portal_download_index_includes_structure_source(self) -> None:
        tsv = portal_index_tsv(
            [
                {
                    "entry_name": "TEST_HUMAN",
                    "status": "ok",
                    "has_alphafold_structure": True,
                    "structure_source": "AlphaFold 3 model",
                }
            ]
        )

        header, row = tsv.splitlines()[:2]
        self.assertIn("structure_source", header.split("\t"))
        self.assertIn("AlphaFold 3 model", row)


def _sequence_from_pdb(path: Path) -> str:
    names = parse_alphafold_pdb(path)["residue_names"]
    mapping = {"MET": "M", "GLY": "G"}
    return "".join(mapping[names[index]] for index in sorted(names))


def _mini_pdb() -> str:
    return (
        "ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00 90.00           N  \n"
        "ATOM      2  CA  MET A   1       1.000   0.000   0.000  1.00 91.00           C  \n"
        "ATOM      3  N   GLY A   2       2.000   0.000   0.000  1.00 80.00           N  \n"
        "ATOM      4  CA  GLY A   2       3.000   0.000   0.000  1.00 81.00           C  \n"
        "END\n"
    )


def _mini_cif() -> str:
    return """data_test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . MET A 1 1 ? 0.0 0.0 0.0 1.00 90.0 ? 1 MET A N 1
ATOM 2 C CA . MET A 1 1 ? 1.0 0.0 0.0 1.00 91.0 ? 1 MET A CA 1
ATOM 3 N N . GLY A 1 2 ? 2.0 0.0 0.0 1.00 80.0 ? 2 GLY A N 1
ATOM 4 C CA . GLY A 1 2 ? 3.0 0.0 0.0 1.00 81.0 ? 2 GLY A CA 1
#
"""


def _write_catalogue_zip(path: Path, *, accession: str, score: float) -> None:
    base = f"openantigen_af3_test_{accession.lower()}"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("terms_of_use.md", "AlphaFold Server terms")
        zf.writestr(
            f"{base}/fold_openantigen_af3_test_{accession.lower()}_job_request.json",
            json.dumps(
                [
                    {
                        "name": f"OpenAntigen_AF3_TEST_{accession}",
                        "sequences": [{"proteinChain": {"sequence": "MG", "count": 1}}],
                    }
                ]
            ),
        )
        zf.writestr(f"{base}/fold_openantigen_af3_test_{accession.lower()}_model_0.cif", _mini_cif())
        zf.writestr(
            f"{base}/fold_openantigen_af3_test_{accession.lower()}_summary_confidences_0.json",
            json.dumps({"ranking_score": score, "ptm": 0.1, "fraction_disordered": 0.0}),
        )
        zf.writestr(
            f"{base}/fold_openantigen_af3_test_{accession.lower()}_full_data_0.json",
            json.dumps({"pae": [[0.0, 1.0], [1.0, 0.0]]}),
        )


if __name__ == "__main__":
    unittest.main()
