"""Build offline human/mouse portals and their public API bundle from synthetic proteins."""
from __future__ import annotations

import json
import csv
from pathlib import Path
import sys

from agdesign2.portal import build_portal
from agdesign2.public_api_data import build_public_api_bundle


def build(root: Path, *, warning_heavy: bool = False, remote: bool = False) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    summaries: dict[str, Path] = {}
    for species, taxon, accession in (("human", 9606, "P00001"), ("mouse", 10090, "P00002")):
        batch = root / species
        batch.mkdir()
        sequence = "ACDEFGHIKLMNPQRSTVWY" * 6
        if species == "mouse":
            sequence = "M" + sequence[1:]
        entry_name = f"SYNTHETIC_{species.upper()}"
        report = {
            "target": {
                "accession": accession, "entry_name": entry_name,
                "gene_symbol": "SYNTHETIC", "protein_name": "Synthetic integration fixture",
                "taxon_id": taxon, "organism": "Homo sapiens" if species == "human" else "Mus musculus",
                "sequence": sequence, "gene_synonyms": ["FIXTURE_ALIAS"],
                **({"canonical_isoform_id": "P00001-1"} if species == "human" else {}),
                "identifier_source": "human_report_mouse_homolog" if species == "mouse" else "uniprot",
                "sequence_source": "ortholog_reference_table.mouse_refseq_sequence" if species == "mouse" else "uniprot",
                **({"refseq_accession": "NP_000001.1"} if species == "mouse" else {}),
            },
            "topology": {"topology_class": "secreted", "source": "synthetic_fixture"},
            "ectodomain": {"start": 25, "end": 90, "label": "Secreted region", "source": "synthetic_fixture"},
            "construct_details": [
                {"name": name, "start": start, "end": end, "length": end - start + 1,
                 "sequence": sequence[start - 1:end], "classification": "soluble_ecd",
                 "score": 0.5, "rationale": "Synthetic fixture; no experimental claim.",
                 "warnings": ["Synthetic record for software tests."], "evidence": ["synthetic_fixture"]}
                for name, start, end in (("full_ectodomain", 25, 90), ("strict_fixture", 25, 44))
            ],
            "construct_recommendations": [], "cross_reactivity_hits": [],
            "family_context": {"source": "synthetic_fixture", "family_names": ["Synthetic family"], "members": []},
            "assembly_requirements": [], "experimental_constructs": [], "interpro_annotations": [],
            "ptms": [{"type": "CARBOHYD", "start": 28, "end": 28, "description": "Synthetic annotation."}],
        }
        report_path = batch / f"synthetic_{species}_report.json"
        if warning_heavy:
            report["construct_details"][0]["warnings"] = ["!" * 12000 for _ in range(40)]
        report_path.write_text(json.dumps(report), encoding="utf-8")
        summary_path = batch / "batch_summary.json"
        summary_path.write_text(json.dumps([{
            "query": entry_name, "resolved_entry_name": entry_name, "status": "ok", "batch_index": 1,
            "json_report": str(report_path), "topology_bucket": "Secreted", "pubtator_count": 0,
            **({"mouse_refseq_accession": "NP_000001.1", "mouse_sequence_source": "ortholog_reference_table.mouse_refseq_sequence"} if species == "mouse" else {}),
        }]), encoding="utf-8")
        summaries[species] = summary_path
        if species == "human":
            associations = [{"entry_name": entry_name, "gene_symbol": "SYNTHETIC", "uniprot_accession": accession, "disease_id": f"EFO_SYNTHETIC_{number:04d}", "disease_name": f"Synthetic association {number}", "score": ((number * 17 + 90) % 101) / 100, "direct_score": "", "indirect_score": ((number * 17 + 90) % 101) / 100, "open_targets_url": f"https://platform.opentargets.org/disease/EFO_SYNTHETIC_{number:04d}"} for number in range(105)]
            with (batch / "open_targets_disease_associations.tsv").open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, tuple(associations[0]), delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(associations)
            (batch / "open_targets_disease_associations.json").write_text(json.dumps({"targets": [{"entry_name": entry_name, "gene_symbol": "SYNTHETIC", "uniprot_accession": accession, "association_count_downloaded": 105, "top_disease_id": associations[0]["disease_id"], "top_disease_name": associations[0]["disease_name"], "top_disease_score": associations[0]["score"]}], "associations": associations}), encoding="utf-8")
        build_portal(
            summary_path, output_dir=root / "site" / ("mouse" if species == "mouse" else ""),
            refresh_reports=False, refresh_stale_reports=False, fetch_literature=False,
            portal_title="OpenAntigens Mouse" if species == "mouse" else "OpenAntigens",
            theme=species, include_disease_context=species == "human",
            sibling_link=("OpenAntigens Human", "../index.html") if species == "mouse" else ("OpenAntigens Mouse", "mouse/index.html"),
            detail_sibling_link=("OpenAntigens Human", "../../index.html") if species == "mouse" else ("OpenAntigens Mouse", "../mouse/index.html"),
        )
    manifest = build_public_api_bundle(
        summaries["human"], summaries["mouse"], root / "bundle",
        public_site_dir=root / "site", software_commit="0" * 40, snapshot_id="synthetic-integration-fixture",
    )
    from agdesign2.public_api_transport import build_runtime_data

    build_runtime_data(root / "bundle", root / "runtime")
    if remote:
        from agdesign2.public_api_query import build_query_runtime
        result = build_query_runtime(
            root / "bundle", root / "remote-runtime", root / "remote-static",
            schema_path=Path(__file__).resolve().parents[1] / "services/public-api/src/query-schema.json",
        )
        (root / "remote-seal.json").write_text(json.dumps(result["seal"], sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def build_remote(root: Path, *, warning_heavy: bool = False) -> Path:
    return build(root, warning_heavy=warning_heavy, remote=True)


if __name__ == "__main__":
    print(build(Path(sys.argv[1]).resolve(), warning_heavy="--warning-heavy" in sys.argv[2:], remote="--remote" in sys.argv[2:]))
