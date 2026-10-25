"""Build offline human/mouse portals and their public API bundle from synthetic proteins."""
from __future__ import annotations

import json
from pathlib import Path
import sys

from agdesign2.portal import build_portal
from agdesign2.public_api_data import build_public_api_bundle


def build(root: Path) -> Path:
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
        report_path.write_text(json.dumps(report), encoding="utf-8")
        summary_path = batch / "batch_summary.json"
        summary_path.write_text(json.dumps([{
            "query": entry_name, "resolved_entry_name": entry_name, "status": "ok", "batch_index": 1,
            "json_report": str(report_path), "topology_bucket": "Secreted", "pubtator_count": 0,
            **({"mouse_refseq_accession": "NP_000001.1", "mouse_sequence_source": "ortholog_reference_table.mouse_refseq_sequence"} if species == "mouse" else {}),
        }]), encoding="utf-8")
        summaries[species] = summary_path
        build_portal(
            summary_path, output_dir=root / "site" / ("mouse" if species == "mouse" else ""),
            refresh_reports=False, refresh_stale_reports=False, fetch_literature=False,
            portal_title="OpenAntigens Mouse" if species == "mouse" else "OpenAntigens",
            theme=species, include_disease_context=species == "human",
            sibling_link=("OpenAntigens Human", "../index.html") if species == "mouse" else ("OpenAntigens Mouse", "mouse/index.html"),
            detail_sibling_link=("OpenAntigens Human", "../../index.html") if species == "mouse" else ("OpenAntigens Mouse", "../mouse/index.html"),
        )
    return build_public_api_bundle(
        summaries["human"], summaries["mouse"], root / "bundle",
        public_site_dir=root / "site", software_commit="0" * 40, snapshot_id="synthetic-integration-fixture",
    )


if __name__ == "__main__":
    print(build(Path(sys.argv[1]).resolve()))
