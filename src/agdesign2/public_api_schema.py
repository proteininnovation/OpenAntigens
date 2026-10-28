"""Closed schemas for the public release projection; all validation is offline."""
from __future__ import annotations

from typing import Any

from jsonschema import Draft202012Validator

SCHEMA_VERSION = 1
SPECIES = {"human": 9606, "mouse": 10090}
SECTIONS = ("topology", "domains", "ptms", "orthologs", "paralogs", "family", "assembly", "structures", "similarity", "literature", "disease")
STATES = ("available", "no_results", "unavailable", "not_applicable", "error", "unknown")
AVAILABILITY = ("available", "unavailable", "unknown")
CATEGORIES = ("full_region", "pdb", "domain_repeat", "strict", "lenient", "multipass_full_length", "multipass_trimmed")
STATUSES = ("ok", "error", "pending", "running", "started", "skipped_existing", "unknown")
TRACKS = ("secreted", "gpi", "single_pass", "multipass", "other", "unknown")


def obj(fields: dict[str, Any]) -> dict[str, Any]:
    return {"type": "object", "properties": fields, "required": list(fields), "additionalProperties": False}


def array(item: dict[str, Any]) -> dict[str, Any]:
    return {"type": "array", "items": item}


def nullable(schema: dict[str, Any]) -> dict[str, Any]:
    return {"anyOf": [schema, {"type": "null"}]}


def ref(name: str) -> dict[str, Any]:
    return {"$ref": f"#/$defs/{name}"}


def mapping(item: dict[str, Any]) -> dict[str, Any]:
    return {"type": "object", "additionalProperties": item}


def enum(values: tuple[str, ...]) -> dict[str, Any]:
    return {"type": "string", "enum": list(values)}


TEXT = {"type": "string", "maxLength": 16384}
ID = {"type": "string", "maxLength": 256}
PATH = {"type": "string", "maxLength": 1024}
NUMBER = {"type": "number"}
INTEGER = {"type": "integer"}
COUNT = {"type": "integer", "minimum": 0}
BOOL = {"type": "boolean"}
HASH = {"type": "string", "pattern": "^[a-f0-9]{64}$"}
N = nullable
S = N(TEXT)
I = N(INTEGER)
F = N(NUMBER)
STRINGS = array(TEXT)
STATE = enum(STATES)

DEFS: dict[str, dict[str, Any]] = {}
DEFS["region"] = obj({"start": I, "end": I, "label": S, "source": S, "confidence": F})
DEFS["design_region"] = obj({"start": INTEGER, "end": INTEGER, "label": S, "source": S})
DEFS["feature"] = obj({"type": S, "start": I, "end": I, "description": S, "source": S})
DEFS["furin"] = obj({"start": I, "end": I, "motif": S, "suggestions": STRINGS, "curated_annotation": S})
DEFS["cysteine"] = obj({"position": I, "paired_with": I, "surface_exposed": N(BOOL), "warning": S, "partner_outside_range": I, "closest_unpaired_cysteine": I, "closest_unpaired_cysteine_distance": F, "closest_unpaired_cysteine_distance_atom": S})
DEFS["cysteine"]["properties"].update({"analysis_start": I, "analysis_end": I})
DEFS["cysteine"]["required"] = list(DEFS["cysteine"]["properties"])
DEFS["processing_record"] = obj({"kind": enum(("feature", "furin", "cysteine")), "feature": N(ref("feature")), "furin": N(ref("furin")), "cysteine": N(ref("cysteine"))})
DEFS["domain"] = obj({"accession": S, "name": S, "type": S, "source_database": S, "start": I, "end": I, "representative": N(BOOL), "integrated_accession": S, "integrated_name": S})
DEFS["topology"] = obj({"topology_class": S, "source": S, **{key: array(ref("region")) for key in ("extracellular_regions", "transmembrane_regions", "intramembrane_regions", "cytoplasmic_regions")}, **{key: N(ref("region")) for key in ("major_extracellular_region", "membrane_core_region", "cytoplasmic_tail_region")}})
DEFS["ortholog"] = obj({"species": S, "accession": S, "entry_name": S, "ectodomain_start": I, "ectodomain_end": I, "identity": F, "coverage": F, "extracellular_surface_identity": F, "extracellular_surface_aligned_positions": I, "available": N(BOOL), "notes": STRINGS, "selection_status": S, "candidates": S, "evidence_count": I, "accession_source": S})
DEFS["family_member"] = obj({"gene_symbol": S, "gene_name": S, "accession": S, "entry_name": S, "sequence_length": I, "ectodomain_start": I, "ectodomain_end": I, "ectodomain_length": I, "notes": STRINGS})
DEFS["family"] = obj({"gene_symbol": S, "source": S, "family_names": STRINGS, "members": array(ref("family_member"))})
DEFS["assembly"] = obj({"summary": S, "obligatory": N(BOOL), "confidence": S, "source": S, "classification": S, "partners": STRINGS, "rationale": STRINGS, "complex_portal_ids": STRINGS})
DEFS["complex_participant"] = obj({"identifier": S, "name": S, "description": S, "stoichiometry": S, "bio_role": S, "interactor_type": S, "binding_regions": STRINGS})
DEFS["complex"] = obj({"complex_ac": S, "name": S, "species": S, "predicted_complex": N(BOOL), "evidence_code": S, "evidence_description": S, "confidence_score": I, "complex_assemblies": STRINGS, "properties": STRINGS, "participants": array(ref("complex_participant"))})
DEFS["complex_lookup"] = obj({"species": S, "status": S, "source": S, "reason": S})
DEFS["assembly_record"] = obj({"kind": enum(("requirement", "complex", "lookup")), "requirement": N(ref("assembly")), "complex": N(ref("complex")), "lookup": N(ref("complex_lookup"))})
DEFS["structure"] = obj({"kind": enum(("pdb", "predicted_model")), "pdb_id": S, "method": S, "resolution": S, "chains": array(ref("region")), "source": S, "live_asset_url": S, "sha256": N(HASH)})
DEFS["similarity"] = obj({"subject_id": S, "description": S, "species": S, "identity": F, "coverage": F, "alignment_length": I, "evalue": F, "bitscore": F, "query_start": I, "query_end": I, "subject_start": I, "subject_end": I, "alignment_source": S, "extracellular_surface_identity": F, "extracellular_surface_aligned_positions": I, "extracellular_surface_matches": I, "scope": enum(("design_region", "full_length"))})
DEFS["literature"] = obj({"count": I, "query_url": S, "source": S})
DEFS["disease"] = obj({"disease_id": S, "disease_name": S, "score": F, "direct_score": F, "indirect_score": F, "open_targets_url": S, "source": S})
DEFS["disease_coverage"] = obj({"state": enum(("available", "unknown", "not_applicable")), "selection": {"const": "portal_top_100"}, "source_count": N(COUNT), "projected_count": {"type": "integer", "minimum": 0, "maximum": 100}, "limit": {"const": 100}, "source_reference": obj({"state": enum(("live", "unknown", "not_applicable")), "url": S, "sha256": N(HASH)})})
DEFS["disease_coverage"]["properties"]["source_count"]["description"] = "Count of rows for this target in the recorded source file. Zero does not prove that a source query ran."
DEFS["metrics"] = obj({key: F for key in ("mean_plddt", "min_plddt", "max_plddt", "structured_fraction", "mean_intra_pae", "max_intra_pae")})
DEFS["hazards"] = obj({"cysteines": array(ref("cysteine")), "furin_sites": array(ref("furin")), "ptms": array(ref("feature")), "ligand_interactions": array(ref("feature")), "complete_domains": array(ref("domain")), "clipped_domains": array(ref("domain"))})
DEFS["summary"] = obj({"target_id": N(ID), "record_id": HASH, "species": enum(tuple(SPECIES)), "taxon_id": {"type": "integer", "enum": list(SPECIES.values())}, "identifier_source": N(enum(("uniprot", "refseq"))), "identifier_accession": S, "entry_name": S, "gene_symbol": S, "protein_name": S, "aliases": STRINGS, "sequence_accession": S, "sequence_source": S, "sequence_sha256": N(HASH), "sequence_length": N(COUNT), "status": enum(STATUSES), "status_raw": S, "track": enum(TRACKS), "track_raw": S, "analyzed_topology": S, "analyzed_topology_source": S, "design_region": N(ref("design_region")), "construct_count": N(COUNT), "structure_state": enum(AVAILABILITY), "construct_state": enum(AVAILABILITY), "section_availability": obj({key: STATE for key in SECTIONS}), "disease_ids": STRINGS, "family_names": STRINGS, "live_report_url": S, "search_text": TEXT})
DEFS["summary"]["properties"].update({"canonical_isoform_id": S, "portal_construct_count": N(COUNT)})
DEFS["summary"]["properties"]["unsupported_variants"] = array(obj({"name": S, "classification": S, "state": {"const": "unavailable"}, "reason": TEXT}))
DEFS["summary"]["properties"]["disease_coverage"] = ref("disease_coverage")
DEFS["summary"]["required"] = list(DEFS["summary"]["properties"])
SEARCH_FIELDS = ("target_id", "record_id", "species", "taxon_id", "search_text", "track", "status", "structure_state", "construct_state")
DEFS["search_row"] = obj({**{key: DEFS["summary"]["properties"][key] for key in SEARCH_FIELDS}, "summary_path": PATH})
DEFS["construct"] = obj({"construct_id": HASH, "target_id": ID, "name": S, "category": enum(CATEGORIES), "classification": S, "start": I, "end": I, "length": I, "sequence_basis": enum(("reference_slice", "recorded_variant", "unavailable")), "sequence_sha256": N(HASH), "sequence_accession": S, "sequence_source": S, "export_state": STATE, "export_reason": S, "score": F, "rationale": S, "warnings": STRINGS, "evidence": STRINGS, "pdb_id": S, "structural_metrics": N(ref("metrics"))})
DEFS["construct"]["properties"]["hazards"] = ref("hazards")
DEFS["construct"]["required"] = list(DEFS["construct"]["properties"])
DEFS["sequence"] = obj({"target_id": ID, "construct_id": N(HASH), "species": enum(tuple(SPECIES)), "taxon_id": INTEGER, "identifier_accession": S, "identifier_source": N(enum(("uniprot", "refseq"))), "sequence_accession": S, "sequence_source": S, "sequence_basis": enum(("reference", "reference_slice", "recorded_variant")), "start": INTEGER, "end": INTEGER, "length": COUNT, "sequence": {"type": "string", "maxLength": 1000000, "pattern": "^[A-Z*]+$"}, "sequence_sha256": HASH, "warnings": STRINGS, "evidence": STRINGS})
DEFS["descriptor"] = obj({"state": STATE, "reason": S, "total": COUNT, "pages": array(PATH)})
DEFS["target"] = obj({"summary": ref("summary"), "sections": obj({key: ref("descriptor") for key in SECTIONS}), "constructs": ref("descriptor"), "reference_sequence_path": N(PATH), "construct_sequence_paths": mapping(N(PATH))})
DEFS["author"] = obj({"given": S, "family": S, "orcid": S})
DEFS["citation"] = obj({"title": S, "authors": array(ref("author")), "doi": S, "version": S, "journal": S, "date": S, "type": S})
DEFS["provenance"] = obj({"software_commit": ID, "snapshot_id": ID, "source_versions": mapping(S), "source_dates": mapping(S), "input_sha256": mapping(HASH), "notes": STRINGS})
DEFS["provenance"]["properties"]["projection_methods"] = obj({"disease": obj({"selection": {"const": "portal_top_100"}, "species": {"const": "human"}, "limit": {"const": 100}, "order": {"const": "published_source_order"}})})
DEFS["file"] = obj({"sha256": HASH, "bytes": COUNT, "media_type": ID, "encoding": N(enum(("gzip",))), "schema": S})
DEFS["download"] = obj({"name": ID, "url": PATH, "path": PATH, "media_type": ID, "bytes": COUNT, "sha256": HASH, "encoding": N(enum(("gzip",))), "row_count": N(COUNT), "schema": S, "license": TEXT})
DEFS["species_manifest"] = obj({"taxon_id": INTEGER, "index_path": PATH, "target_paths": mapping(PATH), "downloads": array(ref("download")), "counts": obj({"targets": COUNT, "resolved_targets": COUNT, "constructs": COUNT, "status": mapping(COUNT)}), "capabilities": obj({"sections": array(enum(SECTIONS)), "disease": BOOL}), "records_download": ref("download")})
DEFS["core"] = obj({"schema_version": {"const": 1}, "provenance": ref("provenance"), "payload_sha256": mapping(HASH), "asset_sha256": mapping(HASH)})
DEFS["release_manifest"] = obj({"release_id": HASH, "schema_version": {"const": 1}, "manifest_url": PATH, "provenance": ref("provenance"), "citation": ref("citation"), "licenses": obj({"software": TEXT, "annotations": TEXT, "third_party": TEXT}), "species": obj({key: ref("species_manifest") for key in SPECIES}), "files": mapping(ref("file")), "core": ref("core")})
DEFS["root_manifest"] = obj({"schema_version": {"const": 1}, "current_release_id": HASH, "releases": array(obj({"release_id": HASH, "manifest_path": PATH, "sha256": HASH}))})
DEFS["retired_release"] = obj({"release_id": HASH, "retired_by_release_id": HASH, "archive_url": S})
DEFS["root_manifest"]["properties"]["retired_releases"] = {**array(ref("retired_release")), "maxItems": 1000}


def envelope(data: dict[str, Any]) -> dict[str, Any]:
    return obj({"release_id": HASH, "schema_version": {"const": 1}, "species": enum(tuple(SPECIES)), "data": data})


DEFS["index"] = envelope(array(ref("search_row")))
DEFS["target_shard"] = envelope(ref("target"))
DEFS["construct_page"] = envelope(array(ref("construct")))
DEFS["sequence_shard"] = envelope(ref("sequence"))
SECTION_TYPES = {"topology": "topology", "domains": "domain", "ptms": "processing_record", "orthologs": "ortholog", "paralogs": "family_member", "family": "family", "assembly": "assembly_record", "structures": "structure", "similarity": "similarity", "literature": "literature", "disease": "disease"}
for section, record_type in SECTION_TYPES.items():
    DEFS[f"section_{section}"] = envelope(obj({"target_id": ID, "section": {"const": section}, "state": STATE, "reason": S, "items": array(ref(record_type))}))
DEFS["section_disease"]["properties"]["data"]["properties"]["coverage"] = ref("disease_coverage")
DEFS["section_disease"]["properties"]["data"]["required"].append("coverage")

PUBLIC_SCHEMA = {"$schema": "https://json-schema.org/draft/2020-12/schema", "$id": "https://openantigens.org/schemas/public-api-v1.json", "x-openantigens-fasta-version": 2, "x-openantigens-index-version": 2, "x-openantigens-disease-selection": "portal_top_100", "$defs": DEFS}


class PublicDataError(ValueError):
    """A public release cannot be represented or verified safely."""


def validate_schema(name: str, value: Any) -> None:
    if name not in DEFS:
        raise PublicDataError("Unknown public schema")
    schema = {**PUBLIC_SCHEMA, "$ref": f"#/$defs/{name}"}
    error = next(Draft202012Validator(schema).iter_errors(value), None)
    if error is not None:
        path = ".".join(map(str, error.absolute_path))
        raise PublicDataError(f"Invalid public {name} field at {path or '<root>'}")
