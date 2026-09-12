from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Callable

from .config import AnalysisConfig
from .models import AnalysisNote, CanonicalFamily, FamilyContext, FamilyMember, Region, TargetRecord
from .pipeline import AntigenAnalyzer
from .sequence_utils import slice_sequence


ModuleHandler = Callable[[AntigenAnalyzer, dict[str, Any]], bool]


def refresh_report_modules(
    summary_path: str | Path,
    *,
    modules: list[str],
    queries: list[str] | None = None,
    only_zero_cross_reactivity: bool = False,
    verbose: bool = False,
    analyzer: AntigenAnalyzer | None = None,
) -> tuple[Path, list[dict[str, Any]], int]:
    summary_file = Path(summary_path)
    batch_dir = summary_file.parent
    results = json.loads(summary_file.read_text(encoding="utf-8"))
    selected_queries = {query.strip() for query in (queries or []) if query and query.strip()}
    normalized_modules = [_normalize_module_name(name) for name in modules]
    tool = analyzer or AntigenAnalyzer(
        config=_module_refresh_config(
            batch_dir=batch_dir,
            verbose=verbose,
            enable_complex_portal="complex-portal" in normalized_modules,
        )
    )
    if normalized_modules == ["cross-reactivity"] and callable(getattr(tool.blast_client, "search_many", None)):
        return _refresh_cross_reactivity_in_bulk(
            summary_file=summary_file,
            batch_dir=batch_dir,
            results=results,
            selected_queries=selected_queries,
            only_zero_cross_reactivity=only_zero_cross_reactivity,
            verbose=verbose,
            analyzer=tool,
        )
    handlers = [_MODULE_HANDLERS[name] for name in normalized_modules]

    updated = 0
    for item in results:
        query = str(item.get("query") or "").strip()
        if selected_queries and query not in selected_queries and str(item.get("resolved_entry_name") or "").strip() not in selected_queries:
            continue
        report_path = _resolve_existing_path(item.get("json_report"), batch_dir)
        if report_path is None or not report_path.exists():
            continue
        try:
            payload = json.loads(report_path.read_text(encoding="utf-8"))
        except Exception as exc:
            item["status"] = "error"
            item["error"] = str(exc)
            summary_file.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
            continue
        if only_zero_cross_reactivity and "cross-reactivity" in normalized_modules:
            if not _has_zero_relevant_cross_reactivity(payload):
                continue
        if verbose:
            print(f"[agdesign2-modules] {query or report_path.name}: {', '.join(normalized_modules)}", flush=True)
        try:
            # Determine "changed" from the actual serialized payload, not the
            # handler return values (several handlers unconditionally return
            # True, which rewrote every report and re-serialized the whole
            # summary per item — O(reports^2) writes even when nothing changed).
            serialized_before = json.dumps(payload, indent=2, sort_keys=True)
            for handler in handlers:
                handler(tool, payload)
            serialized_after = json.dumps(payload, indent=2, sort_keys=True)
            if serialized_after != serialized_before:
                report_path.write_text(serialized_after + "\n", encoding="utf-8")
                _refresh_complex_portal_markdown(
                    report_path=report_path,
                    markdown_path=_resolve_existing_path(item.get("markdown_report"), batch_dir),
                    payload=payload,
                    modules=normalized_modules,
                )
                item["status"] = "ok"
                item["error"] = None
                updated += 1
        except Exception as exc:
            item["status"] = "error"
            item["error"] = str(exc)
    summary_file.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    return summary_file, results, updated


def _refresh_cross_reactivity_in_bulk(
    *,
    summary_file: Path,
    batch_dir: Path,
    results: list[dict[str, Any]],
    selected_queries: set[str],
    only_zero_cross_reactivity: bool,
    verbose: bool,
    analyzer: AntigenAnalyzer,
) -> tuple[Path, list[dict[str, Any]], int]:
    records: list[dict[str, Any]] = []
    blast_queries: list[dict[str, str]] = []

    for index, item in enumerate(results):
        query = str(item.get("query") or "").strip()
        resolved = str(item.get("resolved_entry_name") or "").strip()
        if selected_queries and query not in selected_queries and resolved not in selected_queries:
            continue
        report_path = _resolve_existing_path(item.get("json_report"), batch_dir)
        if report_path is None or not report_path.exists():
            continue
        try:
            payload = json.loads(report_path.read_text(encoding="utf-8"))
        except Exception as exc:
            item["status"] = "error"
            item["error"] = str(exc)
            continue
        if only_zero_cross_reactivity and not _has_zero_relevant_cross_reactivity(payload):
            continue

        target = payload.get("target") or {}
        sequence = str(target.get("sequence") or "")
        topology = payload.get("topology") or {}
        topology_class = str(topology.get("topology_class") or "")
        is_multipass = topology_class.startswith("multipass")
        run_full_length_cross_reactivity = is_multipass or payload.get("ectodomain") is None
        extracellular_surface_positions = _payload_extracellular_surface_positions(payload)
        record: dict[str, Any] = {
            "item": item,
            "payload": payload,
            "path": report_path,
            "query": query or resolved or report_path.name,
            "ecto_task_id": None,
            "ecto_sequence": "",
            "full_task_id": None,
            "full_sequence": "",
            "extracellular_surface_positions": extracellular_surface_positions,
        }

        if not sequence:
            payload["cross_reactivity_hits"] = []
            if run_full_length_cross_reactivity:
                payload["full_length_cross_reactivity_hits"] = []
            records.append(record)
            continue

        ectodomain = payload.get("ectodomain") or {}
        if ectodomain and ectodomain.get("start") is not None and ectodomain.get("end") is not None:
            query_sequence = slice_sequence(sequence, int(ectodomain["start"]), int(ectodomain["end"]))
            task_id = f"r{index}_ecto"
            record["ecto_task_id"] = task_id
            record["ecto_sequence"] = query_sequence
            blast_queries.append(
                {
                    "id": task_id,
                    "sequence": query_sequence,
                    "target_accession": str(target.get("accession") or ""),
                    "target_entry_name": str(target.get("entry_name") or ""),
                }
            )
        else:
            payload["cross_reactivity_hits"] = []

        if run_full_length_cross_reactivity:
            task_id = f"r{index}_full"
            record["full_task_id"] = task_id
            record["full_sequence"] = sequence
            blast_queries.append(
                {
                    "id": task_id,
                    "sequence": sequence,
                    "target_accession": str(target.get("accession") or ""),
                    "target_entry_name": str(target.get("entry_name") or ""),
                }
            )
        records.append(record)

    if verbose:
        print(
            f"[agdesign2-modules] bulk cross-reactivity: {len(blast_queries)} BLAST queries for {len(records)} reports",
            flush=True,
        )

    try:
        hits_by_query = analyzer.blast_client.search_many(blast_queries) if blast_queries else {}
    except Exception as exc:
        if verbose:
            print(f"[agdesign2-modules] bulk BLAST failed; falling back to per-report BLAST: {exc}", flush=True)
        return _refresh_cross_reactivity_sequential_fallback(
            summary_file=summary_file,
            results=results,
            records=records,
            analyzer=analyzer,
        )

    updated = 0
    for record in records:
        payload = record["payload"]
        ecto_task_id = record.get("ecto_task_id")
        if ecto_task_id:
            hits = hits_by_query.get(str(ecto_task_id), [])
            hits = analyzer._canonicalize_macaca_cross_reactivity_hits(str(record.get("ecto_sequence") or ""), hits)
            ectodomain = payload.get("ectodomain") or {}
            hits = analyzer._annotate_blast_hits_with_surface_identity(
                hits,
                query_region_start=int(ectodomain.get("start") or 1),
                extracellular_surface_positions=set(record.get("extracellular_surface_positions") or set()),
            )
            payload["cross_reactivity_hits"] = [hit.to_dict() for hit in hits]
            _remove_stale_cross_reactivity_failure_notes(payload)
        full_task_id = record.get("full_task_id")
        if full_task_id:
            hits = hits_by_query.get(str(full_task_id), [])
            hits = analyzer._canonicalize_macaca_cross_reactivity_hits(str(record.get("full_sequence") or ""), hits)
            hits = analyzer._annotate_blast_hits_with_surface_identity(
                hits,
                query_region_start=1,
                extracellular_surface_positions=set(record.get("extracellular_surface_positions") or set()),
            )
            payload["full_length_cross_reactivity_hits"] = [hit.to_dict() for hit in hits]
            _remove_stale_cross_reactivity_failure_notes(payload)
        record["path"].write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        record["item"]["status"] = "ok"
        record["item"]["error"] = None
        updated += 1

    summary_file.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    return summary_file, results, updated


def _refresh_cross_reactivity_sequential_fallback(
    *,
    summary_file: Path,
    results: list[dict[str, Any]],
    records: list[dict[str, Any]],
    analyzer: AntigenAnalyzer,
) -> tuple[Path, list[dict[str, Any]], int]:
    updated = 0
    for record in records:
        try:
            changed = _refresh_cross_reactivity(analyzer, record["payload"])
            if changed:
                record["path"].write_text(json.dumps(record["payload"], indent=2, sort_keys=True) + "\n", encoding="utf-8")
                record["item"]["status"] = "ok"
                record["item"]["error"] = None
                updated += 1
        except Exception as exc:
            record["item"]["status"] = "error"
            record["item"]["error"] = str(exc)
            notes = record["payload"].setdefault("notes", [])
            if isinstance(notes, list):
                notes.append(
                    {
                        "severity": "warning",
                        "message": f"BLAST sequence similarity refresh failed: {exc}",
                        "source": "BLAST",
                    }
                )
            record["path"].write_text(
                json.dumps(record["payload"], indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
    summary_file.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    return summary_file, results, updated


def _normalize_module_name(name: str) -> str:
    normalized = str(name or "").strip().lower().replace("_", "-")
    if normalized not in _MODULE_HANDLERS:
        valid = ", ".join(sorted(_MODULE_HANDLERS))
        raise ValueError(f"Unknown module '{name}'. Expected one of: {valid}.")
    return normalized


def _module_refresh_config(*, batch_dir: Path, verbose: bool, enable_complex_portal: bool) -> AnalysisConfig:
    snapshot_root = batch_dir.parent.parent if batch_dir.name == "surfy_batch" and batch_dir.parent.name == "outputs" else None
    cache_dir = snapshot_root / "cache" if snapshot_root is not None else Path(".agdesign2/cache")
    data_dir = snapshot_root / "data" if snapshot_root is not None else Path(".agdesign2/data")
    return AnalysisConfig(
        cache_dir=cache_dir,
        data_dir=data_dir,
        blast_db_dir=data_dir / "blastdb",
        ortholog_fasta_dir=data_dir / "proteomes",
        ortholog_blast_db_dir=data_dir / "ortholog_blastdb",
        precomputed_ortholog_table_path=batch_dir / "ortholog_reference_table.tsv",
        precomputed_family_alignment_index_path=batch_dir / "family_alignments" / "family_alignment_index.tsv",
        precomputed_paralog_index_path=batch_dir / "paralogs" / "paralog_index.tsv",
        verbose_progress=verbose,
        generate_assets=False,
        render_structure_images=False,
        render_quality_plots=False,
        enable_complex_portal=enable_complex_portal,
    )


def _refresh_cross_reactivity(analyzer: AntigenAnalyzer, payload: dict[str, Any]) -> bool:
    target = payload.get("target") or {}
    ectodomain = payload.get("ectodomain") or {}
    topology = payload.get("topology") or {}
    sequence = str(target.get("sequence") or "")
    topology_class = str(topology.get("topology_class") or "")
    is_multipass = topology_class.startswith("multipass")
    run_full_length_cross_reactivity = is_multipass or payload.get("ectodomain") is None
    if not sequence:
        payload["cross_reactivity_hits"] = []
        if run_full_length_cross_reactivity:
            payload["full_length_cross_reactivity_hits"] = []
        return True
    if ectodomain and ectodomain.get("start") is not None and ectodomain.get("end") is not None:
        query_sequence = slice_sequence(sequence, int(ectodomain["start"]), int(ectodomain["end"]))
        hits = analyzer.blast_client.search(
            query_sequence,
            target_accession=str(target.get("accession") or ""),
            target_entry_name=str(target.get("entry_name") or ""),
        )
        hits = analyzer._canonicalize_macaca_cross_reactivity_hits(query_sequence, hits)
        hits = analyzer._annotate_blast_hits_with_surface_identity(
            hits,
            query_region_start=int(ectodomain.get("start") or 1),
            extracellular_surface_positions=_payload_extracellular_surface_positions(payload),
        )
        payload["cross_reactivity_hits"] = [hit.to_dict() for hit in hits]
        _remove_stale_cross_reactivity_failure_notes(payload)
    else:
        payload["cross_reactivity_hits"] = []
    if run_full_length_cross_reactivity:
        full_length_hits = analyzer.blast_client.search(
            sequence,
            target_accession=str(target.get("accession") or ""),
            target_entry_name=str(target.get("entry_name") or ""),
        )
        full_length_hits = analyzer._canonicalize_macaca_cross_reactivity_hits(sequence, full_length_hits)
        full_length_hits = analyzer._annotate_blast_hits_with_surface_identity(
            full_length_hits,
            query_region_start=1,
            extracellular_surface_positions=_payload_extracellular_surface_positions(payload),
        )
        payload["full_length_cross_reactivity_hits"] = [hit.to_dict() for hit in full_length_hits]
        _remove_stale_cross_reactivity_failure_notes(payload)
    return True


def _remove_stale_cross_reactivity_failure_notes(payload: dict[str, Any]) -> None:
    notes = payload.get("notes")
    if not isinstance(notes, list):
        return
    payload["notes"] = [
        note
        for note in notes
        if not (
            isinstance(note, dict)
            and str(note.get("source") or "").upper() == "BLAST"
            and any(
                marker in str(note.get("message") or "").lower()
                for marker in (
                    "cross-reactivity blast search failed",
                    "blast sequence similarity search was not run",
                    "blast sequence similarity refresh failed",
                )
            )
        )
    ]


def _has_zero_relevant_cross_reactivity(payload: dict[str, Any]) -> bool:
    topology = payload.get("topology") or {}
    topology_class = str(topology.get("topology_class") or "")
    if topology_class.startswith("multipass") or payload.get("ectodomain") is None:
        return len(payload.get("full_length_cross_reactivity_hits") or []) == 0
    return len(payload.get("cross_reactivity_hits") or []) == 0


def _payload_extracellular_surface_positions(payload: dict[str, Any]) -> set[int]:
    positions: set[int] = set()
    for annotation in payload.get("residue_annotations") or []:
        if isinstance(annotation, dict):
            position = annotation.get("position")
            topology_location = annotation.get("topology_location")
            surface_exposed = annotation.get("surface_exposed")
        else:
            position = getattr(annotation, "position", None)
            topology_location = getattr(annotation, "topology_location", None)
            surface_exposed = getattr(annotation, "surface_exposed", None)
        if topology_location == "extracellular" and surface_exposed is True:
            try:
                positions.add(int(position))
            except (TypeError, ValueError):
                continue
    return positions


def _refresh_family_context(analyzer: AntigenAnalyzer, payload: dict[str, Any]) -> bool:
    target = _payload_to_target(payload)
    precomputed_ortholog = analyzer.precomputed_store.find_ortholog_record(target=target)
    canonical_family = analyzer.precomputed_store.canonical_family_from_record(precomputed_ortholog)
    if canonical_family is None:
        annotations = analyzer.interpro_client.fetch_annotations(target.accession)
        canonical_family = analyzer._select_canonical_family(
            annotations,
            sequence_length=len(target.sequence),
            protein_name=target.protein_name,
            gene_symbol=target.gene_symbol,
        )
    family_context = analyzer.precomputed_store.load_paralog_context(
        target=target,
        max_members=analyzer.config.max_paralog_context_members,
    )
    if family_context is None:
        family_context = analyzer.precomputed_store.load_family_context(
            target=target,
            canonical_family=canonical_family,
            max_members=analyzer.config.max_family_context_members,
        )
    if family_context is None and not analyzer._should_skip_live_family_context_fallback(
        precomputed_ortholog=precomputed_ortholog,
        canonical_family=canonical_family,
    ):
        family_context = analyzer.hgnc_client.fetch_family_context(target.gene_symbol)
    payload["canonical_family"] = canonical_family.to_dict() if canonical_family is not None else None
    payload["family_context"] = family_context.to_dict() if family_context is not None else None
    payload["notes"] = _rewrite_family_context_notes(
        payload.get("notes"),
        family_context=family_context,
        skipped_live_fallback=family_context is None
        and analyzer._should_skip_live_family_context_fallback(
            precomputed_ortholog=precomputed_ortholog,
            canonical_family=canonical_family,
        ),
    )
    return True


def _refresh_target_metadata(analyzer: AntigenAnalyzer, payload: dict[str, Any]) -> bool:
    target_payload = payload.get("target")
    if not isinstance(target_payload, dict):
        return False
    accession = str(target_payload.get("accession") or "").strip()
    entry_name = str(target_payload.get("entry_name") or "").strip()
    query = accession or entry_name
    if not query:
        return False

    entry = analyzer.uniprot_client._fetch_entry(query)
    refreshed = analyzer.uniprot_client._target_from_entry(entry)
    updates = refreshed.to_dict()
    changed = False
    for key in (
        "accession",
        "entry_name",
        "gene_symbol",
        "protein_name",
        "organism",
        "taxon_id",
        "canonical_isoform_id",
        "alternative_names",
        "gene_synonyms",
    ):
        value = updates.get(key)
        if target_payload.get(key) != value:
            target_payload[key] = value
            changed = True
    return changed


def _refresh_complex_portal(analyzer: AntigenAnalyzer, payload: dict[str, Any]) -> bool:
    """Refresh cached Complex Portal context without changing unrelated report fields."""
    target = _payload_to_target(payload)
    notes = _analysis_notes_without_complex_portal_failures(payload.get("notes"))
    homolog_context = _complex_portal_homolog_context(analyzer, payload, target)
    lookup_status: dict[str, dict[str, str]] = {}
    complexes = analyzer._fetch_complex_portal_context(
        target=target,
        homolog_context=homolog_context,
        notes=notes,
        lookup_status=lookup_status,
    )
    payload["complex_portal_complexes"] = [complex_item.to_dict() for complex_item in complexes]
    payload["complex_portal_lookup"] = lookup_status

    if analyzer._is_canonical_integrin_chain(target.gene_symbol):
        entry = analyzer.uniprot_client._fetch_entry(target.accession or target.entry_name)
        requirements = analyzer._collect_assembly_requirements(
            entry=entry,
            gene_symbol=target.gene_symbol,
            family_context=_payload_to_family_context(payload.get("family_context")),
            complex_portal_complexes=complexes,
            target_accession=target.accession,
        )
        payload["assembly_requirements"] = [requirement.to_dict() for requirement in requirements]
    else:
        payload["assembly_requirements"] = []
    payload["notes"] = [note.to_dict() for note in notes]
    return True


def _analysis_notes_without_complex_portal_failures(notes_payload: Any) -> list[AnalysisNote]:
    notes: list[AnalysisNote] = []
    for raw_note in notes_payload if isinstance(notes_payload, list) else []:
        if not isinstance(raw_note, dict):
            continue
        source = raw_note.get("source")
        message = str(raw_note.get("message") or "")
        if str(source or "").lower() == "complex portal" and "lookup failed:" in message.lower():
            continue
        notes.append(
            AnalysisNote(
                severity=str(raw_note.get("severity") or "info"),
                message=message,
                source=str(source) if source is not None else None,
            )
        )
    return notes


def _complex_portal_homolog_context(
    analyzer: AntigenAnalyzer,
    payload: dict[str, Any],
    target: TargetRecord,
) -> dict[str, dict[str, Any]]:
    if not analyzer._is_canonical_integrin_chain(target.gene_symbol):
        return {}
    precomputed_store = getattr(analyzer, "precomputed_store", None)
    if precomputed_store is not None:
        record = precomputed_store.find_ortholog_record(target=target)
        if record is not None:
            ectodomain = _payload_to_region(payload.get("ectodomain"))
            _, _, context = precomputed_store.build_species_matches(
                target=target,
                ectodomain=ectodomain,
                record=record,
            )
            return context
    for match in payload.get("species_name_matches") or []:
        if not isinstance(match, dict) or str(match.get("species") or "") != "mouse":
            continue
        accession = str(match.get("accession") or "").strip()
        if not accession:
            continue
        entry = analyzer.uniprot_client._fetch_entry(accession)
        return {"mouse": {"target": analyzer.uniprot_client._target_from_entry(entry)}}
    return {}


def _payload_to_region(payload: Any) -> Region | None:
    if not isinstance(payload, dict):
        return None
    try:
        return Region(
            start=int(payload["start"]),
            end=int(payload["end"]),
            label=str(payload.get("label") or ""),
            source=str(payload.get("source") or ""),
            confidence=payload.get("confidence"),
            metadata=dict(payload.get("metadata") or {}),
        )
    except (KeyError, TypeError, ValueError):
        return None


def _payload_to_family_context(payload: Any) -> FamilyContext | None:
    if not isinstance(payload, dict):
        return None
    members: list[FamilyMember] = []
    for raw_member in payload.get("members") or []:
        if not isinstance(raw_member, dict):
            continue
        members.append(
            FamilyMember(
                gene_symbol=str(raw_member.get("gene_symbol") or ""),
                gene_name=str(raw_member.get("gene_name") or ""),
                accession=raw_member.get("accession"),
                entry_name=raw_member.get("entry_name"),
                sequence_length=raw_member.get("sequence_length"),
                ectodomain_start=raw_member.get("ectodomain_start"),
                ectodomain_end=raw_member.get("ectodomain_end"),
                ectodomain_length=raw_member.get("ectodomain_length"),
                notes=list(raw_member.get("notes") or []),
            )
        )
    return FamilyContext(
        gene_symbol=str(payload.get("gene_symbol") or ""),
        source=str(payload.get("source") or ""),
        family_names=list(payload.get("family_names") or []),
        members=members,
        identity_matrix_labels=list(payload.get("identity_matrix_labels") or []),
        identity_matrix=list(payload.get("identity_matrix") or []),
        coverage_matrix_labels=list(payload.get("coverage_matrix_labels") or []),
        coverage_matrix=list(payload.get("coverage_matrix") or []),
        metadata=dict(payload.get("metadata") or {}),
    )


def _refresh_complex_portal_markdown(
    *,
    report_path: Path,
    markdown_path: Path | None,
    payload: dict[str, Any],
    modules: list[str],
) -> None:
    if "complex-portal" not in modules:
        return
    destination = markdown_path or report_path.with_suffix(".md")
    existing = destination.read_text(encoding="utf-8") if destination.exists() else ""
    replacement = _complex_portal_markdown_section(payload)
    pattern = r"\n## Complex Portal Context\n.*?(?=\n## |\Z)"
    if re.search(pattern, existing, flags=re.DOTALL):
        updated = re.sub(pattern, "\n" + replacement.rstrip(), existing, count=1, flags=re.DOTALL)
    else:
        updated = existing.rstrip() + "\n\n" + replacement
    destination.write_text(updated.rstrip() + "\n", encoding="utf-8")


def _complex_portal_markdown_section(payload: dict[str, Any]) -> str:
    lookup = payload.get("complex_portal_lookup") or {}
    human = lookup.get("human") if isinstance(lookup, dict) else None
    status = str((human or {}).get("status") or "not_recorded")
    status_text = {
        "ok": "records found.",
        "no_hits": "no curated records matched this target.",
        "error": "failed for this release.",
        "disabled": "disabled for this release.",
    }.get(status, "not recorded for this report.")
    lines = ["## Complex Portal Context", "", f"- Human lookup: {status_text}"]
    complexes = [item for item in (payload.get("complex_portal_complexes") or []) if isinstance(item, dict)]
    if complexes:
        lines.extend(["", "Curated Complex Portal records:", ""])
        for item in complexes:
            participants = ", ".join(
                str(participant.get("name") or participant.get("identifier") or "")
                for participant in (item.get("participants") or [])
                if isinstance(participant, dict)
            )
            lines.append(
                f"- `{item.get('complex_ac') or ''}` {item.get('name') or ''}; "
                f"{item.get('species') or 'n/a'}; participants: {participants or 'n/a'}"
            )
    return "\n".join(lines)


def _payload_to_target(payload: dict[str, Any]) -> TargetRecord:
    target = payload.get("target") or {}
    return TargetRecord(
        accession=str(target.get("accession") or ""),
        entry_name=str(target.get("entry_name") or ""),
        gene_symbol=target.get("gene_symbol"),
        protein_name=str(target.get("protein_name") or target.get("entry_name") or ""),
        organism=str(target.get("organism") or ""),
        taxon_id=target.get("taxon_id"),
        sequence=str(target.get("sequence") or ""),
        canonical_isoform_id=target.get("canonical_isoform_id"),
        alternative_names=list(target.get("alternative_names") or []),
        gene_synonyms=list(target.get("gene_synonyms") or []),
    )


def _resolve_existing_path(path_text: str | None, batch_dir: Path) -> Path | None:
    if not path_text:
        return None
    path = Path(path_text)
    candidates = [
        path,
        Path.cwd() / path,
        batch_dir / path,
        batch_dir / path.name,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return None


def _rewrite_family_context_notes(
    notes_payload: Any,
    *,
    family_context: Any,
    skipped_live_fallback: bool,
) -> list[dict[str, Any]]:
    notes = [note for note in (notes_payload if isinstance(notes_payload, list) else []) if isinstance(note, dict)]
    filtered = []
    for note in notes:
        message = str(note.get("message") or "")
        if (
            "family context loaded from precomputed family alignments" in message.lower()
            or "skipping live family-context fallback because this target is covered by precomputed reference tables" in message.lower()
            or "family context lookup failed:" in message.lower()
        ):
            continue
        filtered.append(note)
    if family_context is not None:
        source = str(getattr(family_context, "source", "") or "")
        if source.lower().startswith("precomputed"):
            message = (
                "Family context loaded from the precomputed paralog reference set."
                if "paralog" in source.lower()
                else "Family context loaded from precomputed family alignments."
            )
            filtered.append(
                {
                    "severity": "info",
                    "message": message,
                    "source": "precomputed",
                }
            )
    elif skipped_live_fallback:
        filtered.append(
            {
                "severity": "info",
                "message": "Skipping live family-context fallback because this target is covered by precomputed reference tables.",
                "source": "precomputed",
            }
        )
    return filtered


_MODULE_HANDLERS: dict[str, ModuleHandler] = {
    "complex-portal": _refresh_complex_portal,
    "cross-reactivity": _refresh_cross_reactivity,
    "family-context": _refresh_family_context,
    "target-metadata": _refresh_target_metadata,
}
