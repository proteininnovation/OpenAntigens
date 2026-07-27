from __future__ import annotations

import re

from .models import (
    AnalysisReport,
    AssemblyRequirement,
    BlastHit,
    ComplexPortalComplex,
    ConstructDetail,
    ConstructSuggestion,
    GPCREngineeringVariant,
    HomologyRecord,
)


def render_markdown_report(report: AnalysisReport) -> str:
    lines: list[str] = []
    lines.append(f"# {report.target.entry_name} Antigen Design Report")
    lines.append("")
    lines.append("## Target")
    lines.append("")
    lines.append(f"- Query: `{report.resolution.get('query', '')}`")
    lines.append(f"- Resolved accession: `{report.target.accession}`")
    lines.append(f"- Entry name: `{report.target.entry_name}`")
    lines.append(f"- Gene symbol: `{report.target.gene_symbol or 'n/a'}`")
    lines.append(f"- Protein: {report.target.protein_name}")
    lines.append(f"- Organism: {report.target.organism}")
    lines.append(f"- Sequence length: {len(report.target.sequence)} aa")
    lines.append("")
    lines.append("## Ectodomain")
    lines.append("")
    if report.ectodomain:
        lines.append(
            f"- Boundary: {report.ectodomain.start}-{report.ectodomain.end} "
            f"({report.ectodomain.end - report.ectodomain.start + 1} aa)"
        )
        lines.append(f"- Assignment: {report.ectodomain.label}")
        lines.append(f"- Source: {report.ectodomain.source}")
        if report.ectodomain.confidence is not None:
            lines.append(f"- Confidence: {report.ectodomain.confidence:.2f}")
    else:
        lines.append("- Boundary: not determined")
    lines.append("")
    lines.append("## Construct Recommendations")
    lines.append("")
    if report.construct_recommendations:
        soluble_constructs, membrane_constructs = _partition_constructs(report.construct_recommendations)
        if soluble_constructs:
            lines.append("### Soluble ECD Constructs")
            lines.append("")
            lines.extend(_render_constructs_table(soluble_constructs))
        if membrane_constructs:
            if soluble_constructs:
                lines.append("")
            lines.append("### Native Membrane-Expression Constructs")
            lines.append("")
            lines.extend(_render_constructs_table(membrane_constructs))
    else:
        lines.append("- No construct recommendations available.")
    if report.construct_details:
        lines.append("")
        lines.append("Detailed construct sections appear later in this report.")
    lines.append("")
    lines.append("## Structural Support")
    lines.append("")
    if report.structural_regions:
        lines.append("AlphaFold regions:")
        lines.append("")
        lines.extend(_render_structural_regions_table(report))
        split_diagnostics = _collect_split_diagnostics(report)
        if split_diagnostics:
            lines.append("")
            lines.append("PAE split diagnostics:")
            lines.append("")
            lines.extend(_render_split_diagnostics_table(split_diagnostics))
    else:
        lines.append("- No AlphaFold structural regions available.")
    if report.experimental_constructs:
        lines.append("")
        lines.append("PDB-backed constructs:")
        lines.append("")
        lines.extend(_render_pdb_constructs_table(report))
    if report.gpcr_annotation:
        lines.append("")
        lines.append("## GPCR Annotation")
        lines.append("")
        lines.extend(_render_gpcr_annotation(report))
    lines.append("")
    lines.append("## InterPro / Pfam Domain Architecture")
    lines.append("")
    if report.interpro_annotations:
        lines.extend(_render_domain_architecture_table(report))
    else:
        lines.append("- No InterPro/Pfam annotations available.")
    lines.append("")
    lines.append("## Canonical Family")
    lines.append("")
    if report.canonical_family:
        lines.extend(_render_canonical_family_table(report))
    else:
        lines.append("- No canonical InterPro family assignment available.")
    lines.append("")
    lines.append("## Ligand Interaction Regions")
    lines.append("")
    if report.ligand_interactions:
        lines.extend(_render_ligand_table(report))
    else:
        lines.append("- No curated ectodomain ligand-interaction annotations found.")
    lines.append("")
    lines.append("## Sequence Liabilities")
    lines.append("")
    lines.append("Post-translational modifications and processing annotations:")
    lines.append("")
    if report.ptms:
        lines.extend(_render_ptm_table(report.ptms))
    else:
        lines.append("- No curated UniProt PTM annotations available.")
    lines.append("")
    if report.furin_sites:
        lines.append("Furin sites:")
        for site in report.furin_sites:
            suggestions = ", ".join(site.suggestions) if site.suggestions else "none"
            lines.append(f"- {site.start}-{site.end} `{site.motif}`; suggested edits: {suggestions}")
    else:
        lines.append("- No furin motifs found in the ectodomain.")
    if report.cysteine_analysis:
        lines.append("")
        lines.append("Cysteines:")
        lines.append("")
        lines.extend(_render_cysteine_table(report))
    else:
        lines.append("- No cysteine analysis available.")
    lines.append("")
    lines.append("## Homology")
    lines.append("")
    if report.ectodomain_homology:
        lines.extend(_render_homology_table(report.ectodomain_homology))
    else:
        lines.append("- No homolog records available.")
    lines.append("")
    lines.append("## Family Context")
    lines.append("")
    if report.family_context and report.family_context.family_names:
        lines.extend(_render_family_context_tables(report))
    elif report.family_context:
        lines.append("- Family context retrieved, but no family names were returned.")
    else:
        lines.append("- No family context available.")
    lines.append("")
    lines.append("## Assembly / Partner Requirements")
    lines.append("")
    if report.assembly_requirements:
        lines.extend(_render_assembly_requirements_summary(report.assembly_requirements))
    else:
        lines.append("- No curated obligatory partner requirement identified.")
    lines.append("")
    lines.append("## Complex Portal Context")
    lines.append("")
    lines.extend(_render_complex_portal_lookup_status(report.complex_portal_lookup))
    if report.complex_portal_complexes:
        lines.append("")
        lines.append("Curated Complex Portal records:")
        lines.append("")
        lines.extend(_render_complex_portal_complexes_table(report.complex_portal_complexes))
    lines.append("")
    if report.gpcr_engineering_variants:
        lines.append("## GPCR Expression and Stabilization Strategy")
        lines.append("")
        lines.extend(_render_gpcr_engineering_variants(report.gpcr_engineering_variants))
        lines.append("")
    lines.append("## Advanced Membrane Engineering Suggestions")
    lines.append("")
    if report.advanced_membrane_suggestions:
        lines.extend(_render_advanced_membrane_suggestions(report.advanced_membrane_suggestions))
    else:
        lines.append("- No advanced membrane-engineering suggestions available.")
    lines.append("")
    lines.append("## Cross-Reactivity Hits")
    lines.append("")
    if report.cross_reactivity_hits:
        lines.extend(_render_blast_hits_table(report.cross_reactivity_hits))
    else:
        lines.append("- No ranked BLAST hits available.")
    if report.notes:
        lines.append("")
        lines.append("## Notes")
        lines.append("")
        for note in report.notes:
            source = f" ({note.source})" if note.source else ""
            lines.append(f"- {note.severity.upper()}{source}: {note.message}")
    if report.construct_details:
        lines.append("")
        lines.append("## Construct Details")
        lines.append("")
        soluble_details, membrane_details = _partition_constructs(report.construct_details)
        if soluble_details:
            lines.append("### Soluble ECD Constructs")
            lines.append("")
            for construct in soluble_details:
                lines.extend(_render_construct_detail(construct))
        if membrane_details:
            if soluble_details:
                lines.append("")
            lines.append("### Native Membrane-Expression Constructs")
            lines.append("")
            for construct in membrane_details:
                lines.extend(_render_construct_detail(construct))
    lines.append("")
    return "\n".join(lines)


def _render_constructs_table(constructs: list[ConstructSuggestion]) -> list[str]:
    rows = []
    for construct in constructs:
        rows.append(
            [
                _construct_name_md(construct.name, construct.pdb_id),
                str(construct.start),
                str(construct.end),
                str(construct.end - construct.start + 1),
                f"{construct.score:.1f}",
                construct.classification or "",
                construct.rationale,
            ]
        )
    return _render_table(
        ["Construct", "Start", "End", "Length", "Score", "Classification", "Rationale"],
        rows,
    )


def _partition_constructs(constructs):
    soluble = []
    membrane = []
    for construct in constructs:
        classification = getattr(construct, "classification", None)
        if classification == "membrane_expression":
            membrane.append(construct)
        else:
            soluble.append(construct)
    return soluble, membrane


def _render_structural_regions_table(report: AnalysisReport) -> list[str]:
    rows = []
    for region in report.structural_regions:
        rows.append(
            [
                region.label,
                str(region.metadata.get("stringency", "")),
                str(region.start),
                str(region.end),
                str(region.end - region.start + 1),
                region.source,
                f"{region.confidence:.2f}" if region.confidence is not None else "",
                (
                    f"{float(region.metadata.get('mean_intra_pae')):.2f}"
                    if region.metadata.get("mean_intra_pae") is not None
                    else ""
                ),
                ", ".join(region.metadata.get("included_interpro_domains", [])),
            ]
        )
    return _render_table(
        [
            "Region",
            "Stringency",
            "Start",
            "End",
            "Length",
            "Source",
            "Confidence",
            "Mean Intra-domain PAE",
            "Included InterPro Domains",
        ],
        rows,
    )


def _render_split_diagnostics_table(split_diagnostics: list[dict[str, float | int | None]]) -> list[str]:
    rows = []
    for split in split_diagnostics:
        rows.append(
            [
                f"{split['parent_start']}-{split['parent_end']}",
                str(split["boundary_after"]),
                f"{split['score']:.2f}",
                f"{split['inter_block_pae']:.2f}",
                f"{split['intra_block_pae']:.2f}",
                f"{split['separation']:.2f}",
                f"{split['boundary_pae']:.2f}",
                f"{split['linker_plddt']:.2f}" if split["linker_plddt"] is not None else "",
            ]
        )
    return _render_table(
        [
            "Parent Region",
            "Boundary After",
            "Score",
            "Inter-block PAE",
            "Intra-block PAE",
            "Separation",
            "Local Boundary PAE",
            "Linker pLDDT",
        ],
        rows,
    )


def _collect_split_diagnostics(report: AnalysisReport) -> list[dict[str, float | int | None]]:
    unique: dict[tuple[int, int, int], dict[str, float | int | None]] = {}
    for region in report.structural_regions:
        if region.label != "Predicted domain":
            continue
        for split in region.metadata.get("split_history", []):
            boundary_after = split.get("boundary_after")
            parent_start = split.get("parent_start")
            parent_end = split.get("parent_end")
            if boundary_after is None or parent_start is None or parent_end is None:
                continue
            key = (int(parent_start), int(boundary_after), int(parent_end))
            unique[key] = split
    return [unique[key] for key in sorted(unique)]


def _render_homology_table(records: list[HomologyRecord]) -> list[str]:
    rows = []
    for record in records:
        if not record.available:
            rows.append([record.species, "unavailable", "", "", "", "", " ".join(record.notes)])
            continue
        rows.append(
            [
                record.species,
                f"`{record.entry_name}`" if record.entry_name else "",
                str(record.ectodomain_start) if record.ectodomain_start is not None else "",
                str(record.ectodomain_end) if record.ectodomain_end is not None else "",
                f"{record.identity:.2f}%" if record.identity is not None else "",
                f"{record.coverage:.2f}%" if record.coverage is not None else "",
                " ".join(record.notes),
            ]
        )
    return _render_table(
        ["Species", "Entry", "Ectodomain Start", "Ectodomain End", "Identity", "Coverage", "Notes"],
        rows,
    )


def _render_blast_hits_table(hits: list[BlastHit]) -> list[str]:
    rows = []
    for index, hit in enumerate(hits, start=1):
        gene_symbol = _extract_blast_gene_symbol(hit)
        rows.append(
            [
                str(index),
                f"`{hit.subject_id}`",
                gene_symbol or "",
                hit.species or "unknown species",
                f"{hit.bitscore:.1f}",
                f"{hit.identity:.1f}%",
                f"{hit.coverage:.1f}%",
                f"{hit.evalue:.2e}",
            ]
        )
    return _render_table(["Rank", "Hit", "Gene", "Species", "Bit score", "Identity", "Coverage", "E-value"], rows)


def _render_assembly_requirements_table(requirements: list[AssemblyRequirement]) -> list[str]:
    rows = []
    for requirement in requirements:
        rows.append(
            [
                requirement.classification or "",
                "yes" if requirement.obligatory else "no",
                requirement.confidence,
                requirement.summary,
                ", ".join(requirement.partners) if requirement.partners else "",
                requirement.source,
                " ".join(requirement.rationale),
            ]
        )
    return _render_table(
        ["Classification", "Obligatory", "Confidence", "Summary", "Named Partner(s)", "Source", "Rationale"],
        rows,
    )


def _render_complex_portal_lookup_status(lookup: dict[str, dict[str, str]]) -> list[str]:
    human = lookup.get("human") if isinstance(lookup, dict) else None
    status = str((human or {}).get("status") or "not_recorded")
    if status == "ok":
        return ["- Human lookup: records found."]
    if status == "no_hits":
        return ["- Human lookup: no curated records matched this target."]
    if status == "error":
        return ["- Human lookup: failed for this release."]
    if status == "disabled":
        return ["- Human lookup: disabled for this release."]
    return ["- Human lookup: not recorded for this report."]


def _render_complex_portal_complexes_table(complexes: list[ComplexPortalComplex]) -> list[str]:
    rows = []
    for complex_item in complexes:
        partner_names = [
            participant.name
            for participant in complex_item.participants
            if (participant.interactor_type or "").lower() == "protein"
        ]
        rows.append(
            [
                f"`{complex_item.complex_ac}` {complex_item.name}",
                str(complex_item.confidence_score) if complex_item.confidence_score is not None else "",
                complex_item.evidence_code or "",
                ", ".join(complex_item.complex_assemblies),
                ", ".join(partner_names[:6]),
            ]
        )
    return _render_table(
        ["Complex", "Stars", "Evidence", "Assembly", "Protein Participants"],
        rows,
    )


def _render_assembly_requirements_summary(requirements: list[AssemblyRequirement]) -> list[str]:
    summarized = [
        requirement
        for requirement in requirements
        if requirement.classification != "interaction_context"
    ]
    omitted_count = len(requirements) - len(summarized)
    if not summarized:
        lines = ["- No curated obligatory or conditional assembly requirement identified."]
        if omitted_count:
            lines.append(f"- Additional interaction-context annotations omitted from summary: {omitted_count}")
        return lines
    lines = _render_assembly_requirements_table(summarized)
    if omitted_count:
        lines.append("")
        lines.append(f"- Additional interaction-context annotations omitted from summary: {omitted_count}")
    return lines


def _render_advanced_membrane_suggestions(suggestions) -> list[str]:
    lines: list[str] = []
    for suggestion in suggestions:
        lines.append(f"### {suggestion.title}")
        lines.append("")
        lines.append(f"- Category: {suggestion.category}")
        if suggestion.start is not None and suggestion.end is not None:
            lines.append(f"- Region: {suggestion.start}-{suggestion.end}")
        lines.append(f"- Summary: {suggestion.summary}")
        if suggestion.rationale:
            lines.append("- Rationale: " + " ".join(suggestion.rationale))
        if suggestion.suggested_actions:
            lines.append("- Suggested actions:")
            for action in suggestion.suggested_actions:
                lines.append(f"- Action: {action}")
        if suggestion.evidence:
            lines.append("- Evidence: " + "; ".join(suggestion.evidence))
        if suggestion.warnings:
            lines.append("- Warnings: " + "; ".join(suggestion.warnings))
        lines.append("")
    return lines


def _render_gpcr_engineering_variants(variants: list[GPCREngineeringVariant]) -> list[str]:
    lines: list[str] = []
    lines.extend(
        [
            "| Variant | Category | Region | Cassette | Sequence role | Length |",
            "|---|---|---:|---|---|---:|",
        ]
    )
    for variant in variants:
        region = (
            f"{variant.start}-{variant.end}"
            if variant.start is not None and variant.end is not None
            else "n/a"
        )
        length = str(variant.engineered_length) if variant.engineered_length is not None else ""
        lines.append(
            "| {name} | {category} | {region} | {cassette} | {role} | {length} |".format(
                name=_md_code(variant.name),
                category=variant.category,
                region=region,
                cassette=variant.cassette_name or variant.cassette_id or "n/a",
                role=variant.sequence_role,
                length=length,
            )
        )
    lines.append("")
    for variant in variants:
        lines.append(f"### `{variant.name}`")
        lines.append("")
        lines.append(f"- Strategy: {variant.strategy}")
        lines.append(f"- Summary: {variant.summary}")
        if variant.start is not None and variant.end is not None:
            lines.append(f"- Focus region: {variant.start}-{variant.end}")
        if variant.replaced_start is not None and variant.replaced_end is not None:
            lines.append(f"- Replaced receptor residues: {variant.replaced_start}-{variant.replaced_end}")
        if variant.n_linker is not None or variant.c_linker is not None:
            lines.append(f"- Linkers: N `{variant.n_linker or ''}`; C `{variant.c_linker or ''}`")
        if variant.gpcr_generic_range:
            lines.append(f"- GPCRdb generic range: {variant.gpcr_generic_range}")
        if variant.evidence:
            lines.append("- Evidence: " + "; ".join(variant.evidence))
        if variant.suggested_actions:
            lines.append("- Suggested actions: " + "; ".join(variant.suggested_actions))
        if variant.warnings:
            lines.append("- Warnings: " + "; ".join(variant.warnings))
        if variant.sequence:
            lines.append("")
            lines.append("```fasta")
            lines.append(f">{variant.name}")
            lines.extend(_wrap_sequence(variant.sequence))
            lines.append("```")
        lines.append("")
    return lines


def _md_code(text: str) -> str:
    return "`" + text.replace("`", "") + "`"


def _construct_name_md(name: str, pdb_id: str | None) -> str:
    """Code-formatted construct name, linked to RCSB when PDB-backed."""
    code = _md_code(name)
    if pdb_id:
        return f"[{code}](https://www.rcsb.org/structure/{pdb_id.strip().upper()})"
    return code


def _extract_blast_gene_symbol(hit: BlastHit) -> str | None:
    for text in (hit.description, hit.subject_id):
        match = re.search(r"\bGN=([A-Za-z0-9_-]+)\b", text)
        if match:
            return match.group(1)
    parts = hit.subject_id.split("|")
    if len(parts) >= 3 and "_" in parts[2]:
        return parts[2].split("_", 1)[0]
    return None


def _render_construct_detail(construct: ConstructDetail) -> list[str]:
    lines: list[str] = []
    lines.append(f"### {_construct_name_md(construct.name, construct.pdb_id)}")
    lines.append("")
    lines.append(f"- human: `{construct.human_entry_name or 'HUMAN'}` {construct.start}-{construct.end}")
    lines.append(f"- Human length: {construct.length} aa")
    lines.append(f"- Score: {construct.score:.1f}")
    if construct.classification:
        lines.append(f"- Classification: {construct.classification}")
    lines.append(f"- Rationale: {construct.rationale}")
    if construct.domain_summary:
        lines.append(f"- Architecture summary: {construct.domain_summary}")
    if construct.complete_domain_annotations:
        lines.append(
            "- Complete domain hits: "
            + "; ".join(_format_domain_annotation(annotation) for annotation in construct.complete_domain_annotations)
        )
    if construct.clipped_domain_annotations:
        lines.append(
            "- Clipped domain hits: "
            + "; ".join(_format_domain_annotation(annotation) for annotation in construct.clipped_domain_annotations)
        )
    if construct.structural_metrics:
        metrics = construct.structural_metrics
        if metrics.mean_plddt is not None:
            lines.append(
                f"- AlphaFold quality: mean pLDDT {metrics.mean_plddt:.2f}; "
                f"min {metrics.min_plddt:.2f}; max {metrics.max_plddt:.2f}"
            )
        if metrics.structured_fraction is not None:
            lines.append(
                f"- Structured coverage: {metrics.structured_fraction * 100:.2f}% residues at or above pLDDT 70"
            )
        if metrics.mean_intra_pae is not None:
            lines.append(
                f"- Intra-construct PAE: mean {metrics.mean_intra_pae:.2f}; max {metrics.max_intra_pae:.2f}"
            )
    if construct.evidence:
        lines.append(f"- Evidence: {', '.join(construct.evidence)}")
    if construct.warnings:
        lines.append(f"- Warnings: {', '.join(construct.warnings)}")
    if construct.cysteine_analysis:
        unpaired = [finding.position for finding in construct.cysteine_analysis if finding.paired_with is None]
        if unpaired:
            lines.append(
                "- Cysteine warning: "
                + ", ".join(f"C{position}" for position in unpaired)
                + " unpaired in this construct"
            )
    if construct.ligand_interactions:
        lines.append(
            "- Ligand-interaction overlap: "
            + "; ".join(
                f"{feature.start}-{feature.end} {feature.metadata.get('interaction_summary') or feature.description or feature.type}"
                for feature in construct.ligand_interactions
            )
        )
    if construct.ptms:
        lines.append(
            "- Included PTMs: "
            + "; ".join(
                f"{_ptm_label(feature)} {feature.metadata.get('ptm_category') or feature.type}: {feature.description or feature.type}"
                for feature in construct.ptms
            )
        )
    if construct.split_diagnostics:
        lines.append("- PAE split diagnostics:")
        lines.append("")
        lines.extend(
            _render_table(
                [
                    "Parent Region",
                    "Boundary After",
                    "Score",
                    "Inter-block PAE",
                    "Intra-block PAE",
                    "Separation",
                    "Local Boundary PAE",
                    "Linker pLDDT",
                ],
                [
                    [
                        f"{split.parent_start}-{split.parent_end}",
                        str(split.boundary_after),
                        f"{split.score:.2f}",
                        f"{split.inter_block_pae:.2f}",
                        f"{split.intra_block_pae:.2f}",
                        f"{split.separation:.2f}",
                        f"{split.boundary_pae:.2f}",
                        f"{split.linker_plddt:.2f}" if split.linker_plddt is not None else "",
                    ]
                    for split in construct.split_diagnostics
                ],
            )
        )
    else:
        lines.append("- PAE split diagnostics: no confident internal boundary identified for this construct")
    lines.append("")
    lines.append("Human sequence:")
    lines.append("")
    lines.append("```text")
    lines.extend(_wrap_sequence(construct.sequence))
    lines.append("```")
    lines.append("")
    if construct.homologs:
        lines.append("Homolog equivalents:")
        lines.append("")
        for homolog in construct.homologs:
            if homolog.available:
                line = f"- {homolog.species}: `{homolog.entry_name}` {homolog.start}-{homolog.end}"
                if homolog.identity_to_human is not None:
                    line += f"; identity to human {homolog.identity_to_human:.2f}%"
                lines.append(line)
                lines.append("")
                lines.append("```text")
                lines.extend(_wrap_sequence(homolog.sequence or ""))
                lines.append("```")
                if homolog.notes:
                    lines.append(f"Notes: {' '.join(homolog.notes)}")
            else:
                note_text = " ".join(homolog.notes) if homolog.notes else "Unavailable."
                lines.append(f"- {homolog.species}: unavailable. {note_text}".strip())
            lines.append("")
    if construct.structure_image:
        lines.append(f"![{construct.name} structure]({construct.structure_image})")
        lines.append("")
    if construct.quality_plot:
        lines.append(f"![{construct.name} quality plot]({construct.quality_plot})")
        lines.append("")
    return lines


def _wrap_sequence(sequence: str, width: int = 80) -> list[str]:
    if not sequence:
        return ["<no sequence available>"]
    return [sequence[index : index + width] for index in range(0, len(sequence), width)]


def _render_family_members(report: AnalysisReport) -> list[str]:
    family = report.family_context
    if family is None:
        return []
    rows = []
    for member in family.members:
        rows.append(
            [
                member.gene_symbol,
                f"`{member.accession}`" if member.accession else "",
                f"`{member.entry_name}`" if member.entry_name else "",
                str(member.ectodomain_start) if member.ectodomain_start is not None else "",
                str(member.ectodomain_end) if member.ectodomain_end is not None else "",
                str(member.ectodomain_length) if member.ectodomain_length is not None else "",
                " ".join(member.notes),
            ]
        )
    return _render_table(
        ["Gene", "Accession", "Entry", "Ectodomain Start", "Ectodomain End", "Ectodomain Length", "Notes"],
        rows,
    )


def _render_pdb_constructs_table(report: AnalysisReport) -> list[str]:
    deduped_constructs = _dedupe_experimental_constructs(report)
    rows = []
    for construct in deduped_constructs:
        chain_text = ", ".join(f"{chain.label} {chain.start}-{chain.end}" for chain in construct.chains)
        rows.append(
            [
                f"`{construct.pdb_id}`",
                chain_text,
                construct.method or "unknown method",
                construct.resolution or "",
            ]
        )
    return _render_table(["PDB", "Chains", "Method", "Resolution"], rows)


def _render_ligand_table(report: AnalysisReport) -> list[str]:
    rows = []
    for feature in report.ligand_interactions:
        rows.append(
            [
                str(feature.start),
                str(feature.end),
                feature.type,
                feature.metadata.get("interaction_summary") or feature.description or feature.type,
            ]
        )
    return _render_table(["Start", "End", "Type", "Summary"], rows)


def _render_ptm_table(ptms: list[Feature]) -> list[str]:
    rows = []
    for feature in ptms:
        rows.append(
            [
                str(feature.start),
                str(feature.end),
                feature.type,
                str(feature.metadata.get("ptm_category") or ""),
                feature.description or "",
            ]
        )
    return _render_table(["Start", "End", "Type", "Category", "Description"], rows)


def _ptm_label(feature: Feature) -> str:
    if feature.type.upper() == "DISULFID":
        return f"C{feature.start}-C{feature.end}"
    return f"{feature.start}-{feature.end}"


def _render_domain_architecture_table(report: AnalysisReport) -> list[str]:
    rows = []
    for annotation in report.interpro_annotations:
        rows.append(
            [
                str(annotation.start),
                str(annotation.end),
                str(annotation.end - annotation.start + 1),
                annotation.source_database,
                f"`{annotation.accession}`",
                annotation.type,
                annotation.name,
                "yes" if annotation.representative else "",
                (
                    f"`{annotation.integrated_accession}` {annotation.integrated_name or ''}".strip()
                    if annotation.integrated_accession
                    else ""
                ),
            ]
        )
    return _render_table(
        ["Start", "End", "Length", "Source", "Accession", "Type", "Name", "Representative", "Integrated"],
        rows,
    )


def _render_cysteine_table(report: AnalysisReport) -> list[str]:
    by_position = {finding.position: finding for finding in report.cysteine_analysis}
    seen: set[int] = set()
    rows = []
    for position in sorted(by_position):
        if position in seen:
            continue
        finding = by_position[position]
        if finding.paired_with is not None and finding.paired_with in by_position:
            partner = by_position[finding.paired_with]
            seen.add(position)
            seen.add(partner.position)
            rows.append(
                [
                    f"C{min(position, partner.position)}-C{max(position, partner.position)}",
                    "paired",
                    f"{_exposure_text(finding.surface_exposed)} / {_exposure_text(partner.surface_exposed)}",
                    "",
                    "; ".join(
                        warning
                        for warning in [finding.warning, partner.warning]
                        if warning
                    ),
                ]
            )
        else:
            seen.add(position)
            rows.append(
                [
                    f"C{position}",
                    "unpaired",
                    _exposure_text(finding.surface_exposed),
                    _closest_unpaired_cysteine_text(finding),
                    finding.warning or "",
                ]
            )
    return _render_table(["Residue(s)", "Status", "Exposure", "Closest unpaired cysteine", "Warning"], rows)


def _closest_unpaired_cysteine_text(finding: CysteineFinding) -> str:
    if finding.paired_with is not None:
        return ""
    if finding.closest_unpaired_cysteine is None or finding.closest_unpaired_cysteine_distance is None:
        return "none detected"
    atom = finding.closest_unpaired_cysteine_distance_atom or "distance"
    return f"C{finding.closest_unpaired_cysteine}, {finding.closest_unpaired_cysteine_distance:.2f} A ({atom})"


def _render_family_context_tables(report: AnalysisReport) -> list[str]:
    family = report.family_context
    if family is None:
        return []
    lines = _render_table(
        ["Family Group", "Source", "Member Count"],
        [[", ".join(family.family_names), family.source, str(len(family.members))]],
    )
    if family.members:
        lines.append("")
        lines.append("Family members:")
        lines.append("")
        lines.extend(_render_family_members(report))
    if family.identity_matrix_labels and family.identity_matrix:
        lines.append("")
        lines.append("Family ectodomain directional identity matrix (%):")
        lines.append("")
        lines.append("Rows are normalized to the row protein ectodomain length, so A->B can differ from B->A.")
        lines.append("")
        lines.extend(_render_matrix_table(family.identity_matrix_labels, family.identity_matrix))
    if family.coverage_matrix_labels and family.coverage_matrix:
        lines.append("")
        lines.append("Family ectodomain directional coverage matrix (%):")
        lines.append("")
        lines.append("Rows report what fraction of the row ectodomain is aligned to the column ectodomain.")
        lines.append("")
        lines.extend(_render_matrix_table(family.coverage_matrix_labels, family.coverage_matrix))
    return lines


def _render_canonical_family_table(report: AnalysisReport) -> list[str]:
    family = report.canonical_family
    if family is None:
        return []
    rows = [[
        family.accession,
        family.name,
        family.source_database,
        str(family.fragment_count),
        f"{family.start}-{family.end}" if family.start is not None and family.end is not None else "n/a",
        f"{family.coverage_fraction:.2f}%" if family.coverage_fraction is not None else "n/a",
    ]]
    lines = _render_table(
        ["Accession", "Name", "Source", "Fragments", "Span", "Coverage"],
        rows,
    )
    if family.notes:
        lines.append("")
        for note in family.notes:
            lines.append(f"- {note}")
    return lines


def _render_gpcr_annotation(report: AnalysisReport) -> list[str]:
    annotation = report.gpcr_annotation
    if annotation is None:
        return []
    lines = [
        f"- GPCRdb entry: `{annotation.entry_name}`",
        f"- GPCRdb name: {annotation.name or 'n/a'}",
        f"- Family: {' > '.join(annotation.family_path) if annotation.family_path else annotation.family_name or 'n/a'}",
        f"- Numbering scheme: {annotation.residue_numbering_scheme or 'n/a'}",
        f"- URL: {annotation.url or 'n/a'}",
        "",
        "Segments:",
        "",
    ]
    lines.extend(
        _render_table(
            ["Segment", "Boundary", "Generic start", "Generic end"],
            [
                [
                    segment.name,
                    f"{segment.start}-{segment.end}",
                    segment.generic_start or "",
                    segment.generic_end or "",
                ]
                for segment in annotation.segments
            ],
        )
        or ["- No GPCRdb segment mapping available."]
    )
    lines.append("")
    lines.append("Conserved motifs:")
    lines.append("")
    lines.extend(
        _render_table(
            ["Motif", "Residues", "Generic numbers", "Sequence", "Status"],
            [
                [
                    motif.name,
                    ", ".join(str(item) for item in motif.positions),
                    ", ".join(motif.generic_numbers),
                    motif.sequence or "",
                    motif.status,
                ]
                for motif in annotation.conserved_motifs
            ],
        )
        or ["- No conserved motif positions detected from GPCRdb numbering."]
    )
    return lines


def _dedupe_experimental_constructs(report: AnalysisReport):
    deduped = []
    seen: set[tuple[tuple[int, int], ...]] = set()
    for construct in report.experimental_constructs:
        signature = tuple(sorted((chain.start, chain.end) for chain in construct.chains))
        if signature in seen:
            continue
        seen.add(signature)
        deduped.append(construct)
    return deduped


def _render_table(headers: list[str], rows: list[list[str]]) -> list[str]:
    if not rows:
        return []
    safe_headers = [_escape_table_cell(header) for header in headers]
    lines = [
        "| " + " | ".join(safe_headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        safe_row = [_escape_table_cell(str(cell)) for cell in row]
        lines.append("| " + " | ".join(safe_row) + " |")
    return lines


def _escape_table_cell(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ").strip()


def _exposure_text(surface_exposed: bool | None) -> str:
    if surface_exposed is True:
        return "surface-exposed"
    if surface_exposed is False:
        return "buried"
    return "unknown exposure"


def _render_matrix_table(labels: list[str], matrix: list[list[float | None]]) -> list[str]:
    header = ["| Protein |", *[f" {label} |" for label in labels]]
    separator = ["| --- |", *[" ---: |" for _ in labels]]
    lines = ["".join(header), "".join(separator)]
    for label, row in zip(labels, matrix, strict=True):
        values = [f" {value:.2f} |" if value is not None else " n/a |" for value in row]
        lines.append("".join([f"| {label} |", *values]))
    return lines


def _format_domain_annotation(annotation) -> str:
    return (
        f"{annotation.source_database} {annotation.accession} {annotation.name} "
        f"{annotation.start}-{annotation.end}"
    )
